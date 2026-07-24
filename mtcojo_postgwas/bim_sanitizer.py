"""
bim_sanitizer.py - PLINK Reference Panel Validator & Sanitizer for mtcojo_postgwas

Responsibilities:
  1. Detect variant ID format in the PLINK BIM file (rsID or CHROM_POS_REF_ALT).
  2. Detect variant ID format in the LD score files (*.l2.ldscore.gz, SNP column).
  3. If BIM format ≠ LD score format → ABORT. IDs must match for GCTA LDSC to work.
  4. Convert non-standard 7-column BIM → standard 6-column BIM using the detected
     format column. IDs are NEVER changed — only the column layout is fixed.
  5. Symlink .bed and .fam files alongside the sanitized BIM.

Returns:
  sanitize_bim() → (bfile_prefix, detected_format)
    bfile_prefix   : Path prefix of the (possibly sanitized) PLINK panel.
    detected_format: 'rsid' or 'chrpos' — used by cli.py to drive VCF conversion.
"""

import os
import re
import gzip
import glob
import polars as pl

from .logger import get_logger, step_banner, log_pass, log_warn, log_info, abort

log = get_logger()
_RSID_RE = re.compile(r"^rs\d+$")


# ─────────────────────────────────────────────────────────────────────────────
# ID format detection helpers
# ─────────────────────────────────────────────────────────────────────────────

def _rsid_fraction(ids: list) -> float:
    """Return the fraction of IDs that match the rs\\d+ pattern."""
    if not ids:
        return 0.0
    return sum(1 for i in ids if _RSID_RE.match(str(i))) / len(ids)


def detect_bim_id_format(bim_path: str, n_sample: int = 2000) -> str:
    """
    Detect whether BIM variant IDs are rsID or CHROM_POS_REF_ALT.

    For a 7-column BIM (extra rsID column), samples BOTH candidate ID columns
    and uses the one with the higher rsID fraction.

    Returns: 'rsid' or 'chrpos'
    """
    ids_col1, ids_col2 = [], []
    with open(bim_path, "r") as fh:
        for i, line in enumerate(fh):
            if i >= n_sample:
                break
            parts = line.strip().split()
            if len(parts) >= 2:
                ids_col1.append(parts[1])
            if len(parts) >= 3:
                ids_col2.append(parts[2])

    frac1 = _rsid_fraction(ids_col1)
    frac2 = _rsid_fraction(ids_col2)
    n_cols = len(line.strip().split()) if ids_col1 else 6

    if n_cols == 7:
        # 7-col: col1 = CHRPOS-style ID, col2 = rsID
        log_info(log, f"7-col BIM: col[1] rsID fraction={frac1:.1%}, col[2] rsID fraction={frac2:.1%}")
        fmt = "rsid" if frac2 >= 0.5 else "chrpos"
        log_info(log, f"→ Using {'col[2]' if frac2 >= 0.5 else 'col[1]'} as canonical ID column")
    else:
        # Standard 6-col BIM
        fmt = "rsid" if frac1 >= 0.5 else "chrpos"

    log_info(log, f"BIM ID format detected: {fmt.upper()}  (rsID fraction = {max(frac1, frac2):.1%})")
    return fmt


def detect_ldscore_id_format(ld_dir: str, n_sample: int = 500) -> str:
    """
    Read the SNP column of the first available *.l2.ldscore.gz file
    to detect whether LD scores use rsID or CHROM_POS_REF_ALT variant IDs.

    Returns: 'rsid', 'chrpos', or 'unknown'
    """
    pattern = os.path.join(ld_dir.rstrip("/"), "*.l2.ldscore.gz")
    files   = sorted(glob.glob(pattern))
    if not files:
        log_warn(log, f"No *.l2.ldscore.gz files found in: {ld_dir}  (skipping LD format check)")
        return "unknown"

    ids = []
    try:
        with gzip.open(files[0], "rt") as fh:
            header = fh.readline().strip().split("\t")
            snp_col = header.index("SNP") if "SNP" in header else 1
            for i, line in enumerate(fh):
                if i >= n_sample:
                    break
                parts = line.strip().split("\t")
                if len(parts) > snp_col:
                    ids.append(parts[snp_col])
    except Exception as e:
        log_warn(log, f"Could not read LD score file {files[0]}: {e}")
        return "unknown"

    frac = _rsid_fraction(ids)
    fmt  = "rsid" if frac >= 0.5 else "chrpos"
    log_info(log, f"LD score file: {os.path.basename(files[0])}")
    log_info(log, f"LD score ID format detected: {fmt.upper()}  (rsID fraction = {frac:.1%})")
    return fmt


# ─────────────────────────────────────────────────────────────────────────────
# Main entry point
# ─────────────────────────────────────────────────────────────────────────────

def sanitize_bim(
    bfile: str,
    out_dir: str,
    ld_dir: str = None
) -> tuple:
    """
    Validate BIM ID format against LD score IDs, sanitize 7-column BIM if needed,
    and return the final panel prefix and detected format.

    Rules:
      - BIM IDs and LD score IDs MUST be the same format. If not → ABORT.
      - IDs are NEVER modified. Only the column layout is fixed (7-col → 6-col).
      - The detected format is returned for use in VCF variant ID construction.

    Args:
        bfile   : PLINK reference panel prefix.
        out_dir : Directory for sanitized output files.
        ld_dir  : LD score files directory (for format cross-check).

    Returns:
        (bfile_prefix: str, detected_format: str)
    """
    step_banner(log, "Validating Reference Panel & LD Score ID Formats", step=1, total=5)

    bim_path = f"{bfile}.bim"

    # ── File existence checks ────────────────────────────────────────────────
    for ext, label in [(".bim", "BIM"), (".bed", "BED"), (".fam", "FAM")]:
        fpath = f"{bfile}{ext}"
        if not os.path.exists(fpath):
            abort(log,
                  f"{label} file not found: {fpath}\n"
                  f"  Ensure --bfile points to a valid PLINK binary panel prefix.")

    log_info(log, f"PLINK panel : {bfile}")

    # ── Detect BIM column count & ID format ──────────────────────────────────
    with open(bim_path) as fh:
        first_line = fh.readline().strip().split()
    n_cols = len(first_line)
    log_info(log, f"BIM columns : {n_cols}")

    bim_fmt = detect_bim_id_format(bim_path)

    # ── Detect LD score ID format ────────────────────────────────────────────
    ld_fmt = "unknown"
    if ld_dir:
        ld_fmt = detect_ldscore_id_format(ld_dir)

    # ── Cross-validate BIM vs LD score format ────────────────────────────────
    if ld_fmt != "unknown" and bim_fmt != ld_fmt:
        abort(log,
              f"ID format mismatch between BIM and LD score files!\n"
              f"  BIM format      : {bim_fmt.upper()}\n"
              f"  LD score format : {ld_fmt.upper()}\n"
              f"  BIM file        : {bim_path}\n"
              f"  LD score dir    : {ld_dir}\n\n"
              f"  Both PLINK BIM and LD score files must use the same variant ID\n"
              f"  format (either both rsID or both CHROM_POS_REF_ALT). GCTA's\n"
              f"  LDSC regression requires IDs to match across these files.\n\n"
              f"  Fix: Ensure you are using matched BIM + LD score resources.")

    if ld_fmt == "unknown":
        log_warn(log, "LD score format check skipped — no LD score directory provided or no files found.")
    else:
        log_pass(log, f"BIM and LD score ID formats match: {bim_fmt.upper()}")

    # ── If already 6-column, return as-is ────────────────────────────────────
    if n_cols == 6:
        log_pass(log, "BIM file is standard 6-column format — no layout fix needed.")
        return bfile, bim_fmt

    # ── Sanitize 7-column → 6-column BIM (preserve original IDs) ────────────
    out_prefix = os.path.join(out_dir, f"{os.path.basename(bfile)}_6col")
    out_bim    = f"{out_prefix}.bim"

    if os.path.exists(out_bim):
        log_pass(log, f"Sanitized 6-col BIM already exists — reusing: {out_bim}")
    else:
        log_info(log, "Converting 7-column BIM → 6-column BIM (IDs are NOT changed) ...")
        df = pl.read_csv(bim_path, separator="\t", has_header=False)
        cols = df.columns

        # For 7-col BIM: col[1]=CHRPOS, col[2]=rsID
        # Pick the column that matches the detected format
        if bim_fmt == "rsid":
            # Use col[2] (rsID column); fall back to col[1] if "."
            id_col = (
                pl.when((pl.col(cols[2]) != ".") & pl.col(cols[2]).is_not_null())
                .then(pl.col(cols[2]))
                .otherwise(pl.col(cols[1]))
            )
        else:
            # Use col[1] (CHRPOS column)
            id_col = pl.col(cols[1])

        df_6col = df.select([
            pl.col(cols[0]),       # CHR
            id_col.alias("SNP"),   # Variant ID (original, unchanged)
            pl.col(cols[3]),       # cM
            pl.col(cols[4]),       # BP
            pl.col(cols[5]),       # A1
            pl.col(cols[6])        # A2
        ])
        df_6col.write_csv(out_bim, separator="\t", include_header=False)
        log_pass(log, f"Sanitized BIM written: {out_bim}  ({len(df_6col):,} variants)")

    # ── Symlink .bed and .fam ─────────────────────────────────────────────────
    for ext in [".bed", ".fam"]:
        src = f"{bfile}{ext}"
        dst = f"{out_prefix}{ext}"
        if not os.path.exists(dst):
            os.symlink(src, dst)
            log_info(log, f"Symlinked {ext}: {dst}")

    log_pass(log, f"Reference panel ready: {out_prefix}")
    return out_prefix, bim_fmt
