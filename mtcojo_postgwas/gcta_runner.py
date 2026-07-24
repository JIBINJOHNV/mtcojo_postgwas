"""
gcta_runner.py - GCTA mtCOJO Process Launcher for mtcojo_postgwas

Constructs command arguments, logs them, and executes GCTA gcta64
for multi-trait conditional and joint analysis (mtCOJO).
"""

import subprocess
from typing import Optional

from .logger import get_logger, step_banner, log_pass, log_info, abort, log_cmd_script

log = get_logger()



def run_gcta_mtcojo(
    gcta64_bin: str,
    bfile: str,
    list_file: str,
    out_prefix: str,
    ref_ld_chr: Optional[str] = None,
    w_ld_chr: Optional[str] = None,
    mbfile: Optional[str] = None,
    gwas_thresh: Optional[float] = 1e-5,
    clump_r2: Optional[float] = None,
    heidi_thresh: Optional[float] = None,
    gsmr_snp_min: Optional[int] = None,
    diff_freq: Optional[float] = None,
    mtcojo_bxy: Optional[str] = None,
    gsmr_ld_fdr: Optional[float] = None
) -> None:
    """
    Execute GCTA mtCOJO multi-trait conditional analysis.

    Args:
        gcta64_bin  : Path to gcta64 binary.
        bfile       : PLINK reference panel prefix.
        list_file   : .mtcojo.list file mapping traits to .ma files.
        out_prefix  : Prefix for GCTA output files.
        ref_ld_chr  : LD score files directory.
        w_ld_chr    : LD score weight directory (defaults to ref_ld_chr).
        mbfile      : File listing multiple PLINK panels.
        gwas_thresh : Index SNP p-value threshold (default 1e-5).
        clump_r2    : LD r² threshold for clumping.
        heidi_thresh: p-value threshold for HEIDI outlier test.
        gsmr_snp_min: Minimum index SNPs for GSMR.
        diff_freq   : Max allele frequency difference vs reference.
        mtcojo_bxy  : Pre-computed bxy causal effect file.
        gsmr_ld_fdr : FDR threshold for LD filtering.
    """
    import os
    step_banner(log, "GCTA mtCOJO Analysis", step=5, total=6)

    # ── Validate gcta64 binary ────────────────────────────────────────────────
    if not os.path.exists(gcta64_bin):
        abort(log,
              f"gcta64 binary not found: {gcta64_bin}\n"
              f"  Pass --gcta64 /path/to/gcta64 or activate the conda environment.")

    # ── Validate inputs ───────────────────────────────────────────────────────
    for label, path in [("mtcojo-list", list_file), ("PLINK BIM", f"{bfile}.bim" if bfile else None)]:
        if path and not os.path.exists(path):
            abort(log, f"Required input file not found [{label}]: {path}")

    if ref_ld_chr and not os.path.isdir(ref_ld_chr):
        abort(log, f"LD score directory not found: {ref_ld_chr}")

    # ── Build command ─────────────────────────────────────────────────────────
    ref_args = ["--mbfile", mbfile] if mbfile else ["--bfile", bfile]
    ld = ref_ld_chr
    w_ld = w_ld_chr or ld

    cmd = [gcta64_bin] + ref_args + ["--mtcojo-file", list_file, "--out", out_prefix]

    if ld:
        cmd += ["--ref-ld-chr", ld]
    if w_ld:
        cmd += ["--w-ld-chr", w_ld]

    optional_flags = [
        ("--gwas-thresh", gwas_thresh),
        ("--clump-r2", clump_r2),
        ("--heidi-thresh", heidi_thresh),
        ("--gsmr-snp-min", gsmr_snp_min),
        ("--diff-freq", diff_freq),
        ("--mtcojo-bxy", mtcojo_bxy),
        ("--gsmr-ld-fdr", gsmr_ld_fdr),
    ]
    for flag, val in optional_flags:
        if val is not None:
            cmd += [flag, str(val)]

    # ── Log full command ──────────────────────────────────────────────────────
    cmd_str = " ".join(cmd)
    log_cmd_script(out_prefix, "GCTA mtCOJO Analysis", cmd_str)
    log_info(log, "GCTA command:")
    log_info(log, "  " + " \\\n      ".join(cmd))

    # ── Execute ───────────────────────────────────────────────────────────────
    log_info(log, "Executing GCTA mtCOJO ... (this may take several minutes)")
    try:
        subprocess.run(cmd, check=True)
    except subprocess.CalledProcessError as e:
        abort(log,
              f"GCTA mtCOJO failed with exit code {e.returncode}.\n"
              f"  Check the GCTA log file: {out_prefix}.log\n"
              f"  Common causes: mismatched variant IDs between .ma files and BIM, "
              f"duplicated SNP IDs in .ma files, insufficient significant SNPs, "
              f"or incompatible LD score files.")

    # ── Validate output ───────────────────────────────────────────────────────
    cma_file = f"{out_prefix}.mtcojo.cma"
    if not os.path.exists(cma_file):
        abort(log,
              f"Expected GCTA output not found: {cma_file}\n"
              f"  Check {out_prefix}.log for GCTA error messages.")

    log_pass(log, f"GCTA mtCOJO completed. Output prefix: {out_prefix}")
    log_info(log, "Output files:")
    for suffix in [".mtcojo.cma", ".mtcojo.list", ".badsnps", ".pleio_snps", ".log"]:
        p = f"{out_prefix}{suffix}"
        if os.path.exists(p):
            size = os.path.getsize(p)
            log_info(log, f"  {suffix:<22} {size:>10,} bytes  → {p}")
