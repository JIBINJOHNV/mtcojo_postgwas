"""
vcf_converter.py - Single-Pass VCF Extraction & GCTA .ma Exporter

Extracts all summary statistic fields from a GWAS VCF file in ONE bcftools query pass
and writes a GCTA-format .ma file.

Variant ID format is automatically determined by the caller (cli.py) based on the
detected format of the PLINK BIM and LD score files — ensuring consistency across
all inputs without requiring the user to specify it manually.

VCF FORMAT fields extracted:
  %ID       → rsID (if annotated)
  %CHROM    → chromosome
  %POS      → base-pair position
  %REF / %ALT → reference / effect allele
  [%AF]     → allele frequency
  [%ES]     → effect size (beta)
  [%SE]     → standard error
  [%LP]     → −log₁₀(p-value)  →  p = 10^(−LP)
  [%NEF]    → effective sample size (N_eff)
  [%SI]     → imputation INFO score
"""

import subprocess
import io
import os
import polars as pl

from mtcojo_postgwas.core.logger import get_logger, step_banner, log_pass, log_warn, log_info, abort, log_cmd_script

log = get_logger()

_MIN_OVERLAP_FRAC = 0.01   # abort if <1% of VCF variants match BIM


def convert_vcf_single_pass(
    vcf_path: str,
    ma_path: str,
    bcftools_bin: str = "bcftools",
    id_format: str = "rsid",   # 'rsid' or 'chrpos' — auto-detected from BIM/LD scores
    bim_ids: set = None        # optional set for overlap validation
) -> pl.DataFrame:
    """
    Extract VCF fields in a single bcftools query pass.

    Args:
        vcf_path    : Path to GWAS summary statistics VCF (.vcf or .vcf.gz).
        ma_path     : Output path for the GCTA .ma file.
        bcftools_bin: Path to bcftools binary.
        id_format   : 'rsid' or 'chrpos' — set automatically from BIM/LD score detection.
        bim_ids     : Set of BIM variant IDs for overlap validation.

    Returns:
        pl.DataFrame with columns [SNP, CHR, POS, SI] for PostGWAS in-memory join.
    """
    step_banner(log, f"VCF Conversion  →  {os.path.basename(ma_path)}")
    log_info(log, f"Input VCF    : {vcf_path}")
    log_info(log, f"Output .ma   : {ma_path}")
    log_info(log, f"ID format    : {id_format.upper()}  (auto-detected from BIM/LD scores)")

    # ── bcftools sanity check ─────────────────────────────────────────────────
    try:
        ver = subprocess.run(
            [bcftools_bin, "--version"], stdout=subprocess.PIPE, stderr=subprocess.PIPE
        ).stdout.decode().split("\n")[0]
        log_info(log, f"bcftools     : {ver}")
    except FileNotFoundError:
        abort(log, f"bcftools not found: {bcftools_bin}\n  Pass --bcftools /path/to/bcftools")

    if not os.path.exists(vcf_path):
        abort(log, f"VCF file not found: {vcf_path}")

    # ── Single-pass bcftools query ────────────────────────────────────────────
    log_info(log, "Running bcftools query (single pass, all fields) ...")
    cmd = [
        bcftools_bin, "query",
        "-f", "%ID\t%CHROM\t%POS\t%REF\t%ALT\t[%AF]\t[%ES]\t[%SE]\t[%LP]\t[%NEF]\t[%SI]\n",
        vcf_path
    ]
    out_prefix_base = os.path.splitext(ma_path)[0]
    log_cmd_script(out_prefix_base, f"VCF Conversion (bcftools query) [{os.path.basename(ma_path)}]", f"{' '.join(cmd)} > {ma_path}")
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                             text=True, check=True)
    except subprocess.CalledProcessError as e:
        abort(log,
              f"bcftools query failed for: {vcf_path}\n"
              f"  stderr: {e.stderr.strip()}")

    if not res.stdout.strip():
        abort(log, f"bcftools returned empty output for: {vcf_path}\n"
                   f"  Ensure the VCF contains FORMAT fields: ES, SE, LP, NEF (or SS), SI.")

    # ── Parse into Polars ─────────────────────────────────────────────────────
    df = pl.read_csv(
        io.StringIO(res.stdout),
        separator="\t", has_header=False,
        new_columns=["ID", "CHR", "POS", "REF", "ALT", "AF", "ES", "SE", "LP", "NEF", "SI"],
        infer_schema_length=10000
    )
    log_info(log, f"Raw variants extracted   : {len(df):,}")

    # ── Build variant IDs matching BIM/LD format ──────────────────────────────
    use_rsid = (id_format == "rsid")
    log_info(log, f"Building variant IDs in {'rsID' if use_rsid else 'CHROM_POS_REF_ALT'} format ...")

    chrpos_id = pl.concat_str(
        [pl.col("CHR"), pl.col("POS").cast(str), pl.col("REF"), pl.col("ALT")],
        separator="_"
    )

    if use_rsid:
        snp_expr = (
            pl.when((pl.col("ID") != ".") & pl.col("ID").is_not_null())
            .then(pl.col("ID"))
            .otherwise(chrpos_id)   # fallback if rsID missing
        )
    else:
        snp_expr = chrpos_id

    df = df.with_columns([
        snp_expr.alias("SNP"),
        pl.col("ALT").alias("A1"),
        pl.col("REF").alias("A2"),
        pl.col("AF").cast(pl.Float64, strict=False).alias("freq"),
        pl.col("ES").cast(pl.Float64, strict=False).alias("b"),
        pl.col("SE").cast(pl.Float64, strict=False).alias("se"),
        (10.0 ** (-pl.col("LP").cast(pl.Float64, strict=False))).clip(1e-300, 1.0).alias("p"),
        pl.col("NEF").cast(pl.Float64, strict=False).alias("N"),
        pl.col("SI").cast(pl.Float64, strict=False)
    ])

    # ── Drop incomplete rows ──────────────────────────────────────────────────
    n_before = len(df)
    df = df.drop_nulls(subset=["SNP", "A1", "A2", "freq", "b", "se", "p", "N"])
    n_dropped = n_before - len(df)
    if n_dropped:
        log_warn(log, f"{n_dropped:,} variants removed — missing critical fields (freq/b/se/p/N)")
    log_info(log, f"Variants after QC        : {len(df):,}")

    # GCTA mtCOJO is designed for biallelic SNP summary statistics. Remove
    # indels and non-standard alleles before de-duplicating rsIDs.
    allele_ok = (
        pl.col("A1").str.to_uppercase().is_in(["A", "C", "G", "T"]) &
        pl.col("A2").str.to_uppercase().is_in(["A", "C", "G", "T"])
    )
    indel_df = df.filter(~allele_ok)
    if len(indel_df) > 0:
        indel_path = f"{os.path.splitext(ma_path)[0]}.indels_removed.tsv"
        indel_df.select(["SNP", "CHR", "POS", "A1", "A2", "freq", "b", "se", "p", "N", "SI"]).write_csv(
            indel_path,
            separator="\t",
        )
        df = df.filter(allele_ok)
        log_warn(log, f"{len(indel_df):,} indel/non-ACGT allele rows removed before GCTA export. Details: {indel_path}")
        log_info(log, f"Variants after SNP-only filter: {len(df):,}")

    # GCTA requires each SNP ID to appear only once in every .ma file. Public
    # GWAS VCFs can contain repeated rsIDs for multi-allelic or duplicated records.
    dup_counts = df.group_by("SNP").len().filter(pl.col("len") > 1)
    if len(dup_counts) > 0:
        n_dup_rows = int(dup_counts["len"].sum() - len(dup_counts))
        dup_path = f"{os.path.splitext(ma_path)[0]}.duplicate_snps.tsv"
        dup_details = df.join(dup_counts.select("SNP"), on="SNP", how="inner")
        dup_details.select(["SNP", "CHR", "POS", "A1", "A2", "freq", "b", "se", "p", "N", "SI"]).write_csv(
            dup_path,
            separator="\t",
        )
        df = (
            df
            .sort(["SNP", "p", "N", "SI"], descending=[False, False, True, True], nulls_last=True)
            .unique(subset=["SNP"], keep="first", maintain_order=True)
        )
        log_warn(
            log,
            f"{n_dup_rows:,} duplicate SNP rows removed across {len(dup_counts):,} SNP IDs "
            f"(kept lowest p-value, then highest N/SI). Details: {dup_path}"
        )
        log_info(log, f"Variants after SNP de-dup  : {len(df):,}")

    # ── BIM overlap validation ────────────────────────────────────────────────
    if bim_ids:
        vcf_ids    = set(df["SNP"].to_list())
        overlap    = vcf_ids & bim_ids
        frac       = len(overlap) / len(vcf_ids) if vcf_ids else 0.0
        log_info(log,
                 f"Overlap check  →  VCF:{len(vcf_ids):,}  BIM:{len(bim_ids):,}  "
                 f"Shared:{len(overlap):,}  ({frac:.1%})")
        if frac < _MIN_OVERLAP_FRAC:
            abort(log,
                  f"Critically low variant ID overlap: {frac:.2%}  ({len(overlap):,} SNPs)\n"
                  f"  VCF format   : {id_format.upper()}\n"
                  f"  VCF IDs (5)  : {list(vcf_ids)[:5]}\n"
                  f"  BIM IDs (5)  : {list(bim_ids)[:5]}\n"
                  f"  Both VCF and BIM must carry the same variant ID type.\n"
                  f"  The pipeline auto-detects the BIM format — check that your VCF\n"
                  f"  is annotated with {'rsIDs' if use_rsid else 'CHROM_POS_REF_ALT IDs'}.")
        elif frac < 0.30:
            log_warn(log, f"Low overlap ({frac:.1%}) — verify genome build and ID annotation.")
        else:
            log_pass(log, f"Variant ID overlap: {frac:.1%}  ({len(overlap):,} shared SNPs)  ✔")

    # ── Write GCTA .ma ────────────────────────────────────────────────────────
    df.select(["SNP", "A1", "A2", "freq", "b", "se", "p", "N"]).write_csv(ma_path, separator="\t")
    log_pass(log, f"GCTA .ma written: {ma_path}  ({len(df):,} variants)")

    coords_df = df.select(["SNP", "CHR", "POS", "SI"])
    coords_path = f"{os.path.splitext(ma_path)[0]}.coords.tsv"
    coords_df.write_csv(coords_path, separator="\t")
    log_info(log, f"Coordinate cache written: {coords_path}")

    return coords_df
