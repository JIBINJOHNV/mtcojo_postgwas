"""
cli.py - Command Line Interface for mtcojo_postgwas

Pipeline execution order:
  [1] Validate inputs (manifest, tool paths)
  [2] Detect variant ID format from BIM + LD score files  →  auto-select rsID or CHRPOS
      Abort if BIM and LD score formats differ
  [3] Sanitize 7-col BIM → 6-col BIM (IDs unchanged)
  [4] Load BIM variant IDs into memory for overlap validation
  [5] Convert each GWAS VCF → GCTA .ma  (single bcftools pass, ID format from step 2)
  [6] Run GCTA mtCOJO
  [7] Optional: PostGWAS Docker harmonisation  (--run-postgwas)
        Requires: --defaults, --resource-folder
        Validates these UPFRONT before any analysis starts.
  [8] Optional: LDSC genetic correlation        (--run-ldsc)
        Independent of PostGWAS — uses .ma files + .mtcojo.cma directly.
        Pass --run-postgwas separately if harmonisation output is also needed.
"""

from __future__ import annotations

import os
import sys
import time
import argparse
import glob

from .core.logger import (
    setup_logger, get_logger, print_logo,
    step_banner, log_pass, log_warn, log_info, abort, summary_table
)

DEF_BCFTOOLS           = "bcftools"
DEF_GCTA64             = "gcta64"
DEF_DOCKER_IMAGE       = "jibinjv/postgwas:1.4"
DEF_RESOURCE_FOLDER    = None
DEF_HARMONISATION_YAML = None


def _pl():
    import polars as pl

    return pl


def probability(value: str) -> float:
    """Parse a probability in the closed interval [0, 1]."""
    try:
        number = float(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Expected a numeric value, received: {value!r}") from exc
    if not 0.0 <= number <= 1.0:
        raise argparse.ArgumentTypeError(f"Value must be between 0 and 1, received: {number}")
    return number


def positive_integer(value: str) -> int:
    """Parse an integer greater than zero."""
    try:
        number = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"Expected an integer, received: {value!r}") from exc
    if number < 1:
        raise argparse.ArgumentTypeError(f"Value must be at least 1, received: {number}")
    return number


class HelpFormatter(argparse.RawTextHelpFormatter):
    """Readable CLI help with defaults and extra space between sections."""

    def start_section(self, heading):
        if self._current_section.heading is not None:
            self._add_item(lambda: "\n", [])
        super().start_section(heading)

    def _get_help_string(self, action):
        help_text = action.help or ""
        if "%(default)" in help_text or action.default is argparse.SUPPRESS:
            return help_text
        if not action.option_strings or action.required:
            return help_text
        default = "not set" if action.default is None else action.default
        if default == "":
            default = "auto-detect"
        return f"{help_text} (default: {default})"


def _load_bim_ids(bim_path: str) -> set:
    """Load all variant IDs from a 6-column BIM file into a Python set."""
    ids = set()
    with open(bim_path) as fh:
        for line in fh:
            p = line.strip().split()
            if len(p) >= 2:
                ids.add(p[1])
    return ids


def _nonempty(path: str) -> bool:
    return bool(path) and os.path.exists(path) and os.path.getsize(path) > 0


def _load_coords_cache(ma_path: str) -> pl.DataFrame:
    pl = _pl()
    coords_path = f"{os.path.splitext(ma_path)[0]}.coords.tsv"
    if _nonempty(coords_path):
        return pl.read_csv(coords_path, separator="\t")
    return pl.DataFrame(schema={"SNP": pl.String, "CHR": pl.Int64, "POS": pl.Int64, "SI": pl.Float64})


def _load_target_coords_for_resume(ma_path: str, merged_tsv: str) -> pl.DataFrame:
    pl = _pl()
    coords_df = _load_coords_cache(ma_path)
    if len(coords_df) > 0:
        return coords_df
    if _nonempty(merged_tsv):
        return (
            pl.scan_csv(merged_tsv, separator="\t")
            .select([
                pl.col("SNP").cast(pl.String),
                pl.col("CHR").cast(pl.Int64, strict=False),
                pl.col("BP").cast(pl.Int64, strict=False).alias("POS"),
                pl.lit(None, dtype=pl.Float64).alias("SI"),
            ])
            .collect()
        )
    return coords_df


def _postgwas_outputs_exist(sdd_dir: str, out_name: str) -> bool:
    out_folder = os.path.join(os.path.abspath(sdd_dir), "03_harmonised_output")
    patterns = [
        os.path.join(out_folder, f"{out_name}*"),
        os.path.join(out_folder, "**", f"{out_name}*"),
    ]
    for pattern in patterns:
        for path in glob.glob(pattern, recursive=True):
            if os.path.isfile(path) and _nonempty(path):
                return True
    return False


def main():
    p = argparse.ArgumentParser(
        prog="mtcojo-postgwas",
        description=(
            "GCTA mtCOJO pipeline for GWAS VCF summary statistics.\n"
            "Converts GWAS VCF files to GCTA .ma, runs mtCOJO, optionally runs\n"
            "PostGWAS harmonisation and LDSC h2/rg, then builds a standalone HTML report.\n"
            "Reruns are resumable by default: completed non-empty outputs are reused\n"
            "unless --force is set."
        ),
        formatter_class=HelpFormatter
    )

    req = p.add_argument_group("Required")
    req.add_argument("-m", "--manifest", required=True,
                     metavar="CSV",
                     help=("Trait manifest CSV with at least two rows:\n"
                           "  row 1 = target trait; rows 2+ = covariates.\n"
                           "Required columns: sample_id,file_path. Optional columns for LDSC liability scale: "
                           "sample_prevalence,population_prevalence."))

    out = p.add_argument_group("Output and Resume")
    out.add_argument("-o", "--out",     default="scz_bip_mtcojo", metavar="PREFIX",
                     help="Output prefix used for GCTA, LDSC, PostGWAS, and final report files.")
    out.add_argument("-d", "--out-dir", metavar="DIR",
                     help="Output directory. Step subdirectories are created inside this directory.")
    out.add_argument("--force",         action="store_true",
                     help="Re-run all stages even when expected non-empty outputs already exist.")
    out.add_argument("--force-report",  action="store_true",
                     help="Regenerate the final standalone HTML report even when it already exists.")

    ref = p.add_argument_group("Reference and Variant-ID Inputs")
    genotype_ref = ref.add_mutually_exclusive_group(required=True)
    genotype_ref.add_argument("-b", "--bfile",   metavar="PREFIX",
                              help=("PLINK reference panel prefix for GCTA --bfile; expects PREFIX.bed/.bim/.fam.\n"
                                    "Used for allele-frequency checks, LD clumping, and LD among GSMR instruments."))
    genotype_ref.add_argument("--mbfile",        metavar="FILE",
                              help=("GCTA --mbfile input listing multiple PLINK binary reference panel prefixes,\n"
                                    "typically one per chromosome. Used instead of --bfile."))
    ref.add_argument("-l", "--ref-ld-chr", "--ld-dir", dest="ref_ld_chr", metavar="DIR_OR_PREFIX",
                     help="LD score reference path passed to GCTA --ref-ld-chr and used for ID-format checks.")
    ref.add_argument("-w", "--w-ld-chr", metavar="DIR_OR_PREFIX",
                     help="LD score weights path passed to GCTA --w-ld-chr. If omitted, uses --ref-ld-chr.")
    ref.add_argument("--no-bim-fix", action="store_true",
                     help="Skip automatic 7-column BIM sanitization. ID-format detection still runs.")

    gcta = p.add_argument_group("GCTA mtCOJO Options")
    gcta.add_argument("--gwas-thresh",  type=probability, default=1e-5, metavar="P",
                      help=("GCTA --gwas-thresh value: p-value threshold used to select candidate SNP instruments\n"
                            "for each conditioning trait before LD clumping."))
    gcta.add_argument("--clump-r2",     type=probability, metavar="R2",
                      help=("Optional GCTA --clump-r2 value: LD r-squared threshold for clumping candidate instruments.\n"
                            "If omitted, this flag is not passed to GCTA."))
    gcta.add_argument("--heidi-thresh", type=probability, metavar="P",
                      help=("Optional GCTA --heidi-thresh value: p-value threshold for HEIDI outlier filtering\n"
                            "of SNP instruments. Set to 0 to disable HEIDI in GCTA versions that support that behavior.\n"
                            "If omitted, GCTA uses its internal default."))
    gcta.add_argument("--gsmr-snp-min", type=positive_integer, metavar="N",
                      help=("Optional GCTA --gsmr-snp-min value: minimum number of valid, quasi-independent\n"
                            "SNP instruments required after clumping/HEIDI filtering."))
    gcta.add_argument("--diff-freq",    type=probability, metavar="FRAC",
                      help=("Optional GCTA --diff-freq value: maximum absolute effect-allele-frequency difference\n"
                            "between GWAS summary statistics and the genotype reference."))
    gcta.add_argument("--mtcojo-bxy",   metavar="FILE",
                      help="Optional precomputed GCTA --mtcojo-bxy causal-effect file.")
    gcta.add_argument("--gsmr-ld-fdr",  type=probability, metavar="FDR",
                      help=("Optional GCTA --gsmr-ld-fdr value: FDR threshold for residual LD/correlation among\n"
                            "SNP instruments. If omitted, this flag is not passed to GCTA."))

    tools = p.add_argument_group("Tool Paths")
    tools.add_argument("--bcftools", default=DEF_BCFTOOLS, metavar="PATH",
                       help="Path to bcftools executable used for a single-pass VCF FORMAT query.")
    tools.add_argument("--gcta64",   default=DEF_GCTA64, metavar="PATH",
                       help="Path to GCTA executable used for --mtcojo.")

    pg = p.add_argument_group("PostGWAS Harmonisation  (enable with --run-postgwas)")
    pg.add_argument("--run-postgwas",    action="store_true",
	                                         help="Run optional PostGWAS Docker harmonisation on the mtCOJO .cma output.")
    pg.add_argument("--sdd-dir",         metavar="DIR",
                                         help="Host directory mounted as /mnt/disks/sdd/ in the PostGWAS Docker container.")
    pg.add_argument("--defaults", "--default-config", dest="defaults",
                                         default=DEF_HARMONISATION_YAML,
                                         metavar="YAML",
	                                         help="PostGWAS harmonisation.yaml defaults/config file.")
    pg.add_argument("--liftover",        default="No", choices=["Yes","No"],
	                                         help="Whether PostGWAS should perform GRCh37-to-GRCh38 liftover.")
    pg.add_argument("--resource-folder", default=DEF_RESOURCE_FOLDER,
                                         metavar="DIR",
	                                         help="PostGWAS/gwas2vcf reference resources directory.")
    pg.add_argument("--docker-image",    default=DEF_DOCKER_IMAGE, metavar="IMAGE",
                    help="Docker image used for PostGWAS harmonisation.")
    pg.add_argument("--docker-platform", default="linux/amd64", metavar="PLATFORM",
                    help="Docker platform passed to docker run.")
    pg.add_argument("--nthreads",        type=positive_integer, default=23, metavar="N",
                    help="Thread count passed to PostGWAS inside Docker.")
    pg.add_argument("--max-mem",         default="50G", metavar="SIZE",
                    help="Maximum memory string passed to PostGWAS inside Docker.")

    ldsc = p.add_argument_group(
        "LDSC Genetic Correlation  (enable with --run-ldsc; independent of --run-postgwas)")
    ldsc.add_argument("--run-ldsc",        action="store_true",
                                           help=("Run LDSC rg analysis directly on .ma files + .mtcojo.cma output.\n"
                                                 "Does not require --run-postgwas; pass both flags when both outputs are needed."))
    ldsc.add_argument("--ldsc-dir",        dest="ldsc_dir",       default="", metavar="DIR",
	                                           help="Path to CBIIT/ldsc clone. Empty means auto-detect from install.sh paths.")
    ldsc.add_argument("--ldsc-ld-dir",     dest="ldsc_ld_dir",
                                           metavar="DIR_OR_PREFIX",
	                                           help="LD score reference path for LDSC --ref-ld-chr/--w-ld-chr. If omitted, uses --ref-ld-chr.")
    ldsc.add_argument("--ldsc-snp-list",   dest="ldsc_snp_list",
                                           metavar="FILE",
	                                           help="HapMap3 SNP include list for LDSC munge_sumstats.py, usually w_hm3.snplist.")
    ldsc.add_argument("--ldsc-n-parallel", dest="ldsc_n_parallel", type=positive_integer, default=4,
                                           metavar="N",
	                                           help="Number of parallel LDSC munge/h2/rg worker processes.")
    ldsc.add_argument("--ldsc-batch-size", dest="ldsc_batch_size",  type=positive_integer, default=10,
                                           metavar="N",
	                                           help="Number of target traits included per ldsc.py --rg batch.")
    ldsc.add_argument("--ldsc-center-z",   dest="ldsc_center_z",    action="store_true",
	                                           help=("Center Z-scores to median 0 before LDSC munge.\n"
                                                     "Useful for regional/test data that fail LDSC mean chi-square checks."))

    args = p.parse_args()
    pl = _pl()

    # ── Validate PostGWAS inputs UPFRONT before any analysis starts ────────────
    # This prevents a late failure after hours of GCTA processing.
    if args.run_postgwas:
        pg_errors = []
        if not args.defaults or not os.path.exists(args.defaults):
            pg_errors.append(
                f"  --defaults / harmonisation.yaml not found: {args.defaults}\n"
                f"    Pass: --defaults /path/to/harmonisation.yaml"
            )
        if not args.resource_folder or not os.path.isdir(args.resource_folder):
            pg_errors.append(
                f"  --resource-folder (gwas2vcf resources) not found: {args.resource_folder}\n"
                f"    Pass: --resource-folder /path/to/postgwas/gwas2vcf/"
            )
        if pg_errors:
            print()
            print("ERROR: --run-postgwas is enabled but required inputs are missing:")
            for e in pg_errors:
                print(e)
            print()
            print("Required for PostGWAS harmonisation:")
            print("  --defaults         Path to harmonisation.yaml config file")
            print("  --resource-folder  Path to gwas2vcf reference resources directory")
            print("  --sdd-dir          Host directory mounted into Docker (auto-created if omitted)")
            print()
            sys.exit(1)

    # ── Determine step total based on active stages ────────────────────────────
    n_steps = 6 + int(args.run_postgwas) + int(args.run_ldsc)

    # ── Setup step-numbered output subdirectories ──────────────────────────────
    out_dir  = os.path.abspath(args.out_dir or ".")
    out_name = args.out or "gcta_mtcojo"

    step0_dir = os.path.join(out_dir, "00_manifest_and_logs")
    step1_dir = os.path.join(out_dir, "01_gcta_ma_conversion")
    step2_dir = os.path.join(out_dir, "02_gcta_mtcojo_results")
    step3_dir = os.path.join(out_dir, "03_postgwas_harmonisation")
    step4_dir = os.path.join(out_dir, "04_ldsc_analysis")
    step5_dir = os.path.join(out_dir, "05_plots_and_tables")
    merged_tsv = os.path.join(step5_dir, "merged_gwas_summary.tsv")

    for d in (step0_dir, step1_dir, step2_dir, step3_dir, step4_dir, step5_dir):
        os.makedirs(d, exist_ok=True)

    out_prefix = os.path.join(step2_dir, out_name)
    log_file   = os.path.join(step0_dir, f"{out_name}.pipeline.log")

    log = setup_logger(log_file=log_file)

    # ── Print logo ────────────────────────────────────────────────────────────
    print_logo()

    log.info(f"  Log file       : {log_file}")
    log.info(f"  Manifest       : {args.manifest}")
    log.info(f"  Output prefix  : {out_prefix}")
    log.info(f"  PostGWAS step  : {'ENABLED' if args.run_postgwas else 'DISABLED'}")
    log.info(f"  LDSC step      : {'ENABLED' if args.run_ldsc else 'DISABLED'}")
    log.info(f"  bcftools       : {args.bcftools}")
    log.info(f"  gcta64         : {args.gcta64}")

    t0 = time.time()
    reran_analysis = False

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 1 — Validate manifest
    # ─────────────────────────────────────────────────────────────────────────
    step_banner(log, "Manifest Validation", step=1, total=n_steps)

    if not os.path.exists(args.manifest):
        abort(log, f"Manifest file not found: {args.manifest}")

    df_man = pl.read_csv(args.manifest)
    missing = {"sample_id", "file_path"} - set(df_man.columns)
    if missing:
        abort(log, f"Manifest missing required columns: {missing}")
    df_man = df_man.filter(
        pl.col("sample_id").is_not_null() & pl.col("file_path").is_not_null()
    )
    if len(df_man) < 2:
        abort(log, "Manifest must contain at least 2 rows (target trait + ≥1 covariate).")

    for r in df_man.iter_rows(named=True):
        if not os.path.exists(r["file_path"]):
            abort(log, f"VCF not found for '{r['sample_id']}': {r['file_path']}")

    traits = df_man["sample_id"].to_list()
    log_pass(log, f"Manifest OK — {len(df_man)} traits: {traits}")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 2 — Detect ID format from BIM + LD scores, validate match,
    #           sanitize 7-col BIM if needed
    # ─────────────────────────────────────────────────────────────────────────
    bfile     = args.bfile
    id_format = "rsid"   # default; will be overridden by BIM/LD detection

    if bfile:
        if args.no_bim_fix:
            log_warn(log, "--no-bim-fix: BIM sanitization skipped. Format detection still runs.")
            from .io.bim_sanitizer import detect_bim_id_format, detect_ldscore_id_format
            step_banner(log, "ID Format Detection (BIM + LD Scores)", step=2, total=6)
            bim_fmt = detect_bim_id_format(f"{bfile}.bim")
            ld_fmt  = detect_ldscore_id_format(args.ref_ld_chr) if args.ref_ld_chr else "unknown"
            if ld_fmt != "unknown" and bim_fmt != ld_fmt:
                abort(log,
                      f"BIM and LD score ID format mismatch!\n"
                      f"  BIM: {bim_fmt.upper()}  LD scores: {ld_fmt.upper()}\n"
                      f"  Both must use the same format.")
            id_format = bim_fmt
        else:
            from .io.bim_sanitizer import sanitize_bim

            bfile, id_format = sanitize_bim(
                bfile   = bfile,
                out_dir = step1_dir,
                ld_dir  = args.ref_ld_chr
            )
        log_pass(log, f"Variant ID format auto-detected: {id_format.upper()}")
        log_info(log,  f"VCF conversion will use {id_format.upper()} variant IDs to match BIM/LD scores")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 3 — Load BIM variant IDs into memory
    # ─────────────────────────────────────────────────────────────────────────
    bim_ids = set()
    if bfile:
        step_banner(log, "Loading BIM Variant IDs", step=3, total=6)
        bim_path = f"{bfile}.bim"
        bim_ids  = _load_bim_ids(bim_path)
        log_pass(log, f"Loaded {len(bim_ids):,} variant IDs from BIM")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 4 — VCF → GCTA .ma conversion (per trait)
    # ─────────────────────────────────────────────────────────────────────────
    rows          = []
    target_coords = None
    n_total       = len(df_man)
    ma_paths      = {}   # {trait_name: ma_path} — used by LDSC

    for idx, r in enumerate(df_man.iter_rows(named=True)):
        trait   = r["sample_id"]
        vcf     = r["file_path"]
        ma_path = os.path.join(step1_dir, f"{trait}.ma")

        step_banner(log,
                    f"VCF Conversion  [{idx+1}/{n_total}]  {trait}",
                    step=4, total=n_steps)

        reuse_ma = not args.force and _nonempty(ma_path)
        if reuse_ma:
            df_coords = _load_coords_cache(ma_path)
            log_pass(log, f"Reusing existing non-empty .ma file: {ma_path}")
            if len(df_coords) == 0:
                log_warn(log, f"No coordinate cache found for {trait}; downstream PostGWAS may use merged report coordinates if available.")
        else:
            from .io.vcf_converter import convert_vcf_single_pass

            df_coords = convert_vcf_single_pass(
                vcf_path    = vcf,
                ma_path     = ma_path,
                bcftools_bin= args.bcftools,
                id_format   = id_format,
                bim_ids     = bim_ids or None
            )
            reran_analysis = True

        if idx == 0:
            target_coords = df_coords

        ma_paths[trait] = ma_path   # record for LDSC

        has_prev = (
            "sample_prevalence" in r
            and r["sample_prevalence"] is not None
            and str(r["sample_prevalence"]).strip().upper() not in ("NA", "")
        )
        prev_str = f"\t{r['sample_prevalence']}\t{r['population_prevalence']}" if has_prev else ""
        rows.append(f"{trait}\t{ma_path}{prev_str}")

    # ── Write .mtcojo.list ────────────────────────────────────────────────────
    list_file = f"{out_prefix}.mtcojo.list"
    list_text = "\n".join(rows) + "\n"
    list_changed = True
    if _nonempty(list_file):
        with open(list_file) as f:
            list_changed = f.read() != list_text
    if args.force or not _nonempty(list_file) or list_changed:
        with open(list_file, "w") as f:
            f.write(list_text)
        log_info(log, f"Trait list written: {list_file}")
        reran_analysis = True
    else:
        log_pass(log, f"Reusing existing mtCOJO trait list: {list_file}")
    for r in rows:
        log_info(log, f"  {r}")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 5 — GCTA mtCOJO
    # ─────────────────────────────────────────────────────────────────────────
    cma_file = f"{out_prefix}.mtcojo.cma"   # used by PostGWAS and LDSC
    if not args.force and not list_changed and _nonempty(cma_file):
        step_banner(log, "GCTA mtCOJO Analysis  [RESUME: OUTPUT EXISTS]", step=5, total=n_steps)
        log_pass(log, f"Reusing existing non-empty mtCOJO output: {cma_file}")
    else:
        from .stages.gcta import run_gcta_mtcojo

        run_gcta_mtcojo(
            gcta64_bin   = args.gcta64,
            bfile        = bfile,
            list_file    = list_file,
            out_prefix   = out_prefix,
            ref_ld_chr   = args.ref_ld_chr,
            w_ld_chr     = args.w_ld_chr,
            mbfile       = args.mbfile,
            gwas_thresh  = args.gwas_thresh,
            clump_r2     = args.clump_r2,
            heidi_thresh = args.heidi_thresh,
            gsmr_snp_min = args.gsmr_snp_min,
            diff_freq    = args.diff_freq,
            mtcojo_bxy   = args.mtcojo_bxy,
            gsmr_ld_fdr  = args.gsmr_ld_fdr
        )
        reran_analysis = True

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 6 — PostGWAS (optional)
    # ─────────────────────────────────────────────────────────────────────────
    if args.run_postgwas:
        sdd_dir = os.path.abspath(args.sdd_dir or os.path.join(step3_dir, "sdd_disk"))
        if not args.force and _postgwas_outputs_exist(sdd_dir, out_name):
            step_banner(log, "PostGWAS Harmonisation  [RESUME: OUTPUT EXISTS]", step=6, total=n_steps)
            log_pass(log, f"Reusing existing PostGWAS harmonised outputs under: {os.path.join(sdd_dir, '03_harmonised_output')}")
        else:
            if target_coords is None or len(target_coords) == 0:
                target_coords = _load_target_coords_for_resume(ma_paths[traits[0]], merged_tsv)
            if target_coords is None or len(target_coords) == 0:
                abort(log, "PostGWAS resume needs target SNP coordinates, but no .coords.tsv or merged_gwas_summary.tsv was found. Re-run without relying on skipped VCF conversion, or use --force.")
            from .stages.postgwas import run_postgwas_harmonisation

            run_postgwas_harmonisation(
                cma_file        = cma_file,
                target_coords   = target_coords,
                out_prefix      = os.path.join(step3_dir, out_name),
                sdd_dir         = sdd_dir,
                defaults_yaml   = args.defaults,
                resource_folder = args.resource_folder,
                docker_image    = args.docker_image,
                docker_platform = args.docker_platform,
                nthreads        = args.nthreads,
                max_mem         = args.max_mem,
                liftover        = args.liftover
            )
            reran_analysis = True
    else:
        step_banner(log, "PostGWAS Harmonisation  [SKIPPED]")
        log_warn(log, "Pass --run-postgwas to enable the PostGWAS Docker harmonisation step.")

    # ─────────────────────────────────────────────────────────────────────────
    # STEP 7 — LDSC (optional)
    # ─────────────────────────────────────────────────────────────────────────
    ldsc_results_csv = None
    ldsc_h2_csv = os.path.join(step4_dir, "ldsc_h2_results.csv")
    if args.run_ldsc:
        # Validate required LDSC inputs
        if not args.ldsc_snp_list:
            abort(log,
                  "--ldsc-snp-list is required when --run-ldsc is set.\n"
                  "  Provide path to HapMap3 SNP list (w_hm3.snplist).\n"
                  "  Download: https://data.broadinstitute.org/alkesgroup/LDSCORE/w_hm3.snplist.bz2")

        # LDSC LD reference directory: prefer --ldsc-ld-dir, else fall back to --ref-ld-chr
        ldsc_ld_dir = args.ldsc_ld_dir or args.ref_ld_chr
        if not ldsc_ld_dir:
            abort(log,
                  "LDSC requires an LD reference directory.\n"
                  "  Pass --ldsc-ld-dir /path/to/eur_w_ld_chr  (or --ref-ld-chr if shared)")

        # Build trait manifest for LDSC: all input traits from manifest
        target_row  = df_man.row(0, named=True)
        trait_manifest_ldsc = []
        for r in df_man.iter_rows(named=True):
            trait_manifest_ldsc.append({
                "name"              : r["sample_id"],
                "ma_path"           : ma_paths[r["sample_id"]],
                "sample_prevalence" : r.get("sample_prevalence", "nan"),
                "pop_prevalence"    : r.get("population_prevalence", "nan"),
            })

        cma_trait_name = f"{traits[0]}_conditioned"

        ldsc_results_csv = os.path.join(step4_dir, "ldsc_results.csv")
        if not args.force and _nonempty(ldsc_results_csv) and _nonempty(ldsc_h2_csv):
            step_banner(log, "LDSC Genetic Correlation  [RESUME: OUTPUT EXISTS]", step=7, total=n_steps)
            log_pass(log, f"Reusing existing LDSC rg results: {ldsc_results_csv}")
            log_pass(log, f"Reusing existing LDSC h2 results: {ldsc_h2_csv}")
        else:
            from .stages.ldsc import run_ldsc_pipeline

            run_ldsc_pipeline(
                trait_manifest   = trait_manifest_ldsc,
                cma_path         = cma_file,
                cma_trait_name   = cma_trait_name,
                cma_sample_prev  = target_row.get("sample_prevalence", "nan"),
                cma_pop_prev     = target_row.get("population_prevalence", "nan"),
                out_prefix       = os.path.join(step4_dir, out_name),
                snp_include_file = args.ldsc_snp_list,
                ld_ref_dir       = ldsc_ld_dir,
                ldsc_dir         = args.ldsc_dir,
                n_parallel       = args.ldsc_n_parallel,
                batch_size       = args.ldsc_batch_size,
                center_z         = args.ldsc_center_z,
            )
            reran_analysis = True
    else:
        step_banner(log, "LDSC Genetic Correlation  [SKIPPED]")
        log_warn(log,
                 "Pass --run-ldsc --ldsc-snp-list /path/to/w_hm3.snplist to enable LDSC rg analysis.\n"
                 "  LDSC runs independently on .ma files + .mtcojo.cma (PostGWAS not required).\n"
                 "  To run PostGWAS harmonisation as well, also pass --run-postgwas with its required arguments.")

    # ─────────────────────────────────────────────────────────────────────────
    # Pipeline summary
    # ─────────────────────────────────────────────────────────────────────────
    elapsed = time.time() - t0
    summary_rows = [
        ("Target trait",       traits[0]),
        ("Covariate trait(s)", ", ".join(traits[1:])),
        ("Variant ID format",  id_format.upper()),
        ("Output directory",   out_dir),
        ("mtCOJO output",      f"{out_prefix}.mtcojo.cma"),
        ("PostGWAS",           "ENABLED" if args.run_postgwas else "SKIPPED"),
        ("LDSC rg",            "ENABLED" if args.run_ldsc else "SKIPPED"),
    ]
    if ldsc_results_csv and os.path.exists(ldsc_results_csv):
        summary_rows.append(("LDSC rg results", ldsc_results_csv))
    if args.run_ldsc and os.path.exists(ldsc_h2_csv):
        summary_rows.append(("LDSC h2 results", ldsc_h2_csv))
    # ── Always generate/update the final HTML report at the end of the pipeline ──
    cli_str = " ".join(sys.argv)
    planned_report = os.path.join(out_dir, f"{out_name}_report.html")
    if not args.force and not args.force_report and not reran_analysis and _nonempty(planned_report):
        report_html = planned_report
        log_pass(log, f"Reusing existing non-empty HTML report: {report_html}")
    else:
        from .reporting.report_generator import build_pipeline_report

        report_html = build_pipeline_report(
            manifest_csv     = args.manifest,
            ma_paths         = ma_paths,
            cma_file         = cma_file,
            ldsc_results_csv = ldsc_results_csv if ldsc_results_csv and os.path.exists(ldsc_results_csv) else None,
            ldsc_h2_csv      = ldsc_h2_csv if args.run_ldsc and os.path.exists(ldsc_h2_csv) else None,
            cli_command      = cli_str,
            config_params    = vars(args),
            out_dir          = out_dir,
            out_prefix       = out_prefix,
            bcftools_bin     = args.bcftools
        )

    if os.path.exists(report_html):
        summary_rows.append(("Final HTML report", report_html))
    summary_rows += [
        ("Pipeline log",       log_file),
        ("Total time",         f"{elapsed:.1f} s"),
    ]
    summary_table(log, summary_rows, title="Pipeline Complete ✔")


if __name__ == "__main__":
    main()
