"""
ldsc_runner.py - LDSC Genetic Correlation and Heritability Module for mtcojo_postgwas

Runs heritability (--h2) and genetic correlation (--rg) analyses using:
  - All input GWAS traits (original summary statistics, from .ma files)
  - The mtCOJO-conditioned target trait (.mtcojo.cma)

LDSC installation: https://github.com/CBIIT/ldsc  (Python 3 fork)
Installed by install.sh into: tools/ldsc/
"""

import os
import sys
import re
import glob
import math
import gzip
import shutil
import base64
import subprocess
import concurrent.futures
import tempfile

import polars as pl

from .logger import get_logger, step_banner, log_pass, log_warn, log_info, abort, log_cmd_script

log = get_logger()


def _load_plotting():
    """Load plotting libraries only when a Python plot fallback is needed."""
    try:
        mpl_config_dir = os.path.join(tempfile.gettempdir(), "mtcojo_postgwas_matplotlib")
        os.makedirs(mpl_config_dir, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", mpl_config_dir)
        os.environ.setdefault("XDG_CACHE_HOME", mpl_config_dir)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        return plt, np
    except ImportError as e:
        raise RuntimeError(
            "Python plotting fallback requires matplotlib and numpy. "
            "Install them in this conda environment or use the R plotting path."
        ) from e


# ─────────────────────────────────────────────────────────────────────────────
# Locate ldsc scripts
# ─────────────────────────────────────────────────────────────────────────────

def _find_ldsc_scripts(ldsc_dir: str) -> tuple:
    """
    Locate munge_sumstats.py and ldsc.py inside the CBIIT/ldsc clone.
    """
    pkg_dir     = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    paths_file  = os.path.join(pkg_dir, "tools", "ldsc_paths.txt")
    cfg         = {}
    if os.path.exists(paths_file):
        with open(paths_file) as f:
            for line in f:
                k, _, v = line.strip().partition("=")
                cfg[k.strip()] = v.strip()

    python_bin  = (
        cfg.get("LDSC_PYTHON")
        or shutil.which("python3")
        or shutil.which("python")
    )

    candidates = [d for d in [
        ldsc_dir,
        cfg.get("LDSC_DIR"),
        os.path.join(pkg_dir, "tools", "ldsc"),
    ] if d]

    munge_script = cfg.get("MUNGE_SCRIPT")
    ldsc_script  = cfg.get("LDSC_SCRIPT")

    for d in candidates:
        if d and os.path.isdir(d):
            m = os.path.join(d, "munge_sumstats.py")
            l = os.path.join(d, "ldsc.py")
            if os.path.exists(m) and os.path.exists(l):
                munge_script = m
                ldsc_script  = l
                break

    missing = []
    if not python_bin:
        missing.append("python3 interpreter")
    if not munge_script or not os.path.exists(munge_script):
        missing.append(f"munge_sumstats.py (looked in: {candidates})")
    if not ldsc_script or not os.path.exists(ldsc_script):
        missing.append(f"ldsc.py (looked in: {candidates})")

    if missing:
        abort(log,
              f"LDSC scripts not found:\n" +
              "\n".join(f"  - {m}" for m in missing) +
              f"\n\n  Run:  bash install.sh\n"
              f"  Or pass:  --ldsc-dir /path/to/CBIIT_ldsc_clone")

    log_info(log, f"  ldsc scripts found:")
    log_info(log, f"    python        : {python_bin}")
    log_info(log, f"    munge_sumstats: {munge_script}")
    log_info(log, f"    ldsc.py       : {ldsc_script}")

    return python_bin, munge_script, ldsc_script


# ─────────────────────────────────────────────────────────────────────────────
# Stage 1 — Summary statistic → munge input TSV
# ─────────────────────────────────────────────────────────────────────────────

def ma_to_munge_input(ma_path: str, out_tsv: str, center_z: bool = False) -> bool:
    """
    Convert a GCTA .ma summary statistics file to munge_sumstats.py input format.
    """
    try:
        df = pl.read_csv(ma_path, separator="\t")
        required = {"SNP", "A1", "A2", "freq", "b", "se", "p", "N"}
        if not required.issubset(df.columns):
            log_warn(log, f"  .ma missing columns: {required - set(df.columns)}  →  {ma_path}")
            return False

        df = (
            df
            .filter(pl.col("se") > 0)
            .with_columns(
                (pl.col("b") / pl.col("se")).alias("EZ")
            )
            .filter(pl.col("EZ").is_finite())
        )
        
        if center_z and len(df) > 0:
            median_ez = df["EZ"].median()
            df = df.with_columns(pl.col("EZ") - median_ez)
            log_info(log, f"  Centered Z-scores for {os.path.basename(ma_path)} (median shifted by {-median_ez:.4f})")

        df = df.select([
            pl.col("SNP").alias("ID"),
            pl.col("A1").alias("ALT"),
            pl.col("A2").alias("REF"),
            pl.col("freq").alias("AF"),
            pl.col("EZ"),
            pl.col("p").alias("P"),
            pl.col("N").alias("NEF"),
        ])
        df.write_csv(out_tsv, separator="\t")
        log_info(log, f"  Munge input: {os.path.basename(out_tsv)}  ({len(df):,} variants)")
        return True
    except Exception as e:
        log_warn(log, f"  Failed to convert .ma → munge input [{ma_path}]: {e}")
        return False


def cma_to_munge_input(cma_path: str, out_tsv: str, use_conditional: bool = True, center_z: bool = False) -> bool:
    """
    Convert a GCTA .mtcojo.cma file to munge_sumstats.py input format.
    """
    try:
        df = pl.read_csv(
            cma_path,
            separator="\t",
            truncate_ragged_lines=True,
            null_values=["nan", "NA", ""]
        )
        required = {"SNP", "A1", "A2", "freq", "b", "se", "p", "N"}
        if not required.issubset(df.columns):
            log_warn(log, f"  .cma missing columns: {required - set(df.columns)}")
            return False

        has_bC = (
            use_conditional
            and {"bC", "bC_se", "bC_pval"}.issubset(df.columns)
            and df.filter(pl.col("bC").is_not_null()).height > 0
        )
        if has_bC:
            b_col, se_col, p_col = "bC", "bC_se", "bC_pval"
            label = "conditional (bC)"
        else:
            b_col, se_col, p_col = "b", "se", "p"
            label = "marginal (b)"
            if use_conditional:
                log_warn(log,
                         "  bC column is all-null (insufficient GW-sig SNPs for GSMR) "
                         "— using marginal b/se/p")

        df = (
            df
            .drop_nulls(subset=[b_col, se_col, p_col])
            .filter(pl.col(se_col).cast(pl.Float64, strict=False) > 0)
            .with_columns(
                (pl.col(b_col).cast(pl.Float64) / pl.col(se_col).cast(pl.Float64)).alias("EZ")
            )
            .filter(pl.col("EZ").is_finite())
        )

        if center_z and len(df) > 0:
            median_ez = df["EZ"].median()
            df = df.with_columns(pl.col("EZ") - median_ez)
            log_info(log, f"  Centered Z-scores for {os.path.basename(cma_path)} (median shifted by {-median_ez:.4f})")

        df = df.select([
            pl.col("SNP").alias("ID"),
            pl.col("A1").alias("ALT"),
            pl.col("A2").alias("REF"),
            pl.col("freq").cast(pl.Float64).alias("AF"),
            pl.col("EZ"),
            pl.col(p_col).cast(pl.Float64).alias("P"),
            pl.col("N").cast(pl.Float64).alias("NEF"),
        ])
        df.write_csv(out_tsv, separator="\t")
        log_info(log,
                 f"  Munge input: {os.path.basename(out_tsv)}  "
                 f"({len(df):,} variants)  [{label}]")
        return True
    except Exception as e:
        log_warn(log, f"  Failed to convert .cma → munge input [{cma_path}]: {e}")
        return False


# ─────────────────────────────────────────────────────────────────────────────
# Helpers
# ─────────────────────────────────────────────────────────────────────────────

def _is_valid_gz(path: str) -> bool:
    try:
        if not os.path.exists(path) or os.path.getsize(path) == 0:
            return False
        with gzip.open(path, "rb") as f:
            f.read(256)
        return True
    except Exception:
        return False


def _run_script(cmd: list, label: str, log_path: str = None, out_prefix: str = None) -> None:
    env = os.environ.copy()
    env["PYTHONNOUSERSITE"] = "1"
    if out_prefix:
        log_cmd_script(out_prefix, f"LDSC Analysis [{label}]", " ".join(cmd))
    result = subprocess.run(
        cmd,
        stdout=subprocess.PIPE,
        stderr=subprocess.STDOUT,
        text=True,
        env=env
    )
    if log_path:
        with open(log_path, "w") as f:
            f.write(result.stdout)
    if result.returncode != 0:
        tail = "\n".join(result.stdout.splitlines()[-25:])
        raise RuntimeError(
            f"[{label}] failed (exit {result.returncode}).\n"
            f"  Output tail:\n{tail}"
        )


# ─────────────────────────────────────────────────────────────────────────────
# Stage 2 — Parallel munge_sumstats.py
# ─────────────────────────────────────────────────────────────────────────────

def parallel_munge_sumstats(
    traits: list,
    ldsc_input_dir: str,
    snp_include_file: str,
    python_bin: str,
    munge_script: str,
    n_parallel: int = 4,
) -> list:
    step_banner(log, "LDSC Stage 2: munge_sumstats.py (parallel)")
    os.makedirs(ldsc_input_dir, exist_ok=True)

    if not os.path.exists(snp_include_file):
        abort(log, f"HapMap3 SNP list not found: {snp_include_file}")

    def munge_one(trait: dict) -> tuple:
        name       = trait["name"]
        in_file    = trait["munge_tsv"]
        out_prefix = os.path.join(ldsc_input_dir, name)
        out_gz     = f"{out_prefix}.sumstats.gz"
        munge_log  = f"{out_prefix}.munge.log"

        if _is_valid_gz(out_gz):
            log_info(log, f"  Skipping {name} — {out_gz} already exists")
            return name, True

        cmd = [
            python_bin, munge_script,
            "--sumstat",        in_file,
            "--N-col",          "NEF",
            "--snp",            "ID",
            "--a1",             "ALT",
            "--a2",             "REF",
            "--p",              "P",
            "--frq",            "AF",
            "--maf-min",        "0.005",
            "--signed-sumstats","EZ,0",
            "--merge-alleles",  snp_include_file,
            "--out",            out_prefix,
        ]
        log_info(log, f"  Munging: {name}")
        try:
            _run_script(cmd, f"munge_{name}", log_path=munge_log)
            if _is_valid_gz(out_gz):
                log_pass(log, f"  Munged: {name}  →  {os.path.basename(out_gz)}")
                return name, True
            else:
                log_warn(log, f"  Munge produced no .sumstats.gz for {name}")
                return name, False
        except RuntimeError as e:
            log_warn(log, f"  Munge failed: {name}\n  {e}")
            return name, False

    success = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as pool:
        futures = {pool.submit(munge_one, t): t["name"] for t in traits}
        for fut in concurrent.futures.as_completed(futures):
            name, ok = fut.result()
            if ok:
                success.append(name)

    log_pass(log, f"Munge complete: {len(success)}/{len(traits)} traits succeeded")
    return success


# ─────────────────────────────────────────────────────────────────────────────
# Stage 3a — Parallel ldsc.py --h2 (heritability)
# ─────────────────────────────────────────────────────────────────────────────

def parallel_ldsc_h2(
    manifest: list,
    ldsc_input_dir: str,
    ldsc_h2_dir: str,
    ld_ref_dir: str,
    python_bin: str,
    ldsc_script: str,
    n_parallel: int = 4,
) -> None:
    step_banner(log, "LDSC Stage 3a: --h2 heritability estimation (parallel)")
    os.makedirs(ldsc_h2_dir, exist_ok=True)

    def run_h2_one(row: dict) -> None:
        name       = row["name"]
        sumstats   = os.path.join(ldsc_input_dir, f"{name}.sumstats.gz")
        out_prefix = os.path.join(ldsc_h2_dir, name)
        h2_log     = f"{out_prefix}.log"

        if os.path.exists(h2_log) and os.path.getsize(h2_log) > 0:
            log_info(log, f"  Skipping heritability for {name} (log exists)")
            return

        cmd = [
            python_bin, ldsc_script,
            "--h2",         sumstats,
            "--ref-ld-chr", f"{ld_ref_dir}/",
            "--w-ld-chr",   f"{ld_ref_dir}/",
            "--out",        out_prefix,
        ]
        
        s_prev = str(row.get("sample_prevalence", "nan"))
        p_prev = str(row.get("pop_prevalence", "nan"))
        if s_prev != "nan" and p_prev != "nan":
            cmd += ["--samp-prev", s_prev, "--pop-prev", p_prev]

        log_info(log, f"  h2: {name}")
        try:
            _run_script(cmd, f"ldsc_h2_{name}", log_path=h2_log)
            log_pass(log, f"  Estimated h2: {name}")
        except RuntimeError as e:
            log_warn(log, f"  Failed h2 estimation for {name}: {e}")

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as pool:
        futures = [pool.submit(run_h2_one, r) for r in manifest]
        concurrent.futures.wait(futures)


def compile_ldsc_h2_results(ldsc_h2_dir: str, out_csv: str) -> pl.DataFrame:
    step_banner(log, "LDSC Stage 4a: Compile Heritability Results")
    log_files = glob.glob(os.path.join(ldsc_h2_dir, "*.log"))
    
    rows = []
    
    p_obs  = re.compile(r"Total Observed [sS]cale h2:\s*([-\d.e]+)\s*\(([-\d.e]+)\)")
    p_liab = re.compile(r"Total Liability [sS]cale h2:\s*([-\d.e]+)\s*\(([-\d.e]+)\)")
    p_int  = re.compile(r"Intercept:\s*([-\d.e]+)\s*\(([-\d.e]+)\)")
    p_lam  = re.compile(r"Lambda GC:\s*([-\d.e]+)")
    p_chi  = re.compile(r"Mean Chi\^2:\s*([-\d.e]+)")
    p_snps = re.compile(r"After merging with.*,\s*(\d+)\s*SNPs remain")

    for lf in sorted(log_files):
        name = os.path.basename(lf).replace(".log", "")
        try:
            with open(lf) as f:
                content = f.read()
            
            obs_match  = p_obs.search(content)
            liab_match = p_liab.search(content)
            int_match  = p_int.search(content)
            lam_match  = p_lam.search(content)
            chi_match  = p_chi.search(content)
            snps_match = p_snps.search(content)
            
            obs_h2, obs_h2_se = (float(obs_match.group(1)), float(obs_match.group(2))) if obs_match else (None, None)
            liab_h2, liab_h2_se = (float(liab_match.group(1)), float(liab_match.group(2))) if liab_match else (None, None)
            intercept, intercept_se = (float(int_match.group(1)), float(int_match.group(2))) if int_match else (None, None)
            lambda_gc = float(lam_match.group(1)) if lam_match else None
            mean_chi2 = float(chi_match.group(1)) if chi_match else None
            n_snps = int(snps_match.group(1)) if snps_match else None
            
            rows.append({
                "trait": name,
                "n_snps": n_snps,
                "h2_obs": obs_h2,
                "h2_obs_se": obs_h2_se,
                "h2_liab": liab_h2,
                "h2_liab_se": liab_h2_se,
                "intercept": intercept,
                "intercept_se": intercept_se,
                "lambda_gc": lambda_gc,
                "mean_chi2": mean_chi2
            })
        except Exception as e:
            log_warn(log, f"  Failed parsing {os.path.basename(lf)}: {e}")
            
    if not rows:
        log_warn(log, "  No heritability results found.")
        return pl.DataFrame()
        
    df = pl.DataFrame(rows)
    df.write_csv(out_csv)
    log_pass(log, f"  Compiled heritability results → {out_csv}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Stage 3b — Parallel ldsc.py --rg
# ─────────────────────────────────────────────────────────────────────────────

def parallel_ldsc_rg(
    manifest: list,
    ldsc_input_dir: str,
    ldsc_results_dir: str,
    ld_ref_dir: str,
    python_bin: str,
    ldsc_script: str,
    n_parallel: int = 4,
    batch_size: int = 10,
) -> None:
    step_banner(log, "LDSC Stage 3b: --rg genetic correlation (all pairwise)")
    os.makedirs(ldsc_results_dir, exist_ok=True)

    if not os.path.isdir(ld_ref_dir):
        abort(log, f"LD reference directory not found: {ld_ref_dir}")

    n_traits = len(manifest)
    n_batches_per_ref = math.ceil((n_traits - 1) / batch_size) if n_traits > 1 else 0
    total_jobs = n_traits * n_batches_per_ref

    def run_rg_batch(ref_row: dict, targets: list, part: int) -> None:
        ref_name   = ref_row["name"]
        ref_file   = os.path.join(ldsc_input_dir, f"{ref_name}.sumstats.gz")
        out_prefix = os.path.join(ldsc_results_dir, f"{ref_name}_batch{part:03d}")
        rg_log     = f"{out_prefix}.log"

        if os.path.exists(rg_log) and os.path.getsize(rg_log) > 0:
            log_info(log, f"  Skipping {ref_name} batch {part} (log exists)")
            return

        target_files = ",".join(
            os.path.join(ldsc_input_dir, f"{t['name']}.sumstats.gz")
            for t in targets
        )
        rg_input = f"{ref_file},{target_files}"

        s_prevs = [str(ref_row.get("sample_prevalence", "nan"))] + \
                  [str(t.get("sample_prevalence", "nan")) for t in targets]
        p_prevs = [str(ref_row.get("pop_prevalence",    "nan"))] + \
                  [str(t.get("pop_prevalence",    "nan")) for t in targets]

        cmd = [
            python_bin, ldsc_script,
            "--rg",         rg_input,
            "--ref-ld-chr", f"{ld_ref_dir}/",
            "--w-ld-chr",   f"{ld_ref_dir}/",
            "--samp-prev",  ",".join(s_prevs),
            "--pop-prev",   ",".join(p_prevs),
            "--out",        out_prefix,
        ]
        log_info(log, f"  rg: {ref_name} vs batch {part}  ({len(targets)} targets)")
        try:
            _run_script(cmd, f"ldsc_rg_{ref_name}_b{part}", log_path=rg_log)
            log_pass(log, f"  Done: {ref_name} batch {part}")
        except RuntimeError as e:
            log_warn(log, f"  Failed: {ref_name} batch {part}: {e}")

    jobs = []
    for ref_row in manifest:
        others = [t for t in manifest if t["name"] != ref_row["name"]]
        for part, start in enumerate(range(0, len(others), batch_size)):
            jobs.append((ref_row, others[start:start + batch_size], part))

    with concurrent.futures.ThreadPoolExecutor(max_workers=n_parallel) as pool:
        futures = [pool.submit(run_rg_batch, r, t, p) for r, t, p in jobs]
        concurrent.futures.wait(futures)

    log_pass(log, f"LDSC rg jobs complete — logs in: {ldsc_results_dir}")


# ─────────────────────────────────────────────────────────────────────────────
# Stage 4b — Compile results
# ─────────────────────────────────────────────────────────────────────────────

def compile_ldsc_results(ldsc_results_dir: str, out_csv: str) -> pl.DataFrame:
    step_banner(log, "LDSC Stage 4b: Compile Correlation Results")

    log_files = glob.glob(os.path.join(ldsc_results_dir, "*.log"))
    all_rows, headers = [], []

    for lf in sorted(log_files):
        try:
            with open(lf) as fh:
                lines = fh.readlines()
            for i, line in enumerate(lines):
                if "gcov_int_se" in line:
                    if not headers:
                        headers = line.split()
                    for data_line in lines[i + 1:]:
                        parts = data_line.split()
                        if not parts or "Summary" in data_line:
                            break
                        if len(parts) == len(headers):
                            all_rows.append(parts)
                    break
        except Exception as e:
            log_warn(log, f"  Could not parse {os.path.basename(lf)}: {e}")

    if not all_rows or not headers:
        log_warn(log, "  No genetic correlation results found in log files.")
        return pl.DataFrame()

    df = pl.DataFrame(all_rows, schema={h: pl.Utf8 for h in headers}, orient="row")

    for col in ["p1", "p2"]:
        if col in df.columns:
            df = df.with_columns(
                pl.col(col).map_elements(
                    lambda x: os.path.basename(x).replace(".sumstats.gz", ""),
                    return_dtype=pl.Utf8
                )
            )

    df = df.filter(pl.col("p1") != "p1").unique()

    for c in [col for col in df.columns if col not in ("p1", "p2")]:
        try:
            df = df.with_columns(pl.col(c).cast(pl.Float64, strict=False))
        except Exception:
            pass

    df = df.sort(["p1", "p2"])
    df.write_csv(out_csv)
    log_pass(log, f"  Compiled {len(df)} genetic correlations → {out_csv}")
    return df


# ─────────────────────────────────────────────────────────────────────────────
# Robust Multi-Trait Plotting Suite (Matplotlib)
# ─────────────────────────────────────────────────────────────────────────────

def make_gwas_visualizations(
    gwas_dfs: dict,
    ldsc_rg_csv: str,
    out_dir: str,
    ma_paths: dict,
    cma_file: str,
    bim_file: str,
) -> dict:
    """
    Generate:
      1. Stacked Multi-Trait Manhattan plot (rectangular fashion for multiple traits).
      2. Multi-Trait Q-Q plot overlaid on the same grid.
      3. Heatmap of genetic correlations (LDSC rg).
    Tries calling plot_gwas.R for CMplot/rMVP, falls back to matplotlib.
    """
    plot_b64s = {"manhattan": "", "qq": "", "heatmap": ""}
    r_ok = False

    # ── Attempt to run CMplot / rMVP via Rscript ─────────────────────────────
    try:
        rscript_bin = shutil.which("Rscript")
        if rscript_bin:
            pkg_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
            r_script = os.path.join(pkg_dir, "mtcojo_postgwas", "plot_gwas.R")
            
            trait_files = []
            trait_names = []
            for name, path in ma_paths.items():
                if os.path.exists(path):
                    trait_files.append(path)
                    trait_names.append(name)
            if os.path.exists(cma_file):
                trait_files.append(cma_file)
                trait_names.append("Conditioned_Target")
                
            if len(trait_files) >= 1 and os.path.exists(bim_file):
                cmd = [
                    rscript_bin, r_script,
                    out_dir,
                    ",".join(trait_files),
                    ",".join(trait_names),
                    bim_file,
                    ldsc_rg_csv if ldsc_rg_csv and os.path.exists(ldsc_rg_csv) else "none"
                ]
                log_info(log, "  [R] Invoking Rscript to run CMplot / rMVP...")
                res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
                if res.returncode == 0:
                    man_png = os.path.join(out_dir, "Rectangular-Manhattan.multi_trait.png")
                    qq_png = os.path.join(out_dir, "QQ-Plot.multi_trait.png")
                    hm_png = os.path.join(out_dir, "ldsc_rg_heatmap.png")
                    
                    if os.path.exists(man_png):
                        with open(man_png, "rb") as f:
                            plot_b64s["manhattan"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
                    if os.path.exists(qq_png):
                        with open(qq_png, "rb") as f:
                            plot_b64s["qq"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
                    if os.path.exists(hm_png):
                        with open(hm_png, "rb") as f:
                            plot_b64s["heatmap"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
                            
                    log_pass(log, "  [R] Successfully generated plots using CMplot/rMVP")
                    r_ok = True
                else:
                    log_warn(log, f"  [R] Rscript failed (exit {res.returncode}):\n{res.stdout}")
    except Exception as e:
        log_warn(log, f"  [R] Failed to invoke R plotting script: {e}")

    if r_ok:
        return plot_b64s

    log_info(log, "  Falling back to Python-based Matplotlib visualizations...")
    try:
        plt, np = _load_plotting()
    except RuntimeError as e:
        log_warn(log, str(e))
        return plot_b64s

    # ── 1. Stacked Manhattan Plot Fallback ────────────────────────────────────
    try:
        active_gwas = {k: v for k, v in gwas_dfs.items() if len(v) > 0}
        if active_gwas:
            n_traits = len(active_gwas)
            fig, axes = plt.subplots(n_traits, 1, figsize=(10, 3 * n_traits), sharex=True)
            if n_traits == 1:
                axes = [axes]
                
            colors = ["#4197d8", "#9e9e9e"]
            
            for ax, (trait, df) in zip(axes, active_gwas.items()):
                df = df.with_columns((-pl.col("P").log(10)).alias("minus_log10_p"))
                df = df.filter(pl.col("minus_log10_p").is_finite() & (pl.col("minus_log10_p") >= 0))
                df = df.sort(["CHR", "BP"])
                
                chrs = df["CHR"].unique().sort()
                for idx, c in enumerate(chrs):
                    sub = df.filter(pl.col("CHR") == c)
                    ax.scatter(sub["BP"] / 1e6, sub["minus_log10_p"], c=colors[idx % len(colors)], s=8, alpha=0.75, edgecolors="none")
                
                ax.axhline(y=-math.log10(1e-5), color="black", linestyle=":", alpha=0.8)
                ax.axhline(y=-math.log10(1e-8), color="red", linestyle="--", alpha=0.8)
                
                sig_snps = df.filter(pl.col("P") < 1e-5)
                if len(sig_snps) > 0:
                    for idx, c in enumerate(chrs):
                        sub_sig = sig_snps.filter(pl.col("CHR") == c)
                        ax.scatter(sub_sig["BP"] / 1e6, sub_sig["minus_log10_p"], c="#d93025", s=15, edgecolors="none")
                
                ax.set_title(f"Manhattan Plot: {trait}", fontsize=11, fontweight="bold", pad=8)
                ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=9)
                ax.grid(True, linestyle=":", alpha=0.4)
            
            axes[-1].set_xlabel("Position (Mb)", fontsize=10)
            plt.tight_layout()
            
            out_man = os.path.join(out_dir, "multi_trait_manhattan.png")
            plt.savefig(out_man, dpi=150)
            plt.close()
            
            with open(out_man, "rb") as f:
                plot_b64s["manhattan"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
    except Exception as e:
        log_warn(log, f"Failed to generate stacked Manhattan plot fallback: {e}")

    # ── 2. Multi-Trait Q-Q Plot Fallback ──────────────────────────────────────
    try:
        active_gwas = {k: v for k, v in gwas_dfs.items() if len(v) > 0}
        if active_gwas:
            plt.figure(figsize=(6.5, 6))
            palette = ["#1a73e8", "#d93025", "#12b886", "#f5a623", "#9061f9"]
            
            max_observed = 0
            for idx, (trait, df) in enumerate(active_gwas.items()):
                p_obs = df["P"].cast(pl.Float64, strict=False).drop_nulls().sort().to_list()
                M = len(p_obs)
                if M == 0:
                    continue
                
                p_exp = [(i + 1) / (M + 1) for i in range(M)]
                log10_exp = -np.log10(p_exp)
                log10_obs = -np.log10(p_obs)
                
                max_observed = max(max_observed, max(log10_obs))
                
                plt.scatter(
                    log10_exp,
                    log10_obs,
                    label=trait,
                    color=palette[idx % len(palette)],
                    s=6,
                    alpha=0.8,
                    edgecolors="none"
                )
            
            max_val = max(10.0, max_observed + 1.0)
            plt.plot([0, max_val], [0, max_val], color="red", linestyle="--", alpha=0.7, label="Expected (null)")
            
            plt.xlim(0, max_val)
            plt.ylim(0, max_val)
            plt.xlabel(r"Expected $-\log_{10}(P)$", fontsize=10)
            plt.ylabel(r"Observed $-\log_{10}(P)$", fontsize=10)
            plt.title("Multi-Trait Q-Q Plot", fontsize=12, fontweight="bold", pad=12)
            plt.grid(True, linestyle=":", alpha=0.4)
            plt.legend(loc="upper left")
            plt.tight_layout()
            
            out_qq = os.path.join(out_dir, "multi_trait_qq.png")
            plt.savefig(out_qq, dpi=150)
            plt.close()
            
            with open(out_qq, "rb") as f:
                plot_b64s["qq"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
    except Exception as e:
        log_warn(log, f"Failed to generate Q-Q plot fallback: {e}")

    # ── 3. LDSC rg Heatmap Fallback ───────────────────────────────────────────
    try:
        if ldsc_rg_csv and os.path.exists(ldsc_rg_csv):
            df_rg = pl.read_csv(ldsc_rg_csv)
            if len(df_rg) > 0:
                traits = sorted(list(set(df_rg["p1"].to_list() + df_rg["p2"].to_list())))
                n_t = len(traits)
                
                mat = np.ones((n_t, n_t))
                for row in df_rg.iter_rows(named=True):
                    i = traits.index(row["p1"])
                    j = traits.index(row["p2"])
                    mat[i, j] = row["rg"]
                    mat[j, i] = row["rg"]
                
                fig, ax = plt.subplots(figsize=(6, 5))
                im = ax.imshow(mat, cmap="coolwarm", vmin=-1.0, vmax=1.0)
                
                for i in range(n_t):
                    for j in range(n_t):
                        ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", 
                                color="white" if abs(mat[i, j]) > 0.5 else "black",
                                fontweight="bold", fontsize=10)
                
                ax.set_xticks(range(n_t))
                ax.set_yticks(range(n_t))
                ax.set_xticklabels(traits, rotation=45, ha="right", fontsize=9)
                ax.set_yticklabels(traits, fontsize=9)
                
                fig.colorbar(im, ax=ax, label="Genetic Correlation (rg)")
                plt.title("LDSC Genetic Correlation Heatmap", fontsize=11, fontweight="bold", pad=12)
                plt.tight_layout()
                
                out_hm = os.path.join(out_dir, "ldsc_rg_heatmap.png")
                plt.savefig(out_hm, dpi=150)
                plt.close()
                
                with open(out_hm, "rb") as f:
                    plot_b64s["heatmap"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
    except Exception as e:
        log_warn(log, f"Failed to generate genetic correlation heatmap fallback: {e}")

    return plot_b64s


# ─────────────────────────────────────────────────────────────────────────────
# HTML Report Generator
# ─────────────────────────────────────────────────────────────────────────────

def get_sig_stats(df: pl.DataFrame, p_col: str) -> dict:
    total = len(df)
    if total == 0:
        return {"total": 0, "sig_005": 0, "pct_005": 0, "sig_1e5": 0, "pct_1e5": 0, "sig_1e8": 0, "pct_1e8": 0}
    try:
        p_series = df[p_col].cast(pl.Float64, strict=False).drop_nulls()
        s_005 = p_series.filter(p_series < 0.05).len()
        s_1e5 = p_series.filter(p_series < 1e-5).len()
        s_1e8 = p_series.filter(p_series < 1e-8).len()
    except Exception:
        s_005, s_1e5, s_1e8 = 0, 0, 0
    return {
        "total": total,
        "sig_005": s_005,
        "pct_005": (s_005 / total) * 100 if total > 0 else 0,
        "sig_1e5": s_1e5,
        "pct_1e5": (s_1e5 / total) * 100 if total > 0 else 0,
        "sig_1e8": s_1e8,
        "pct_1e8": (s_1e8 / total) * 100 if total > 0 else 0,
    }


def generate_html_report(
    out_dir: str,
    out_prefix: str,
    cma_file: str,
    ma_paths: dict,
    bim_file: str,
    pipeline_args: dict = None,
    ldsc_results_csv: str = None,
    ldsc_h2_csv: str = None,
) -> str:
    """
    Generate a beautiful HTML report summarizing all results.
    """
    step_banner(log, "Generating HTML Report")
    report_html = f"{out_prefix}_report.html"

    # ── Load BIM coordinates mapping ──────────────────────────────────────────
    coords = {}
    if os.path.exists(bim_file):
        try:
            df_bim = pl.read_csv(
                bim_file,
                separator="\t",
                has_header=False,
                new_columns=["CHR", "SNP", "GP", "BP", "A1", "A2"]
            )
            for r in df_bim.iter_rows(named=True):
                coords[r["SNP"]] = (int(r["CHR"]), int(r["BP"]))
        except Exception as e:
            log_warn(log, f"  Failed to load PLINK BIM coordinates: {e}")

    # Helper to map coordinates and create df for Manhattan
    def get_manhattan_df(df_in: pl.DataFrame, snp_col: str, p_col: str) -> pl.DataFrame:
        chrs, bps = [], []
        for snp in df_in[snp_col]:
            c_bp = coords.get(snp)
            if c_bp:
                chrs.append(c_bp[0])
                bps.append(c_bp[1])
            else:
                chrs.append(1)  # Fallback
                bps.append(0)
        return df_in.with_columns([
            pl.Series("CHR", chrs),
            pl.Series("BP", bps),
            pl.col(p_col).cast(pl.Float64, strict=False).alias("P")
        ])

    # Build collection of traits for plotting
    plot_dfs = {}
    stats_rows = ""

    # Process original input files (.ma)
    for name, ma_path in ma_paths.items():
        if os.path.exists(ma_path):
            try:
                df = pl.read_csv(ma_path, separator="\t")
                st = get_sig_stats(df, "p")
                stats_rows += f"""
                <tr>
                    <td><strong>{name} (Input GCTA .ma)</strong></td>
                    <td>{st['total']:,}</td>
                    <td>{st['sig_005']:,} ({st['pct_005']:.2f}%)</td>
                    <td>{st['sig_1e5']:,} ({st['pct_1e5']:.4f}%)</td>
                    <td>{st['sig_1e8']:,} ({st['pct_1e8']:.4f}%)</td>
                </tr>
                """
                plot_dfs[name] = get_manhattan_df(df, "SNP", "p")
            except Exception as e:
                log_warn(log, f"  Failed stats collection for input {name}: {e}")

    # Process conditioned target output (.cma)
    if os.path.exists(cma_file):
        try:
            df_cma = pl.read_csv(cma_file, separator="\t", truncate_ragged_lines=True)
            p_col = "bC_pval" if "bC_pval" in df_cma.columns and df_cma["bC_pval"].drop_nulls().len() > 0 else "p"
            st = get_sig_stats(df_cma, p_col)
            label = "Conditioned GCTA output (.mtcojo.cma)"
            stats_rows += f"""
            <tr style="background-color: #f1f3f4;">
                <td><strong>{label}</strong></td>
                <td>{st['total']:,}</td>
                <td>{st['sig_005']:,} ({st['pct_005']:.2f}%)</td>
                <td>{st['sig_1e5']:,} ({st['pct_1e5']:.4f}%)</td>
                <td>{st['sig_1e8']:,} ({st['pct_1e8']:.4f}%)</td>
            </tr>
            """
            plot_dfs["Conditioned Target"] = get_manhattan_df(df_cma, "SNP", p_col)
        except Exception as e:
            log_warn(log, f"  Failed stats collection for GCTA conditioned output: {e}")

    # Generate visual plots (Manhattan stack, Q-Q overlay, Heatmap)
    plots = make_gwas_visualizations(
        gwas_dfs    = plot_dfs,
        ldsc_rg_csv = ldsc_results_csv,
        out_dir     = out_dir,
        ma_paths    = ma_paths,
        cma_file    = cma_file,
        bim_file    = bim_file
    )

    # ── Expandable Parameters section ──────────────────────────────────────────
    params_html = ""
    if pipeline_args:
        params_html += "<h3>Command Line Execution</h3>"
        cmd_used = " ".join(subprocess.list2cmdline([arg]) for arg in sys.argv) if hasattr(sys, "argv") else "mtcojo-postgwas"
        params_html += f"<pre style='background: #333; color: #fff; padding: 12px; border-radius: 6px; overflow-x: auto;'>{cmd_used}</pre>"
        
        params_html += "<h3>Active Settings & Configurations</h3>"
        params_html += "<table><thead><tr><th>Parameter</th><th>Value</th></tr></thead><tbody>"
        for k, v in sorted(pipeline_args.items()):
            params_html += f"<tr><td><code>--{k.replace('_', '-')}</code></td><td>{v}</td></tr>"
        params_html += "</tbody></table>"
    else:
        params_html = "<p>Parameter logs not provided.</p>"

    # ── GCTA top SNPs ──────────────────────────────────────────────────────────
    cma_rows = ""
    if os.path.exists(cma_file):
        try:
            df_cma = pl.read_csv(cma_file, separator="\t", truncate_ragged_lines=True)
            p_col = "bC_pval" if "bC_pval" in df_cma.columns and df_cma["bC_pval"].drop_nulls().len() > 0 else "p"
            df_cma_sig = df_cma.sort(p_col).head(15)
            
            for row in df_cma_sig.iter_rows(named=True):
                # Safely parse conditional items
                bC_val = f"{row.get('bC'):.4f}" if row.get('bC') is not None else "N/A"
                bC_se_val = f"{row.get('bC_se'):.4f}" if row.get('bC_se') is not None else "N/A"
                bC_pval_val = f"{row.get('bC_pval'):.2e}" if row.get('bC_pval') is not None else "N/A"
                
                cma_rows += f"""
                <tr>
                    <td>{row['SNP']}</td>
                    <td>{row['A1']}</td>
                    <td>{row['A2']}</td>
                    <td>{row['freq']:.4f}</td>
                    <td>{row.get('b', 0):.4f}</td>
                    <td>{row.get('se', 0):.4f}</td>
                    <td>{row.get('p', 0):.2e}</td>
                    <td><strong>{bC_val}</strong></td>
                    <td>{bC_se_val}</td>
                    <td><strong>{bC_pval_val}</strong></td>
                </tr>
                """
        except Exception as e:
            cma_rows = f"<tr><td colspan='10'>Error loading mtCOJO results: {e}</td></tr>"
    else:
        cma_rows = "<tr><td colspan='10'>No mtCOJO output file found.</td></tr>"

    # ── LDSC heritability ──────────────────────────────────────────────────────
    h2_rows = ""
    if ldsc_h2_csv and os.path.exists(ldsc_h2_csv):
        try:
            df_h2 = pl.read_csv(ldsc_h2_csv)
            for row in df_h2.iter_rows(named=True):
                obs = f"{row['h2_obs']:.4f} ({row['h2_obs_se']:.4f})" if row.get('h2_obs') is not None else "N/A"
                liab = f"{row['h2_liab']:.4f} ({row['h2_liab_se']:.4f})" if row.get('h2_liab') is not None else "N/A"
                inter = f"{row['intercept']:.4f} ({row['intercept_se']:.4f})" if row.get('intercept') is not None else "N/A"
                n_snps_val = f"{row['n_snps']:,}" if row.get('n_snps') is not None else "0"
                lam = f"{row['lambda_gc']:.4f}" if row.get('lambda_gc') is not None else "N/A"
                chi = f"{row['mean_chi2']:.4f}" if row.get('mean_chi2') is not None else "N/A"
                h2_rows += f"""
                <tr>
                    <td><strong>{row['trait']}</strong></td>
                    <td>{n_snps_val}</td>
                    <td>{obs}</td>
                    <td>{liab}</td>
                    <td>{inter}</td>
                    <td>{lam}</td>
                    <td>{chi}</td>
                </tr>
                """
        except Exception as e:
            h2_rows = f"<tr><td colspan='7'>Error loading heritability results: {e}</td></tr>"
    else:
        h2_rows = "<tr><td colspan='7'>LDSC heritability not run or not compiled.</td></tr>"

    # ── LDSC genetic correlation ───────────────────────────────────────────────
    rg_rows = ""
    if ldsc_results_csv and os.path.exists(ldsc_results_csv):
        try:
            df_rg = pl.read_csv(ldsc_results_csv)
            for row in df_rg.iter_rows(named=True):
                rg_rows += f"""
                <tr>
                    <td>{row['p1']}</td>
                    <td>{row['p2']}</td>
                    <td><strong>{row['rg']:.4f}</strong></td>
                    <td>{row['se']:.4f}</td>
                    <td>{row['z']:.4f}</td>
                    <td>{row['p']:.2e}</td>
                    <td>{row['gcov_int']:.4f} ({row['gcov_int_se']:.4f})</td>
                </tr>
                """
        except Exception as e:
            rg_rows = f"<tr><td colspan='7'>Error loading genetic correlation results: {e}</td></tr>"
    else:
        rg_rows = "<tr><td colspan='7'>LDSC genetic correlation not run or not compiled.</td></tr>"

    # ── Display Visualizations (Manhattan stack, Q-Q overlay, Heatmap side-by-side) ──
    plots_html = ""
    if plots["manhattan"] or plots["qq"] or plots["heatmap"]:
        plots_html += "<div class='vis-grid'>"
        if plots["manhattan"]:
            plots_html += f"""
            <div class="vis-card" style="grid-column: span 2;">
                <h3>Manhattan Plot (Multiple Traits Stacked)</h3>
                <img src="{plots['manhattan']}" alt="Manhattan Plot Stack" style="width:100%; border-radius:6px; border:1px solid #dee2e6;" />
            </div>
            """
        if plots["qq"]:
            plots_html += f"""
            <div class="vis-card">
                <h3>Multi-Trait Q-Q Plot</h3>
                <img src="{plots['qq']}" alt="Multi-Trait Q-Q Plot" style="width:100%; border-radius:6px; border:1px solid #dee2e6;" />
            </div>
            """
        if plots["heatmap"]:
            plots_html += f"""
            <div class="vis-card">
                <h3>Genetic Correlation (rg) Heatmap</h3>
                <img src="{plots['heatmap']}" alt="LDSC rg Heatmap" style="width:100%; border-radius:6px; border:1px solid #dee2e6;" />
            </div>
            """
        plots_html += "</div>"
    else:
        plots_html = "<p>No visualizations generated.</p>"

    html_content = f"""<!DOCTYPE html>
<html>
<head>
    <meta charset="utf-8">
    <title>mtcojo_postgwas Final Report</title>
    <style>
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, "Helvetica Neue", Arial, sans-serif;
            background-color: #f8f9fa;
            color: #212529;
            margin: 0;
            padding: 40px;
        }}
        .container {{
            max-width: 1200px;
            margin: 0 auto;
            background: #ffffff;
            padding: 40px;
            border-radius: 12px;
            box-shadow: 0 4px 20px rgba(0, 0, 0, 0.05);
        }}
        h1 {{
            color: #1a73e8;
            font-size: 2.2rem;
            margin-top: 0;
            border-bottom: 2px solid #e8f0fe;
            padding-bottom: 20px;
        }}
        h2 {{
            color: #3c4043;
            font-size: 1.5rem;
            margin-top: 40px;
            margin-bottom: 20px;
        }}
        table {{
            width: 100%;
            border-collapse: collapse;
            margin-bottom: 30px;
        }}
        th, td {{
            padding: 12px 15px;
            text-align: left;
            border-bottom: 1px solid #dee2e6;
        }}
        th {{
            background-color: #f1f3f4;
            color: #5f6368;
            font-weight: 600;
        }}
        tr:hover {{
            background-color: #f8f9fa;
        }}
        details {{
            background: #f8f9fa;
            padding: 16px;
            border-radius: 8px;
            border: 1px solid #dee2e6;
            margin-bottom: 25px;
        }}
        details summary {{
            font-weight: bold;
            font-size: 1.1rem;
            cursor: pointer;
            color: #1a73e8;
        }}
        .vis-grid {{
            display: grid;
            grid-template-columns: 1fr 1fr;
            gap: 25px;
            margin-top: 20px;
        }}
        .vis-card {{
            background: #ffffff;
            border: 1px solid #dee2e6;
            border-radius: 8px;
            padding: 20px;
            box-shadow: 0 2px 8px rgba(0,0,0,0.03);
        }}
        .vis-card h3 {{
            margin-top: 0;
            margin-bottom: 15px;
            color: #202124;
            font-size: 1.1rem;
            border-bottom: 1px solid #f1f3f4;
            padding-bottom: 8px;
        }}
    </style>
</head>
<body>
    <div class="container">
        <h1>mtcojo_postgwas Final Report</h1>
        <p>Analysis completed successfully. Expand sections below to inspect configurations, statistical counts, and plots.</p>

        <!-- Expandable Execution Parameters -->
        <details>
            <summary>🔧 Pipeline Parameters & Command Used</summary>
            <div style="margin-top: 15px;">
                {params_html}
            </div>
        </details>

        <h2>1. Variant Abundance & Significance Counts</h2>
        <table>
            <thead>
                <tr>
                    <th>GWAS Trait / File</th>
                    <th>Total SNPs</th>
                    <th>Significant at 0.05 level (Pct)</th>
                    <th>Significant at 1e-5 level (Pct)</th>
                    <th>Significant at 1e-8 level (Pct)</th>
                </tr>
            </thead>
            <tbody>
                {stats_rows}
            </tbody>
        </table>

        <h2>2. Visualizations (Manhattan stack, Q-Q, & Correlation Heatmap)</h2>
        {plots_html}

        <h2>3. GCTA mtCOJO Top 15 SNPs</h2>
        <table>
            <thead>
                <tr>
                    <th>SNP</th>
                    <th>A1</th>
                    <th>A2</th>
                    <th>Freq</th>
                    <th>b</th>
                    <th>se</th>
                    <th>p</th>
                    <th>bC (Conditional)</th>
                    <th>bC_se</th>
                    <th>bC_pval</th>
                </tr>
            </thead>
            <tbody>
                {cma_rows}
            </tbody>
        </table>
        <p style="font-size: 0.9rem; color: #5f6368; font-style: italic;">
            * Note: Conditional (bC) values are N/A if GSMR estimation was bypassed (fewer than 10 genome-wide significant markers).
        </p>

        <h2>4. LDSC Heritability Summary</h2>
        <table>
            <thead>
                <tr>
                    <th>Trait</th>
                    <th>SNPs</th>
                    <th>Observed h2 (SE)</th>
                    <th>Liability h2 (SE)</th>
                    <th>Intercept (SE)</th>
                    <th>Lambda GC</th>
                    <th>Mean Chi^2</th>
                </tr>
            </thead>
            <tbody>
                {h2_rows}
            </tbody>
        </table>

        <h2>5. LDSC Genetic Correlation (Pairwise rg)</h2>
        <table>
            <thead>
                <tr>
                    <th>Trait 1</th>
                    <th>Trait 2</th>
                    <th>rg</th>
                    <th>SE</th>
                    <th>Z-score</th>
                    <th>P-value</th>
                    <th>Gen Cov Intercept (SE)</th>
                </tr>
            </thead>
            <tbody>
                {rg_rows}
            </tbody>
        </table>
    </div>
</body>
</html>
"""

    with open(report_html, "w") as f:
        f.write(html_content)

    log_pass(log, f"HTML Report generated successfully → {report_html}")
    return report_html


# ─────────────────────────────────────────────────────────────────────────────
# Top-level orchestrator — called from cli.py
# ─────────────────────────────────────────────────────────────────────────────

def run_ldsc_pipeline(
    trait_manifest: list,
    cma_path: str,
    cma_trait_name: str,
    cma_sample_prev,
    cma_pop_prev,
    out_prefix: str,
    snp_include_file: str,
    ld_ref_dir: str,
    ldsc_dir: str = "",
    n_parallel: int = 4,
    batch_size: int = 10,
    center_z: bool = False,
) -> None:
    """
    Orchestrate the full LDSC heritability and genetic correlation pipeline.
    """
    step_banner(log, "LDSC Genetic Correlation & Heritability Pipeline")

    out_dir          = os.path.dirname(out_prefix)
    munge_input_dir  = os.path.join(out_dir, "ldsc_munge_input")
    ldsc_input_dir   = os.path.join(out_dir, "ldsc_sumstats")
    ldsc_h2_dir      = os.path.join(out_dir, "ldsc_h2")
    ldsc_results_dir = os.path.join(out_dir, "ldsc_results")
    os.makedirs(munge_input_dir, exist_ok=True)

    python_bin, munge_script, ldsc_script = _find_ldsc_scripts(ldsc_dir)

    step_banner(log, "LDSC Stage 1b: Convert Summary Statistics → Munge Input")

    munge_traits = []
    for t in trait_manifest:
        tsv = os.path.join(munge_input_dir, f"{t['name']}_mungeinput.tsv")
        if ma_to_munge_input(t["ma_path"], tsv, center_z=center_z):
            munge_traits.append({
                "name"              : t["name"],
                "munge_tsv"         : tsv,
                "sample_prevalence" : t.get("sample_prevalence", "nan"),
                "pop_prevalence"    : t.get("pop_prevalence",    "nan"),
            })

    if os.path.exists(cma_path):
        cma_tsv = os.path.join(munge_input_dir, f"{cma_trait_name}_mungeinput.tsv")
        if cma_to_munge_input(cma_path, cma_tsv, use_conditional=True, center_z=center_z):
            munge_traits.append({
                "name"              : cma_trait_name,
                "munge_tsv"         : cma_tsv,
                "sample_prevalence" : cma_sample_prev,
                "pop_prevalence"    : cma_pop_prev,
            })
    else:
        log_warn(log, f"mtCOJO .cma not found: {cma_path} — conditioned output excluded from LDSC")

    if len(munge_traits) < 2:
        abort(log, f"Need ≥ 2 valid traits for LDSC (only {len(munge_traits)} prepared).")

    log_pass(log, f"Prepared {len(munge_traits)} munge inputs: {[t['name'] for t in munge_traits]}")

    succeeded = parallel_munge_sumstats(
        traits          = munge_traits,
        ldsc_input_dir  = ldsc_input_dir,
        snp_include_file= snp_include_file,
        python_bin      = python_bin,
        munge_script    = munge_script,
        n_parallel      = n_parallel,
    )

    ldsc_manifest = [t for t in munge_traits if t["name"] in set(succeeded)]
    if len(ldsc_manifest) < 2:
        abort(log, f"Only {len(ldsc_manifest)} trait(s) munged successfully — need ≥ 2.")

    parallel_ldsc_h2(
        manifest       = ldsc_manifest,
        ldsc_input_dir = ldsc_input_dir,
        ldsc_h2_dir    = ldsc_h2_dir,
        ld_ref_dir     = ld_ref_dir,
        python_bin     = python_bin,
        ldsc_script    = ldsc_script,
        n_parallel     = n_parallel,
    )
    ldsc_h2_csv = os.path.join(out_dir, "ldsc_h2_results.csv")
    compile_ldsc_h2_results(ldsc_h2_dir, ldsc_h2_csv)

    parallel_ldsc_rg(
        manifest         = ldsc_manifest,
        ldsc_input_dir   = ldsc_input_dir,
        ldsc_results_dir = ldsc_results_dir,
        ld_ref_dir       = ld_ref_dir,
        python_bin       = python_bin,
        ldsc_script      = ldsc_script,
        n_parallel       = n_parallel,
        batch_size       = batch_size,
    )

    out_csv = os.path.join(out_dir, "ldsc_results.csv")
    df      = compile_ldsc_results(ldsc_results_dir, out_csv)
