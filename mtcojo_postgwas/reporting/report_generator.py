#!/usr/bin/env python3
# =============================================================================
# mtcojo_postgwas/reporting.py
# Standalone Module for Summary Statistics, Visualizations, & HTML Report Generation
# =============================================================================

import os
import sys
import io
import math
import shutil
import base64
import subprocess
import html
import re
import tempfile
import polars as pl

from mtcojo_postgwas.core.logger import get_logger, log_info, log_pass, log_warn, log_fail, log_cmd_script, abort
log = get_logger()


def _load_plotting():
    """Load Python plotting libraries only when fallback plots are generated."""
    try:
        mpl_config_dir = os.path.join(tempfile.gettempdir(), "mtcojo_postgwas_matplotlib")
        os.makedirs(mpl_config_dir, exist_ok=True)
        os.environ.setdefault("MPLCONFIGDIR", mpl_config_dir)
        os.environ.setdefault("XDG_CACHE_HOME", mpl_config_dir)
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import numpy as np
        return matplotlib, plt, np
    except ImportError as e:
        raise RuntimeError(
            "Python plotting requires matplotlib and numpy. "
            "Install them in this conda environment, or rely on available R plot outputs."
        ) from e


def _esc(value) -> str:
    """HTML-escape values interpolated into the standalone report."""
    return html.escape("" if value is None else str(value), quote=True)


def _fmt_num(value, digits: int = 4) -> str:
    try:
        if value is None:
            return "N/A"
        f = float(value)
        if not math.isfinite(f):
            return "N/A"
        return f"{f:.{digits}f}"
    except Exception:
        return _esc(value)


def _fmt_sci(value) -> str:
    try:
        if value is None:
            return "N/A"
        f = float(value)
        if not math.isfinite(f):
            return "N/A"
        return f"{f:.2e}"
    except Exception:
        return _esc(value)


def _read_image_b64(path: str) -> str:
    if not path or not os.path.exists(path):
        return ""
    ext = os.path.splitext(path)[1].lower()
    mime = "jpeg" if ext in (".jpg", ".jpeg") else "png"
    with open(path, "rb") as f:
        return f"data:image/{mime};base64,{base64.b64encode(f.read()).decode('utf-8')}"


def _file_size(path: str) -> str:
    if not path or not os.path.exists(path):
        return "pending"
    size = os.path.getsize(path)
    for unit in ("B", "KB", "MB", "GB"):
        if size < 1024 or unit == "GB":
            return f"{size:.1f} {unit}" if unit != "B" else f"{size} B"
        size /= 1024


def _build_file_inventory(out_dir: str, planned_report: str = None) -> str:
    if not out_dir or not os.path.isdir(out_dir):
        return "<p><em>No generated files found yet.</em></p>"

    folder_descriptions = {
        ".": "Run root: final report and top-level run outputs",
        "00_manifest_and_logs": "Validated manifest, execution logs, and command records",
        "01_gcta_ma_conversion": "Converted GCTA .ma files, reference-panel helpers, and variant QC tables",
        "02_gcta_mtcojo_results": "GCTA mtCOJO input list, log, and conditioned result files",
        "03_postgwas_harmonisation": "Optional PostGWAS/gwas2vcf harmonisation outputs",
        "04_ldsc_analysis": "Optional LDSC munged summary statistics, logs, h2, rg, and heatmap inputs",
        "05_plots_and_tables": "Merged GWAS summary, plots, overlap outputs, and report-ready tables",
    }
    step_order = {name: idx for idx, name in enumerate(folder_descriptions)}

    def describe_file(path: str, rel: str) -> str:
        name = os.path.basename(path)
        if planned_report and os.path.abspath(path) == os.path.abspath(planned_report):
            return "Final standalone HTML report"
        if name.endswith(".pipeline.log") or name.endswith(".log"):
            return "Execution log"
        if name.endswith(".ma"):
            return "GCTA .ma summary statistics"
        if name.endswith(".cma"):
            return "GCTA mtCOJO conditioned output"
        if name.endswith(".csv") or name.endswith(".tsv"):
            return "Tabular output"
        if name.endswith((".png", ".jpg", ".jpeg")):
            return "Report plot image"
        if "duplicate" in name.lower():
            return "Duplicate variant QC table"
        if "indel" in name.lower():
            return "Indel/allele QC table"
        return "Generated pipeline artifact"

    dir_paths = []
    for root, dirs, files in os.walk(out_dir):
        dirs[:] = sorted(d for d in dirs if not d.startswith("."))
        if os.path.basename(root).startswith("."):
            continue
        if files or dirs or os.path.abspath(root) == os.path.abspath(out_dir):
            dir_paths.append(root)

    if not dir_paths:
        return "<p><em>No generated files found yet.</em></p>"

    def dir_sort_key(path: str) -> tuple:
        rel = os.path.relpath(path, out_dir)
        first = "." if rel == "." else rel.split(os.sep, 1)[0]
        return (
            step_order.get(first, len(step_order)),
            rel.count(os.sep),
            rel,
        )

    sections = ""
    for directory in sorted(dir_paths, key=dir_sort_key):
        rel_dir = os.path.relpath(directory, out_dir)
        display_dir = "." if rel_dir == "." else rel_dir
        child_dirs = [
            name for name in sorted(os.listdir(directory))
            if os.path.isdir(os.path.join(directory, name)) and not name.startswith(".")
        ]
        files = [
            name for name in sorted(os.listdir(directory))
            if os.path.isfile(os.path.join(directory, name)) and not name.startswith(".")
        ]
        rows = ""
        for name in child_dirs:
            path = os.path.join(directory, name)
            rel = os.path.relpath(path, out_dir)
            rows += f"""
            <tr>
                <td>Folder</td>
                <td><code class="path-scroll" title="{_esc(path)}">{_esc(name)}/</code></td>
                <td>-</td>
                <td>{_esc(folder_descriptions.get(rel, 'Pipeline output subfolder'))}</td>
            </tr>
            """
        for name in files:
            path = os.path.join(directory, name)
            rel = os.path.relpath(path, out_dir)
            rows += f"""
            <tr>
                <td>File</td>
                <td><code class="path-scroll" title="{_esc(path)}">{_esc(name)}</code></td>
                <td>{_esc(_file_size(path))}</td>
                <td>{_esc(describe_file(path, rel))}</td>
            </tr>
            """
        if not rows:
            rows = '<tr><td colspan="4"><em>No files in this folder yet.</em></td></tr>'

        sections += f"""
        <details class="folder-inventory">
            <summary>{_esc(display_dir)}/ <span style="color:#5f6368; font-weight:400;">({_esc(folder_descriptions.get(display_dir, 'Pipeline output folder'))})</span></summary>
            <table class="data-table compact-table">
                <thead><tr><th>Type</th><th>Name</th><th>Size</th><th>Description</th></tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </details>
        """

    return sections


def _parse_ldsc_h2_log(log_path: str) -> dict:
    if not log_path or not os.path.exists(log_path):
        return {}
    try:
        with open(log_path) as f:
            content = f.read()
    except Exception:
        return {}

    patterns = {
        "h2_obs": r"Total Observed [sS]cale h2:\s*([-\d.eE]+)\s*\(([-\d.eE]+)\)",
        "h2_liab": r"Total Liability [sS]cale h2:\s*([-\d.eE]+)\s*\(([-\d.eE]+)\)",
        "intercept": r"Intercept:\s*([-\d.eE]+)\s*\(([-\d.eE]+)\)",
        "lambda_gc": r"Lambda GC:\s*([-\d.eE]+)",
        "mean_chi2": r"Mean Chi\^2:\s*([-\d.eE]+)",
        "n_snps": r"After merging with.*,\s*(\d+)\s*SNPs remain",
    }
    parsed = {}
    for key, pattern in patterns.items():
        match = re.search(pattern, content)
        if not match:
            continue
        if key in ("h2_obs", "h2_liab", "intercept"):
            parsed[key] = float(match.group(1))
            parsed[f"{key}_se"] = float(match.group(2))
        elif key == "n_snps":
            parsed[key] = int(match.group(1))
        else:
            parsed[key] = float(match.group(1))
    return parsed


def _build_ldsc_rg_failure_html(out_dir: str, ldsc_results_csv: str = None) -> str:
    candidate_dirs = [
        os.path.join(out_dir, "ldsc_results"),
        os.path.join(out_dir, "04_ldsc_analysis", "ldsc_results"),
    ]
    if ldsc_results_csv:
        candidate_dirs.append(os.path.join(os.path.dirname(ldsc_results_csv), "ldsc_results"))

    log_paths = []
    seen = set()
    for directory in candidate_dirs:
        if not os.path.isdir(directory):
            continue
        for name in sorted(os.listdir(directory)):
            if not name.endswith(".log"):
                continue
            path = os.path.join(directory, name)
            abspath = os.path.abspath(path)
            if abspath not in seen:
                seen.add(abspath)
                log_paths.append(path)

    rows = ""
    for path in log_paths:
        try:
            with open(path, encoding="utf-8", errors="replace") as f:
                lines = f.readlines()
        except Exception:
            continue
        content = "".join(lines)
        if "Summary of Genetic Correlation Results" in content and "gcov_int_se" in content:
            continue
        error_lines = [
            line.strip()
            for line in lines
            if (
                "ERROR computing rg" in line or
                "Traceback" in line or
                "TypeError:" in line or
                "ValueError:" in line or
                "FloatingPointError:" in line or
                "WARNING: number of SNPs less" in line
            )
        ]
        if not error_lines and "Analysis finished" in content:
            error_lines = ["LDSC rg log finished without a parsable genetic-correlation results table."]
        if error_lines:
            trimmed = " | ".join(error_lines[:6])
            rows += f"""
            <tr>
                <td><code class="path-scroll" title="{_esc(path)}">{_esc(os.path.relpath(path, out_dir))}</code></td>
                <td>{_esc(trimmed)}</td>
            </tr>
            """

    if not rows:
        return ""

    return f"""
    <div class="alert-box">
        <strong>LDSC rg did not produce a parsable pairwise results table.</strong>
        The heatmap is only generated after at least one rg estimate is compiled.
    </div>
    <table class="data-table compact-table">
        <thead><tr><th>LDSC rg Log</th><th>Detected Issue</th></tr></thead>
        <tbody>{rows}</tbody>
    </table>
    """


# ─────────────────────────────────────────────────────────────────────────────
# 1. Summary Statistics Helper
# ─────────────────────────────────────────────────────────────────────────────

def compute_p_stats(df: pl.DataFrame, p_col: str) -> dict:
    """Compute count and percentage of SNPs at 0.05, 1e-5, and 1e-8 thresholds."""
    total = len(df)
    if total == 0 or p_col not in df.columns:
        return {"total": 0, "sig_005": 0, "pct_005": 0.0, "sig_1e5": 0, "pct_1e5": 0.0, "sig_1e8": 0, "pct_1e8": 0.0}
    try:
        p_series = df[p_col].cast(pl.Float64, strict=False).drop_nulls()
        s_005 = p_series.filter(p_series < 0.05).len()
        s_1e5 = p_series.filter(p_series < 1e-5).len()
        s_1e8 = p_series.filter(p_series < 1e-8).len()
        return {
            "total": total,
            "sig_005": s_005,
            "pct_005": (s_005 / total * 100.0) if total > 0 else 0.0,
            "sig_1e5": s_1e5,
            "pct_1e5": (s_1e5 / total * 100.0) if total > 0 else 0.0,
            "sig_1e8": s_1e8,
            "pct_1e8": (s_1e8 / total * 100.0) if total > 0 else 0.0,
        }
    except Exception:
        return {"total": total, "sig_005": 0, "pct_005": 0.0, "sig_1e5": 0, "pct_1e5": 0.0, "sig_1e8": 0, "pct_1e8": 0.0}


# ─────────────────────────────────────────────────────────────────────────────
# 2. VCF Query Helper (Extract CHR, BP, SNP, P-value from VCF)
# ─────────────────────────────────────────────────────────────────────────────

def extract_vcf_coords_and_p(vcf_path: str, bcftools_bin: str = "bcftools") -> pl.DataFrame:
    """Queries a GWAS VCF file using bcftools to extract CHR, BP, SNP, and P-value."""
    if not os.path.exists(vcf_path):
        return pl.DataFrame(schema={"CHR": pl.Int64, "BP": pl.Int64, "SNP": pl.String, "P": pl.Float64})
    
    cmd = [bcftools_bin, "query", "-f", "%CHROM\t%POS\t%ID\t[%LP]\n", vcf_path]
    try:
        res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True, check=True)
        df = pl.read_csv(
            io.StringIO(res.stdout),
            separator="\t",
            has_header=False,
            new_columns=["CHR", "BP", "SNP", "LP"],
            schema_overrides={"CHR": pl.String, "BP": pl.Int64, "SNP": pl.String, "LP": pl.Float64}
        )
        
        # Clean chromosome (strip 'chr')
        df = df.with_columns(
            pl.col("CHR").str.replace("chr", "", literal=True).cast(pl.Int64, strict=False).alias("CHR")
        )
        
        # Convert LP (-log10 P) to P-value: P = 10^(-LP)
        df = df.with_columns(
            (10.0 ** (-pl.col("LP"))).alias("P")
        )
        return df.filter(pl.col("CHR").is_not_null() & pl.col("BP").is_not_null())
    except Exception as e:
        log_warn(log, f"Failed to extract VCF data via bcftools for {vcf_path}: {e}")
        return pl.DataFrame(schema={"CHR": pl.Int64, "BP": pl.Int64, "SNP": pl.String, "P": pl.Float64})


# ─────────────────────────────────────────────────────────────────────────────
# 3. Multi-Trait Plotting Engine (R CMplot/rMVP with Matplotlib Fallback)
# ─────────────────────────────────────────────────────────────────────────────

def generate_multi_trait_plots(
    merged_tsv: str,
    ldsc_rg_csv: str,
    out_dir: str,
    plot_dfs_python: dict
) -> dict:
    """Generates Manhattan, Q-Q, and Heatmap plots using R plot_gwas.R or Matplotlib."""
    plot_b64s = {
        "manhattan": "", "manhattan_overlaid": "", "qq": "", "qq_combined": "", "heatmap": "",
        "files": {}, "notes": []
    }
    r_ok = False

    # ── Attempt Rscript plot_gwas.R ─────────────────────────────────────────
    try:
        rscript_bin = None
        if "CONDA_PREFIX" in os.environ:
            conda_r = os.path.join(os.environ["CONDA_PREFIX"], "bin", "Rscript")
            if os.path.exists(conda_r):
                rscript_bin = conda_r
        if not rscript_bin:
            rscript_bin = shutil.which("Rscript")

        if rscript_bin and os.path.exists(merged_tsv):
            r_script = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", "plot_gwas.R")
            
            cmd = [
                rscript_bin, r_script,
                out_dir,
                merged_tsv,
                ldsc_rg_csv if ldsc_rg_csv and os.path.exists(ldsc_rg_csv) else "none"
            ]
            log_info(log, "  [R] Invoking Rscript to run CMplot / rMVP...")
            res = subprocess.run(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
            if res.returncode == 0:
                for img_ext in ["jpg", "png"]:
                    man_f  = os.path.join(out_dir, f"Rectangular-Manhattan.multi_trait.{img_ext}")
                    ovr_f  = os.path.join(out_dir, f"Overlaid-Manhattan.multi_trait.{img_ext}")
                    qq_f   = os.path.join(out_dir, f"QQ-Plot.multi_trait.{img_ext}")
                    qqc_f  = os.path.join(out_dir, f"QQ-Plot.combined.{img_ext}")
                    
                    mime = "jpeg" if img_ext == "jpg" else "png"
                    if not plot_b64s["manhattan"] and os.path.exists(man_f):
                        plot_b64s["manhattan"] = _read_image_b64(man_f)
                        plot_b64s["files"]["manhattan"] = man_f
                    if not plot_b64s["manhattan_overlaid"] and os.path.exists(ovr_f):
                        plot_b64s["manhattan_overlaid"] = _read_image_b64(ovr_f)
                        plot_b64s["files"]["manhattan_overlaid"] = ovr_f
                    if not plot_b64s["qq"] and os.path.exists(qq_f):
                        plot_b64s["qq"] = _read_image_b64(qq_f)
                        plot_b64s["files"]["qq"] = qq_f
                    if not plot_b64s["qq_combined"] and os.path.exists(qqc_f):
                        plot_b64s["qq_combined"] = _read_image_b64(qqc_f)
                        plot_b64s["files"]["qq_combined"] = qqc_f
                            
                hm_png = os.path.join(out_dir, "ldsc_rg_heatmap.png")
                if os.path.exists(hm_png):
                    plot_b64s["heatmap"] = _read_image_b64(hm_png)
                    plot_b64s["files"]["heatmap"] = hm_png
                        
                log_pass(log, "  [R] Successfully generated plots using CMplot/rMVP")
                plot_b64s["notes"].append("R plotting completed successfully.")
                r_ok = True
            else:
                log_warn(log, f"  [R] Rscript output:\n{res.stdout}")
                plot_b64s["notes"].append("R plotting failed; Python fallback images were used where needed.")
    except Exception as e:
        log_warn(log, f"  [R] Failed to execute R script: {e}")
        plot_b64s["notes"].append("R plotting failed; Python fallback images were used where needed.")

    # ── Python Matplotlib Fallbacks for Any Missing Image Slots ──────────────
    try:
        matplotlib, plt, np = _load_plotting()
        active_gwas = {k: v for k, v in plot_dfs_python.items() if len(v) > 0}
        if active_gwas:
            clean_gwas = {}
            skipped_tracks = []
            for trait, df in active_gwas.items():
                required_cols = {"P", "CHR", "BP"}
                if not required_cols.issubset(set(df.columns)):
                    skipped_tracks.append(trait)
                    continue
                df_sub = (
                    df
                    .with_columns([
                        pl.col("P").cast(pl.Float64, strict=False).alias("P"),
                        pl.col("CHR").cast(pl.Int64, strict=False).alias("CHR"),
                        pl.col("BP").cast(pl.Float64, strict=False).alias("BP"),
                    ])
                    .filter(
                        pl.col("P").is_not_null() &
                        (pl.col("P") > 0) &
                        pl.col("CHR").is_not_null() &
                        pl.col("BP").is_not_null()
                    )
                    .with_columns(
                        pl.when(pl.col("P") < 1e-300)
                        .then(1e-300)
                        .when(pl.col("P") > 1.0)
                        .then(1.0)
                        .otherwise(pl.col("P"))
                        .alias("P")
                    )
                    .with_columns((-pl.col("P").log(10)).alias("minus_log10_p"))
                    .filter(
                        pl.col("minus_log10_p").is_finite() &
                        (pl.col("minus_log10_p") >= 0) &
                        (pl.col("minus_log10_p") <= 300)
                    )
                    .with_row_index("__plot_row")
                    .filter((pl.col("P") < 0.1) | ((pl.col("__plot_row") % 20) == 0))
                    .drop("__plot_row")
                    .sort(["CHR", "BP"])
                )
                if len(df_sub) == 0:
                    skipped_tracks.append(trait)
                else:
                    clean_gwas[trait] = df_sub

            chr_max_bp = {}
            for df_sub in clean_gwas.values():
                for row in df_sub.group_by("CHR").agg(pl.col("BP").max().alias("max_bp")).iter_rows(named=True):
                    chrom = int(row["CHR"])
                    max_bp = float(row["max_bp"])
                    chr_max_bp[chrom] = max(chr_max_bp.get(chrom, 0.0), max_bp)

            chr_order = sorted(chr_max_bp)
            chr_offsets = {}
            chr_ticks = []
            chr_tick_labels = []
            running_bp = 0.0
            for chrom in chr_order:
                chr_offsets[chrom] = running_bp
                chr_len = chr_max_bp[chrom]
                chr_ticks.append((running_bp + chr_len / 2.0) / 1e6)
                chr_tick_labels.append(str(chrom))
                running_bp += chr_len + 1e6

            if clean_gwas:
                clean_gwas = {
                    trait: df_sub.with_columns(
                        pl.col("CHR")
                        .map_elements(lambda c: chr_offsets.get(int(c), 0.0), return_dtype=pl.Float64)
                        .alias("__chr_offset")
                    ).with_columns((pl.col("BP") + pl.col("__chr_offset")).alias("plot_pos"))
                    for trait, df_sub in clean_gwas.items()
                }

            # 1. Polished custom stacked Manhattan plot. Prefer this over R output.
            if clean_gwas:
                n_traits = len(clean_gwas)
                fig, axes = plt.subplots(
                    n_traits,
                    1,
                    figsize=(14, max(4.2, 3.45 * n_traits)),
                    sharex=True,
                    facecolor="white",
                )
                axes = np.atleast_1d(axes)
                band_colors = ["#3F3F46", "#A1A1AA"]
                trait_gradients = [
                    ["#3F3F46", "#71717A", "#A1A1AA"],
                    ["#0072B2", "#56B4E9", "#A1A1AA"],
                    ["#D55E00", "#E69F00", "#A1A1AA"],
                    ["#009E73", "#56B4E9", "#3F3F46"],
                ]
                max_y = max(float(df["minus_log10_p"].max()) for df in clean_gwas.values())
                y_limit = max(5.4, min(max_y * 1.18 + 0.4, 30.0))
                for trait_idx, (ax, (trait, df_sub)) in enumerate(zip(axes, clean_gwas.items())):
                    ax.set_facecolor("#fbfcfe")
                    for spine in ax.spines.values():
                        spine.set_color("#d0d5dd")
                        spine.set_linewidth(0.8)
                    chrs = df_sub["CHR"].unique().sort()
                    if len(chrs) == 1:
                        cmap = matplotlib.colors.LinearSegmentedColormap.from_list(
                            f"mtcojo_{trait_idx}",
                            trait_gradients[trait_idx % len(trait_gradients)],
                        )
                        ax.scatter(
                            df_sub["plot_pos"] / 1e6,
                            df_sub["minus_log10_p"],
                            c=df_sub["plot_pos"],
                            cmap=cmap,
                            s=0.7,
                            alpha=0.76,
                            edgecolors="none",
                        )
                    else:
                        for idx, chrom in enumerate(chrs):
                            sub = df_sub.filter(pl.col("CHR") == chrom)
                            ax.scatter(
                                sub["plot_pos"] / 1e6,
                                sub["minus_log10_p"],
                                color=band_colors[idx % len(band_colors)],
                                s=0.3,
                                alpha=0.68,
                                edgecolors="none",
                            )
                    ax.axhline(y=-math.log10(1e-5), color="#334155", linestyle=(0, (1, 2)), linewidth=1.2, alpha=0.9)
                    ax.axhline(y=-math.log10(1e-8), color="#BE123C", linestyle="--", linewidth=1.4, alpha=0.9)
                    ax.set_title(trait.replace("_", " "), fontsize=12, fontweight="bold", pad=7, color="#101828")
                    ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=9, color="#344054")
                    ax.set_ylim(0, y_limit)
                    ax.grid(True, axis="y", linestyle=":", linewidth=0.8, color="#d0d5dd", alpha=0.75)
                    ax.grid(True, axis="x", linestyle=":", linewidth=0.6, color="#eaecf0", alpha=0.75)
                    ax.tick_params(axis="both", labelsize=8, colors="#475467")
                    if running_bp > 0:
                        ax.set_xlim(0, running_bp / 1e6)
                axes[-1].set_xticks(chr_ticks)
                axes[-1].set_xticklabels(chr_tick_labels)
                axes[-1].set_xlabel("Chromosome", fontsize=10, color="#344054")
                fig.suptitle("Multi-Trait Manhattan Plot", fontsize=16, fontweight="bold", color="#101828", y=0.995)
                fig.tight_layout(rect=[0, 0, 1, 0.965])
                out_man = os.path.join(out_dir, "multi_trait_manhattan.png")
                fig.savefig(out_man, dpi=200, bbox_inches="tight")
                plt.close(fig)
                plot_b64s["manhattan"] = _read_image_b64(out_man)
                plot_b64s["files"]["manhattan"] = out_man
                plot_b64s["notes"].append("Custom polished Python Manhattan plot used for the stacked tracks; plotting p-values were clamped to [1e-300, 1], and 95% of variants with P between 0.1 and 1 were thinned for memory-efficient rendering.")
            if skipped_tracks:
                plot_b64s["notes"].append(f"Skipped empty Manhattan track(s): {', '.join(skipped_tracks)}.")

            # 2. Overlaid Multi-Trait Manhattan Plot. Prefer this custom version
            # over R output so trait colours and point size remain readable.
            if clean_gwas:
                fig, ax = plt.subplots(figsize=(12.5, 5.3), facecolor="white")
                palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]
                for idx, (trait, df_sub) in enumerate(clean_gwas.items()):
                    ax.scatter(
                        df_sub["plot_pos"] / 1e6,
                        df_sub["minus_log10_p"],
                        label=trait.replace("_", " "),
                        color=palette[idx % len(palette)],
                        s=0.3,
                        alpha=0.32,
                        edgecolors="none",
                        rasterized=True,
                    )
                ax.axhline(y=-math.log10(1e-5), color="#111827", linestyle=(0, (1, 2)), linewidth=1.1, alpha=0.9)
                ax.axhline(y=-math.log10(1e-8), color="#BE123C", linestyle="--", linewidth=1.2, alpha=0.9)
                ax.set_title("Multi-Trait Overlaid Manhattan Plot", fontsize=13, fontweight="bold", color="#101828", pad=10)
                ax.set_xlabel("Chromosome", fontsize=10, color="#344054")
                ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=10, color="#344054")
                if running_bp > 0:
                    ax.set_xlim(0, running_bp / 1e6)
                ax.set_xticks(chr_ticks)
                ax.set_xticklabels(chr_tick_labels)
                ax.grid(True, linestyle=":", linewidth=0.8, color="#d0d5dd", alpha=0.7)
                ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.92, markerscale=2.4)
                ax.tick_params(axis="both", labelsize=8, colors="#475467")
                for spine in ax.spines.values():
                    spine.set_color("#d0d5dd")
                fig.tight_layout()
                out_ovr = os.path.join(out_dir, "multi_trait_overlaid_manhattan.png")
                fig.savefig(out_ovr, dpi=190, bbox_inches="tight")
                plt.close(fig)
                plot_b64s["manhattan_overlaid"] = _read_image_b64(out_ovr)
                plot_b64s["files"]["manhattan_overlaid"] = out_ovr
                plot_b64s["notes"].append("Custom Python overlaid Manhattan plot used with cumulative chromosome positions, smaller points, reduced opacity, and muted phenotype colours.")

            # 3. Multi-Trait 3-Panel / Separate Q-Q Plot (if missing)
            if not plot_b64s["qq"]:
                plt.figure(figsize=(6, 5.5))
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
                    plt.scatter(log10_exp, log10_obs, label=trait, color=palette[idx % len(palette)], s=6, alpha=0.8, edgecolors="none")
                max_val = max(10.0, max_observed + 1.0)
                plt.plot([0, max_val], [0, max_val], color="red", linestyle="--", alpha=0.7, label="Expected (null)")
                plt.xlim(0, max_val)
                plt.ylim(0, max_val)
                plt.xlabel(r"Expected $-\log_{10}(P)$", fontsize=9)
                plt.ylabel(r"Observed $-\log_{10}(P)$", fontsize=9)
                plt.title("Multi-Trait Q-Q Plot", fontsize=11, fontweight="bold", pad=10)
                plt.grid(True, linestyle=":", alpha=0.4)
                plt.legend(loc="upper left", fontsize=8)
                plt.tight_layout()
                out_qq = os.path.join(out_dir, "multi_trait_qq.png")
                plt.savefig(out_qq, dpi=150)
                plt.close()
                plot_b64s["qq"] = _read_image_b64(out_qq)
                plot_b64s["files"]["qq"] = out_qq

            # 4. Multi-Trait Combined Single-Panel Q-Q Plot (if missing)
            if not plot_b64s["qq_combined"]:
                fig, ax = plt.subplots(figsize=(7.2, 6.2), facecolor="white")
                palette = ["#255C99", "#D00000", "#2A9D8F", "#F77F00", "#7B2CBF", "#0081A7"]
                max_observed = 0.0
                max_expected = 0.0
                plotted = 0
                for idx, (trait, df) in enumerate(active_gwas.items()):
                    if "P" not in df.columns:
                        continue
                    p_obs = df["P"].cast(pl.Float64, strict=False).drop_nulls()
                    p_vals = [p for p in p_obs.to_list() if p is not None and math.isfinite(float(p)) and float(p) > 0]
                    p_vals = sorted(p_vals)
                    n_obs = len(p_vals)
                    if n_obs == 0:
                        continue
                    p_exp = np.arange(1, n_obs + 1, dtype=float) / (n_obs + 1)
                    log10_exp = -np.log10(p_exp)
                    log10_obs = -np.log10(np.array(p_vals, dtype=float))
                    max_expected = max(max_expected, float(np.max(log10_exp)))
                    max_observed = max(max_observed, float(np.max(log10_obs)))
                    ax.scatter(
                        log10_exp,
                        log10_obs,
                        label=trait.replace("_", " "),
                        color=palette[idx % len(palette)],
                        s=12,
                        alpha=0.62,
                        edgecolors="none",
                    )
                    plotted += 1
                if plotted > 0:
                    max_val = max(3.0, min(max(max_observed, max_expected) + 0.35, 30.0))
                    ax.plot([0, max_val], [0, max_val], color="#111827", linestyle="--", linewidth=1.2, alpha=0.75, label="Expected")
                    ax.set_xlim(0, max_val)
                    ax.set_ylim(0, max_val)
                    ax.set_xlabel(r"Expected $-\log_{10}(P)$", fontsize=10, color="#344054")
                    ax.set_ylabel(r"Observed $-\log_{10}(P)$", fontsize=10, color="#344054")
                    ax.set_title("Multi-Trait Combined Q-Q Plot", fontsize=13, fontweight="bold", color="#101828", pad=10)
                    ax.grid(True, linestyle=":", linewidth=0.8, color="#d0d5dd", alpha=0.75)
                    ax.legend(loc="upper left", fontsize=8, frameon=True, framealpha=0.92)
                    for spine in ax.spines.values():
                        spine.set_color("#d0d5dd")
                    fig.tight_layout()
                    out_qq_combined = os.path.join(out_dir, "multi_trait_qq_combined.png")
                    fig.savefig(out_qq_combined, dpi=180, bbox_inches="tight")
                    plt.close(fig)
                    plot_b64s["qq_combined"] = _read_image_b64(out_qq_combined)
                    plot_b64s["files"]["qq_combined"] = out_qq_combined
                else:
                    plt.close(fig)

    except Exception as e:
        log_warn(log, f"Failed to generate Matplotlib fallbacks: {e}")
        plot_b64s["notes"].append(str(e))
        return plot_b64s

    # 5. LDSC rg Heatmap. Prefer this custom Python heatmap over R output so
    # labels have enough margin in the standalone HTML report.
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
                fig, ax = plt.subplots(figsize=(7.2, 6.4), facecolor="white")
                im = ax.imshow(mat, cmap="coolwarm", vmin=-1.0, vmax=1.0)
                for i in range(n_t):
                    for j in range(n_t):
                        ax.text(j, i, f"{mat[i, j]:.3f}", ha="center", va="center", 
                                color="white" if abs(mat[i, j]) > 0.5 else "black",
                                fontweight="bold", fontsize=9)
                ax.set_xticks(range(n_t))
                ax.set_yticks(range(n_t))
                ax.set_xticklabels(traits, rotation=45, ha="right", fontsize=9, rotation_mode="anchor")
                ax.set_yticklabels(traits, fontsize=9)
                fig.colorbar(im, ax=ax, label="Genetic Correlation (rg)")
                plt.title("LDSC Genetic Correlation Heatmap", fontsize=11, fontweight="bold", pad=10)
                fig.subplots_adjust(left=0.22, right=0.92, bottom=0.30, top=0.88)
                out_hm = os.path.join(out_dir, "ldsc_rg_heatmap.png")
                plt.savefig(out_hm, dpi=170, bbox_inches="tight", pad_inches=0.18)
                plt.close()
                plot_b64s["heatmap"] = _read_image_b64(out_hm)
                plot_b64s["files"]["heatmap"] = out_hm
    except Exception as e:
        log_warn(log, f"Python fallback plotting failed: {e}")
        plot_b64s["notes"].append(str(e))

    return plot_b64s


# ─────────────────────────────────────────────────────────────────────────────
# 4. Variant Overlap Venn Diagram & UpSet Plot Engine
# ─────────────────────────────────────────────────────────────────────────────

def generate_venn_and_upset_plots(plot_dfs: dict, out_dir: str) -> dict:
    """
    Generates multi-level Venn diagrams and UpSet plots for variant ID intersections
    across Total Variants, P < 0.05, P < 1e-5, and P < 1e-8.
    """
    levels = [
        ("total", "Total Variants", None),
        ("p005",  "Nominal Significance (P < 0.05)", 0.05),
        ("p1e5",  "Suggestive Significance (P < 10⁻⁵)", 1e-5),
        ("p1e8",  "Genome-Wide Significance (P < 10⁻⁸)", 1e-8),
    ]
    try:
        _, plt, np = _load_plotting()
        plotting_available = True
    except RuntimeError as e:
        log_warn(log, str(e))
        plotting_available = False
    venn_available = False
    venn2 = venn3 = None
    if plotting_available:
        try:
            from matplotlib_venn import venn2, venn3
            venn_available = True
        except ImportError:
            log_warn(
                log,
                "matplotlib-venn is not installed; Venn images will be skipped. "
                "Overlap count tables and UpSet plots will still be generated. "
                "Install with: conda install -n mtcojo_postgwas matplotlib-venn"
            )
    
    results = {}
    
    for lvl_key, lvl_title, p_cutoff in levels:
        variant_sets = {}
        for label, df in plot_dfs.items():
            if "SNP" in df.columns and len(df) > 0:
                sub_df = df
                if p_cutoff is not None and "P" in df.columns:
                    sub_df = df.filter(pl.col("P").cast(pl.Float64, strict=False) < p_cutoff)
                snps = set(sub_df["SNP"].drop_nulls().to_list())
                if len(snps) > 0:
                    variant_sets[label] = snps
                    
        size_rows = "".join(
            f"<tr><td>{_esc(label.replace('_', ' '))}</td><td>{len(snps):,}</td></tr>"
            for label, snps in variant_sets.items()
        )
        intersection_rows = ""
        if len(variant_sets) >= 2:
            labels_for_counts = list(variant_sets.keys())[:3]
            for i, label_a in enumerate(labels_for_counts):
                for label_b in labels_for_counts[i + 1:]:
                    n_shared = len(variant_sets[label_a] & variant_sets[label_b])
                    intersection_rows += (
                        f"<tr><td>{_esc(label_a.replace('_', ' '))} / "
                        f"{_esc(label_b.replace('_', ' '))}</td><td>{n_shared:,}</td></tr>"
                    )
            if len(labels_for_counts) >= 3:
                n_shared_all = len(set.intersection(*(variant_sets[k] for k in labels_for_counts)))
                intersection_rows += f"<tr><td>All displayed datasets</td><td>{n_shared_all:,}</td></tr>"
        if not size_rows:
            size_rows = "<tr><td colspan='2'><em>No variants available at this threshold.</em></td></tr>"
        if not intersection_rows:
            intersection_rows = "<tr><td colspan='2'><em>At least two non-empty datasets are required.</em></td></tr>"
        summary_html = f"""
        <details class="overlap-summary">
            <summary>Overlap Counts Table</summary>
            <div class="overlap-table-grid">
                <table class="data-table compact-table">
                    <thead><tr><th>Dataset</th><th>Variants</th></tr></thead>
                    <tbody>{size_rows}</tbody>
                </table>
                <table class="data-table compact-table">
                    <thead><tr><th>Intersection</th><th>Shared Variants</th></tr></thead>
                    <tbody>{intersection_rows}</tbody>
                </table>
            </div>
        </details>
        """

        lvl_res = {"venn": "", "upset": "", "summary": summary_html, "title": lvl_title, "n_sets": len(variant_sets)}
        
        if len(variant_sets) >= 2 and plotting_available:
            # A. Venn Diagram
            if venn_available:
                try:
                    fig, ax = plt.subplots(figsize=(7.8, 6.3), facecolor="white")
                    labels = list(variant_sets.keys())[:3]
                    sets = [variant_sets[k] for k in labels]
                    
                    if len(sets) == 2:
                        venn = venn2(sets, set_labels=[lbl.replace("_", " ") for lbl in labels], ax=ax)
                    elif len(sets) >= 3:
                        venn = venn3(sets[:3], set_labels=[lbl.replace("_", " ") for lbl in labels[:3]], ax=ax)
                    else:
                        venn = None
                    if venn:
                        palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#3F3F46"]
                        patch_idx = 0
                        for patch in venn.patches:
                            if patch:
                                patch.set_facecolor(palette[patch_idx % len(palette)])
                                patch.set_alpha(0.56)
                                patch.set_edgecolor("#ffffff")
                                patch.set_linewidth(1.2)
                                patch_idx += 1
                        for text in list(venn.set_labels or []) + list(venn.subset_labels or []):
                            if text:
                                text.set_fontsize(10)
                                text.set_color("#101828")
                                text.set_fontweight("bold")
                    ax.set_title(lvl_title, fontsize=13, fontweight="bold", color="#101828", pad=12)
                    plt.tight_layout()
                    venn_png = os.path.join(out_dir, f"venn_{lvl_key}.png")
                    plt.savefig(venn_png, dpi=180)
                    plt.close()
                    
                    with open(venn_png, "rb") as f:
                        lvl_res["venn"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
                except Exception as e:
                    log_warn(log, f"Failed Venn [{lvl_key}]: {e}")

            # B. UpSet Plot
            try:
                labels = list(variant_sets.keys())[:3]
                import itertools
                combo_counts = {}
                for r in range(1, len(labels) + 1):
                    for combo in itertools.combinations(labels, r):
                        in_combo = set.intersection(*(variant_sets[k] for k in combo))
                        others = [variant_sets[k] for k in labels if k not in combo]
                        out_combo = set.union(*others) if others else set()
                        exact_snps = in_combo - out_combo
                        if len(exact_snps) > 0:
                            combo_counts[combo] = len(exact_snps)
                            
                sorted_combos = sorted(combo_counts.keys(), key=lambda c: combo_counts[c], reverse=True)
                counts = [combo_counts[c] for c in sorted_combos]
                
                if sorted_combos:
                    fig = plt.figure(figsize=(8.5, 5.5))
                    gs = fig.add_gridspec(2, 1, height_ratios=[3, 1.2], hspace=0.15)
                    ax_bar = fig.add_subplot(gs[0])
                    ax_matrix = fig.add_subplot(gs[1], sharex=ax_bar)
                    
                    x_positions = np.arange(len(sorted_combos))
                    bars = ax_bar.bar(x_positions, counts, color="#1b4965", width=0.45)
                    max_c = max(counts) if counts else 1
                    for x_pos, cnt in zip(x_positions, counts):
                        ax_bar.text(x_pos, cnt + max_c * 0.02, f"{cnt:,}", ha="center", va="bottom", fontsize=8, fontweight="bold")
                    ax_bar.set_ylabel("Intersection Size", fontsize=9, fontweight="bold")
                    ax_bar.grid(True, axis="y", linestyle=":", alpha=0.4)
                    ax_bar.tick_params(labelbottom=False)
                    
                    for y_idx, label_name in enumerate(labels):
                        for x_idx, combo in enumerate(sorted_combos):
                            present = label_name in combo
                            dot_color = "#1b4965" if present else "#d1d5db"
                            dot_size = 90 if present else 35
                            ax_matrix.scatter(x_idx, y_idx, color=dot_color, s=dot_size, zorder=3)
                        for x_idx, combo in enumerate(sorted_combos):
                            active_indices = [labels.index(lbl) for lbl in combo if lbl in labels]
                            if len(active_indices) > 1:
                                ax_matrix.plot([x_idx, x_idx], [min(active_indices), max(active_indices)], color="#1b4965", linewidth=2, zorder=2)
                                
                    ax_matrix.set_yticks(np.arange(len(labels)))
                    ax_matrix.set_yticklabels(labels, fontsize=9, fontweight="bold")
                    ax_matrix.set_ylim(-0.5, len(labels) - 0.5)
                    ax_matrix.grid(True, axis="y", linestyle="--", alpha=0.3)
                    ax_matrix.tick_params(bottom=False, labelbottom=False)
                    
                    upset_png = os.path.join(out_dir, f"upset_{lvl_key}.png")
                    plt.savefig(upset_png, dpi=180, bbox_inches="tight")
                    plt.close()
                    
                    with open(upset_png, "rb") as f:
                        lvl_res["upset"] = f"data:image/png;base64,{base64.b64encode(f.read()).decode('utf-8')}"
            except Exception as e:
                log_warn(log, f"Failed UpSet [{lvl_key}]: {e}")
                
        results[lvl_key] = lvl_res
        
    return results


def _build_conditional_shift_html(summary_tsv: str, detail_tsv: str) -> str:
    if not summary_tsv or not os.path.exists(summary_tsv):
        return "<p><em>Conditional significance shift summary was not generated.</em></p>"

    try:
        summary_df = pl.read_csv(summary_tsv, separator="\t")
    except Exception as e:
        return f"<p><em>Unable to read conditional shift summary: {_esc(e)}</em></p>"

    if len(summary_df) == 0:
        return "<p><em>No conditional significance shift rows were generated.</em></p>"

    summary_df = (
        summary_df
        .with_columns(
            pl.when(pl.col("comparison_trait").str.ends_with("_Target"))
            .then(0)
            .otherwise(1)
            .alias("__trait_order")
        )
        .sort(["__trait_order", "comparison_trait", "threshold"])
        .drop("__trait_order")
    )

    summary_rows = ""
    for row in summary_df.iter_rows(named=True):
        summary_rows += (
            f"<tr><td>{_esc(row.get('comparison_trait'))}</td>"
            f"<td>{_esc(row.get('conditioned_trait'))}</td>"
            f"<td>{_fmt_sci(row.get('threshold'))}</td>"
            f"<td>{_fmt_num(row.get('gwas_significant_n'), 0)}</td>"
            f"<td>{_fmt_num(row.get('conditioned_significant_n'), 0)}</td>"
            f"<td>{_fmt_num(row.get('lost_after_conditioning_n'), 0)}</td>"
            f"<td>{_fmt_num(row.get('gained_after_conditioning_n'), 0)}</td>"
            f"<td>{_fmt_num(row.get('both_significant_n'), 0)}</td></tr>"
        )

    detail_html = "<p><em>No gained/lost SNP examples were found at these thresholds.</em></p>"
    if detail_tsv and os.path.exists(detail_tsv):
        try:
            detail_df = (
                pl.scan_csv(detail_tsv, separator="\t")
                .with_columns(
                    pl.when(pl.col("direction") == "gained_after_conditioning")
                    .then(pl.col("conditioned_p"))
                    .otherwise(pl.col("gwas_p"))
                    .alias("__display_p")
                )
                .with_columns(
                    pl.when(pl.col("comparison_trait").str.ends_with("_Target"))
                    .then(0)
                    .otherwise(1)
                    .alias("__trait_order")
                )
                .sort(["__trait_order", "comparison_trait", "threshold", "direction", "__display_p"])
                .with_columns(
                    pl.col("__display_p")
                    .rank("ordinal")
                    .over(["comparison_trait", "threshold", "direction"])
                    .alias("__group_rank")
                )
                .filter(pl.col("__group_rank") <= 10)
                .drop(["__display_p", "__group_rank", "__trait_order"])
                .limit(80)
                .collect()
            )
            if len(detail_df) > 0:
                detail_rows = ""
                for row in detail_df.iter_rows(named=True):
                    detail_rows += (
                        f"<tr><td>{_esc(row.get('comparison_trait'))}</td>"
                        f"<td>{_fmt_sci(row.get('threshold'))}</td>"
                        f"<td>{_esc(row.get('direction'))}</td>"
                        f"<td>{_esc(row.get('SNP'))}</td>"
                        f"<td>{_esc(row.get('CHR'))}</td>"
                        f"<td>{_esc(row.get('BP'))}</td>"
                        f"<td>{_fmt_sci(row.get('gwas_p'))}</td>"
                        f"<td>{_fmt_sci(row.get('conditioned_p'))}</td></tr>"
                    )
                detail_html = f"""
                <details>
                    <summary>Top gained/lost SNP examples shown in this report</summary>
                    <table class="data-table">
                        <thead>
                            <tr>
                                <th>GWAS Trait</th><th>Threshold</th><th>Direction</th><th>SNP</th>
                                <th>CHR</th><th>BP</th><th>Original GWAS P</th>
                                <th>Conditioned P</th>
                            </tr>
                        </thead>
                        <tbody>{detail_rows}</tbody>
                    </table>
                </details>
                """
        except Exception as e:
            detail_html = f"<p><em>Unable to read conditional shift SNP details: {_esc(e)}</em></p>"

    return f"""
    <p style="color:#5f6368; font-size:14px; margin:12px 0 16px 0;">
        Counts compare each original target/covariate GWAS p-value column against the mtCOJO conditioned target.
        <strong>Lost</strong> means original GWAS P is at or below the threshold but conditioned P is above it;
        <strong>gained</strong> means the reverse.
    </p>
    <table class="data-table" style="margin-top:12px;">
        <thead>
            <tr>
                <th>GWAS Trait</th><th>Conditioned Trait</th><th>Threshold</th>
                <th>Original Significant</th><th>Conditioned Significant</th>
                <th>Lost After Conditioning</th><th>Gained After Conditioning</th><th>Both Significant</th>
            </tr>
        </thead>
        <tbody>{summary_rows}</tbody>
    </table>
    <p style="color:#5f6368; font-size:13px;">
        Full TSV outputs: <code>{_esc(summary_tsv)}</code> and <code>{_esc(detail_tsv)}</code>
    </p>
    {detail_html}
    """


# ─────────────────────────────────────────────────────────────────────────────
# 5. Main HTML Report Generator Entrypoint
# ─────────────────────────────────────────────────────────────────────────────

def build_pipeline_report(
    manifest_csv: str,
    ma_paths: dict,
    cma_file: str,
    ldsc_results_csv: str,
    ldsc_h2_csv: str,
    cli_command: str,
    config_params: dict,
    out_dir: str,
    out_prefix: str,
    bcftools_bin: str = "bcftools"
) -> str:
    """
    Builds a comprehensive HTML report containing:
      - Input VCF stats & percentages (P < 0.05, 1e-5, 1e-8)
      - GCTA Harmonized .ma stats & percentages
      - GCTA mtCOJO .cma stats & percentages
      - Merged wide TSV dataset (CHR, BP, SNP, P-values for all traits)
      - Stacked Manhattan, Q-Q, and LDSC Heatmap plots
      - Collapsible parameter and command logs
    """
    log_info(log, "Generating Standalone HTML Pipeline Report...")
    out_name = os.path.basename(out_prefix)
    if cli_command:
        log_cmd_script(out_prefix, "Pipeline CLI Command", cli_command)

    # Prefer explicit LDSC paths from the active run, but also discover completed
    # LDSC artifacts in the report output directory. This keeps regenerated HTML
    # from showing LDSC as skipped when the CSVs already exist beside the report.
    if not ldsc_results_csv or not os.path.exists(ldsc_results_csv):
        for candidate in (
            os.path.join(out_dir, "ldsc_results.csv"),
            os.path.join(out_dir, "04_ldsc_analysis", "ldsc_results.csv"),
        ):
            if os.path.exists(candidate):
                ldsc_results_csv = candidate
                break
    if not ldsc_h2_csv or not os.path.exists(ldsc_h2_csv):
        for candidate in (
            os.path.join(out_dir, "ldsc_h2_results.csv"),
            os.path.join(out_dir, "04_ldsc_analysis", "ldsc_h2_results.csv"),
        ):
            if os.path.exists(candidate):
                ldsc_h2_csv = candidate
                break
    ldsc_artifacts_available = bool(
        (ldsc_results_csv and os.path.exists(ldsc_results_csv)) or
        (ldsc_h2_csv and os.path.exists(ldsc_h2_csv))
    )

    # Read manifest to locate VCF files
    manifest_df = pl.read_csv(manifest_csv).filter(
        pl.col("sample_id").is_not_null() & pl.col("file_path").is_not_null()
    )
    vcf_dict = dict(zip(manifest_df["sample_id"].to_list(), manifest_df["file_path"].to_list()))
    step5_plots_dir = os.path.join(out_dir, "05_plots_and_tables")
    os.makedirs(step5_plots_dir, exist_ok=True)
    merged_tsv_path = os.path.join(step5_plots_dir, "merged_gwas_summary.tsv")

    stats_rows = ""
    plot_dfs = {}
    master_coords_df = None
    if os.path.exists(merged_tsv_path):
        try:
            previous_merged = pl.read_csv(merged_tsv_path, separator="\t")
            if {"SNP", "CHR", "BP"}.issubset(set(previous_merged.columns)) and len(previous_merged) > 0:
                master_coords_df = previous_merged.select(["SNP", "CHR", "BP"]).unique(subset=["SNP"])
        except Exception as e:
            log_warn(log, f"  Existing merged summary coordinates could not be reused: {e}")
    top_results_html = ""
    cma_top_html = "<p><em>mtCOJO output file was not found.</em></p>"

    def stats_row(label: str, st: dict, background: str) -> str:
        return f"""
        <tr style="background-color: {background};">
            <td><strong>{_esc(label)}</strong></td>
            <td>{st['total']:,}</td>
            <td>{st['sig_005']:,} ({st['pct_005']:.2f}%)</td>
            <td>{st['sig_1e5']:,} ({st['pct_1e5']:.4f}%)</td>
            <td>{st['sig_1e8']:,} ({st['pct_1e8']:.4f}%)</td>
        </tr>
        """

    def significant_table(
        title: str,
        df: pl.DataFrame,
        p_col: str,
        columns: list,
        threshold: float = 0.05,
        export_tsv: str = None,
        max_html_rows: int = 500,
    ) -> str:
        if p_col not in df.columns or len(df) == 0:
            return f"<p><em>No rows available for {_esc(title)}.</em></p>"
        try:
            sort_df = (
                df
                .with_columns(pl.col(p_col).cast(pl.Float64, strict=False).alias("__sort_p"))
                .drop_nulls(subset=["__sort_p"])
                .filter(pl.col("__sort_p") < threshold)
                .sort("__sort_p")
            )
            if export_tsv:
                sort_df.drop("__sort_p").write_csv(export_tsv, separator="\t")
            total_sig = len(sort_df)
            display_df = sort_df.head(max_html_rows)
        except Exception as e:
            return f"<p><em>Unable to build significant-variant table for {_esc(title)}: {_esc(e)}</em></p>"

        if total_sig == 0:
            return f"<p><em>No variants with P &lt; {_fmt_sci(threshold)} available for {_esc(title)}.</em></p>"

        header = "".join(f"<th>{_esc(c)}</th>" for c in columns)
        rows = ""
        for row in display_df.iter_rows(named=True):
            cells = ""
            for col in columns:
                value = row.get(col)
                if col.lower() in ("p", "bc_pval", "bC_pval".lower()):
                    cells += f"<td>{_fmt_sci(value)}</td>"
                elif col in ("freq", "b", "se", "bC", "bC_se", "N"):
                    cells += f"<td>{_fmt_num(value, 4)}</td>"
                else:
                    cells += f"<td>{_esc(value)}</td>"
            rows += f"<tr>{cells}</tr>"

        export_note = ""
        if export_tsv:
            export_note = (
                f'<p style="color:#5f6368; font-size:13px;">'
                f'Full P &lt; {_fmt_sci(threshold)} table: <code>{_esc(export_tsv)}</code></p>'
            )
        preview_note = ""
        if total_sig > len(display_df):
            preview_note = (
                f'<p style="color:#5f6368; font-size:13px;">'
                f'HTML preview shows the {len(display_df):,} most significant variants to keep the report responsive; '
                f'use the TSV for all {total_sig:,} variants.</p>'
            )

        return f"""
        <details open>
            <summary>{_esc(title)} ({total_sig:,} variants with P &lt; {_fmt_sci(threshold)})</summary>
            {preview_note}
            {export_note}
            <table class="data-table">
                <thead><tr>{header}</tr></thead>
                <tbody>{rows}</tbody>
            </table>
        </details>
        """

    # ── A. Process Input VCF Files (Target & Covariates) ─────────────────────
    log_info(log, "  [1/4] Extracting variant stats and coordinates from Input VCF files...")
    target_trait = manifest_df["sample_id"][0]
    
    for idx, row in enumerate(manifest_df.iter_rows(named=True)):
        trait_id = row["sample_id"]
        vcf_path = row["file_path"]
        vcf_df   = extract_vcf_coords_and_p(vcf_path, bcftools_bin)
        st       = compute_p_stats(vcf_df, "P")
        
        is_target = (idx == 0)
        label_text = f"{trait_id} (Target Input VCF)" if is_target else f"{trait_id} (Covariate Input VCF)"
        stats_rows += stats_row(label_text, st, "#f8f9fa")
        
        plot_label = f"{trait_id}_Target" if is_target else f"{trait_id}_Covariate"
        plot_dfs[plot_label] = vcf_df
        
        if len(vcf_df) > 0:
            sub_coords = vcf_df.select(["SNP", "CHR", "BP"])
            if master_coords_df is None:
                master_coords_df = sub_coords
            else:
                master_coords_df = pl.concat([master_coords_df, sub_coords]).unique(subset=["SNP"])

    # ── B. Process converted GCTA .ma files ──────────────────────────────────
    log_info(log, "  [2/4] Extracting stats from converted GCTA .ma files...")
    for idx, row in enumerate(manifest_df.iter_rows(named=True)):
        trait_id = row["sample_id"]
        role = "Target" if idx == 0 else "Covariate"
        ma_path = ma_paths.get(trait_id)
        if not ma_path or not os.path.exists(ma_path):
            stats_rows += stats_row(f"{trait_id} ({role} .ma missing)", compute_p_stats(pl.DataFrame(), "p"), "#fff4e5")
            top_results_html += f"<p><em>{_esc(trait_id)} .ma file not found.</em></p>"
            continue
        try:
            ma_df = pl.read_csv(ma_path, separator="\t")
            st = compute_p_stats(ma_df, "p")
            stats_rows += stats_row(f"{trait_id} ({role} converted .ma)", st, "#f3f8ff")
            top_results_html += significant_table(
                f"{trait_id} converted .ma nominally significant SNPs",
                ma_df,
                "p",
                ["SNP", "A1", "A2", "freq", "b", "se", "p", "N"],
                export_tsv=os.path.join(step5_plots_dir, f"{trait_id}.ma.p_lt_0.05.tsv"),
            )
            plot_label = f"{trait_id}_{role}"
            if master_coords_df is not None:
                plot_dfs[plot_label] = ma_df.join(master_coords_df, on="SNP", how="inner").rename({"p": "P"})
        except Exception as e:
            log_warn(log, f"  Failed to summarize .ma file [{ma_path}]: {e}")
            top_results_html += f"<p><em>Error loading {_esc(trait_id)} .ma: {_esc(e)}</em></p>"

    # ── B. Process GCTA mtCOJO .cma File ─────────────────────────────────────
    log_info(log, "  [3/4] Extracting stats from GCTA mtCOJO .cma file...")
    if os.path.exists(cma_file):
        cma_df = pl.read_csv(cma_file, separator="\t", truncate_ragged_lines=True)
        p_col = "bC_pval" if "bC_pval" in cma_df.columns and cma_df["bC_pval"].drop_nulls().len() > 0 else "p"
        st = compute_p_stats(cma_df, p_col)
        
        label_text = f"{target_trait}_conditioned_output (.cma)"
        stats_rows += stats_row(label_text, st, "#e8f0fe")
        cma_columns = ["SNP", "A1", "A2", "freq", "b", "se", "p"]
        for extra in ["bC", "bC_se", "bC_pval"]:
            if extra in cma_df.columns:
                cma_columns.append(extra)
        if "N" in cma_df.columns:
            cma_columns.append("N")
        cma_top_html = significant_table(
            "mtCOJO conditioned .cma nominally significant SNPs",
            cma_df,
            p_col,
            cma_columns,
            export_tsv=os.path.join(step5_plots_dir, f"{target_trait}.conditioned.p_lt_0.05.tsv"),
        )
        
        plot_label = f"{target_trait}_conditioned"
        if master_coords_df is not None:
            plot_dfs[plot_label] = cma_df.join(master_coords_df, on="SNP", how="inner").rename({p_col: "P"})
        else:
            plot_dfs[plot_label] = cma_df.rename({p_col: "P"})

    # ── Build Merged Wide Summary Stats TSV for CMplot ───────────────────────
    if master_coords_df is not None and len(master_coords_df) > 0:
        merged_wide = master_coords_df
        for trait_label, sub_df in plot_dfs.items():
            if "P" in sub_df.columns and len(sub_df) > 0:
                p_sub = sub_df.select(["SNP", "P"]).rename({"P": trait_label})
                merged_wide = merged_wide.join(p_sub, on="SNP", how="left")
        
        merged_wide.write_csv(merged_tsv_path, separator="\t")
        log_pass(log, f"  Created merged wide summary stats TSV → {merged_tsv_path}")
    else:
        merged_tsv_path = "none"

    conditional_shift_html = "<p><em>Conditional significance shift summary requires a merged GWAS summary table.</em></p>"
    if merged_tsv_path != "none":
        try:
            from mtcojo_postgwas.reporting.conditional_shift_summary import build_conditional_shift_summary

            shift_outputs = build_conditional_shift_summary(
                merged_tsv=merged_tsv_path,
                out_dir=step5_plots_dir,
            )
            conditional_shift_html = _build_conditional_shift_html(
                shift_outputs["summary_tsv"],
                shift_outputs["detail_tsv"],
            )
            log_pass(log, f"  Created conditional significance shift TSVs → {shift_outputs['summary_tsv']}")
        except Exception as e:
            log_warn(log, f"  Conditional significance shift summary failed: {e}")
            conditional_shift_html = f"<p><em>Conditional significance shift summary failed: {_esc(e)}</em></p>"

    # ── Generate Figures ─────────────────────────────────────────────────────
    plots = generate_multi_trait_plots(merged_tsv_path, ldsc_results_csv, step5_plots_dir, plot_dfs)
    overlaps = generate_venn_and_upset_plots(plot_dfs, step5_plots_dir)

    # ── Read LDSC Tables ─────────────────────────────────────────────────────
    rg_failure_html = _build_ldsc_rg_failure_html(out_dir, ldsc_results_csv)
    rg_table_html = rg_failure_html or "<p><em>LDSC genetic correlation analysis was not run or returned no results.</em></p>"
    if ldsc_results_csv and os.path.exists(ldsc_results_csv):
        df_rg = pl.read_csv(ldsc_results_csv)
        if len(df_rg) > 0:
            rg_rows = ""
            for row in df_rg.iter_rows(named=True):
                rg_val = f"{row['rg']:.4f}" if row['rg'] is not None else "N/A"
                se_val = f"{row['se']:.4f}" if row['se'] is not None else "N/A"
                p_val = f"{row['p']:.2e}" if row['p'] is not None else "N/A"
                rg_rows += f"<tr><td>{row['p1']}</td><td>{row['p2']}</td><td>{rg_val}</td><td>{se_val}</td><td>{p_val}</td></tr>"
            rg_table_html = f"""
            <table class="data-table">
                <thead>
                    <tr><th>Trait 1</th><th>Trait 2</th><th>rg</th><th>SE</th><th>P-value</th></tr>
                </thead>
                <tbody>{rg_rows}</tbody>
            </table>
            """
        elif rg_failure_html:
            rg_table_html = rg_failure_html

    h2_table_html = "<p><em>LDSC heritability analysis was not run or returned no results.</em></p>"
    if ldsc_h2_csv and os.path.exists(ldsc_h2_csv):
        df_h2 = pl.read_csv(ldsc_h2_csv)
        if len(df_h2) > 0:
            h2_rows = ""
            h2_log_dirs = [
                os.path.join(out_dir, "ldsc_h2"),
                os.path.join(out_dir, "04_ldsc_analysis", "ldsc_h2"),
                os.path.join(os.path.dirname(ldsc_h2_csv), "ldsc_h2"),
            ]
            for row in df_h2.iter_rows(named=True):
                trait = row.get("trait")
                log_values = {}
                for h2_log_dir in h2_log_dirs:
                    parsed = _parse_ldsc_h2_log(os.path.join(h2_log_dir, f"{trait}.log"))
                    if parsed:
                        log_values = parsed
                        break

                n_snps = row.get("n_snps") if row.get("n_snps") is not None else log_values.get("n_snps")
                h2_obs = row.get("h2_obs") if row.get("h2_obs") is not None else log_values.get("h2_obs")
                h2_obs_se = row.get("h2_obs_se") if row.get("h2_obs_se") is not None else log_values.get("h2_obs_se")
                h2_liab = row.get("h2_liab") if row.get("h2_liab") is not None else log_values.get("h2_liab")
                h2_liab_se = row.get("h2_liab_se") if row.get("h2_liab_se") is not None else log_values.get("h2_liab_se")
                intercept = row.get("intercept") if row.get("intercept") is not None else log_values.get("intercept")
                lambda_gc = row.get("lambda_gc") if row.get("lambda_gc") is not None else row.get("lambda_GC") if "lambda_GC" in row else log_values.get("lambda_gc")
                mean_chi2 = row.get("mean_chi2") if row.get("mean_chi2") is not None else log_values.get("mean_chi2")

                h2_rows += (
                    f"<tr><td>{_esc(trait)}</td><td>{_fmt_num(n_snps, 0)}</td>"
                    f"<td>{_fmt_num(h2_obs)}</td><td>{_fmt_num(h2_obs_se)}</td>"
                    f"<td>{_fmt_num(h2_liab)}</td><td>{_fmt_num(h2_liab_se)}</td>"
                    f"<td>{_fmt_num(intercept)}</td><td>{_fmt_num(lambda_gc)}</td>"
                    f"<td>{_fmt_num(mean_chi2)}</td></tr>"
                )
            h2_table_html = f"""
            <table class="data-table">
                <thead>
                    <tr><th>Trait</th><th>SNPs</th><th>Observed h²</th><th>Obs SE</th><th>Liability h²</th><th>Liab SE</th><th>Intercept</th><th>λ GC</th><th>Mean χ²</th></tr>
                </thead>
                <tbody>{h2_rows}</tbody>
            </table>
            """

    # ── Input Traits & Manifest Summary Table ─────────────────────────────────
    trait_overview_rows = ""
    target_trait = manifest_df["sample_id"][0]
    covariate_traits = manifest_df["sample_id"].to_list()[1:]
    
    for idx, row in enumerate(manifest_df.iter_rows(named=True)):
        trait_id  = row["sample_id"]
        vcf_path  = row["file_path"]
        is_target = (idx == 0)
        role_badge = '<span style="background:#1a73e8; color:white; padding:3px 8px; border-radius:4px; font-weight:600; font-size:12px;">Target Trait</span>' if is_target else '<span style="background:#5f6368; color:white; padding:3px 8px; border-radius:4px; font-weight:600; font-size:12px;">Covariate Trait</span>'
        
        samp_prev = row.get("sample_prevalence")
        pop_prev  = row.get("population_prevalence")
        sp_str = f"{float(samp_prev):.4f}" if samp_prev is not None and str(samp_prev).strip() != "" and str(samp_prev).upper() != "NA" else "N/A (Quantitative)"
        pp_str = f"{float(pop_prev):.4f}" if pop_prev is not None and str(pop_prev).strip() != "" and str(pop_prev).upper() != "NA" else "N/A (Quantitative)"
        
        ma_file = ma_paths.get(trait_id, f"{trait_id}.ma")
        
        trait_overview_rows += f"""
        <tr>
            <td><strong>{trait_id}</strong></td>
            <td>{role_badge}</td>
            <td><code class="path-scroll" title="{vcf_path}">{vcf_path}</code></td>
            <td>{sp_str}</td>
            <td>{pp_str}</td>
        </tr>
        """

    # ── Complete Per-Step Parameters & Configuration Matrix ────────────────────
    step_params = [
        ("Step 1: Manifest & VCF Conversion", "Target Trait", target_trait),
        ("Step 1: Manifest & VCF Conversion", "Covariate Trait(s)", ", ".join(covariate_traits)),
        ("Step 1: Manifest & VCF Conversion", "bcftools Binary Path", config_params.get("bcftools", bcftools_bin)),
        ("Step 1: Manifest & VCF Conversion", "Variant ID Format", config_params.get("id_format", "RSID (Auto-detected)")),
        ("Step 1: Manifest & VCF Conversion", "Min Variant Overlap Cutoff", "1% (0.01)"),
        ("Step 1: Manifest & VCF Conversion", "VCF FORMAT Fields Extracted", "ES (Beta), SE (StdErr), LP (-log10 P), NEF (N_eff), SI (INFO)"),
        
        ("Step 2: GCTA mtCOJO Execution", "gcta64 Binary Path", config_params.get("gcta64", "gcta64")),
        ("Step 2: GCTA mtCOJO Execution", "PLINK Reference Panel (bfile)", config_params.get("bfile", "N/A")),
        ("Step 2: GCTA mtCOJO Execution", "LD Reference Directory (ref_ld_chr)", config_params.get("ref_ld_chr", "N/A")),
        ("Step 2: GCTA mtCOJO Execution", "LD Weight Directory (w_ld_chr)", config_params.get("w_ld_chr", config_params.get("ref_ld_chr", "N/A"))),
        ("Step 2: GCTA mtCOJO Execution", "Index SNP P-value Threshold (gwas_thresh)", config_params.get("gwas_thresh", "1e-05 (Default)")),
        ("Step 2: GCTA mtCOJO Execution", "LD Clumping r² Cutoff (clump_r2)", config_params.get("clump_r2", "Default / Auto")),
        ("Step 2: GCTA mtCOJO Execution", "HEIDI Outlier Cutoff (heidi_thresh)", config_params.get("heidi_thresh", "0.01 (Default)")),
        ("Step 2: GCTA mtCOJO Execution", "Min Index SNPs for GSMR (gsmr_snp_min)", config_params.get("gsmr_snp_min", "10 (Default)")),
        ("Step 2: GCTA mtCOJO Execution", "Max Allele Frequency Diff (diff_freq)", config_params.get("diff_freq", "Default / None")),
        
        ("Step 3: PostGWAS Harmonisation", "PostGWAS Enabled", "YES" if config_params.get("run_postgwas") else "NO (Skipped)"),
        ("Step 3: PostGWAS Harmonisation", "Docker Image Tag", config_params.get("docker_image", "cbiit/postgwas:v1")),
        ("Step 3: PostGWAS Harmonisation", "GRCh37→GRCh38 Liftover", config_params.get("liftover", "No")),
        ("Step 3: PostGWAS Harmonisation", "Docker Threads", config_params.get("nthreads", "23")),
        ("Step 3: PostGWAS Harmonisation", "Docker Max Memory", config_params.get("max_mem", "50G")),
        ("Step 3: PostGWAS Harmonisation", "Config YAML (defaults)", config_params.get("defaults", "harmonisation.yaml")),
        
        ("Step 4: LDSC Heritability & rg", "LDSC Enabled", "YES" if config_params.get("run_ldsc") else "NO (Skipped; existing results detected)" if ldsc_artifacts_available else "NO (Skipped)"),
        ("Step 4: LDSC Heritability & rg", "LDSC Results CSV", ldsc_results_csv if ldsc_results_csv and os.path.exists(ldsc_results_csv) else "Not found"),
        ("Step 4: LDSC Heritability & rg", "LDSC h² CSV", ldsc_h2_csv if ldsc_h2_csv and os.path.exists(ldsc_h2_csv) else "Not found"),
        ("Step 4: LDSC Heritability & rg", "LDSC Directory (CBIIT/ldsc)", config_params.get("ldsc_dir", "Auto-detected")),
        ("Step 4: LDSC Heritability & rg", "HapMap3 SNP List (ldsc_snp_list)", config_params.get("ldsc_snp_list", "w_hm3.snplist")),
        ("Step 4: LDSC Heritability & rg", "Parallel Processes (ldsc_n_parallel)", config_params.get("ldsc_n_parallel", "4")),
        ("Step 4: LDSC Heritability & rg", "Batch Size per Call (ldsc_batch_size)", config_params.get("ldsc_batch_size", "10")),
        ("Step 4: LDSC Heritability & rg", "Center Z-scores (ldsc_center_z)", config_params.get("ldsc_center_z", "False")),
        
        ("Step 5: Visualizations & Plotting", "R Plotting Engine", "rMVP::MVP.Report / CMplot"),
        ("Step 5: Visualizations & Plotting", "Overlaid Manhattan Point Size", "s = 0.3, alpha = 0.32"),
        ("Step 5: Visualizations & Plotting", "Stacked Manhattan Point Size", "s = 0.3 for chromosome bands; s = 0.7 for single-chromosome tracks"),
        ("Step 5: Visualizations & Plotting", "Manhattan Background Thinning", "Retain all P < 0.1; retain 5% of variants with P between 0.1 and 1"),
        ("Step 5: Visualizations & Plotting", "Plot P-value Range", "Clamp plotting p-values to [1e-300, 1]"),
        ("Step 5: Visualizations & Plotting", "Color Palette", "Muted publication GWAS palette: #0072B2, #D55E00, #009E73, #CC79A7, #E69F00, #56B4E9, #000000"),
        ("Step 5: Visualizations & Plotting", "P-value Threshold Lines", "1e-5 (Suggestive, blue dotted) & 1e-8 (Genome-wide, red dashed)"),
        ("Step 5: Visualizations & Plotting", "Q-Q Confidence Envelope", "conf.int = TRUE (Gray shaded 95% CI band)"),
    ]
    
    step_params_matrix = ""
    for step_cat, param_name, param_val in step_params:
        val_str = str(param_val)
        val_cell = f'<code class="path-scroll" title="{val_str}">{val_str}</code>' if len(val_str) > 35 else f'<code>{val_str}</code>'
        step_params_matrix += f"""
        <tr>
            <td style="color:#1b4965; font-weight:600;">{step_cat}</td>
            <td style="font-weight:600;">{param_name}</td>
            <td>{val_cell}</td>
        </tr>
        """

    report_file = os.path.join(out_dir, f"{out_name}_report.html")
    file_inventory_html = _build_file_inventory(out_dir, report_file)

    plot_notes_html = ""
    if plots.get("notes"):
        plot_notes_html = "<ul>" + "".join(f"<li>{_esc(note)}</li>" for note in plots["notes"]) + "</ul>"

    plot_file_rows = ""
    for label, path in plots.get("files", {}).items():
        plot_file_rows += f"""
        <tr>
            <td>{_esc(label.replace('_', ' ').title())}</td>
            <td><code class="path-scroll" title="{_esc(path)}">{_esc(os.path.relpath(path, out_dir))}</code></td>
        </tr>
        """
    plot_files_html = ""
    if plot_file_rows:
        plot_files_html = f"""
        <details>
            <summary>Displayed Plot Files</summary>
            <table class="data-table">
                <thead><tr><th>Plot</th><th>File Embedded in HTML</th></tr></thead>
                <tbody>{plot_file_rows}</tbody>
            </table>
        </details>
        """

    heatmap_card_html = ""
    if plots.get("heatmap"):
        heatmap_card_html = f"""
        <div class="img-card">
            <h3>LDSC Pairwise Genetic Correlation (rg) Heatmap</h3>
            <img src="{plots["heatmap"]}" alt="LDSC Genetic Correlation Heatmap"/>
        </div>
        """
    elif config_params.get("run_ldsc") or ldsc_artifacts_available:
        heatmap_message = (
            "LDSC rg logs were detected, but no pairwise rg estimates were compiled; "
            "the heatmap requires at least one successful rg result. See the LDSC Pairwise Genetic Correlations section below."
            if rg_failure_html else
            "LDSC results were detected, but no LDSC rg heatmap image was generated."
        )
        heatmap_card_html = f"""
        <div class="img-card">
            <h3>LDSC Pairwise Genetic Correlation (rg) Heatmap</h3>
            <p><em>{_esc(heatmap_message)}</em></p>
        </div>
        """

    ldsc_sections_html = ""
    if config_params.get("run_ldsc") or ldsc_artifacts_available:
        ldsc_sections_html = f"""
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">LDSC Heritability (h²)</summary>
                <div style="margin-top:12px;">
                    {h2_table_html}
                </div>
            </details>
        </div>

        <div class="card">
            <details open class="card-details">
                <summary class="section-title">LDSC Pairwise Genetic Correlations (rg)</summary>
                <div style="margin-top:12px;">
                    {heatmap_card_html}
                    {rg_table_html}
                </div>
            </details>
        </div>
        """
    else:
        ldsc_sections_html = """
        <div class="card">
            <details class="card-details">
                <summary class="section-title">LDSC Analysis [Skipped]</summary>
                <p><em>LDSC was not requested for this run, so heritability, rg tables, and rg heatmap are intentionally absent.</em></p>
            </details>
        </div>
        """

    html_content = f"""<!DOCTYPE html>
<html lang="en">
<head>
    <meta charset="UTF-8">
    <meta name="viewport" content="width=device-width, initial-scale=1.0">
    <title>GCTA mtCOJO & LDSC Analysis Dashboard</title>
    <style>
        :root {{
            --primary: #1a73e8;
            --bg-dark: #121212;
            --card-bg: #ffffff;
            --border-color: #e0e0e0;
            --text-main: #202124;
            --text-sub: #5f6368;
        }}
        body {{
            font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", Roboto, Helvetica, Arial, sans-serif;
            background-color: #f8f9fa;
            color: var(--text-main);
            margin: 0;
            padding: 24px;
            line-height: 1.5;
        }}
        .container {{
            max-width: 2200px;
            margin: 0 auto;
        }}
        .header {{
            background: linear-gradient(135deg, #1a73e8 0%, #0d47a1 100%);
            color: white;
            padding: 32px;
            border-radius: 12px;
            box-shadow: 0 4px 12px rgba(0,0,0,0.1);
            margin-bottom: 24px;
        }}
        .header h1 {{ margin: 0 0 8px 0; font-size: 28px; font-weight: 700; }}
        .header p {{ margin: 0; opacity: 0.9; font-size: 15px; }}
        
        .card {{
            background: var(--card-bg);
            border: 1px solid var(--border-color);
            border-radius: 12px;
            padding: 24px;
            margin-bottom: 24px;
            box-shadow: 0 2px 6px rgba(0,0,0,0.04);
        }}
        .card h2 {{
            margin-top: 0;
            font-size: 20px;
            border-bottom: 2px solid #f1f3f4;
            padding-bottom: 12px;
            color: #1a73e8;
        }}
        
        .card-details summary.section-title {{
            font-size: 20px;
            font-weight: 700;
            color: #1a73e8;
            cursor: pointer;
            padding-bottom: 8px;
            outline: none;
            user-select: none;
        }}
        
        /* Flowchart Styling */
        .flowchart-container {{
            display: flex;
            align-items: stretch;
            justify-content: space-between;
            gap: 12px;
            overflow-x: auto;
            padding: 12px 0;
        }}
        .flow-step {{
            flex: 1;
            min-width: 210px;
            background: linear-gradient(135deg, #ffffff 0%, #f8f9fa 100%);
            border: 1px solid #dadce0;
            border-top: 4px solid #1a73e8;
            border-radius: 8px;
            padding: 16px;
            box-shadow: 0 2px 4px rgba(0,0,0,0.04);
        }}
        .step-num {{
            font-size: 11px;
            font-weight: 700;
            text-transform: uppercase;
            color: #1a73e8;
            letter-spacing: 0.5px;
            margin-bottom: 4px;
        }}
        .step-title {{
            font-size: 14px;
            font-weight: 700;
            color: #202124;
            margin-bottom: 8px;
        }}
        .step-desc {{
            font-size: 12px;
            color: #5f6368;
            line-height: 1.4;
        }}
        .flow-arrow {{
            display: flex;
            align-items: center;
            font-size: 22px;
            color: #1a73e8;
            font-weight: bold;
            user-select: none;
        }}
        
        table.data-table {{
            width: 100%;
            border-collapse: collapse;
            margin: 16px 0;
            font-size: 14px;
        }}
        table.data-table th, table.data-table td {{
            padding: 10px 14px;
            text-align: left;
            border-bottom: 1px solid #e0e0e0;
        }}
        table.data-table th {{
            background-color: #f1f3f4;
            font-weight: 600;
            color: #3c4043;
            cursor: pointer;
            position: relative;
            user-select: none;
        }}
        table.data-table th.sort-asc::after {{ content: " ▲"; color: #1a73e8; font-size: 11px; }}
        table.data-table th.sort-desc::after {{ content: " ▼"; color: #1a73e8; font-size: 11px; }}
        table.data-table tr:hover {{ background-color: #f8f9fa; }}
        table.compact-table {{
            margin: 8px 0;
            font-size: 13px;
        }}
        table.compact-table th, table.compact-table td {{
            padding: 7px 10px;
        }}
        .overlap-table-grid {{
            display: grid;
            grid-template-columns: repeat(auto-fit, minmax(280px, 1fr));
            gap: 16px;
            margin-top: 10px;
        }}
        .overlap-summary {{
            background: #ffffff;
            margin-top: 12px;
        }}
        .folder-inventory {{
            background: #ffffff;
            border-color: #e0e0e0;
            margin-top: 12px;
        }}
        .folder-inventory > summary {{
            color: #202124;
        }}
        .table-shell {{
            margin: 16px 0;
            overflow-x: auto;
        }}
        .table-shell table.data-table {{
            margin: 0;
        }}
        .table-controls {{
            display: flex;
            flex-wrap: wrap;
            align-items: center;
            justify-content: space-between;
            gap: 10px;
            margin: 16px 0 8px 0;
            color: #5f6368;
            font-size: 13px;
        }}
        .table-controls label {{
            display: inline-flex;
            align-items: center;
            gap: 6px;
        }}
        .table-controls select,
        .table-controls button {{
            border: 1px solid #dadce0;
            background: #ffffff;
            color: #202124;
            border-radius: 6px;
            padding: 5px 8px;
            font-size: 13px;
        }}
        .table-controls button {{
            cursor: pointer;
            font-weight: 600;
        }}
        .table-controls button:disabled {{
            color: #9aa0a6;
            cursor: default;
            background: #f1f3f4;
        }}
        .table-pager {{
            display: inline-flex;
            align-items: center;
            gap: 8px;
        }}
        .alert-box {{
            background: #fff4e5;
            border: 1px solid #f4b860;
            border-left: 4px solid #f77f00;
            border-radius: 6px;
            color: #5f370e;
            padding: 12px 14px;
            margin: 12px 0;
            font-size: 14px;
        }}
        
        details {{
            background-color: #f8f9fa;
            border: 1px solid #dadce0;
            border-radius: 8px;
            padding: 12px 16px;
            margin-bottom: 12px;
        }}
        summary {{
            font-weight: 600;
            cursor: pointer;
            color: #1a73e8;
        }}
        code, pre {{
            background-color: #272822;
            color: #f8f8f2;
            padding: 12px;
            border-radius: 6px;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
            font-size: 13px;
            white-space: pre-wrap;
            word-break: break-all;
            overflow-x: auto;
            display: block;
        }}
        code.path-scroll {{
            display: inline-block;
            max-width: 260px;
            overflow-x: auto;
            white-space: nowrap;
            vertical-align: middle;
            background-color: #272822;
            color: #a6e22e;
            padding: 3px 8px;
            border-radius: 4px;
            font-size: 12px;
            font-family: "SFMono-Regular", Consolas, "Liberation Mono", Menlo, monospace;
        }}
        code.path-scroll::-webkit-scrollbar {{
            height: 4px;
        }}
        code.path-scroll::-webkit-scrollbar-thumb {{
            background-color: #666;
            border-radius: 2px;
        }}
        
        .img-grid {{
            display: grid;
            grid-template-columns: 1fr;
            gap: 24px;
        }}
        .img-card {{
            text-align: center;
            background: #fafafa;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            padding: 16px;
        }}
        .img-card img {{
            max-width: 100%;
            max-height: 480px;
            width: auto;
            height: auto;
            object-fit: contain;
            border-radius: 4px;
            margin: 0 auto;
            display: block;
        }}
        .img-card-wide img {{
            max-width: 100%;
            max-height: 600px;
            width: auto;
            height: auto;
            object-fit: contain;
            border-radius: 4px;
            margin: 0 auto;
            display: block;
        }}
        .img-card-fullwide {{
            text-align: center;
            background: #ffffff;
            border: 1px solid #e0e0e0;
            border-radius: 8px;
            padding: 20px;
            overflow-x: auto;
            width: 100%;
            box-sizing: border-box;
        }}
        .img-card-fullwide img {{
            max-width: 100% !important;
            height: auto !important;
            max-height: 700px !important;
            object-fit: contain !important;
            border-radius: 4px;
            margin: 0 auto;
            display: block;
        }}
    </style>
</head>
<body>
    <div class="container">
        <div class="header">
            <h1>GCTA mtCOJO & LDSC Pipeline Dashboard</h1>
            <p>Target Trait Output: <code>{out_prefix}.mtcojo.cma</code></p>
        </div>

        <!-- 1a. Input Traits & Manifest Overview -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Input Datasets &amp; Trait Roles</summary>
                <table class="data-table" style="margin-top:12px;">
                    <thead>
                        <tr>
                            <th>Sample ID</th>
                            <th>Pipeline Role</th>
                            <th>Input VCF File Path</th>
                            <th>Sample Prevalence (K)</th>
                            <th>Population Prevalence (P)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {trait_overview_rows}
                    </tbody>
                </table>
            </details>
        </div>

        <!-- 1b. Pipeline Analytical Workflow & Execution Flowchart -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Pipeline Analytical Workflow &amp; Execution Flowchart</summary>
                <div class="flowchart-container" style="margin-top: 16px;">
                    <div class="flow-step">
                        <div class="step-num">Step 1</div>
                        <div class="step-title">Manifest &amp; VCF Extraction</div>
                        <div class="step-desc">
                            Reads target &amp; covariate GWAS summary stats in VCF format using <code>bcftools</code>. Extracts effect sizes (&beta;), StdErrors (SE), -log<sub>10</sub>(P), and aligns variant IDs.
                        </div>
                    </div>
                    <div class="flow-arrow">➔</div>
                    <div class="flow-step">
                        <div class="step-num">Step 2</div>
                        <div class="step-title">GCTA .ma Format Conversion</div>
                        <div class="step-desc">
                            Converts extracted VCF statistics into GCTA <code>.ma</code> format (SNP, A1, A2, freq, b, se, p, N) for all target and covariate traits in a single streaming pass.
                        </div>
                    </div>
                    <div class="flow-arrow">➔</div>
                    <div class="flow-step">
                        <div class="step-num">Step 3</div>
                        <div class="step-title">GCTA mtCOJO Conditioning</div>
                        <div class="step-desc">
                            Executes <code>gcta64 --mtcojo</code> with PLINK reference panel &amp; LD scores. Performs multi-trait conditioning, joint effect calculation, &amp; HEIDI outlier filtering.
                        </div>
                    </div>
                    <div class="flow-arrow">➔</div>
                    <div class="flow-step">
                        <div class="step-num">Step 4</div>
                        <div class="step-title">PostGWAS &amp; LDSC rg Analysis</div>
                        <div class="step-desc">
                            Executes optional Docker harmonisation &amp; runs LD Score Regression (<code>ldsc.py</code>) to compute SNP heritability (h²) and pairwise genetic correlations (rg).
                        </div>
                    </div>
                    <div class="flow-arrow">➔</div>
                    <div class="flow-step">
                        <div class="step-num">Step 5</div>
                        <div class="step-title">Multi-Trait Visualization Suite</div>
                        <div class="step-desc">
                            Invokes R (<code>rMVP</code>/<code>CMplot</code>) &amp; Python to render Stacked &amp; Overlaid Manhattan plots, Q-Q plots, LDSC Heatmaps, Venn Diagrams, and UpSet plots.
                        </div>
                    </div>
                </div>
            </details>
        </div>

        <!-- 1c. Output Files -->
        <div class="card">
            <details class="card-details">
                <summary class="section-title">Generated Output Files</summary>
                <p style="color:#5f6368; font-size:14px; margin: 12px 0 16px 0;">
                    Files detected inside the root output directory:
                    <code style="background:#e8f0fe; color:#1a73e8; padding:3px 6px; border-radius:4px; font-weight:600;">{out_dir}</code>
                </p>
                {file_inventory_html}
            </details>
        </div>

        <!-- 1c. Executed CLI & Parameters Matrix -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Pipeline Parameters &amp; Command Execution</summary>
                <div style="margin-top:12px;">
                    <details>
                        <summary>Click to view exact command-line string executed</summary>
                        <pre>{cli_command}</pre>
                    </details>
                    <details>
                        <summary>Click to view complete per-step parameters &amp; configuration matrix (including defaults)</summary>
                        <table class="data-table">
                            <thead>
                                <tr>
                                    <th>Pipeline Stage</th>
                                    <th>Parameter Name</th>
                                    <th>Value / Configuration</th>
                                </tr>
                            </thead>
                            <tbody>
                                {step_params_matrix}
                            </tbody>
                        </table>
                    </details>
                </div>
            </details>
        </div>

        <!-- 2. Variant Counts & Significance Percentages Table -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Variant Counts &amp; Threshold Significance Summary</summary>
                <table class="data-table" style="margin-top:12px;">
                    <thead>
                        <tr>
                            <th>Stage / Dataset</th>
                            <th>Total Variants</th>
                            <th>P &lt; 0.05 (Count &amp; %)</th>
                            <th>P &lt; 10⁻⁵ (Count &amp; %)</th>
                            <th>P &lt; 10⁻⁸ (Count &amp; %)</th>
                        </tr>
                    </thead>
                    <tbody>
                        {stats_rows}
                    </tbody>
                </table>
            </details>
        </div>

        <!-- 2b. Top Result Tables -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Top Association Results</summary>
                <h3>Converted .ma Files</h3>
                {top_results_html}
                <h3>mtCOJO Conditioned Output</h3>
                {cma_top_html}
            </details>
        </div>

        <!-- 2c. Conditional Significance Shift Summary -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Conditional Significance Shift Summary</summary>
                {conditional_shift_html}
            </details>
        </div>

        <!-- 3. Multi-Level Variant Overlap & Intersection Suite -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Multi-Level Variant Overlap &amp; Intersection Visualizations</summary>
                <p style="color:#5f6368; font-size:14px; margin: 12px 0 16px 0;">
                    Venn Diagrams and UpSet Plots displaying variant ID overlaps across 4 significance levels: 
                    <strong>Total Variants</strong> (shown below), and expandable sections for <strong>P &lt; 0.05</strong>, <strong>P &lt; 10⁻⁵</strong>, and <strong>P &lt; 10⁻⁸</strong>.
                </p>

                <!-- Level 1: Total Variants (Open by default) -->
                <details open style="background: white; border: 1px solid #dadce0; margin-bottom: 16px;">
                    <summary style="font-size: 15px; font-weight: 700; color: #1b4965;">1. Total Variants Overlap (All Variants)</summary>
                    <div class="img-grid" style="grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); margin-top: 12px;">
                        {f'''<div class="img-card">
                            <h4>Venn Diagram: Total Variants</h4>
                            <img src="{overlaps['total']['venn']}" alt="Venn Diagram: Total Variants"/>
                        </div>''' if overlaps['total']['venn'] else '<div class="img-card"><h4>Venn Diagram: Total Variants</h4><p><em>Venn image was not generated; overlap counts are shown below.</em></p></div>'}
                        {f'''<div class="img-card">
                            <h4>UpSet Plot: Total Variants</h4>
                            <img src="{overlaps['total']['upset']}" alt="UpSet Plot: Total Variants"/>
                        </div>''' if overlaps['total']['upset'] else ''}
                    </div>
                    {overlaps['total'].get('summary', '')}
                </details>

                <!-- Level 2: P < 0.05 (Expandable) -->
                <details style="background: white; border: 1px solid #dadce0; margin-bottom: 16px;">
                    <summary style="font-size: 15px; font-weight: 700; color: #1b4965;">2. Nominal Significance (P &lt; 0.05) Overlap</summary>
                    <div class="img-grid" style="grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); margin-top: 12px;">
                        {f'''<div class="img-card">
                            <h4>Venn Diagram: P &lt; 0.05</h4>
                            <img src="{overlaps['p005']['venn']}" alt="Venn Diagram: P < 0.05"/>
                        </div>''' if overlaps['p005']['venn'] else '<p><em>No overlapping variants at P &lt; 0.05 threshold.</em></p>'}
                        {f'''<div class="img-card">
                            <h4>UpSet Plot: P &lt; 0.05</h4>
                            <img src="{overlaps['p005']['upset']}" alt="UpSet Plot: P < 0.05"/>
                        </div>''' if overlaps['p005']['upset'] else ''}
                    </div>
                    {overlaps['p005'].get('summary', '')}
                </details>

                <!-- Level 3: P < 1e-5 (Expandable) -->
                <details style="background: white; border: 1px solid #dadce0; margin-bottom: 16px;">
                    <summary style="font-size: 15px; font-weight: 700; color: #1b4965;">3. Suggestive Significance (P &lt; 10⁻⁵) Overlap</summary>
                    <div class="img-grid" style="grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); margin-top: 12px;">
                        {f'''<div class="img-card">
                            <h4>Venn Diagram: P &lt; 10⁻⁵</h4>
                            <img src="{overlaps['p1e5']['venn']}" alt="Venn Diagram: P < 1e-5"/>
                        </div>''' if overlaps['p1e5']['venn'] else '<p><em>No overlapping variants at P &lt; 10⁻⁵ threshold.</em></p>'}
                        {f'''<div class="img-card">
                            <h4>UpSet Plot: P &lt; 10⁻⁵</h4>
                            <img src="{overlaps['p1e5']['upset']}" alt="UpSet Plot: P < 1e-5"/>
                        </div>''' if overlaps['p1e5']['upset'] else ''}
                    </div>
                    {overlaps['p1e5'].get('summary', '')}
                </details>

                <!-- Level 4: P < 1e-8 (Expandable) -->
                <details style="background: white; border: 1px solid #dadce0;">
                    <summary style="font-size: 15px; font-weight: 700; color: #1b4965;">4. Genome-Wide Significance (P &lt; 10⁻⁸) Overlap</summary>
                    <div class="img-grid" style="grid-template-columns: repeat(auto-fit, minmax(360px, 1fr)); margin-top: 12px;">
                        {f'''<div class="img-card">
                            <h4>Venn Diagram: P &lt; 10⁻⁸</h4>
                            <img src="{overlaps['p1e8']['venn']}" alt="Venn Diagram: P < 1e-8"/>
                        </div>''' if overlaps['p1e8']['venn'] else '<p><em>No overlapping variants at P &lt; 10⁻⁸ threshold.</em></p>'}
                        {f'''<div class="img-card">
                            <h4>UpSet Plot: P &lt; 10⁻⁸</h4>
                            <img src="{overlaps['p1e8']['upset']}" alt="UpSet Plot: P < 1e-8"/>
                        </div>''' if overlaps['p1e8']['upset'] else ''}
                    </div>
                    {overlaps['p1e8'].get('summary', '')}
                </details>
            </details>
        </div>

        <!-- 4. Visualizations Suite -->
        <div class="card">
            <details open class="card-details">
                <summary class="section-title">Multi-Trait Visualizations (rMVP Graphics Suite)</summary>
                {plot_notes_html}
                {plot_files_html}
                <div class="img-grid" style="margin-top:16px;">
                    <div class="img-card-fullwide">
                        <h3>Multi-Trait Stacked Manhattan Plot (Rectangular Tracks)</h3>
                        {f'<img src="{plots["manhattan"]}" alt="Multi-Trait Stacked Manhattan Plot"/>' if plots["manhattan"] else '<p><em>Stacked Manhattan plot image not available.</em></p>'}
                    </div>
                    {f'''<div class="img-card img-card-wide">
                        <h3>Multi-Trait Overlaid Manhattan Plot (Combined Single Track)</h3>
                        <img src="{plots["manhattan_overlaid"]}" alt="Multi-Trait Overlaid Manhattan Plot"/>
                    </div>''' if plots["manhattan_overlaid"] else ''}
                    <div class="img-card">
                        <h3>Multi-Trait 3-Panel Q-Q Plot (Separate Panel Tracks)</h3>
                        {f'<img src="{plots["qq"]}" alt="Multi-Trait 3-Panel Q-Q Plot"/>' if plots["qq"] else '<p><em>Q-Q plot image not available.</em></p>'}
                    </div>
                    {f'''<div class="img-card">
                        <h3>Multi-Trait Combined Q-Q Plot (All Traits Overlaid Single Track)</h3>
                        <img src="{plots["qq_combined"]}" alt="Multi-Trait Combined Q-Q Plot"/>
                    </div>''' if plots["qq_combined"] else ''}
                </div>
            </details>
        </div>

        {ldsc_sections_html}
    </div>
    <script>
        (function () {{
            function parseCell(value) {{
                var text = (value || "").trim();
                if (!text) return "";
                var numeric = Number(text.replace(/,/g, ""));
                return Number.isFinite(numeric) ? numeric : text.toLowerCase();
            }}

            function initTable(table, index) {{
                if (!table.tBodies.length || table.dataset.interactiveReady === "1") return;
                var tbody = table.tBodies[0];
                var originalRows = Array.from(tbody.rows);
                if (originalRows.length <= 10) return;
                table.dataset.interactiveReady = "1";

                var state = {{
                    rows: originalRows.slice(),
                    page: 1,
                    pageSize: 10,
                    sortColumn: null,
                    sortDirection: "asc"
                }};

                var controls = document.createElement("div");
                controls.className = "table-controls";
                controls.innerHTML =
                    '<label>Rows <select aria-label="Rows per page">' +
                    '<option value="10">10</option><option value="20">20</option>' +
                    '<option value="50">50</option><option value="all">All</option>' +
                    '</select></label>' +
                    '<span class="table-info"></span>' +
                    '<span class="table-pager">' +
                    '<button type="button" data-action="prev">Previous</button>' +
                    '<span class="table-page"></span>' +
                    '<button type="button" data-action="next">Next</button>' +
                    '</span>';

                var shell = document.createElement("div");
                shell.className = "table-shell";
                table.parentNode.insertBefore(controls, table);
                table.parentNode.insertBefore(shell, table);
                shell.appendChild(table);

                var select = controls.querySelector("select");
                var info = controls.querySelector(".table-info");
                var pageText = controls.querySelector(".table-page");
                var prev = controls.querySelector('[data-action="prev"]');
                var next = controls.querySelector('[data-action="next"]');

                function pageCount() {{
                    return state.pageSize === "all" ? 1 : Math.max(1, Math.ceil(state.rows.length / state.pageSize));
                }}

                function render() {{
                    var pages = pageCount();
                    state.page = Math.min(Math.max(1, state.page), pages);
                    var start = state.pageSize === "all" ? 0 : (state.page - 1) * state.pageSize;
                    var end = state.pageSize === "all" ? state.rows.length : start + state.pageSize;
                    var visible = state.rows.slice(start, end);

                    tbody.replaceChildren.apply(tbody, visible);
                    info.textContent = "Showing " + (visible.length ? start + 1 : 0) + "-" + Math.min(end, state.rows.length) + " of " + state.rows.length;
                    pageText.textContent = "Page " + state.page + " / " + pages;
                    prev.disabled = state.page <= 1 || state.pageSize === "all";
                    next.disabled = state.page >= pages || state.pageSize === "all";
                }}

                Array.from(table.tHead ? table.tHead.rows[0].cells : []).forEach(function (th, colIndex) {{
                    th.title = "Click to sort";
                    th.addEventListener("click", function () {{
                        var direction = state.sortColumn === colIndex && state.sortDirection === "asc" ? "desc" : "asc";
                        state.sortColumn = colIndex;
                        state.sortDirection = direction;
                        Array.from(th.parentNode.cells).forEach(function (cell) {{
                            cell.classList.remove("sort-asc", "sort-desc");
                        }});
                        th.classList.add(direction === "asc" ? "sort-asc" : "sort-desc");
                        state.rows.sort(function (a, b) {{
                            var av = parseCell(a.cells[colIndex] ? a.cells[colIndex].textContent : "");
                            var bv = parseCell(b.cells[colIndex] ? b.cells[colIndex].textContent : "");
                            if (av < bv) return direction === "asc" ? -1 : 1;
                            if (av > bv) return direction === "asc" ? 1 : -1;
                            return 0;
                        }});
                        state.page = 1;
                        render();
                    }});
                }});

                select.addEventListener("change", function () {{
                    state.pageSize = select.value === "all" ? "all" : Number(select.value);
                    state.page = 1;
                    render();
                }});
                prev.addEventListener("click", function () {{ state.page -= 1; render(); }});
                next.addEventListener("click", function () {{ state.page += 1; render(); }});
                render();
            }}

            document.addEventListener("DOMContentLoaded", function () {{
                document.querySelectorAll("table.data-table").forEach(initTable);
            }});
        }})();
    </script>
</body>
</html>
"""

    with open(report_file, "w", encoding="utf-8") as f:
        f.write(html_content)

    log_pass(log, f"  HTML Pipeline Report successfully written → {report_file}")
    return report_file
