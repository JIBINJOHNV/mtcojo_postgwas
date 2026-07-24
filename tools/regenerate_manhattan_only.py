#!/usr/bin/env python3
"""Regenerate lightweight Manhattan PNGs from merged_gwas_summary.tsv."""

from __future__ import annotations

import argparse
import base64
import math
import re
from pathlib import Path

import polars as pl


def _clamped_p_expr(col: str) -> pl.Expr:
    p = pl.col(col).cast(pl.Float64, strict=False)
    return (
        pl.when(p.is_null() | (p <= 0))
        .then(None)
        .when(p < 1e-300)
        .then(1e-300)
        .when(p > 1.0)
        .then(1.0)
        .otherwise(p)
        .alias(col)
    )


def _replace_html_image(report_html: Path, alt_text: str, image_path: Path) -> None:
    if not report_html.exists():
        return
    b64 = base64.b64encode(image_path.read_bytes()).decode("ascii")
    src = f"data:image/png;base64,{b64}"
    html = report_html.read_text(errors="ignore")
    pattern = re.compile(rf'(<img src=")data:image/[^"]+("[^>]*alt="{re.escape(alt_text)}")')
    html_new, count = pattern.subn(rf"\1{src}\2", html, count=1)
    if count:
        report_html.write_text(html_new)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--merged-tsv", required=True)
    parser.add_argument("--report-html")
    args = parser.parse_args()

    merged_tsv = Path(args.merged_tsv)
    out_dir = merged_tsv.parent
    report_html = Path(args.report_html) if args.report_html else None

    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import numpy as np

    header = merged_tsv.open().readline().strip().split("\t")
    p_cols = [c for c in header if c not in {"SNP", "CHR", "BP"}]
    if not p_cols:
        raise SystemExit("No trait p-value columns found.")

    base_lf = (
        pl.scan_csv(str(merged_tsv), separator="\t")
        .select(["CHR", "BP", *p_cols])
        .with_columns(
            [
                pl.col("CHR").cast(pl.Int64, strict=False),
                pl.col("BP").cast(pl.Float64, strict=False),
                *[_clamped_p_expr(c) for c in p_cols],
            ]
        )
        .filter(pl.col("CHR").is_not_null() & pl.col("BP").is_not_null())
    )

    chr_stats = (
        base_lf.group_by("CHR")
        .agg(pl.col("BP").max().alias("max_bp"))
        .collect()
        .sort("CHR")
    )
    chr_offsets = {}
    chr_ticks = []
    chr_tick_labels = []
    running_bp = 0.0
    for row in chr_stats.iter_rows(named=True):
        chrom = int(row["CHR"])
        chr_len = float(row["max_bp"])
        chr_offsets[chrom] = running_bp
        chr_ticks.append((running_bp + chr_len / 2.0) / 1e6)
        chr_tick_labels.append(str(chrom))
        running_bp += chr_len + 1e6

    keep_expr = (pl.col("__row") % 20) == 0
    for col in p_cols:
        keep_expr = keep_expr | (pl.col(col) < 0.1)

    plot_df = (
        base_lf.with_row_index("__row")
        .filter(keep_expr)
        .with_columns(
            pl.col("CHR")
            .map_elements(lambda c: chr_offsets.get(int(c), 0.0), return_dtype=pl.Float64)
            .alias("__chr_offset")
        )
        .with_columns((pl.col("BP") + pl.col("__chr_offset")).alias("plot_pos"))
        .collect()
        .sort(["CHR", "BP"])
    )

    max_y = 0.0
    for col in p_cols:
        vals = plot_df[col].drop_nulls()
        if len(vals):
            max_y = max(max_y, float((-vals.log(10)).max()))
    y_limit = max(5.4, min(max_y * 1.18 + 0.4, 30.0))

    band_colors = ["#3F3F46", "#A1A1AA"]
    palette = ["#0072B2", "#D55E00", "#009E73", "#CC79A7", "#E69F00", "#56B4E9", "#000000"]

    fig, axes = plt.subplots(len(p_cols), 1, figsize=(14, max(4.2, 3.45 * len(p_cols))), sharex=True, facecolor="white")
    axes = np.atleast_1d(axes)
    for trait_idx, (ax, trait) in enumerate(zip(axes, p_cols)):
        ax.set_facecolor("#fbfcfe")
        y = -plot_df[trait].log(10)
        df_trait = plot_df.with_columns(y.alias("minus_log10_p")).filter(pl.col("minus_log10_p").is_not_null())
        for idx, chrom in enumerate(chr_tick_labels):
            sub = df_trait.filter(pl.col("CHR") == int(chrom))
            ax.scatter(
                sub["plot_pos"] / 1e6,
                sub["minus_log10_p"],
                color=band_colors[idx % len(band_colors)],
                s=0.3,
                alpha=0.68,
                edgecolors="none",
                rasterized=True,
            )
        ax.axhline(y=-math.log10(1e-5), color="#334155", linestyle=(0, (1, 2)), linewidth=1.2, alpha=0.9)
        ax.axhline(y=-math.log10(1e-8), color="#BE123C", linestyle="--", linewidth=1.4, alpha=0.9)
        ax.set_title(trait.replace("_", " "), fontsize=12, fontweight="bold", pad=7, color="#101828")
        ax.set_ylabel(r"$-\log_{10}(P)$", fontsize=9, color="#344054")
        ax.set_ylim(0, y_limit)
        ax.set_xlim(0, running_bp / 1e6)
        ax.grid(True, axis="y", linestyle=":", linewidth=0.8, color="#d0d5dd", alpha=0.75)
        ax.tick_params(axis="both", labelsize=8, colors="#475467")
    axes[-1].set_xticks(chr_ticks)
    axes[-1].set_xticklabels(chr_tick_labels)
    axes[-1].set_xlabel("Chromosome", fontsize=10, color="#344054")
    fig.suptitle("Multi-Trait Manhattan Plot", fontsize=16, fontweight="bold", color="#101828", y=0.995)
    fig.tight_layout(rect=[0, 0, 1, 0.965])
    stacked = out_dir / "multi_trait_manhattan.png"
    fig.savefig(stacked, dpi=200, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(12.5, 5.3), facecolor="white")
    for idx, trait in enumerate(p_cols):
        y = -plot_df[trait].log(10)
        df_trait = plot_df.with_columns(y.alias("minus_log10_p")).filter(pl.col("minus_log10_p").is_not_null())
        ax.scatter(
            df_trait["plot_pos"] / 1e6,
            df_trait["minus_log10_p"],
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
    ax.set_xlim(0, running_bp / 1e6)
    ax.set_ylim(0, y_limit)
    ax.set_xticks(chr_ticks)
    ax.set_xticklabels(chr_tick_labels)
    ax.grid(True, linestyle=":", linewidth=0.8, color="#d0d5dd", alpha=0.7)
    ax.legend(loc="upper right", fontsize=8, frameon=True, framealpha=0.92, markerscale=4)
    ax.tick_params(axis="both", labelsize=8, colors="#475467")
    fig.tight_layout()
    overlaid = out_dir / "multi_trait_overlaid_manhattan.png"
    fig.savefig(overlaid, dpi=190, bbox_inches="tight")
    plt.close(fig)

    if report_html:
        _replace_html_image(report_html, "Multi-Trait Stacked Manhattan Plot", stacked)
        _replace_html_image(report_html, "Multi-Trait Overlaid Manhattan Plot", overlaid)

    print(f"Retained {len(plot_df):,} plotted rows from {merged_tsv}")
    print(stacked)
    print(overlaid)


if __name__ == "__main__":
    main()
