#!/usr/bin/env python3
"""Summarize SNPs whose significance changes after mtCOJO conditioning."""

import argparse
import os
from typing import Dict, Iterable, List, Optional

import polars as pl


DEFAULT_THRESHOLDS = (1e-5, 1e-8)
ID_COLUMNS = ("SNP", "CHR", "BP")


def _schema_names(lf: pl.LazyFrame) -> List[str]:
    try:
        return list(lf.collect_schema().names())
    except AttributeError:
        return list(lf.schema.keys())


def _safe_label(value: float) -> str:
    return f"{value:.0e}".replace("+", "")


def _pick_conditioned_column(columns: Iterable[str], requested: Optional[str] = None) -> str:
    names = list(columns)
    if requested:
        if requested not in names:
            raise ValueError(f"Conditioned p-value column not found: {requested}")
        return requested

    conditioned = [c for c in names if "conditioned" in c.lower()]
    if not conditioned:
        raise ValueError("No conditioned p-value column found in merged GWAS summary.")
    return conditioned[0]


def build_conditional_shift_summary(
    merged_tsv: str,
    out_dir: str,
    thresholds: Iterable[float] = DEFAULT_THRESHOLDS,
    conditioned_col: Optional[str] = None,
) -> Dict[str, str]:
    """Write summary/detail TSVs for SNPs gained or lost after conditioning.

    A lost SNP is significant in the original GWAS trait but not in the
    conditioned target. A gained SNP is the opposite.
    """
    if not merged_tsv or not os.path.exists(merged_tsv):
        raise FileNotFoundError(f"Merged GWAS summary not found: {merged_tsv}")

    os.makedirs(out_dir, exist_ok=True)
    summary_tsv = os.path.join(out_dir, "conditional_significance_shift_summary.tsv")
    detail_tsv = os.path.join(out_dir, "conditional_significance_shift_snps.tsv")

    lf = pl.scan_csv(merged_tsv, separator="\t")
    columns = _schema_names(lf)
    conditioned = _pick_conditioned_column(columns, conditioned_col)
    p_columns = [c for c in columns if c not in ID_COLUMNS]
    comparison_traits = [c for c in p_columns if c != conditioned and "conditioned" not in c.lower()]

    summary_rows = []
    detail_frames = []
    for trait in comparison_traits:
        valid = (
            lf.select(
                [c for c in ID_COLUMNS if c in columns]
                + [
                    pl.col(trait).cast(pl.Float64, strict=False).alias("gwas_p"),
                    pl.col(conditioned).cast(pl.Float64, strict=False).alias("conditioned_p"),
                ]
            )
            .filter(
                pl.col("gwas_p").is_not_null()
                & pl.col("conditioned_p").is_not_null()
                & (pl.col("gwas_p") > 0)
                & (pl.col("conditioned_p") > 0)
            )
        )

        for threshold in thresholds:
            lost = (pl.col("gwas_p") <= threshold) & (pl.col("conditioned_p") > threshold)
            gained = (pl.col("gwas_p") > threshold) & (pl.col("conditioned_p") <= threshold)
            both = (pl.col("gwas_p") <= threshold) & (pl.col("conditioned_p") <= threshold)
            counts = valid.select(
                [
                    pl.len().alias("valid_variant_n"),
                    (pl.col("gwas_p") <= threshold).sum().alias("gwas_significant_n"),
                    (pl.col("conditioned_p") <= threshold).sum().alias("conditioned_significant_n"),
                    lost.sum().alias("lost_after_conditioning_n"),
                    gained.sum().alias("gained_after_conditioning_n"),
                    both.sum().alias("both_significant_n"),
                ]
            ).collect()
            row = counts.row(0, named=True)
            summary_rows.append(
                {
                    "comparison_trait": trait,
                    "conditioned_trait": conditioned,
                    "threshold": threshold,
                    **row,
                }
            )

            for direction, predicate in (
                ("lost_after_conditioning", lost),
                ("gained_after_conditioning", gained),
            ):
                detail_frames.append(
                    valid.filter(predicate)
                    .with_columns(
                        [
                            pl.lit(trait).alias("comparison_trait"),
                            pl.lit(conditioned).alias("conditioned_trait"),
                            pl.lit(threshold).alias("threshold"),
                            pl.lit(_safe_label(threshold)).alias("threshold_label"),
                            pl.lit(direction).alias("direction"),
                            (pl.col("conditioned_p") / pl.col("gwas_p")).alias("conditioned_to_gwas_p_ratio"),
                        ]
                    )
                    .select(
                        [
                            "comparison_trait",
                            "conditioned_trait",
                            "threshold",
                            "threshold_label",
                            "direction",
                            *[c for c in ID_COLUMNS if c in columns],
                            "gwas_p",
                            "conditioned_p",
                            "conditioned_to_gwas_p_ratio",
                        ]
                    )
                )

    pl.DataFrame(summary_rows).write_csv(summary_tsv, separator="\t")
    if detail_frames:
        pl.concat(detail_frames, how="vertical_relaxed").collect().write_csv(detail_tsv, separator="\t")
    else:
        pl.DataFrame(
            schema={
                "comparison_trait": pl.String,
                "conditioned_trait": pl.String,
                "threshold": pl.Float64,
                "threshold_label": pl.String,
                "direction": pl.String,
                "SNP": pl.String,
                "CHR": pl.String,
                "BP": pl.String,
                "gwas_p": pl.Float64,
                "conditioned_p": pl.Float64,
                "conditioned_to_gwas_p_ratio": pl.Float64,
            }
        ).write_csv(detail_tsv, separator="\t")

    return {
        "summary_tsv": summary_tsv,
        "detail_tsv": detail_tsv,
        "conditioned_trait": conditioned,
    }


def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "Create TSV tables of SNPs that lose or gain threshold significance "
            "after mtCOJO conditioning."
        )
    )
    parser.add_argument("--merged-tsv", required=True, help="Path to 05_plots_and_tables/merged_gwas_summary.tsv.")
    parser.add_argument("--out-dir", required=True, help="Directory where conditional shift TSV files will be written.")
    parser.add_argument("--conditioned-col", default=None, help="Conditioned p-value column. Default: first column containing 'conditioned'.")
    parser.add_argument(
        "--threshold",
        dest="thresholds",
        action="append",
        type=float,
        default=None,
        help="Significance threshold to compare. May be repeated. Default: 1e-5 and 1e-8.",
    )
    args = parser.parse_args()
    result = build_conditional_shift_summary(
        merged_tsv=args.merged_tsv,
        out_dir=args.out_dir,
        thresholds=args.thresholds or DEFAULT_THRESHOLDS,
        conditioned_col=args.conditioned_col,
    )
    print(result["summary_tsv"])
    print(result["detail_tsv"])


if __name__ == "__main__":
    main()
