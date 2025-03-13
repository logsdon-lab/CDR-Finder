import os
import argparse
import itertools
import polars as pl

IDENT_CUTOFF = 97.5
IDENT_INCREMENT = 0.25
IDENT_COLORS = [
    "#4b3991",
    "#2974af",
    "#4a9da8",
    "#57b894",
    "#9dd893",
    "#e1f686",
    "#ffffb2",
    "#fdda79",
    "#fb9e4f",
    "#ee5634",
    "#c9273e",
    "#8a0033",
]

IDENT_RANGE_ENDS = tuple(
    # Increment by IDENT_INCREMENT from IDENT_CUTOFF up to 100%
    i + IDENT_CUTOFF
    for i in itertools.accumulate(IDENT_INCREMENT for _ in range(10))
)
IDENT_RANGE = [
    (0, 90),
    (90, 97.5),
    *zip((IDENT_CUTOFF, *IDENT_RANGE_ENDS[:-1]), IDENT_RANGE_ENDS),
]

IDENT_COLOR_RANGE: dict[tuple[float, float], str] = dict(zip(IDENT_RANGE, IDENT_COLORS))
BED9_COLS = [
    "chrom",
    "chrom_st",
    "chrom_end",
    "name",
    "score",
    "strand",
    "thick_chrom_st",
    "thick_chrom_end",
    "color",
]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--ctg", required=True, help="Contig name to filter for.")
    ap.add_argument(
        "--bed_binned_freq",
        required=True,
        help="Bedfile for binned methylation frequency.",
    )
    ap.add_argument("--bed_cdr", required=True, help="Bedfile for CDR calls.")
    ap.add_argument(
        "--bed_rm", required=True, help="Bedfile for repeatmasker annotations."
    )
    ap.add_argument("--bed_ident", required=True, help="Bedfile for self-identity.")
    ap.add_argument("--outdir", required=True, help="Output directory.")

    args = ap.parse_args()

    outdir = os.path.join(args.outdir, args.ctg)
    os.makedirs(outdir, exist_ok=True)

    df_binned_freq = pl.read_csv(
        args.bed_binned_freq,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "chrom_st", "chrom_end", "avg_methyl", "cov"],
    ).with_columns(methyl_cov=pl.col("cov") * (pl.col("avg_methyl") / 100))
    df_cdr = pl.read_csv(
        args.bed_cdr,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "chrom_st", "chrom_end"],
    )
    hsat_repeats = [
        "SAR",
        "HSAT",
        "(CATTC)n",
        "(GAATG)n",
    ]
    df_rm = (
        pl.read_csv(
            args.bed_rm,
            separator="\t",
            has_header=False,
            new_columns=["chrom", "chrom_st", "chrom_end", "rtype", "rclass", "strand"],
        )
        .with_columns(
            strand=pl.when(pl.col("strand") == "C")
            .then(pl.lit("-"))
            .otherwise(pl.lit("+")),
            name=pl.when(pl.col("rtype") == "ALR/Alpha")
            .then(pl.lit("α-satellite"))
            .when(pl.col("rtype").str.contains_any(hsat_repeats))
            .then(pl.lit("Human satellite"))
            .otherwise(pl.lit("Other satellite")),
            color=pl.when(pl.col("rtype") == "ALR/Alpha")
            .then(pl.lit("#522758"))
            .when(pl.col("rtype").str.contains_any(hsat_repeats))
            .then(pl.lit("#84bac6"))
            .otherwise(pl.lit("#808080")),
        )
        .select("chrom", "chrom_st", "chrom_end", "name", "strand", "color")
    )
    # Annotate intervals between RM intervals as non-satellite.
    df_uniq_sequence = (
        df_rm.group_by(["chrom"])
        .agg(chrom_st=pl.col("chrom_end"), chrom_end=pl.col("chrom_st").shift(-1))
        .explode(["chrom_st", "chrom_end"])
        .with_columns(
            strand=pl.lit("+"),
            name=pl.lit("Non-satellite"),
            color=pl.lit("#FFFFFF"),
        )
        .select("chrom", "chrom_st", "chrom_end", "name", "strand", "color")
        .drop_nulls()
    )
    df_rm_all = pl.concat([df_rm, df_uniq_sequence]).sort(by=["chrom", "chrom_st"])
    # cdr_bed
    df_cdr.filter(pl.col("chrom") == args.ctg).write_csv(
        os.path.join(outdir, "cdr.bed"), separator="\t", include_header=False
    )
    df_binned_freq.filter(pl.col("chrom") == args.ctg).select(
        "chrom", "chrom_st", "chrom_end", "cov"
    ).write_csv(os.path.join(outdir, "cov.bed"), separator="\t", include_header=False)
    df_binned_freq.filter(pl.col("chrom") == args.ctg).select(
        "chrom", "chrom_st", "chrom_end", "methyl_cov"
    ).write_csv(
        os.path.join(outdir, "methyl_cov.bed"), separator="\t", include_header=False
    )
    df_binned_freq.filter(pl.col("chrom") == args.ctg).select(
        "chrom", "chrom_st", "chrom_end", "avg_methyl"
    ).write_csv(
        os.path.join(outdir, "methyl_avg_cov.bed"), separator="\t", include_header=False
    )
    df_rm_all.filter(pl.col("chrom") == args.ctg).with_columns(
        thick_chrom_st=pl.col("chrom_st"),
        thick_chrom_end=pl.col("chrom_end"),
        score=pl.lit(0),
    ).select(BED9_COLS).write_csv(
        os.path.join(outdir, "rm.bed"), separator="\t", include_header=False
    )

    df_self_ident = pl.read_csv(
        args.bed_ident,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "chrom_st", "chrom_end", "score"],
    )
    # Build expr
    expr_colors = None
    for (rng_st, rng_end), color in IDENT_COLOR_RANGE.items():
        rng_name = f"{rng_st}% - {rng_end}%"
        color_name = pl.lit(f"{color}_{rng_name}")
        expr_in_rng = pl.col("score").is_between(rng_st, rng_end)
        if isinstance(expr_colors, pl.Expr):
            expr_colors = expr_colors.when(expr_in_rng).then(color_name)
        else:
            expr_colors = pl.when(expr_in_rng).then(color_name)

    df_self_ident = df_self_ident.with_columns(
        strand=pl.lit("+"),
        thick_chrom_st=pl.col("chrom_st"),
        thick_chrom_end=pl.col("chrom_end"),
        color_name=expr_colors.otherwise(pl.lit("white_none"))
        .str.split_exact("_", 1)
        .struct.rename_fields(["color", "name"]),
    ).unnest("color_name")

    df_self_ident.select(BED9_COLS).write_csv(
        os.path.join(outdir, "ident.bed"), separator="\t", include_header=False
    )


if __name__ == "__main__":
    raise SystemExit(main())
