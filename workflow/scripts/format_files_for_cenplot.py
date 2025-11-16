import os
import argparse
import polars as pl


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

    cdr_bed = os.path.join(outdir, "cdr.bed")
    try:
        df_cdr = pl.read_csv(
            args.bed_cdr,
            separator="\t",
            has_header=False,
            columns=[0, 1, 2],
            new_columns=["chrom", "chrom_st", "chrom_end"],
        )
        # cdr_bed
        df_cdr.filter(pl.col("chrom") == args.ctg).write_csv(
            cdr_bed, separator="\t", include_header=False
        )
    except pl.exceptions.NoDataError:
        _ = open(cdr_bed, "wt")

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
            new_columns=["chrom", "chrom_st", "chrom_end", "name"],
        )
        .with_columns(
            strand=pl.lit("+"),
            name=pl.when(pl.col("name") == "ALR/Alpha")
            .then(pl.lit("α-satellite"))
            .when(pl.col("name").str.contains_any(hsat_repeats))
            .then(pl.lit("Human satellite"))
            .otherwise(pl.lit("Other satellite")),
            color=pl.when(pl.col("name") == "ALR/Alpha")
            .then(pl.lit("#522758"))
            .when(pl.col("name").str.contains_any(hsat_repeats))
            .then(pl.lit("#84bac6"))
            .otherwise(pl.lit("#808080")),
        )
        .select("chrom", "chrom_st", "chrom_end", "name", "strand", "color")
        .sort(by=["chrom", "chrom_st"])
    )
    # Annotate intervals between RM intervals as non-satellite.
    df_uniq_sequence = (
        df_rm.group_by(["chrom"], maintain_order=True)
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
    df_rm_all = (
        pl.concat([df_rm, df_uniq_sequence])
        .sort(by=["chrom", "chrom_st"])
        .filter(pl.col("chrom_st") != pl.col("chrom_end"))
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
    ).select(
        "chrom",
        "chrom_st",
        "chrom_end",
        "name",
        0,
        "strand",
        "thick_chrom_st",
        "thick_chrom_end",
        "color",
    ).write_csv(os.path.join(outdir, "rm.bed"), separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
