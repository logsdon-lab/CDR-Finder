import sys
import argparse
import polars as pl


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        help="Input repeatmasker file.",
        required=True,
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-r",
        "--rename_key",
        help="Rename key for repeatmasker.",
        required=True,
        type=argparse.FileType("rb"),
    )
    ap.add_argument(
        "-o",
        "--outfile",
        help="Repeatmasker bed.",
        type=argparse.FileType("wt"),
        default=sys.stdout,
    )
    args = ap.parse_args()

    infile = args.infile
    outfile = args.outfile
    rename_key = args.rename_key

    df_infile = pl.read_csv(
        infile,
        separator="\t",
        truncate_ragged_lines=True,
        skip_rows=3,
        has_header=False,
        columns=[4, 5, 6, 9],
        new_columns=["chrom", "chrom_st", "chrom_end", "rtype"],
    )
    df_rename_key = (
        pl.read_csv(
            rename_key,
            separator="\t",
            has_header=False,
            new_columns=["old", "new"],
        )
        .with_columns(
            mtch_chrom=pl.col("new").str.extract_groups(
                r"^(?<ctg>.+):(?<ctg_st>\d+)-(?<ctg_end>\d+)$"
            )
        )
        .unnest("mtch_chrom")
        .cast({"ctg_st": pl.Int64, "ctg_end": pl.Int64})
    )

    # Rename back
    df_outfile = (
        df_infile.join(df_rename_key, left_on="chrom", right_on="old", how="left")
        .with_columns(
            chrom=pl.col("ctg"),
            chrom_st=pl.col("chrom_st") + pl.col("ctg_st"),
            chrom_end=pl.col("chrom_end") + pl.col("ctg_st"),
        )
        .drop("ctg", "ctg_st", "ctg_end", "new")
    )

    df_outfile.write_csv(outfile, separator="\t", include_header=False)


if __name__ == "__main__":
    raise SystemExit(main())
