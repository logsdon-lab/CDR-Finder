import sys
import math
import argparse
from collections import defaultdict

import polars as pl
from scipy import signal
from intervaltree import Interval, IntervalTree


def get_interval(df: pl.DataFrame, interval: Interval):
    return df.filter((pl.col("st") >= interval.begin) & (pl.col("end") <= interval.end))


def main():
    ap = argparse.ArgumentParser(description="CDR finder.")
    ap.add_argument(
        "-i",
        "--infile",
        default=sys.stdin,
        type=argparse.FileType("rb"),
        required=True,
        help="Average 5mC methylation signal as 4-column bedfile.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="CDR regions as 3-column bedfile.",
    )
    ap.add_argument(
        "--bp_merge", type=int, default=None, help="Base pairs to merge CDRs."
    )
    ap.add_argument(
        "--thr_height_perc_valley",
        type=float,
        default=0.5,
        help="Threshold percent of the median methylation percentage needed as the minimal height/prominence of a valley from the median. Larger values filter for deeper valleys.",
    )
    ap.add_argument(
        "--thr_prom_perc_valley",
        type=float,
        default=0.5,
        help="Threshold percent of the median methylation percentage needed as the minimal prominence of a valley from the median. Larger values filter for prominent valleys.",
    )
    ap.add_argument(
        "--bp_edge",
        type=int,
        default=5_000,
        help="Bases to look on both edges of cdr to determine effective height.",
    )

    args = ap.parse_args()
    df = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "avg"],
    )

    cdr_intervals: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    for chrom, df_chr_methyl in df.group_by(["chrom"]):
        chrom: str = chrom[0]

        # Group adjacent, contiguous intervals.
        df_chr_methyl_adj_groups = (
            df_chr_methyl.with_columns(brk=pl.col("end") == pl.col("st").shift(-1))
            .fill_null(True)
            .with_columns(pl.col("brk").rle_id())
            # Group contiguous intervals.
            .with_columns(
                pl.when(pl.col("brk") % 2 == 0)
                .then(pl.col("brk") + 1)
                .otherwise(pl.col("brk"))
            )
            .partition_by("brk")
        )
        cdr_prom_thr = (
            df_chr_methyl["avg"].median() * args.thr_prom_perc_valley
            if args.thr_prom_perc_valley
            else None
        )
        cdr_height_thr = df_chr_methyl["avg"].median() * args.thr_height_perc_valley

        # Find peaks within the signal per group.
        for df_chr_methyl_adj_grp in df_chr_methyl_adj_groups:
            df_chr_methyl_adj_grp = df_chr_methyl_adj_grp.with_row_index()

            # Require valley has prominence of some percentage of median methyl signal.
            # Invert for peaks.
            _, peak_info = signal.find_peaks(
                -df_chr_methyl_adj_grp["avg"], width=1, prominence=cdr_prom_thr
            )

            grp_cdr_intervals = []
            for cdr_st_idx, cdr_end_idx in zip(
                peak_info["left_ips"], peak_info["right_ips"]
            ):
                # Convert approx indices to indices
                cdr_st = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.floor(cdr_st_idx)
                ).row(0, named=True)["st"]
                cdr_end = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.ceil(cdr_end_idx)
                ).row(0, named=True)["end"]

                grp_cdr_intervals.append(Interval(cdr_st, cdr_end))

            for interval in grp_cdr_intervals:
                cdr_st, cdr_end = interval.begin, interval.end

                # Get left and right side of CDR.
                df_cdr = get_interval(df_chr_methyl_adj_grp, interval)
                df_cdr_left = get_interval(
                    df_chr_methyl_adj_grp, Interval(cdr_st - args.bp_edge, cdr_st)
                )
                df_cdr_right = get_interval(
                    df_chr_methyl_adj_grp, Interval(cdr_end, cdr_end + args.bp_edge)
                )

                cdr_low = df_cdr["avg"].min()
                cdr_right_max = df_cdr_right["avg"].max()
                cdr_left_max = df_cdr_left["avg"].max()
                cdr_edge_height = max(
                    cdr_right_max if cdr_right_max else 0,
                    cdr_left_max if cdr_left_max else 0,
                )

                # Calculate the height of this CDR looking at edges.
                cdr_height = max(0, cdr_edge_height - cdr_low)

                # Add merge distance bp.
                if args.bp_merge:
                    cdr_st = cdr_st - args.bp_merge
                    cdr_end = cdr_end + args.bp_merge

                if cdr_height >= cdr_height_thr:
                    cdr_intervals[chrom].add(Interval(cdr_st, cdr_end))

    # Merge overlaps and output.
    for chrom, cdrs in cdr_intervals.items():
        if args.bp_merge:
            starting_intervals = len(cdrs)
            cdrs.merge_overlaps()
            print(
                f"Merged {starting_intervals - len(cdrs)} intervals in {chrom}.",
                file=sys.stderr,
            )

        for cdr in cdrs.iter():
            cdr_st, cdr_end = cdr.begin, cdr.end
            if args.bp_merge:
                cdr_st += args.bp_merge
                cdr_end -= args.bp_merge

            args.outfile.write(f"{chrom}\t{cdr_st}\t{cdr_end}\n")


if __name__ == "__main__":
    raise SystemExit(main())
