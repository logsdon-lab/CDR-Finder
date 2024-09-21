import sys
import math
import argparse
from collections import defaultdict

import polars as pl
from scipy import signal
from intervaltree import Interval, IntervalTree


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
        default=0.15,
        help="Threshold percent of the median methylation percentage needed as the minimal height/prominence of a valley from the median. Larger values filter for deeper valleys.",
    )
    ap.add_argument(
        "--thr_prom_perc_valley",
        type=float,
        default=0.3,
        help="Threshold percent of the median methylation percentage needed as the minimal prominence of a valley from the median. Larger values filter for deeper valleys and remove smaller valleys at edges.",
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

        thr_valley_prom = df_chr_methyl["avg"].median() * args.thr_prom_perc_valley

        # Find peaks within the signal per group.
        for df_chr_methyl_adj_grp in df_chr_methyl_adj_groups:
            df_chr_methyl_adj_grp = df_chr_methyl_adj_grp.with_row_index()

            # Require valley has prominence of some percentage of median methyl signal.
            # Invert for peaks.
            _, peak_info = signal.find_peaks(
                -df_chr_methyl_adj_grp["avg"], width=1, prominence=thr_valley_prom
            )

            grp_cdr_intervals = []
            for i, (cdr_st_idx, cdr_end_idx) in enumerate(
                zip(peak_info["left_ips"], peak_info["right_ips"])
            ):
                # Convert approx indices to indices
                cdr_st = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.floor(cdr_st_idx)
                ).row(0, named=True)["st"]
                cdr_end = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.ceil(cdr_end_idx)
                ).row(0, named=True)["end"]

                # Add merge distance bp.
                grp_cdr_intervals.append(
                    Interval(cdr_st, cdr_end, -peak_info["width_heights"][i])
                )

            # Calculate median of group removing valleys
            exprs_valley_intervals = [
                (pl.col("st") >= i.begin) & (pl.col("end") <= i.end)
                for i in grp_cdr_intervals
            ]
            try:
                median_noncdr = (
                    df_chr_methyl_adj_grp.filter(
                        ~pl.any_horizontal(exprs_valley_intervals)
                    )
                    .get_column("avg")
                    .median()
                )
                cdr_height_thr = median_noncdr * args.thr_height_perc_valley
            except Exception:
                # No CDRs in this region.
                continue

            for interval in grp_cdr_intervals:
                cdr_st, cdr_end, cdr_low = interval.begin, interval.end, interval.data
                # Calculate the height of this CDR from the median.
                cdr_height_from_median = max(0, median_noncdr - cdr_low)

                if args.bp_merge:
                    cdr_st = cdr_st - args.bp_merge
                    cdr_end = cdr_end + args.bp_merge

                if cdr_height_from_median >= cdr_height_thr:
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
