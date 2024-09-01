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
        "--window_median_filt",
        type=int,
        default=None,
        help="Window in number of rows to apply median filter to denoise.",
    )
    ap.add_argument(
        "--bp_merge", type=int, default=None, help="Base pairs to merge CDRs."
    )
    ap.add_argument(
        "--thr_quantile_valley",
        type=float,
        default=0.1,
        help="Threshold quantile to filter low confidence CDRs. ex 0.1 allows valleys with average value less than the 10% percentile.",
    )
    ap.add_argument(
        "--thr_prominence_perc_valley",
        type=float,
        default=0.2,
        help="Threshold percent of the maximum methylation percentage as the minimal prominence of a valley to filter low confidence CDRs.",
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
        thr_valley_prom = df_chr_methyl["avg"].max() * args.thr_prominence_perc_valley
        thr_valley_quantile = df_chr_methyl["avg"].quantile(args.thr_quantile_valley)

        # Per group apply optionally apply a median filter to denoise.
        # And find peaks within the signal.
        for df_chr_methyl_adj_grp in df_chr_methyl_adj_groups:
            df_chr_methyl_adj_grp = df_chr_methyl_adj_grp.with_row_index()

            methyl_signal = df_chr_methyl_adj_grp["avg"]
            if args.window_median_filt:
                methyl_signal = signal.medfilt(methyl_signal, args.window_median_filt)

            # Require valley has prominence of some percentage of max methyl signal.
            # Invert for peaks.
            peaks, peak_info = signal.find_peaks(
                -methyl_signal, width=1, prominence=thr_valley_prom
            )

            df_filtered_regions = (
                df_chr_methyl_adj_grp.filter(pl.col("index").is_in(peaks))
                .rename({"index": "index_original"})
                .with_row_index()
                # Filter peaks less than desired quantile.
                .filter(pl.col("avg") < thr_valley_quantile)
            )
            filtered_region_idxs = set(df_filtered_regions["index"])

            for i, (cdr_st_idx, cdr_end_idx) in enumerate(
                zip(peak_info["left_ips"], peak_info["right_ips"])
            ):
                if i not in filtered_region_idxs:
                    continue

                # Convert approx indices to indices
                cdr_st = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.floor(cdr_st_idx)
                ).row(0, named=True)["st"]
                cdr_end = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.floor(cdr_end_idx)
                ).row(0, named=True)["end"]

                if args.bp_merge:
                    cdr_st = cdr_st - args.bp_merge
                    cdr_end = cdr_end + args.bp_merge

                # Add merge distance bp.
                cdr_intervals[chrom].add(Interval(cdr_st, cdr_end))

    # Merge overlaps and output.
    for chrom, cdrs in cdr_intervals.items():
        if args.bp_merge:
            starting_intervals = len(cdrs)
            cdrs.merge_overlaps()
            print(
                f"Merged {starting_intervals - len(cdrs)} intervals.", file=sys.stderr
            )

        for cdr in cdrs.iter():
            cdr_st, cdr_end = cdr.begin, cdr.end
            if args.bp_merge:
                cdr_st += args.bp_merge
                cdr_end -= args.bp_merge

            args.outfile.write(f"{chrom}\t{cdr_st}\t{cdr_end}\n")


if __name__ == "__main__":
    raise SystemExit(main())
