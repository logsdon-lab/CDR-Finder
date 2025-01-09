import sys
import math
import json
import argparse
import statistics
from collections import defaultdict
from typing import Iterable

import polars as pl
from scipy import signal
from intervaltree import Interval, IntervalTree


def get_interval(
    df: pl.DataFrame,
    interval: Interval,
    ignore_intervals: Iterable[Interval] | None = None,
) -> pl.DataFrame:
    df_res = df.filter(
        (pl.col("st") >= interval.begin) & (pl.col("end") <= interval.end)
    )
    if ignore_intervals:
        df_res = df_res.with_columns(
            ignore=pl.when(
                pl.any_horizontal(
                    (pl.col("st") >= interval.begin) & (pl.col("end") <= interval.end)
                    for interval in ignore_intervals
                )
            )
            .then(pl.lit(True))
            .otherwise(pl.lit(False))
        )
    else:
        df_res = df_res.with_columns(ignore=pl.lit(False))

    return df_res


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
        default=None,
        help="Threshold percent of the median methylation percentage needed as the minimal prominence of a valley from the median. Larger values filter for prominent valleys.",
    )
    ap.add_argument(
        "--bp_edge",
        type=int,
        default=5_000,
        help="Bases to look on both edges of cdr to determine effective height.",
    )
    ap.add_argument(
        "--extend_edges_std",
        type=int,
        default=None,
        help="Extend edges of CDR until the mean signal with the given stdev is reached. Adds smaller, less prominent CDRs if missed. A value of 0 is the mean while +1/-1 is one stdev above/below the mean.",
    )
    ap.add_argument(
        "--edge_height_heuristic",
        type=str,
        choices=["min", "max", "avg"],
        default="min",
        help="Heuristic used to determine edge height of CDR.",
    )
    ap.add_argument(
        "--override_chrom_params",
        type=str,
        default=None,
        help="JSON file with params per chromosome name to override default settings. Each parameter name should match.",
    )

    args = ap.parse_args()
    df = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "avg"],
    )
    if args.override_chrom_params:
        with open(args.override_chrom_params, "rt") as fh:
            override_params = json.load(fh)
    else:
        override_params = {}

    cdr_intervals: defaultdict[str, IntervalTree] = defaultdict(IntervalTree)
    for chrom, df_chr_methyl in df.group_by(["chrom"]):
        chrom: str = chrom[0]

        # Get override params.
        thr_prom_perc_valley = override_params.get(chrom, {}).get(
            "thr_prom_perc_valley", args.thr_prom_perc_valley
        )
        thr_height_perc_valley = override_params.get(chrom, {}).get(
            "thr_height_perc_valley", args.thr_height_perc_valley
        )
        bp_edge = override_params.get(chrom, {}).get("bp_edge", args.bp_edge)
        edge_height_heuristic = override_params.get(chrom, {}).get(
            "edge_height_heuristic", args.edge_height_heuristic
        )
        extend_edges_std = override_params.get(chrom, {}).get(
            "extend_edges_std", args.extend_edges_std
        )
        bp_merge = override_params.get(chrom, {}).get("bp_merge", args.bp_merge)

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
        avg_methyl_median = df_chr_methyl["avg"].median()
        avg_methyl_mean = df_chr_methyl["avg"].mean()
        avg_methyl_std = df_chr_methyl["avg"].std()
        cdr_prom_thr = (
            avg_methyl_median * thr_prom_perc_valley if thr_prom_perc_valley else None
        )
        cdr_height_thr = avg_methyl_median * thr_height_perc_valley
        print(
            f"Using CDR height threshold of {cdr_height_thr} and prominence threshold of {cdr_prom_thr} for {chrom}.",
            file=sys.stderr,
        )

        # Find peaks within the signal per group.
        for df_chr_methyl_adj_grp in df_chr_methyl_adj_groups:
            df_chr_methyl_adj_grp = df_chr_methyl_adj_grp.with_row_index()

            # Require valley has prominence of some percentage of median methyl signal.
            # Invert for peaks.
            _, peak_info = signal.find_peaks(
                -df_chr_methyl_adj_grp["avg"], width=1, prominence=cdr_prom_thr
            )

            grp_cdr_intervals: set[Interval] = set()
            for cdr_st_idx, cdr_end_idx, cdr_prom in zip(
                peak_info["left_ips"], peak_info["right_ips"], peak_info["prominences"]
            ):
                # Convert approx indices to indices
                cdr_st = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.floor(cdr_st_idx)
                ).row(0, named=True)["st"]
                cdr_end = df_chr_methyl_adj_grp.filter(
                    pl.col("index") == math.ceil(cdr_end_idx)
                ).row(0, named=True)["end"]

                grp_cdr_intervals.add(Interval(cdr_st, cdr_end, cdr_prom))

            for interval in grp_cdr_intervals:
                cdr_st, cdr_end, cdr_prom = interval.begin, interval.end, interval.data
                ignore_intervals = grp_cdr_intervals.difference([interval])
                df_cdr = get_interval(df_chr_methyl_adj_grp, interval)

                interval_cdr_left = Interval(cdr_st - bp_edge, cdr_st)
                interval_cdr_right = Interval(cdr_end, cdr_end + bp_edge)

                # Get left and right side of CDR.
                # Subtract intervals if overlapping bp edge region.
                # Set ignored intervals on sides of CDR to average methylation median.
                # This does not affect other calls and is just to look at valley in isolation.
                df_cdr_left = get_interval(
                    df_chr_methyl_adj_grp, interval_cdr_left, ignore_intervals
                ).with_columns(
                    avg=pl.when(pl.col("ignore"))
                    .then(avg_methyl_median)
                    .otherwise(pl.col("avg"))
                )
                df_cdr_right = get_interval(
                    df_chr_methyl_adj_grp, interval_cdr_right, ignore_intervals
                ).with_columns(
                    avg=pl.when(pl.col("ignore"))
                    .then(avg_methyl_median)
                    .otherwise(pl.col("avg"))
                )

                cdr_low = df_cdr["avg"].min()
                cdr_right_median = df_cdr_right["avg"].median()
                cdr_left_median = df_cdr_left["avg"].median()

                # If empty, use median.
                edge_heights = [
                    cdr_right_median
                    if cdr_right_median
                    else df_chr_methyl_adj_grp["avg"].median(),
                    cdr_left_median
                    if cdr_left_median
                    else df_chr_methyl_adj_grp["avg"].median(),
                ]
                if edge_height_heuristic == "min":
                    cdr_edge_height = min(edge_heights)
                elif edge_height_heuristic == "max":
                    cdr_edge_height = max(edge_heights)
                else:
                    cdr_edge_height = statistics.mean(edge_heights)

                # Calculate the height of this CDR looking at edges.
                cdr_height = cdr_edge_height - cdr_low

                # Ignore CDR if less than height.
                if cdr_height < cdr_height_thr:
                    continue

                print(
                    f"Found CDR at {chrom}:{interval.begin}-{interval.end} with height of {cdr_height} and prominence {cdr_prom}.",
                    file=sys.stderr,
                )

                if isinstance(extend_edges_std, int):
                    edge_thr = avg_methyl_mean + (extend_edges_std * avg_methyl_std)
                    try:
                        cdr_st = df_cdr_left.filter(
                            (pl.col("st") < cdr_st) & (pl.col("avg") >= edge_thr)
                        ).row(-1, named=True)["st"]
                    except pl.exceptions.OutOfBoundsError:
                        cdr_st = cdr_st
                    try:
                        cdr_end = df_cdr_right.filter(
                            (pl.col("end") > cdr_end) & (pl.col("avg") >= edge_thr)
                        ).row(0, named=True)["end"]
                    except pl.exceptions.OutOfBoundsError:
                        cdr_end = cdr_end

                # Add merge distance bp.
                if bp_merge:
                    cdr_st = cdr_st - bp_merge
                    cdr_end = cdr_end + bp_merge

                cdr_intervals[chrom].add(Interval(cdr_st, cdr_end))

    # Merge overlaps and output.
    for chrom, cdrs in cdr_intervals.items():
        if bp_merge:
            starting_intervals = len(cdrs)
            cdrs.merge_overlaps()
            print(
                f"Merged {starting_intervals - len(cdrs)} intervals in {chrom}.",
                file=sys.stderr,
            )

        for cdr in cdrs.iter():
            cdr_st, cdr_end = cdr.begin, cdr.end
            if bp_merge:
                cdr_st += bp_merge
                cdr_end -= bp_merge

            args.outfile.write(f"{chrom}\t{cdr_st}\t{cdr_end}\n")


if __name__ == "__main__":
    raise SystemExit(main())
