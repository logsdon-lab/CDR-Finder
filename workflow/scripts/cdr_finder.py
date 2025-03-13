import sys
import json
import argparse
from collections import deque
from typing import Callable, Iterable

import polars as pl
import numpy as np
from scipy.stats import chi2
from scipy.spatial import distance
from scipy.linalg import inv
from scipy.signal import peak_prominences
from intervaltree import Interval, IntervalTree


def fn_cmp_def(_itv_1: Interval, _itv_2: Interval) -> bool:
    return True


def fn_merge_itv_def(itv_1: Interval, itv_2: Interval) -> Interval:
    return Interval(begin=itv_1.begin, end=itv_2.end, data=None)


def merge_itvs(
    itvs: Iterable[Interval],
    dst: int = 1,
    fn_cmp: Callable[[Interval, Interval], bool] | None = None,
    fn_merge_itv: Callable[[Interval, Interval], Interval] | None = None,
) -> list[Interval]:
    if not fn_cmp:
        fn_cmp = fn_cmp_def
    if not fn_merge_itv:
        fn_merge_itv = fn_merge_itv_def

    final_itvs = []
    sorted_itvs = deque(sorted(itvs))
    while sorted_itvs:
        try:
            itv_1 = sorted_itvs.popleft()
        except IndexError:
            break
        try:
            itv_2 = sorted_itvs.popleft()
        except IndexError:
            final_itvs.append(itv_1)
            break
        dst_between = itv_2.begin - itv_1.end
        passes_cmp = fn_cmp(itv_1, itv_2)
        if dst_between <= dst and passes_cmp:
            sorted_itvs.appendleft(fn_merge_itv(itv_1, itv_2))
        else:
            final_itvs.append(itv_1)
            sorted_itvs.appendleft(itv_2)

    return final_itvs


def detect_outlier_regions(
    df: pl.DataFrame, cols: Iterable[str] = ("avg", "prom"), p_val: float = 0.05
) -> IntervalTree:
    """
    Use Mahalanobis Distance to detect outlier regions based on prominence and average methylation of peaks.
    * multivariate equivalent of the Euclidean distance.

    Adapted idea from https://www.sciencedirect.com/science/article/pii/S1746809419303404

    > Mahalanobis distance (dm) was employed to detect outliers in the multivariate dataset.
    > dm enables a powerful statistical technique to determine how likely (or unlikely)
    > an outlier lies at a specific distance away (or closer) from the center of the distribution.

    """
    itree = IntervalTree()
    # Use average methylation and prominence of each peak.
    df_coords = df.select("st", "end")
    vals = df.select(cols).to_numpy()
    means = np.mean(vals, axis=0)
    try:
        inv_cov = inv(np.cov(vals.T))
    except Exception:
        return itree
    # https://stackoverflow.com/a/66856690
    dst = [distance.mahalanobis(row, means, inv_cov) for row in vals]
    # https://www.machinelearningplus.com/statistics/mahalanobis-distance/
    # Filter coords that are not significant assuming distances follow chi-square distribution.
    dof = len(cols) - 1
    df_coords = df_coords.with_columns(p_val=1 - chi2.cdf(dst, dof))
    for itv in df_coords.filter(pl.col("p_val") < p_val).iter_rows(named=True):
        itree.add(Interval(itv["st"], itv["end"], itv["p_val"]))
    return itree


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
        "--thr_pvalue", type=float, default=0.05, help="Required p-value to call CDR."
    )
    ap.add_argument(
        "--thr_zscore",
        type=float,
        default=-1.8,
        help="Number of z-scores to find/extend CDRs. Should be negative.",
    )
    ap.add_argument("--thr_cdr", type=int, default=1, help="Minimum CDR size.")
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

    # Human centromeres generally have a plateau with some number of dips.
    #   _________     __
    #  /         \   /  \_
    # /           \_/     \
    for chrom, df_chr_methyl in df.group_by(["chrom"]):
        chrom: str = chrom[0]

        # Get override params.
        bp_merge = override_params.get(chrom, {}).get("bp_merge", args.bp_merge)
        if not bp_merge:
            bp_merge = 0

        thr_zscore = override_params.get(chrom, {}).get("thr_zscore", args.thr_zscore)
        if thr_zscore > 0:
            raise ValueError(f"z-score threshold for {chrom} is non-negative.")

        thr_cdr = override_params.get(chrom, {}).get("thr_cdr", args.thr_cdr)
        thr_pvalue = override_params.get(chrom, {}).get("thr_pvalue", args.thr_pvalue)

        df_chr_methyl = df_chr_methyl.with_row_index()

        # Calculate prominences.
        proms, _, _ = peak_prominences(-df_chr_methyl["avg"], df_chr_methyl["index"])
        df_chr_methyl = df_chr_methyl.with_columns(prom=pl.Series(proms))

        # Group adjacent, contiguous intervals.
        # Split regions with no methylation due to N's.
        # Allow a jump of ~8000 bp (rough size of LINE element)
        df_chr_methyl_adj_groups = (
            df_chr_methyl.filter(pl.col("avg") != 0.0)
            .with_columns(
                brk=(pl.col("st").shift(-1) - pl.col("end").add(1)).abs().le(8000)
            )
            .fill_null(True)
            .with_columns(pl.col("brk").rle_id())
            .with_columns(
                pl.when(pl.col("brk") % 2 == 0)
                .then(pl.col("brk") - 1)
                .otherwise(pl.col("brk"))
            )
            .partition_by("brk")
        )
        itree_methyl = IntervalTree(
            Interval(row["st"], row["end"], row["avg"])
            for row in df_chr_methyl.iter_rows(named=True)
        )
        # # Use std of entire region above mean.
        chr_methyl_avg_std = df_chr_methyl["avg"].mean() ** 0.5

        for df_grp in df_chr_methyl_adj_groups:
            grp_st, grp_end = df_grp["st"].first(), df_grp["end"].max()

            # Calculate significant regions within each group using Mahalanobis Distance.
            # This only finds the major peaks.
            # We use the z-scores method below to find entire regions before intersecting them with these intervals.
            outlier_itvs = detect_outlier_regions(
                df_grp,
                cols=["prom", "avg"],
                p_val=thr_pvalue,
            )
            print(
                f"Detected {len(outlier_itvs)} outlier intervals at p-value of {thr_pvalue} for {chrom}:{grp_st}-{grp_end}.",
                file=sys.stderr,
            )

            # Calculate z-score across subregion.
            # A z-score threshold will find bins with abnormally low methylation.
            # - We assume the average methyl percentage is normally distributed.
            df_grp = df_grp.with_columns(
                zscore=(pl.col("avg") - pl.col("avg").mean()) / chr_methyl_avg_std
            )

            # Then filter for expected negative zscore.
            df_peaks = df_grp.filter(pl.col("zscore") < thr_zscore)

            print(
                f"Filtered {df_grp.shape[0] - df_peaks.shape[0]} dips below {thr_zscore=} for {chrom}:{grp_st}-{grp_end}.",
                file=sys.stderr,
            )

            cdrs = [
                Interval(row["st"], row["end"] + 1, (row["avg"], row["prom"]))
                for row in df_peaks.iter_rows(named=True)
            ]

            def check_in_between(itv_1: Interval, itv_2: Interval) -> bool:
                df_between = df_grp.filter(
                    (pl.col("st") >= itv_1.end) & (pl.col("end") <= itv_2.begin)
                )
                mean_between = df_between["avg"].mean()
                if not mean_between:
                    return True
                return mean_between <= df_grp["avg"].mean()

            # Extend CDRs making sure that we don't extend past the mean between two intervals.
            merged_cdrs = merge_itvs(
                cdrs,
                dst=bp_merge,
                fn_cmp=check_in_between,
                fn_merge_itv=lambda x, y: (
                    Interval(
                        x.begin,
                        y.end,
                        ((x.data[0] + y.data[0]) / 2, max(x.data[1], y.data[1])),
                    )
                ),
            )

            num_merged = len(cdrs) - len(merged_cdrs)
            if num_merged:
                print(
                    f"Merged {len(cdrs) - len(merged_cdrs)} CDRs for {chrom}:{grp_st}-{grp_end}.",
                    file=sys.stderr,
                )

            for cdr in sorted(merged_cdrs):
                if not outlier_itvs.overlap(cdr):
                    print(
                        f"Omitted {chrom}:{cdr.begin}-{cdr.end} not overlapping significant interval.",
                        file=sys.stderr,
                    )
                    continue
                cdr_st, cdr_end, (cdr_avg, cdr_prom) = cdr.begin, cdr.end, cdr.data

                edge_len = 1_000_000
                ovl_left_edge = itree_methyl.overlap(cdr_st - edge_len, cdr_st)
                ovl_right_edge = itree_methyl.overlap(cdr_end, cdr_end + edge_len)
                left_edge_median = np.median([itv.data for itv in ovl_left_edge])
                right_edge_median = np.median([itv.data for itv in ovl_right_edge])
                left_edge_len = sum(itv.length() for itv in ovl_left_edge)
                right_edge_len = sum(itv.length() for itv in ovl_right_edge)

                print(
                    cdr,
                    [
                        left_edge_median * (left_edge_len / edge_len),
                        right_edge_median * (right_edge_len / edge_len),
                    ],
                )
                length = cdr.length()

                # # CDR should have prominence greater than average dip height.
                # if cdr_avg > cdr_prom:
                #     continue

                if length < thr_cdr:
                    print(
                        f"Omitted {chrom}:{cdr.begin}-{cdr.end} with {length=} below {thr_cdr=}.",
                        file=sys.stderr,
                    )
                    continue
                row = [
                    chrom,
                    str(cdr_st),
                    str(cdr_end),
                    "cdr",
                    # str(cdr_avg),
                    f"{cdr_avg},{cdr_prom}",
                    "+",
                    str(cdr_st),
                    str(cdr_end),
                    "255,0,0",
                ]
                args.outfile.write("\t".join(row) + "\n")


if __name__ == "__main__":
    raise SystemExit(main())
