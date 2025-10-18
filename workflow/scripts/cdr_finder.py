import os
import sys
import math
import json
import argparse
import statistics
from collections import defaultdict
from typing import Iterable

import matplotlib.axes
import polars as pl
import matplotlib.pyplot as plt
import matplotlib.patheffects as PathEffects

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



def group_by_dst(df: pl.DataFrame, dst: int, group_name: str) -> pl.DataFrame:
    try:
        df = df.drop("index")
    except pl.exceptions.ColumnNotFoundError:
        pass
    return (
        df.with_columns(
            # c1  st1 (end1)
            # c1 (st2) end2
            dst_behind=(pl.col("st") - pl.col("end").shift(1)).fill_null(0),
            dst_ahead=(pl.col("st").shift(-1) - pl.col("end")).fill_null(0),
        )
        .with_row_index()
        .with_columns(
            **{
                # Group HOR units based on distance.
                group_name: pl.when(pl.col("dst_behind").le(dst))
                # We assign 0 if within merge dst.
                .then(pl.lit(0))
                # Otherwise, give unique index.
                .otherwise(pl.col("index") + 1)
                # Then create run-length ID to group on.
                # Contiguous rows within distance will be grouped together.
                .rle_id()
            },
        )
        .with_columns(
            # Adjust groups in scenarios where should join group ahead or behind but given unique group.
            # B:64617 A:52416 G:1
            # B:52416 A:1357  G:2 <- This should be group 3.
            # B:1357  A:1358  G:3
            pl.when(pl.col("dst_behind").le(dst) & pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name))
            .when(pl.col("dst_behind").le(dst))
            .then(pl.col(group_name).shift(1))
            .when(pl.col("dst_ahead").le(dst))
            .then(pl.col(group_name).shift(-1))
            .otherwise(pl.col(group_name))
        )
    )


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
        help="Bases to look on both edges of cdr to determine relative height.",
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
        "--baseline_avg_methyl",
        type=float,
        default=0.4,
        help=" ".join(
            [
                "Average methylation baseline used to scale thresholds in cases of low average methylation.",
                "If average methylation is below baseline, multiplies threshold by the ratio in methylation. (baseline / avg_methyl)",
                "Will only increase threshold.",
                "Reduces false positives when low average methylation coverage.",
            ]
        ),
    )
    ap.add_argument(
        "--override_chrom_params",
        type=str,
        default=None,
        help="JSON file with params per chromosome name to override default settings. Each parameter name should match.",
    )
    ap.add_argument(
        "--output_plot_dir",
        help="Output debug plots per chrom.",
        type=str,
        default=None,
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

    output_plot_dir = args.output_plot_dir
    if output_plot_dir:
        os.makedirs(output_plot_dir, exist_ok=True)

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
        baseline_avg_methyl = override_params.get(chrom, {}).get(
            "baseline_avg_methyl", args.baseline_avg_methyl * 100
        )

        # Group adjacent, contiguous intervals.
        df_chr_methyl_adj_groups = group_by_dst(df_chr_methyl, 1, "brk").partition_by("brk")
        avg_methyl_median = df_chr_methyl["avg"].median()
        avg_methyl_mean = df_chr_methyl["avg"].mean()
        # Threshold scaling factor. Will always be at least 1.
        # Reduces false positives when mean avg methylation is low.
        thr_scaling_factor = max(
            baseline_avg_methyl / avg_methyl_mean if baseline_avg_methyl != 0.0 else 1, 1
        )
        avg_methyl_std = df_chr_methyl["avg"].std()
        cdr_prom_thr = (
            (avg_methyl_median * thr_prom_perc_valley) * thr_scaling_factor
            if thr_prom_perc_valley
            else None
        )
        cdr_height_thr = (
            avg_methyl_median * thr_height_perc_valley
        ) * thr_scaling_factor
        print(
            f"Using CDR height threshold of {cdr_height_thr} and prominence threshold of {cdr_prom_thr} for {chrom}.",
            file=sys.stderr,
        )

        # Find peaks within the signal per group.
        for df_chr_methyl_adj_grp in df_chr_methyl_adj_groups:
            df_chr_methyl_adj_grp = df_chr_methyl_adj_grp.drop("index").with_row_index()

            # Require valley has prominence of some percentage of median methyl signal.
            # Invert for peaks.
            # wlen applies to both edges (+/- n/2 indices) to find prominence. Lower prevents overestimate of peak prom.
            window = (df_chr_methyl_adj_grp["end"] - df_chr_methyl_adj_grp["st"]).mode()[0]
            wlen = (bp_edge // window) * 2
            _, peak_info = signal.find_peaks(
                -df_chr_methyl_adj_grp["avg"], width=1, prominence=cdr_prom_thr, wlen=wlen
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
                if not cdr_low:
                    continue

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

                cdr_intervals[chrom].chop(cdr_st, cdr_end)

                # Add merge distance bp.
                if bp_merge:
                    cdr_st = cdr_st - bp_merge
                    cdr_end = cdr_end + bp_merge
                # Trim overlaps.
                itv_cdr = Interval(cdr_st, cdr_end, cdr_height)
                cdr_intervals[chrom].add(itv_cdr)

        if output_plot_dir:
            itvs_cdr = cdr_intervals[chrom]
            fig, ax = plt.subplots(figsize=(16, 4), layout="constrained")
            ax: matplotlib.axes.Axes
            ax.fill_between(
                x=df_chr_methyl["st"],
                y1=df_chr_methyl["avg"],
                color="black",
            )
            ax.set_xlim(
                left=df_chr_methyl["st"].min(),
                right=df_chr_methyl["st"].max(),
            )
            lines = [
                (avg_methyl_median, "Median", "black"),
                (cdr_prom_thr, "Prominence threshold", "red"),
                (cdr_height_thr, "Height threshold", "blue"),
            ]
            _, xmax = ax.get_xlim()
            for value, label, color in lines:
                value = round(value)
                ax.axhline(round(value), label=label, linestyle="dotted", color=color)
                txt = ax.text(
                    xmax,
                    value,
                    str(value),
                    ha="center",
                    va="center",
                    color=color,
                    fontsize="small",
                )
                txt.set_path_effects(
                    [PathEffects.withStroke(linewidth=2, foreground="w")]
                )

            ax.set_xlabel("Position (bp)")
            ax.set_ylabel("Average CpG methylation (%)")
            ax.set_ylim(0.0, 100.0)
            for cdr in itvs_cdr.iter():
                midpt = cdr.begin + ((cdr.end - cdr.begin) / 2)
                # https://osxastrotricks.wordpress.com/2014/12/02/add-border-around-text-with-matplotlib/
                txt = ax.text(
                    midpt,
                    avg_methyl_median + 5,
                    round(cdr.data),
                    color="black",
                    fontsize="small",
                    ha="center",
                    va="center",
                )
                txt.set_path_effects(
                    [PathEffects.withStroke(linewidth=2, foreground="w")]
                )
                ax.axvspan(cdr.begin, cdr.end, color="red", alpha=0.5, label="CDR")

            handles, labels = ax.get_legend_handles_labels()
            labels_handles = dict(zip(labels, handles))

            for spine in ["top", "right"]:
                ax.spines[spine].set_visible(False)

            # No scientific notation
            ax.ticklabel_format(useOffset=False, style="plain")
            ax.legend(
                labels=labels_handles.keys(),
                handles=labels_handles.values(),
                loc="center left",
                bbox_to_anchor=(1.0, 1.0),
            )

            output_plot = os.path.join(output_plot_dir, f"{chrom}.png")
            fig.savefig(output_plot, bbox_inches="tight", dpi=600)
            plt.close()

    # Merge overlaps and output.
    for chrom, cdrs in cdr_intervals.items():
        if bp_merge:
            starting_intervals = len(cdrs)
            cdrs.merge_overlaps(data_reducer=lambda x, y: max(x, y))
            print(
                f"Merged {starting_intervals - len(cdrs)} intervals in {chrom}.",
                file=sys.stderr,
            )

        for cdr in sorted(cdrs.iter()):
            cdr_st, cdr_end, cdr_height = cdr.begin, cdr.end, cdr.data
            if bp_merge:
                cdr_st += bp_merge
                cdr_end -= bp_merge
            row = [
                chrom,
                cdr_st,
                cdr_end,
                "cdr",
                cdr_height,
                ".",
                cdr_st,
                cdr_end,
                "0,0,0",
            ]
            args.outfile.write(f"{'\t'.join(str(e) for e in row)}\n")


if __name__ == "__main__":
    raise SystemExit(main())
