import os
import sys
import json
import argparse
from collections import deque
from typing import TYPE_CHECKING, Callable, Iterable, Any

import polars as pl
import matplotlib.pyplot as plt
from matplotlib.axes import Axes
from intervaltree import Interval, IntervalTree

if TYPE_CHECKING:
    SubArgumentParser = argparse._SubParsersAction[argparse.ArgumentParser]
else:
    SubArgumentParser = Any


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


def add_detect_cli(subparser: SubArgumentParser) -> None:
    ap = subparser.add_parser("detect", description="Detect CDRs.")
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
        "-b", "--bp_merge", type=int, default=None, help="Base pairs to merge CDRs."
    )
    ap.add_argument(
        "-g",
        "--bp_group",
        type=int,
        default=8_000,
        help="Base pairs to group methylated regions to evaluate CDRs.",
    )
    ap.add_argument(
        "-m",
        "--thr_methyl_score",
        type=float,
        default=-3.0,
        help="Number of adjusted z-scores to find dip regions. ",
    )
    ap.add_argument(
        "-d",
        "--thr_dmethyl_score",
        type=float,
        default=1.0,
        help="Number of z-scores of the derivative of the average to find steep edges of CDRs. Higher values find steeper regions.",
    )
    ap.add_argument(
        "-l",
        "--thr_plen_score",
        type=float,
        default=1.0,
        help="Number of z-scores of the CDR length. Higher filters for larger CDRs.",
    )
    ap.add_argument("--thr_cdr", type=int, default=1, help="Minimum CDR size.")
    ap.add_argument(
        "-p",
        "--output_plot_dir",
        type=str,
        default=None,
        help="Output debug plots to directory.",
    )
    ap.add_argument(
        "--override_chrom_params",
        type=str,
        default=None,
        help="JSON file with params per chromosome name to override default settings. Each parameter name should match.",
    )


def add_plot_hist_cli(subparser: SubArgumentParser) -> None:
    ap = subparser.add_parser(
        "hist",
        description="Plot histogram of data for both edges and methyl average. Use to set proper scores for detect command.",
    )
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
        "--output_dir",
        required=True,
        type=str,
        help="Output directory for plots.",
    )


def main():
    ap = argparse.ArgumentParser(description="CDR finder.")
    sub_ap = ap.add_subparsers(dest="cmd")
    add_detect_cli(sub_ap)
    add_plot_hist_cli(sub_ap)

    args = ap.parse_args()

    if args.cmd == "detect":
        return detect(args)
    elif args.cmd == "hist":
        return plot_hist(args)
    else:
        raise ValueError(f"Not a valid command ({args.cmd})")


def add_metric_ndev_lines(
    ax: Axes, ax_top: Axes, middle: float, ndev: float, itv_bounds: Interval
):
    line_kwargs = {"color": "black", "linestyle": "dotted"}
    valid_lines = [(middle, 0)]
    # Draw middle line.
    ax.axvline(middle, **line_kwargs)

    # Draw remaining lines if within range.
    n = 1
    while True:
        lower = middle - (ndev * n)
        upper = middle + (ndev * n)
        lower_valid = itv_bounds.contains_point(lower)
        upper_valid = itv_bounds.contains_point(upper)

        if not lower_valid and not upper_valid:
            break
        if lower_valid:
            ax.axvline(lower, **line_kwargs)
            valid_lines.append((lower, f"-{n}"))
        if upper_valid:
            ax.axvline(upper, **line_kwargs)
            valid_lines.append((upper, f"+{n}"))
        n += 1
    # Clear labels.
    ax_top.set_xticks([], [])
    xticklocs, xticklbls = zip(*valid_lines)
    ax_top.set_xticks(xticklocs, xticklbls)
    ax_top.set_xlim(ax.get_xlim())


def plot_hist_avg(df_chr_methyl: pl.DataFrame, chrom: str, outfile_avg: str):
    itv_avg = Interval(df_chr_methyl["avg"].min(), df_chr_methyl["avg"].max())
    median_avg = df_chr_methyl.get_column("avg").median()
    # https://www.ibm.com/docs/en/cognos-analytics/12.0.0?topic=terms-modified-z-score
    mstd = 1.486 * (df_chr_methyl.get_column("avg") - median_avg).abs().median()

    # Plot methyl average percent.
    fig, ax = plt.subplots()
    ax_top = ax.twiny()
    ax.hist(df_chr_methyl["avg"])
    ax.set_title(f"Distribution of mean CpG methylation ({chrom})")
    ax.set_ylabel("Count")
    ax.set_xlabel("Mean CpG methylation (%)")

    add_metric_ndev_lines(ax, ax_top, median_avg, mstd, itv_avg)
    ax_top.set_xlabel("Number of modified z-scores")

    fig.savefig(outfile_avg, dpi=600, bbox_inches="tight")
    fig.clear()


def plot_hist(args: argparse.Namespace):
    os.makedirs(args.output_dir, exist_ok=True)

    df = pl.read_csv(
        args.infile,
        separator="\t",
        has_header=False,
        new_columns=["chrom", "st", "end", "avg"],
    )

    for chrom, df_chr_methyl in df.group_by(["chrom"]):
        chrom: str = chrom[0]

        print(f"Plotting {chrom}.", file=sys.stderr)

        outfile_avg = os.path.join(args.output_dir, f"{chrom}_avg.png")

        plot_hist_avg(df_chr_methyl, chrom, outfile_avg)


def group_df(df: pl.DataFrame, dst: int) -> pl.DataFrame:
    return (
        df.with_columns(
            # c1  st1 (end1)
            # c1 (st2) end2
            dst_behind=(pl.col("st") - pl.col("end").shift(1)).fill_null(0),
            dst_ahead=(pl.col("st").shift(-1) - pl.col("end")).fill_null(0),
        )
        .with_row_index()
        .with_columns(
            # Group rows based on distance.
            grp=pl.when(pl.col("dst_behind").le(dst))
            # We assign 0 if within merge dst.
            .then(pl.lit(0))
            # Otherwise, give unique index.
            .otherwise(pl.col("index") + 1)
            # Then create run-length ID to group on.
            # Contiguous rows within distance will be grouped together.
            .rle_id()
        )
        .with_columns(
            # Adjust groups in scenarios where should join group ahead or behind but given unique group.
            pl.when(pl.col("dst_behind").le(dst) & pl.col("dst_ahead").le(dst))
            .then(pl.col("grp"))
            .when(pl.col("dst_behind").le(dst))
            .then(pl.col("grp").shift(1))
            .when(pl.col("dst_ahead").le(dst))
            .then(pl.col("grp").shift(-1))
            .otherwise(pl.col("grp"))
        )
    )


def detect(args: argparse.Namespace):
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
        min_st, max_end = df_chr_methyl["st"].min(), df_chr_methyl["end"].max()
        chrom_coords = f"{chrom}:{min_st}-{max_end}"

        # Get override params.
        bp_merge = override_params.get(chrom, {}).get("bp_merge", args.bp_merge)
        if not bp_merge:
            bp_merge = 0
        bp_group = override_params.get(chrom, {}).get("bp_group", args.bp_group)

        thr_cdr = override_params.get(chrom, {}).get("thr_cdr", args.thr_cdr)

        thr_methyl_score = override_params.get(chrom, {}).get(
            "thr_methyl_score", args.thr_methyl_score
        )
        if thr_methyl_score > 0:
            raise ValueError(f"Adjusted z-score threshold for {chrom} is non-negative.")

        thr_dmethyl_score = override_params.get(chrom, {}).get(
            "thr_dmethyl_score", args.thr_dmethyl_score
        )

        thr_plen_score = override_params.get(chrom, {}).get(
            "thr_plen_score", args.thr_plen_score
        )

        if args.output_plot_dir:
            os.makedirs(args.output_plot_dir, exist_ok=True)

            fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(12, 6))
            fig.suptitle(chrom_coords)
            fig.supxlabel("Position")

            ax_1: Axes = axes[0]
            ax_2: Axes = axes[1]
            ax_3: Axes = axes[2]
            labels = [
                "Mean\nCpG methylation (%)",
                "Mean\nCpG methylation\nz-scores",
                "Change in mean\nCpG methylation\nz-scores",
            ]
            for ax, lbl in zip(axes, labels):
                ax: Axes
                ax.set_xlim(min_st, max_end)
                ax.set_ylabel(lbl)

            # Show median.
            ax_1.axhline(df_chr_methyl["avg"].median(), linestyle="dotted")

        # Group adjacent, contiguous intervals.
        df_chr_methyl = group_df(df_chr_methyl, bp_group)

        grp_lowest_methyl = (
            df_chr_methyl.filter(pl.col("avg") != 0)
            .group_by(["grp"])
            .agg(min_dip=pl.col("avg").min())
            .filter(pl.col("min_dip") == pl.col("min_dip").min())
            .row(0, named=True)
            .get("grp")
        )
        df_grp = df_chr_methyl.filter(pl.col("grp") == grp_lowest_methyl)
        st_grp, end_grp = df_grp["st"].min(), df_grp["end"].max()

        grp_chrom = f"{chrom}:{st_grp}-{end_grp}"

        # Calculate MAD-score across subregion.
        # - A MAD-score threshold will find bins with abnormally low methylation but will resist outliers.
        mad = (df_grp["avg"] - df_grp["avg"].median()).abs().median()

        df_edges = (
            df_grp.with_columns(
                # First derivative
                d_avg=pl.col("avg").diff() / 2.0
            )
            # Group derivatives in one direction and within distance
            .with_columns(d_grp=pl.col("d_avg").lt(0).rle_id())
            .group_by(["d_grp"])
            # Calculate derivative across all oriented derivatives. Uses sum rule.
            # https://www.geogebra.org/m/wtfhdbf3
            # And get intervals.
            .agg(
                st=pl.col("st").min(),
                end=pl.col("end").max() + 1,
                d_avg=pl.col("d_avg").sum(),
            )
            .with_columns(
                # Second derivative.
                pl.col("d_avg").diff() / 2.0
            )
            .drop_nulls()
            # Calculate z-score
            .with_columns(
                zscore_d_avg=(pl.col("d_avg") - pl.col("d_avg").mean())
                / pl.col("d_avg").std()
            )
        )
        # Add all edges found to intervaltree to query.
        itree_edges = IntervalTree(
            [
                Interval(row["st"], row["end"], row["zscore_d_avg"])
                for row in df_edges.iter_rows(named=True)
            ]
        )

        if args.output_plot_dir:
            for itv in itree_edges:
                itv: Interval
                ax_3.bar(
                    x=itv.begin,
                    height=itv.data if itv.data else 0.0,
                    width=itv.length(),
                    color="blue",
                    alpha=0.5,
                    align="edge",
                )

        df_grp = df_grp.with_columns(
            # https://www.statology.org/modified-z-score/
            mscore_avg=((pl.col("avg") - pl.col("avg").median()) * 0.6745) / mad,
        )
        # Then filter for expected negative mscore and peak length score.
        df_peaks = (
            df_grp.with_columns(
                is_dip=pl.when(pl.col("mscore_avg") < thr_methyl_score)
                .then(pl.col("index").max() + 1)
                .otherwise(pl.col("index"))
                .rle_id()
            )
            .group_by(["is_dip"])
            .agg(
                pl.col("st").min(),
                pl.col("end").max(),
                pl.col("mscore_avg").median(),
            )
            .with_columns(plen=pl.col("end") - pl.col("st"))
            .with_columns(
                zscore_plen=(pl.col("plen") - pl.col("plen").mean())
                / pl.col("plen").std(),
            )
            .filter(pl.col("mscore_avg") < thr_methyl_score)
        )

        print(
            f"Filtered {df_grp.shape[0] - df_peaks.shape[0]} dips below {thr_methyl_score=} and above {thr_plen_score=} for {grp_chrom}.",
            file=sys.stderr,
        )

        cdrs = [
            Interval(row["st"], row["end"] + 1, row["mscore_avg"])
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
                    (x.data + y.data) / 2,
                )
            ),
        )

        num_merged = len(cdrs) - len(merged_cdrs)
        if num_merged:
            print(
                f"Merged {len(cdrs) - len(merged_cdrs)} CDRs for {grp_chrom}.",
                file=sys.stderr,
            )

        if args.output_plot_dir:
            ax_1.bar(
                x=df_grp["st"],
                height=df_grp["avg"],
                width=df_grp["end"] - df_grp["st"],
                color="black",
                align="edge",
            )
            ax_2.bar(
                x=df_grp["st"],
                height=df_grp["mscore_avg"],
                width=df_grp["end"] - df_grp["st"],
                color="orange",
                alpha=0.5,
                align="edge",
            )

        for cdr in sorted(merged_cdrs):
            cdr_st, cdr_end, cdr_mscore_avg = cdr.begin, cdr.end, cdr.data
            cdr_length = cdr.length()

            if cdr_length < thr_cdr:
                print(
                    f"Filtered {chrom}:{cdr.begin}-{cdr.end} with {cdr_length=} below {thr_cdr=}.",
                    file=sys.stderr,
                )
                continue

            # Merge overlap edges by two windows if edges jagged.
            ovl_edges = sorted(itree_edges.overlap(cdr.begin, cdr.end))
            if not ovl_edges:
                print(
                    f"Filtered {chrom}:{cdr.begin}-{cdr.end} not overlapping 2 edge intervals.",
                    file=sys.stderr,
                )
                continue

            # # All edges sloped in the same direction is wrong (x).
            # #  x _ o _ x
            # #   / \_/ \
            # #  /       \
            # elif all((itv.data * ovl_edges[0].data) > 0 for itv in ovl_edges[1:]):
            #     print(
            #         f"Filtered {chrom}:{cdr.begin}-{cdr.end} where edges are sloped in only one direction. ({ovl_edges})",
            #         file=sys.stderr,
            #     )
            #     continue

            # Check that among valid edges at least one negative and positive value is present.
            # min \_/ max
            valid_edges = [
                itv for itv in ovl_edges if abs(itv.data) > thr_dmethyl_score
            ]
            if not valid_edges:
                print(
                    f"Filtered {chrom}:{cdr.begin}-{cdr.end} where no edge ({ovl_edges}) is below ({thr_dmethyl_score=})",
                    file=sys.stderr,
                )
                continue

            min_edge, max_edge = (
                min(valid_edges, key=lambda x: x.data),
                max(valid_edges, key=lambda x: x.data),
            )
            # Must have negative and positive edge and must cover the majority of the CDR.
            cdr_edges_valid = min_edge.data < 0 and max_edge.data > 0
            if not cdr_edges_valid:
                print(
                    f"Filtered {chrom}:{cdr.begin}-{cdr.end} missing a valid descending and ascending CDR edge. ({ovl_edges})",
                    file=sys.stderr,
                )
                continue

            row = [
                chrom,
                str(cdr_st),
                str(cdr_end),
                "cdr",
                str(cdr_mscore_avg),
                ".",
                str(cdr_st),
                str(cdr_end),
                "255,0,0",
            ]
            if args.output_plot_dir:
                ax_1.axvspan(cdr_st, cdr_end, alpha=0.34, color="red")

            args.outfile.write("\t".join(row) + "\n")

        if args.output_plot_dir:
            # Reset ticks so in 1 z-score increments.
            ymin_2, ymax_2 = ax_2.get_ylim()
            ymin_3, ymax_3 = ax_3.get_ylim()
            ymin_2, ymax_2 = int(ymin_2), int(ymax_2)
            ymin_3, ymax_3 = int(ymin_3), int(ymax_3)
            yticks_2 = (*range(ymin_2, 0), *range(0, ymax_2))
            yticks_3 = (*range(ymin_3, 0), *range(0, ymax_3))
            ax_2.set_yticks(yticks_2, [str(tk) for tk in yticks_2])
            ax_3.set_yticks(yticks_3, [str(tk) for tk in yticks_3])
            fig.savefig(
                os.path.join(args.output_plot_dir, f"{chrom_coords}.png"),
                bbox_inches="tight",
            )
            plt.close(fig)


if __name__ == "__main__":
    raise SystemExit(main())
