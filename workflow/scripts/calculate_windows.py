import sys
import math
import argparse
import statistics
import itertools
import multiprocessing as mp

from collections import defaultdict
from typing import Generator
from intervaltree import IntervalTree, Interval


METH_BED_COLS = [
    "chrom",
    "start",
    "end",
    "mod_base_name",
    "score",
    "strand",
    "n1",
    "n2",
    "n3",
    "coverage",
    "freq",
    "nmod",
    "ncanon",
    "naltmod",
    "ndel",
    "nfail",
    "ndiff",
    "nnocall",
]


def bin_freq(
    chrom: str, start: int, stop: int, meth_intervals: IntervalTree, window_size: int
) -> list[tuple[str, int, int, float, int]]:
    sys.stderr.write(
        f"Calculating average methylation frequency for {chrom}:{start}-{stop}\n"
    )

    interval_overlap = IntervalTree(meth_intervals.overlap(start, stop))

    # Removed slide and just kept window size as effectively the same.
    # ex. seq=[0:10_001] and slide = 5000 and window size = 10_000 -> window=(0, 10_000)
    # Why have slide when can just have window size?
    # Window size only altered index so would produce unexpected behavior.
    # ex. seq=[0:10_001] and slide = 2000 and window size = 5000 -> window=[0]
    window_ends = tuple(range(start, stop + window_size, window_size))
    window_intervals = [
        Interval(window_ends[i], window_ends[i + 1])
        for i in range(len(window_ends) - 1)
    ]

    averaged_intervals = []
    for window in window_intervals:
        window_overlap: set[Interval] = interval_overlap.overlap(window)

        if len(window_overlap) == 0:
            continue

        # Use absolute start/stop
        bin_start, bin_stop = window.begin, window.end
        perc_methyl, coverage = zip(*(i.data for i in window_overlap))
        avg_meth = round(statistics.mean(perc_methyl), 2)
        avg_cov = round(statistics.mean(coverage), 2)

        # Removed index.
        averaged_intervals.append((chrom, bin_start, bin_stop, avg_meth, avg_cov))

    return averaged_intervals


def average_methyl_freq_windows(
    target_bed: str, methylation_tsv: str, processes: int, window_size: int
) -> Generator[tuple[str, int, int, float, int], None, None]:
    """
    Average methylation frequency across a given window size and target region.
    """
    with open(target_bed, "rt") as bed_fh:
        target_regions = []
        for line in bed_fh.readlines():
            chrom, start, stop, *_ = line.strip().split("\t")
            start, stop = int(start), int(stop)
            target_regions.append((chrom, start, stop))

    meth_intervals = defaultdict(IntervalTree)
    with open(methylation_tsv, "rt") as meth_fh:
        for i, line in enumerate(meth_fh):
            line = line.strip().split("\t")
            try:
                line_info = dict(zip(METH_BED_COLS, line))

                # Check for nans and set to zero.
                freq = float(line_info["freq"])
                if math.isnan(freq):
                    freq = 0.0

                interval = Interval(
                    int(line_info["start"]),
                    int(line_info["end"]),
                    (freq, int(line_info["coverage"])),
                )
                meth_intervals[line_info["chrom"]].add(interval)
            except ValueError:
                sys.stderr.write(f"Line {i} in {methylation_tsv} is invalid: {line}\n")
                continue

    sys.stderr.write(
        f"Using {processes} process(es) for {len(target_regions)} regions.\n"
    )
    if processes == 1:
        results = []
        for chrom, start, stop in target_regions:
            results.append(
                bin_freq(chrom, start, stop, meth_intervals[chrom], window_size)
            )
    else:
        with mp.Pool(processes) as pool:
            results = pool.starmap(
                bin_freq,
                [
                    (chrom, start, stop, meth_intervals[chrom], window_size)
                    for chrom, start, stop in target_regions
                ],
            )

    yield from itertools.chain(*results)


def main() -> None:
    ap = argparse.ArgumentParser(description="")
    ap.add_argument(
        "-b",
        "--target_bed",
        required=True,
        type=str,
        help="Regions to average methylation over.",
    )
    ap.add_argument(
        "-m",
        "--methylation_tsv",
        required=True,
        type=str,
        help="Methylation bed from modbam2bed.",
    )
    ap.add_argument(
        "-o",
        "--binned_freq",
        help="Output binned frequency.",
        default=sys.stdout,
        type=argparse.FileType("wt"),
    )
    ap.add_argument(
        "-p",
        "--processes",
        help="Number of processes to run each row in --target_bed.",
        default=1,
        type=int,
    )
    ap.add_argument(
        "--window_size",
        help="Window size to bin and average methylation frequency.",
        type=int,
        default=5_000,
    )

    args = ap.parse_args()

    for chrom, start, stop, avg_meth, avg_cov in average_methyl_freq_windows(
        args.target_bed, args.methylation_tsv, args.processes, args.window_size
    ):
        args.binned_freq.write(f"{chrom}\t{start}\t{stop}\t{avg_meth}\t{avg_cov}\n")


if __name__ == "__main__":
    raise SystemExit(main())
