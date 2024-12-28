import sys
import csv
import argparse
from statistics import mean
from typing import TextIO
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=argparse.FileType("rt"),
        help="Input self-identity alignment bedfile.",
    )
    ap.add_argument(
        "-o",
        "--outfile",
        default=sys.stdout,
        type=argparse.FileType("wt"),
        help="Input self-identity alignment bedfile.",
    )
    ap.add_argument(
        "-b",
        "--n_bins",
        default=5,
        type=int,
        help="Number of bins to calculate average sequence identity over.",
    )
    ap.add_argument(
        "--ignore_bands",
        default=2,
        type=int,
        help="Number of bands ignored along self-identity diagonal.",
    )
    args = ap.parse_args()

    infile: TextIO = args.infile
    outfile: TextIO = args.outfile
    n_bins = args.n_bins
    ignore_bands = args.ignore_bands

    header = next(infile).strip()
    reader = csv.DictReader(infile, delimiter="\t", fieldnames=header.split("\t"))
    window_size = None
    ctg = None

    # Use dictionary to avoid sparse mtx.
    aln_mtx = defaultdict(dict)
    for row in reader:
        if not window_size:
            window_size = int(row["query_end"]) - int(row["query_start"])
        if not ctg:
            ctg = row["#query_name"]
        # Convert position to indices.
        x = int(row["query_start"]) // window_size
        y = int(row["reference_start"]) // window_size
        aln_mtx[x][y] = float(row["perID_by_events"])

    st_idxs = list(aln_mtx.keys())
    for st_idx in st_idxs:
        st = st_idx * window_size + 1
        end = st + window_size - 1
        band_end_idx = st_idx + n_bins
        # Within the alignment matrix with a n_bins of 5 and ignore_bands of 2:
        # - '*' is the calculated aln band
        # - '+' is self aln.
        # 4 * * *   +
        # 3 * *   +
        # 2 *   +
        # 1   +
        # 0 +
        #   0 1 2 3 4
        mean_ident = mean(
            aln_mtx[x].get(y, 0.0)
            for x in range(st_idx, band_end_idx)
            for y in range(x + ignore_bands, band_end_idx)
        )
        print(f"{ctg}\t{st}\t{end}\t{mean_ident}", file=outfile)


if __name__ == "__main__":
    raise SystemExit(main())
