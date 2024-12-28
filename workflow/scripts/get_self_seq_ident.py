import os
import sys
import pysam
import argparse
import itertools

from statistics import mean
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor

from estimate_identity import convertMatrixToBed, createSelfMatrix
from read_fasta import generateKmersFromFasta


def get_self_seq_ident(
    seq_id: str,
    seq: str,
    outdir: str,
    window: int,
    delta: float,
    kmer_size: int,
    ident_thr: float,
    modimizer: int,
    n_bins: int,
    ignore_bands: int,
) -> None:
    print(f"Generating self sequence identity for {seq_id}.", file=sys.stderr)
    kmers = [kmer_hash for kmer_hash in generateKmersFromFasta(seq, kmer_size)]
    mtx = createSelfMatrix(kmers, window, delta, kmer_size, ident_thr, False, modimizer)
    bed = convertMatrixToBed(mtx, window, ident_thr, seq_id, seq_id, True)

    print(
        f"Converting 2D self sequence identity matrix to 1D for {seq_id}.",
        file=sys.stderr,
    )
    # Use dictionary to avoid sparse mtx.
    aln_mtx = defaultdict(dict)
    # Use islice to avoid copying.
    for line in itertools.islice(bed, 1, len(bed) + 1):
        q_ctg, q_st, q_end, r_ctg, r_st, r_end, perID = line
        # Convert position to indices.
        x = q_st // window
        y = r_st // window
        aln_mtx[x][y] = perID

    outfile = os.path.join(outdir, f"{seq_id}.bed")
    st_idxs = list(aln_mtx.keys())

    print(
        f"Writing 1D self sequence identity array for {seq_id} to {outfile}",
        file=sys.stderr,
    )
    with open(outfile, "wt") as fh:
        for st_idx in st_idxs:
            st = st_idx * window + 1
            end = st + window - 1
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
            fh.write(f"{seq_id}\t{st}\t{end}\t{mean_ident}\n")

    return None


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-i",
        "--infile",
        required=True,
        type=str,
        help="Input fasta.",
    )
    ap.add_argument(
        "-o",
        "--outdir",
        type=str,
        required=True,
        help="Output directory for 1D self-identity alignment bedfile by contig.",
    )
    ap.add_argument(
        "-p", "--processes", default=4, type=int, help="Number of processes."
    )
    ap.add_argument(
        "-t", "--ident_thr", default=0.86, type=float, help="Identity threshold."
    )
    ap.add_argument("-w", "--window", default=5000, type=int, help="Window size.")
    ap.add_argument("-k", "--kmer_size", default=21, type=int, help="K-mer size.")
    ap.add_argument(
        "-d",
        "--delta",
        default=0.5,
        type=float,
        help="Fraction of neighboring partition to include in identity estimation. Must be between 0 and 1, use > 0.5 is not recommended.",
    )
    ap.add_argument(
        "-m",
        "--modimizer",
        default=1000,
        type=int,
        help="Modimizer sketch size. A lower value will reduce the number of modimizers, but will increase performance. Must be less than --window.",
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
    ident_thr = args.ident_thr
    kmer_size = args.kmer_size
    window = args.window
    delta = args.delta
    modimizer = args.modimizer
    n_bins = args.n_bins
    ignore_bands = args.ignore_bands
    processes = args.processes
    outdir = args.outdir

    os.makedirs(outdir, exist_ok=True)

    seq = pysam.FastaFile(args.infile)
    with ProcessPoolExecutor(max_workers=processes) as pool:
        _ = pool.map(
            get_self_seq_ident,
            *zip(
                *[
                    (
                        seq_id,
                        seq.fetch(seq_id),
                        outdir,
                        window,
                        delta,
                        kmer_size,
                        ident_thr,
                        modimizer,
                        n_bins,
                        ignore_bands,
                    )
                    for seq_id in seq.references
                ]
            ),
        )


if __name__ == "__main__":
    raise SystemExit(main())
