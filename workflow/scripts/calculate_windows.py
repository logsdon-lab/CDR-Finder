import sys
import argparse
import pandas as pd
import numpy as np


def getOverlap(alin_rec):
    a = [alin_rec["start"], alin_rec["end"]]
    b = [alin_rec["bin_start"], alin_rec["bin_stop"]]
    return max(0, min(a[1], b[1]) - max(a[0], b[0]))


def main():
    ap = argparse.ArgumentParser(description="")
    ap.add_argument("-b", "--target_bed", required=True, type=str, help="")
    ap.add_argument("-m", "--methylation_tsv", required=True, type=str, help="")
    ap.add_argument("-o", "--binned_freq", help="", default=sys.stdout, type=argparse.FileType("wt"))
    ap.add_argument("--window_size", help="", type=int, default=5_000)
    ap.add_argument("--slide", help="", type=int, default=5_000)

    args = ap.parse_args()

    bed_df = pd.read_csv(
        args.target_bed, header=None, sep="\t", names=["chrom", "start", "stop"]
    )
    meth_df = pd.read_csv(
        args.methylation_tsv,
        header=None,
        sep="\t",
        names=[
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
            "ncanon",
            "nmod",
            "nfilt",
            "nnocall",
            "naltmod",
        ],
    )
    keep_window = pd.DataFrame()

    for index in bed_df.index:
        start = bed_df.at[index, "start"]
        stop = bed_df.at[index, "stop"]
        chrom = bed_df.at[index, "chrom"]
        chrom_meth_df = meth_df.loc[meth_df["chrom"] == chrom].copy()
        seq = np.arange(
            int(start), int(stop) + int(args.slide), int(args.slide)
        )
        window_size = (int(int(args.window_size) / int(args.slide))) + 1
        window_avg = pd.DataFrame(columns=["Chr", "Bin", "Freq"])

        for i in range(len(seq) - window_size + 1):
            window = seq[i : i + window_size]
            window_meth_pre = chrom_meth_df.copy()
            RM_start = window[0] - start
            RM_stop = window[-1] - start
            window_meth_pre["bin_start"] = window[0]
            window_meth_pre["bin_stop"] = window[-1]
            print(RM_start, RM_stop)
            if len(window_meth_pre) == 0:
                continue
            window_meth_pre["overlap"] = window_meth_pre.apply(getOverlap, axis=1)
            window_meth = window_meth_pre.loc[window_meth_pre["overlap"] != 0]
            if len(window_meth) == 0:
                continue
            avg_meth = window_meth["freq"].mean()
            new_row = pd.DataFrame({
                    "Chr": [str(chrom)],
                    "Bin": [str(RM_start) + "-" + str(RM_stop)],
                    "Freq": [avg_meth]
            })
            window_avg = pd.concat([window_avg, new_row], ignore_index=True)

        window_avg = window_avg.round({"Freq": 2})
        keep_window = pd.concat([keep_window, window_avg])

    keep_window.reset_index(inplace=True)

    with open(args.binned_freq, "w+") as outfile:
        for idx_2 in keep_window.index:
            pos1 = str(keep_window.at[idx_2, "Bin"]).split("-")[0]
            pos2 = str(keep_window.at[idx_2, "Bin"]).split("-")[1]
            freq = str(keep_window.at[idx_2, "Freq"])
            chrom = str(keep_window.at[idx_2, "Chr"])
            outfile.write("\t".join([chrom, pos1, pos2, freq, f"{idx_2}\n"]))


def __main__():
    raise SystemExit(main())
