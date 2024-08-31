import polars as pl
from scipy import signal


def main():
    df = pl.read_csv(
        "../../test/cdr/input/HG00731_intersect.bed",
        separator="\t",
        has_header=False,
        new_columns=["chr", "st", "end", "avg", "cov", "idx"],
    )

    for chr, df_chr_methyl in df.group_by(["chr"]):
        chr = chr[0]
        methyl_signal = df_chr_methyl["avg"]
        smoothed_methyl_signal = signal.medfilt(methyl_signal, 3)
        valley_prom = smoothed_methyl_signal.mean() * 0.33
        peaks, peak_info = signal.find_peaks(
            -smoothed_methyl_signal, width=1, prominence=valley_prom
        )

        # plt.plot(methyl_signal)
        # ax = plt.gca()
        # for l, r in zip(peak_info["left_ips"], peak_info["right_ips"]):
        #     ax.axvspan(l, r, color="red", alpha=0.5)
        # plt.savefig(f"{chr}.png")
        # plt.close()


if __name__ == "__main__":
    raise SystemExit(main())
