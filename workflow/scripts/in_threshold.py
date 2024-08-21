import argparse
import pandas as pd


def main():

    
    # read in intersection bin bed
    intersect_bed = pd.read_csv(
        input.intersect_bed,
        header=None,
        sep="\t",
        names=["chrom", "start", "stop", "freq", "idx"],
        index_col="idx",
    )
    with open(input.target_bed) as infile, open(output.final_call, "w+") as outfile:
        # Evaluate over all windows in target bed        
        for line in infile:
            # Empty dataframe for windows
            window_avg_all = pd.DataFrame(columns=["Bin", "Freq"])
            # Parse target bed
            chrom, start = line.split("\t")[0], line.split("\t")[1]
            # Subset to chromosome
            intersect_bed_sub = intersect_bed.loc[intersect_bed["chrom"] == chrom]
            # Convert intersect bed to expected format for df
            for index in intersect_bed_sub.index:
                chrom_int = intersect_bed_sub.at[index, "chrom"]
                start_int = int(intersect_bed_sub.at[index, "start"]) + int(start)
                stop_int = int(intersect_bed_sub.at[index, "stop"]) + int(start)
                freq = intersect_bed_sub.at[index, "freq"]
                new_row2 = pd.DataFrame({
                    "Bin": [str(start_int) + "-" + str(stop_int)],
                    "pos": [start_int],
                    "end": [stop_int],
                    "Freq": [float(freq)],
                    "idx_2": [index]
                })
                window_avg_all = pd.concat([window_avg_all, new_row2], ignore_index=True)

            window_avg_all["idx_2"] = window_avg_all["idx_2"].astype(int)
            window_avg_all = window_avg_all.set_index("idx_2")
            # Establish threshold, stddev, max
            median_freq = np.median(window_avg_all["Freq"])
            stddev_freq = np.std(window_avg_all["Freq"])
            max_freq = max(window_avg_all["Freq"])
            mean_freq = np.mean(window_avg_all["Freq"])
            max_thresh = max_freq - 1 * stddev_freq
            # Subset df to those windows under the threshold for instances where mean represents high methylation (large centromeres)
            if mean_freq > LOW_THRESHOLD:
                threshold = mean_freq - 1 * stddev_freq
                window_avg_under = window_avg_all.loc[
                    window_avg_all["Freq"] <= threshold
                ].copy()
            # If methylation mean is low (small centromeres with more windows in CDR than not)
            else:
                threshold = mean_freq + 1 * (stddev_freq/2)
                window_avg_under = window_avg_all.loc[
                    ~(window_avg_all["Freq"] > threshold)
                ].copy()
            s = pd.Series(window_avg_under.index.tolist())
            # Puts windows into consecutive ranges
            groups = (
                s.groupby(s.diff().ne(1).cumsum())
                .apply(
                    lambda x: [x.iloc[0], x.iloc[-1]]
                    if len(x) >= 2
                    else [x.iloc[0]]
                )
                .tolist()
            )
            groups = [sub for sub in groups if len(sub) > 1]
            # Begin checking for content of gaps
            check_df = pd.DataFrame()
            # Make df containing only windows which are consecutive and under threshold
            for y in groups:
                check_df = pd.concat([check_df, window_avg_all.loc[y[0] : y[1]]])
            # get iterable of indexes
            initial_index_list = check_df.index
            # Begin checking for gaps in the ranges
            for i, index in enumerate(initial_index_list[:-1]):
                # If next window found,continue
                if index + 1 in initial_index_list:
                    continue
                # If not, make df containing only the windows between the gaps
                else:
                    gap_df = window_avg_all.loc[
                        initial_index_list[i] + 1 : initial_index_list[i + 1] - 1
                    ]
                    # If one window is within one stddev of the max, do not add gap
                    if (
                        len(gap_df.loc[gap_df["Freq"] > (max_freq - stddev_freq)])
                        == 0
                    ):
                        check_df = pd.concat([check_df, gap_df])
                    else:
                        continue
            s = pd.Series(sorted(check_df.index.tolist()))
            # Puts windows into consecutive ranges, again
            groups = (
                s.groupby(s.diff().ne(1).cumsum())
                .apply(
                    lambda x: [x.iloc[0], x.iloc[-1]]
                    if len(x) >= 2
                    else [x.iloc[0]]
                )
                .tolist()
            )
            groups = [sub for sub in groups if len(sub) > 1]
            # Writes output
            window_df_out = pd.DataFrame()
            for y in groups:
                pos1 = int(check_df.at[y[0], "Bin"].split("-")[0])
                pos2 = int(check_df.at[y[1], "Bin"].split("-")[1])
                # left_check and right_check are checking for flanking peaks of high methylation frequency
                left_check = window_avg_all.loc[(window_avg_all['end'] < pos1 ) & (window_avg_all['end'] >= pos1-EDGE_SEARCH )].copy()
                right_check = window_avg_all.loc[(window_avg_all['end'] > pos2 ) & (window_avg_all['end'] <= pos2+EDGE_SEARCH )].copy()
                # If neighbors are found on either side label as high confidence
                if np.max(right_check['Freq']) > max_thresh and np.max(left_check['Freq']) > max_thresh:
                    confidence = "high_confidence"
                # Else assign and label low confidence
                else:
                    confidence = "low_confidence-neighbor_peaks"
                # Label region with size threshold confidence
                if pos2 - pos1 < int(REPORT_THRESHOLD):
                    confidence += "+low_confidence-size"
                # Write to output bed file
                outfile.write(f"{chrom}\t{pos1}\t{pos2}\t{confidence}\n")
                # Add to combined DF
                window_df_out = pd.concat([window_df_out, check_df.loc[y[0]:y[1]]])
            window_df_out.to_csv(output.windows_call, sep='\t', index=False)


if __name__ == "__main__":
    raise SystemExit(main())
