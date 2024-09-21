# Changes
* Added basic test files with Git LFS
* Added plot cdr R script.
* Added output, log, and benchmark dir.
* Remove singularity for rules.
* Rewrote `calculate_windows.py` as unbearably slow.
    * Use `intervaltree` library for faster overlap detection.
    * Remove `pandas` and `numpy`.
    * Remove slide argument as unnecessary and leads to undefined behavior

    ```bash
    time python workflow/scripts/calculate_windows_test.py         --target_bed test/bed/CHM13_cen_500kbp.bed         --methylation_tsv results/bed/CHM13_subset.bed         --window_size 5000      -p 1 > /dev/null
    ```
    ```
    real    0m12.428s
    user    0m10.869s
    sys     0m0.208s
    ```

    ```bash
    time python tests/scripts/calculate_windows.py         --target_bed test/bed/CHM13_cen_500kbp.bed         --methylation_tsv results/bed/CHM13_subset.bed         --window_size 5000         --slide 5000  > /dev/null
    ```
    ```
    real    22m59.106s
    user    22m34.406s
    sys     0m21.116s
    ```
* Use `modkit` instead of `modbam2bed`.
    * `modkit` produces identical results and is actively supported by ONT.
* Rework CDR detection script to not use custom valley detection algorithm and rely on pre-existing well-tested libraries. Using the [`scipy.signals`](https://docs.scipy.org/doc/scipy/reference/signal.html) library, we can achieve similar results out of the box.
    * Remove confidence.
    * Changed parameters so more robust thresholds. Rather than `low_threshold` of `39`, the only param is:
        * `thr_height_perc_valley`
            * Threshold percent of the average methylation percentage needed as the minimal height of a valley from the median.
            * ex. 0.1 allows valleys with a height from the median equal to 10% of the median.
            * *How deep should this valley be in the context of the median methylation in this region?*
        * Also, use prominence of a valley which helps in removing low-confidence CDRs.
            * See https://en.wikipedia.org/wiki/Topographic_prominence.
