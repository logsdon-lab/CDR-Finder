# Changes
* Added basic test files with Git LFS
* Added plot cdr R script.
* Added output, log, and benchmark dir.
* Remove singularity for rules.
* Rewrote calculate_windows.py as unbearably slow.
    * Use intervaltrees for faster overlap detection.
    * Remove pandas and numpy.
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
