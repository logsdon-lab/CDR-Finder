import pytest
from test.helpers.integration import run_integration_test


OVERWRITE_OUTPUT: bool = False


@pytest.mark.parametrize(
    [
        "target_bed",
        "methylation_tsv",
        "window_size",
        "exp_bed",
    ],
    [
        (
            "test/calculate_windows/input/CHM13_cen_500kbp.bed",
            "test/calculate_windows/input/CHM13_methyl.bed",
            5000,
            "test/calculate_windows/expected/expected_bin_freq.bed",
        )
    ],
)
def test_calculate_windows(
    target_bed: str,
    methylation_tsv: str,
    window_size: int,
    exp_bed: str,
):
    run_integration_test(
        "python",
        "workflow/scripts/calculate_windows.py",
        "--target_bed",
        target_bed,
        "--methylation_tsv",
        methylation_tsv,
        "--window_size",
        str(window_size),
        expected_output=exp_bed,
        overwrite_output=OVERWRITE_OUTPUT,
    )
