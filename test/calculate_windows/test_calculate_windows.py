import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    [
        "target_bed",
        "methylation_tsv",
        "window_size",
        "exp_bed",
        # The output produced by the original script.
        "exp_og_bed",
    ],
    [
        (
            "test/calculate_windows/input/CHM13_cen_500kbp.bed",
            "test/calculate_windows/input/CHM13_methyl.bed",
            5000,
            "test/calculate_windows/expected/expected_bin_freq.bed",
            "test/calculate_windows/expected/original_bin_freq.bed",
        )
    ],
)
def test_calculate_windows(
    target_bed: str,
    methylation_tsv: str,
    window_size: int,
    exp_bed: str,
    exp_og_bed: str,
):
    with open(exp_bed) as fh, open(exp_og_bed) as ofh:
        exp_regions = sorted(
            enumerate(line.strip().split("\t") for line in fh.readlines()),
            key=lambda x: x[1],
        )
        exp_og_regions = sorted(line.strip().split("\t") for line in ofh.readlines())

    # The command output equals the expected output.
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
    )

    # The expected output equals the original output.
    for (i, exp_region), exp_og_region in zip(exp_regions, exp_og_regions):
        *coords, freq, _ = exp_region
        *ecoords, efreq, ei = exp_og_region
        freq, efreq = float(freq), float(efreq)
        assert str(i) == ei and coords == ecoords
        # Allow some approximation for floating point.
        assert freq == pytest.approx(efreq, 0.5)
