import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["in_bed", "exp_bed", "thr_quantile", "thr_prom", "args"],
    [
        (
            "test/cdr/input/CHM13_intersect.bed",
            "test/cdr/expected/CHM13_cdr.bed",
            0.1,
            0.2,
            tuple(),
        ),
        (
            "test/cdr/input/HG00731_intersect.bed",
            "test/cdr/expected/HG00731_cdr.bed",
            0.1,
            0.2,
            tuple(),
        ),
        (
            "test/cdr/input/HG00731_intersect.bed",
            "test/cdr/expected/HG00731_cdr_merged.bed",
            0.1,
            0.2,
            tuple(["--bp_merge", str(10_000)]),
        ),
    ],
)
def test_cdr_finder(
    in_bed: str, thr_quantile: int, thr_prom: int, exp_bed: str, args: tuple[str] | None
):
    args = [] if not args else args
    run_integration_test(
        "python",
        "workflow/scripts/cdr_finder.py",
        "--infile",
        in_bed,
        "--thr_quantile_valley",
        str(thr_quantile),
        "--thr_prominence_perc_valley",
        str(thr_prom),
        *args,
        expected_output=exp_bed,
    )
