import pytest
from test.helpers.integration import run_integration_test


@pytest.mark.parametrize(
    ["in_bed", "exp_bed", "thr_height", "thr_prom", "args"],
    [
        (
            "test/cdr/input/CHM13_intersect.bed",
            "test/cdr/expected/CHM13_cdr.bed",
            0.5,
            0.5,
            tuple(["--bp_edge", str(500_000)]),
        ),
        (
            "test/cdr/input/HG00731_intersect.bed",
            "test/cdr/expected/HG00731_cdr.bed",
            0.5,
            0.5,
            tuple(["--bp_edge", str(500_000)]),
        ),
        (
            "test/cdr/input/HG00731_intersect.bed",
            "test/cdr/expected/HG00731_cdr_ext.bed",
            0.5,
            0.5,
            tuple(["--bp_edge", str(500_000), "--extend_edges_std", str(0)]),
        ),
        (
            "test/cdr/input/NA19331_intersect.bed",
            "test/cdr/expected/NA19331_cdr.bed",
            0.5,
            0.5,
            tuple(["--bp_edge", str(500_000)]),
        ),
    ],
)
def test_cdr_finder(
    in_bed: str,
    thr_height: float,
    thr_prom: float,
    exp_bed: str,
    args: tuple[str] | None,
):
    args = [] if not args else args
    run_integration_test(
        "python",
        "workflow/scripts/cdr_finder.py",
        "--infile",
        in_bed,
        "--thr_height_perc_valley",
        str(thr_height),
        "--thr_prom_perc_valley",
        str(thr_prom),
        *args,
        expected_output=exp_bed,
    )
