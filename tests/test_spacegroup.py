import pytest
from xtalkit.spacegroup import (
    wyckoff_positions,
    crystal_system,
    default_cell_params,
    sg_name,
    WyckoffInfo,
)


def test_sg_name_f43m():
    """Space group 216 should be F-43m."""
    name = sg_name(216)
    assert "F" in name and "43" in name


def test_wyckoff_positions_f43m():
    """F-43m has 9 Wyckoff positions: 4a,4b,4c,4d,16e,24f,24g,48h,96i.

    The general position is 96i (mult 96 = 24 point-group ops x 4 F-centering);
    the legacy hand-coded DB omitted it. 24g is (x,1/4,1/4), not (x,0,1/2).
    """
    positions = wyckoff_positions(216)
    letters = {p.letter for p in positions}
    assert letters == {"4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h", "96i"}
    total_multiplicity = sum(p.multiplicity for p in positions)
    assert total_multiplicity > 0


def test_all_230_space_groups_have_wyckoff_data():
    """Every space group 1-230 now has Wyckoff data (bundled dataset)."""
    for n in range(1, 231):
        positions = wyckoff_positions(n)
        assert len(positions) > 0, f"SG {n} has no Wyckoff positions"
        for p in positions:
            assert p.multiplicity >= 1
            assert p.coordinates  # non-empty


def test_wyckoff_position_data():
    """Each WyckoffInfo has letter, multiplicity, site_symmetry, coordinates."""
    positions = wyckoff_positions(216)
    for p in positions:
        assert isinstance(p.letter, str)
        assert isinstance(p.multiplicity, int)
        assert isinstance(p.site_symmetry, str)
        assert isinstance(p.coordinates, str)


def test_crystal_system_cubic():
    """SG 216 should be cubic."""
    system = crystal_system(216)
    assert "cubic" in system.lower()


def test_crystal_system_monoclinic():
    """SG 14 should be monoclinic."""
    system = crystal_system(14)
    assert "monoclinic" in system.lower()


def test_default_cell_params_cubic():
    """Cubic systems get a=b=c."""
    params = default_cell_params(216)
    assert "a" in params and "b" in params and "c" in params
    assert params["alpha"] == 90 and params["beta"] == 90 and params["gamma"] == 90


def test_default_cell_params_returns_floats():
    """Cell params dict values should be floats."""
    params = default_cell_params(216)
    for k in ("a", "b", "c", "alpha", "beta", "gamma"):
        assert isinstance(params[k], float), f"{k} is not float"


def test_invalid_sg_number():
    """SG number out of range should raise ValueError."""
    with pytest.raises(ValueError):
        wyckoff_positions(0)
    with pytest.raises(ValueError):
        wyckoff_positions(231)


def test_sg_1_wyckoff():
    """SG 1 (P1) has exactly 1 Wyckoff position: 1a."""
    positions = wyckoff_positions(1)
    assert len(positions) == 1
    assert positions[0].letter == "1a"


def test_sg_2_wyckoff():
    """SG 2 (P-1) has Wyckoff positions."""
    positions = wyckoff_positions(2)
    letters = {p.letter for p in positions}
    assert len(positions) > 0
    assert "1a" in letters
