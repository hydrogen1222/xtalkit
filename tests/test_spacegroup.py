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
    """F-43m has 8 Wyckoff positions: 4a,4b,4c,4d,16e,24f,24g,48h."""
    positions = wyckoff_positions(216)
    letters = {p.letter for p in positions}
    assert letters == {"4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h"}
    total_multiplicity = sum(p.multiplicity for p in positions)
    assert total_multiplicity > 0


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
