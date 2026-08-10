"""Tests for xtalkit.skeleton — pure Wyckoff skeleton generation."""

import os
import tempfile

import pytest
import gemmi
from xtalkit.skeleton import generate


def test_generate_creates_cif():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")


def test_generate_valid_cif():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        doc = gemmi.cif.read_file(out + ".cif")
        assert len(doc) > 0


def test_generate_contains_dummy_atoms_only():
    """Skeleton should have no real atoms, only WYCK_* dummy atoms."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a"],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        assert "WYCK_4a" in content


def test_generate_writes_each_selected_position_once():
    """Regression: dummies must not be present in both structure and overlay."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        block = gemmi.cif.read_file(out + ".cif").sole_block()
        assert list(block.find_values("_atom_site_label")) == [
            "WYCK_4a", "WYCK_4c",
        ]


def test_generate_all_wyckoffs():
    """Should generate all 9 Wyckoff positions for SG 216."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=[
                "4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h", "96i",
            ],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        for letter in [
            "4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h", "96i",
        ]:
            assert f"WYCK_{letter}" in content


def test_generate_custom_cell():
    """Custom cell parameters should appear in output."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a"],
            cell_params={"a": 10.5, "b": 10.5, "c": 10.5,
                         "alpha": 90, "beta": 90, "gamma": 90},
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        assert "10.5" in content or "10.500" in content


def test_generate_multiple_formats():
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a"],
            cell_params=None,
            element_map=None,
            formats=["cif", "xyz"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")
        assert os.path.exists(out + ".xyz")
        lines = open(out + ".xyz").read().splitlines()
        assert int(lines[0]) == 4  # SG 216, Wyckoff 4a full orbit
        assert "symmetry-expanded" in lines[1]


def test_generate_custom_element_map():
    """Custom element_map should assign specified elements to Wyckoff letters."""
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "SG216_skeleton")
        generate(
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            cell_params=None,
            element_map={"4a": "He", "4c": "Ne"},
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        assert "He" in content
        assert "Ne" in content


def test_generate_invalid_sg():
    with pytest.raises(ValueError):
        generate(
            sg_number=0,
            wyckoff_letters=["4a"],
            cell_params=None,
            element_map=None,
            formats=["cif"],
            output_base="out",
        )


def test_generate_rejects_invalid_element_cell_and_format(tmp_path):
    base = str(tmp_path / "bad")
    with pytest.raises(ValueError, match="Invalid dummy element"):
        generate(216, ["4a"], None, {"4a": "NoSuch"}, ["cif"], base)
    with pytest.raises(ValueError, match="cell lengths"):
        generate(
            216, ["4a"],
            {"a": 0, "b": 5, "c": 5, "alpha": 90, "beta": 90, "gamma": 90},
            None, ["cif"], base,
        )
    with pytest.raises(ValueError, match="unsupported output"):
        generate(216, ["4a"], None, None, ["pdb"], base)
    with pytest.raises(ValueError, match="unit-cell metric"):
        generate(
            216, ["4a"],
            {"a": 5, "b": 6, "c": 5,
             "alpha": 90, "beta": 90, "gamma": 90},
            None, ["cif"], base,
        )
