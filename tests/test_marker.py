import os
import tempfile

import pytest
import gemmi
from xtalkit.marker import mark


@pytest.fixture
def simple_cif():
    return os.path.join(os.path.dirname(__file__), "fixtures", "simple.cif")


def test_mark_overlay_mode(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        # Check output file exists
        assert os.path.exists(out + ".cif")
        # Check CIF can be re-parsed
        doc = gemmi.cif.read_file(out + ".cif")
        assert len(doc) > 0
        # Should contain both original and dummy atoms
        content = open(out + ".cif").read()
        assert "Li" in content  # original
        assert "WYCK_4a" in content  # dummy


def test_mark_replace_mode(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a"],
            mode="replace",
            tolerance=0.5,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")


def test_mark_all_wyckoff_positions(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h"],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        # All 8 wyckoff positions should have dummy atoms
        for letter in ["WYCK_4a", "WYCK_4b", "WYCK_4c", "WYCK_4d",
                       "WYCK_16e", "WYCK_24f", "WYCK_24g", "WYCK_48h"]:
            assert letter in content, f"Missing {letter}"


def test_mark_custom_element_map(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            mode="overlay",
            tolerance=0.5,
            element_map={"4a": "He", "4c": "Ne"},
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        assert "He" in content
        assert "Ne" in content


def test_mark_all_output_formats(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif", "vesta", "xyz"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")
        assert os.path.exists(out + ".vesta")
        assert os.path.exists(out + ".xyz")


def test_mark_returns_output_paths(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        result = mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a"],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        assert out + ".cif" in result


def test_mark_invalid_sg(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        with pytest.raises(ValueError):
            mark(
                cif_path=simple_cif,
                sg_number=999,
                wyckoff_letters=["4a"],
                mode="overlay",
                tolerance=0.5,
                element_map=None,
                formats=["cif"],
                output_base=out,
            )


def test_mark_invalid_wyckoff_letter(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        with pytest.raises(ValueError):
            mark(
                cif_path=simple_cif,
                sg_number=216,
                wyckoff_letters=["99z"],
                mode="overlay",
                tolerance=0.5,
                element_map=None,
                formats=["cif"],
                output_base=out,
            )
