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
        mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            mode="overlay",
            tolerance=0.01,
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
        mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a"],
            mode="replace",
            tolerance=0.01,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")
        block = gemmi.cif.read_file(out + ".cif").sole_block()
        labels = list(block.find_values("_atom_site_label"))
        types = list(block.find_values("_atom_site_type_symbol"))
        occupancies = list(block.find_values("_atom_site_occupancy"))
        assert labels == ["WYCK_4a", "P2"]
        assert types == ["Xe", "P"]
        assert [float(value) for value in occupancies] == [1.0, 1.0]
        assert block.find_value("_xtalkit_output_purpose") == "visualization_only"


def test_mark_all_wyckoff_positions(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=[
                "4a", "4b", "4c", "4d", "16e", "24f", "24g", "48h", "96i",
            ],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif"],
            output_base=out,
        )
        content = open(out + ".cif").read()
        # All 9 Wyckoff positions should have dummy atoms.
        for letter in ["WYCK_4a", "WYCK_4b", "WYCK_4c", "WYCK_4d",
                       "WYCK_16e", "WYCK_24f", "WYCK_24g", "WYCK_48h",
                       "WYCK_96i"]:
            assert letter in content, f"Missing {letter}"


def test_mark_custom_element_map(simple_cif):
    with tempfile.TemporaryDirectory() as tmp:
        out = os.path.join(tmp, "simple_WYCK")
        mark(
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
        mark(
            cif_path=simple_cif,
            sg_number=216,
            wyckoff_letters=["4a", "4c"],
            mode="overlay",
            tolerance=0.5,
            element_map=None,
            formats=["cif", "xyz"],
            output_base=out,
        )
        assert os.path.exists(out + ".cif")
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


def test_mark_rejects_requested_spacegroup_mismatch(tmp_path):
    source = os.path.join(
        os.path.dirname(__file__), "fixtures", "disordered_binary.cif"
    )
    with pytest.raises(ValueError, match="disagrees with CIF-declared"):
        mark(
            cif_path=source,
            sg_number=2,
            wyckoff_letters=["1a"],
            mode="overlay",
            tolerance=0.01,
            element_map=None,
            formats=["cif"],
            output_base=str(tmp_path / "bad"),
        )


def test_mark_rejects_spacegroup_inconsistent_cell_metric(tmp_path):
    source = tmp_path / "bad_metric.cif"
    source.write_text("""data_bad
_cell_length_a 5
_cell_length_b 6
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F -4 3 m'
_symmetry_Int_Tables_number 216
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
Li1 Li 0 0 0
""")
    with pytest.raises(ValueError, match="unit-cell metric"):
        mark(
            str(source), 216, ["4a"], "overlay", 0.01, None, ["cif"],
            str(tmp_path / "bad"),
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


def test_mark_preserves_partial_occupancy_and_exact_default_position(tmp_path):
    source = os.path.join(
        os.path.dirname(__file__), "fixtures", "disordered_binary.cif"
    )
    out = str(tmp_path / "marked")
    mark(
        cif_path=source,
        sg_number=1,
        wyckoff_letters=["1a"],
        mode="overlay",
        tolerance=0.01,
        element_map={"1a": "He"},
        formats=["cif"],
        output_base=out,
    )
    block = gemmi.cif.read_file(out + ".cif").sole_block()
    labels = list(block.find_values("_atom_site_label"))
    occupancies = list(block.find_values("_atom_site_occupancy"))
    rows = dict(zip(labels, occupancies))
    assert float(rows["Au1"]) == pytest.approx(0.5)
    assert float(rows["Cu1"]) == pytest.approx(0.5)
    marker_index = labels.index("WYCK_1a")
    for tag in ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z"):
        assert float(list(block.find_values(tag))[marker_index]) == pytest.approx(
            [2 ** 0.5 / 10, 3 ** 0.5 / 10, 5 ** 0.5 / 10][
                ("_atom_site_fract_x", "_atom_site_fract_y", "_atom_site_fract_z").index(tag)
            ],
            abs=1e-6,
        )


def test_mark_replace_preserves_refined_parameter_and_occupancy(tmp_path):
    source = tmp_path / "parameterized.cif"
    source.write_text("""data_parameterized
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'F -4 3 m'
_symmetry_Int_Tables_number 216
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li 0.31 0.31 0.31 0.625
""")
    out = str(tmp_path / "replaced")
    mark(
        str(source), 216, ["16e"], "replace", 0.01,
        {"16e": "He"}, ["cif"], out,
    )
    block = gemmi.cif.read_file(out + ".cif").sole_block()
    assert list(block.find_values("_atom_site_label")) == ["WYCK_16e"]
    assert list(block.find_values("_atom_site_type_symbol")) == ["He"]
    assert float(block.find_value("_atom_site_fract_x")) == pytest.approx(0.31)
    assert float(block.find_value("_atom_site_fract_y")) == pytest.approx(0.31)
    assert float(block.find_value("_atom_site_fract_z")) == pytest.approx(0.31)
    assert float(block.find_value("_atom_site_occupancy")) == pytest.approx(0.625)


def test_mark_nonzero_offset_is_explicitly_visual_only(simple_cif, tmp_path):
    with pytest.warns(RuntimeWarning, match="visual-only"):
        mark(
            simple_cif, 216, ["4a"], "overlay", 0.01, None, ["cif"],
            str(tmp_path / "offset"), offset=0.02,
        )


def test_mark_rejects_negative_tolerance(simple_cif, tmp_path):
    with pytest.raises(ValueError, match="tolerance"):
        mark(
            simple_cif, 216, ["4a"], "overlay", -0.1, None, ["cif"],
            str(tmp_path / "bad"),
        )


def test_mark_rejects_empty_or_duplicate_requests(simple_cif, tmp_path):
    with pytest.raises(ValueError, match="at least one Wyckoff"):
        mark(
            simple_cif, 216, [], "overlay", 0.01, None, ["cif"],
            str(tmp_path / "empty"),
        )
    with pytest.raises(ValueError, match="duplicates"):
        mark(
            simple_cif, 216, ["4a", "4a"], "overlay", 0.01, None, ["cif"],
            str(tmp_path / "duplicate"),
        )
