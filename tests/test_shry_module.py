import json
import os
import shutil
import subprocess
import sys

import pytest

from xtalkit.enumeration.shry_backend import ShryBackend, ShryResult
from xtalkit.enumeration.shry_count import count_shry_structures
from xtalkit.enumeration.degeneracy import compute_degeneracy, raw_combination_count
from xtalkit.enumeration.shry_enum import enumerate_with_shry
from xtalkit.enumeration.shry_prepare import prepare_shry_input
from xtalkit.enumeration.shry_verify import verify_shry_outputs
from xtalkit.enumeration.supercell_backend import parse_supercell_count
from xtalkit.enumeration.postprocess import (
    rank_by_shortest_distance,
    write_tblite_inputs,
    write_slurm_array,
)


P1_PARTIAL = """data_partial
_cell_length_a 4.1
_cell_length_b 4.2
_cell_length_c 4.3
_cell_angle_alpha 71
_cell_angle_beta 82
_cell_angle_gamma 83
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li 0.11 0.23 0.37 1/2
Li2 Li 0.41 0.67 0.19 1/2
Li3 Li 0.83 0.05 0.51 1/2
"""


P1_ORDERED_WITH_X = """data_ordered
_cell_length_a 4.1
_cell_length_b 4.2
_cell_length_c 4.3
_cell_angle_alpha 71
_cell_angle_beta 82
_cell_angle_gamma 83
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li 0 0 0 1
X1 X 0.5 0.5 0.5 1
"""


class FakeShryBackend:
    """Mimics real SHRY 1.1.x CLI behaviour for unit tests.

    Real SHRY has no --output-dir/--write-cif flags: it writes
    ``shry-<base>-<ts>/sliceN/*.cif`` (ordered configs) and
    ``shry-<base>-<ts>/<base>-<scaling>.cif`` (modified parent) into its CWD.
    --mod-only writes only the modified parent; --count-only writes nothing.
    """

    def __init__(self, tmp_path=None, count=7):
        self.tmp_path = tmp_path
        self.count = count
        self.commands = []

    def run(self, args, cwd=None):
        self.commands.append(args)
        if "--count-only" in args:
            return ShryResult(
                ["shry", *args], 0,
                f"inequivalent structures: {self.count}\n", "",
            )
        if "--version" in args:
            return ShryResult(["shry", *args], 0, "SHRY 1.1.8\n", "")
        if "--mod-only" in args:
            if cwd:
                d = os.path.join(cwd, "shry-fake")
                os.makedirs(d, exist_ok=True)
                with open(os.path.join(d, "fake-1-1-1.cif"), "w", encoding="utf-8") as f:
                    f.write(P1_PARTIAL)
            return ShryResult(["shry", *args], 0, "ok\n", "")
        # Full enum: write one ordered config into slice0/.
        if cwd:
            raw = os.path.join(cwd, "shry-fake", "slice0")
            os.makedirs(raw, exist_ok=True)
            with open(os.path.join(raw, "conf_000001.cif"), "w", encoding="utf-8") as f:
                f.write(P1_ORDERED_WITH_X)
        return ShryResult(
            ["shry", *args], 0,
            "Total number of combinations is 1\ngenerated 1\n", "",
        )

    def version(self):
        return "SHRY 1.1.8"


def test_shry_prepare_fills_vacancy_and_writes_manifest(tmp_path):
    inp = tmp_path / "partial.cif"
    out = tmp_path / "ready.cif"
    inp.write_text(P1_PARTIAL)

    result = prepare_shry_input(
        str(inp),
        str(out),
        parent_spacegroup=1,
        scaling_matrix=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        target_formula="Li3",
    )

    text = out.read_text()
    assert "Li1_X" in text
    assert " X " in text
    manifest = json.loads((tmp_path / "ready.cif.manifest.json").read_text())
    assert manifest["parent_spacegroup_detected"] == 1
    assert manifest["target_formula_after_removing_vacancy"] == "Li3"
    assert result["output_cif"] == str(out)


def test_shry_count_uses_backend_and_records_atol(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)

    result = count_shry_structures(
        str(cif),
        scaling_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]],
        atol=2e-5,
        out_json=str(tmp_path / "count.json"),
        backend=FakeShryBackend(),
    )

    assert result["count_only_result"] == 7
    assert result["atol"] == 2e-5
    assert (tmp_path / "count.json").exists()


def test_shry_count_rejects_zero_configurations(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    with pytest.raises(ValueError, match="at least one configuration"):
        count_shry_structures(
            str(cif), backend=FakeShryBackend(count=0),
            out_json=str(tmp_path / "count.json"),
        )
    assert not (tmp_path / "count.json").exists()


def test_shry_workflows_reject_invalid_matrix_and_tolerances(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    with pytest.raises(ValueError, match="determinant"):
        count_shry_structures(
            str(cif), scaling_matrix=[[0, 0, 0], [0, 1, 0], [0, 0, 1]],
            backend=FakeShryBackend(),
        )
    with pytest.raises(ValueError, match="symprec"):
        count_shry_structures(
            str(cif), symprec=float("nan"), backend=FakeShryBackend(),
        )
    with pytest.raises(ValueError, match="angle_tolerance"):
        count_shry_structures(
            str(cif), angle_tolerance=-2, backend=FakeShryBackend(),
        )
    with pytest.raises(ValueError, match="finite integer"):
        count_shry_structures(
            str(cif), scaling_matrix=[[float("inf"), 0, 0], [0, 1, 0], [0, 0, 1]],
            backend=FakeShryBackend(),
        )
    with pytest.raises(ValueError, match="exactly 3x3"):
        count_shry_structures(
            str(cif), scaling_matrix="bad", backend=FakeShryBackend(),
        )


def test_shry_enum_streams_raw_outputs_and_removes_vacancy(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"

    result = enumerate_with_shry(
        str(cif),
        str(out),
        expect_count=1,
        target_formula="Li1",
        write_cif_output=True,
        write_poscar=True,
        write_degeneracy=True,
        backend=FakeShryBackend(tmp_path),
    )

    assert result["generated_clean_count"] == 1
    clean = (out / "clean_cif" / "conf_000001.cif").read_text()
    assert " X " not in clean
    assert (out / "poscar" / "POSCAR_000001").exists()
    line = json.loads((out / "manifest.jsonl").read_text().splitlines()[0])
    assert line["degeneracy"] == 1


def test_shry_enum_rejects_degeneracy_for_nonunit_supercell(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"

    with pytest.raises(ValueError, match="unit-determinant"):
        enumerate_with_shry(
            str(cif), str(out),
            scaling_matrix=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
            expect_count=1,
            write_cif_output=True,
            write_degeneracy=True,
            backend=FakeShryBackend(tmp_path),
        )
    assert not out.exists()


def test_shry_enum_rejects_real_element_as_vacancy(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    with pytest.raises(ValueError, match="must not collide"):
        enumerate_with_shry(
            str(cif), str(tmp_path / "enum"), expect_count=1,
            remove_vacancy="Li", backend=FakeShryBackend(tmp_path),
        )


def test_shry_enum_rejects_degeneracy_when_declared_group_differs(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    (tmp_path / "ready.cif.orbit_grouping.json").write_text(json.dumps({
        "parent_spacegroup_detected": 2,
        "site_orbits": [],
        "hall_symbol": "-P 1",
    }))
    with pytest.raises(ValueError, match="declared 1, detected 2"):
        enumerate_with_shry(
            str(cif), str(tmp_path / "enum"), expect_count=1,
            write_cif_output=True, write_degeneracy=True,
            backend=FakeShryBackend(tmp_path),
        )
    failure = json.loads((tmp_path / "enum" / "manifest.json").read_text())
    assert failure["failure_stage"] == "input_audit_metadata"


def test_degeneracy_matches_across_periodic_boundary(tmp_path):
    parent = tmp_path / "parent.cif"
    structure = tmp_path / "structure.cif"
    p_minus_one = P1_ORDERED_WITH_X.replace(
        "'P 1'", "'P -1'",
    ).replace(
        "_symmetry_Int_Tables_number 1",
        "_symmetry_Int_Tables_number 2",
    ).replace(
        "1 'x,y,z'",
        "1 'x,y,z'\n2 '-x,-y,-z'",
    ).replace(
        "Li1 Li 0 0 0 1",
        "Li1 Li 0.99999 0 0 1",
    )
    parent.write_text(p_minus_one)
    structure.write_text(p_minus_one)

    assert compute_degeneracy(str(structure), str(parent), tol=1e-4) == 1


def test_degeneracy_rejects_incomplete_declared_group(tmp_path):
    parent = tmp_path / "parent.cif"
    structure = tmp_path / "structure.cif"
    incomplete = P1_ORDERED_WITH_X.replace("'P 1'", "'P -1'").replace(
        "_symmetry_Int_Tables_number 1", "_symmetry_Int_Tables_number 2",
    )
    parent.write_text(incomplete)
    structure.write_text(P1_ORDERED_WITH_X)

    with pytest.raises(ValueError, match="declared space group has order 2"):
        compute_degeneracy(str(structure), str(parent))


def test_raw_combination_count_is_exact_multinomial(tmp_path):
    parent = tmp_path / "parent.cif"
    parent.write_text(P1_PARTIAL.replace(
        "Li1 Li 0.11 0.23 0.37 1/2",
        "Li1 Li 0.11 0.23 0.37 1/2\nX1 X 0.11 0.23 0.37 1/2",
    ))
    site_orbits = [{
        "orbit_id": 1,
        "representative": ["0.110000", "0.230000", "0.370000"],
        "multiplicity": 2,
    }]
    assert raw_combination_count(str(parent), site_orbits) == 2


def test_shry_enum_tracks_outputs_that_are_actually_written(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum_no_files"

    result = enumerate_with_shry(
        str(cif),
        str(out),
        strict=False,
        expect_count=None,
        target_formula="Li1",
        write_cif_output=False,
        write_poscar=False,
        backend=FakeShryBackend(tmp_path),
    )

    assert result["generated_raw_count"] == 1
    assert result["generated_clean_count"] == 0
    assert result["generated_poscar_count"] == 0


def test_shry_enum_rejects_nonempty_output_directory(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    out.mkdir()
    (out / "old-result.cif").write_text("stale")

    with pytest.raises(FileExistsError, match="not empty"):
        enumerate_with_shry(
            str(cif), str(out), expect_count=1,
            backend=FakeShryBackend(tmp_path),
        )


def test_shry_enum_count_mismatch_writes_failure_manifest(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum_mismatch"

    with pytest.raises(ValueError, match="generated 1 structures, expected 2"):
        enumerate_with_shry(
            str(cif), str(out), expect_count=2,
            target_formula="Li1", backend=FakeShryBackend(tmp_path),
        )

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["failure_reason"] == "generated_count_mismatch"
    assert manifest["generated_raw_count"] == 1


def test_shry_enum_postprocess_failure_is_auditable(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum_bad_formula"

    with pytest.raises(ValueError, match="formula mismatch"):
        enumerate_with_shry(
            str(cif), str(out), expect_count=1, target_formula="Na1",
            write_cif_output=True, backend=FakeShryBackend(tmp_path),
        )

    manifest = json.loads((out / "manifest.json").read_text())
    assert manifest["status"] == "failed"
    assert manifest["failure_stage"] == "postprocess"
    assert manifest["generated_raw_count"] == 1
    assert manifest["generated_clean_count"] == 0


def test_shry_poscar_preserves_nonorthogonal_cell(tmp_path):
    """POSCAR lattice vectors must preserve all six CIF cell parameters."""
    pytest.importorskip("pymatgen")
    from pymatgen.core import Structure
    from xtalkit.enumeration.cif_io import AtomRow, CifData
    from xtalkit.enumeration.shry_enum import _write_poscar

    data = CifData(
        "triclinic",
        {"a": "4.1", "b": "4.2", "c": "4.3",
         "alpha": "71", "beta": "82", "gamma": "83"},
        "P 1", 1, ["x,y,z"],
        [AtomRow("Li1", "Li", "0.1", "0.2", "0.3", "1")],
    )
    path = tmp_path / "POSCAR"
    _write_poscar(data, str(path))
    structure = Structure.from_file(path)
    assert structure.lattice.abc == pytest.approx((4.1, 4.2, 4.3), abs=1e-8)
    assert structure.lattice.angles == pytest.approx((71, 82, 83), abs=1e-8)


def test_shry_verify_checks_count_formula_and_dedup(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    enumerate_with_shry(
        str(cif),
        str(out),
        expect_count=1,
        target_formula="Li1",
        write_cif_output=True,
        backend=FakeShryBackend(tmp_path),
    )

    result = verify_shry_outputs(
        str(out),
        check_count=True,
        check_formula=True,
        check_dedup=True,
        target_formula="Li1",
    )

    assert result["clean_count"] == 1
    assert result["count_ok"] is True
    formula_rows = (out / "checks" / "formula_check.jsonl").read_text().splitlines()
    assert len(formula_rows) == 1
    assert isinstance(json.loads(formula_rows[0]), dict)


def test_shry_verify_symprec_scan_and_degeneracy(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    backend = FakeShryBackend(tmp_path, count=1)
    enumerate_with_shry(
        str(cif),
        str(out),
        expect_count=1,
        target_formula="Li1",
        write_cif_output=True,
        backend=backend,
    )

    result = verify_shry_outputs(
        str(out),
        check_count=True,
        target_formula="Li1",
        symprec_list=[1e-3, 1e-2],
        backend=backend,
        check_degeneracy=True,
    )

    assert [row["count"] for row in result["symprec_scan"]] == [1, 1]
    assert result["degeneracy_sum"] == 1


def test_shry_verify_records_all_failures_before_raising(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    enumerate_with_shry(
        str(cif), str(out), expect_count=1, target_formula="Li1",
        write_cif_output=True, backend=FakeShryBackend(tmp_path),
    )
    manifest_path = out / "manifest.json"
    manifest = json.loads(manifest_path.read_text())
    manifest["count_only_result"] = 2
    manifest_path.write_text(json.dumps(manifest))

    with pytest.raises(ValueError, match="SHRY verification failed"):
        verify_shry_outputs(
            str(out), check_count=True, check_formula=True,
            check_dedup=False, target_formula="Na1",
        )

    report = json.loads((out / "checks" / "verify.json").read_text())
    assert report["ok"] is False
    assert {failure["check"] for failure in report["failures"]} == {
        "count", "formula",
    }


def test_shry_verify_never_accepts_an_empty_clean_set(tmp_path):
    clean = tmp_path / "clean_cif"
    clean.mkdir(parents=True)
    with pytest.raises(ValueError, match="outputs"):
        verify_shry_outputs(
            str(tmp_path), check_count=False, check_formula=False,
            check_dedup=False,
        )
    report = json.loads((tmp_path / "checks" / "verify.json").read_text())
    assert report["ok"] is False


def test_supercell_count_parser():
    assert parse_supercell_count("inequivalent structures: 42") == 42
    assert parse_supercell_count("7") == 7


def test_dedup_two_stage_detects_confirmed_duplicate(tmp_path):
    """Candidate bucketing is followed by StructureMatcher within buckets."""
    from xtalkit.enumeration.fingerprint import verify_no_duplicates

    base = """data_c
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li {x} {y} {z} 1
"""
    p1 = tmp_path / "a.cif"
    p1.write_text(base.format(x=0.1, y=0.2, z=0.3))
    p2 = tmp_path / "b.cif"
    p2.write_text(base.format(x=0.1, y=0.2, z=0.3))  # identical
    # One-site P1 structures at different absolute origins are translation-
    # equivalent and must not be called distinct by a strict verifier.
    p3 = tmp_path / "c.cif"
    p3.write_text(base.format(x=0.4, y=0.5, z=0.6))

    result = verify_no_duplicates([str(p1), str(p2), str(p3)])
    assert result["has_duplicates"] is True
    assert len(result["duplicates"]) == 1
    group = result["duplicates"][0]["groups"][0]
    assert {os.path.basename(g) for g in group} == {"a.cif", "b.cif", "c.cif"}


def test_dedup_does_not_gate_on_origin_dependent_fingerprint(tmp_path):
    """Translation-equivalent structures must reach StructureMatcher."""
    from xtalkit.enumeration.fingerprint import (
        structure_fingerprint,
        verify_no_duplicates,
    )
    template = """data_t
_cell_length_a 5
_cell_length_b 5
_cell_length_c 5
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li {x1} 0.2 0.3 1
Li2 Li {x2} 0.2 0.3 1
"""
    a = tmp_path / "origin_a.cif"
    b = tmp_path / "origin_b.cif"
    a.write_text(template.format(x1=0.1, x2=0.4))
    b.write_text(template.format(x1=0.3, x2=0.6))
    assert structure_fingerprint(str(a)) != structure_fingerprint(str(b))
    result = verify_no_duplicates([str(a), str(b)])
    assert result["has_duplicates"] is True


def test_dedup_reports_direct_pairs_for_nontransitive_matcher():
    from xtalkit.enumeration.fingerprint import _match_graph

    class NonTransitiveMatcher:
        def fit(self, left, right):
            return abs(left - right) <= 1

    groups, pairs = _match_graph([0, 1, 2], NonTransitiveMatcher())
    assert groups == [[0, 1, 2]]
    assert pairs == [(0, 1), (1, 2)]


def test_symmetry_audit_rejects_split_orbit_and_symmetrize_passes(tmp_path):
    """A P1-expanded CIF that splits one orbit across labels is rejected.

    Declaring the *detected* group exposes the orbit-split; --symmetrize must
    repair it (the documented escape hatch).
    """
    p1_expanded = """data_p1exp
_cell_length_a 4
_cell_length_b 4
_cell_length_c 4
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
Li1 Li 0.1 0.5 0.5 1/2
Li2 Li 0.9 0.5 0.5 1/2
"""
    inp = tmp_path / "p1exp.cif"
    inp.write_text(p1_expanded)

    # Declared group matches detected (123) but labels split the orbit -> error.
    try:
        prepare_shry_input(str(inp), str(tmp_path / "out.cif"),
                           parent_spacegroup=123, scaling_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    except ValueError as exc:
        assert "does not match spglib symmetry orbits" in str(exc)
    else:
        raise AssertionError("expected orbit-split error")

    # --symmetrize repairs it.
    result = prepare_shry_input(str(inp), str(tmp_path / "out.cif"),
                                parent_spacegroup=123, symmetrize=True,
                                scaling_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert result["site_orbits"] is not None


def test_symmetry_audit_signature_includes_exact_occupancies():
    from xtalkit.enumeration.cif_io import AtomRow, CifData
    from xtalkit.enumeration.symmetry_audit import _expand_full_cell

    atoms = [
        AtomRow("A_Li", "Li", "0.1", "0.2", "0.3", "1/2"),
        AtomRow("A_X", "X", "0.1", "0.2", "0.3", "1/2"),
        AtomRow("B_Li", "Li", "0.4", "0.2", "0.3", "1/4"),
        AtomRow("B_X", "X", "0.4", "0.2", "0.3", "3/4"),
    ]
    data = CifData(
        "occ", {"a": "4", "b": "4", "c": "4", "alpha": "90",
                "beta": "90", "gamma": "90"},
        "P 1", 1, ["x,y,z"], atoms,
    )
    groups = {}
    for atom in atoms:
        groups.setdefault(atom.frac_key, []).append(atom)

    full, _ = _expand_full_cell(data, groups)
    assert len(full) == 2
    assert full[0].species_number != full[1].species_number


def test_symmetry_audit_preserves_original_label_underscores():
    from xtalkit.enumeration.cif_io import AtomRow
    from xtalkit.enumeration.symmetry_audit import (
        _base_label,
        _coord_group_base_labels,
    )

    assert _base_label("Li_A1", "Li") == "Li_A1"
    assert _base_label("M1_Ge", "Ge") == "M1"
    assert _base_label("Li1_X", "X") == "Li1"
    assert _coord_group_base_labels([
        AtomRow("site_Li", "Li", "0", "0", "0", "1"),
    ]) == {"site_Li"}
    assert _coord_group_base_labels([
        AtomRow("M1_Ge", "Ge", "0", "0", "0", "1/2"),
        AtomRow("M1_P", "P", "0", "0", "0", "1/2"),
    ]) == {"M1"}


def test_prepare_drops_zero_occupancy_and_rejects_negative_occupancy():
    from xtalkit.enumeration.cif_io import AtomRow
    from xtalkit.enumeration.shry_prepare import _fill_and_validate_vacancies

    zero = AtomRow("Li0", "Li", "0", "0", "0", "0")
    rows = _fill_and_validate_vacancies([zero], "X")
    assert [(row.type_symbol, row.occupancy) for row in rows] == [("X", "1")]

    negative = AtomRow("Li1", "Li", "0", "0", "0", "-1/2")
    with pytest.raises(ValueError, match="between 0 and 1"):
        _fill_and_validate_vacancies([negative], "X")


def test_occupancy_overrides_reject_duplicate_declarations():
    from xtalkit.enumeration.occupancy import parse_set_occupancy

    with pytest.raises(ValueError, match="duplicate"):
        parse_set_occupancy("Li1:1/2 Li1:3/4")


@pytest.mark.parametrize("formula", ["", "Li0", "Xx2", "Li-1"])
def test_formula_parser_rejects_nonphysical_formulas(formula):
    from xtalkit.enumeration.formula import parse_formula

    with pytest.raises(ValueError):
        parse_formula(formula)


def test_shry_output_iterator_uses_natural_numeric_order(tmp_path):
    from xtalkit.enumeration.shry_enum import iter_shry_outputs

    for slice_name, files in {
        "slice10": ["conf_10.cif"],
        "slice2": ["conf_11.cif", "conf_2.cif"],
    }.items():
        directory = tmp_path / "run" / slice_name
        directory.mkdir(parents=True)
        for name in files:
            (directory / name).write_text("data_x")
    names = [
        f"{os.path.basename(os.path.dirname(path))}/{os.path.basename(path)}"
        for path in iter_shry_outputs(str(tmp_path))
    ]
    assert names == ["slice2/conf_2.cif", "slice2/conf_11.cif", "slice10/conf_10.cif"]


def test_shortest_distance_rejects_residual_vacancy(tmp_path):
    from xtalkit.enumeration.postprocess import shortest_distance

    cif = tmp_path / "with_vacancy.cif"
    cif.write_text(P1_ORDERED_WITH_X)
    with pytest.raises(ValueError, match="unknown chemical element"):
        shortest_distance(str(cif))


def test_prepare_checks_target_formula_before_writing(tmp_path):
    inp = tmp_path / "partial.cif"
    out = tmp_path / "ready.cif"
    inp.write_text(P1_PARTIAL)
    with pytest.raises(ValueError, match="target formula mismatch"):
        prepare_shry_input(
            str(inp), str(out), target_formula="Li4",
            scaling_matrix=[[2, 0, 0], [0, 1, 0], [0, 0, 1]],
        )
    assert not out.exists()


def test_cif_reader_rejects_spacegroup_name_number_disagreement(tmp_path):
    from xtalkit.enumeration.cif_io import read_cif

    cif = tmp_path / "bad_sg.cif"
    cif.write_text(P1_PARTIAL.replace("'P 1'", "'P -1'"))
    with pytest.raises(ValueError, match="name and International Tables number disagree"):
        read_cif(str(cif))


def test_mod_only_composition_uses_expanded_wyckoff_multiplicities(tmp_path):
    from pymatgen.core import Structure
    from pymatgen.io.cif import CifWriter
    from xtalkit.builder import AtomSite, build_structure
    from xtalkit.enumeration.shry_enum import _composition_diff
    from xtalkit.exporter import write_structure_cif

    cell = {"a": 5.46, "b": 5.46, "c": 5.46,
            "alpha": 90, "beta": 90, "gamma": 90}
    source = tmp_path / "fluorite_asymmetric.cif"
    write_structure_cif(build_structure(225, cell, [
        AtomSite("Ca", "4a", [], 1),
        AtomSite("F", "8c", [], 1),
    ]), str(source))

    expanded_structure = Structure.from_file(source)
    for index, site in enumerate(expanded_structure, 1):
        site.label = f"{site.specie.symbol}{index}"
    expanded = tmp_path / "fluorite_p1.cif"
    CifWriter(expanded_structure, symprec=None).write_file(expanded)
    assert _composition_diff(str(source), str(expanded)) == ""

    expanded_structure.remove_sites([len(expanded_structure) - 1])
    altered = tmp_path / "altered.cif"
    CifWriter(expanded_structure, symprec=None).write_file(altered)
    assert _composition_diff(str(source), str(altered))


def test_shry_cif_reader_defaults_missing_occupancy_to_one(tmp_path):
    from xtalkit.enumeration.cif_io import read_cif

    cif = tmp_path / "no_occ.cif"
    cif.write_text(P1_ORDERED_WITH_X.replace("_atom_site_occupancy\n", "").replace(
        "Li1 Li 0 0 0 1", "Li1 Li 0 0 0"
    ).replace("X1 X 0.5 0.5 0.5 1", "X1 X 0.5 0.5 0.5"))
    assert [atom.occupancy for atom in read_cif(str(cif)).atoms] == ["1", "1"]


def test_shry_cif_reader_accepts_esd_before_exponent(tmp_path):
    from xtalkit.enumeration.cif_io import read_cif

    cif = tmp_path / "esd.cif"
    cif.write_text(P1_ORDERED_WITH_X.replace(
        "Li1 Li 0 0 0 1", "Li1 Li 1.2(3)e-2 0 0 1"
    ))
    assert float(read_cif(str(cif)).atoms[0].x) == pytest.approx(0.012)


def test_symmetry_audit_rejects_conflicting_symmetry_overlap():
    from xtalkit.enumeration.cif_io import AtomRow, CifData
    from xtalkit.enumeration.symmetry_audit import _expand_full_cell

    atoms = [
        AtomRow("A_Li", "Li", "0.1", "0.2", "0.3", "1/2"),
        AtomRow("A_X", "X", "0.1", "0.2", "0.3", "1/2"),
        AtomRow("B_Li", "Li", "0.9", "0.2", "0.3", "1/4"),
        AtomRow("B_X", "X", "0.9", "0.2", "0.3", "3/4"),
    ]
    data = CifData(
        "overlap", {"a": "4", "b": "4", "c": "4", "alpha": "90",
                    "beta": "90", "gamma": "90"},
        "P m", None, ["x,y,z", "-x,y,z"], atoms,
    )
    groups = {}
    for atom in atoms:
        groups.setdefault(atom.frac_key, []).append(atom)

    with pytest.raises(ValueError, match="different occupancy compositions"):
        _expand_full_cell(data, groups)


def test_shry_postprocess_writes_rankings_and_jobs(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    enumerate_with_shry(
        str(cif),
        str(out),
        expect_count=1,
        target_formula="Li1",
        write_cif_output=True,
        backend=FakeShryBackend(tmp_path),
    )

    rows = rank_by_shortest_distance(str(out))
    assert len(rows) == 1
    assert rows[0]["shortest_distance"] is None
    tblite_jobs = write_tblite_inputs(str(out))
    assert (out / "tblite" / "conf_000001" / "POSCAR").exists()
    tblite_script = (out / "tblite" / "conf_000001" / "run.sh").read_text()
    assert "tblite run" in tblite_script
    assert "--charge 0 --spin 0 POSCAR" in tblite_script
    assert len(tblite_jobs) == 1
    script = write_slurm_array(str(out))
    assert os.path.exists(script)


def test_job_generation_rejects_invalid_spin(tmp_path):
    clean = tmp_path / "clean_cif"
    clean.mkdir()
    (clean / "conf.cif").write_text(P1_ORDERED_WITH_X.replace(
        "X1 X 0.5 0.5 0.5 1\n", ""
    ))
    with pytest.raises(ValueError, match="non-negative"):
        write_tblite_inputs(str(tmp_path), uhf=-1)


@pytest.mark.skipif(shutil.which("tblite") is None, reason="tblite CLI not installed")
def test_real_tblite_generated_periodic_job_smoke(tmp_path):
    """Run one generated periodic H2 job with the installed tblite binary."""
    clean = tmp_path / "clean_cif"
    clean.mkdir()
    (clean / "H2.cif").write_text("""data_h2
_cell_length_a 10
_cell_length_b 10
_cell_length_c 10
_cell_angle_alpha 90
_cell_angle_beta 90
_cell_angle_gamma 90
_symmetry_space_group_name_H-M 'P 1'
_symmetry_Int_Tables_number 1
loop_
_symmetry_equiv_pos_site_id
_symmetry_equiv_pos_as_xyz
1 'x,y,z'
loop_
_atom_site_label
_atom_site_type_symbol
_atom_site_fract_x
_atom_site_fract_y
_atom_site_fract_z
_atom_site_occupancy
H1 H 0.463 0.5 0.5 1
H2 H 0.537 0.5 0.5 1
""")
    [job] = write_tblite_inputs(str(tmp_path), charge=0, uhf=0)
    env = os.environ.copy()
    env.update({
        "OMP_NUM_THREADS": "1",
        "OPENBLAS_NUM_THREADS": "1",
        "MKL_NUM_THREADS": "1",
    })
    result = subprocess.run(
        ["./run.sh"], cwd=job, env=env, capture_output=True, text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stderr
    output = os.path.join(job, "tblite.out")
    with open(output, encoding="utf-8") as handle:
        assert "total energy" in handle.read().lower()


def test_cli_shry_help():
    result = subprocess.run(
        [sys.executable, "-m", "xtalkit.cli", "shry", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "prepare" in result.stdout
    assert "count" in result.stdout


def test_cli_shry_postprocess_help():
    result = subprocess.run(
        [sys.executable, "-m", "xtalkit.cli", "shry", "postprocess", "--help"],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert "--write-tblite" in result.stdout


@pytest.mark.skipif(shutil.which("shry") is None, reason="real SHRY CLI not installed")
def test_real_shry_count_and_enumeration_smoke(tmp_path, monkeypatch):
    """Exercise the installed isolated SHRY backend on a two-site supercell."""
    for name in ("OMP_NUM_THREADS", "OPENBLAS_NUM_THREADS", "MKL_NUM_THREADS"):
        monkeypatch.setenv(name, "1")
    source = os.path.join(
        os.path.dirname(__file__), "fixtures", "disordered_binary.cif"
    )
    matrix = [[2, 0, 0], [0, 1, 0], [0, 0, 1]]
    backend = ShryBackend()
    assert backend.version() == "1.1.8"

    count = count_shry_structures(
        source,
        scaling_matrix=matrix,
        out_json=str(tmp_path / "count.json"),
        backend=backend,
    )
    assert count["count_only_result"] == 1
    assert count["total_combinations"] == 2

    result = enumerate_with_shry(
        source,
        str(tmp_path / "enum"),
        scaling_matrix=matrix,
        expect_count=1,
        target_formula="Au1Cu1",
        write_cif_output=True,
        write_poscar=True,
        backend=backend,
    )
    assert result["generated_raw_count"] == 1
    assert result["generated_clean_count"] == 1
    assert result["generated_poscar_count"] == 1
    assert result["backend_version"] == "1.1.8"
