import json
import os
import subprocess
import sys

from xtalkit.enumeration.shry_backend import ShryResult
from xtalkit.enumeration.shry_count import count_shry_structures
from xtalkit.enumeration.shry_enum import enumerate_with_shry
from xtalkit.enumeration.shry_prepare import prepare_shry_input
from xtalkit.enumeration.shry_verify import verify_shry_outputs
from xtalkit.enumeration.supercell_backend import parse_supercell_count
from xtalkit.enumeration.postprocess import (
    rank_by_shortest_distance,
    write_tblite_inputs,
    write_cp2k_inputs,
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

    def __init__(self, tmp_path=None):
        self.tmp_path = tmp_path
        self.commands = []

    def run(self, args, cwd=None):
        self.commands.append(args)
        if "--count-only" in args:
            return ShryResult(["shry", *args], 0, "inequivalent structures: 7\n", "")
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
        return ShryResult(["shry", *args], 0, "generated 1\n", "")

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


def test_shry_verify_symprec_scan_and_degeneracy(tmp_path):
    cif = tmp_path / "ready.cif"
    cif.write_text(P1_PARTIAL)
    out = tmp_path / "enum"
    backend = FakeShryBackend(tmp_path)
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

    assert [row["count"] for row in result["symprec_scan"]] == [7, 7]
    assert result["degeneracy_sum"] == 1


def test_supercell_count_parser():
    assert parse_supercell_count("inequivalent structures: 42") == 42
    assert parse_supercell_count("7") == 7


def test_dedup_two_stage_detects_confirmed_duplicate(tmp_path):
    """Plan §10: fingerprint bucketing then StructureMatcher within buckets."""
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
    p1 = tmp_path / "a.cif"; p1.write_text(base.format(x=0.1, y=0.2, z=0.3))
    p2 = tmp_path / "b.cif"; p2.write_text(base.format(x=0.1, y=0.2, z=0.3))  # identical
    p3 = tmp_path / "c.cif"; p3.write_text(base.format(x=0.4, y=0.5, z=0.6))  # distinct

    result = verify_no_duplicates([str(p1), str(p2), str(p3)])
    assert result["has_duplicates"] is True
    assert len(result["duplicates"]) == 1
    group = result["duplicates"][0]["groups"][0]
    assert {os.path.basename(g) for g in group} == {"a.cif", "b.cif"}


def test_symmetry_audit_rejects_split_orbit_and_symmetrize_passes(tmp_path):
    """Plan §3.4 / case_007: a P1-expanded CIF splits one orbit across labels.

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
        raise AssertionError("expected §3.4 orbit-split error")

    # --symmetrize repairs it.
    result = prepare_shry_input(str(inp), str(tmp_path / "out.cif"),
                                parent_spacegroup=123, symmetrize=True,
                                scaling_matrix=[[1, 0, 0], [0, 1, 0], [0, 0, 1]])
    assert result["site_orbits"] is not None


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
    cp2k_jobs = write_cp2k_inputs(str(out))
    assert len(tblite_jobs) == 1
    assert len(cp2k_jobs) == 1
    script = write_slurm_array(str(out), job_kind="tblite")
    assert os.path.exists(script)


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
