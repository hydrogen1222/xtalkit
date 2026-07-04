"""Tests for xtalkit.builder — building CIFs from refinement parameters."""

import os
import tempfile

import gemmi
import pytest

from xtalkit.builder import (
    AtomSite, FracAtom, free_params, eval_coord, find_wyckoff,
    build_structure, build_structure_frac, stoichiometry, stoichiometry_frac,
    format_formula, validate_cell, validate_atoms, validate_atoms_frac,
    detect_wyckoff,
)
from xtalkit.exporter import write_structure_cif, write_structure_xyz


# --- free_params / eval_coord ---------------------------------------------

@pytest.mark.parametrize("coord,expected", [
    ("x,-x,z", ["x", "z"]),
    ("x,x,x", ["x"]),
    ("1/3,2/3,z", ["z"]),
    ("x,2x,z", ["x", "z"]),
    ("x,1/4+x,7/8", ["x"]),
    ("0,0,0", []),
    ("1/2,1/2,1/2", []),
])
def test_free_params(coord, expected):
    assert free_params(coord) == expected


@pytest.mark.parametrize("coord,vals,expected", [
    ("x,-x,z", [0.3, 0.4], (0.3, 0.7, 0.4)),
    ("x,x,x", [0.25], (0.25, 0.25, 0.25)),
    ("1/3,2/3,z", [0.3], (1/3, 2/3, 0.3)),
    ("x,2x,z", [0.2, 0.5], (0.2, 0.4, 0.5)),
    ("x,1/4+x,7/8", [0.3], (0.3, 0.55, 0.875)),
    ("0,0,0", [], (0.0, 0.0, 0.0)),
    ("1/2,1/2,1/2", [], (0.5, 0.5, 0.5)),
])
def test_eval_coord(coord, vals, expected):
    got = eval_coord(coord, vals)
    assert got == pytest.approx(expected, abs=1e-4)


def test_eval_coord_wrong_arg_count():
    with pytest.raises(ValueError, match="expects 1 free"):
        eval_coord("x,x,x", [0.2, 0.3])


def test_eval_coord_wraps_to_unit_cell():
    """Negative results wrap into [0, 1)."""
    assert eval_coord("x,-x,z", [0.3, 0.4]) == pytest.approx((0.3, 0.7, 0.4), abs=1e-4)
    assert 0.0 <= eval_coord("x,-x,z", [0.1, 0.1])[0] < 1.0


# --- find_wyckoff ----------------------------------------------------------

def test_find_wyckoff_valid():
    wp = find_wyckoff(216, "16e")
    assert wp.letter == "16e"
    assert wp.coordinates == "x,x,x"


def test_find_wyckoff_invalid_label():
    with pytest.raises(ValueError, match="Invalid Wyckoff label"):
        find_wyckoff(216, "99z")


# --- build_structure -------------------------------------------------------

CUBIC = {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}


def test_build_nacl():
    """NaCl (Fm-3m, 225): Na 4a, Cl 4b."""
    atoms = [AtomSite("Na", "4a", [], 1.0), AtomSite("Cl", "4b", [], 1.0)]
    s = build_structure(225, CUBIC, atoms)
    rows = [(a.name, a.element.name, round(a.occ, 2)) for m in s for c in m for r in c for a in r]
    assert rows == [("Na1", "Na", 1.0), ("Cl1", "Cl", 1.0)]
    assert s.spacegroup_hm  # SG set


def test_build_with_free_param():
    """F-43m Li on 16e (x,x,x) with x=0.3."""
    s = build_structure(216, CUBIC, [AtomSite("Li", "16e", [0.3], 1.0)])
    atom = next(a for m in s for c in m for r in c for a in r)
    frac = s.cell.fractionalize(atom.pos)
    assert (frac.x, frac.y, frac.z) == pytest.approx((0.3, 0.3, 0.3), abs=1e-4)


def test_build_hexagonal():
    """P63/m (176): Mg 2c (1/3,2/3,1/4) + H 4f (1/3,2/3,z=0.3)."""
    cell = {"a": 3.2, "b": 3.2, "c": 5.1, "alpha": 90.0, "beta": 90.0, "gamma": 120.0}
    atoms = [AtomSite("Mg", "2c", [], 1.0), AtomSite("H", "4f", [0.3], 1.0)]
    s = build_structure(176, cell, atoms)
    fracs = [s.cell.fractionalize(a.pos) for m in s for c in m for r in c for a in r]
    # gemmi.Fractional has no __eq__; compare component-wise.
    assert (fracs[0].x, fracs[0].y, fracs[0].z) == pytest.approx((1/3, 2/3, 0.25), abs=1e-3)
    assert (fracs[1].x, fracs[1].y, fracs[1].z) == pytest.approx((1/3, 2/3, 0.3), abs=1e-3)


def test_build_partial_occupancy():
    """Disordered Li0.5/Cu0.5 on 4a (Fm-3m) — two atoms same site."""
    atoms = [AtomSite("Li", "4a", [], 0.5), AtomSite("Cu", "4a", [], 0.5)]
    s = build_structure(225, CUBIC, atoms)
    occs = sorted(round(a.occ, 2) for m in s for c in m for r in c for a in r)
    assert occs == [0.5, 0.5]


def test_build_bad_element():
    with pytest.raises(ValueError, match="Unknown element"):
        build_structure(225, CUBIC, [AtomSite("Xx", "4a", [], 1.0)])


def test_build_bad_free_count():
    """16e (dof 1) given 2 free values is an error."""
    with pytest.raises(ValueError, match="expects 1 free"):
        build_structure(216, CUBIC, [AtomSite("Li", "16e", [0.3, 0.4], 1.0)])


# --- stoichiometry / formula ----------------------------------------------

def test_stoichiometry_nacl():
    atoms = [AtomSite("Na", "4a", [], 1.0), AtomSite("Cl", "4b", [], 1.0)]
    assert stoichiometry(225, atoms) == {"Na": 4.0, "Cl": 4.0}


def test_formula_normalizes_by_gcd():
    """Na4Cl4 -> NaCl (GCD 4)."""
    assert format_formula({"Na": 4.0, "Cl": 4.0}) == "NaCl"


def test_formula_no_gcd():
    """Li6PS5Cl stays Li6PS5Cl (GCD 1)."""
    assert format_formula({"Li": 6.0, "P": 1.0, "S": 5.0, "Cl": 1.0}) == "Li6PS5Cl"


def test_formula_partial_occupancy():
    """Partial occupancy -> fractional, unnormalized."""
    assert format_formula({"Li": 2.0, "Cu": 2.0}) == "LiCu"  # integers, GCD 2
    # 16e mult 16 x occ 0.4 -> 6.4 (fractional)
    assert format_formula({"Li": 6.4}) == "Li6.4"


# --- validation ------------------------------------------------------------

def test_validate_cell_cubic_a_not_equal_b():
    """Cubic cell with a != b should warn."""
    cell = {"a": 5.0, "b": 6.0, "c": 5.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    warnings = validate_cell(225, cell)
    assert any("a should equal b" in w for w in warnings)


def test_validate_cell_clean():
    """A consistent cubic cell produces no warnings."""
    assert validate_cell(225, CUBIC) == []


def test_validate_atoms_over_occupancy():
    """Two species summing > 1.0 on one site is flagged."""
    atoms = [AtomSite("Li", "4a", [], 0.6), AtomSite("Cu", "4a", [], 0.6)]
    warnings = validate_atoms(225, atoms)
    assert any("exceeds 1.0" in w for w in warnings)


def test_validate_atoms_partial_ok():
    """Partial occupancy < 1.0 is intentional (disorder), not flagged."""
    atoms = [AtomSite("Li", "4a", [], 0.5), AtomSite("Cu", "4a", [], 0.5)]
    assert validate_atoms(225, atoms) == []


# --- CIF / XYZ round-trip --------------------------------------------------

def test_cif_round_trips_through_gemmi():
    """Built CIF reads back with the right cell, SG, and asymmetric-unit atoms."""
    atoms = [AtomSite("Na", "4a", [], 1.0), AtomSite("Cl", "4b", [], 1.0)]
    s = build_structure(225, {"a": 5.64, "b": 5.64, "c": 5.64,
                              "alpha": 90.0, "beta": 90.0, "gamma": 90.0}, atoms)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "NaCl.cif")
        write_structure_cif(s, path)
        block = gemmi.cif.read_file(path).sole_block()
        assert float(block.find_value("_cell_length_a")) == pytest.approx(5.64)
        assert "F m -3 m" in block.find_value("_symmetry_space_group_name_H-M")
        labels = list(block.find_values("_atom_site_label"))
        assert labels == ["Na1", "Cl1"]
        occ = list(block.find_values("_atom_site_occupancy"))
        assert occ == ["1.0000", "1.0000"]


def test_xyz_expands_to_full_cell():
    """XYZ writer expands the asymmetric unit to the full conventional cell."""
    atoms = [AtomSite("Na", "4a", [], 1.0), AtomSite("Cl", "4b", [], 1.0)]
    cell = {"a": 5.64, "b": 5.64, "c": 5.64, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure(225, cell, atoms)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "NaCl.xyz")
        write_structure_xyz(s, 225, path)
        lines = open(path).read().splitlines()
        # Fm-3m conventional NaCl: 4 Na + 4 Cl = 8 atoms
        assert int(lines[0]) == 8


def test_cif_read_back_by_pymatgen():
    """The built CIF is parseable by pymatgen and expands to the right count."""
    pytest.importorskip("pymatgen")
    from pymatgen.core import Structure
    atoms = [AtomSite("Na", "4a", [], 1.0), AtomSite("Cl", "4b", [], 1.0)]
    cell = {"a": 5.64, "b": 5.64, "c": 5.64, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure(225, cell, atoms)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "NaCl.cif")
        write_structure_cif(s, path)
        ps = Structure.from_file(path)
        assert ps.composition.reduced_formula == "NaCl"
        assert len(ps) == 8  # 4 Na + 4 Cl


# --- fractional-coordinate input mode -------------------------------------

def test_detect_wyckoff_canonical_and_noncanonical():
    """detect_wyckoff identifies the orbit from any representative.

    SG 137 4d is canonically (1/2,0,1/4+z); (0,1/2,0.9446) is the 4-fold-
    rotated equivalent (non-canonical). Both should detect as 4d, mult 4.
    """
    assert detect_wyckoff(137, 0.5, 0.0, 0.9446) == ("4d", 4)
    assert detect_wyckoff(137, 0.0, 0.5, 0.9446) == ("4d", 4)  # non-canonical
    assert detect_wyckoff(137, 0.0, 0.0, 0.5) == ("2b", 2)
    assert detect_wyckoff(137, 0.2463, 0.2463, 0.0) == ("8f", 8)


def test_build_structure_frac_places_atoms_at_given_coords():
    """FracAtom places the atom at exactly its fractional coordinates."""
    cell = {"a": 5.0, "b": 5.0, "c": 5.0, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure_frac(225, cell, [FracAtom("Na", 0.0, 0.0, 0.0, 1.0)])
    atom = next(a for m in s for c in m for r in c for a in r)
    frac = s.cell.fractionalize(atom.pos)
    assert (frac.x, frac.y, frac.z) == pytest.approx((0.0, 0.0, 0.0), abs=1e-4)


def test_stoichiometry_frac_uses_orbit_multiplicity():
    """stoichiometry_frac computes multiplicity by orbit expansion.

    Na on 4a (Fm-3m, mult 4) at occ 1 -> 4 Na per cell.
    """
    cell = {"a": 5.64, "b": 5.64, "c": 5.64, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    atoms = [FracAtom("Na", 0.0, 0.0, 0.0, 1.0), FracAtom("Cl", 0.5, 0.5, 0.5, 1.0)]
    assert stoichiometry_frac(225, atoms) == {"Na": 4.0, "Cl": 4.0}


def test_stoichiometry_frac_partial_occupancy():
    """Partial occupancy on a non-canonical representative counts correctly."""
    # Li0.5/Cu0.5 on 4a (Fm-3m): mult 4 * 0.5 = 2 each.
    atoms = [FracAtom("Li", 0.0, 0.0, 0.0, 0.5),
             FracAtom("Cu", 0.0, 0.0, 0.0, 0.5)]
    counts = stoichiometry_frac(225, atoms)
    assert counts == {"Li": 2.0, "Cu": 2.0}


def test_validate_atoms_frac_over_occupancy():
    atoms = [FracAtom("Li", 0.0, 0.0, 0.0, 0.6),
             FracAtom("Cu", 0.0, 0.0, 0.0, 0.6)]
    warnings = validate_atoms_frac(225, atoms)
    assert any("exceeds 1.0" in w for w in warnings)


def test_frac_cif_round_trip_lgps_li_site():
    """A disordered Li site built from fractional coords reads back correctly."""
    # LGPS-like: Li on 16h (partial) in SG 137, non-canonical rep.
    cell = {"a": 8.69, "b": 8.69, "c": 12.60, "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    atoms = [FracAtom("Li", 0.2563, 0.2718, 0.1832, 0.691)]
    s = build_structure_frac(137, cell, atoms)
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "Li.cif")
        write_structure_cif(s, path)
        block = gemmi.cif.read_file(path).sole_block()
        occ = list(block.find_values("_atom_site_occupancy"))
        assert occ == ["0.6910"]
