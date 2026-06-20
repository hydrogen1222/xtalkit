import pytest
import gemmi
from xtalkit.matcher import match_atoms
from xtalkit.spacegroup import wyckoff_positions


def _make_cubic_structure(atoms_info: list[tuple[str, tuple[float,float,float]]]) -> gemmi.Structure:
    """Build a gemmi.Structure with cubic cell a=5.0."""
    s = gemmi.Structure()
    cell = gemmi.UnitCell(5.0, 5.0, 5.0, 90, 90, 90)
    s.cell = cell
    # Space group 216 F-43m
    s.spacegroup_hm = gemmi.SpaceGroup(216).hm
    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")
    for i, (element, frac) in enumerate(atoms_info):
        atom = gemmi.Atom()
        atom.name = f"{element}{i+1}"
        atom.element = gemmi.Element(element)
        # Convert fractional to Cartesian via cell
        atom.pos = cell.orthogonalize(gemmi.Fractional(*frac))
        residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    s.add_model(model)
    return s


def test_match_atoms_exact_position():
    """Atom at (0,0,0) should match 4a in F-43m."""
    structure = _make_cubic_structure([("Li", (0.0, 0.0, 0.0))])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.5)
    assert "Li1" in mapping
    assert mapping["Li1"] == "4a"


def test_match_atoms_near_position():
    """Atom near (0.25,0.25,0.25) should match 4c in F-43m (0.25,0.25,0.25)."""
    structure = _make_cubic_structure([("P", (0.251, 0.249, 0.250))])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.5)
    assert "P1" in mapping
    assert mapping["P1"] == "4c"


def test_match_atoms_no_match_beyond_tolerance():
    """Atom at (0.0, 0.0, 0.0) with tiny tolerance 0.01 on a SG with no
    Wyckoff at exactly (0,0,0) — but SG 216 does have 4a at (0,0,0).
    Let's use a position far from any Wyckoff: (0.99, 0.99, 0.99)."""
    structure = _make_cubic_structure([("Li", (0.99, 0.99, 0.99))])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.01)
    assert "Li1" not in mapping


def test_match_atoms_multiple_atoms():
    """Multiple atoms should each get matched."""
    structure = _make_cubic_structure([
        ("Li", (0.0, 0.0, 0.0)),
        ("P", (0.25, 0.25, 0.25)),
    ])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.5)
    assert mapping["Li1"] == "4a"
    assert mapping["P2"] == "4c"


def test_match_atoms_empty_structure():
    """Empty structure returns empty dict."""
    structure = _make_cubic_structure([])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.5)
    assert mapping == {}


def test_match_atoms_symmetry_equivalent():
    """Atom at (0.5,0.5,0.5) should match 4b in F-43m.
    This tests that symmetry-equivalent positions are generated correctly."""
    structure = _make_cubic_structure([("Na", (0.5, 0.5, 0.5))])
    wyckoffs = wyckoff_positions(216)
    mapping = match_atoms(structure, wyckoffs, tolerance=0.5)
    assert "Na1" in mapping
    assert mapping["Na1"] == "4b"
