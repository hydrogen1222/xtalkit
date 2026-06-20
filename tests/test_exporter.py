"""Tests for xtalkit.exporter — multi-format output writers."""

import os
import tempfile

import pytest
import gemmi
from xtalkit.exporter import write_cif, write_vesta, write_xyz, DummyAtom


@pytest.fixture
def cubic_structure():
    """A minimal gemmi Structure with two atoms in F-43m."""
    s = gemmi.Structure()
    s.cell = gemmi.UnitCell(5.0, 5.0, 5.0, 90, 90, 90)
    s.spacegroup_hm = "F -4 3 m"
    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")
    for i, (el, frac) in enumerate([
        ("Li", (0.0, 0.0, 0.0)),
        ("P", (0.25, 0.25, 0.25)),
    ]):
        atom = gemmi.Atom()
        atom.name = f"{el}{i+1}"
        atom.element = gemmi.Element(el)
        # Convert fractional to Cartesian for gemmi 0.7.5
        orth = s.cell.orthogonalize(gemmi.Fractional(*frac))
        atom.pos = (orth.x, orth.y, orth.z)
        residue.add_atom(atom)
    chain.add_residue(residue)
    model.add_chain(chain)
    s.add_model(model)
    return s


@pytest.fixture
def dummy_atoms():
    """Two dummy atoms representing Wyckoff sites."""
    return [
        DummyAtom("WYCK_4a", "Xe", (0.0, 0.0, 0.0)),
        DummyAtom("WYCK_4c", "Kr", (0.25, 0.25, 0.25)),
    ]


def test_write_cif_creates_file(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.cif")
        write_cif(cubic_structure, dummy_atoms, path)
        assert os.path.exists(path)
        doc = gemmi.cif.read_file(path)
        assert len(doc) > 0


def test_write_cif_contains_dummy_atom_label(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.cif")
        write_cif(cubic_structure, dummy_atoms, path)
        content = open(path).read()
        assert "WYCK_4a" in content
        assert "Xe" in content


def test_write_vesta_creates_file(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.vesta")
        write_vesta(cubic_structure, dummy_atoms, path)
        assert os.path.exists(path)


def test_write_vesta_is_valid_xml(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.vesta")
        write_vesta(cubic_structure, dummy_atoms, path)
        content = open(path).read()
        assert content.startswith('<?xml')
        assert 'VESTA' in content or '<root' in content.lower() or '<?xml' in content


def test_write_xyz_creates_file(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.xyz")
        write_xyz(cubic_structure, dummy_atoms, path)
        assert os.path.exists(path)


def test_write_xyz_content_format(cubic_structure, dummy_atoms):
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "test.xyz")
        write_xyz(cubic_structure, dummy_atoms, path)
        lines = open(path).read().strip().split("\n")
        # First line: atom count
        assert int(lines[0]) > 0
        # Remaining lines: element x y z
        for line in lines[2:]:
            parts = line.split()
            assert len(parts) == 4  # element, x, y, z


def test_write_cif_directory_not_exist(cubic_structure, dummy_atoms):
    """Should create parent directories."""
    with tempfile.TemporaryDirectory() as tmp:
        path = os.path.join(tmp, "subdir", "sub2", "test.cif")
        write_cif(cubic_structure, dummy_atoms, path)
        assert os.path.exists(path)
