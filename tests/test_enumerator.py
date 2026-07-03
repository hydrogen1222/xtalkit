"""Tests for xtalkit.enumerator — symmetry-inequivalent configuration enumeration."""

import os
import tempfile

import pytest

from xtalkit._env import setup_for_enumlib
from xtalkit.enumerator import enumerate_structures, _has_pymatgen


HAS_PYMATGEN = _has_pymatgen()
skip_no_pymatgen = pytest.mark.skipif(
    not HAS_PYMATGEN,
    reason="pymatgen not installed in this environment",
)


@pytest.fixture
def disordered_cif():
    return os.path.join(os.path.dirname(__file__), "fixtures", "disordered_binary.cif")


def test_setup_idempotent():
    """Calling setup_for_enumlib twice should not raise."""
    setup_for_enumlib()
    setup_for_enumlib()
    # If we reach here, the function is idempotent
    assert True


def test_setup_adds_pathext_on_windows():
    """After setup, PATHEXT should include .X (on Windows)."""
    import sys, os
    setup_for_enumlib()
    if sys.platform == "win32":
        assert ".X" in os.environ.get("PATHEXT", "").upper()
        assert ".PY" in os.environ.get("PATHEXT", "").upper()


@skip_no_pymatgen
def test_enumerate_simple_binary(disordered_cif):
    """Real enumlib run on Au0.5/Cu0.5 — should find >=1 ordered configs."""
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        paths = enumerate_structures(
            cif_path=disordered_cif,
            min_cell_size=1,
            max_cell_size=2,
            symm_prec=0.1,
            vacancy_symbol="X",
            output_dir=out_dir,
        )
        assert len(paths) >= 1
        for p in paths:
            assert os.path.exists(p)
            assert p.endswith(".cif")


@skip_no_pymatgen
def test_enumerate_output_cifs_parse(disordered_cif):
    """Each output CIF must be re-parseable by gemmi."""
    import gemmi
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        paths = enumerate_structures(
            cif_path=disordered_cif,
            min_cell_size=1,
            max_cell_size=2,
            symm_prec=0.1,
            output_dir=out_dir,
        )
        for p in paths:
            doc = gemmi.cif.read_file(p)
            block = doc.sole_block()
            # Should have cell params
            assert float(block.find_value("_cell_length_a")) > 0
            # Should have atom sites
            atoms = list(block.find_values("_atom_site_label"))
            assert len(atoms) >= 1


@skip_no_pymatgen
def test_enumerate_max_structures_limits(disordered_cif):
    """--max-structures should limit the number of output files."""
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        paths = enumerate_structures(
            cif_path=disordered_cif,
            min_cell_size=1,
            max_cell_size=2,
            output_dir=out_dir,
            max_structures=1,
        )
        assert len(paths) == 1


@skip_no_pymatgen
def test_enumerate_xyz_format(disordered_cif):
    """--format xyz must actually write parseable XYZ files.

    pymatgen's Structure.to() does not support fmt='xyz', so enumerate has
    to route xyz output through pymatgen.io.xyz.XYZ. This test guards that
    path (a latent bug: the CLI accepted '--format xyz' but it crashed
    before the fix).
    """
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        paths = enumerate_structures(
            cif_path=disordered_cif,
            min_cell_size=1,
            max_cell_size=2,
            output_dir=out_dir,
            format="xyz",
        )
        assert len(paths) >= 1
        for p in paths:
            assert p.endswith(".xyz")
            assert os.path.exists(p)
            lines = open(p).read().strip().split("\n")
            # Line 1: atom count; line 2: comment; line 3+: "el x y z"
            assert int(lines[0]) > 0
            for line in lines[2:]:
                parts = line.split()
                assert len(parts) == 4


@skip_no_pymatgen
def test_enumerate_nonexistent_cif_raises():
    """Missing CIF should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        enumerate_structures(cif_path="/nonexistent/path.cif")


def test_enumerate_missing_pymatgen_message(monkeypatch):
    """When pymatgen is unavailable, raise RuntimeError with install instructions."""
    # Force _has_pymatgen to return False
    import xtalkit.enumerator as enum_mod
    monkeypatch.setattr(enum_mod, "_has_pymatgen", lambda: False)
    # Need a CIF that exists so we don't trip FileNotFoundError first
    cif = os.path.join(os.path.dirname(__file__), "fixtures", "disordered_binary.cif")
    with pytest.raises(RuntimeError) as exc_info:
        enum_mod.enumerate_structures(cif_path=cif)
    msg = str(exc_info.value)
    assert "pymatgen" in msg.lower()
    assert "uv" in msg.lower()


@skip_no_pymatgen
def test_partial_occupancy_augmentation():
    """Pure-Python test: a 0.5-occupancy site should be augmented with DummySpecies."""
    from pymatgen.core import Structure, Lattice, DummySpecies
    lat = Lattice.cubic(4.0)
    struct = Structure(lat, [{"Au": 0.5, "Cu": 0.5}], [[0, 0, 0]])
    # Manually augment using the function under test
    from xtalkit.enumerator import _augment_partial_occupancy
    n = _augment_partial_occupancy(struct, "X")
    # The mixed-species site should NOT be augmented (already mixed)
    assert n == 0

    # Now test a single-species partial-occupancy site
    struct2 = Structure(lat, ["Li"], [[0, 0, 0]])
    # Set occupancy manually
    from pymatgen.core import Composition
    struct2[0] = {"Li": 0.5}
    n2 = _augment_partial_occupancy(struct2, "X")
    assert n2 == 1
    # The site should now have Li + X
    site = struct2[0]
    assert len(site.species) == 2
    assert DummySpecies("X") in site.species
