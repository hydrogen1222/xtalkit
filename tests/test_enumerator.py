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
def test_enumerate_custom_vacancy_symbol_is_removed(tmp_path):
    """A non-default vacancy symbol should not leak into output structures."""
    from pymatgen.core import Structure
    from xtalkit.builder import build_structure_frac, FracAtom
    from xtalkit.exporter import write_structure_cif

    cell = {"a": 4.0, "b": 4.0, "c": 4.0,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure_frac(225, cell, [FracAtom("Li", 0, 0, 0, 0.5)])
    cif = str(tmp_path / "li_vacancy.cif")
    write_structure_cif(s, cif)

    paths = enumerate_structures(
        cif_path=cif,
        max_cell_size=2,
        output_dir=str(tmp_path / "out"),
        vacancy_symbol="Xa",
    )

    assert paths
    for path in paths:
        symbols = {sp.symbol for sp in Structure.from_file(path).composition}
        assert "Xa" not in symbols


@skip_no_pymatgen
def test_enum_input_keeps_distinct_li_vacancy_orbits_separate(tmp_path):
    """Same-species disorder on different Wyckoff orbits needs distinct labels."""
    from pymatgen.core import Structure
    from xtalkit.builder import build_structure_frac, FracAtom
    from xtalkit.exporter import write_structure_cif
    from xtalkit.enumerator import (
        _augment_partial_occupancy,
        _gen_input_file_preserve_site_constraints,
    )

    cell = {"a": 8.69407, "b": 8.69407, "c": 12.5994,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    built = build_structure_frac(137, cell, [
        FracAtom("Li", 0.2563, 0.2718, 0.1832, 0.6875),  # 16h: 11 Li + 5 X
        FracAtom("Li", 0.0, 0.5, 0.9446, 1.0),
        FracAtom("Li", 0.2463, 0.2463, 0.0, 0.625),     # 8f: 5 Li + 3 X
        FracAtom("Ge", 0.0, 0.5, 0.6907, 0.5),
        FracAtom("P", 0.0, 0.5, 0.6907, 0.5),
        FracAtom("P", 0.0, 0.0, 0.5, 1.0),
        FracAtom("S", 0.0, 0.1843, 0.4103, 1.0),
        FracAtom("S", 0.0, 0.2991, 0.0950, 1.0),
        FracAtom("S", 0.0, 0.6990, 0.7914, 1.0),
    ])
    cif = str(tmp_path / "two_li_orbits.cif")
    write_structure_cif(built, cif)
    struct = Structure.from_file(cif).get_primitive_structure()
    _augment_partial_occupancy(struct, "X")

    setup_for_enumlib()
    from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
    adaptor = EnumlibAdaptor(struct, min_cell_size=1, max_cell_size=1, symm_prec=0.1)

    prev = os.getcwd()
    os.chdir(tmp_path)
    try:
        _gen_input_file_preserve_site_constraints(adaptor)
        text = (tmp_path / "struct_enum.in").read_text()
    finally:
        os.chdir(prev)

    assert [str(sp) for sp in adaptor.index_species] == [
        "Li", "X0+", "Li", "X0+", "Ge", "P",
    ]
    assert text.count(" 0/1") == 16
    assert text.count(" 2/3") == 8
    assert text.count(" 4/5") == 4
    assert "110 110 280" in text  # 11 Li on the 16h orbit
    assert "50 50 280" in text    # 5 vacancies on 16h and 5 Li on 8f
    assert "30 30 280" in text    # 3 vacancies on the 8f orbit
    assert "20 20 280" in text    # 2 Ge and 2 P on the 4d orbit


@skip_no_pymatgen
def test_enumerate_nonexistent_cif_raises():
    """Missing CIF should raise FileNotFoundError."""
    with pytest.raises(FileNotFoundError):
        enumerate_structures(cif_path="/nonexistent/path.cif")


@skip_no_pymatgen
def test_enumerate_fully_ordered_cif_returns_empty(tmp_path):
    """A fully ordered CIF has no occupational disorder to enumerate."""
    from xtalkit.builder import build_structure_frac, FracAtom
    from xtalkit.exporter import write_structure_cif

    cell = {"a": 5.64, "b": 5.64, "c": 5.64,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    nacl = build_structure_frac(225, cell, [
        FracAtom("Na", 0.0, 0.0, 0.0, 1.0),
        FracAtom("Cl", 0.5, 0.5, 0.5, 1.0),
    ])
    cif = str(tmp_path / "nacl.cif")
    write_structure_cif(nacl, cif)

    assert enumerate_structures(
        cif_path=cif,
        output_dir=str(tmp_path / "out"),
    ) == []


def test_enumerate_rejects_invalid_parameters(disordered_cif):
    """Bad API parameters should fail before enumlib is launched."""
    bad_kwargs = [
        {"min_cell_size": 0},
        {"min_cell_size": 3, "max_cell_size": 2},
        {"max_structures": 0},
        {"jobs": -1},
        {"batch_size": 0},
        {"format": "poscar"},
        {"vacancy_symbol": ""},
    ]
    for kwargs in bad_kwargs:
        with pytest.raises(ValueError):
            enumerate_structures(cif_path=disordered_cif, **kwargs)


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


def _reduced_formulas(out_dir):
    """Sorted reduced formulas of every CIF in out_dir (set comparison helper)."""
    import glob
    from pymatgen.core import Structure
    formulas = []
    for p in sorted(glob.glob(os.path.join(out_dir, "*.cif"))):
        formulas.append(Structure.from_file(p).composition.reduced_formula)
    return sorted(formulas)


def _pymatgen_reference_formulas(cif_path, min_cell=1, max_cell=2):
    """Run pymatgen's EnumlibAdaptor.run() path directly for comparison."""
    from pymatgen.core import Structure
    from pymatgen.command_line.enumlib_caller import EnumlibAdaptor
    from xtalkit._env import setup_for_enumlib
    from xtalkit.enumerator import _augment_partial_occupancy
    setup_for_enumlib()
    struct = Structure.from_file(cif_path).get_primitive_structure()
    _augment_partial_occupancy(struct, "X")
    adaptor = EnumlibAdaptor(struct, min_cell_size=min_cell,
                             max_cell_size=max_cell, symm_prec=0.1)
    with tempfile.TemporaryDirectory() as scratch:
        prev = os.getcwd()
        os.chdir(scratch)
        try:
            adaptor._gen_input_file()
            n = adaptor._run_multienum()
            ref = adaptor._get_structures(n) if n > 0 else []
        finally:
            os.chdir(prev)
    return sorted(s.composition.reduced_formula for s in ref)


@skip_no_pymatgen
def test_streaming_matches_pymatgen_adaptor(disordered_cif):
    """Streamed batched output is the same SET of structures as pymatgen's
    EnumlibAdaptor. This validates the vasp->Structure parsing replication."""
    ref = _pymatgen_reference_formulas(disordered_cif, max_cell=2)
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        enumerate_structures(
            cif_path=disordered_cif, min_cell_size=1, max_cell_size=2,
            output_dir=out_dir, batch_size=1,  # force many small batches
        )
        assert _reduced_formulas(out_dir) == ref


@skip_no_pymatgen
def test_batch_size_does_not_change_result(disordered_cif):
    """The set of structures is independent of batch_size."""
    with tempfile.TemporaryDirectory() as tmp:
        out_a = os.path.join(tmp, "a")
        out_b = os.path.join(tmp, "b")
        enumerate_structures(cif_path=disordered_cif, max_cell_size=2,
                             output_dir=out_a, batch_size=1)
        enumerate_structures(cif_path=disordered_cif, max_cell_size=2,
                             output_dir=out_b, batch_size=256)
        assert _reduced_formulas(out_a) == _reduced_formulas(out_b)


@skip_no_pymatgen
def test_parallel_matches_serial(disordered_cif):
    """--jobs=2 produces the same set of structures as --jobs=1."""
    with tempfile.TemporaryDirectory() as tmp:
        out_s = os.path.join(tmp, "serial")
        out_p = os.path.join(tmp, "parallel")
        enumerate_structures(cif_path=disordered_cif, max_cell_size=2,
                             output_dir=out_s, jobs=1)
        enumerate_structures(cif_path=disordered_cif, max_cell_size=2,
                             output_dir=out_p, jobs=2)
        assert _reduced_formulas(out_s) == _reduced_formulas(out_p)


@skip_no_pymatgen
def test_max_structures_does_not_exceed_nor_invent(disordered_cif):
    """--max-structures caps output; a cap larger than the natural count
    returns the natural count (no padding, no truncation beyond the cap)."""
    ref_count = len(_pymatgen_reference_formulas(disordered_cif, max_cell=2))
    with tempfile.TemporaryDirectory() as tmp:
        out_dir = os.path.join(tmp, "out")
        paths = enumerate_structures(
            cif_path=disordered_cif, max_cell_size=2,
            output_dir=out_dir, max_structures=ref_count + 100,
        )
        assert len(paths) == ref_count


@skip_no_pymatgen
def test_non_integerizable_species_flags_bad_occupancy():
    """_non_integerizable_species flags 0.56 (1*0.56=0.56, never integer)."""
    from pymatgen.core import Structure, Lattice
    from xtalkit.enumerator import _non_integerizable_species, _augment_partial_occupancy
    struct = Structure(Lattice.cubic(4.0), ["Li"], [[0, 0, 0]])
    struct[0] = {"Li": 0.56}
    _augment_partial_occupancy(struct, "X")
    bad = _non_integerizable_species(struct, 0.1, 1, 2)
    assert len(bad) == 1
    assert bad[0]["species"] == "Li"
    assert bad[0]["count"] == 0.56


@skip_no_pymatgen
def test_non_integerizable_species_passes_clean_occupancy():
    """0.5 on a mult-1 site is integerizable at cell_size 2 (2*0.5=1)."""
    from pymatgen.core import Structure, Lattice
    from xtalkit.enumerator import _non_integerizable_species, _augment_partial_occupancy
    struct = Structure(Lattice.cubic(4.0), ["Li"], [[0, 0, 0]])
    struct[0] = {"Li": 0.5}
    _augment_partial_occupancy(struct, "X")
    assert _non_integerizable_species(struct, 0.1, 1, 2) == []


@skip_no_pymatgen
def test_enumerate_refuses_non_integer_occupancy(tmp_path):
    """enumerate_structures refuses a non-integer-occupancy CIF before enum.x
    runs (instead of letting enumlib exhaust memory and crash)."""
    from xtalkit.builder import build_structure_frac, FracAtom
    from xtalkit.exporter import write_structure_cif
    cell = {"a": 4.0, "b": 4.0, "c": 4.0,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure_frac(225, cell, [FracAtom("Li", 0, 0, 0, 0.56)])
    bad_cif = str(tmp_path / "bad.cif")
    write_structure_cif(s, bad_cif)
    with pytest.raises(RuntimeError, match="non-integer stoichiometry"):
        enumerate_structures(
            cif_path=bad_cif, max_cell_size=1,
            output_dir=str(tmp_path / "out"),
        )


@skip_no_pymatgen
def test_enumerate_skip_preflight_bypasses_check(disordered_cif, tmp_path):
    """--skip-preflight bypasses the integerizability check (escape hatch)."""
    # A clean CIF passes the check anyway; this just confirms the flag is wired
    # through and doesn't raise for a non-integer case when skipped. Use a bad
    # CIF and confirm skip_preflight lets enumlib run (returns 0 cleanly here
    # because 0.56 can't integerize, but the pre-flight is what we're skipping).
    from xtalkit.builder import build_structure_frac, FracAtom
    from xtalkit.exporter import write_structure_cif
    cell = {"a": 4.0, "b": 4.0, "c": 4.0,
            "alpha": 90.0, "beta": 90.0, "gamma": 90.0}
    s = build_structure_frac(225, cell, [FracAtom("Li", 0, 0, 0, 0.56)])
    bad_cif = str(tmp_path / "bad.cif")
    write_structure_cif(s, bad_cif)
    # With skip_preflight, the pre-flight RuntimeError is NOT raised; enumlib
    # then runs and returns 0 (caught as RuntimeError "0 structures").
    with pytest.raises(RuntimeError, match="0 structures"):
        enumerate_structures(
            cif_path=bad_cif, max_cell_size=1,
            output_dir=str(tmp_path / "out"), skip_preflight=True,
        )
