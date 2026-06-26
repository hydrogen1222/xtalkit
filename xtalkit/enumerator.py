"""Enumerate symmetry-inequivalent ordered configurations of a disordered CIF.

Wraps pymatgen's ``EnumlibAdaptor`` (which shells out to the Fortran
``enum.x`` and ``makestr.x`` binaries from msg-byu/enumlib). Requires a
conda env with pymatgen 2023.x and self-compiled enumlib binaries — see
README "Enumeration" section for setup. Lazy-imports pymatgen so that the
rest of xtalkit works without it.
"""

from __future__ import annotations

import os
import sys

from xtalkit._env import setup_for_enumlib


def _has_pymatgen() -> bool:
    """Check whether pymatgen is importable. Must NOT import enumlib_caller
    here — that module's module-level code calls ``which("enum.x")`` at
    import time, which would cache a None result before setup_for_enumlib()
    has a chance to patch PATHEXT."""
    try:
        import pymatgen  # noqa: F401
        return True
    except ImportError:
        return False


def _check_enumlib_binaries() -> str | None:
    """Return error message if enum.x or makestr.x not found, else None."""
    from shutil import which
    missing = []
    if not which("enum.x"):
        missing.append("enum.x")
    if not which("makestr.x"):
        missing.append("makestr.x")
    if missing:
        return (
            f"enumlib binaries not found in PATH: {', '.join(missing)}. "
            "Compile them from https://github.com/msg-byu/enumlib "
            "(see README 'Enumeration' section)."
        )
    return None


def _augment_partial_occupancy(struct, vacancy_symbol: str) -> int:
    """Replace partial-occupancy single-species sites with explicit
    species + DummySpecies(vacancy). Returns count of augmented sites.

    pymatgen's EnumlibAdaptor needs an explicit vacancy species to enumerate
    orderings of a partially-occupied site.
    """
    from pymatgen.core import DummySpecies
    augmented = 0
    for i, site in enumerate(struct):
        if len(site.species) != 1:
            continue
        sp, occ = next(iter(site.species.items()))
        if abs(occ - 1.0) <= 0.01:
            continue
        new_species = {sp: occ, DummySpecies(vacancy_symbol): 1.0 - occ}
        struct[i] = new_species
        augmented += 1
    return augmented


def enumerate_structures(
    cif_path: str,
    min_cell_size: int = 1,
    max_cell_size: int = 2,
    symm_prec: float = 0.1,
    vacancy_symbol: str = "X",
    output_dir: str | None = None,
    max_structures: int | None = None,
    timeout: float | None = None,
    format: str = "cif",
) -> list[str]:
    """Enumerate ordered configurations of a disordered CIF.

    Returns list of output file paths. Raises RuntimeError if pymatgen or
    enumlib binaries are unavailable. Raises the underlying pymatgen
    EnumError if enumeration fails (e.g., non-integer occupancy with too
    small max_cell_size).
    """
    resolved = os.path.abspath(cif_path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"CIF not found: {resolved}")

    if not _has_pymatgen():
        raise RuntimeError(
            "pymatgen is required for `xtalkit enumerate` but is not installed.\n"
            "Install in a conda env:\n"
            "  conda create -n xtalkit -c conda-forge python=3.11 m2w64-gcc-fortran make git\n"
            "  conda activate xtalkit\n"
            "  conda install -c conda-forge pymatgen==2023.5.31\n"
            "Then compile enumlib (see README 'Enumeration' section)."
        )

    setup_for_enumlib()

    bin_err = _check_enumlib_binaries()
    if bin_err:
        raise RuntimeError(bin_err)

    from pymatgen.core import Structure
    from pymatgen.command_line.enumlib_caller import EnumlibAdaptor, EnumError

    struct = Structure.from_file(resolved)
    # Convert to primitive cell before enumeration. enumlib does not do this
    # itself, and a conventional cell on an F/I-centered space group can
    # have 2×–4× the candidate sites of the primitive cell — enough to
    # overflow enumlib's tree_class array (e.g. F-43m 48h × 48 candidates
    # → C(48,24) overflows; the primitive cell has only 12 candidates →
    # C(12,6) = 924, which is exactly the regime that produces the ~48
    # symmetry-inequivalent Li orderings reported in the argyrodite
    # literature).
    struct = struct.get_primitive_structure()
    n_augmented = _augment_partial_occupancy(struct, vacancy_symbol)

    basename = os.path.splitext(os.path.basename(resolved))[0]
    out_dir = output_dir or f"{basename}_enum"
    os.makedirs(out_dir, exist_ok=True)

    try:
        adaptor = EnumlibAdaptor(
            struct,
            min_cell_size=min_cell_size,
            max_cell_size=max_cell_size,
            symm_prec=symm_prec,
            timeout=timeout,
        )
        adaptor.run()
    except EnumError as e:
        raise RuntimeError(
            f"enumlib returned 0 structures: {e}\n"
            f"Likely cause: non-integer stoichiometry at max_cell_size={max_cell_size}.\n"
            "Try a larger --max-cell-size, or provide a 'clean' parent CIF "
            "with integer stoichiometry (e.g., round 0.56 -> 0.5)."
        ) from e

    structures = adaptor.structures
    if max_structures is not None:
        structures = structures[:max_structures]

    paths: list[str] = []
    for i, s in enumerate(structures):
        path = os.path.join(out_dir, f"{basename}_{i:03d}.{format}")
        s.to(filename=path)
        paths.append(path)

    return paths
