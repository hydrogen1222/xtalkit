"""Enumerate symmetry-inequivalent ordered configurations of a disordered CIF.

Wraps pymatgen's ``EnumlibAdaptor`` (which shells out to the Fortran
``enum.x`` and ``makestr.x`` binaries from msg-byu/enumlib). Requires the
``enumerate`` uv extra (``uv sync --extra enumerate``) plus source-compiled
enumlib binaries (``scripts/build_enumlib.sh``) — see the README
"Enumeration setup" section. Lazy-imports pymatgen so the rest of xtalkit
works without it.
"""

from __future__ import annotations

import os
import tempfile

from xtalkit._env import setup_for_enumlib


def _has_pymatgen() -> bool:
    """Check whether pymatgen is importable. Must NOT import enumlib_caller
    here — that module's module-level code calls ``which("enum.x")`` at
    import time, which would cache a None result before setup_for_enumlib()
    has had a chance to put the binaries on PATH."""
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
            "Compile them with: bash scripts/build_enumlib.sh "
            "(or set XTALKIT_ENUMLIB_BIN to a directory containing enum.x and "
            "makestr.x). See README 'Enumeration setup'."
        )
    return None


def _augment_partial_occupancy(struct, vacancy_symbol: str) -> int:
    """Replace partial-occupancy sites with an explicit vacancy species.

    enumlib's EnumlibAdaptor needs an explicit vacancy species to enumerate
    orderings of a partially-occupied site. Handles both single- and
    multi-species partial-occupancy sites by adding ``DummySpecies`` for the
    shortfall:

    - ``{Li: 0.5}``            -> ``{Li: 0.5, X: 0.5}``
    - ``{Au: 0.3, Cu: 0.3}``   -> ``{Au: 0.3, Cu: 0.3, X: 0.4}``

    Sites whose species sum to 1.0 (fully occupied, e.g. ``{Au:0.5,Cu:0.5}``)
    are left untouched. Returns the count of augmented sites.
    """
    from pymatgen.core import DummySpecies
    augmented = 0
    for i, site in enumerate(struct):
        species = site.species
        total = sum(species.values())
        if abs(total - 1.0) <= 0.01:
            continue
        new_species = dict(species)
        new_species[DummySpecies(vacancy_symbol)] = 1.0 - total
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
    enumlib binaries are unavailable, if enumeration times out, or if
    enumlib returns 0 structures (e.g. non-integer occupancy with too small
    a max_cell_size).
    """
    resolved = os.path.abspath(cif_path)
    if not os.path.exists(resolved):
        raise FileNotFoundError(f"CIF not found: {resolved}")

    if not _has_pymatgen():
        raise RuntimeError(
            "pymatgen is required for `xtalkit enumerate` but is not installed.\n"
            "Install the optional extra:  uv sync --extra enumerate\n"
            "Then compile enumlib binaries: bash scripts/build_enumlib.sh\n"
            "See README 'Enumeration setup'."
        )

    setup_for_enumlib()

    bin_err = _check_enumlib_binaries()
    if bin_err:
        raise RuntimeError(bin_err)

    from pymatgen.core import Structure
    try:
        from pymatgen.command_line.enumlib_caller import (
            EnumlibAdaptor, EnumError)
    except ImportError as e:
        raise RuntimeError(
            "pymatgen is installed but its enumlib_caller module is missing "
            "(removed in some versions). Reinstall with a version that "
            "includes it, e.g.: uv pip install 'pymatgen>=2024.5'."
        ) from e

    struct = Structure.from_file(resolved)
    # Convert to primitive cell before enumeration. enumlib does not do this
    # itself, and a conventional cell on an F/I-centered space group can
    # have 2x-4x the candidate sites of the primitive cell — enough to
    # overflow enumlib's tree_class array (e.g. F-43m 48h x 48 candidates
    # -> C(48,24) overflows; the primitive cell has only 12 candidates ->
    # C(12,6) = 924, which is exactly the regime that produces the ~48
    # symmetry-inequivalent Li orderings reported in the argyrodite
    # literature).
    struct = struct.get_primitive_structure()
    n_augmented = _augment_partial_occupancy(struct, vacancy_symbol)

    basename = os.path.splitext(os.path.basename(resolved))[0]
    # Resolve to an absolute path *before* changing directory below, so the
    # output lands in the user's CWD rather than the throwaway scratch dir.
    out_dir = os.path.abspath(output_dir or f"{basename}_enum")
    os.makedirs(out_dir, exist_ok=True)

    adaptor = EnumlibAdaptor(
        struct,
        min_cell_size=min_cell_size,
        max_cell_size=max_cell_size,
        symm_prec=symm_prec,
        timeout=timeout,
    )
    try:
        # enumlib_caller writes struct_enum.in/.out and vasp.* to the CWD;
        # run inside a throwaway directory so the user's CWD stays clean.
        prev_cwd = os.getcwd()
        with tempfile.TemporaryDirectory() as scratch:
            os.chdir(scratch)
            try:
                adaptor.run()
            finally:
                os.chdir(prev_cwd)
    except EnumError as e:
        raise RuntimeError(
            f"enumlib returned 0 structures: {e}\n"
            f"Likely cause: non-integer stoichiometry at max_cell_size={max_cell_size}.\n"
            "Try a larger --max-cell-size, or provide a 'clean' parent CIF "
            "with integer stoichiometry (e.g., round 0.56 -> 0.5)."
        ) from e
    except TimeoutError as e:
        raise RuntimeError(
            f"enumlib exceeded the {timeout}-minute timeout.\n"
            "Increase --timeout, reduce --max-cell-size, or use a cleaner "
            "parent CIF."
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
