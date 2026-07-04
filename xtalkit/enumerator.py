"""Enumerate symmetry-inequivalent ordered configurations of a disordered CIF.

Wraps pymatgen's ``EnumlibAdaptor`` (which shells out to the Fortran
``enum.x`` and ``makestr.x`` binaries from msg-byu/enumlib). Requires the
``enumerate`` uv extra (``uv sync --extra enumerate``) plus source-compiled
enumlib binaries (``scripts/build_enumlib.sh``) — see the README
"Enumeration setup" section. Lazy-imports pymatgen so the rest of xtalkit
works without it.

Memory/disk/CPU notes
---------------------
``enum.x`` (the enumeration core) is single-threaded Fortran and cannot be
parallelised here — it is often the dominant cost. What we *can* and do
improve vs. calling ``EnumlibAdaptor.run()`` directly:

* **Streaming** — ``EnumlibAdaptor`` runs ``makestr.x`` for *all* N structures
  at once, then holds all N parsed ``Structure`` objects in memory. We instead
  process structures in batches: ``makestr.x`` is invoked per batch (disjoint
  index ranges), each batch is parsed and written to disk, and the batch's
  ``vasp.*`` intermediates are deleted before the next batch. Memory peaks at
  one batch, not N structures.
* **Cap at generation** — ``--max-structures`` is passed to the ``makestr.x``
  range, so the limit actually reduces work (the default wrapper generates all
  N then slices).
* **Parallel structure generation** — with ``--jobs > 1``, batches run in
  parallel worker processes (``makestr.x`` + parse + write CIF). This
  parallelises the post-enumeration phase only; ``enum.x`` stays serial.
* **Ramdisk scratch** — ``--scratch-dir /dev/shm`` puts ``struct_enum.out``
  and the ``vasp.*`` intermediates on tmpfs to avoid disk I/O.

This relies on two pymatgen-internal methods (``_gen_input_file`` and
``_run_multienum``) and the ``index_species``/``ordered_sites`` attributes
they populate; the vasp→Structure parsing is copied faithfully from
``EnumlibAdaptor._get_structures`` so the output matches.
"""

from __future__ import annotations

import glob
import os
import re
import subprocess
import tempfile
from concurrent.futures import ProcessPoolExecutor

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


def _non_integerizable_species(struct, symm_prec: float, min_cell: int,
                               max_cell: int) -> list[dict]:
    """Find species whose cell count is not integerizable.

    For each symmetry-distinct site, the count of a species in a kx supercell
    is ``multiplicity * k * occupancy``. enumlib requires this to be an integer
    for some k in [min_cell, max_cell]. If no such k exists, enumlib cannot
    place an integer number of that species and may run away consuming memory
    (it does not always fail cleanly — see the crash this guards against).

    Returns a list of ``{species, mult, occ, count, suggested}`` dicts for the
    offending REAL species (vacancies are derived, so not reported). Each
    ``suggested`` is the nearest valid fraction at ``min_cell`` (i.e.
    ``round(mult*min_cell*occ) / (mult*min_cell)``).
    """
    from pymatgen.core import DummySpecies
    from pymatgen.symmetry.analyzer import SpacegroupAnalyzer
    try:
        sym = SpacegroupAnalyzer(struct, symprec=symm_prec).get_symmetrized_structure()
    except Exception:
        return []  # don't block on a symmetry-analysis failure

    bad: list[dict] = []
    seen: set[tuple] = set()
    for group in sym.equivalent_sites:
        mult = len(group)
        rep = group[0]
        for sp, occ in rep.species.items():
            if isinstance(sp, DummySpecies):
                continue  # vacancy is derived from the real species' occupancy
            key = (str(sp), mult, round(float(occ), 6))
            if key in seen:
                continue
            seen.add(key)
            counts = [mult * k * float(occ) for k in range(min_cell, max_cell + 1)]
            if not any(abs(c - round(c)) < 0.01 for c in counts):
                j = max(0, min(round(mult * min_cell * float(occ)),
                               mult * min_cell))
                suggested = j / (mult * min_cell)
                bad.append({
                    "species": str(sp), "mult": mult, "occ": float(occ),
                    "count": mult * float(occ), "suggested": suggested,
                })
    return bad


def _makestr_range_args(struct_enum_out: str, start: int, end_exclusive: int) -> list[str]:
    """Build the range arguments for makestr.x for indices [start, end_exclusive).

    enumlib's Fortran ``makestr.x`` takes (file, start, end) with end INCLUSIVE
    (pymatgen calls it with ``0 num_structs-1``). The Python fallback
    ``makeStr.py`` takes ``-input file 1 N`` (1-indexed inclusive).
    """
    from shutil import which
    makestr = which("makestr.x") or which("makeStr.x") or which("makeStr.py")
    if makestr and ".py" in makestr:
        # 1-indexed inclusive: [start, end_exclusive) -> (start+1, end_exclusive)
        return [makestr, "-input", struct_enum_out,
                str(start + 1), str(end_exclusive)]
    # Fortran: 0-indexed inclusive: [start, end_exclusive) -> (start, end_exclusive-1)
    return [makestr, struct_enum_out, str(start), str(end_exclusive - 1)]


def _vasp_numeric_key(path: str) -> int:
    """Sort key for vasp.* files by their numeric suffix (vasp.2 before vasp.10)."""
    try:
        return int(os.path.basename(path).split(".")[-1])
    except ValueError:
        return 0


def _parse_vasp_to_structure(data: str, index_species, ordered_sites) -> "Structure":
    """Parse one vasp-format string into a Structure.

    Faithful copy of pymatgen's ``EnumlibAdaptor._get_structures`` per-file
    logic: regex-fix the POSCAR, map the enumerated lattice back to the
    original ordered sites' lattice (supercell construction), drop vacancies
    (species X). ``index_species`` and ``ordered_sites`` come from the
    adaptor (set by ``_gen_input_file``).
    """
    import numpy as np
    from pymatgen.core import PeriodicSite, Structure
    from pymatgen.io.vasp.inputs import Poscar

    data = re.sub(r"scale factor", "1", data)
    data = re.sub(r"(\d+)-(\d+)", r"\1 -\2", data)
    poscar = Poscar.from_str(data, index_species)
    sub_structure = poscar.structure
    new_latt = sub_structure.lattice

    disordered_site_properties: dict = {}
    sites: list = []

    if len(ordered_sites) > 0:
        original_latt = ordered_sites[0].lattice
        site_properties: dict = {}
        for site in ordered_sites:
            for k, v in site.properties.items():
                disordered_site_properties[k] = None
                if k in site_properties:
                    site_properties[k].append(v)
                else:
                    site_properties[k] = [v]
        ordered_structure = Structure(
            original_latt,
            [site.species for site in ordered_sites],
            [site.frac_coords for site in ordered_sites],
            site_properties=site_properties,
        )
        inv_org_latt = np.linalg.inv(original_latt.matrix)
        transformation = np.dot(new_latt.matrix, inv_org_latt)
        transformation = [[round(cell) for cell in row] for row in transformation]
        struct = ordered_structure * transformation
        sites.extend([site.to_unit_cell() for site in struct])
        super_latt = sites[-1].lattice
    else:
        super_latt = new_latt

    for site in sub_structure:
        if site.specie.symbol != "X":  # exclude vacancies
            sites.append(
                PeriodicSite(
                    site.species,
                    site.frac_coords,
                    super_latt,
                    to_unit_cell=True,
                    properties=disordered_site_properties,
                )
            )
    return Structure.from_sites(sorted(sites))


def _process_chunk(args) -> list[tuple[int, str]]:
    """Worker: generate + parse + write one batch of structures.

    Args:
        tuple of (struct_enum_out_abs, start, end_exclusive, index_species,
                  ordered_sites, out_dir, basename, fmt)

    Returns a list of (global_index, output_path), in index order.
    """
    (struct_enum_out, start, end_exclusive, index_species, ordered_sites,
     out_dir, basename, fmt) = args

    paths: list[tuple[int, str]] = []
    # Per-batch temp dir so vasp.* files from parallel batches never collide.
    with tempfile.TemporaryDirectory(prefix="xtalkit_batch_") as batch_dir:
        cmd = _makestr_range_args(struct_enum_out, start, end_exclusive)
        subprocess.run(cmd, cwd=batch_dir, check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.PIPE)
        vasp_files = sorted(
            glob.glob(os.path.join(batch_dir, "vasp.*")),
            key=_vasp_numeric_key,
        )
        for i, vf in enumerate(vasp_files):
            with open(vf, encoding="utf-8") as f:
                data = f.read()
            struct = _parse_vasp_to_structure(data, index_species, ordered_sites)
            idx = start + i
            path = os.path.join(out_dir, f"{basename}_{idx:03d}.{fmt}")
            if fmt == "xyz":
                from pymatgen.io.xyz import XYZ
                XYZ(struct).write_file(path)
            else:
                struct.to(filename=path, fmt="cif")
            paths.append((idx, path))
    return paths


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
    jobs: int = 1,
    batch_size: int = 256,
    scratch_dir: str | None = None,
    skip_preflight: bool = False,
) -> list[str]:
    """Enumerate ordered configurations of a disordered CIF.

    Returns list of output file paths. Raises RuntimeError if pymatgen or
    enumlib binaries are unavailable, if enumeration times out, or if
    enumlib returns 0 structures (e.g. non-integer occupancy with too small
    a max_cell_size).

    Args:
        jobs: number of parallel worker processes for the makestr+parse+write
            phase. 1 = serial streaming (default, lowest memory). 0 = auto
            (``os.cpu_count()``). Each worker loads its own copy of pymatgen,
            so memory scales with ``jobs``.
        batch_size: structures per batch. Larger = fewer ``makestr.x``
            invocations but higher peak memory.
        scratch_dir: directory for the enumlib scratch (``struct_enum.in/out``
            and per-batch ``vasp.*`` intermediates). Defaults to the system
            temp; ``/dev/shm`` puts them on tmpfs (faster, but limited to
            ~half of RAM).
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

    # Pre-flight: refuse non-integerizable stoichiometry BEFORE running enum.x.
    # enumlib cannot place a non-integer count of a species; without this guard
    # it tends to run away and exhaust memory (crashing the system) rather than
    # failing cleanly. Catch it here with actionable guidance instead.
    if not skip_preflight:
        bad = _non_integerizable_species(struct, symm_prec, min_cell_size, max_cell_size)
        if bad:
            lines = [
                "Cannot enumerate: non-integer stoichiometry "
                "(enumlib would likely run away and exhaust memory).",
                "",
                "  species  mult  occ      mult*occ   nearest valid (at "
                f"cell_size {min_cell_size})",
            ]
            for b in bad:
                lines.append(
                    f"  {b['species']:<7} {b['mult']:<5} {b['occ']:.4f}   "
                    f"{b['count']:.4f}    -> {b['suggested']:.4f}"
                )
            lines.append("")
            lines.append(
                "enumlib needs each species' count (mult * occ * cell_size) to "
                f"be an integer for some cell_size in [{min_cell_size}, {max_cell_size}]."
            )
            lines.append(
                "Fix: rebuild the CIF with rounded occupancies (e.g. "
                "`xtalkit build --atom-frac ...` using the 'nearest valid' "
                "values above), then re-run enumerate. See README 'When "
                "enumlib returns 0 structures'."
            )
            raise RuntimeError("\n".join(lines))

    basename = os.path.splitext(os.path.basename(resolved))[0]
    out_dir = os.path.abspath(output_dir or f"{basename}_enum")
    os.makedirs(out_dir, exist_ok=True)

    adaptor = EnumlibAdaptor(
        struct,
        min_cell_size=min_cell_size,
        max_cell_size=max_cell_size,
        symm_prec=symm_prec,
        timeout=timeout,
    )

    # Run enum.x ourselves (via the adaptor's helpers) so we keep struct_enum.out
    # and can drive makestr.x in streamed batches. enum.x writes struct_enum.in
    # and struct_enum.out to the CWD, so run inside a scratch dir.
    scratch = tempfile.TemporaryDirectory(prefix="xtalkit_enum_", dir=scratch_dir)
    try:
        prev_cwd = os.getcwd()
        os.chdir(scratch.name)
        try:
            adaptor._gen_input_file()          # writes struct_enum.in
            num_structs = adaptor._run_multienum()  # runs enum.x -> struct_enum.out
        finally:
            os.chdir(prev_cwd)

        if num_structs <= 0:
            raise EnumError("Unable to enumerate structure.")

        struct_enum_out = os.path.abspath(
            os.path.join(scratch.name, "struct_enum.out"))
        index_species = adaptor.index_species
        ordered_sites = adaptor.ordered_sites

        # Cap generation at max_structures (the limit now actually reduces work).
        effective_n = min(num_structs, max_structures) if max_structures else num_structs

        # Disjoint [start, end) batches.
        chunks = [(s, min(s + batch_size, effective_n))
                  for s in range(0, effective_n, batch_size)]
        chunk_args = [
            (struct_enum_out, s, e, index_species, ordered_sites,
             out_dir, basename, format)
            for (s, e) in chunks
        ]

        # jobs: 1 = serial; 0 = auto (cpu_count); >1 = that many workers.
        n_workers = (jobs if jobs > 1
                     else (os.cpu_count() or 1 if jobs == 0 else 1))
        if n_workers > 1 and len(chunks) > 1:
            with ProcessPoolExecutor(max_workers=min(n_workers, len(chunks))) as ex:
                results = list(ex.map(_process_chunk, chunk_args))
        else:
            results = [_process_chunk(a) for a in chunk_args]

        # ex.map preserves input order; each batch returns index-ordered paths,
        # so the flattened list is already in global-index order.
        paths = [p for batch in results for (_, p) in batch]
        return paths
    except EnumError as e:
        raise RuntimeError(
            f"enumlib returned 0 structures: {e}\n"
            f"Likely cause: non-integer stoichiometry at max_cell_size={max_cell_size}.\n"
            "Try a larger --max-cell-size, or provide a 'clean' parent CIF "
            "with integer stoichiometry (e.g., round 0.56 -> 0.5)."
        ) from e
    except TimeoutError as e:
        raise RuntimeError(
            f"enum.x exceeded the {timeout}-minute timeout.\n"
            "Increase --timeout, reduce --max-cell-size, or use a cleaner "
            "parent CIF. (enum.x is single-threaded; --jobs does not speed it up.)"
        ) from e
    finally:
        scratch.cleanup()
