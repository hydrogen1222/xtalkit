"""Strict SHRY generation and streaming post-processing."""

from __future__ import annotations

import glob
import os
import shlex
import tempfile
from dataclasses import replace
from fractions import Fraction
from functools import reduce
from math import gcd, lcm

from xtalkit.enumeration.cif_io import CifData, copy_file, read_cif, write_cif
from xtalkit.enumeration.degeneracy import compute_degeneracy
from xtalkit.enumeration.formula import (
    assert_formula,
    counts_from_atom_rows,
    element_symbol,
    real_element_symbol,
    validate_vacancy_symbol,
)
from xtalkit.enumeration.manifest import (
    append_jsonl,
    base_manifest,
    pip_freeze,
    sha256_file,
    write_json,
)
from xtalkit.enumeration.shry_backend import (
    ShryBackend,
    matrix_args,
    parse_total_combinations,
    scaling_determinant,
    validate_scaling_matrix,
    validate_tolerances,
)


def enumerate_with_shry(
    cif_path: str,
    output_dir: str,
    scaling_matrix: list[list[int]] | None = None,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
    atol: float = 1e-5,
    expect_count: int | None = None,
    remove_vacancy: str = "X",
    target_formula: str | None = None,
    write_cif_output: bool = True,
    write_poscar: bool = False,
    write_degeneracy: bool = False,
    dir_size: int = 10000,
    strict: bool = True,
    symmetrize: bool = False,
    backend: ShryBackend | None = None,
) -> dict:
    """Run exhaustive SHRY enumeration and stream clean outputs."""
    if strict and expect_count is None:
        raise ValueError("--expect-count is required in strict SHRY enum mode")
    if expect_count is not None and expect_count <= 0:
        raise ValueError("expect_count must be a positive integer")
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)
    if dir_size < 1:
        raise ValueError("dir_size must be >= 1")
    validate_tolerances(symprec, angle_tolerance, atol)
    remove_vacancy = validate_vacancy_symbol(remove_vacancy)
    if write_degeneracy and not write_cif_output:
        raise ValueError("write_degeneracy requires write_cif_output=True")

    matrix = validate_scaling_matrix(
        scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    )
    if write_degeneracy and scaling_determinant(matrix) != 1:
        raise ValueError(
            "write_degeneracy is only rigorous for a unit-determinant scaling "
            "matrix; disable it for supercell enumeration"
        )
    backend = backend or ShryBackend()
    paths = _make_layout(output_dir)
    try:
        copy_file(cif_path, os.path.join(paths["input"], "shry_ready.cif"))
    except Exception as exc:
        _write_failure_manifest(
            paths, cif_path, matrix, symprec, angle_tolerance, atol,
            strict, symmetrize, remove_vacancy, target_formula, expect_count,
            backend, "input_copy", exc,
        )
        raise
    # Carry the prepare-stage orbit grouping sidecar (if any) into the workflow
    # so verify can write checks/orbit_grouping.json.
    orbit_sidecar = f"{cif_path}.orbit_grouping.json"
    site_orbits = None
    hall_symbol = None
    parent_sg_detected = None
    try:
        if os.path.exists(orbit_sidecar):
            copy_file(orbit_sidecar, os.path.join(
                paths["input"], "shry_ready.cif.orbit_grouping.json"))
            import json as _json
            with open(orbit_sidecar, encoding="utf-8") as _f:
                _orbit_data = _json.load(_f)
            site_orbits = _orbit_data.get("site_orbits")
            hall_symbol = _orbit_data.get("hall_symbol")
            parent_sg_detected = _orbit_data.get("parent_spacegroup_detected")
            if write_degeneracy and parent_sg_detected is not None:
                declared = read_cif(cif_path).spacegroup_number
                if declared != int(parent_sg_detected):
                    raise ValueError(
                        "write_degeneracy requires the prepared CIF's declared "
                        "space group to equal the spglib-detected enumeration "
                        f"group; declared {declared}, detected {parent_sg_detected}"
                    )
        _write_pip_freeze(
            backend, os.path.join(paths["input"], "pip_freeze.txt"))
    except Exception as exc:
        _write_failure_manifest(
            paths, cif_path, matrix, symprec, angle_tolerance, atol,
            strict, symmetrize, remove_vacancy, target_formula, expect_count,
            backend, "input_audit_metadata", exc,
        )
        raise

    try:
        mod_audit = _run_mod_only_audit(
            backend, cif_path, paths["input"], matrix, symprec, angle_tolerance,
            atol, symmetrize=symmetrize, strict=strict)
    except Exception as exc:
        _write_failure_manifest(
            paths, cif_path, matrix, symprec, angle_tolerance, atol,
            strict, symmetrize, remove_vacancy, target_formula, expect_count,
            backend, "mod_only_audit", exc,
        )
        raise

    args = [
        os.path.abspath(cif_path),
        "--scaling-matrix",
        *matrix_args(matrix),
        "--symprec",
        str(symprec),
        "--angle-tolerance",
        str(angle_tolerance),
        "--atol",
        str(atol),
        "--dir-size",
        str(dir_size),
        "--disable-progressbar",
    ]
    if symmetrize:
        args.append("--symmetrize")
    # SHRY has no --output-dir/--write-cif/--write-poscar flags. It always
    # writes CIFs (to shry-<basename>-<timestamp>/slice*/ in the CWD) and has
    # no POSCAR writer at all. We run it with cwd=raw so its output lands in
    # paths["raw"], then build clean CIFs/POSCARs ourselves in post-processing.
    try:
        result = backend.run(args, cwd=paths["raw"])
    except Exception as exc:
        _write_failure_manifest(
            paths, cif_path, matrix, symprec, angle_tolerance, atol,
            strict, symmetrize, remove_vacancy, target_formula, expect_count,
            backend, "shry_enumeration", exc, mod_audit=mod_audit,
        )
        raise
    with open(os.path.join(paths["input"], "command.txt"), "w", encoding="utf-8") as f:
        f.write(shlex.join(result.command) + "\n")

    try:
        generated = _postprocess_raw_outputs(
            paths=paths,
            target_formula=target_formula,
            remove_vacancy=remove_vacancy,
            write_cif_output=write_cif_output,
            write_poscar_output=write_poscar,
            parent_cif=cif_path,
            write_degeneracy=write_degeneracy,
        )
    except Exception as exc:
        _write_failure_manifest(
            paths, cif_path, matrix, symprec, angle_tolerance, atol,
            strict, symmetrize, remove_vacancy, target_formula, expect_count,
            backend, "postprocess", exc, command=result.command,
            mod_audit=mod_audit,
        )
        raise
    count_matches = (
        expect_count is None or generated["generated_raw_count"] == expect_count
    )

    manifest = base_manifest(
        mode=("strict_exhaustive_success" if strict else "shry_enum")
        if count_matches else "shry_enum_count_mismatch",
        status="success" if count_matches else "failed",
        failure_reason=None if count_matches else "generated_count_mismatch",
        input_cif=os.path.abspath(cif_path),
        input_sha256=sha256_file(cif_path),
        supercell_matrix=matrix,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        atol=atol,
        symmetrize_flag=bool(symmetrize),
        vacancy_symbol=remove_vacancy,
        target_formula_after_removing_vacancy=target_formula,
        count_only_result=expect_count,
        total_combinations=parse_total_combinations(result.stdout, result.stderr),
        generated_raw_count=generated["generated_raw_count"],
        generated_clean_count=generated["generated_clean_count"],
        generated_poscar_count=generated["generated_poscar_count"],
        write_cif_output=bool(write_cif_output),
        write_poscar_output=bool(write_poscar),
        degeneracy_written=bool(write_degeneracy),
        shry_command_line=shlex.join(result.command),
        mod_only_audit=mod_audit,
        backend_version=backend.version(),
        site_orbits=site_orbits,
        hall_symbol=hall_symbol,
        parent_spacegroup_detected=parent_sg_detected,
    )
    write_json(os.path.join(paths["root"], "manifest.json"), manifest)
    if not count_matches:
        raise ValueError(
            "SHRY generated "
            f"{generated['generated_raw_count']} structures, expected {expect_count}; "
            f"failure manifest written to {os.path.join(paths['root'], 'manifest.json')}"
        )
    return manifest


def _write_failure_manifest(
    paths: dict[str, str],
    cif_path: str,
    matrix: list[list[int]],
    symprec: float,
    angle_tolerance: float,
    atol: float,
    strict: bool,
    symmetrize: bool,
    vacancy_symbol: str,
    target_formula: str | None,
    expect_count: int | None,
    backend,
    stage: str,
    error: Exception,
    command: list[str] | None = None,
    mod_audit: dict | None = None,
) -> None:
    """Best-effort audit artifact for an interrupted or failed enum stage.

    This helper deliberately swallows its own diagnostic-write failures so the
    scientifically relevant original exception is never replaced.
    """
    manifest_path = os.path.join(paths["root"], "manifest.json")
    if os.path.exists(manifest_path):
        return
    try:
        try:
            backend_version = backend.version()
        except Exception:  # noqa: BLE001 - failure reporting must be best effort
            backend_version = None
        raw_count = sum(1 for _ in iter_shry_outputs(paths["raw"]))
        clean_count = len(glob.glob(os.path.join(paths["clean"], "*.cif")))
        poscar_count = len(glob.glob(os.path.join(paths["poscar"], "POSCAR*")))
        payload = base_manifest(
            mode="strict_exhaustive_failed" if strict else "shry_enum_failed",
            status="failed",
            failure_stage=stage,
            failure_type=type(error).__name__,
            failure_reason=str(error),
            input_cif=os.path.abspath(cif_path),
            input_sha256=sha256_file(cif_path),
            supercell_matrix=matrix,
            symprec=symprec,
            angle_tolerance=angle_tolerance,
            atol=atol,
            symmetrize_flag=bool(symmetrize),
            vacancy_symbol=vacancy_symbol,
            target_formula_after_removing_vacancy=target_formula,
            count_only_result=expect_count,
            generated_raw_count=raw_count,
            generated_clean_count=clean_count,
            generated_poscar_count=poscar_count,
            shry_command_line=shlex.join(command) if command else None,
            mod_only_audit=mod_audit,
            backend_version=backend_version,
        )
        write_json(manifest_path, payload)
    except Exception:  # noqa: BLE001 - never mask the original workflow error
        return


def _write_pip_freeze(backend, path: str) -> None:
    """Write the SHRY isolated environment's package snapshot.

    Prefers the backend's own ``pip_freeze`` (which recovers the SHRY
    interpreter from the ``shry`` shebang); falls back to the manifest helper
    (xtalkit's env) for backends that don't expose it (e.g. test fakes).
    """
    backend_freeze = getattr(backend, "pip_freeze", None)
    if callable(backend_freeze):
        backend_freeze(path)
        return
    pip_freeze(path)


def iter_shry_outputs(raw_dir: str):
    """Yield SHRY-generated ordered-configuration CIFs in stable order.

    SHRY writes ``shry-<basename>-<timestamp>/sliceN/<formula>_<i>_<w>.cif``
    plus a top-level ``<basename>-<scaling>.cif`` that is the *disordered*
    modified parent (still partial occupancy). Only the ``slice*/`` CIFs are
    ordered configurations — the parent CIF must be skipped or post-processing
    crashes on its residual partial occupancy.
    """
    for root, dirs, files in os.walk(raw_dir):
        dirs.sort(key=_natural_sort_key)
        if not os.path.basename(root).startswith("slice"):
            continue
        for name in sorted(files, key=_natural_sort_key):
            if name.lower().endswith(".cif"):
                yield os.path.join(root, name)


def _natural_sort_key(value: str) -> tuple:
    """Sort SHRY ``slice2``/``slice10`` and numbered CIFs numerically."""
    return tuple(
        int(part) if part.isdigit() else part.casefold()
        for part in __import__("re").split(r"(\d+)", value)
    )


def _make_layout(output_dir: str) -> dict[str, str]:
    root = os.path.abspath(output_dir)
    if os.path.exists(root):
        if not os.path.isdir(root):
            raise NotADirectoryError(root)
        if os.listdir(root):
            raise FileExistsError(
                f"output directory is not empty: {root}. Choose a new or empty "
                "directory so files from separate SHRY runs cannot be mixed."
            )
    else:
        os.makedirs(root)
    paths = {
        "root": root,
        "input": os.path.join(root, "input"),
        "raw": os.path.join(root, "raw_shry"),
        "clean": os.path.join(root, "clean_cif"),
        "poscar": os.path.join(root, "poscar"),
        "checks": os.path.join(root, "checks"),
    }
    for name, path in paths.items():
        if name != "root":
            os.makedirs(path)
    return paths


def _run_mod_only_audit(
    backend: ShryBackend,
    cif_path: str,
    input_dir: str,
    matrix: list[list[int]],
    symprec: float,
    angle_tolerance: float,
    atol: float,
    symmetrize: bool,
    strict: bool,
) -> dict:
    """Run ``shry --mod-only`` and verify SHRY did not alter the input chemistry.

    SHRY writes the modified structure to ``shry-<basename>-<ts>/<basename>-<scaling>.cif``
    in its CWD (never to stdout), so we run it in a throwaway directory, locate
    that file, and compare its reduced composition to the prepared input. A
    text/coordinate diff is too fragile — SHRY re-symmetrizes (e.g. P1 -> Pm-3m)
    and re-emits via pymatgen, so even an unchanged structure rarely matches
    byte-for-byte. Composition catches the real danger: SHRY snapping
    occupancies or dropping/adding a species.
    """
    args = [
        os.path.abspath(cif_path),
        "--mod-only",
        "--scaling-matrix",
        *matrix_args(matrix),
        "--symprec",
        str(symprec),
        "--angle-tolerance",
        str(angle_tolerance),
        "--atol",
        str(atol),
        "--disable-progressbar",
    ]
    if symmetrize:
        args.append("--symmetrize")
    with tempfile.TemporaryDirectory(prefix="xtalkit_modonly_") as tmp:
        result = backend.run(args, cwd=tmp)
        candidates = sorted(glob.glob(os.path.join(tmp, "shry-*", "*.cif")))
        if not candidates:
            raise RuntimeError(
                "SHRY --mod-only did not write a modified CIF.\n"
                f"stdout:\n{result.stdout}\nstderr:\n{result.stderr}"
            )
        mod_path = os.path.join(input_dir, "shry_mod_only.cif")
        copy_file(candidates[0], mod_path)
    diff = _composition_diff(cif_path, mod_path)
    if strict and diff:
        raise ValueError(
            "SHRY --mod-only altered the prepared CIF chemistry "
            f"({diff}). Inspect {mod_path} before enumerating."
        )
    return {"command": shlex.join(result.command), "mod_only_cif": mod_path,
            "diff": diff}


def _reduced_composition(cif_path: str) -> dict[str, int]:
    """Reduce the expanded physical cell to a normalized element ratio.

    Parsing with pymatgen expands CIF symmetry before composition is counted.
    This is essential when species occupy Wyckoff positions with different
    multiplicities: summing asymmetric-unit rows alone is not a composition.
    Clearing rationalized occupancy denominators and dividing by the GCD keeps
    the comparison invariant to conventional-versus-primitive cell choice.
    """
    try:
        from pymatgen.core import Structure
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "pymatgen is required for the SHRY mod-only chemistry audit"
        ) from exc

    structure = Structure.from_file(cif_path)
    comp: dict[str, Fraction] = {}
    for site in structure:
        for species, occupancy in site.species.items():
            symbol = element_symbol(getattr(species, "symbol", str(species)))
            rational = Fraction(str(float(occupancy))).limit_denominator(1_000_000)
            comp[symbol] = comp.get(symbol, Fraction(0)) + rational
    if not comp:
        return {}
    denom_lcm = reduce(lcm, (v.denominator for v in comp.values()))
    ints = {k: int(v * denom_lcm) for k, v in comp.items() if v != 0}
    if not ints:
        return {}
    g = reduce(gcd, ints.values())
    return {k: v // g for k, v in ints.items()}


def _composition_diff(input_cif: str, modified_cif: str) -> str:
    """Return '' if chemistries match, else a one-line description."""
    a = _reduced_composition(input_cif)
    b = _reduced_composition(modified_cif)
    if a == b:
        return ""
    return f"input {a} != modified {b}"


def _postprocess_raw_outputs(
    paths: dict[str, str],
    target_formula: str | None,
    remove_vacancy: str,
    write_cif_output: bool,
    write_poscar_output: bool,
    parent_cif: str,
    write_degeneracy: bool,
) -> dict[str, int]:
    manifest_jsonl = os.path.join(paths["root"], "manifest.jsonl")
    raw_count = 0
    clean_count = 0
    poscar_count = 0
    for raw_file in iter_shry_outputs(paths["raw"]):
        raw_count += 1
        raw = read_cif(raw_file)
        # SHRY emits type_symbols carrying oxidation states (Li+, S2-, ...).
        # Normalize to bare element symbols and drop the vacancy species so
        # clean CIFs/POSCARs and formula checks use plain chemistry.
        vac = element_symbol(remove_vacancy)
        clean_atoms = []
        for atom in raw.atoms:
            if element_symbol(atom.type_symbol) == vac:
                continue
            clean_atoms.append(replace(
                atom, type_symbol=real_element_symbol(atom.type_symbol),
            ))
        clean = CifData(
            raw.block_name,
            raw.cell,
            raw.spacegroup_name,
            raw.spacegroup_number,
            raw.symops,
            clean_atoms,
        )
        formula = assert_formula(clean.atoms, target_formula, remove_vacancy)
        clean_path = os.path.join(paths["clean"], f"conf_{raw_count:06d}.cif")
        poscar_path = os.path.join(paths["poscar"], f"POSCAR_{raw_count:06d}")
        if write_cif_output:
            write_cif(clean, clean_path)
            clean_count += 1
        else:
            clean_path = None
        if write_poscar_output:
            _write_poscar(clean, poscar_path)
            poscar_count += 1
        else:
            poscar_path = None
        degeneracy = compute_degeneracy(clean_path, parent_cif, remove_vacancy) \
            if write_degeneracy and clean_path else None
        append_jsonl(manifest_jsonl, {
            "index": raw_count,
            "raw_file": raw_file,
            "clean_cif": clean_path,
            "poscar": poscar_path,
            "formula": formula,
            "removed_vacancy_count": len(raw.atoms) - len(clean.atoms),
            "degeneracy": degeneracy,
            "status": "ok",
        })
    return {
        "generated_raw_count": raw_count,
        "generated_clean_count": clean_count,
        "generated_poscar_count": poscar_count,
    }


def _write_poscar(data: CifData, path: str) -> None:
    """Write a POSCAR using the complete triclinic lattice geometry."""
    import gemmi

    counts_from_atom_rows(data.atoms)
    cell = data.cell
    species_order: list[str] = []
    grouped: dict[str, list] = {}
    for atom in data.atoms:
        species = real_element_symbol(atom.type_symbol)
        if species not in grouped:
            species_order.append(species)
            grouped[species] = []
        grouped[species].append(atom)
    with open(path, "w", encoding="utf-8") as f:
        f.write("xtalkit shry\n1.0\n")
        unit_cell = gemmi.UnitCell(
            float(cell["a"]), float(cell["b"]), float(cell["c"]),
            float(cell["alpha"]), float(cell["beta"]), float(cell["gamma"]),
        )
        for frac in (
            gemmi.Fractional(1, 0, 0),
            gemmi.Fractional(0, 1, 0),
            gemmi.Fractional(0, 0, 1),
        ):
            vector = unit_cell.orthogonalize(frac)
            f.write(f"{vector.x:.10f} {vector.y:.10f} {vector.z:.10f}\n")
        f.write(" ".join(species_order) + "\n")
        f.write(" ".join(str(len(grouped[sp])) for sp in species_order) + "\n")
        f.write("Direct\n")
        for sp in species_order:
            for atom in grouped[sp]:
                f.write(f"{float(atom.x):.10f} {float(atom.y):.10f} {float(atom.z):.10f}\n")
