"""Strict SHRY generation and streaming post-processing."""

from __future__ import annotations

import glob
import os
import tempfile
from dataclasses import replace
from fractions import Fraction
from functools import reduce
from math import gcd, lcm

from xtalkit.enumeration.cif_io import CifData, copy_file, read_cif, write_cif
from xtalkit.enumeration.degeneracy import compute_degeneracy
from xtalkit.enumeration.formula import assert_formula, element_symbol
from xtalkit.enumeration.manifest import (
    append_jsonl,
    base_manifest,
    pip_freeze,
    sha256_file,
    write_json,
)
from xtalkit.enumeration.occupancy import parse_fraction
from xtalkit.enumeration.shry_backend import ShryBackend, matrix_args


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
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)
    if dir_size < 1:
        raise ValueError("dir_size must be >= 1")
    if write_degeneracy and not write_cif_output:
        raise ValueError("write_degeneracy requires write_cif_output=True")

    matrix = scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    backend = backend or ShryBackend()
    paths = _make_layout(output_dir)
    copy_file(cif_path, os.path.join(paths["input"], "shry_ready.cif"))
    # Carry the prepare-stage orbit grouping sidecar (if any) into the workflow
    # so verify can write checks/orbit_grouping.json (plan §7).
    orbit_sidecar = f"{cif_path}.orbit_grouping.json"
    site_orbits = None
    hall_symbol = None
    parent_sg_detected = None
    if os.path.exists(orbit_sidecar):
        copy_file(orbit_sidecar, os.path.join(paths["input"], "shry_ready.cif.orbit_grouping.json"))
        import json as _json
        with open(orbit_sidecar, encoding="utf-8") as _f:
            _orbit_data = _json.load(_f)
        site_orbits = _orbit_data.get("site_orbits")
        hall_symbol = _orbit_data.get("hall_symbol")
        parent_sg_detected = _orbit_data.get("parent_spacegroup_detected")
    _write_pip_freeze(backend, os.path.join(paths["input"], "pip_freeze.txt"))

    mod_audit = _run_mod_only_audit(
        backend, cif_path, paths["input"], matrix, symprec, angle_tolerance,
        atol, symmetrize=symmetrize, strict=strict)

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
    result = backend.run(args, cwd=paths["raw"])
    with open(os.path.join(paths["input"], "command.txt"), "w", encoding="utf-8") as f:
        f.write(" ".join(result.command) + "\n")

    generated = _postprocess_raw_outputs(
        paths=paths,
        target_formula=target_formula,
        remove_vacancy=remove_vacancy,
        write_cif_output=write_cif_output,
        write_poscar_output=write_poscar,
        parent_cif=cif_path,
        write_degeneracy=write_degeneracy,
    )
    if expect_count is not None and generated != expect_count:
        raise ValueError(
            f"SHRY generated {generated} structures, expected {expect_count}"
        )

    manifest = base_manifest(
        mode="strict_exhaustive_success" if strict else "shry_enum",
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
        generated_raw_count=generated,
        generated_clean_count=generated,
        degeneracy_written=bool(write_degeneracy),
        shry_command_line=" ".join(result.command),
        mod_only_audit=mod_audit,
        backend_version=backend.version(),
        site_orbits=site_orbits,
        hall_symbol=hall_symbol,
        parent_spacegroup_detected=parent_sg_detected,
    )
    write_json(os.path.join(output_dir, "manifest.json"), manifest)
    return manifest


def _write_pip_freeze(backend, path: str) -> None:
    """Write the SHRY isolated env's pip freeze (plan §5/§8).

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
        dirs.sort()
        if not os.path.basename(root).startswith("slice"):
            continue
        for name in sorted(files):
            if name.lower().endswith(".cif"):
                yield os.path.join(root, name)


def _make_layout(output_dir: str) -> dict[str, str]:
    paths = {
        "root": output_dir,
        "input": os.path.join(output_dir, "input"),
        "raw": os.path.join(output_dir, "raw_shry"),
        "clean": os.path.join(output_dir, "clean_cif"),
        "poscar": os.path.join(output_dir, "poscar"),
        "checks": os.path.join(output_dir, "checks"),
    }
    for path in paths.values():
        os.makedirs(path, exist_ok=True)
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
    return {"command": " ".join(result.command), "mod_only_cif": mod_path,
            "diff": diff}


def _reduced_composition(rows) -> dict[str, int]:
    """Reduce atom rows to a GCD-normalized {element: count} for comparison.

    Robust to cell choice (conventional vs primitive) and to whether a partial
    site is split into two rows (Li 0.5 + X 0.5) or kept as one: it sums
    occupancies per species, clears denominators, and divides by the GCD.
    Species labels are normalized to bare element symbols so that ``Li1+``
    (refinement-style oxidation notation in the input) and ``Li+`` (pymatgen's
    normalized form in SHRY's modified CIF) compare equal — the audit checks
    chemistry, not notation.
    """
    comp: dict[str, Fraction] = {}
    for r in rows:
        occ = parse_fraction(r.occupancy)
        el = element_symbol(r.type_symbol)
        comp[el] = comp.get(el, Fraction(0)) + occ
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
    a = _reduced_composition(read_cif(input_cif).atoms)
    b = _reduced_composition(read_cif(modified_cif).atoms)
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
) -> int:
    manifest_jsonl = os.path.join(paths["root"], "manifest.jsonl")
    count = 0
    for raw_file in iter_shry_outputs(paths["raw"]):
        count += 1
        raw = read_cif(raw_file)
        # SHRY emits type_symbols carrying oxidation states (Li+, S2-, ...).
        # Normalize to bare element symbols and drop the vacancy species so
        # clean CIFs/POSCARs and formula checks use plain chemistry.
        vac = element_symbol(remove_vacancy)
        clean_atoms = [
            replace(a, type_symbol=element_symbol(a.type_symbol))
            for a in raw.atoms
            if element_symbol(a.type_symbol) != vac
        ]
        clean = CifData(
            raw.block_name,
            raw.cell,
            raw.spacegroup_name,
            raw.spacegroup_number,
            raw.symops,
            clean_atoms,
        )
        formula = assert_formula(clean.atoms, target_formula, remove_vacancy)
        clean_path = os.path.join(paths["clean"], f"conf_{count:06d}.cif")
        poscar_path = os.path.join(paths["poscar"], f"POSCAR_{count:06d}")
        if write_cif_output:
            write_cif(clean, clean_path)
        else:
            clean_path = None
        if write_poscar_output:
            _write_poscar(clean, poscar_path)
        else:
            poscar_path = None
        degeneracy = compute_degeneracy(clean_path, parent_cif, remove_vacancy) \
            if write_degeneracy and clean_path else None
        append_jsonl(manifest_jsonl, {
            "index": count,
            "raw_file": raw_file,
            "clean_cif": clean_path,
            "poscar": poscar_path,
            "formula": formula,
            "removed_vacancy_count": len(raw.atoms) - len(clean.atoms),
            "degeneracy": degeneracy,
            "status": "ok",
        })
    return count


def _write_poscar(data: CifData, path: str) -> None:
    cell = data.cell
    species_order: list[str] = []
    grouped: dict[str, list] = {}
    for atom in data.atoms:
        if atom.type_symbol not in grouped:
            species_order.append(atom.type_symbol)
            grouped[atom.type_symbol] = []
        grouped[atom.type_symbol].append(atom)
    with open(path, "w", encoding="utf-8") as f:
        f.write("xtalkit shry\n1.0\n")
        f.write(f"{float(cell['a']):.10f} 0 0\n")
        f.write(f"0 {float(cell['b']):.10f} 0\n")
        f.write(f"0 0 {float(cell['c']):.10f}\n")
        f.write(" ".join(species_order) + "\n")
        f.write(" ".join(str(len(grouped[sp])) for sp in species_order) + "\n")
        f.write("Direct\n")
        for sp in species_order:
            for atom in grouped[sp]:
                f.write(f"{float(atom.x):.10f} {float(atom.y):.10f} {float(atom.z):.10f}\n")
