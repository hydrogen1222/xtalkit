"""Strict SHRY generation and streaming post-processing."""

from __future__ import annotations

import difflib
import os

from xtalkit.enumeration.cif_io import CifData, copy_file, read_cif, write_cif
from xtalkit.enumeration.degeneracy import compute_degeneracy
from xtalkit.enumeration.formula import assert_formula
from xtalkit.enumeration.manifest import (
    append_jsonl,
    base_manifest,
    pip_freeze,
    sha256_file,
    write_json,
)
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
    pip_freeze(os.path.join(paths["input"], "pip_freeze.txt"))

    mod_audit = _run_mod_only_audit(
        backend, cif_path, paths["input"], matrix, symprec, angle_tolerance,
        atol, symmetrize=symmetrize, strict=strict)

    args = [
        cif_path,
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
        "--output-dir",
        paths["raw"],
        "--disable-progressbar",
    ]
    if write_cif_output:
        args.append("--write-cif")
    if write_poscar:
        args.append("--write-poscar")
    if symmetrize:
        args.append("--symmetrize")

    result = backend.run(args)
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
    )
    write_json(os.path.join(output_dir, "manifest.json"), manifest)
    return manifest


def iter_shry_outputs(raw_dir: str):
    """Yield SHRY-generated CIF files in stable order."""
    for root, dirs, files in os.walk(raw_dir):
        dirs.sort()
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
    args = [
        cif_path,
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
    result = backend.run(args)
    mod_path = os.path.join(input_dir, "shry_mod_only.cif")
    if "data_" in result.stdout:
        with open(mod_path, "w", encoding="utf-8") as f:
            f.write(result.stdout)
        diff = _text_diff(cif_path, mod_path)
        if strict and diff:
            raise ValueError(
                "SHRY --mod-only changed the prepared CIF; inspect "
                f"{mod_path} before enumerating."
            )
        return {"command": " ".join(result.command), "mod_only_cif": mod_path,
                "diff": diff}
    return {"command": " ".join(result.command), "mod_only_cif": None,
            "diff": []}


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
        clean = CifData(
            raw.block_name,
            raw.cell,
            raw.spacegroup_name,
            raw.spacegroup_number,
            raw.symops,
            [a for a in raw.atoms if a.type_symbol != remove_vacancy],
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


def _text_diff(a: str, b: str) -> list[str]:
    with open(a, encoding="utf-8") as fa, open(b, encoding="utf-8") as fb:
        left = [line.rstrip() for line in fa]
        right = [line.rstrip() for line in fb]
    return list(difflib.unified_diff(left, right, fromfile=a, tofile=b, lineterm=""))
