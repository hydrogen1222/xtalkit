"""Verification checks for SHRY output directories."""

from __future__ import annotations

import json
import os

from xtalkit.enumeration.degeneracy import degeneracy_sum
from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.fingerprint import structure_fingerprint
from xtalkit.enumeration.formula import assert_formula
from xtalkit.enumeration.manifest import write_json
from xtalkit.enumeration.shry_backend import ShryBackend
from xtalkit.enumeration.shry_count import count_shry_structures
from xtalkit.enumeration.supercell_backend import run_supercell_count


def verify_shry_outputs(
    output_dir: str,
    check_count: bool = True,
    check_formula: bool = True,
    check_dedup: bool = True,
    target_formula: str | None = None,
    vacancy_symbol: str = "X",
    symprec_list: list[float] | None = None,
    cross_backend: str | None = None,
    backend: ShryBackend | None = None,
    check_degeneracy: bool = False,
) -> dict:
    """Verify count, formula, residual vacancies, and duplicate buckets."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    if not os.path.isdir(clean_dir):
        raise FileNotFoundError(f"clean_cif directory not found: {clean_dir}")
    files = _cif_files(clean_dir)

    result = {
        "output_dir": os.path.abspath(output_dir),
        "clean_count": len(files),
        "count_ok": None,
        "formula_errors": [],
        "vacancy_errors": [],
        "duplicate_buckets": {},
        "symprec_scan": [],
        "cross_backend": None,
        "degeneracy_sum": None,
    }
    manifest = _manifest(output_dir)
    if check_count:
        expected = manifest.get("count_only_result") if manifest else None
        result["expected_count"] = expected
        result["count_ok"] = expected is None or expected == len(files)
        if expected is not None and expected != len(files):
            raise ValueError(f"count mismatch: got {len(files)}, expected {expected}")

    if check_formula:
        for path in files:
            data = read_cif(path)
            if any(a.type_symbol == vacancy_symbol for a in data.atoms):
                result["vacancy_errors"].append(path)
            try:
                assert_formula(data.atoms, target_formula, vacancy_symbol)
            except ValueError as exc:
                result["formula_errors"].append({"path": path, "error": str(exc)})
        if result["vacancy_errors"] or result["formula_errors"]:
            raise ValueError("formula/vacancy verification failed")

    if check_dedup:
        buckets: dict[str, list[str]] = {}
        for path in files:
            buckets.setdefault(structure_fingerprint(path), []).append(path)
        result["duplicate_buckets"] = {
            fp: paths for fp, paths in buckets.items() if len(paths) > 1
        }
        if result["duplicate_buckets"]:
            raise ValueError(
                f"possible duplicate structures in "
                f"{len(result['duplicate_buckets'])} fingerprint bucket(s)"
            )

    if symprec_list:
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        angle = manifest.get("angle_tolerance", 5.0) if manifest else 5.0
        atol = manifest.get("atol", 1e-5) if manifest else 1e-5
        for symprec in symprec_list:
            payload = count_shry_structures(
                input_cif,
                scaling_matrix=matrix,
                symprec=float(symprec),
                angle_tolerance=float(angle),
                atol=float(atol),
                out_json=os.path.join(output_dir, "checks", f"count_symprec_{symprec}.json"),
                backend=backend,
            )
            result["symprec_scan"].append({
                "symprec": float(symprec),
                "count": payload["count_only_result"],
            })

    if cross_backend:
        if cross_backend != "supercell":
            raise ValueError(f"unknown cross backend {cross_backend!r}")
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        count = run_supercell_count(input_cif, scaling_matrix=matrix)
        result["cross_backend"] = {"backend": "supercell", "count": count}
        expected = result.get("expected_count")
        if expected is not None and count != expected:
            raise ValueError(
                f"supercell count mismatch: got {count}, expected {expected}"
            )

    if check_degeneracy:
        input_cif = _manifest_input_cif(output_dir, manifest)
        total = degeneracy_sum(files, input_cif, vacancy_symbol=vacancy_symbol)
        result["degeneracy_sum"] = total

    checks_dir = os.path.join(output_dir, "checks")
    os.makedirs(checks_dir, exist_ok=True)
    write_json(os.path.join(checks_dir, "verify.json"), result)
    return result


def _cif_files(path: str) -> list[str]:
    return sorted(
        os.path.join(path, name)
        for name in os.listdir(path)
        if name.lower().endswith(".cif")
    )


def _expected_count(output_dir: str) -> int | None:
    data = _manifest(output_dir)
    return data.get("count_only_result") if data else None


def _manifest(output_dir: str) -> dict:
    manifest = os.path.join(output_dir, "manifest.json")
    if not os.path.exists(manifest):
        return {}
    with open(manifest, encoding="utf-8") as f:
        return json.load(f)


def _manifest_input_cif(output_dir: str, manifest: dict) -> str:
    path = manifest.get("input_cif")
    if path and os.path.exists(path):
        return path
    fallback = os.path.join(output_dir, "input", "shry_ready.cif")
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError("Cannot locate SHRY input CIF from manifest")
