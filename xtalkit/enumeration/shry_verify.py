"""Verification checks for SHRY output directories (plan §4.4, §7, §10)."""

from __future__ import annotations

import json
import os

from xtalkit.enumeration.degeneracy import degeneracy_sum
from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.fingerprint import verify_no_duplicates
from xtalkit.enumeration.formula import assert_formula
from xtalkit.enumeration.manifest import append_jsonl, write_json
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
    """Verify count, formula, residual vacancies, and duplicates.

    Per plan §7 each check is written to its own file under ``checks/``:
    ``count.json``, ``formula_check.jsonl`` (one row per structure),
    ``dedup_check.json``, ``symprec_scan.json``, ``orbit_grouping.json``;
    ``verify.json`` is the overall summary.
    """
    clean_dir = os.path.join(output_dir, "clean_cif")
    if not os.path.isdir(clean_dir):
        raise FileNotFoundError(f"clean_cif directory not found: {clean_dir}")
    files = _cif_files(clean_dir)
    checks_dir = os.path.join(output_dir, "checks")
    os.makedirs(checks_dir, exist_ok=True)

    manifest = _manifest(output_dir)
    result = {
        "output_dir": os.path.abspath(output_dir),
        "clean_count": len(files),
        "count_ok": None,
        "formula_errors": [],
        "vacancy_errors": [],
        "dedup": None,
        "symprec_scan": [],
        "cross_backend": None,
        "degeneracy_sum": None,
    }

    if check_count:
        expected = manifest.get("count_only_result") if manifest else None
        count_result = {
            "expected": expected,
            "actual": len(files),
            "ok": expected is None or expected == len(files),
        }
        result["expected_count"] = expected
        result["count_ok"] = count_result["ok"]
        write_json(os.path.join(checks_dir, "count.json"), count_result)
        if expected is not None and expected != len(files):
            raise ValueError(f"count mismatch: got {len(files)}, expected {expected}")

    if check_formula:
        rows = []
        for path in files:
            data = read_cif(path)
            has_vacancy = any(a.type_symbol == vacancy_symbol for a in data.atoms)
            if has_vacancy:
                result["vacancy_errors"].append(path)
            try:
                formula = assert_formula(data.atoms, target_formula, vacancy_symbol)
                rows.append({"path": path, "formula": formula, "ok": True})
            except ValueError as exc:
                rows.append({"path": path, "formula": None, "ok": False,
                             "error": str(exc)})
                result["formula_errors"].append({"path": path, "error": str(exc)})
        append_jsonl(os.path.join(checks_dir, "formula_check.jsonl"), rows)
        if result["vacancy_errors"] or result["formula_errors"]:
            raise ValueError("formula/vacancy verification failed")

    if check_dedup:
        dedup = verify_no_duplicates(files)
        result["dedup"] = dedup
        write_json(os.path.join(checks_dir, "dedup_check.json"), dedup)
        if dedup["has_duplicates"]:
            raise ValueError(
                f"StructureMatcher-confirmed duplicates in "
                f"{len(dedup['duplicates'])} bucket(s); see dedup_check.json"
            )

    if symprec_list:
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        angle = manifest.get("angle_tolerance", 5.0) if manifest else 5.0
        atol = manifest.get("atol", 1e-5) if manifest else 1e-5
        scan = []
        for symprec in symprec_list:
            payload = count_shry_structures(
                input_cif,
                scaling_matrix=matrix,
                symprec=float(symprec),
                angle_tolerance=float(angle),
                atol=float(atol),
                out_json=os.path.join(checks_dir, f"count_symprec_{symprec}.json"),
                backend=backend,
            )
            scan.append({"symprec": float(symprec),
                         "count": payload["count_only_result"]})
        result["symprec_scan"] = scan
        write_json(os.path.join(checks_dir, "symprec_scan.json"), {"results": scan})

    if cross_backend:
        if cross_backend != "supercell":
            raise ValueError(f"unknown cross backend {cross_backend!r}")
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        count = run_supercell_count(input_cif, scaling_matrix=matrix)
        cross = {"backend": "supercell", "count": count}
        expected = result.get("expected_count")
        if expected is not None and count != expected:
            raise ValueError(
                f"supercell count mismatch: got {count}, expected {expected}")
        result["cross_backend"] = cross
        write_json(os.path.join(checks_dir, "cross_backend.json"), cross)

    if check_degeneracy:
        input_cif = _manifest_input_cif(output_dir, manifest)
        total = degeneracy_sum(files, input_cif, vacancy_symbol=vacancy_symbol)
        result["degeneracy_sum"] = total
        write_json(os.path.join(checks_dir, "degeneracy.json"),
                   {"degeneracy_sum": total})

    # Orbit grouping record (plan §7 checks/orbit_grouping.json): prefer the
    # manifest's site_orbits, else the prepare sidecar if present.
    orbits = manifest.get("site_orbits") if manifest else None
    if orbits is None:
        sidecar = os.path.join(output_dir, "input", "shry_ready.cif.orbit_grouping.json")
        if os.path.exists(sidecar):
            with open(sidecar, encoding="utf-8") as f:
                orbits = json.load(f).get("site_orbits")
    if orbits is not None:
        write_json(os.path.join(checks_dir, "orbit_grouping.json"),
                   {"site_orbits": orbits})

    write_json(os.path.join(checks_dir, "verify.json"), result)
    return result


def _cif_files(path: str) -> list[str]:
    return sorted(
        os.path.join(path, name)
        for name in os.listdir(path)
        if name.lower().endswith(".cif")
    )


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
