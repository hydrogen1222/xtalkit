"""Verification checks for SHRY output directories."""

from __future__ import annotations

import json
import os

from xtalkit.enumeration.degeneracy import degeneracy_sum, raw_combination_count
from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.fingerprint import verify_no_duplicates
from xtalkit.enumeration.formula import (
    assert_formula,
    element_symbol,
    validate_vacancy_symbol,
)
from xtalkit.enumeration.manifest import write_json, write_jsonl
from xtalkit.enumeration.shry_backend import (
    ShryBackend,
    scaling_determinant,
    validate_scaling_matrix,
)
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

    Each check is written to its own file under ``checks/``:
    ``count.json``, ``formula_check.jsonl`` (one row per structure),
    ``dedup_check.json``, ``symprec_scan.json``, ``orbit_grouping.json``;
    ``verify.json`` is the overall summary.
    """
    clean_dir = os.path.join(output_dir, "clean_cif")
    vacancy_symbol = validate_vacancy_symbol(vacancy_symbol)
    if not os.path.isdir(clean_dir):
        raise FileNotFoundError(f"clean_cif directory not found: {clean_dir}")
    files = _cif_files(clean_dir)
    checks_dir = os.path.join(output_dir, "checks")
    os.makedirs(checks_dir, exist_ok=True)

    manifest = _manifest(output_dir)
    failures: list[dict[str, str]] = []
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
        "failures": failures,
        "ok": None,
    }
    if not files:
        failures.append({
            "check": "outputs",
            "error": "no clean CIF configurations found",
        })

    if check_count:
        expected = manifest.get("count_only_result") if manifest else None
        count_ok = expected is not None and expected == len(files)
        count_result = {
            "expected": expected,
            "actual": len(files),
            "ok": count_ok,
        }
        if expected is None:
            count_result["error"] = "manifest has no count_only_result"
            failures.append({
                "check": "count",
                "error": "count check requested but manifest has no expected count",
            })
        elif not count_ok:
            failures.append({
                "check": "count",
                "error": f"got {len(files)}, expected {expected}",
            })
        result["expected_count"] = expected
        result["count_ok"] = count_ok
        write_json(os.path.join(checks_dir, "count.json"), count_result)

    if check_formula:
        rows = []
        vacancy_element = element_symbol(vacancy_symbol)
        for path in files:
            data = read_cif(path)
            has_vacancy = any(
                element_symbol(a.type_symbol) == vacancy_element
                for a in data.atoms
            )
            if has_vacancy:
                result["vacancy_errors"].append(path)
            try:
                formula = assert_formula(data.atoms, target_formula, vacancy_symbol)
                rows.append({
                    "path": path,
                    "formula": formula,
                    "has_vacancy": has_vacancy,
                    "ok": not has_vacancy,
                })
            except ValueError as exc:
                rows.append({
                    "path": path,
                    "formula": None,
                    "has_vacancy": has_vacancy,
                    "ok": False,
                    "error": str(exc),
                })
                result["formula_errors"].append({"path": path, "error": str(exc)})
        write_jsonl(os.path.join(checks_dir, "formula_check.jsonl"), rows)
        if result["vacancy_errors"] or result["formula_errors"]:
            failures.append({
                "check": "formula",
                "error": (
                    f"{len(result['formula_errors'])} formula error(s), "
                    f"{len(result['vacancy_errors'])} residual-vacancy file(s)"
                ),
            })

    if check_dedup:
        try:
            dedup = verify_no_duplicates(files)
            result["dedup"] = dedup
            write_json(os.path.join(checks_dir, "dedup_check.json"), dedup)
            if dedup["has_duplicates"]:
                failures.append({
                    "check": "dedup",
                    "error": (
                        "StructureMatcher-confirmed duplicates in "
                        f"{len(dedup['duplicates'])} safe bucket(s)"
                    ),
                })
        except Exception as exc:  # noqa: BLE001 - preserve all audit results
            payload = {"ok": False, "error": str(exc)}
            result["dedup"] = payload
            write_json(os.path.join(checks_dir, "dedup_check.json"), payload)
            failures.append({"check": "dedup", "error": str(exc)})

    if symprec_list:
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        angle = manifest.get("angle_tolerance", 5.0) if manifest else 5.0
        atol = manifest.get("atol", 1e-5) if manifest else 1e-5
        scan = []
        for symprec in symprec_list:
            try:
                payload = count_shry_structures(
                    input_cif,
                    scaling_matrix=matrix,
                    symprec=float(symprec),
                    angle_tolerance=float(angle),
                    atol=float(atol),
                    out_json=os.path.join(checks_dir, f"count_symprec_{symprec}.json"),
                    backend=backend,
                )
                scan.append({
                    "symprec": float(symprec),
                    "count": payload["count_only_result"],
                    "ok": True,
                })
            except Exception as exc:  # noqa: BLE001 - report entire scan
                scan.append({
                    "symprec": float(symprec),
                    "count": None,
                    "ok": False,
                    "error": str(exc),
                })
        result["symprec_scan"] = scan
        counts = [row["count"] for row in scan if row["ok"]]
        stable = len(counts) == len(scan) and len(set(counts)) <= 1
        expected = manifest.get("count_only_result") if manifest else None
        matches_expected = (
            expected is None
            or (len(counts) == len(scan) and all(count == expected for count in counts))
        )
        scan_payload = {
            "results": scan,
            "stable": stable,
            "expected": expected,
            "matches_expected": matches_expected,
        }
        write_json(os.path.join(checks_dir, "symprec_scan.json"), scan_payload)
        if not stable:
            failures.append({
                "check": "symprec_scan",
                "error": "SHRY counts are not stable across requested tolerances",
            })
        if not matches_expected:
            failures.append({
                "check": "symprec_scan",
                "error": f"one or more scan counts differ from expected {expected}",
            })

    if cross_backend:
        if cross_backend != "supercell":
            raise ValueError(f"unknown cross backend {cross_backend!r}")
        input_cif = _manifest_input_cif(output_dir, manifest)
        matrix = manifest.get("supercell_matrix") if manifest else None
        try:
            count = run_supercell_count(input_cif, scaling_matrix=matrix)
            expected = manifest.get("count_only_result") if manifest else None
            cross = {
                "backend": "supercell",
                "count": count,
                "expected": expected,
                "ok": expected is None or count == expected,
            }
            if not cross["ok"]:
                failures.append({
                    "check": "cross_backend",
                    "error": f"supercell count {count} != expected {expected}",
                })
        except Exception as exc:  # noqa: BLE001 - preserve audit artifact
            cross = {"backend": "supercell", "count": None,
                     "ok": False, "error": str(exc)}
            failures.append({"check": "cross_backend", "error": str(exc)})
        result["cross_backend"] = cross
        write_json(os.path.join(checks_dir, "cross_backend.json"), cross)

    if check_degeneracy:
        matrix = validate_scaling_matrix(
            manifest.get("supercell_matrix")
            if manifest and manifest.get("supercell_matrix")
            else [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
        )
        determinant = scaling_determinant(matrix)
        expected_total = manifest.get("total_combinations") if manifest else None
        if determinant != 1:
            degeneracy_payload = {
                "degeneracy_sum": None,
                "supercell_determinant": determinant,
                "expected_total_combinations": expected_total,
                "ok": False,
                "error": (
                    "degeneracy calculation is not implemented rigorously for "
                    "non-unit supercells"
                ),
            }
            failures.append({
                "check": "degeneracy",
                "error": degeneracy_payload["error"],
            })
        else:
            try:
                input_cif = _manifest_input_cif(output_dir, manifest)
                declared = read_cif(input_cif).spacegroup_number
                detected = manifest.get("parent_spacegroup_detected") \
                    if manifest else None
                if detected is not None and declared != int(detected):
                    raise ValueError(
                        "degeneracy requires the prepared CIF's declared space "
                        "group to equal the spglib-detected enumeration group; "
                        f"declared {declared}, detected {detected}"
                    )
                if expected_total is None and manifest:
                    site_orbits = manifest.get("site_orbits")
                    if site_orbits:
                        expected_total = raw_combination_count(
                            input_cif, site_orbits,
                        )
                if expected_total is None:
                    raise ValueError(
                        "cannot verify degeneracy completeness: neither the "
                        "SHRY total-combination count nor audited site-orbit "
                        "metadata is available"
                    )
                total = degeneracy_sum(files, input_cif, vacancy_symbol=vacancy_symbol)
                result["degeneracy_sum"] = total
                degeneracy_payload = {
                    "degeneracy_sum": total,
                    "supercell_determinant": determinant,
                    "expected_total_combinations": expected_total,
                    "ok": expected_total is None or total == expected_total,
                }
                if not degeneracy_payload["ok"]:
                    failures.append({
                        "check": "degeneracy",
                        "error": f"degeneracy sum {total} != expected {expected_total}",
                    })
            except Exception as exc:  # noqa: BLE001 - preserve audit artifact
                degeneracy_payload = {
                    "degeneracy_sum": None,
                    "supercell_determinant": determinant,
                    "expected_total_combinations": expected_total,
                    "ok": False,
                    "error": str(exc),
                }
                failures.append({"check": "degeneracy", "error": str(exc)})
        write_json(os.path.join(checks_dir, "degeneracy.json"), degeneracy_payload)

    # Orbit grouping record: prefer the
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

    result["ok"] = not failures
    verify_path = os.path.join(checks_dir, "verify.json")
    write_json(verify_path, result)
    if failures:
        summary = "; ".join(
            f"{failure['check']}: {failure['error']}" for failure in failures
        )
        raise ValueError(f"SHRY verification failed ({summary}); see {verify_path}")
    return result


def _cif_files(path: str) -> list[str]:
    return sorted(
        os.path.join(path, name)
        for name in os.listdir(path)
        if name.lower().endswith(".cif")
        and os.path.isfile(os.path.join(path, name))
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
