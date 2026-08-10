"""High-throughput post-processing helpers for SHRY outputs."""

from __future__ import annotations

import math
import os
import shlex

import gemmi

from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.formula import (
    counts_from_atom_rows,
    real_element_symbol,
)
from xtalkit.enumeration.manifest import write_json
from xtalkit.enumeration.shry_verify import _cif_files


def shortest_distance(cif_path: str, species: tuple[str, str] | None = None) -> float | None:
    """Return the shortest periodic pair distance in Angstrom."""
    data = read_cif(cif_path)
    counts_from_atom_rows(data.atoms)
    normalized_species = (
        tuple(sorted(real_element_symbol(value) for value in species))
        if species is not None else None
    )
    normalized_atoms = [
        (atom, real_element_symbol(atom.type_symbol)) for atom in data.atoms
    ]
    cell = gemmi.UnitCell(
        float(data.cell["a"]), float(data.cell["b"]), float(data.cell["c"]),
        float(data.cell["alpha"]), float(data.cell["beta"]), float(data.cell["gamma"]),
    )
    best: float | None = None
    for i, (a, a_species) in enumerate(normalized_atoms):
        for b, b_species in normalized_atoms[i + 1:]:
            pair = tuple(sorted((a_species, b_species)))
            if normalized_species is not None and pair != normalized_species:
                continue
            d = _periodic_distance(cell, (float(a.x), float(a.y), float(a.z)),
                                   (float(b.x), float(b.y), float(b.z)))
            if best is None or d < best:
                best = d
    return best


def rank_by_shortest_distance(
    output_dir: str,
    species: tuple[str, str] | None = None,
) -> list[dict]:
    """Rank clean CIFs by shortest distance."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    rows = []
    for path in _required_cif_files(clean_dir):
        rows.append({"path": path, "shortest_distance": shortest_distance(path, species)})
    rows.sort(key=lambda r: math.inf if r["shortest_distance"] is None else r["shortest_distance"])
    checks = os.path.join(output_dir, "checks")
    os.makedirs(checks, exist_ok=True)
    write_json(os.path.join(checks, "shortest_distance.json"), {"results": rows})
    return rows


def ewald_energy(cif_path: str, charges: dict[str, float]) -> float:
    """Compute a charge-neutral ordered Ewald energy with shared safeguards."""
    try:
        from pymatgen.core import Structure
    except ImportError as exc:
        raise RuntimeError("pymatgen is required for Ewald scoring") from exc
    from xtalkit.ewald import compute_ewald_energy

    struct = Structure.from_file(cif_path)
    return compute_ewald_energy(struct, charges=charges)


def rank_by_ewald(output_dir: str, charges: dict[str, float]) -> list[dict]:
    """Rank clean CIFs by Ewald energy."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    rows = [{"path": path, "ewald_energy": ewald_energy(path, charges)}
            for path in _required_cif_files(clean_dir)]
    rows.sort(key=lambda r: r["ewald_energy"])
    checks = os.path.join(output_dir, "checks")
    os.makedirs(checks, exist_ok=True)
    write_json(os.path.join(checks, "ewald.json"), {"charges": charges, "results": rows})
    return rows


def write_tblite_inputs(output_dir: str, charge: int = 0, uhf: int = 0) -> list[str]:
    """Write periodic tblite jobs using POSCAR geometry and current CLI flags."""
    charge = _exact_integer(charge, "charge")
    uhf = _exact_integer(uhf, "spin")
    if uhf < 0:
        raise ValueError("spin (number of unpaired electrons) must be non-negative")
    clean_dir = os.path.join(output_dir, "clean_cif")
    cifs = _required_cif_files(clean_dir)
    root = os.path.join(output_dir, "tblite")
    _make_empty_job_root(root)
    paths = []
    for i, cif in enumerate(cifs, 1):
        job = os.path.join(root, f"conf_{i:06d}")
        os.makedirs(job, exist_ok=True)
        poscar = os.path.join(job, "POSCAR")
        from xtalkit.enumeration.shry_enum import _write_poscar
        _write_poscar(read_cif(cif), poscar)
        run = os.path.join(job, "run.sh")
        with open(run, "w", encoding="utf-8") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(
                f"tblite run --no-restart --charge {charge} "
                f"--spin {uhf} POSCAR > tblite.out\n"
            )
        os.chmod(run, 0o755)
        paths.append(job)
    return paths


def write_slurm_array(output_dir: str) -> str:
    """Write a Slurm array script for generated tblite job folders."""
    job_root = os.path.join(output_dir, "tblite")
    jobs = sorted(d for d in os.listdir(job_root)
                  if os.path.isdir(os.path.join(job_root, d)))
    if not jobs:
        raise ValueError(f"no tblite jobs found in {job_root}")
    script = os.path.join(output_dir, "run_tblite_array.slurm")
    with open(script, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write(f"#SBATCH --array=1-{len(jobs)}\n")
        f.write("#SBATCH --job-name=xtalkit\nset -euo pipefail\n")
        f.write(f"jobs=({ ' '.join(jobs) })\n")
        f.write("job=${jobs[$SLURM_ARRAY_TASK_ID-1]}\n")
        f.write(f"cd {shlex.quote(job_root)}/\"$job\"\n./run.sh\n")
    return script


def _periodic_distance(cell, a, b) -> float:
    return cell.find_nearest_pbc_image(
        gemmi.Fractional(*a), gemmi.Fractional(*b)
    ).dist()


def _make_empty_job_root(path: str) -> None:
    """Create a job root without mixing files from separate generations."""
    if os.path.exists(path):
        if not os.path.isdir(path):
            raise NotADirectoryError(path)
        if os.listdir(path):
            raise FileExistsError(
                f"job output directory is not empty: {path}. Choose a clean "
                "workflow directory or move the previous jobs first."
            )
    else:
        os.makedirs(path)


def _required_cif_files(path: str) -> list[str]:
    if not os.path.isdir(path):
        raise FileNotFoundError(f"clean CIF directory not found: {path}")
    files = _cif_files(path)
    if not files:
        raise ValueError(f"no clean CIF configurations found in {path}")
    return files


def _exact_integer(value, name: str) -> int:
    if isinstance(value, bool):
        raise ValueError(f"{name} must be an integer")
    try:
        numeric = float(value)
        integer = int(value)
    except (TypeError, ValueError, OverflowError) as exc:
        raise ValueError(f"{name} must be an integer") from exc
    if not math.isfinite(numeric) or numeric != integer:
        raise ValueError(f"{name} must be an integer")
    return integer
