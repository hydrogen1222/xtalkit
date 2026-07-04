"""High-throughput post-processing helpers for SHRY outputs."""

from __future__ import annotations

import json
import math
import os

import gemmi

from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.manifest import write_json
from xtalkit.enumeration.shry_verify import _cif_files


def shortest_distance(cif_path: str, species: tuple[str, str] | None = None) -> float | None:
    """Return the shortest periodic pair distance in Angstrom."""
    data = read_cif(cif_path)
    cell = gemmi.UnitCell(
        float(data.cell["a"]), float(data.cell["b"]), float(data.cell["c"]),
        float(data.cell["alpha"]), float(data.cell["beta"]), float(data.cell["gamma"]),
    )
    best: float | None = None
    atoms = data.atoms
    for i, a in enumerate(atoms):
        for b in atoms[i + 1:]:
            if species and tuple(sorted((a.type_symbol, b.type_symbol))) != tuple(sorted(species)):
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
    for path in _cif_files(clean_dir):
        rows.append({"path": path, "shortest_distance": shortest_distance(path, species)})
    rows.sort(key=lambda r: math.inf if r["shortest_distance"] is None else r["shortest_distance"])
    checks = os.path.join(output_dir, "checks")
    os.makedirs(checks, exist_ok=True)
    write_json(os.path.join(checks, "shortest_distance.json"), {"results": rows})
    return rows


def ewald_energy(cif_path: str, charges: dict[str, float]) -> float:
    """Compute Ewald energy using pymatgen, with explicit oxidation states."""
    try:
        from pymatgen.analysis.ewald import EwaldSummation
        from pymatgen.core import Structure
    except ImportError as exc:
        raise RuntimeError("pymatgen is required for Ewald scoring") from exc
    struct = Structure.from_file(cif_path)
    missing = sorted({sp.symbol for sp in struct.composition} - set(charges))
    if missing:
        raise ValueError(f"missing charges for species: {missing}")
    struct.add_oxidation_state_by_element(charges)
    return float(EwaldSummation(struct).total_energy)


def rank_by_ewald(output_dir: str, charges: dict[str, float]) -> list[dict]:
    """Rank clean CIFs by Ewald energy."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    rows = [{"path": path, "ewald_energy": ewald_energy(path, charges)}
            for path in _cif_files(clean_dir)]
    rows.sort(key=lambda r: r["ewald_energy"])
    checks = os.path.join(output_dir, "checks")
    os.makedirs(checks, exist_ok=True)
    write_json(os.path.join(checks, "ewald.json"), {"charges": charges, "results": rows})
    return rows


def write_tblite_inputs(output_dir: str, charge: int = 0, uhf: int = 0) -> list[str]:
    """Write simple tblite job folders with XYZ inputs."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    root = os.path.join(output_dir, "tblite")
    os.makedirs(root, exist_ok=True)
    paths = []
    for i, cif in enumerate(_cif_files(clean_dir), 1):
        job = os.path.join(root, f"conf_{i:06d}")
        os.makedirs(job, exist_ok=True)
        xyz = os.path.join(job, "structure.xyz")
        _write_xyz_from_cif(cif, xyz)
        run = os.path.join(job, "run.sh")
        with open(run, "w", encoding="utf-8") as f:
            f.write("#!/usr/bin/env bash\nset -euo pipefail\n")
            f.write(f"tblite structure.xyz --chrg {charge} --uhf {uhf} > tblite.out\n")
        os.chmod(run, 0o755)
        paths.append(job)
    return paths


def write_cp2k_inputs(output_dir: str, template: str | None = None) -> list[str]:
    """Write CP2K job folders with structure CIF and minimal input."""
    clean_dir = os.path.join(output_dir, "clean_cif")
    root = os.path.join(output_dir, "cp2k")
    os.makedirs(root, exist_ok=True)
    template_text = None
    if template:
        with open(template, encoding="utf-8") as f:
            template_text = f.read()
    paths = []
    for i, cif in enumerate(_cif_files(clean_dir), 1):
        job = os.path.join(root, f"conf_{i:06d}")
        os.makedirs(job, exist_ok=True)
        local_cif = os.path.join(job, "structure.cif")
        with open(cif, encoding="utf-8") as src, open(local_cif, "w", encoding="utf-8") as dst:
            dst.write(src.read())
        inp = os.path.join(job, "cp2k.inp")
        with open(inp, "w", encoding="utf-8") as f:
            f.write(template_text or _default_cp2k_input())
        paths.append(job)
    return paths


def write_slurm_array(output_dir: str, job_kind: str = "tblite") -> str:
    """Write a simple Slurm array script for generated job folders."""
    job_root = os.path.join(output_dir, job_kind)
    jobs = sorted(d for d in os.listdir(job_root)
                  if os.path.isdir(os.path.join(job_root, d)))
    if not jobs:
        raise ValueError(f"no {job_kind} jobs found in {job_root}")
    script = os.path.join(output_dir, f"run_{job_kind}_array.slurm")
    command = "./run.sh" if job_kind == "tblite" else "cp2k -i cp2k.inp -o cp2k.out"
    with open(script, "w", encoding="utf-8") as f:
        f.write("#!/usr/bin/env bash\n")
        f.write(f"#SBATCH --array=1-{len(jobs)}\n")
        f.write("#SBATCH --job-name=xtalkit\nset -euo pipefail\n")
        f.write(f"jobs=({ ' '.join(jobs) })\n")
        f.write("job=${jobs[$SLURM_ARRAY_TASK_ID-1]}\n")
        f.write(f"cd {job_root}/$job\n{command}\n")
    return script


def _periodic_distance(cell, a, b) -> float:
    diff = [a[i] - b[i] for i in range(3)]
    diff = [x - round(x) for x in diff]
    pos = cell.orthogonalize(gemmi.Fractional(diff[0], diff[1], diff[2]))
    return math.sqrt(pos.x ** 2 + pos.y ** 2 + pos.z ** 2)


def _write_xyz_from_cif(cif_path: str, xyz_path: str) -> None:
    data = read_cif(cif_path)
    cell = gemmi.UnitCell(
        float(data.cell["a"]), float(data.cell["b"]), float(data.cell["c"]),
        float(data.cell["alpha"]), float(data.cell["beta"]), float(data.cell["gamma"]),
    )
    lines = [str(len(data.atoms)), "xtalkit shry"]
    for atom in data.atoms:
        pos = cell.orthogonalize(gemmi.Fractional(float(atom.x), float(atom.y), float(atom.z)))
        lines.append(f"{atom.type_symbol:<3} {pos.x:.10f} {pos.y:.10f} {pos.z:.10f}")
    with open(xyz_path, "w", encoding="utf-8") as f:
        f.write("\n".join(lines) + "\n")


def _default_cp2k_input() -> str:
    return """&GLOBAL
  PROJECT xtalkit
  RUN_TYPE ENERGY
&END GLOBAL
&FORCE_EVAL
  METHOD Quickstep
  &SUBSYS
    ! Replace this minimal template with your production CP2K setup.
  &END SUBSYS
&END FORCE_EVAL
"""
