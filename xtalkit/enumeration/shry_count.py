"""SHRY count-only workflow."""

from __future__ import annotations

import os

from xtalkit.enumeration.manifest import base_manifest, sha256_file, write_json
from xtalkit.enumeration.shry_backend import ShryBackend, matrix_args, parse_count_output


def count_shry_structures(
    cif_path: str,
    scaling_matrix: list[list[int]] | None = None,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
    atol: float = 1e-5,
    strict: bool = True,
    symmetrize: bool = False,
    out_json: str | None = None,
    backend: ShryBackend | None = None,
) -> dict:
    """Run ``shry --count-only`` and write a count manifest."""
    if not os.path.exists(cif_path):
        raise FileNotFoundError(cif_path)
    matrix = scaling_matrix or [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    backend = backend or ShryBackend()
    args = [
        cif_path,
        "--count-only",
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
    count = parse_count_output(result.stdout, result.stderr)
    payload = base_manifest(
        mode="strict_count_only" if strict else "count_only",
        input_cif=os.path.abspath(cif_path),
        input_sha256=sha256_file(cif_path),
        supercell_matrix=matrix,
        symprec=symprec,
        angle_tolerance=angle_tolerance,
        atol=atol,
        symmetrize_flag=bool(symmetrize),
        count_only_result=count,
        shry_command_line=" ".join(result.command),
        backend_version=backend.version(),
    )
    path = out_json or "count.json"
    write_json(path, payload)
    return payload
