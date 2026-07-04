"""Cheap structure fingerprints for deduplication buckets + two-stage dedup.

Plan §10 mandates a two-level duplicate check instead of all-pairs
``StructureMatcher`` (which is O(N²) and infeasible for thousands of
structures):

  1. **Fingerprint hash bucketing (O(N)).** A cheap composition + quantized
     fractional-coordinate hash assigns each structure to a bucket.
  2. **``StructureMatcher`` only within buckets.** A bucket with more than one
     member is refined with pymatgen's ``StructureMatcher``: true duplicates are
     flagged, while a mere fingerprint collision (different structures that
     happen to hash equal) is recorded but not fatal.
"""

from __future__ import annotations

import hashlib

from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.formula import counts_from_atom_rows, formula_from_counts


def structure_fingerprint(cif_path: str, quantum: float = 1e-4) -> str:
    """Return a quantized composition + fractional-coordinate hash."""
    data = read_cif(cif_path)
    formula = formula_from_counts(counts_from_atom_rows(data.atoms))
    atoms = []
    for atom in data.atoms:
        q = tuple(round(float(v) / quantum) for v in (atom.x, atom.y, atom.z))
        atoms.append((atom.type_symbol, q))
    payload = repr((formula, sorted(atoms))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def verify_no_duplicates(
    structure_paths: list[str],
    stol: float = 0.3,
    quantum: float = 1e-4,
) -> dict:
    """Two-stage dedup (plan §10).

    Returns a dict with the bucket counts, the StructureMatcher-confirmed
    duplicate groups (fatal), and the fingerprint-only collisions (non-fatal).
    """
    buckets: dict[str, list[str]] = {}
    for path in structure_paths:
        buckets.setdefault(structure_fingerprint(path, quantum=quantum), []).append(path)

    multi = {fp: paths for fp, paths in buckets.items() if len(paths) > 1}
    duplicates: list[dict] = []
    collisions: list[dict] = []

    if multi:
        try:
            from pymatgen.analysis.structure_matcher import StructureMatcher
            from pymatgen.core import Structure
        except ImportError as exc:  # pragma: no cover - depends on env
            raise RuntimeError(
                "pymatgen is required for the §10 second-stage dedup. "
                "Install it with `uv sync --extra enumerate`."
            ) from exc

        matcher = StructureMatcher(primitive_cell=False, scale=False, stol=stol)
        for fp, paths in multi.items():
            structs = [Structure.from_file(p) for p in paths]
            groups = _match_groups(structs, matcher)
            dup_groups = [
                [paths[i] for i in grp] for grp in groups if len(grp) > 1
            ]
            if dup_groups:
                duplicates.append({"fingerprint": fp, "groups": dup_groups})
            else:
                collisions.append({"fingerprint": fp, "paths": paths})

    return {
        "total": len(structure_paths),
        "num_buckets": len(buckets),
        "multi_bucket_count": len(multi),
        "duplicates": duplicates,
        "fingerprint_collisions": collisions,
        "has_duplicates": len(duplicates) > 0,
    }


def _match_groups(structs: list, matcher) -> list[list[int]]:
    """Group structure indices by StructureMatcher equivalence (within bucket)."""
    assigned: list[int | None] = [None] * len(structs)
    groups: list[list[int]] = []
    for i, s in enumerate(structs):
        if assigned[i] is not None:
            continue
        group = [i]
        assigned[i] = len(groups)
        for j in range(i + 1, len(structs)):
            if assigned[j] is None and matcher.fit(s, structs[j]):
                group.append(j)
                assigned[j] = len(groups)
        groups.append(group)
    return groups
