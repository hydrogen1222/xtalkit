"""Safe coarse bucketing plus ``StructureMatcher`` duplicate verification.

The old implementation used an absolute fractional-coordinate hash as a hard
gate before ``StructureMatcher``.  That hash changes under an origin shift or
symmetry operation, so genuine duplicates could land in different buckets and
never be compared.  A strict verifier must prefer a larger candidate bucket
over a false negative.

We therefore bucket only by invariants that every duplicate must share under
the configured matcher: element counts and number of sites.  The coordinate
fingerprint remains useful diagnostic metadata, but never decides that two
structures are different.  ``StructureMatcher`` compares the complete safe
bucket and is the sole authority for equivalence.
"""

from __future__ import annotations

import hashlib

from xtalkit.enumeration.cif_io import read_cif
from xtalkit.enumeration.formula import counts_from_atom_rows, formula_from_counts


def structure_fingerprint(cif_path: str, quantum: float = 1e-4) -> str:
    """Return a fast diagnostic hash of literal quantized coordinates.

    This hash is intentionally *not* used as a no-duplicate proof because it
    is not invariant to origin shifts or crystallographic basis changes.
    """
    data = read_cif(cif_path)
    formula = formula_from_counts(counts_from_atom_rows(data.atoms))
    atoms = []
    for atom in data.atoms:
        q = tuple(round(float(v) / quantum) for v in (atom.x, atom.y, atom.z))
        atoms.append((atom.type_symbol, q))
    payload = repr((formula, sorted(atoms))).encode("utf-8")
    return hashlib.sha256(payload).hexdigest()


def _safe_bucket_key(structure) -> tuple[tuple[tuple[str, float], ...], int]:
    """Return a no-false-negative coarse key for the strict matcher settings."""
    counts = structure.composition.element_composition.as_dict()
    composition = tuple(sorted(
        (str(element), round(float(count), 12))
        for element, count in counts.items()
    ))
    return composition, len(structure)


def verify_no_duplicates(
    structure_paths: list[str],
    stol: float = 0.3,
    quantum: float = 1e-4,
) -> dict:
    """Conservatively find all ``StructureMatcher``-confirmed duplicates.

    Literal coordinate fingerprints are reported for traceability but are not
    used to exclude comparisons.
    """
    try:
        from pymatgen.analysis.structure_matcher import StructureMatcher
        from pymatgen.core import Structure
    except ImportError as exc:  # pragma: no cover - depends on env
        raise RuntimeError(
            "pymatgen is required for structure duplicate verification. "
            "Install it with `uv sync --extra enumerate`."
        ) from exc

    buckets: dict[tuple, list[tuple[str, object]]] = {}
    for path in structure_paths:
        structure = Structure.from_file(path)
        buckets.setdefault(_safe_bucket_key(structure), []).append((path, structure))

    multi = {key: entries for key, entries in buckets.items() if len(entries) > 1}
    duplicates: list[dict] = []
    matcher = StructureMatcher(primitive_cell=False, scale=False, stol=stol)
    for key, entries in multi.items():
        paths = [entry[0] for entry in entries]
        structs = [entry[1] for entry in entries]
        groups, confirmed_pairs = _match_graph(structs, matcher)
        dup_groups = [
            [paths[i] for i in grp] for grp in groups if len(grp) > 1
        ]
        if dup_groups:
            duplicates.append({
                "safe_bucket": repr(key),
                "groups": dup_groups,
                "confirmed_pairs": [
                    [paths[i], paths[j]] for i, j in confirmed_pairs
                ],
                "literal_fingerprints": {
                    path: structure_fingerprint(path, quantum=quantum)
                    for path in paths
                },
            })

    return {
        "total": len(structure_paths),
        "num_buckets": len(buckets),
        "multi_bucket_count": len(multi),
        "duplicates": duplicates,
        "fingerprint_collisions": [],
        "bucket_strategy": "element_counts_and_site_count",
        "group_semantics": (
            "connected components of directly StructureMatcher-confirmed pairs; "
            "approximate matching need not be transitive, so consult "
            "confirmed_pairs for pairwise claims"
        ),
        "has_duplicates": len(duplicates) > 0,
    }


def _match_groups(structs: list, matcher) -> list[list[int]]:
    """Return connected components of all pairwise matcher-confirmed edges.

    Matcher tolerance can make approximate matching non-transitive. Comparing
    only each group's first representative can therefore miss a confirmed
    B-C duplicate when A matches B but not C. The strict verifier evaluates
    every pair in a safe bucket and reports graph components.
    """
    return _match_graph(structs, matcher)[0]


def _match_graph(structs: list, matcher) -> tuple[list[list[int]], list[tuple[int, int]]]:
    """Return matcher-edge components and every directly confirmed pair."""
    parents = list(range(len(structs)))
    confirmed_pairs: list[tuple[int, int]] = []

    def find(index: int) -> int:
        while parents[index] != index:
            parents[index] = parents[parents[index]]
            index = parents[index]
        return index

    def union(left: int, right: int) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            parents[right_root] = left_root

    for i, structure in enumerate(structs):
        for j in range(i + 1, len(structs)):
            if matcher.fit(structure, structs[j]):
                confirmed_pairs.append((i, j))
                union(i, j)

    components: dict[int, list[int]] = {}
    for index in range(len(structs)):
        components.setdefault(find(index), []).append(index)
    return list(components.values()), confirmed_pairs
