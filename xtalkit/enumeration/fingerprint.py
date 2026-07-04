"""Cheap structure fingerprints for deduplication buckets."""

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
