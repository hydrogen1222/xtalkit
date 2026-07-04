"""Degeneracy checks for ordered configurations."""

from __future__ import annotations

from xtalkit.enumeration.cif_io import AtomRow, read_cif


def compute_degeneracy(
    structure_cif: str,
    parent_cif: str,
    vacancy_symbol: str = "X",
    tol: float = 1e-4,
) -> int:
    """Return orbit degeneracy ``|G| / |stabilizer|`` for one structure."""
    import gemmi

    parent = read_cif(parent_cif)
    structure = read_cif(structure_cif)
    if not parent.symops:
        raise ValueError("parent CIF has no symmetry operations")
    ops = [gemmi.Op(op) for op in parent.symops]
    original = _site_set(structure.atoms, vacancy_symbol, tol)
    stabilizer = 0
    for op in ops:
        transformed = _site_set(_transform_atoms(structure.atoms, op), vacancy_symbol, tol)
        if transformed == original:
            stabilizer += 1
    if stabilizer == 0:
        raise ValueError("no stabilizer operations found; check parent symmetry/tolerance")
    return len(ops) // stabilizer


def degeneracy_sum(
    structure_cifs: list[str],
    parent_cif: str,
    vacancy_symbol: str = "X",
    tol: float = 1e-4,
) -> int:
    """Sum degeneracies over a generated structure set."""
    return sum(compute_degeneracy(path, parent_cif, vacancy_symbol, tol)
               for path in structure_cifs)


def _transform_atoms(atoms: list[AtomRow], op) -> list[AtomRow]:
    out = []
    for atom in atoms:
        c = op.apply_to_xyz([float(atom.x), float(atom.y), float(atom.z)])
        out.append(AtomRow(
            atom.label,
            atom.type_symbol,
            str(float(c[0]) % 1.0),
            str(float(c[1]) % 1.0),
            str(float(c[2]) % 1.0),
            atom.occupancy,
        ))
    return out


def _site_set(atoms: list[AtomRow], vacancy_symbol: str, tol: float) -> tuple:
    q = max(tol, 1e-12)
    rows = []
    for atom in atoms:
        if atom.type_symbol == vacancy_symbol:
            continue
        coords = tuple(round((float(v) % 1.0) / q) for v in (atom.x, atom.y, atom.z))
        rows.append((atom.type_symbol, coords))
    return tuple(sorted(rows))
