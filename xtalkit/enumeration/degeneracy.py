"""Degeneracy checks for ordered configurations."""

from __future__ import annotations

import math
from collections import defaultdict
from fractions import Fraction
from math import factorial

from xtalkit.enumeration.cif_io import AtomRow, read_cif
from xtalkit.enumeration.formula import element_symbol


def compute_degeneracy(
    structure_cif: str,
    parent_cif: str,
    vacancy_symbol: str = "X",
    tol: float = 1e-4,
) -> int:
    """Return orbit degeneracy ``|G| / |stabilizer|`` for one unit-cell structure.

    Atom configurations are compared as species-labelled periodic multisets,
    using a bipartite match rather than coordinate binning.  The parent CIF's
    operations are normalized and checked for identity, closure, and agreement
    with the declared space-group order before the orbit--stabilizer theorem is
    applied.
    """
    import gemmi

    if not math.isfinite(tol) or tol <= 0:
        raise ValueError("degeneracy tolerance must be a finite positive number")
    parent = read_cif(parent_cif)
    structure = read_cif(structure_cif)
    _require_same_unit_cell(parent, structure)
    if not parent.symops:
        raise ValueError("parent CIF has no symmetry operations")
    ops = _validated_group_operations(parent, gemmi)
    original = _configuration(structure.atoms, vacancy_symbol, tol)
    stabilizer = 0
    for op in ops:
        transformed = _configuration(
            _transform_atoms(structure.atoms, op), vacancy_symbol, tol,
        )
        if _configurations_match(original, transformed, tol):
            stabilizer += 1
    if stabilizer == 0:
        raise ValueError("no stabilizer operations found; check parent symmetry/tolerance")
    if len(ops) % stabilizer:
        raise ValueError(
            "stabilizer size does not divide parent group order; the structure "
            "or tolerance does not define a consistent group action"
        )
    return len(ops) // stabilizer


def _require_same_unit_cell(parent, structure) -> None:
    """Require the unit-cell action to refer to the same physical lattice."""
    for name in ("a", "b", "c"):
        expected = float(parent.cell[name])
        actual = float(structure.cell[name])
        if not math.isclose(actual, expected, rel_tol=1e-8, abs_tol=1e-8):
            raise ValueError(
                f"ordered and parent CIF unit cells differ at {name}: "
                f"{actual} != {expected}"
            )
    for name in ("alpha", "beta", "gamma"):
        expected = float(parent.cell[name])
        actual = float(structure.cell[name])
        if not math.isclose(actual, expected, rel_tol=0.0, abs_tol=1e-7):
            raise ValueError(
                f"ordered and parent CIF unit cells differ at {name}: "
                f"{actual} != {expected}"
            )


def degeneracy_sum(
    structure_cifs: list[str],
    parent_cif: str,
    vacancy_symbol: str = "X",
    tol: float = 1e-4,
) -> int:
    """Sum degeneracies over a generated structure set."""
    return sum(compute_degeneracy(path, parent_cif, vacancy_symbol, tol)
               for path in structure_cifs)


def raw_combination_count(parent_cif: str, site_orbits: list[dict]) -> int:
    """Return the exact number of labelled occupational configurations.

    For each independently substituted orbit this is the multinomial
    ``n! / product(n_i!)``; independent orbit factors are multiplied.  The
    input must be a prepared CIF whose representative occupancies are exactly
    integerizable in the unit cell.
    """
    parent = read_cif(parent_cif)
    atoms_by_coord: dict[tuple[str, str, str], list[AtomRow]] = {}
    for atom in parent.atoms:
        atoms_by_coord.setdefault(atom.frac_key, []).append(atom)
    total = 1
    for record in site_orbits:
        representative = tuple(record["representative"])
        atoms = atoms_by_coord.get(representative)
        if not atoms:
            raise ValueError(
                f"orbit {record.get('orbit_id')} representative "
                f"{representative} is absent from the prepared CIF"
            )
        multiplicity = int(record["multiplicity"])
        occupancies: dict[str, Fraction] = defaultdict(Fraction)
        for atom in atoms:
            occupancies[element_symbol(atom.type_symbol)] += Fraction(
                str(atom.occupancy)
            )
        counts = {species: occupancy * multiplicity
                  for species, occupancy in occupancies.items()}
        if any(count.denominator != 1 for count in counts.values()):
            raise ValueError(
                f"orbit {record.get('orbit_id')} is not integerizable in the "
                f"unit cell: {counts}"
            )
        integer_counts = [count.numerator for count in counts.values()]
        if sum(integer_counts) != multiplicity:
            raise ValueError(
                f"orbit {record.get('orbit_id')} species counts sum to "
                f"{sum(integer_counts)}, expected multiplicity {multiplicity}"
            )
        factor = factorial(multiplicity)
        for count in integer_counts:
            factor //= factorial(count)
        total *= factor
    return total


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


def _configuration(
    atoms: list[AtomRow], vacancy_symbol: str, tol: float,
) -> dict[str, tuple[tuple[float, float, float], ...]]:
    """Return a validated species-labelled periodic multiset."""
    vacancy = element_symbol(vacancy_symbol)
    rows: dict[str, list[tuple[float, float, float]]] = defaultdict(list)
    for atom in atoms:
        species = element_symbol(atom.type_symbol)
        if species == vacancy:
            continue
        if Fraction(str(atom.occupancy)) != 1:
            raise ValueError(
                "degeneracy requires an ordered structure with unit occupancy; "
                f"got {atom.label} {atom.type_symbol} occ={atom.occupancy}"
            )
        coords = tuple(float(v) % 1.0 for v in (atom.x, atom.y, atom.z))
        rows[species].append(coords)

    flattened = [(species, coord) for species, coords in rows.items() for coord in coords]
    for index, (species_a, coord_a) in enumerate(flattened):
        for species_b, coord_b in flattened[index + 1:]:
            if _periodic_close(coord_a, coord_b, tol):
                raise ValueError(
                    "ordered structure has overlapping atom sites within the "
                    f"degeneracy tolerance: {species_a} and {species_b}"
                )
    return {species: tuple(coords) for species, coords in rows.items()}


def _configurations_match(
    first: dict[str, tuple[tuple[float, float, float], ...]],
    second: dict[str, tuple[tuple[float, float, float], ...]],
    tol: float,
) -> bool:
    """Test equality of periodic multisets with a complete bipartite matching."""
    if first.keys() != second.keys():
        return False
    for species, left in first.items():
        right = second[species]
        if len(left) != len(right):
            return False
        edges = [
            [j for j, candidate in enumerate(right)
             if _periodic_close(point, candidate, tol)]
            for point in left
        ]
        matched_left = [-1] * len(right)

        def augment(i: int, seen: set[int]) -> bool:
            for j in edges[i]:
                if j in seen:
                    continue
                seen.add(j)
                if matched_left[j] < 0 or augment(matched_left[j], seen):
                    matched_left[j] = i
                    return True
            return False

        if any(not augment(i, set()) for i in range(len(left))):
            return False
    return True


def _periodic_close(first, second, tol: float) -> bool:
    return all(
        min(abs(a - b), 1.0 - abs(a - b)) <= tol
        for a, b in zip(first, second)
    )


def _validated_group_operations(parent, gemmi) -> list:
    """Normalize CIF operations and verify that they form the declared group."""
    unique = {}
    for expression in parent.symops:
        try:
            op = gemmi.Op(expression).wrap()
        except Exception as exc:
            raise ValueError(f"invalid parent symmetry operation {expression!r}") from exc
        unique.setdefault(_op_key(op), op)
    ops = list(unique.values())

    expected_order = _declared_group_order(parent, gemmi)
    if expected_order is not None and len(ops) != expected_order:
        raise ValueError(
            f"parent CIF contains {len(ops)} unique symmetry operations, but "
            f"the declared space group has order {expected_order}"
        )

    keys = set(unique)
    identity = _op_key(gemmi.Op("x,y,z"))
    if identity not in keys:
        raise ValueError("parent symmetry operations do not contain the identity")
    for left in ops:
        for right in ops:
            if _op_key((left * right).wrap()) not in keys:
                raise ValueError(
                    "parent symmetry operations are not closed under composition"
                )
    return ops


def _declared_group_order(parent, gemmi) -> int | None:
    group = None
    if parent.spacegroup_name:
        try:
            group = gemmi.SpaceGroup(str(parent.spacegroup_name).strip("'\""))
        except Exception as exc:
            raise ValueError(
                f"invalid declared space-group name {parent.spacegroup_name!r}"
            ) from exc
    if parent.spacegroup_number is not None:
        if not 1 <= parent.spacegroup_number <= 230:
            raise ValueError(
                f"invalid declared space-group number {parent.spacegroup_number}"
            )
        numbered = gemmi.find_spacegroup_by_number(parent.spacegroup_number)
        if numbered is None:
            raise ValueError(
                f"unknown declared space-group number {parent.spacegroup_number}"
            )
        if group is not None and group.number != numbered.number:
            raise ValueError(
                "declared space-group name and International Tables number disagree"
            )
        group = group or numbered
    return len(group.operations()) if group is not None else None


def _op_key(op) -> tuple[int, ...]:
    wrapped = op.wrap()
    return tuple(value for row in wrapped.rot for value in row) + tuple(wrapped.tran)
