"""Symmetry and site-grouping checks for strict SHRY input."""

from __future__ import annotations

from dataclasses import dataclass

import gemmi

from xtalkit.enumeration.cif_io import AtomRow, CifData


@dataclass
class OrbitRecord:
    orbit_id: int
    labels: list[str]
    species: list[str]
    representative: tuple[str, str, str]
    multiplicity: int

    def as_dict(self) -> dict:
        return {
            "orbit_id": self.orbit_id,
            "labels": self.labels,
            "species": self.species,
            "representative": list(self.representative),
            "multiplicity": self.multiplicity,
        }


def audit_symmetry_grouping(
    data: CifData,
    expected_spacegroup: int | None = None,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
    symmetrize: bool = False,
) -> dict:
    """Check space group declaration and record symmetry orbits.

    The first implementation intentionally uses Gemmi's declared space-group
    operations instead of rewriting the CIF with a new detected group. This
    keeps prepare deterministic and preserves the user's crystallographic
    setting; strict mode still refuses an explicit parent-spacegroup mismatch.
    """
    del symprec, angle_tolerance  # recorded by caller; Gemmi ops are exact here.

    if expected_spacegroup is not None and data.spacegroup_number != expected_spacegroup:
        raise ValueError(
            f"Detected/declared space group {data.spacegroup_number} != "
            f"declared parent space group {expected_spacegroup}. "
            "The enumeration group defines the enumeration itself."
        )
    if not data.symops:
        raise ValueError(
            "CIF carries no symmetry operations; SHRY strict mode refuses P1-like "
            "or symmetry-degraded input."
        )

    records = _orbit_records(data)
    _check_label_does_not_span_orbits(records, symmetrize=symmetrize)
    return {
        "parent_spacegroup_detected": data.spacegroup_number,
        "spacegroup_name": data.spacegroup_name,
        "hall_symbol": None,
        "symmetrize_flag": bool(symmetrize),
        "site_orbits": [r.as_dict() for r in records],
    }


def _orbit_records(data: CifData) -> list[OrbitRecord]:
    coord_groups: dict[tuple[str, str, str], list[AtomRow]] = {}
    for atom in data.atoms:
        coord_groups.setdefault(atom.frac_key, []).append(atom)

    ops = [gemmi.Op(op) for op in data.symops]
    seen: set[tuple[str, str, str]] = set()
    records: list[OrbitRecord] = []
    for key, atoms in coord_groups.items():
        if key in seen:
            continue
        orbit = _orbit_keys(key, ops)
        matching_keys = sorted(k for k in coord_groups if k in orbit)
        seen.update(matching_keys)
        orbit_atoms = [atom for k in matching_keys for atom in coord_groups[k]]
        records.append(OrbitRecord(
            orbit_id=len(records) + 1,
            labels=[a.label for a in orbit_atoms],
            species=sorted({a.type_symbol for a in orbit_atoms}),
            representative=key,
            multiplicity=len(orbit),
        ))
    return records


def _orbit_keys(key: tuple[str, str, str], ops: list[gemmi.Op]) -> set[tuple[str, str, str]]:
    xyz = [float(key[0]), float(key[1]), float(key[2])]
    out: set[tuple[str, str, str]] = set()
    for op in ops:
        c = op.apply_to_xyz(xyz)
        out.add(tuple(f"{(float(v) % 1.0):.6f}" for v in c))
    return out


def _check_label_does_not_span_orbits(records: list[OrbitRecord],
                                      symmetrize: bool) -> None:
    label_to_orbits: dict[str, set[int]] = {}
    for record in records:
        for label in record.labels:
            label_to_orbits.setdefault(label, set()).add(record.orbit_id)
    bad = {label: sorted(orbits) for label, orbits in label_to_orbits.items()
           if len(orbits) > 1}
    if bad and not symmetrize:
        detail = "; ".join(f"{label}: orbits {orbits}" for label, orbits in bad.items())
        raise ValueError(
            "CIF site-label grouping does not match symmetry orbits. "
            f"{detail}. Fix labels or rerun with --symmetrize."
        )
