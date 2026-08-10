"""Symmetry and site-grouping checks for strict SHRY input.

The parent space group is detected *independently* with spglib — not merely
read from the CIF declaration — by expanding the asymmetric unit to the full
cell (using the CIF's declared symmetry operations) and assigning each site a
synthetic species derived from its disorder composition signature. This makes
two sites symmetry-equivalent iff they sit at symmetry-related positions *and*
carry the same disorder composition, which is exactly the equivalence SHRY
relies on. The detected group is then compared against the user-declared
``--parent-spacegroup``, and the CIF's ``_atom_site_label`` grouping is checked
against the spglib orbit grouping: an orbit that spans several CIF coord-groups
would silently change SHRY's enumeration group and must be refused (or repaired
with ``--symmetrize``).
"""

from __future__ import annotations

from collections import defaultdict
from dataclasses import dataclass

import gemmi

from xtalkit.enumeration.cif_io import AtomRow, CifData
from xtalkit.enumeration.formula import element_symbol
from xtalkit.enumeration.occupancy import fraction_string, parse_fraction


@dataclass
class OrbitRecord:
    orbit_id: int
    labels: list[str]
    species: list[str]
    base_labels: list[str]
    representative: tuple[str, str, str]
    multiplicity: int
    wyckoff: str | None = None
    coord_group_count: int = 1
    split: bool = False

    def as_dict(self) -> dict:
        return {
            "orbit_id": self.orbit_id,
            "labels": self.labels,
            "species": self.species,
            "base_labels": self.base_labels,
            "representative": list(self.representative),
            "multiplicity": self.multiplicity,
            "wyckoff": self.wyckoff,
            "coord_group_count": self.coord_group_count,
            "split": self.split,
        }


@dataclass
class _FullAtom:
    """One atom of the expanded full cell, remembering its AU origin."""

    x: float
    y: float
    z: float
    species_number: int
    coord_keys: set[tuple[str, str, str]]


def _field(dataset, name):
    """Read a field from a spglib dataset (attribute in >=2.x, key in <2.x)."""
    if hasattr(dataset, name):
        return getattr(dataset, name)
    return dataset.get(name) if isinstance(dataset, dict) else None


def _detect_dataset(cell, symprec: float, angle_tolerance: float):
    """Run spglib.get_symmetry_dataset using the modern (non-deprecated) path.

    spglib >=2.7 deprecates silent ``None``-on-failure (``_throw=False``),
    emitting a ``DeprecationWarning``. Passing ``_throw=True`` makes spglib
    raise ``SpglibError`` instead — the documented modern path — which we
    catch as "no group detected". Older spglib builds without the kwarg fall
    back to the classic call.
    """
    import spglib
    try:
        try:
            return spglib.get_symmetry_dataset(
                cell, symprec=symprec, angle_tolerance=angle_tolerance,
                _throw=True,
            )
        except TypeError:
            return spglib.get_symmetry_dataset(
                cell, symprec=symprec, angle_tolerance=angle_tolerance,
            )
    except spglib.SpglibError:
        return None


def audit_symmetry_grouping(
    data: CifData,
    expected_spacegroup: int | None = None,
    symprec: float = 0.01,
    angle_tolerance: float = 5.0,
    symmetrize: bool = False,
) -> dict:
    """Detect the parent space group with spglib and audit label/orbit grouping.

    Detect the parent group from the full cell, compare it with the requested
    parent, record its Wyckoff orbits, and verify that CIF coordinate groups do
    not split a detected orbit. ``symmetrize`` explicitly delegates repair of
    an orbit split to SHRY downstream.
    """
    try:
        __import__("spglib")
    except ImportError as exc:  # pragma: no cover - depends on env
        raise RuntimeError(
            "spglib is required for the strict symmetry audit. "
            "Install it with `uv sync --extra enumerate`."
        ) from exc

    coord_to_atoms: dict[tuple[str, str, str], list[AtomRow]] = {}
    for atom in data.atoms:
        coord_to_atoms.setdefault(atom.frac_key, []).append(atom)

    full_atoms, lattice = _expand_full_cell(data, coord_to_atoms)
    if not full_atoms:
        raise ValueError("no atoms to audit after expanding the CIF")

    cell = (
        lattice,
        [[a.x, a.y, a.z] for a in full_atoms],
        [a.species_number for a in full_atoms],
    )
    dataset = _detect_dataset(cell, symprec, angle_tolerance)
    if dataset is None:
        raise ValueError(
            "spglib could not detect a space group; the input is too distorted "
            f"at symprec={symprec}, angle_tolerance={angle_tolerance}"
        )

    detected = int(_field(dataset, "number"))
    hall = _field(dataset, "hall")
    international = _field(dataset, "international")

    if expected_spacegroup is not None and detected != expected_spacegroup:
        raise ValueError(
            f"Detected space group {international} ({detected}) != declared "
            f"parent space group {expected_spacegroup}. "
            "The enumeration group defines the enumeration itself. "
            "Inspect the input symmetry and use the detected group only after "
            "confirming that it is the intended parent structure."
        )

    records = _orbit_records(full_atoms, dataset, coord_to_atoms)
    _check_label_orbit_consistency(records, symmetrize=symmetrize)

    return {
        "parent_spacegroup_detected": detected,
        "spacegroup_name": international,
        "spacegroup_declared": data.spacegroup_number,
        "hall_symbol": hall,
        "symmetrize_flag": bool(symmetrize),
        "symprec": symprec,
        "angle_tolerance": angle_tolerance,
        "site_orbits": [r.as_dict() for r in records],
    }


def _expand_full_cell(
    data: CifData,
    coord_to_atoms: dict[tuple[str, str, str], list[AtomRow]],
) -> tuple[list[_FullAtom], list[list[float]]]:
    """Expand the CIF asymmetric unit to the full cell.

    Sites sharing a fractional coordinate (a disorder group, e.g. Li + X, or
    Ge + P) are treated as one crystallographic site and assigned a synthetic
    species number from their composition signature, so spglib recognises them
    as a single site with one disorder composition.
    """
    sig_to_num: dict[tuple[tuple[str, str], ...], int] = {}
    au_sites: list[tuple[tuple[str, str, str], int]] = []
    for key, grp in coord_to_atoms.items():
        sig = _composition_signature(grp)
        if sig not in sig_to_num:
            sig_to_num[sig] = len(sig_to_num) + 1
        au_sites.append((key, sig_to_num[sig]))

    # No symmetry loop => the listed atoms are the full cell already (P1-expanded).
    ops = ([gemmi.Op(op) for op in data.symops]
           if data.symops else [gemmi.Op("x,y,z")])
    seen: dict[tuple[float, float, float], _FullAtom] = {}
    full: list[_FullAtom] = []
    for (key, num) in au_sites:
        x, y, z = (float(v) for v in key)
        for op in ops:
            cx, cy, cz = (float(v) % 1.0 for v in op.apply_to_xyz([x, y, z]))
            rk = tuple(round(value, 10) % 1.0 for value in (cx, cy, cz))
            if rk in seen:
                existing = seen[rk]
                if existing.species_number != num:
                    raise ValueError(
                        "symmetry expansion maps crystallographic sites with "
                        "different occupancy compositions onto the same "
                        f"fractional coordinate {rk}: "
                        f"{sorted(existing.coord_keys)} and {key}"
                    )
                existing.coord_keys.add(key)
                continue
            atom = _FullAtom(cx, cy, cz, num, {key})
            seen[rk] = atom
            full.append(atom)

    lattice = _lattice_matrix(data.cell)
    return full, lattice


def _composition_signature(atoms: list[AtomRow]) -> tuple[tuple[str, str], ...]:
    """Return an exact species/occupancy signature for one disordered site."""
    totals = defaultdict(lambda: 0)
    for atom in atoms:
        totals[element_symbol(atom.type_symbol)] += parse_fraction(atom.occupancy)
    return tuple(
        (symbol, fraction_string(occupancy))
        for symbol, occupancy in sorted(totals.items())
        if occupancy != 0
    )


def _lattice_matrix(cell: dict[str, str]) -> list[list[float]]:
    uc = gemmi.UnitCell(
        float(cell["a"]), float(cell["b"]), float(cell["c"]),
        float(cell["alpha"]), float(cell["beta"]), float(cell["gamma"]),
    )
    a = uc.orthogonalize(gemmi.Fractional(1, 0, 0))
    b = uc.orthogonalize(gemmi.Fractional(0, 1, 0))
    c = uc.orthogonalize(gemmi.Fractional(0, 0, 1))
    return [[a.x, a.y, a.z], [b.x, b.y, b.z], [c.x, c.y, c.z]]


def _orbit_records(
    full_atoms: list[_FullAtom],
    dataset,
    coord_to_atoms: dict[tuple[str, str, str], list[AtomRow]],
) -> list[OrbitRecord]:
    """Group full-cell atoms into orbits and map each back to its AU sites."""
    equiv = _field(dataset, "equivalent_atoms")
    wyckoffs = _field(dataset, "wyckoffs")
    if wyckoffs is None:
        wyckoffs = ["?"] * len(full_atoms)
    by_rep: dict[int, list[int]] = defaultdict(list)
    for i, rep in enumerate(equiv):
        by_rep[int(rep)].append(i)

    records: list[OrbitRecord] = []
    for orbit_id, (rep, idxs) in enumerate(sorted(by_rep.items()), 1):
        # Distinct CIF coord-groups that feed this orbit.
        coord_keys: list[tuple[str, str, str]] = []
        seen_keys: set[tuple[str, str, str]] = set()
        for i in idxs:
            for key in sorted(full_atoms[i].coord_keys):
                if key not in seen_keys:
                    seen_keys.add(key)
                    coord_keys.append(key)
        rep_key = coord_keys[0]
        labels = sorted({a.label for key in coord_keys for a in coord_to_atoms[key]})
        species = sorted({
            element_symbol(atom.type_symbol)
            for key in coord_keys
            for atom in coord_to_atoms[key]
        })
        # Base-label-set per coord-group (strip the _<species> suffix added on
        # vacancy fill / split, so Li1 and Li1_X share base "Li1"). An orbit is
        # "split" only if several coord-groups carry *different* base labels.
        coord_group_bases: list[frozenset[str]] = []
        for key in coord_keys:
            bases = _coord_group_base_labels(coord_to_atoms[key])
            coord_group_bases.append(bases)
        split = len(set(coord_group_bases)) > 1
        records.append(OrbitRecord(
            orbit_id=orbit_id,
            labels=labels,
            species=species,
            base_labels=sorted(set().union(*coord_group_bases)),
            representative=rep_key,
            multiplicity=len(idxs),
            wyckoff=wyckoffs[idxs[0]] if idxs else None,
            coord_group_count=len(coord_keys),
            split=split,
        ))
    return records


def _base_label(label: str, type_symbol: str) -> str:
    """Strip only the species suffix added by SHRY preparation.

    Arbitrary underscores in an original CIF label are preserved; splitting at
    the first underscore would incorrectly merge labels such as ``Li_A1`` and
    ``Li_A2`` into one SHRY disorder group.
    """
    suffix = f"_{element_symbol(type_symbol)}"
    return label[:-len(suffix)] if label.endswith(suffix) else label


def _coord_group_base_labels(atoms: list[AtomRow]) -> frozenset[str]:
    """Infer preparation-generated suffixes only when the whole site agrees.

    A lone original label such as ``site_Li`` is not evidence that ``_Li`` was
    generated.  For a true split/fill group (``M1_Ge``/``M1_P`` or
    ``Li1``/``Li1_X``), all rows reduce to the same base and can be grouped.
    """
    labels = [atom.label for atom in atoms]
    if len(labels) < 2:
        return frozenset(labels)
    candidates = [_base_label(atom.label, atom.type_symbol) for atom in atoms]
    if len(set(candidates)) == 1:
        return frozenset(candidates)
    return frozenset(labels)


def _check_label_orbit_consistency(records: list[OrbitRecord],
                                  symmetrize: bool) -> None:
    """Verify that CIF coordinate groups match the spglib orbit grouping.

    Two failure modes, both only fatal without ``--symmetrize``:
      * an orbit spans several CIF coord-groups with *different* base labels
        (the P1-expanded / mislabelled case — SHRY would split one Wyckoff
        orbit into independent groups). An orbit whose members share one label
        (a full-cell CIF) is fine.
      * a single CIF label spans several orbits (a label reused on inequivalent
        sites).
    """
    label_to_orbits: dict[str, set[int]] = defaultdict(set)
    for record in records:
        if record.split:
            _raise_or_warn(symmetrize,
                f"orbit {record.orbit_id} (Wyckoff {record.wyckoff}, mult "
                f"{record.multiplicity}) spans {record.coord_group_count} CIF "
                f"coord-groups with different labels {record.labels} — one "
                "Wyckoff orbit is split across several labels")
        for label in record.base_labels:
            label_to_orbits[label].add(record.orbit_id)

    bad = {label: sorted(orbits) for label, orbits in label_to_orbits.items()
           if len(orbits) > 1}
    if bad and not symmetrize:
        detail = "; ".join(f"{label}: orbits {orbits}" for label, orbits in bad.items())
        raise ValueError(
            "CIF site-label grouping does not match spglib symmetry orbits. "
            f"{detail}. Fix labels or rerun with --symmetrize."
        )


def _raise_or_warn(symmetrize: bool, detail: str) -> None:
    if symmetrize:
        return
    raise ValueError(
        "CIF site-label grouping does not match spglib symmetry orbits. "
        f"{detail}. Fix labels or rerun with --symmetrize."
    )
