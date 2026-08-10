"""Build a crystal-structure CIF from XRD-refinement parameters.

Given a space group, unit-cell parameters, and a list of atomic sites (each an
element on a Wyckoff position with its free-coordinate values and occupancy),
assemble a standards-compliant CIF. The asymmetric-unit representatives are
written together with the space group's symmetry operations; the reader
(VESTA, gemmi, pymatgen) expands them to the full unit cell.

This is the runtime counterpart to ``scripts/build_wyckoff_db.py``: the
coordinate templates stored in ``data/wyckoff.json`` use the canonical format
produced there (e.g. ``"x,-x,z"``, ``"1/3,2/3,z"``, ``"x,2x,z"``, ``"1/4+x"``),
and ``eval_coord`` below parses exactly that format.

Crystallographic input model
---------------------------
The independent inputs are:

* space group number (1-230) — the crystal system is *derived* from it and is
  used only to sanity-check the cell parameters, never asked separately;
* unit-cell parameters (a, b, c, alpha, beta, gamma);
* per atom: element, Wyckoff label (e.g. ``"16e"``), the free-coordinate
  values in order of appearance in the position's template, and occupancy
  (defaults to 1.0; partial/mixed occupancy is supported and yields a
  disordered CIF that ``xtalkit enumerate`` can later order).
"""

from __future__ import annotations

import math
import re
from collections import namedtuple

import gemmi

from xtalkit.spacegroup import wyckoff_positions, crystal_system
from xtalkit.utils import evaluate_coord_template, unit_cell_from_parameters

# A single atomic site: element symbol, full Wyckoff label (e.g. "16e"),
# the free-coordinate values in template order, and occupancy in (0, 1].
AtomSite = namedtuple("AtomSite", ["element", "wyckoff", "free", "occ"])

# A site specified by direct fractional coordinates (refinement-table style):
# element, fractional x/y/z, and occupancy. The Wyckoff orbit is detected.
FracAtom = namedtuple("FracAtom", ["element", "x", "y", "z", "occ"])


def _validated_element(symbol: str) -> gemmi.Element:
    """Return an exact real element, rejecting Gemmi's prefix fallbacks."""
    if not isinstance(symbol, str) or not re.fullmatch(r"[A-Z][a-z]?", symbol):
        raise ValueError(f"Unknown element {symbol!r} — use an exact element symbol.")
    try:
        element = gemmi.Element(symbol)
    except (RuntimeError, ValueError, TypeError) as exc:
        raise ValueError(f"Unknown element {symbol!r}: {exc}") from exc
    if element.atomic_number == 0:
        raise ValueError(f"Unknown element {symbol!r} — not a real element.")
    return element


def _require_consistent_cell(sg_number: int, cell_params: dict) -> gemmi.UnitCell:
    """Construct a valid cell and enforce the reference setting's metric."""
    cell = unit_cell_from_parameters(cell_params)
    issues = validate_cell(sg_number, cell_params)
    if issues:
        raise ValueError(
            "unit cell is inconsistent with the selected space group: "
            + "; ".join(issues)
        )
    return cell


def free_params(coord_str: str) -> list[str]:
    """Distinct free-parameter names (x/y/z) in order of first appearance.

    ``"x,-x,z"`` -> ``["x", "z"]``; ``"x,x,x"`` -> ``["x"]``;
    ``"1/3,2/3,z"`` -> ``["z"]``; ``"0,0,0"`` -> ``[]``.
    """
    seen: list[str] = []
    for ch in coord_str:
        if ch in "xyz" and ch not in seen:
            seen.append(ch)
    return seen


def _wrap01(x: float) -> float:
    """Wrap a coordinate into [0, 1)."""
    return x - math.floor(x)


def _min_image(d: float) -> float:
    """Wrap a difference into [-0.5, 0.5]."""
    return d - round(d)


def eval_coord(coord_str: str, values: list[float]) -> tuple[float, float, float]:
    """Evaluate a coordinate template at the given free-param values.

    Parses the canonical format from ``data/wyckoff.json``: each comma-
    separated axis is a sum of signed terms; a term is a rational constant or
    ``[coefficient]variable`` (e.g. ``x``, ``-x``, ``2x``, ``1/2-x``,
    ``1/4+x``). Values map to variables (x/y/z) in order of first appearance.
    Results are wrapped to [0, 1).
    """
    if len(values) != len(free_params(coord_str)):
        raise ValueError(
            f"coordinate {coord_str!r} expects "
            f"{len(free_params(coord_str))} free value(s), got {len(values)}."
        )
    name_to_val: dict[str, float] = {}
    for ch in coord_str:
        if ch in "xyz" and ch not in name_to_val:
            name_to_val[ch] = values[len(name_to_val)]

    return evaluate_coord_template(coord_str, name_to_val)


def find_wyckoff(sg_number: int, label: str):
    """Return the WyckoffInfo for ``label`` (e.g. "16e") in the space group.

    Raises ValueError if the label is not a valid Wyckoff position.
    """
    for wp in wyckoff_positions(sg_number):
        if wp.letter == label:
            return wp
    valid = ", ".join(wp.letter for wp in wyckoff_positions(sg_number))
    raise ValueError(
        f"Invalid Wyckoff label {label!r} for SG #{sg_number}. Valid: {valid}"
    )


def build_structure(
    sg_number: int,
    cell_params: dict,
    atoms: list[AtomSite],
) -> gemmi.Structure:
    """Assemble a gemmi.Structure (asymmetric unit) from refinement parameters.

    Each atom is placed at its Wyckoff representative (free values substituted
    into the coordinate template). The structure carries the cell and space
    group so the CIF writer can emit symmetry operations.

    Raises ValueError on invalid element, Wyckoff label, or free-value count.
    """
    sg = gemmi.SpaceGroup(sg_number)
    structure = gemmi.Structure()
    structure.cell = _require_consistent_cell(sg_number, cell_params)
    structure.spacegroup_hm = sg.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    elem_counter: dict[str, int] = {}
    for site in atoms:
        occ = float(site.occ)
        if not math.isfinite(occ) or not 0.0 < occ <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} on {site.wyckoff} must be in "
                f"(0, 1], got {site.occ!r}"
            )
        wp = find_wyckoff(sg_number, site.wyckoff)
        rep = eval_coord(wp.coordinates, list(site.free))
        element = _validated_element(site.element)

        elem_counter[site.element] = elem_counter.get(site.element, 0) + 1
        atom = gemmi.Atom()
        atom.name = f"{site.element}{elem_counter[site.element]}"
        atom.element = element
        atom.occ = occ
        frac = gemmi.Fractional(*rep)
        atom.pos = structure.cell.orthogonalize(frac)
        residue.add_atom(atom)

    occupancy_issues = validate_atoms(sg_number, atoms)
    if occupancy_issues:
        raise ValueError("invalid site occupancy: " + "; ".join(occupancy_issues))

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    return structure


def stoichiometry(sg_number: int, atoms: list[AtomSite]) -> dict[str, float]:
    """Per-element atom count per unit cell = sum(multiplicity * occupancy).

    Mixed/partial occupancy is handled naturally: two elements on the same
    Wyckoff site each contribute multiplicity * their occupancy.
    """
    counts: dict[str, float] = {}
    for site in atoms:
        _validated_element(site.element)
        occupancy = float(site.occ)
        if not math.isfinite(occupancy) or not 0.0 < occupancy <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} on {site.wyckoff} must be in "
                f"(0, 1], got {site.occ!r}"
            )
        wp = find_wyckoff(sg_number, site.wyckoff)
        counts[site.element] = counts.get(site.element, 0.0) \
            + wp.multiplicity * occupancy
    return counts


def format_formula(counts: dict[str, float]) -> str:
    """Format a composition dict as a formula string.

    Integer counts are normalized by their GCD (so Na4Cl4 -> NaCl); fractional
    counts (from partial occupancy) retain up to 8 significant decimal digits,
    in insertion order.
    """
    # Preserve enough digits for refined fractional occupancies while keeping
    # a compact display formula; treat only numerical-noise residuals as ints.
    numeric_counts = {}
    for element, count in counts.items():
        _validated_element(element)
        numeric = float(count)
        if not math.isfinite(numeric) or numeric <= 0:
            raise ValueError(
                f"formula count for {element} must be finite and positive"
            )
        numeric_counts[element] = numeric
    if not numeric_counts:
        raise ValueError("cannot format an empty chemical formula")
    int_counts = {el: int(round(c)) for el, c in numeric_counts.items()}
    # Never round a small positive fractional count down to zero. Counts near
    # a positive integer may differ only by accumulated floating-point noise.
    all_integer = all(
        int_counts[el] >= 1
        and abs(numeric_counts[el] - int_counts[el])
        <= 1e-10 * max(1.0, abs(numeric_counts[el]))
        for el in numeric_counts
    )

    if all_integer and int_counts:
        g = _gcd_all(int_counts.values())
        if g > 1:
            int_counts = {el: c // g for el, c in int_counts.items()}

    parts = []
    for el, c in int_counts.items():
        if all_integer:
            parts.append(el if c == 1 else f"{el}{c}")
        else:
            parts.append(f"{el}{numeric_counts[el]:.8g}")
    return "".join(parts)


def _gcd_all(values):
    from math import gcd
    from functools import reduce
    return reduce(gcd, values, 0)


# --- fractional-coordinate input mode -------------------------------------
#
# Refinement tables give final fractional coordinates, often as non-canonical
# representatives of a Wyckoff orbit (e.g. (0, 1/2, z) instead of the canonical
# (1/2, 0, 1/4+z) for SG 137 4d — they are related by a 4-fold rotation and
# expand to the same orbit). Converting those to the canonical free-parameter
# form by hand is error-prone, so `--atom-frac` accepts fractional coordinates
# directly and detects the orbit.

def _orbit_size(sg_number: int, rep: tuple[float, float, float]) -> int:
    """Number of distinct symmetry-equivalent points of ``rep`` (= multiplicity)."""
    sg = gemmi.SpaceGroup(sg_number)
    seen: list[tuple[float, float, float]] = []
    for op in sg.operations():
        c = op.apply_to_xyz([rep[0], rep[1], rep[2]])
        f = (c[0] - math.floor(c[0]), c[1] - math.floor(c[1]),
             c[2] - math.floor(c[2]))
        if not any(max(abs(_min_image(f[i] - s[i])) for i in range(3)) < 1e-3
                   for s in seen):
            seen.append(f)
    return len(seen)


def _point_in_orbit(sg_number: int, a, b) -> bool:
    """True if point ``b`` lies in the symmetry orbit of point ``a``."""
    sg = gemmi.SpaceGroup(sg_number)
    for op in sg.operations():
        c = op.apply_to_xyz([a[0], a[1], a[2]])
        f = (c[0] - math.floor(c[0]), c[1] - math.floor(c[1]),
             c[2] - math.floor(c[2]))
        if max(abs(_min_image(f[i] - b[i])) for i in range(3)) < 1e-3:
            return True
    return False


def detect_wyckoff(
    sg_number: int,
    x: float,
    y: float,
    z: float,
    tolerance: float = 1e-3,
) -> tuple[str, int]:
    """Return (label, multiplicity) of the Wyckoff orbit containing (x, y, z).

    The point may be any representative (canonical or not). For each Wyckoff
    position we ask: does some symmetry image of the point match the position's
    fixed-offset template (with the free variables solved consistently)? The
    first match wins. Returns ("?", multiplicity) if no bundled position matches.
    """
    rep = (_wrap01(x), _wrap01(y), _wrap01(z))
    from xtalkit.matcher import detect_wyckoff_position
    match = detect_wyckoff_position(
        rep,
        gemmi.SpaceGroup(sg_number),
        wyckoff_positions(sg_number),
        tolerance=tolerance,
    )
    if match is not None:
        return match.letter, match.multiplicity
    return "?", _orbit_size(sg_number, rep)


def build_structure_frac(sg_number: int, cell_params: dict,
                         atoms: list[FracAtom]) -> gemmi.Structure:
    """Assemble a gemmi.Structure from direct fractional coordinates.

    Each atom is placed at its given (x, y, z); the structure carries the cell
    and space group so the CIF writer emits symmetry operations. Raises
    ValueError on an unknown element.
    """
    sg = gemmi.SpaceGroup(sg_number)
    structure = gemmi.Structure()
    structure.cell = _require_consistent_cell(sg_number, cell_params)
    structure.spacegroup_hm = sg.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    elem_counter: dict[str, int] = {}
    for site in atoms:
        coords = (float(site.x), float(site.y), float(site.z))
        if not all(math.isfinite(value) for value in coords):
            raise ValueError(
                f"fractional coordinates for {site.element} must be finite"
            )
        occ = float(site.occ)
        if not math.isfinite(occ) or not 0.0 < occ <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} must be in (0, 1], "
                f"got {site.occ!r}"
            )
        element = _validated_element(site.element)
        elem_counter[site.element] = elem_counter.get(site.element, 0) + 1
        atom = gemmi.Atom()
        atom.name = f"{site.element}{elem_counter[site.element]}"
        atom.element = element
        atom.occ = occ
        frac = gemmi.Fractional(*(_wrap01(value) for value in coords))
        atom.pos = structure.cell.orthogonalize(frac)
        residue.add_atom(atom)

    occupancy_issues = validate_atoms_frac(sg_number, atoms)
    if occupancy_issues:
        raise ValueError("invalid site occupancy: " + "; ".join(occupancy_issues))

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    return structure


def stoichiometry_frac(sg_number: int, atoms: list[FracAtom]) -> dict[str, float]:
    """Per-element cell count = sum(multiplicity * occupancy), multiplicity
    computed by expanding each fractional coordinate under the SG operations."""
    counts: dict[str, float] = {}
    for site in atoms:
        _validated_element(site.element)
        coords = tuple(float(value) for value in (site.x, site.y, site.z))
        if not all(math.isfinite(value) for value in coords):
            raise ValueError(
                f"fractional coordinates for {site.element} must be finite"
            )
        occupancy = float(site.occ)
        if not math.isfinite(occupancy) or not 0.0 < occupancy <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} must be in (0, 1], "
                f"got {site.occ!r}"
            )
        mult = _orbit_size(sg_number, tuple(_wrap01(value) for value in coords))
        counts[site.element] = counts.get(site.element, 0.0) + mult * occupancy
    return counts


def validate_atoms_frac(sg_number: int, atoms: list[FracAtom]) -> list[str]:
    """Warn when symmetry-equivalent site rows sum to occupancy above one."""
    sites = []
    for site in atoms:
        coords = tuple(float(value) for value in (site.x, site.y, site.z))
        occupancy = float(site.occ)
        if not all(math.isfinite(value) for value in coords):
            raise ValueError(
                f"fractional coordinates for {site.element} must be finite"
            )
        if not math.isfinite(occupancy) or not 0.0 < occupancy <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} must be in (0, 1], "
                f"got {site.occ!r}"
            )
        sites.append((tuple(_wrap01(value) for value in coords), occupancy))
    groups = _group_site_occupancies(sg_number, sites)
    warnings = []
    for representative, total in groups:
        if total > 1.0 + 1e-10:
            warnings.append(
                f"site {tuple(round(v, 4) for v in representative)}: total "
                f"occupancy {total:g} exceeds 1.0 "
                f"(unphysical for disorder)."
            )
    return warnings


# Cell-parameter constraints per crystal system (the independent axes/angles
# that must be equal or 90/120 deg). Used only to warn, never to override.
_CELL_CONSTRAINTS = {
    "triclinic": (),
    "monoclinic": (("gamma", 90.0), ("alpha", 90.0)),  # unique axis b
    "orthorhombic": (("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0)),
    "tetragonal": (("a", "b"), ("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0)),
    "trigonal": (("a", "b"), ("alpha", 90.0), ("beta", 90.0), ("gamma", 120.0)),
    "hexagonal": (("a", "b"), ("alpha", 90.0), ("beta", 90.0), ("gamma", 120.0)),
    "cubic": (("a", "b"), ("a", "c"), ("alpha", 90.0), ("beta", 90.0), ("gamma", 90.0)),
}


def validate_cell(sg_number: int, cell_params: dict) -> list[str]:
    """Return warning strings for cell params inconsistent with the crystal system."""
    system = crystal_system(sg_number)
    constraints = _CELL_CONSTRAINTS.get(system, ())
    warnings = []
    for c in constraints:
        if isinstance(c, tuple) and len(c) == 2 and isinstance(c[1], float):
            name, expected = c
            if abs(float(cell_params[name]) - expected) > 1e-3:
                warnings.append(
                    f"{system}: {name} should be {expected:g} but is "
                    f"{cell_params[name]:g}."
                )
        elif isinstance(c, tuple) and len(c) == 2:
            n1, n2 = c
            if abs(float(cell_params[n1]) - float(cell_params[n2])) > 1e-3:
                warnings.append(
                    f"{system}: {n1} should equal {n2} "
                    f"({cell_params[n1]:g} vs {cell_params[n2]:g})."
                )
    return warnings


def validate_atoms(sg_number: int, atoms: list[AtomSite]) -> list[str]:
    """Return warning strings for occupancy issues.

    Atoms sharing a representative coordinate (mixed occupancy on one site)
    must have occupancies summing to <= 1.0. Sums below 1.0 are intentional
    partial occupancy (disorder — the feature's whole point) and are NOT
    flagged; only unphysical sums > 1.0 are warned about.
    """
    sites = []
    for site in atoms:
        occupancy = float(site.occ)
        if not math.isfinite(occupancy) or not 0.0 < occupancy <= 1.0:
            raise ValueError(
                f"occupancy for {site.element} on {site.wyckoff} must be in "
                f"(0, 1], got {site.occ!r}"
            )
        wp = find_wyckoff(sg_number, site.wyckoff)
        rep = eval_coord(wp.coordinates, list(site.free))
        sites.append((rep, occupancy))
    warnings = []
    for representative, total in _group_site_occupancies(sg_number, sites):
        if total > 1.0 + 1e-10:
            warnings.append(
                f"site {tuple(round(v, 4) for v in representative)}: total "
                f"occupancy {total:g} exceeds 1.0 "
                f"(unphysical for disorder)."
            )
    return warnings


def _group_site_occupancies(
    sg_number: int,
    sites: list[tuple[tuple[float, float, float], float]],
) -> list[list]:
    """Combine rows that denote the same orbit, not only identical triples."""
    groups: list[list] = []
    for representative, occupancy in sites:
        for group in groups:
            if _point_in_orbit(sg_number, group[0], representative):
                group[1] += occupancy
                break
        else:
            groups.append([representative, occupancy])
    return groups
