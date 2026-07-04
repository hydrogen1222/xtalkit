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
from collections import namedtuple

import gemmi

from xtalkit.spacegroup import wyckoff_positions, crystal_system

# A single atomic site: element symbol, full Wyckoff label (e.g. "16e"),
# the free-coordinate values in template order, and occupancy in (0, 1].
AtomSite = namedtuple("AtomSite", ["element", "wyckoff", "free", "occ"])

# A site specified by direct fractional coordinates (refinement-table style):
# element, fractional x/y/z, and occupancy. The Wyckoff orbit is detected.
FracAtom = namedtuple("FracAtom", ["element", "x", "y", "z", "occ"])


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

    out: list[float] = []
    for part in coord_str.split(","):
        out.append(_wrap01(_eval_axis(part, name_to_val)))
    return (out[0], out[1], out[2])


def _eval_axis(part: str, name_to_val: dict[str, float]) -> float:
    """Evaluate one comma-separated axis expression."""
    # Split into signed terms, keeping the leading sign of each.
    terms: list[str] = []
    buf = ""
    for ch in part:
        if ch in "+-" and buf:
            terms.append(buf)
            buf = ch
        else:
            buf += ch
    if buf:
        terms.append(buf)
    if not terms:
        terms = [part]

    total = 0.0
    for term in terms:
        term = term.strip()
        if not term:
            continue
        sign = -1.0 if term.startswith("-") else 1.0
        if term[0] in "+-":
            term = term[1:]
        # Trailing variable letter, if any.
        var = None
        if term and term[-1] in "xyz":
            var = term[-1]
            term = term[:-1]
        if term:  # numeric coefficient or constant
            if "/" in term:
                num, den = term.split("/")
                coef = float(num) / float(den)
            else:
                coef = float(term)
        else:
            coef = 1.0
        if var is not None:
            coef *= name_to_val[var]
        total += sign * coef
    return total


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
    structure.cell = gemmi.UnitCell(
        cell_params["a"], cell_params["b"], cell_params["c"],
        cell_params["alpha"], cell_params["beta"], cell_params["gamma"],
    )
    structure.spacegroup_hm = sg.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    elem_counter: dict[str, int] = {}
    for site in atoms:
        wp = find_wyckoff(sg_number, site.wyckoff)
        rep = eval_coord(wp.coordinates, list(site.free))
        # Validate element. gemmi.Element maps unknown symbols to a
        # placeholder (atomic_number 0) instead of raising, so check the
        # atomic number to reject non-elements like "Xx".
        try:
            element = gemmi.Element(site.element)
        except (RuntimeError, ValueError, TypeError) as e:
            raise ValueError(f"Unknown element {site.element!r}: {e}") from e
        if element.atomic_number == 0:
            raise ValueError(
                f"Unknown element {site.element!r} — not a real element."
            )

        elem_counter[site.element] = elem_counter.get(site.element, 0) + 1
        atom = gemmi.Atom()
        atom.name = f"{site.element}{elem_counter[site.element]}"
        atom.element = element
        atom.occ = float(site.occ)
        frac = gemmi.Fractional(*rep)
        atom.pos = structure.cell.orthogonalize(frac)
        residue.add_atom(atom)

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
        wp = find_wyckoff(sg_number, site.wyckoff)
        counts[site.element] = counts.get(site.element, 0.0) \
            + wp.multiplicity * float(site.occ)
    return counts


def format_formula(counts: dict[str, float]) -> str:
    """Format a composition dict as a formula string.

    Integer counts are normalized by their GCD (so Na4Cl4 -> NaCl); fractional
    counts (from partial occupancy) are shown rounded to 2 decimals, unsorted
    by count, in insertion order.
    """
    # Round to 2 decimals; treat near-integers as integers.
    rounded = {el: round(c, 2) for el, c in counts.items()}
    int_counts = {el: int(round(c)) for el, c in rounded.items()}
    all_integer = all(abs(rounded[el] - int_counts[el]) < 1e-6 for el in rounded)

    if all_integer and int_counts:
        g = _gcd_all(int_counts.values())
        if g > 1:
            int_counts = {el: c // g for el, c in int_counts.items()}

    parts = []
    for el, c in int_counts.items():
        if c == 1:
            parts.append(el)
        elif all_integer:
            parts.append(f"{el}{c}")
        else:
            parts.append(f"{el}{rounded[el]:g}")
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
        if not any(max(abs(f[i] - s[i]) for i in range(3)) < 1e-3
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
        if max(abs(f[i] - b[i]) for i in range(3)) < 1e-3:
            return True
    return False


def _parse_template_axes(template: str):
    """Parse a coordinate template into per-axis (const, var|None, coef).

    ``"1/4+x,3/4-y,1/2"`` -> ``[(0.25,'x',1.0),(0.75,'y',-1.0),(0.5,None,1.0)]``.
    """
    import re
    axes = []
    for ax in template.split(","):
        ax = ax.strip()
        terms = re.findall(r"[+-]?[^+-]+", ax)
        const = 0.0
        var = None
        coef = 1.0
        for t in terms:
            t = t.strip()
            if not t:
                continue
            if t[-1] in "xyz":
                var = t[-1]
                num = t[:-1]
                if num in ("", "+"):
                    coef = 1.0
                elif num == "-":
                    coef = -1.0
                elif "/" in num:
                    n, d = num.split("/")
                    coef = float(n) / float(d)
                else:
                    coef = float(num)
            else:
                if "/" in t:
                    n, d = t.split("/")
                    const += float(n) / float(d)
                else:
                    const += float(t)
        const = const - math.floor(const)
        axes.append((const, var, coef))
    return axes


def detect_wyckoff(sg_number: int, x: float, y: float, z: float
                   ) -> tuple[str, int]:
    """Return (label, multiplicity) of the Wyckoff orbit containing (x, y, z).

    The point may be any representative (canonical or not). For each Wyckoff
    position we ask: does some symmetry image of the point match the position's
    fixed-offset template (with the free variables solved consistently)? The
    first match wins. Returns ("?", multiplicity) if no bundled position matches.
    """
    import gemmi
    sg = gemmi.SpaceGroup(sg_number)
    rep = (_wrap01(x), _wrap01(y), _wrap01(z))
    # Precompute the symmetry images of the point once.
    images = []
    for op in sg.operations():
        c = op.apply_to_xyz([rep[0], rep[1], rep[2]])
        images.append((c[0] - math.floor(c[0]),
                       c[1] - math.floor(c[1]),
                       c[2] - math.floor(c[2])))

    for wp in wyckoff_positions(sg_number):
        axes = _parse_template_axes(wp.coordinates)
        for img in images:
            vals: dict[str, float] = {}
            ok = True
            for (const, var, coef), v in zip(axes, img):
                if var is None:
                    if min(abs(v - const), abs(v - const - 1), abs(v - const + 1)) > 1e-3:
                        ok = False
                        break
                else:
                    fv = (v - const) / coef
                    fv = fv - math.floor(fv)
                    if var in vals and min(abs(vals[var] - fv),
                                           abs(vals[var] - fv - 1),
                                           abs(vals[var] - fv + 1)) > 1e-3:
                        ok = False
                        break
                    vals[var] = fv
            if ok and len(vals) == len({a[1] for a in axes if a[1]}):
                return wp.letter, wp.multiplicity
    return "?", _orbit_size(sg_number, rep)


# Irrational sentinel values for orbit detection — never land on a rational
# special position (mirrors the build script's convention).
_SENTINELS = [0.1414213562, 0.1732050808, 0.2236067977]


def build_structure_frac(sg_number: int, cell_params: dict,
                         atoms: list[FracAtom]) -> gemmi.Structure:
    """Assemble a gemmi.Structure from direct fractional coordinates.

    Each atom is placed at its given (x, y, z); the structure carries the cell
    and space group so the CIF writer emits symmetry operations. Raises
    ValueError on an unknown element.
    """
    sg = gemmi.SpaceGroup(sg_number)
    structure = gemmi.Structure()
    structure.cell = gemmi.UnitCell(
        cell_params["a"], cell_params["b"], cell_params["c"],
        cell_params["alpha"], cell_params["beta"], cell_params["gamma"],
    )
    structure.spacegroup_hm = sg.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    elem_counter: dict[str, int] = {}
    for site in atoms:
        try:
            element = gemmi.Element(site.element)
        except (RuntimeError, ValueError, TypeError) as e:
            raise ValueError(f"Unknown element {site.element!r}: {e}") from e
        if element.atomic_number == 0:
            raise ValueError(
                f"Unknown element {site.element!r} — not a real element."
            )
        elem_counter[site.element] = elem_counter.get(site.element, 0) + 1
        atom = gemmi.Atom()
        atom.name = f"{site.element}{elem_counter[site.element]}"
        atom.element = element
        atom.occ = float(site.occ)
        frac = gemmi.Fractional(_wrap01(site.x), _wrap01(site.y), _wrap01(site.z))
        atom.pos = structure.cell.orthogonalize(frac)
        residue.add_atom(atom)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)
    return structure


def stoichiometry_frac(sg_number: int, atoms: list[FracAtom]) -> dict[str, float]:
    """Per-element cell count = sum(multiplicity * occupancy), multiplicity
    computed by expanding each fractional coordinate under the SG operations."""
    counts: dict[str, float] = {}
    for site in atoms:
        mult = _orbit_size(sg_number, (_wrap01(site.x), _wrap01(site.y), _wrap01(site.z)))
        counts[site.element] = counts.get(site.element, 0.0) + mult * float(site.occ)
    return counts


def validate_atoms_frac(sg_number: int, atoms: list[FracAtom]) -> list[str]:
    """Warn when atoms sharing a site have occupancies summing > 1.0."""
    warnings = []
    occ_by_site: dict[tuple, float] = {}
    for site in atoms:
        key = (round(_wrap01(site.x), 4), round(_wrap01(site.y), 4),
               round(_wrap01(site.z), 4))
        occ_by_site[key] = occ_by_site.get(key, 0.0) + float(site.occ)
    for key, total in occ_by_site.items():
        if total > 1.0 + 1e-3:
            warnings.append(
                f"site {key}: total occupancy {total:g} exceeds 1.0 "
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
    warnings = []
    occ_by_site: dict[tuple, float] = {}
    for site in atoms:
        wp = find_wyckoff(sg_number, site.wyckoff)
        rep = eval_coord(wp.coordinates, list(site.free))
        key = (round(rep[0], 4), round(rep[1], 4), round(rep[2], 4))
        occ_by_site[key] = occ_by_site.get(key, 0.0) + float(site.occ)
    for key, total in occ_by_site.items():
        if total > 1.0 + 1e-3:
            warnings.append(
                f"site {key}: total occupancy {total:g} exceeds 1.0 "
                f"(unphysical for disorder)."
            )
    return warnings
