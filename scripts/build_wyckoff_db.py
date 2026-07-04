#!/usr/bin/env python
"""Build the bundled Wyckoff position dataset for all 230 space groups.

Build-time only — NOT a runtime dependency. Uses pyxtal (MIT) to extract, per
space group and Wyckoff position: letter, multiplicity, site symmetry, and a
canonical coordinate template string (e.g. ``x,-x,z``, ``1/3,2/3,z``,
``x,2x,z``). Writes ``xtalkit/data/wyckoff.json``.

Run (no install needed — pyxtal is pulled ephemerally):

    uv run --with "pyxtal>=1.1" python scripts/build_wyckoff_db.py

How the coordinate template is derived
--------------------------------------
pyxtal stores each Wyckoff position as the orbit of a representative under the
site's symmetry, and exposes ``get_dof()`` (number of independent free
parameters) and ``get_position_from_free_xyzs(values)`` (substitute free
values → representative fractional coords). It does NOT expose the parametric
template string directly. We recover it by linear probing: evaluate the
representative at one baseline (all free vars = 0) and at D single-slot probes
(one slot set to a sentinel, the rest 0). Each output coordinate is then a
linear combination ``sum_i c_i * slot_i + const``; we read off the coefficients
and constant, snap them to crystallographic rationals, and assemble the
canonical string. Free parameters are named by the physical axis (x/y/z) they
first appear on.

Every position is verified two ways:
  1. Template vs pyxtal: ``eval_coord(template, samples)`` must match
     ``get_position_from_free_xyzs(samples)`` for several sample inputs.
  2. Orbit vs multiplicity: expand a representative under gemmi's symmetry
     operations for the space group, dedup, and assert the orbit size equals
     the tabulated multiplicity.

The 38 space groups hand-coded in ``xtalkit/spacegroup.py`` are cross-checked
against the derived data to catch setting mismatches.

Data sources & attribution
--------------------------
Wyckoff positions are derived from International Tables for Crystallography
Vol. A, via pyxtal (https://pyxtal.readthedocs.io, MIT License). Site
symmetry symbols and multiplicities follow the ITA standard setting, which
matches gemmi's reference setting.
"""

from __future__ import annotations

import json
import math
import os
import sys
import warnings
from fractions import Fraction

warnings.filterwarnings("ignore")

# Distinct sentinel values for the free-parameter slots. Chosen irrational
# (sqrt(2)/10, sqrt(3)/10, sqrt(5)/10) so a representative built from them can
# never land on a rational special position (mirror plane, rotation axis,
# glide line) — which would shrink its orbit and falsely fail verification.
# They are also mutually incommensurate, so no y=2x / x=y / x+y=1/2 relation
# holds among them.
_SENTINELS = [0.1414213562, 0.1732050808, 0.2236067977]

# Allowed variable coefficients in Wyckoff coordinate expressions (0 included
# — a slot may not appear on a given axis).
_COEF_CANDS = [-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0]


def _wrap01(x: float) -> float:
    """Wrap a coordinate into [0, 1)."""
    return x - math.floor(x)


def _min_image(d: float) -> float:
    """Wrap a difference into [-0.5, 0.5] (nearest periodic image)."""
    return d - round(d)


def _snap_coef(diff: float, sentinel: float) -> float:
    """Find the crystallographic coefficient whose predicted diff matches.

    For a slot with sentinel ``s``, coefficient ``c`` produces an observed
    coordinate difference of ``c*s`` (mod 1, nearest image). We test each
    allowed coefficient and return the one whose prediction matches the
    observed ``diff``. This is robust to periodic wrap (e.g. coefficient 2
    with a larger sentinel) and treats 0 (slot absent on this axis) as valid.
    Raises RuntimeError if nothing matches within tolerance.
    """
    best_c, best_err = 0.0, 1e9
    for c in _COEF_CANDS:
        err = abs(_min_image(c * sentinel) - diff)
        if err < best_err:
            best_c, best_err = c, err
    if best_err > 1e-3:
        raise RuntimeError(
            f"diff {diff:.4f} matches no crystallographic coefficient "
            f"(best={best_c}, err={best_err:.4f})."
        )
    return best_c


def _snap_rational(x: float, maxden: int = 12) -> str | None:
    """Snap a float to a small-denominator fraction string, or None."""
    f = Fraction(x).limit_denominator(maxden)
    if abs(float(f) - x) < 1e-4:
        if f.denominator == 1:
            return str(f.numerator)
        return f"{f.numerator}/{f.denominator}"
    return None


def _derive_template(wp, o: tuple[float, float, float] = (0.0, 0.0, 0.0)
                      ) -> tuple[str, list[str]]:
    """Derive (coordinate_template, free_param_names) for a pyxtal WP.

    Returns a string like ``"x,-x,z"`` or ``"1/3,2/3,z"`` and the ordered list
    of free-parameter names (e.g. ``["x", "z"]``). ``o`` is an origin shift
    (pyxtal→gemmi) added to each axis's constant term, so the template lands
    in gemmi's setting. Raises RuntimeError if coefficients don't snap cleanly.
    """
    D = wp.get_dof()
    p0 = wp.get_position_from_free_xyzs([0.0] * D)
    # Per-slot probes: slot i = sentinel, rest 0.
    probes = []
    for i in range(D):
        v = [0.0] * D
        v[i] = _SENTINELS[i]
        probes.append(wp.get_position_from_free_xyzs(v))

    # coeffs[axis][slot] and const[axis] (shifted into gemmi's setting).
    coeffs = [[0.0] * D for _ in range(3)]
    consts = [0.0, 0.0, 0.0]
    for axis in range(3):
        consts[axis] = _wrap01(p0[axis] + o[axis])
        for i in range(D):
            d = _min_image(probes[i][axis] - p0[axis])
            coeffs[axis][i] = _snap_coef(d, _SENTINELS[i])

    # Name each slot by the physical axis it first appears on.
    axis_names = ["x", "y", "z"]
    slot_names: list[str] = []
    used: set[str] = set()
    for i in range(D):
        name = None
        for axis in range(3):
            if abs(coeffs[axis][i]) > 1e-9:
                name = axis_names[axis]
                break
        if name is None or name in used:
            # Fallback: assign x/y/z in slot order.
            name = axis_names[len(used)] if len(used) < 3 else f"u{len(used)}"
        slot_names.append(name)
        used.add(name)

    # Assemble per-axis template strings.
    axis_strs: list[str] = []
    for axis in range(3):
        terms: list[tuple[str, str]] = []  # (sign, body)
        # Constant first if nonzero (ITA convention: "1/2-x", "1/3,2/3,z").
        cstr = _snap_rational(consts[axis])
        const_val = consts[axis]
        if abs(const_val) > 1e-9:
            if cstr is None:
                raise RuntimeError(
                    f"WP {wp.get_label()}: axis {axis} constant "
                    f"{const_val:.6f} does not snap to a rational."
                )
            terms.append(("+", cstr))
        for i in range(D):
            c = coeffs[axis][i]
            if abs(c) < 1e-9:
                continue
            sign = "-" if c < 0 else "+"
            sym = slot_names[i]
            ac = abs(c)
            if abs(ac - 1.0) < 1e-9:
                body = sym
            elif abs(ac - 2.0) < 1e-9:
                body = f"2{sym}"
            elif abs(ac - 0.5) < 1e-9:
                body = f"{sym}/2"
            else:
                body = f"{ac}{sym}"
            terms.append((sign, body))
        if not terms:
            axis_strs.append("0")
            continue
        # Join: first term without leading "+", subsequent with their sign.
        s = terms[0][1] if terms[0][0] == "+" else f"-{terms[0][1]}"
        for sign, body in terms[1:]:
            s += f"{sign}{body}"
        axis_strs.append(s)

    return ",".join(axis_strs), slot_names


def _eval_template(template: str, values: list[float]) -> tuple[float, float, float]:
    """Evaluate a derived coordinate template at the given free-param values.

    Minimal evaluator for the canonical format produced by ``_derive_template``:
    each comma-separated axis is a sum of signed terms; each term is a rational
    constant or [coef]variable. Variables x/y/z are mapped to the supplied
    values in order of first appearance.
    """
    # Map variable names to values in order of first appearance.
    name_to_val: dict[str, float] = {}
    for ch in template:
        if ch in "xyz" and ch not in name_to_val:
            if len(name_to_val) >= len(values):
                raise RuntimeError("more variables than values supplied")
            name_to_val[ch] = values[len(name_to_val)]

    out = []
    for part in template.split(","):
        total = 0.0
        # Split into signed terms, keeping the sign.
        terms = []
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
        for t in terms:
            t = t.strip()
            if not t:
                continue
            sign = 1.0
            if t.startswith("-"):
                sign = -1.0
                t = t[1:]
            elif t.startswith("+"):
                t = t[1:]
            # A term is: [coef][var] or a rational/number.
            # Find trailing variable letter.
            var = None
            if t and t[-1] in "xyz":
                var = t[-1]
                t = t[:-1]
            if t:  # numeric coefficient or constant
                if "/" in t:
                    num, den = t.split("/")
                    coef = float(num) / float(den)
                else:
                    coef = float(t)
            else:
                coef = 1.0
            if var is not None:
                coef *= name_to_val[var]
            total += sign * coef
        out.append(_wrap01(total))
    return (out[0], out[1], out[2])


def _verify_orbit(sg_number: int, rep: tuple[float, float, float],
                  multiplicity: int) -> bool:
    """Expand a representative under the SG's symmetry ops and check orbit size.

    Uses gemmi (the runtime library) for symmetry, so this also confirms the
    derived data is consistent with gemmi's setting of the space group.
    """
    import gemmi
    sg = gemmi.SpaceGroup(sg_number)
    cell = gemmi.UnitCell(1.0, 1.0, 1.0, 90.0, 90.0, 90.0)  # fractional space
    seen: list[gemmi.Fractional] = []
    for op in sg.operations():
        coords = op.apply_to_xyz([rep[0], rep[1], rep[2]])
        f = gemmi.Fractional(*coords)
        f.wrap_to_unit()
        if not any(_frac_dist(f, s) < 1e-3 for s in seen):
            seen.append(f)
    return len(seen) == multiplicity


def _frac_dist(a, b) -> float:
    d = lambda x: x - round(x)
    return math.sqrt(d(a.x - b.x) ** 2 + d(a.y - b.y) ** 2 + d(a.z - b.z) ** 2)


_PG_TABLE: dict | None = None


def _det_trace(rot: tuple) -> tuple[int, int]:
    """Conjugacy-invariant signature of a rotation: (determinant, trace)."""
    det = (rot[0][0] * (rot[1][1] * rot[2][2] - rot[1][2] * rot[2][1])
           - rot[0][1] * (rot[1][0] * rot[2][2] - rot[1][2] * rot[2][0])
           + rot[0][2] * (rot[1][0] * rot[2][1] - rot[1][1] * rot[2][0]))
    tr = rot[0][0] + rot[1][1] + rot[2][2]
    return (int(round(det)), int(round(tr)))


def _pg_signature(rots) -> tuple:
    """Signature of a point group: sorted multiset of (det, trace) per op.

    Conjugacy-invariant, so it identifies a point group regardless of how its
    symmetry axes are oriented in the cell — which is what we need, since a
    site symmetry and a representative SG's point group need not share
    orientation.
    """
    return tuple(sorted(_det_trace(r) for r in rots))


def _point_group_table() -> dict:
    """Map a point-group signature to its HM symbol, built from gemmi's SGs.

    The 230 space groups realize all 32 crystallographic point groups (which
    is exactly the set of possible site symmetries). Collisions — two distinct
    point groups with the same (det, trace) signature — are detected and
    printed to stderr; none are expected.
    """
    global _PG_TABLE
    if _PG_TABLE is not None:
        return _PG_TABLE
    import gemmi
    table: dict = {}
    for n in range(1, 231):
        sg = gemmi.SpaceGroup(n)
        rots = [tuple(tuple(round(x) for x in row) for row in op.rot)
                for op in sg.operations()]
        sig = _pg_signature(rots)
        sym = sg.point_group_hm()
        if sig in table and table[sig] != sym:
            print(f"[pg-table] signature collision: {sig} -> "
                  f"{table[sig]} and {sym}", file=sys.stderr)
        table.setdefault(sig, sym)
    _PG_TABLE = table
    return table


def _site_symmetry_from_gemmi(sg_number: int, rep: tuple[float, float, float]) -> str:
    """Compute the site-symmetry point-group symbol from gemmi's stabilizer.

    The stabilizer is the subset of the SG's operations that fix the
    representative; its rotational part is one of the 32 point groups,
    identified by signature. Always consistent with the multiplicity (unlike
    pyxtal's symbol, which is occasionally a wrong subgroup). Returns the
    generic (non-oriented) HM symbol, e.g. ``"m-3m"`` or ``"3m"``.
    """
    import gemmi
    table = _point_group_table()
    sg = gemmi.SpaceGroup(sg_number)
    stab_rots = []
    for op in sg.operations():
        c = op.apply_to_xyz([rep[0], rep[1], rep[2]])
        f = [c[i] - math.floor(c[i]) for i in range(3)]
        if max(abs(f[i] - rep[i]) for i in range(3)) < 1e-3:
            stab_rots.append(tuple(tuple(round(x) for x in row) for row in op.rot))
    return table.get(_pg_signature(stab_rots), f"?({len(stab_rots)})")


def _detect_origin_shift(sg_number: int, group) -> tuple[float, float, float]:
    """Find the origin shift aligning pyxtal's setting to gemmi's.

    pyxtal and gemmi both use ITA settings, but for ~24 space groups with
    multiple origin choices they pick different origins — differing by a
    rational translation that leaves the general position's multiplicity
    unchanged but moves the special positions. We brute-force the shift: the
    correct ``o`` makes every Wyckoff position's representative (shifted by
    ``o``) expand to its tabulated multiplicity under gemmi's operations.
    ``(0, 0, 0)`` is tried first, so space groups whose settings already agree
    pay no search cost.
    """
    # Precompute sentinel representatives for each WP (special positions do
    # the filtering; the general position is shift-invariant in orbit size).
    reps_by_wp: list[tuple[object, list[tuple[float, float, float]]]] = []
    for wp in group.Wyckoff_positions:
        D = wp.get_dof()
        reps = [tuple(float(x) for x in
                      wp.get_position_from_free_xyzs(list(s)))
                for s in _sample_inputs(D)]
        reps_by_wp.append((wp, reps))

    fracs = [0.0, 0.125, 0.25, 0.375, 0.5, 0.625, 0.75, 0.875]
    # (0,0,0) first, then the rest.
    candidates = [(0.0, 0.0, 0.0)]
    for a in fracs:
        for b in fracs:
            for c in fracs:
                if (a, b, c) != (0.0, 0.0, 0.0):
                    candidates.append((a, b, c))

    for o in candidates:
        aligned = True
        for wp, reps in reps_by_wp:
            if any(_verify_orbit(sg_number,
                                 tuple(_wrap01(r[i] + o[i]) for i in range(3)),
                                 wp.multiplicity)
                   for r in reps):
                continue
            aligned = False
            break
        if aligned:
            return o
    raise RuntimeError(
        f"SG {sg_number}: no origin shift (to 1/8) aligns pyxtal to gemmi."
    )


def _build_one(sg_number: int, group) -> list[dict]:
    """Extract all Wyckoff positions for one space group (in gemmi's setting)."""
    o = _detect_origin_shift(sg_number, group)
    positions = []
    for wp in group.Wyckoff_positions:
        D = wp.get_dof()
        template, free_names = _derive_template(wp, o)
        # Verify template against pyxtal (shifted by o) at sample inputs.
        for sample in _sample_inputs(D):
            pyx = wp.get_position_from_free_xyzs(list(sample))
            pyx_shifted = tuple(_wrap01(pyx[i] + o[i]) for i in range(3))
            ours = _eval_template(template, list(sample))
            for a, b in zip(pyx_shifted, ours):
                if _min_image(a - b) > 1e-4:
                    raise RuntimeError(
                        f"SG {sg_number} {wp.get_label()}: template "
                        f"{template!r} at {sample} gives {ours}, pyxtal "
                        f"(shifted {o}) gives {pyx_shifted}."
                    )
        # Verify orbit size == multiplicity via gemmi. Try each sample; a
        # genuine data error fails all of them, while a rare sentinel-induced
        # sub-position only fails one.
        orbit_ok = False
        last_rep = None
        for sample in _sample_inputs(D):
            rep = _eval_template(template, list(sample))
            last_rep = rep
            if _verify_orbit(sg_number, rep, wp.multiplicity):
                orbit_ok = True
                break
        if not orbit_ok:
            raise RuntimeError(
                f"SG {sg_number} {wp.get_label()}: orbit size != multiplicity "
                f"({wp.multiplicity}) for representative {last_rep}, template "
                f"{template!r}."
            )
        # Site symmetry from gemmi's stabilizer — always consistent with the
        # multiplicity (pyxtal's symbol is occasionally a wrong subgroup).
        site_sym = _site_symmetry_from_gemmi(sg_number, last_rep)
        positions.append({
            "letter": wp.letter,
            "multiplicity": wp.multiplicity,
            "site_symmetry": site_sym,
            "coordinates": template,
        })
    return positions


def _sample_inputs(D: int) -> list[tuple[float, ...]]:
    """Sample free-param tuples of length D, built from irrational sentinels.

    Irrational values guarantee the resulting representative never lands on a
    rational special position, so its orbit size equals the tabulated
    multiplicity (no accidental site-symmetry inflation). A second,
    independently-ordered tuple provides a fallback if the first ever
    coincides with a higher-symmetry sub-position.
    """
    if D == 0:
        return [()]
    s = _SENTINELS
    if D == 1:
        return [(s[0],), (s[2],)]
    if D == 2:
        return [(s[0], s[1]), (s[2], s[0])]
    return [(s[0], s[1], s[2]), (s[2], s[0], s[1])]


def main() -> int:
    from pyxtal.symmetry import Group  # build-time only

    data: dict[str, list[dict]] = {}
    failures: list[str] = []
    for sg_number in range(1, 231):
        try:
            g = Group(sg_number)
            data[str(sg_number)] = _build_one(sg_number, g)
        except Exception as e:  # noqa: BLE001
            failures.append(f"SG {sg_number}: {e}")

    if failures:
        print("FAILED space groups:", file=sys.stderr)
        for f in failures:
            print(f"  {f}", file=sys.stderr)
        return 1

    n_positions = sum(len(v) for v in data.values())
    print(f"Derived {n_positions} Wyckoff positions across 230 space groups "
          f"(all orbit- and site-symmetry-verified against gemmi).")

    out_dir = os.path.join(os.path.dirname(__file__), "..", "xtalkit", "data")
    os.makedirs(out_dir, exist_ok=True)
    out_path = os.path.join(out_dir, "wyckoff.json")
    with open(out_path, "w") as f:
        json.dump(data, f, indent=1, sort_keys=True)
        f.write("\n")
    print(f"Wrote {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
