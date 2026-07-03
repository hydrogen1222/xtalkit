"""Match atoms in a structure to their Wyckoff positions."""

import math

import gemmi
from xtalkit.spacegroup import WyckoffInfo


def _frac_dist(a: gemmi.Fractional, b: gemmi.Fractional) -> float:
    """Minimum-image distance in fractional coordinate space.

    Each component is wrapped into [-0.5, 0.5] so that an atom just inside
    one cell face (e.g. 0.99) is treated as adjacent to a Wyckoff site on
    the opposite face (e.g. 0.0), which is the same physical point under
    periodic boundary conditions. Without this, boundary atoms get
    mis-assigned to the wrong Wyckoff position.
    """
    def _wrap(d: float) -> float:
        return d - round(d)  # nearest image, in [-0.5, 0.5]
    return math.sqrt(
        _wrap(a.x - b.x) ** 2 + _wrap(a.y - b.y) ** 2 + _wrap(a.z - b.z) ** 2
    )


def _apply_symmetry(
    formula_coords: str, sg: gemmi.SpaceGroup
) -> list[gemmi.Fractional]:
    """Generate all symmetry-equivalent positions for a fractional coordinate string.

    The coordinate string uses gemmi format, e.g. 'x,x,z' or '0,0,0'.
    Special positions with literal fractions (0, 1/2, etc.) are resolved
    by evaluating the formula at a representative x,y,z then applying
    symmetry operations.
    """
    parts = formula_coords.split(",")

    try:
        # Try literal numeric/fraction parsing first, e.g. "0,0,0" or "1/2,1/2,1/2"

        def _parse_part(s: str) -> float:
            s = s.strip()
            if "/" in s:
                num, den = s.split("/")
                return float(num) / float(den)
            return float(s)

        base = gemmi.Fractional(*[_parse_part(p) for p in parts])
    except (ValueError, ZeroDivisionError):
        # Parameterized coordinate like 'x,x,z' -- use 0.3 as representative
        rep = (
            formula_coords.replace("x", "0.3")
            .replace("y", "0.3")
            .replace("z", "0.3")
        )
        parts = rep.split(",")
        base = gemmi.Fractional(*[_parse_part(p) for p in parts])

    # Apply all symmetry operations
    equivalents: list[gemmi.Fractional] = []
    for op in sg.operations():
        coords = op.apply_to_xyz([base.x, base.y, base.z])
        equiv = gemmi.Fractional(*coords)
        equiv.wrap_to_unit()
        equivalents.append(equiv)

    # Deduplicate within tolerance
    unique: list[gemmi.Fractional] = []
    for eq in equivalents:
        if not any(_frac_dist(eq, ue) < 0.001 for ue in unique):
            unique.append(eq)
    return unique


def match_atoms(
    structure: gemmi.Structure,
    wyckoffs: list[WyckoffInfo],
    tolerance: float = 0.5,
) -> dict[str, str]:
    """Match each atom in the structure to its nearest Wyckoff position.

    For each atom, compute the distance to every symmetry-equivalent
    Wyckoff site. If the minimum distance is within tolerance, the atom
    is considered to occupy that Wyckoff position.

    Args:
        structure: gemmi.Structure with atoms (Cartesian positions stored
                   in atom.pos) and a valid cell for fractional conversion.
        wyckoffs: List of Wyckoff positions for the space group.
        tolerance: Maximum distance (in fractional coordinate units)
                   for a match. Default 0.5.

    Returns:
        dict mapping atom label (str) to Wyckoff letter (str).
        Atoms with no match within tolerance are excluded.
    """
    sg_str = structure.spacegroup_hm
    if sg_str:
        sg = gemmi.SpaceGroup(sg_str)
    else:
        # Fallback: use P1 (identity only) if no space group set
        sg = gemmi.SpaceGroup(1)

    cell = structure.cell
    if cell is None:
        return {}

    # Pre-compute all Wyckoff equivalent sites
    wyckoff_sites: dict[str, list[gemmi.Fractional]] = {}
    for w in wyckoffs:
        wyckoff_sites[w.letter] = _apply_symmetry(w.coordinates, sg)

    mapping: dict[str, str] = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_frac = cell.fractionalize(atom.pos)
                    best_letter: str | None = None
                    best_dist = float("inf")

                    for letter, sites in wyckoff_sites.items():
                        for site in sites:
                            d = _frac_dist(atom_frac, site)
                            if d < best_dist:
                                best_dist = d
                                best_letter = letter

                    if best_dist <= tolerance and best_letter is not None:
                        mapping[atom.name] = best_letter

    return mapping
