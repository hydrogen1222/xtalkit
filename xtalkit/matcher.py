"""Match atoms in a structure to their Wyckoff positions."""

import math

import gemmi
from xtalkit.spacegroup import WyckoffInfo
from xtalkit.utils import parse_linear_axis, representative_coord


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
    base = gemmi.Fractional(*representative_coord(formula_coords))

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


def _axis_distance(a: float, b: float) -> float:
    """Absolute minimum-image separation of two fractional axis values."""
    return abs((a - b) - round(a - b))


def _variable_candidates(
    observed: float,
    constant: float,
    coefficient: float,
    tolerance: float,
) -> list[float]:
    """Solve ``constant + coefficient*x == observed (mod 1)`` for x in [0,1).

    Coefficients such as 2 have several periodic branches.  Explicitly
    retaining those branches avoids the false negatives caused by simply
    dividing and taking modulo one.
    """
    if abs(coefficient) < 1e-15:
        return []
    candidates: list[float] = []
    delta = observed - constant
    # Derive every periodic image that can put the variable in [0, 1). This is
    # exact for arbitrary finite linear coefficients and avoids a hidden fixed
    # image bound for future database expressions.
    lower = min(0.0, coefficient) - delta
    upper = max(0.0, coefficient) - delta
    first_image = math.ceil(lower - 1e-12)
    last_image = math.floor(upper + 1e-12)
    for image in range(first_image, last_image + 1):
        raw = (delta + image) / coefficient
        if raw < -1e-9 or raw >= 1.0 + 1e-9:
            continue
        candidate = raw % 1.0
        predicted = constant + coefficient * candidate
        if _axis_distance(predicted, observed) <= tolerance and not any(
            _axis_distance(candidate, old) <= 1e-10 for old in candidates
        ):
            candidates.append(candidate)
    return candidates


def _template_matches(
    point: tuple[float, float, float],
    coordinate_template: str,
    tolerance: float,
) -> bool:
    """Return whether one point satisfies a canonical Wyckoff template."""
    axes = coordinate_template.split(",")
    if len(axes) != 3:
        raise ValueError(
            f"Wyckoff coordinate must have three axes, got {coordinate_template!r}"
        )

    equations: dict[str, list[tuple[float, float, float]]] = {}
    for expression, observed in zip(axes, point):
        constant, coefficients = parse_linear_axis(expression)
        if len(coefficients) > 1:
            raise ValueError(
                "runtime Wyckoff matcher requires at most one free variable "
                f"per axis, got {expression!r}"
            )
        if not coefficients:
            if _axis_distance(observed, constant) > tolerance:
                return False
            continue
        variable, coefficient = next(iter(coefficients.items()))
        equations.setdefault(variable, []).append(
            (observed, constant, coefficient)
        )

    for variable_equations in equations.values():
        observed, constant, coefficient = variable_equations[0]
        candidates = _variable_candidates(
            observed, constant, coefficient, tolerance
        )
        if not any(
            all(
                _axis_distance(c + k * candidate, value) <= tolerance
                for value, c, k in variable_equations
            )
            for candidate in candidates
        ):
            return False
    return True


def detect_wyckoff_position(
    point: tuple[float, float, float],
    sg: gemmi.SpaceGroup,
    wyckoffs: list[WyckoffInfo],
    tolerance: float = 1e-3,
) -> WyckoffInfo | None:
    """Identify the most-special Wyckoff orbit containing ``point``.

    Positions must be supplied in increasing multiplicity order (as returned
    by :func:`xtalkit.spacegroup.wyckoff_positions`).  Testing every symmetry
    image permits non-canonical asymmetric-unit representatives while solving
    the template's free variables exactly up to ``tolerance``.
    """
    if not math.isfinite(tolerance) or tolerance < 0:
        raise ValueError("matching tolerance must be finite and non-negative")
    numeric_point = tuple(float(value) for value in point)
    if not all(math.isfinite(value) for value in numeric_point):
        raise ValueError("fractional coordinates must be finite")
    wrapped = tuple(value % 1.0 for value in numeric_point)
    images: list[tuple[float, float, float]] = []
    for op in sg.operations():
        transformed = op.apply_to_xyz(list(wrapped))
        image = tuple(float(value) % 1.0 for value in transformed)
        if not any(max(_axis_distance(a, b) for a, b in zip(image, old)) < 1e-10
                   for old in images):
            images.append(image)

    for wyckoff in wyckoffs:
        if any(_template_matches(image, wyckoff.coordinates, tolerance)
               for image in images):
            return wyckoff
    return None


def match_atoms(
    structure: gemmi.Structure,
    wyckoffs: list[WyckoffInfo],
    tolerance: float = 0.01,
) -> dict[str, str]:
    """Match each atom in the structure to its nearest Wyckoff position.

    For each atom, compute the distance to every symmetry-equivalent
    Wyckoff site. If the minimum distance is within tolerance, the atom
    is considered to occupy that Wyckoff position.

    Args:
        structure: gemmi.Structure with atoms (Cartesian positions stored
                   in atom.pos) and a valid cell for fractional conversion.
        wyckoffs: List of Wyckoff positions for the space group.
        tolerance: Maximum residual in the Wyckoff coordinate constraints,
                   in fractional units. Default 0.01.

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

    mapping: dict[str, str] = {}

    for model in structure:
        for chain in model:
            for residue in chain:
                for atom in residue:
                    atom_frac = cell.fractionalize(atom.pos)
                    match = detect_wyckoff_position(
                        (atom_frac.x, atom_frac.y, atom_frac.z),
                        sg,
                        wyckoffs,
                        tolerance=tolerance,
                    )
                    if match is not None:
                        mapping[atom.name] = match.letter

    return mapping
