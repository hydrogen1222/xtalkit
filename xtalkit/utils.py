"""Shared utility functions for xtalkit.

The bundled Wyckoff database uses a small, deliberately restricted linear
coordinate language (for example ``x``, ``-x``, ``2x``, ``1/4+x`` and
``1/2-y``).  All runtime consumers must evaluate that language identically;
otherwise a database entry can be valid for ``build`` but crash ``mark`` or
``skeleton``.  The helpers below are the single implementation of that
language.
"""

import math
import re
from fractions import Fraction

import gemmi

_DUMMY_ELEMENTS = ["Xe", "Kr", "Rn", "Ar", "Ne", "He"]

# Generic, mutually incommensurate values used to choose one representative
# from a parameterised Wyckoff position.  Rational values such as 0.3 can put
# x/y/z accidentally on a higher-symmetry sub-position (notably when x == y).
_GENERIC_VARIABLES = {
    "x": math.sqrt(2.0) / 10.0,
    "y": math.sqrt(3.0) / 10.0,
    "z": math.sqrt(5.0) / 10.0,
}


def _number(value: str) -> float:
    """Parse a decimal/integer/fraction used in a Wyckoff expression."""
    try:
        return float(Fraction(value))
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError(f"invalid numeric term {value!r}") from exc


def parse_linear_axis(expression: str) -> tuple[float, dict[str, float]]:
    """Parse one canonical Wyckoff axis into ``constant + coefficients``.

    The accepted grammar is a sum of signed rational constants and linear
    variable terms.  Besides the forms currently present in the bundled
    database, ``x/2`` and ``1/2x`` are supported so regenerated databases with
    half-coefficients remain runtime-compatible.  Products of variables and
    arbitrary Python expressions are intentionally rejected.
    """
    compact = "".join(str(expression).split())
    if not compact:
        raise ValueError("empty Wyckoff coordinate expression")

    constant = 0.0
    coefficients: dict[str, float] = {}
    terms = re.findall(r"[+-]?[^+-]+", compact)
    if not terms or "".join(terms) != compact:
        raise ValueError(f"invalid Wyckoff coordinate expression {expression!r}")

    for signed in terms:
        sign = -1.0 if signed.startswith("-") else 1.0
        term = signed[1:] if signed[:1] in "+-" else signed
        variables = [ch for ch in term if ch in "xyz"]
        if not variables:
            constant += sign * _number(term)
            continue
        if len(variables) != 1:
            raise ValueError(
                f"non-linear or multi-variable term {signed!r} in {expression!r}"
            )

        variable = variables[0]
        before, after = term.split(variable, 1)
        coefficient = _number(before) if before else 1.0
        if after:
            if not after.startswith("/") or len(after) == 1:
                raise ValueError(
                    f"invalid variable term {signed!r} in {expression!r}"
                )
            denominator = _number(after[1:])
            if denominator == 0:
                raise ValueError(
                    f"zero denominator in variable term {signed!r} "
                    f"in {expression!r}"
                )
            coefficient /= denominator
        coefficients[variable] = coefficients.get(variable, 0.0) \
            + sign * coefficient

    coefficients = {
        name: value for name, value in coefficients.items()
        if abs(value) > 1e-15
    }
    return constant, coefficients


def evaluate_linear_axis(expression: str, variables: dict[str, float]) -> float:
    """Evaluate one parsed Wyckoff axis without periodic wrapping."""
    constant, coefficients = parse_linear_axis(expression)
    missing = sorted(set(coefficients) - set(variables))
    if missing:
        raise ValueError(
            f"coordinate {expression!r} requires value(s) for {', '.join(missing)}"
        )
    numeric_variables = {name: float(variables[name]) for name in coefficients}
    if not all(math.isfinite(value) for value in numeric_variables.values()):
        raise ValueError(f"non-finite variable in coordinate {expression!r}")
    result = constant + sum(coeff * numeric_variables[name]
                            for name, coeff in coefficients.items())
    if not math.isfinite(result):
        raise ValueError(f"non-finite result for coordinate {expression!r}")
    return result


def evaluate_coord_template(
    coordinate_template: str,
    variables: dict[str, float],
) -> tuple[float, float, float]:
    """Evaluate a three-axis Wyckoff template and wrap it into ``[0, 1)``."""
    axes = coordinate_template.split(",")
    if len(axes) != 3:
        raise ValueError(
            f"Wyckoff coordinate must have three axes, got {coordinate_template!r}"
        )
    values = [evaluate_linear_axis(axis, variables) % 1.0 for axis in axes]
    return values[0], values[1], values[2]


def representative_coord(coordinate_template: str) -> tuple[float, float, float]:
    """Return a generic (non-special) representative for a Wyckoff template."""
    return evaluate_coord_template(coordinate_template, _GENERIC_VARIABLES)


_CIF_ESD_RE = re.compile(r"(?<=\d)\(\d+\)(?=(?:[eEdD][+-]?\d+)?$)")


def _cif_float(value: str, field: str) -> float:
    """Parse a CIF number, accepting a trailing estimated standard deviation."""
    text = str(value).strip()
    if text in {"", ".", "?"}:
        raise ValueError(f"missing numeric CIF value for {field}")
    text = _CIF_ESD_RE.sub("", text).replace("D", "E").replace("d", "e")
    try:
        result = float(text)
    except ValueError as exc:
        raise ValueError(f"invalid CIF number {value!r} for {field}") from exc
    if not math.isfinite(result):
        raise ValueError(f"non-finite CIF number {value!r} for {field}")
    return result


def _element_from_cif(raw: str | None, label: str) -> gemmi.Element:
    """Resolve a real element without silently substituting a dummy species."""
    candidates = [raw or "", label]
    for candidate in candidates:
        match = re.match(r"[A-Z][a-z]?", str(candidate).strip())
        if not match:
            continue
        try:
            element = gemmi.Element(match.group(0))
        except (RuntimeError, ValueError):
            continue
        if element.atomic_number > 0:
            return element
    raise ValueError(
        f"could not determine a valid element for CIF atom {label!r} "
        f"from type symbol {raw!r}"
    )


def unit_cell_from_parameters(cell_params: dict) -> gemmi.UnitCell:
    """Validate six cell parameters and construct a non-degenerate unit cell."""
    names = ("a", "b", "c", "alpha", "beta", "gamma")
    missing = [name for name in names if name not in cell_params]
    if missing:
        raise ValueError(f"missing cell parameter(s): {', '.join(missing)}")
    values = {name: float(cell_params[name]) for name in names}
    if not all(math.isfinite(value) for value in values.values()):
        raise ValueError("cell parameters must all be finite")
    if any(values[name] <= 0.0 for name in ("a", "b", "c")):
        raise ValueError("cell lengths a, b and c must be positive")
    if any(not 0.0 < values[name] < 180.0
           for name in ("alpha", "beta", "gamma")):
        raise ValueError("cell angles alpha, beta and gamma must lie in (0, 180)")
    cell = gemmi.UnitCell(*(values[name] for name in names))
    if not math.isfinite(cell.volume) or cell.volume <= 1e-10:
        raise ValueError("cell parameters define a degenerate unit cell")
    return cell


def read_cif_structure(path: str, sg_number: int) -> gemmi.Structure:
    """Read a CIF file and build a gemmi.Structure with the given space group.

    Uses gemmi.cif.read_file for robust parsing. Handles both standard CIF
    (_cell_length_a, _atom_site_fract_x) and mmCIF (_cell.length_a,
    _atom_site.Cartn_x) tag conventions. Falls back to gemmi.read_structure
    if neither format is detected.
    """
    doc = gemmi.cif.read_file(path)
    block = doc.sole_block()
    requested_group = gemmi.SpaceGroup(sg_number)

    declared_number = None
    for tag in ("_symmetry_Int_Tables_number", "_space_group_IT_number",
                "_symmetry.Int_Tables_number"):
        raw_number = block.find_value(tag)
        if raw_number and str(raw_number).strip() not in {".", "?"}:
            try:
                declared_number = int(str(raw_number).strip())
            except ValueError as exc:
                raise ValueError(
                    f"invalid CIF space-group number {raw_number!r}"
                ) from exc
            break
    declared_name = None
    for tag in ("_symmetry_space_group_name_H-M", "_space_group_name_H-M_alt",
                "_symmetry.space_group_name_H-M"):
        raw_name = block.find_value(tag)
        if raw_name and str(raw_name).strip() not in {".", "?"}:
            declared_name = str(raw_name).strip(" '\"")
            break
    if declared_name:
        try:
            named_number = gemmi.SpaceGroup(declared_name).number
        except Exception as exc:
            raise ValueError(
                f"invalid CIF space-group name {declared_name!r}"
            ) from exc
        if declared_number is not None and named_number != declared_number:
            raise ValueError(
                "CIF space-group name and number disagree: "
                f"{declared_name!r} is #{named_number}, not #{declared_number}"
            )
        declared_number = named_number
    if declared_number is not None and declared_number != requested_group.number:
        raise ValueError(
            f"requested space group #{requested_group.number} disagrees with "
            f"CIF-declared space group #{declared_number}"
        )

    # Try standard CIF cell tags first, then mmCIF
    def _cell_val(*names: str) -> float | None:
        for name in names:
            v = block.find_value(name)
            if v is not None:
                return _cif_float(v, name)
        return None

    a = _cell_val("_cell_length_a", "_cell.length_a")
    b = _cell_val("_cell_length_b", "_cell.length_b")
    c = _cell_val("_cell_length_c", "_cell.length_c")
    alpha = _cell_val("_cell_angle_alpha", "_cell.angle_alpha")
    beta = _cell_val("_cell_angle_beta", "_cell.angle_beta")
    gamma = _cell_val("_cell_angle_gamma", "_cell.angle_gamma")

    if None in (a, b, c, alpha, beta, gamma):
        raise ValueError(f"Could not read cell parameters from CIF: {path}")

    cell = unit_cell_from_parameters({
        "a": a, "b": b, "c": c,
        "alpha": alpha, "beta": beta, "gamma": gamma,
    })

    # Try to read atoms from standard CIF loop first, then mmCIF
    labels = list(block.find_values("_atom_site_label"))
    if not labels:
        labels = list(block.find_values("_atom_site.label_atom_id"))

    type_col = list(block.find_values("_atom_site_type_symbol"))
    if not type_col:
        type_col = list(block.find_values("_atom_site.type_symbol"))

    occupancies = list(block.find_values("_atom_site_occupancy"))
    if not occupancies:
        occupancies = list(block.find_values("_atom_site.occupancy"))

    fx = list(block.find_values("_atom_site_fract_x"))
    fy = list(block.find_values("_atom_site_fract_y"))
    fz = list(block.find_values("_atom_site_fract_z"))

    # mmCIF uses Cartesian; convert to fractional
    cart_x = list(block.find_values("_atom_site.Cartn_x"))
    cart_y = list(block.find_values("_atom_site.Cartn_y"))
    cart_z = list(block.find_values("_atom_site.Cartn_z"))

    use_cartesian = False
    if not fx:
        fx, fy, fz = cart_x, cart_y, cart_z
        use_cartesian = True

    n_atoms = len(labels)
    if n_atoms == 0:
        raise ValueError("No atom sites found in CIF file")
    if len(set(labels)) != n_atoms:
        raise ValueError("CIF atom-site labels must be unique")
    if not (len(fx) == len(fy) == len(fz) == n_atoms):
        coordinate_kind = "Cartesian" if use_cartesian else "fractional"
        raise ValueError(
            f"incomplete {coordinate_kind} coordinates in CIF: expected "
            f"{n_atoms} rows, got x/y/z={len(fx)}/{len(fy)}/{len(fz)}"
        )

    # Build structure
    structure = gemmi.Structure()
    structure.name = block.name
    structure.cell = cell
    structure.spacegroup_hm = requested_group.hm

    model = gemmi.Model(0)
    chain = gemmi.Chain("A")
    residue = gemmi.Residue()
    residue.name = "UNK"
    residue.seqid = gemmi.SeqId("1")

    for i in range(n_atoms):
        atom = gemmi.Atom()
        atom.name = labels[i] if i < len(labels) else f"X{i}"
        raw_element = type_col[i] if i < len(type_col) else None
        atom.element = _element_from_cif(raw_element, atom.name)
        atom.occ = (
            _cif_float(occupancies[i], f"occupancy for {atom.name}")
            if i < len(occupancies) else 1.0
        )
        if not 0.0 <= atom.occ <= 1.0:
            raise ValueError(
                f"occupancy for {atom.name} must be in [0, 1], got {atom.occ}"
            )

        if use_cartesian:
            atom.pos = gemmi.Position(
                _cif_float(fx[i], f"Cartn_x for {atom.name}"),
                _cif_float(fy[i], f"Cartn_y for {atom.name}"),
                _cif_float(fz[i], f"Cartn_z for {atom.name}"),
            )
        else:
            frac = gemmi.Fractional(
                _cif_float(fx[i], f"fract_x for {atom.name}"),
                _cif_float(fy[i], f"fract_y for {atom.name}"),
                _cif_float(fz[i], f"fract_z for {atom.name}"),
            )
            atom.pos = cell.orthogonalize(frac)

        residue.add_atom(atom)

    chain.add_residue(residue)
    model.add_chain(chain)
    structure.add_model(model)

    return structure


def parse_element_map(value: str) -> dict[str, str]:
    """Parse an element override map like '4a:Xe,16e:Kr' into a dict.

    Each entry must be '<letter>:<element>'. Entries are validated and a
    ValueError with a clear, user-facing message is raised on malformed input
    (missing colon, empty letter/element, etc.) so callers can surface it
    cleanly instead of crashing on tuple-unpacking.
    """
    mapping: dict[str, str] = {}
    for pair in value.split(","):
        pair = pair.strip()
        if not pair:
            continue
        if ":" not in pair:
            raise ValueError(
                f"Invalid element map entry '{pair}' — expected "
                f"'<letter>:<element>' (e.g. '4a:Xe,16e:Kr')."
            )
        letter, elem = pair.split(":", 1)
        letter = letter.strip()
        elem = elem.strip()
        if not letter or not elem:
            raise ValueError(
                f"Invalid element map entry '{pair}' — both Wyckoff letter "
                f"and element must be non-empty (e.g. '4a:Xe')."
            )
        if letter in mapping:
            raise ValueError(f"Duplicate element map entry for {letter}")
        mapping[letter] = elem
    if not mapping:
        raise ValueError(
            "Empty element map — expected '<letter>:<element>' entries "
            "(e.g. '4a:Xe,16e:Kr')."
        )
    return mapping


def assign_dummy_elements(
    wyckoff_letters: list[str],
    element_map: dict[str, str] | None,
) -> dict[str, str]:
    """Assign dummy elements to Wyckoff letters.

    Priority: Xe -> Kr -> Rn -> Ar -> Ne -> He (cycling if needed).

    When element_map is provided, every requested letter must have an entry.
    """
    if element_map is not None:
        extra = sorted(set(element_map) - set(wyckoff_letters))
        if extra:
            raise ValueError(
                f"Element map contains unrequested Wyckoff letter(s): "
                f"{', '.join(extra)}"
            )
        for letter in wyckoff_letters:
            if letter not in element_map:
                raise ValueError(f"No element assigned for Wyckoff letter {letter}")
            raw_element = element_map[letter]
            if not re.fullmatch(r"[A-Z][a-z]?", raw_element):
                raise ValueError(
                    f"Invalid dummy element {raw_element!r} for {letter}"
                )
            try:
                element = gemmi.Element(raw_element)
            except (RuntimeError, ValueError) as exc:
                raise ValueError(
                    f"Invalid dummy element {raw_element!r} for {letter}"
                ) from exc
            if element.atomic_number == 0:
                raise ValueError(
                    f"Invalid dummy element {raw_element!r} for {letter}"
                )
        return dict(element_map)

    assignment = {}
    for i, letter in enumerate(sorted(
        wyckoff_letters,
        key=lambda w: (int("".join(c for c in w if c.isdigit()) or 0), w),
    )):
        assignment[letter] = _DUMMY_ELEMENTS[i % len(_DUMMY_ELEMENTS)]
    return assignment


def parse_coord(s: str, variable_default: float = 0.3) -> float:
    """Evaluate one coordinate expression using one common variable value.

    Kept as a small public compatibility helper.  New code that evaluates a
    complete three-axis Wyckoff template should use :func:`representative_coord`
    so x/y/z receive distinct generic values and cannot accidentally land on a
    special sub-position.
    """
    variables = {name: float(variable_default) for name in "xyz"}
    return evaluate_linear_axis(s, variables)
