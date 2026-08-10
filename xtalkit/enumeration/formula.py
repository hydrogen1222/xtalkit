"""Formula helpers for SHRY post-processing."""

from __future__ import annotations

import re
from collections import OrderedDict
from fractions import Fraction

import gemmi


_ORDERED_TYPE_RE = re.compile(
    r"(?P<element>[A-Z][a-z]?)(?:(?:\d+(?:\.\d+)?)?[+-])?"
)


def _validate_element_symbol(symbol: str) -> str:
    """Return a canonical real-element symbol or fail closed."""
    if not isinstance(symbol, str) or not re.fullmatch(r"[A-Z][a-z]?", symbol):
        raise ValueError(f"invalid chemical element symbol {symbol!r}")
    try:
        element = gemmi.Element(symbol)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"unknown chemical element {symbol!r}") from exc
    if element.atomic_number <= 0:
        raise ValueError(f"unknown chemical element {symbol!r}")
    return element.name


def real_element_symbol(type_symbol: str) -> str:
    """Normalize an ordered CIF type symbol, including standard charge suffixes.

    Examples accepted here include ``Li``, ``Li+``, ``Li1+``, and ``S2-``.
    Composite/disordered tokens and arbitrary suffixes are deliberately refused:
    ordered outputs must identify one unambiguous real element per atom row.
    """
    text = str(type_symbol).strip()
    match = _ORDERED_TYPE_RE.fullmatch(text)
    if match is None:
        raise ValueError(f"invalid ordered atom type symbol {type_symbol!r}")
    return _validate_element_symbol(match.group("element"))


def parse_formula(formula: str) -> dict[str, int]:
    """Parse a simple formula such as ``Li20Ge2P4S24``."""
    if not formula:
        raise ValueError("formula must not be empty")
    parts = re.findall(r"([A-Z][a-z]?)([0-9]*)", formula)
    if not parts or "".join(el + n for el, n in parts) != formula:
        raise ValueError(f"cannot parse formula {formula!r}")
    out: dict[str, int] = {}
    for el, n in parts:
        el = _validate_element_symbol(el)
        count = int(n or "1")
        if count <= 0:
            raise ValueError(f"formula count for {el} must be positive, got {count}")
        out[el] = out.get(el, 0) + count
    return out


def formula_from_counts(counts: dict[str, int | Fraction | float]) -> str:
    """Format counts without reducing by GCD; preserves insertion order."""
    parts = []
    for el, count in counts.items():
        el = _validate_element_symbol(el)
        if isinstance(count, bool):
            raise ValueError(f"invalid formula count for {el}: {count!r}")
        if isinstance(count, Fraction):
            if count.denominator != 1:
                raise ValueError(f"non-integer formula count for {el}: {count}")
            n = count.numerator
        else:
            n = int(round(float(count)))
            if abs(float(count) - n) > 1e-6:
                raise ValueError(f"non-integer formula count for {el}: {count}")
        if n <= 0:
            raise ValueError(f"formula count for {el} must be positive, got {n}")
        parts.append(el if n == 1 else f"{el}{n}")
    if not parts:
        raise ValueError("cannot format an empty chemical formula")
    return "".join(parts)


def counts_from_atom_rows(rows, vacancy_symbol: str = "X") -> OrderedDict[str, int]:
    """Count integer atoms in an ordered atom-row list, skipping vacancies.

    Species labels are normalized to bare element symbols so that
    ``Li1+``/``Li+``/``Li`` (refinement vs pymatgen notation) all count as Li.
    """
    counts: OrderedDict[str, int] = OrderedDict()
    vac = element_symbol(vacancy_symbol)
    for row in rows:
        el = element_symbol(row.type_symbol)
        if el == vac:
            continue
        el = real_element_symbol(row.type_symbol)
        occ = Fraction(str(row.occupancy))
        if occ != 1:
            raise ValueError(
                f"ordered output still has partial occupancy: {row.label} "
                f"{row.type_symbol} occ={row.occupancy}"
            )
        counts[el] = counts.get(el, 0) + 1
    return counts


def element_symbol(type_symbol: str) -> str:
    """Extract the leading element symbol from a CIF type_symbol.

    Handles ``Li``, ``Li1+``, ``Li+``, ``S2-``, ``Cl1-``, ``X`` (vacancy),
    and pymatgen composite tokens like ``Li0.5X0.5`` (takes the first element).
    """
    s = str(type_symbol).strip()
    m = re.match(r"[A-Z][a-z]?", s)
    return m.group(0) if m else s


def validate_vacancy_symbol(symbol: str) -> str:
    """Validate a pymatgen-compatible non-element vacancy pseudo-species."""
    if not isinstance(symbol, str) or not re.fullmatch(r"[A-Z][a-z]?", symbol):
        raise ValueError(
            "vacancy symbol must be a one- or two-letter pseudo-species such "
            "as X, Xa, or Q"
        )
    try:
        from pymatgen.core import DummySpecies
    except ImportError as exc:  # pragma: no cover - optional dependency
        raise RuntimeError(
            "pymatgen is required to validate vacancy species; "
            "run 'uv sync --extra enumerate'"
        ) from exc
    try:
        dummy = DummySpecies(symbol)
    except (TypeError, ValueError) as exc:
        raise ValueError(
            f"invalid vacancy pseudo-species {symbol!r}; it must not collide "
            "with a real chemical element"
        ) from exc
    return dummy.symbol


def assert_formula(rows, target_formula: str | None,
                   vacancy_symbol: str = "X") -> str:
    """Return formula and optionally check it against a target formula."""
    counts = counts_from_atom_rows(rows, vacancy_symbol=vacancy_symbol)
    formula = formula_from_counts(counts)
    if target_formula is not None and parse_formula(formula) != parse_formula(target_formula):
        raise ValueError(f"formula mismatch: got {formula}, expected {target_formula}")
    return formula
