"""Formula helpers for SHRY post-processing."""

from __future__ import annotations

import re
from collections import OrderedDict
from fractions import Fraction


def parse_formula(formula: str) -> dict[str, int]:
    """Parse a simple formula such as ``Li20Ge2P4S24``."""
    if not formula:
        return {}
    parts = re.findall(r"([A-Z][a-z]?)([0-9]*)", formula)
    if not parts or "".join(el + n for el, n in parts) != formula:
        raise ValueError(f"cannot parse formula {formula!r}")
    out: dict[str, int] = {}
    for el, n in parts:
        out[el] = out.get(el, 0) + int(n or "1")
    return out


def formula_from_counts(counts: dict[str, int | Fraction | float]) -> str:
    """Format counts without reducing by GCD; preserves insertion order."""
    parts = []
    for el, count in counts.items():
        if isinstance(count, Fraction):
            if count.denominator != 1:
                raise ValueError(f"non-integer formula count for {el}: {count}")
            n = count.numerator
        else:
            n = int(round(float(count)))
            if abs(float(count) - n) > 1e-6:
                raise ValueError(f"non-integer formula count for {el}: {count}")
        parts.append(el if n == 1 else f"{el}{n}")
    return "".join(parts)


def counts_from_atom_rows(rows, vacancy_symbol: str = "X") -> OrderedDict[str, int]:
    """Count integer atoms in an ordered atom-row list, skipping vacancies."""
    counts: OrderedDict[str, int] = OrderedDict()
    for row in rows:
        if row.type_symbol == vacancy_symbol:
            continue
        occ = Fraction(str(row.occupancy))
        if occ != 1:
            raise ValueError(
                f"ordered output still has partial occupancy: {row.label} "
                f"{row.type_symbol} occ={row.occupancy}"
            )
        counts[row.type_symbol] = counts.get(row.type_symbol, 0) + 1
    return counts


def assert_formula(rows, target_formula: str | None,
                   vacancy_symbol: str = "X") -> str:
    """Return formula and optionally check it against a target formula."""
    counts = counts_from_atom_rows(rows, vacancy_symbol=vacancy_symbol)
    formula = formula_from_counts(counts)
    if target_formula and parse_formula(formula) != parse_formula(target_formula):
        raise ValueError(f"formula mismatch: got {formula}, expected {target_formula}")
    return formula
