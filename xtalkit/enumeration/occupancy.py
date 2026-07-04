"""Fraction-based occupancy parsing and validation."""

from __future__ import annotations

import re
from fractions import Fraction


def parse_fraction(value: str | float | int) -> Fraction:
    """Parse an occupancy as an exact Fraction."""
    if isinstance(value, Fraction):
        return value
    if isinstance(value, int):
        return Fraction(value, 1)
    if isinstance(value, float):
        return Fraction(str(value)).limit_denominator(1000000)
    text = str(value).strip()
    if not text:
        raise ValueError("empty occupancy value")
    return Fraction(text)


def fraction_string(value: Fraction) -> str:
    """Return a compact exact string for a Fraction."""
    return str(value.numerator) if value.denominator == 1 else str(value)


def parse_set_occupancy(spec: str | None) -> dict[str, dict[str, Fraction]]:
    """Parse ``--set-occupancy`` declarations.

    Supported tokens:
      - ``Li1:1/2`` sets the existing row/group occupancy to 1/2.
      - ``M1:Ge1/2P1/2`` declares a mixed site composition.

    The parser intentionally does not guess labels or nearby fractions.
    """
    if not spec:
        return {}
    overrides: dict[str, dict[str, Fraction]] = {}
    for token in spec.split():
        if ":" not in token:
            raise ValueError(f"--set-occupancy token {token!r} lacks ':'")
        label, value = token.split(":", 1)
        label = label.strip()
        if not label:
            raise ValueError(f"--set-occupancy token {token!r} has empty label")
        overrides[label] = _parse_occupancy_value(value)
    return overrides


def _parse_occupancy_value(value: str) -> dict[str, Fraction]:
    value = value.strip()
    if not value:
        raise ValueError("empty --set-occupancy value")
    try:
        return {"__self__": parse_fraction(value)}
    except ValueError:
        pass

    parts = re.findall(r"([A-Z][a-z]?)([0-9]+(?:/[0-9]+)?|(?:0?\.\d+)|1(?:\.0+)?)", value)
    if not parts or "".join(el + occ for el, occ in parts) != value:
        raise ValueError(
            f"cannot parse occupancy composition {value!r}; "
            "use forms like '1/2' or 'Ge1/2P1/2'"
        )
    return {el: parse_fraction(occ) for el, occ in parts}


def check_integerizable(n_sites: int, occupancies: dict[str, Fraction],
                        group_label: str) -> None:
    """Raise if all occupancies cannot be represented by integer site counts."""
    for species, occ in occupancies.items():
        count = occ * n_sites
        if count.denominator != 1:
            raise ValueError(
                f"occupancy for {group_label}:{species} is not integerizable: "
                f"{n_sites} * {fraction_string(occ)} = {count}. "
                "Declare a commensurate --set-occupancy value or use a larger supercell."
            )
