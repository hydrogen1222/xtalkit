"""Shared utility functions for xtalkit."""

_DUMMY_ELEMENTS = ["Xe", "Kr", "Rn", "Ar", "Ne", "He"]


def assign_dummy_elements(
    wyckoff_letters: list[str],
    element_map: dict[str, str] | None,
) -> dict[str, str]:
    """Assign dummy elements to Wyckoff letters.

    Priority: Xe -> Kr -> Rn -> Ar -> Ne -> He (cycling if needed).

    When element_map is provided, every requested letter must have an entry.
    """
    if element_map is not None:
        for letter in wyckoff_letters:
            if letter not in element_map:
                raise ValueError(f"No element assigned for Wyckoff letter {letter}")
        return dict(element_map)

    assignment = {}
    for i, letter in enumerate(sorted(
        wyckoff_letters,
        key=lambda w: (int("".join(c for c in w if c.isdigit()) or 0), w),
    )):
        assignment[letter] = _DUMMY_ELEMENTS[i % len(_DUMMY_ELEMENTS)]
    return assignment


def parse_coord(s: str, variable_default: float = 0.3) -> float:
    """Parse a coordinate expression like '0', '1/4', '0.25' to float.

    Variable expressions ('x', 'y', 'z') are resolved to a representative
    value (default 0.3) for dummy atom placement.
    """
    s = s.strip()
    if s in ("x", "y", "z"):
        return variable_default
    if "/" in s:
        num, den = s.split("/")
        return float(num) / float(den)
    return float(s)
