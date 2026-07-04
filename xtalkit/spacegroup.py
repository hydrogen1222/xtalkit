"""Space group data queries via Gemmi.

Wyckoff positions for all 230 space groups are loaded from the bundled
``data/wyckoff.json`` dataset, which is generated at build time by
``scripts/build_wyckoff_db.py`` (derived from International Tables for
Crystallography Vol. A via pyxtal, MIT-licensed, and verified against gemmi's
symmetry operations — see the build script for details). The dataset uses
gemmi's reference setting, so coordinates are consistent with
``gemmi.SpaceGroup`` used elsewhere in the package.

Each Wyckoff position is stored as ``(letter, multiplicity, site_symmetry,
coordinates)`` where ``letter`` is the raw Wyckoff letter (e.g. ``"a"``); the
full label (``"4a"``) is assembled as ``f"{multiplicity}{letter}"`` in
``wyckoff_positions``.
"""

import json
import os
from collections import namedtuple

import gemmi

WyckoffInfo = namedtuple(
    "WyckoffInfo", ["letter", "multiplicity", "site_symmetry", "coordinates"]
)

# Default cell parameters per crystal system — rough templates used by
# `skeleton` and `build` when the user does not supply real values. Users
# should override with actual refined values for any serious work.
_DEFAULT_CELLS = {
    gemmi.CrystalSystem.Triclinic: (5.0, 6.0, 7.0, 80.0, 90.0, 100.0),
    gemmi.CrystalSystem.Monoclinic: (5.0, 6.0, 7.0, 90.0, 110.0, 90.0),
    gemmi.CrystalSystem.Orthorhombic: (5.0, 6.0, 7.0, 90.0, 90.0, 90.0),
    gemmi.CrystalSystem.Tetragonal: (5.0, 5.0, 7.0, 90.0, 90.0, 90.0),
    gemmi.CrystalSystem.Trigonal: (5.0, 5.0, 8.0, 90.0, 90.0, 120.0),
    gemmi.CrystalSystem.Hexagonal: (5.0, 5.0, 8.0, 90.0, 90.0, 120.0),
    gemmi.CrystalSystem.Cubic: (5.0, 5.0, 5.0, 90.0, 90.0, 90.0),
}

_DATA_PATH = os.path.join(os.path.dirname(__file__), "data", "wyckoff.json")
# Lazy-loaded: {sg_number: [(letter, multiplicity, site_symmetry, coordinates), ...]}.
_WYCKOFF_DB: dict[int, list[tuple[str, int, str, str]]] | None = None


def _load_wyckoff() -> dict[int, list[tuple[str, int, str, str]]]:
    """Load and cache the bundled Wyckoff dataset."""
    global _WYCKOFF_DB
    if _WYCKOFF_DB is None:
        with open(_DATA_PATH, encoding="utf-8") as f:
            raw = json.load(f)
        _WYCKOFF_DB = {
            int(sg): [(p["letter"], p["multiplicity"],
                       p["site_symmetry"], p["coordinates"]) for p in positions]
            for sg, positions in raw.items()
        }
    return _WYCKOFF_DB


def _get_sg(sg_number: int) -> gemmi.SpaceGroup:
    """Validate and return space group object."""
    if not 1 <= sg_number <= 230:
        raise ValueError(f"Space group number must be 1-230, got {sg_number}")
    return gemmi.SpaceGroup(sg_number)


def sg_name(sg_number: int) -> str:
    """Return the Hermann-Mauguin symbol for the space group."""
    sg = _get_sg(sg_number)
    return sg.short_name()


def crystal_system(sg_number: int) -> str:
    """Return the crystal system name (lowercase)."""
    sg = _get_sg(sg_number)
    return sg.crystal_system_str()


def wyckoff_positions(sg_number: int) -> list[WyckoffInfo]:
    """Return all Wyckoff positions for the given space group.

    Args:
        sg_number: ITA space group number (1-230).

    Returns:
        List of WyckoffInfo namedtuples with letter (full label, e.g.
        ``"4a"``), multiplicity, site_symmetry, and coordinates string,
        sorted by multiplicity then letter.
    """
    _get_sg(sg_number)  # validate number

    db = _load_wyckoff()
    raw = db.get(sg_number)
    if raw is None:
        raise NotImplementedError(
            f"Wyckoff data not available for space group {sg_number}"
        )

    result = []
    for letter_code, multiplicity, site_symmetry, coordinates in raw:
        result.append(WyckoffInfo(
            letter=f"{multiplicity}{letter_code}",
            multiplicity=multiplicity,
            site_symmetry=site_symmetry,
            coordinates=coordinates,
        ))
    # Sort by multiplicity first, then letter
    result.sort(key=lambda w: (
        int("".join(c for c in w.letter if c.isdigit()) or 0),
        w.letter,
    ))
    return result


def default_cell_params(sg_number: int) -> dict[str, float]:
    """Return default cell parameters for a space group.

    Uses crystal-system-based defaults. These are rough templates -
    users should override with real values via --cell for accurate work.

    Returns:
        dict with keys a, b, c, alpha, beta, gamma (all floats).
    """
    sg = _get_sg(sg_number)
    cs = sg.crystal_system
    a, b, c, alpha, beta, gamma = _DEFAULT_CELLS.get(
        cs, (5.0, 5.0, 5.0, 90.0, 90.0, 90.0)
    )
    return {"a": a, "b": b, "c": c, "alpha": alpha, "beta": beta, "gamma": gamma}
