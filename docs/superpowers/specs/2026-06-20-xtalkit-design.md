# xtalkit — Crystal Wyckoff Toolkit Design

**Date:** 2026-06-20
**Status:** approved

## Overview

xtalkit is a CLI + TUI toolkit for marking Wyckoff positions in crystal structures with dummy atoms, so they can be visually inspected in VESTA. It helps users learn and understand crystal structures, Wyckoff site occupancy, and international space groups.

Primary use cases:

1. **Mark a CIF file**: take a real structure (e.g., Li6PS5Cl), let user specify which Wyckoff positions to highlight, place dummy atoms at those positions, export for VESTA.
2. **Generate a Wyckoff skeleton**: given only a space group number, generate a "pure Wyckoff" structure (no real atoms, only dummy atoms at Wyckoff positions) as a reference template.
3. **Query space group info**: list Wyckoff positions, site symmetries, coordinates for a given space group.

## Design decisions

| Decision | Choice | Rationale |
|----------|--------|-----------|
| Name | xtalkit | crystal + toolkit, short |
| Form | CLI + TUI (rich) | Multiwfn-style menus + pure CLI passthrough |
| Package manager | uv | user preference |
| Core library | Gemmi | built-in 230 space groups with Wyckoff data, fast CIF I/O, CCP4-maintained |
| Dummy atoms | rare elements (Xe→Kr→Rn→Ar→Ne→He) | never appear in real CIFs, VESTA assigns distinct colors |
| Matching tolerance | 0.5 fractional default (configurable) | loose enough for experimental data from Materials Project |
| Output formats | .cif, .xyz | .cif is primary; .xyz as auxiliary |
| Space group data | Gemmi offline + `fetch` for updates | offline-first, network optional |
| Cell params (skeleton) | built-in typical values per SG + user override | each SG has one preset from a known material |

## Architecture

```
xtalkit/
├── xtalkit/
│   ├── __init__.py
│   ├── cli.py          # argparse CLI entry + TUI fallback (no args)
│   ├── tui.py          # Multiwfn-style interactive menus
│   ├── spacegroup.py   # space group queries via Gemmi + fetch
│   ├── marker.py       # core: Wyckoff marking orchestration
│   ├── skeleton.py     # pure Wyckoff skeleton generation
│   ├── exporter.py     # multi-format output (cif/xyz)
│   └── matcher.py      # atom coordinate → Wyckoff position matching
├── data/
│   └── cell_params.json  # typical cell params for all 230 SGs
├── tests/
│   └── fixtures/          # test CIF files
├── pyproject.toml
└── README.md
```

### Dependencies between modules

```
cli → tui → (marker → matcher + exporter, skeleton, spacegroup)
```

### Subcommands

```
xtalkit mark      # Mark Wyckoff positions in a CIF file (core)
xtalkit skeleton  # Generate pure Wyckoff skeleton (no real atoms)
xtalkit info      # Query space group Wyckoff information
xtalkit fetch     # Online sync of space group data
```

No subcommand launches the TUI.

## Data flow: `mark` command

```
Input CIF + space group number + Wyckoff selection
        │
        ▼
  [spacegroup.py]  Query Gemmi for Wyckoff positions of that SG
        │
        ▼
  [matcher.py]  For each atom in CIF, find nearest Wyckoff position
        │  Generate equivalent positions via symmetry operations
        │  Distance ≤ tolerance → match
        │  Output: {atom_label: wyckoff_letter} mapping
        ▼
  [marker.py]  Orchestrator
        │
        ├── overlay: keep real atoms + add dummy atoms at Wyckoff sites
        ├── replace: replace matched real atoms with dummy atoms
        │
        ▼
  [exporter.py]  Write .cif / .xyz
        │  Dummy atoms: WYCK_<letter> label, unique rare element per letter
        ▼
    Output files
```

## Data flow: `skeleton` command

```
Space group number + cell params (default or user) + Wyckoff selection
        │
        ▼
  [spacegroup.py]  Get all Wyckoff positions
        │
        ▼
  [skeleton.py]  Generate structure with only dummy atoms
        │  One dummy atom per Wyckoff position (at representative coordinate)
        │  Apply symmetry to generate all equivalent positions
        ▼
  [exporter.py]  Output
```

## Dummy atom assignment

Rare elements in priority order: Xe(54) → Kr(36) → Rn(86) → Ar(18) → Ne(10) → He(2)

- Wyckoff letters sorted alphabetically, assigned elements in order
- Cycle through the 6 elements if a SG has many Wyckoff positions
- CIF label format: `WYCK_<letter>` (e.g., `WYCK_4a`, `WYCK_16e`)
- User can override with `--map 4a:Xe,16e:Kr`

## Wyckoff selection

- User always specifies which Wyckoff letters to mark: `--wyckoff 4a,24f`
- `--wyckoff all` to mark all Wyckoff positions in the space group
- Invalid letters → error with list of valid letters for that SG

## Matching logic (matcher.py)

1. For each Wyckoff position, generate all equivalent sites via the SG symmetry operations
2. For each atom in the CIF, compute distance to every Wyckoff equivalent site
3. If minimum distance ≤ tolerance (default 0.5 fractional units), the atom occupies that Wyckoff position
4. Return mapping: which atoms are at which Wyckoff positions
5. Warn if some atoms have no match (possible SG mismatch or poor data)

## TUI design

Two-level menu structure. Main menu:

```
  ╔══════════════════════════════════════════╗
  ║        xtalkit · Crystal Wyckoff Toolkit║
  ╠══════════════════════════════════════════╣
  ║  [1] Mark CIF    — Mark Wyckoff in a structure    ║
  ║  [2] Skeleton    — Generate pure Wyckoff skeleton  ║
  ║  [3] Query SG    — View space group information    ║
  ║  [4] Fetch DB    — Update database online          ║
  ║  [0] Exit                               ║
  ╚══════════════════════════════════════════╝
```

### Mark workflow (select 1)

```
  > CIF file path: ./Li6PS5Cl.cif
    ✓ found (relative path resolved)

  > Space group number: 216

    Wyckoff positions in F-43m (No. 216):
    4a  4b  4c  4d  16e  24f  24g  48h

  > Wyckoff positions to mark (comma-separated, or 'all'): 4a,24f

  > Mode: [1] Overlay  [2] Replace  > 1

  > Output format: [1] cif  [2] xyz  [3] all  > 3

  > Tolerance (default 0.5, fractional): [Enter for default]

  ✓ Done. Saved to:
    ./Li6PS5Cl_WYCK.cif
    ./Li6PS5Cl_WYCK.xyz
```

### Skeleton workflow (select 2)

```
  > Space group number: 216

  > Wyckoff positions (comma-separated, or 'all'): all

  > Cell parameters:
    Default for SG #216: a=5.905 b=5.905 c=5.905 α=90 β=90 γ=90
    [1] Use default  [2] Enter manually  > 1

  > Output format: [1] cif  [2] xyz  [3] all  > 1

  ✓ Done. Saved to: ./SG216_skeleton.cif
```

### Info workflow (select 3)

```
  > Space group number: 216

  Space Group #216: F-43m
  Crystal System: cubic
  Wyckoff Positions:
    Label    Site Sym    Coordinates
    4a       -4          0,0,0
    4b       -4          1/2,1/2,1/2
    4c       -4          1/4,1/4,1/4
    4d       -4          3/4,3/4,3/4
    16e      .3m         x,x,x
    24f      ..m         x,0,0
    24g      ..m         1/4,1/4,z
    48h      1..         x,x,z
```

## Output formats

| Format | Implementation | Notes |
|--------|---------------|-------|
| .cif | Gemmi native CIF writer | Full symmetry info preserved |
| .xyz | Plain text writer | Element + xyz, loses cell info; auxiliary only |

Naming: `{original_name}_WYCK.{ext}` for mark, `SG{num}_skeleton.{ext}` for skeleton.

## Error handling

| Scenario | Behavior |
|----------|----------|
| CIF file not found | Red error, allow re-input (TUI) or exit with message (CLI) |
| Invalid SG number (not 1–230) | Immediate error, prompt re-input |
| Invalid Wyckoff letter for SG | Show valid letters, prompt re-input |
| No atoms matched within tolerance | Warning, continue processing |
| Gemmi not installed / data corrupt | Startup check, clear install instructions |

## Configuration

- `--tol` / `TOLERANCE` env var: matching tolerance in fractional coordinate units (default 0.5)
- `--map 4a:Xe,16e:Rn`: override dummy element assignment
- `--cell a b c α β γ`: override cell parameters for skeleton

## Testing

- **Unit tests** (pytest): matcher logic, spacegroup queries, exporter format, dummy element assignment
- **Integration tests**: process real CIF files (Li6PS5Cl, NaCl, Si), verify output CIFs parseable by Gemmi
- **Snapshot tests**: fixed inputs → fixed output comparison (prevent regressions)
- Test fixtures in `tests/fixtures/`

## Dependencies

- `gemmi` — space group data + CIF I/O
- `rich` — TUI table/color formatting
- `argparse` — CLI (stdlib)
- `pytest` — testing (dev)
