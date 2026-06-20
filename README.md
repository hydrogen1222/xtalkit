# xtalkit — Crystal Wyckoff Toolkit

Mark Wyckoff positions in crystal structures with dummy atoms for intuitive visualization in [VESTA](https://jp-minerals.org/vesta/en/).

## Installation

Requires Python 3.10+.

```bash
git clone <repo-url> && cd xtalkit
uv sync
uv pip install -e .
```

Verify:

```bash
xtalkit --version
# xtalkit 0.1.0
```

## Quick Start

Mark two Wyckoff positions in a CIF file:

```bash
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f
# Output: Li6PS5Cl_WYCK.cif
```

Open the result in VESTA — dummy atoms (Xe, Kr, etc.) mark the Wyckoff sites with distinct colors.

---

# Usage

xtalkit has two interfaces:

| Interface | How to launch | Best for |
|-----------|--------------|----------|
| **CLI** | `xtalkit <command> ...` | Scripting, batch processing, one-liners |
| **TUI** | `xtalkit` (no arguments) | Interactive exploration, guided workflows |

---

## CLI Reference

```
xtalkit <command> [options]
```

Four subcommands:

| Command | Purpose |
|---------|---------|
| `mark` | Mark Wyckoff positions in a CIF file |
| `skeleton` | Generate a pure Wyckoff skeleton (no real atoms) |
| `info` | Query Wyckoff positions for a space group |
| `fetch` | Verify space group database integrity |

---

### `mark` — Mark Wyckoff Positions in a CIF

```
xtalkit mark <input.cif> --sg <N> --wyckoff <letters> [options]
```

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `input.cif` | Path to the CIF file (relative paths supported) |
| `--sg N` | Space group number, 1–230 |
| `--wyckoff L` | Wyckoff letters to mark, comma-separated (e.g. `4a,24f`) or `all` for every position |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--mode overlay` | overlay | `overlay`: keep real atoms + add dummies. `replace`: swap matched real atoms with dummies |
| `--tol 0.5` | 0.5 | Matching tolerance in fractional coordinates |
| `--map 4a:Xe,16e:Kr` | (auto) | Override dummy element assignment per Wyckoff letter |
| `--format cif` | cif | Output format(s): `cif`, `vesta`, `xyz`, or comma-separated `cif,vesta,xyz` |
| `-o base` | `{name}_WYCK` | Output base path (extension added per format) |

**Examples:**

```bash
# Mark 4a and 24f in a cubic F-43m structure (overlay mode)
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f

# Mark all 8 Wyckoff positions in all three output formats
xtalkit mark structure.cif --sg 216 --wyckoff all --format cif,vesta,xyz

# Replace mode: swap real atoms at 4a with Xe dummy atoms
xtalkit mark structure.cif --sg 216 --wyckoff 4a --mode replace

# Custom element mapping with tight tolerance
xtalkit mark NaCl.cif --sg 225 --wyckoff 4a,4b --map 4a:He,4b:Ne --tol 0.01

# Specify output path explicitly
xtalkit mark input.cif --sg 216 --wyckoff 4a -o ./output/marked
# Produces: ./output/marked.cif
```

**Output files:**

| Format | File | Notes |
|--------|------|-------|
| CIF | `{name}_WYCK.cif` | Full CIF with dummy atoms appended. Open in VESTA directly. |
| VESTA | `{name}_WYCK.vesta` | VESTA-native XML with cell + atom sites. |
| XYZ | `{name}_WYCK.xyz` | Simple Cartesian XYZ (loses cell info — VESTA will prompt for cell). |

---

### `skeleton` — Generate Wyckoff Skeleton

Generate a structure containing **only** dummy atoms at Wyckoff positions — no real atoms. Useful as a reference template.

```
xtalkit skeleton --sg <N> --wyckoff <letters> [options]
```

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `--sg N` | Space group number, 1–230 |
| `--wyckoff L` | Wyckoff letters, comma-separated or `all` |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--cell "a b c α β γ"` | System-based default | Custom cell parameters |
| `--map 4a:Xe,...` | (auto) | Element override |
| `--format cif` | cif | `cif`, `vesta`, `xyz`, or `cif,vesta,xyz` |
| `-o base` | `SG{N}_skeleton` | Output base path |

**Examples:**

```bash
# Skeleton for F-43m with default cubic cell (a=b=c=5.0)
xtalkit skeleton --sg 216 --wyckoff all

# Skeleton with real cell parameters for Li6PS5Cl
xtalkit skeleton --sg 216 --wyckoff 4a,4c,16e,24f \
    --cell "9.85 9.85 9.85 90 90 90"

# Skeleton for P2₁/c (monoclinic) with custom cell
xtalkit skeleton --sg 14 --wyckoff 2a,4e \
    --cell "5.5 6.3 8.1 90 108.5 90" --format cif,vesta
```

**Default cell parameters by crystal system:**

| System | a | b | c | α | β | γ |
|-----------|---|---|---|---|---|---|
| Triclinic | 5.0 | 6.0 | 7.0 | 80° | 90° | 100° |
| Monoclinic | 5.0 | 6.0 | 7.0 | 90° | 110° | 90° |
| Orthorhombic | 5.0 | 6.0 | 7.0 | 90° | 90° | 90° |
| Tetragonal | 5.0 | 5.0 | 7.0 | 90° | 90° | 90° |
| Trigonal | 5.0 | 5.0 | 8.0 | 90° | 90° | 120° |
| Hexagonal | 5.0 | 5.0 | 8.0 | 90° | 90° | 120° |
| Cubic | 5.0 | 5.0 | 5.0 | 90° | 90° | 90° |

**Always override with real cell parameters for accurate work.**

---

### `info` — Query Space Group Information

```
xtalkit info --sg <N>
```

**Example:**

```
$ xtalkit info --sg 216

Space Group #216: F-43m
Crystal System: cubic
Default cell: a=5.0 b=5.0 c=5.0 α=90 β=90 γ=90

Wyckoff Positions (8):
  Letter   Mult   Site Sym   Coordinates
  4a       4      -4         0,0,0
  4b       4      -4         1/2,1/2,1/2
  4c       4      -4         1/4,1/4,1/4
  4d       4      -4         3/4,3/4,3/4
  16e      16     .3m        x,x,x
  24f      24     2..        1/4,0,0
  24g      24     .2.        1/4,1/4,1/4
  48h      48     1          1/4,1/4,1/4
```

---

### `fetch` — Verify Database

```
xtalkit fetch
```

Verifies that all space group data is intact. Output:

```
✓ Space group data intact (230/230 OK)
```

---

## TUI (Interactive Mode)

Launch with no arguments:

```bash
xtalkit
```

```
╔═══════════════════════════════════════════╗
║      xtalkit · Crystal Wyckoff Toolkit   ║
╠═══════════════════════════════════════════╣
║  [1] Mark CIF    — Mark Wyckoff in a     ║
║                     structure            ║
║  [2] Skeleton    — Generate pure Wyckoff ║
║                     skeleton             ║
║  [3] Query SG    — View space group      ║
║                     information          ║
║  [4] Fetch DB    — Verify database       ║
║                     online               ║
║  [0] Exit                                ║
╚═══════════════════════════════════════════╝
```

### TUI Workflow Example — Mark CIF

```
Input your choice: 1

═══════════════════════════════════════════════
  Mark Wyckoff Positions in CIF

  CIF file path: ./Li6PS5Cl.cif
    ✓ found D:\structures\Li6PS5Cl.cif

  Space group number: 216

    Space Group #216: F-43m
    Available Wyckoff positions: 4a  4b  4c  4d  16e  24f  24g  48h

  Wyckoff positions to mark (comma-separated, or 'all'): 4a,24f

  Mode: [1] Overlay  [2] Replace  > 1

  Output format: [1] cif  [2] vesta  [3] xyz  [4] all  > 4

  Tolerance in Å (default 0.5): [Enter for default]

  Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip: [Enter]

  Output base path [default: D:\structures\Li6PS5Cl_WYCK]: [Enter]

  ✓ Saved to: D:\structures\Li6PS5Cl_WYCK.cif,
               D:\structures\Li6PS5Cl_WYCK.vesta,
               D:\structures\Li6PS5Cl_WYCK.xyz

  Press Enter to continue...
```

---

## Dummy Atom System

xtalkit uses **rare/inert elements** that almost never appear in real CIF files:

| Priority | Element | Z | VESTA Color |
|----------|---------|---|-------------|
| 1 | Xe (Xenon) | 54 | Silver |
| 2 | Kr (Krypton) | 36 | Light gray |
| 3 | Rn (Radon) | 86 | Pink |
| 4 | Ar (Argon) | 18 | Light blue |
| 5 | Ne (Neon) | 10 | Light green |
| 6 | He (Helium) | 2 | White |

Wyckoff letters are sorted alphabetically and assigned elements in priority order, cycling if the space group has more than 6 Wyckoff positions.

**CIF label format:** `WYCK_<letter>` (e.g., `WYCK_4a`, `WYCK_16e`).

**Custom mapping:**

```bash
xtalkit mark file.cif --sg 216 --wyckoff 4a,4c --map 4a:He,4c:Ne
```

---

## Overlay vs. Replace Modes

### Overlay (default)

```
Before:  Li at (0,0,0)   P at (0.25,0.25,0.25)
After:   Li at (0,0,0)   P at (0.25,0.25,0.25)
         WYCK_4a (Xe) at (0,0,0)
         WYCK_4c (Kr) at (0.25,0.25,0.25)
         ... more dummy atoms for other requested Wyckoff positions
```

Real atoms are preserved. Dummy atoms are added on top. In VESTA, you see both.

### Replace

```
Before:  Li at (0,0,0)   P at (0.25,0.25,0.25)
After:   WYCK_4a (Xe) at (0,0,0)    ← Li replaced (matched 4a)
         P at (0.25,0.25,0.25)      ← P stayed (didn't match 4c within tolerance)
```

Atoms that match a requested Wyckoff position (within tolerance) are **replaced** by dummy atoms. Other atoms are left untouched.

---

## Matching Tolerance

The tolerance (`--tol`, default 0.5) controls how closely an atom's coordinates must match a Wyckoff position's theoretical coordinates to be considered occupying that position.

- **Larger values** (0.5–1.0): lenient matching, suitable for experimental structures with slight coordinate deviations
- **Smaller values** (0.01–0.1): strict matching, only atoms very close to ideal Wyckoff positions are considered

The tolerance is applied in **fractional coordinate space** (not angstroms). For cubic cells with roughly equal axes, 0.5 in fractional space ≈ 0.5 × a in real distance along each axis.

---

## Workflow Recipes

### Recipe 1: Study Li₆PS₅Cl (SG 216, F-43m) Wyckoff occupancy

```bash
# Step 1: Download CIF from Materials Project
# (you already have the file)

# Step 2: See what Wyckoff positions exist in F-43m
xtalkit info --sg 216

# Step 3: Mark ALL Wyckoff positions, overlay mode, all formats
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff all --format cif,vesta,xyz

# Step 4: Open Li6PS5Cl_WYCK.vesta in VESTA
# → All 8 Wyckoff positions are now visible as colored dummy atoms
# → Real atoms (Li, P, S, Cl) are still shown
# → You can toggle atoms in VESTA to compare
```

### Recipe 2: Create a Wyckoff reference skeleton

```bash
# Generate a skeleton for F-43m with real cell parameters
xtalkit skeleton --sg 216 --wyckoff all \
    --cell "9.85 9.85 9.85 90 90 90" \
    --format vesta

# Open SG216_skeleton.vesta in VESTA
# → See exactly where each Wyckoff position sits in the unit cell
# → No real atoms — pure reference template
```

### Recipe 3: Check which atoms occupy specific Wyckoff positions

```bash
# Mark only the Wyckoff positions you care about
xtalkit mark structure.cif --sg 225 --wyckoff 4a,8c --mode replace

# In replace mode, atoms at 4a and 8c are swapped with dummy atoms
# → Instantly see in VESTA: "are there atoms at these positions?"
```

### Recipe 4: Batch process multiple structures

```bash
# All .cif files in a directory, same space group
for f in *.cif; do
    xtalkit mark "$f" --sg 216 --wyckoff all -o "${f%.cif}_WYCK"
done
```

---

## Supported Space Groups

Wyckoff position data is currently available for 38 space groups:

| Range | Crystal System | Count |
|-------|---------------|-------|
| 1–2 | Triclinic | 2 |
| 195–230 | Cubic | 36 |

Unsupported space groups raise `NotImplementedError` with a clear message. Expanding to all 230 space groups is planned.

**Most common battery materials (cubic SGs) are fully supported.**

---

## Development

```bash
uv sync          # Install dependencies
uv run pytest    # Run all tests (64 tests)
```

### Project structure

```
xtalkit/
├── xtalkit/
│   ├── __init__.py      # Package, version
│   ├── cli.py           # argparse CLI + 4 subcommands
│   ├── tui.py           # rich-based interactive TUI
│   ├── spacegroup.py    # Gemmi space group queries
│   ├── matcher.py       # Atom → Wyckoff position matching
│   ├── marker.py        # Core: mark Wyckoff in CIF
│   ├── skeleton.py      # Pure Wyckoff skeleton generation
│   ├── exporter.py      # .cif / .vesta / .xyz writers
│   └── utils.py         # Shared helpers
├── tests/
│   ├── fixtures/
│   │   └── simple.cif   # Test CIF (F-43m, Li + P)
│   ├── test_spacegroup.py
│   ├── test_matcher.py
│   ├── test_exporter.py
│   ├── test_marker.py
│   ├── test_skeleton.py
│   ├── test_cli.py
│   ├── test_tui.py
│   └── test_integration.py
├── docs/superpowers/
│   ├── specs/2026-06-20-xtalkit-design.md
│   └── plans/2026-06-20-xtalkit.md
├── pyproject.toml
└── README.md
```

---

## Dependencies

| Package | Purpose |
|---------|---------|
| [gemmi](https://gemmi.readthedocs.io/) | Space group data, CIF I/O |
| [rich](https://rich.readthedocs.io/) | TUI formatting (tables, panels, colors) |
| pytest (dev) | Test framework |
