**English** | [中文](README.zh-CN.md)

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

Five subcommands:

| Command | Purpose |
|---------|---------|
| `mark` | Mark Wyckoff positions in a CIF file |
| `skeleton` | Generate a pure Wyckoff skeleton (no real atoms) |
| `info` | Query Wyckoff positions for a space group |
| `fetch` | Verify space group database integrity |
| `enumerate` | Enumerate symmetry-inequivalent ordered configurations (requires pymatgen + enumlib) |

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

### `enumerate` — Enumerate Ordered Configurations

Enumerate all symmetry-inequivalent ordered configurations of a disordered CIF using [pymatgen](https://pymatgen.org/) + [enumlib](https://github.com/msg-byu/enumlib) (Hart-Forcade algorithm). This subcommand is **opt-in** — it requires the `enumerate` uv extra and source-compiled enumlib binaries (see [Enumeration setup](#enumeration-setup) below).

```
xtalkit enumerate <input.cif> [options]
```

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `input.cif` | Path to a CIF with partial/disordered occupancy (e.g. Au0.5/Cu0.5) |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--min-cell-size N` | 1 | Minimum supercell size to enumerate |
| `--max-cell-size N` | 2 | Maximum supercell size to enumerate |
| `--symm-prec TOL` | 0.1 | Symmetry tolerance for spacegroup analysis |
| `--vacancy-symbol S` | `X` | DummySpecies symbol for vacancies |
| `--output-dir DIR` | `{name}_enum/` | Output directory for enumerated CIFs |
| `--max-structures N` | unlimited | Cap the number of output files |
| `--timeout MIN` | none | Timeout in minutes for the enumlib subprocess |
| `--format F` | cif | Output format (`cif`; `xyz` reserved for future) |

**Examples:**

```bash
# Enumerate a 50/50 Au/Cu binary in a 1×–2× supercell
xtalkit enumerate AuCu_disordered.cif --max-cell-size 2
# Output: AuCu_disordered_enum/AuCu_disordered_000.cif, _001.cif, _002.cif

# Limit to first 5 structures
xtalkit enumerate disordered.cif --max-cell-size 3 --max-structures 5

# Custom output directory
xtalkit enumerate parent.cif --output-dir ./runs/exp01/
```

**How it works:**

1. Reads the CIF via `pymatgen.core.Structure.from_file` (handles partial occupancy + mmCIF natively)
2. **Converts to the primitive cell** (`Structure.get_primitive_structure()`). This is essential: enumlib does not auto-detect the primitive cell, and a conventional F-centered/I-centered cell can have 2×–4× the candidate sites, enough to overflow enumlib's tree_class array (e.g. F-43m 48h with 48 candidates → C(48,24) overflows; the primitive cell has only 12 candidates → C(12,6) = 924, exactly the regime that reproduces literature results).
3. Augments any partial-occupancy site (single- or multi-species) with an explicit `DummySpecies(vacancy_symbol)` to represent the vacancy (e.g. `Li0.5` → `Li0.5 + X0.5`; `Au0.3/Cu0.3` → `Au0.3 + Cu0.3 + X0.4`)
4. Calls `pymatgen.command_line.enumlib_caller.EnumlibAdaptor` which shells out to the Fortran `enum.x` and `makestr.x` binaries
5. Writes each symmetry-inequivalent ordered configuration as `<basename>_<NNN>.cif` (in the primitive cell)

#### How CIF partial occupancy works

In a CIF, the same crystallographic coordinate can be listed on **multiple atom_site rows**, each with its own `_atom_site_occupancy` value. A row whose occupancy < 1 means "this atom is here only that fraction of the time". Two cases:

**Case A — mixed occupancy** (sum of rows at one site = 1.0): describes *occupational disorder* — the site is always occupied, but by different species with given probabilities. The diffraction experiment sees a spatial/time average.

```
Au1 Au 0.0 0.0 0.0 0.5      ← 50% Au
Cu1 Cu 0.0 0.0 0.0 0.5      ← 50% Cu  (same coords, sums to 1.0 → site is full)
```

**Case B — single partial occupancy** (one row, occ < 1.0): describes a *true vacancy*. The remaining (1 − occ) is empty space.

```
Li1 Li 0.3148 0.018 0.6852 0.56    ← 56% Li, 44% vacancy (implicit)
```

xtalkit's `enumerate` handles both:
- **Case A** (multi-species mix, occupancies sum to 1.0): enumlib enumerates which species goes where directly — no augmentation needed. (If a multi-species site sums to < 1, xtalkit adds the shortfall as a vacancy species first.)
- **Case B** (single species + true vacancy): xtalkit auto-augments with an explicit `DummySpecies("X")` at `1 − occ`, turning it into Case A internally (`Li0.56 → Li0.56 + X0.44`), then enumlib enumerates Li-vs-vacancy orderings.

#### Why `Li6PS5Cl_clean.cif` exists (vs `EntryWithCollCode418490.cif`)

The argyrodite CIF from the literature (`EntryWithCollCode418490.cif`) reports Li on the 48h Wyckoff position with occupancy **0.56** (= Li6.72PS5Cl stoichiometry):

```
Li1 Li1+ 48 h 0.3148(19) 0.018(4) 0.6852(19) 0.104(14) 0.56(6) 0
```

`0.56 = 14/25` cannot be integerized in any practical supercell (would need 25× supercell = 15625 unit cells). enumlib returns 0 structures.

The "clean" version (`Li6PS5Cl_clean.cif`) replaces **only** that one occupancy value, rounding 0.56 → 0.5 (= Li6PS5Cl stoichiometry):

```
Li1 Li1+ 48 h 0.3148(19) 0.018(4) 0.6852(19) 0.104(14) 0.5 0
```

With 0.5 (= 1/2), a 1× supercell already has integer Li count (12 Li sites in the primitive cell × 0.5 = 6 Li), so enumlib runs cleanly and produces **48 symmetry-inequivalent ordered configurations** — matching the literature's result. Everything else in the CIF (cell, P/S/Cl positions, space group) is unchanged.

This rounding is a standard approximation in the argyrodite literature: the real material is Li6.72PS5Cl with positional Li disorder, but computational studies enumerate Li6PS5Cl (a nearby rational stoichiometry) because it's tractable. If you need the true Li6.72 stoichiometry, you'd need a much larger supercell or a different enumeration approach.

**Known limitations:**

- **Non-integer stoichiometry**: A site with occupancy 0.56 (= 14/25) cannot be integerized in any small supercell. enumlib will return 0 structures; xtalkit surfaces this with a clear error suggesting either a larger `--max-cell-size` or a "clean" parent CIF with rational occupancy (e.g. round 0.56 → 0.5). For Li6.72PS5Cl (Li occ 0.56), prepare a Li6PS5Cl parent CIF with Li occ 0.5 — this reproduces the ~48 Li orderings from the argyrodite literature.
- **Platform note**: `scripts/build_enumlib.sh` targets Linux/macOS (system `gfortran`). Windows users can compile under [WSL](https://learn.microsoft.com/windows/wsl/), or follow the legacy conda + `m2w64-gcc-fortran` path (place `enum.x`/`makestr.x` in the env's `Library/mingw-w64/bin`).

---

#### Enumeration setup

xtalkit's core (`mark`, `skeleton`, `info`, `fetch`) does **not** require pymatgen. Only `enumerate` does, and it's an opt-in uv extra — the core stays lightweight (gemmi + rich only). No conda, no root.

**Step 1 — Install the `enumerate` extra:**

```bash
uv sync --extra enumerate
```

This adds `pymatgen` (>=2024.5). Current PyPI builds of pymatgen still ship `pymatgen.command_line.enumlib_caller.EnumlibAdaptor`, so there is no need to pin an old version. (The `2023.5.31` pin seen in older docs was a conda-forge-specific workaround; uv uses PyPI, which is unaffected.)

**Step 2 — Compile the enumlib binaries (one-time, no root):**

```bash
bash scripts/build_enumlib.sh
```

This clones [msg-byu/enumlib](https://github.com/msg-byu/enumlib) (plus its `symlib` submodule), compiles `enum.x` and `makestr.x` with the system `gfortran`, and installs them to `~/.local/share/xtalkit/bin/`. Requires `gfortran`, `git`, and `make` (e.g. `sudo apt install gfortran make git`). Override the install location with `XTALKIT_ENUMLIB_BIN`.

xtalkit finds these binaries automatically — **no manual PATH setup needed**. `xtalkit._env.setup_for_enumlib()` prepends the install directory to `PATH` before `enumlib_caller`'s import-time `which("enum.x")` runs. It also checks `$XTALKIT_ENUMLIB_BIN`, and falls back to an in-repo `enumlib_src/enumlib/src/` clone for development.

**Step 3 — Verify:**

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif --max-cell-size 2
# Expect: 3 ordered CIFs from Au0.5/Cu0.5
uv run pytest tests/test_enumerator.py -v
# Expect: 8 passed
```

On Windows, `xtalkit/_env.py` additionally applies three runtime workarounds (these do not run on Linux/macOS):

1. Appends `.X` and `.PY` to `PATHEXT` so `shutil.which("enum.x")` finds the binary
2. Calls `os.add_dll_directory(env/Library/bin)` so scipy's native extension loads
3. Monkey-patches `shutil.which` to return absolute paths (Windows otherwise returns `.\makestr.x`, which `subprocess.Popen` cannot launch)

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
║  [5] Enumerate   — Enumerate ordered     ║
║                     configurations       ║
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
uv sync                     # Install core dependencies
uv sync --extra enumerate   # Also enable `enumerate` (pulls pymatgen)
uv run pytest               # Run all tests (5 skip without the enumerate extra)
```

### Project structure

```
xtalkit/
├── xtalkit/
│   ├── __init__.py      # Package, version
│   ├── cli.py           # argparse CLI + 5 subcommands
│   ├── tui.py           # rich-based interactive TUI
│   ├── spacegroup.py    # Gemmi space group queries
│   ├── matcher.py       # Atom → Wyckoff position matching
│   ├── marker.py        # Core: mark Wyckoff in CIF
│   ├── skeleton.py      # Pure Wyckoff skeleton generation
│   ├── exporter.py      # .cif / .vesta / .xyz writers
│   ├── enumerator.py    # enumlib wrapper (lazy pymatgen import)
│   ├── _env.py          # enumlib binary discovery + Windows env fixes
│   └── utils.py         # Shared helpers
├── tests/
│   ├── fixtures/
│   │   ├── simple.cif         # Test CIF (F-43m, Li + P)
│   │   └── disordered_binary.cif  # Au0.5/Cu0.5 for enumerate tests
│   ├── test_spacegroup.py
│   ├── test_matcher.py
│   ├── test_exporter.py
│   ├── test_marker.py
│   ├── test_skeleton.py
│   ├── test_enumerator.py     # enumlib integration (skips without pymatgen)
│   ├── test_cli.py
│   ├── test_tui.py
│   └── test_integration.py
├── docs/superpowers/
│   ├── specs/2026-06-20-xtalkit-design.md
│   └── plans/2026-06-20-xtalkit.md
├── scripts/
│   └── build_enumlib.sh    # Compile enumlib (enum.x, makestr.x) from source
├── pyproject.toml
└── README.md
```

---

## Dependencies

| Package | Purpose | Required? |
|---------|---------|-----------|
| [gemmi](https://gemmi.readthedocs.io/) | Space group data, CIF I/O | Yes |
| [rich](https://rich.readthedocs.io/) | TUI formatting (tables, panels, colors) | Yes |
| [pymatgen](https://pymatgen.org/) >=2024.5 | enumlib wrapper for `enumerate` | `enumerate` extra only |
| [enumlib](https://github.com/msg-byu/enumlib) | Symmetry-inequivalent configuration enumeration (Fortran) | Only for `enumerate` |
| pytest (dev) | Test framework | Yes |

The `enumerate` subcommand lazy-imports pymatgen, so the core toolkit works without it. See [Enumeration setup](#enumeration-setup) for the `uv sync --extra enumerate` + `build_enumlib.sh` path.
