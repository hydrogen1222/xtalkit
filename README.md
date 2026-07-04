**English** | [中文](README.zh-CN.md)

# xtalkit — Crystal Wyckoff Toolkit

A command-line toolkit for **periodic crystals of any kind** — not limited to solid-state electrolytes. (The maintainer works on argyrodite / LGPS-type Li conductors, but xtalkit handles alloys, minerals, ceramics, MOFs — anything with a space group.) It does four jobs:

- **Visualize Wyckoff positions** — mark crystallographic sites with dummy atoms for [VESTA](https://jp-minerals.org/vesta/en/).
- **Build CIFs from refinement data** — assemble a standards-compliant CIF from a space group, unit cell, and site list.
- **Enumerate occupational disorder** — generate all symmetry-inequivalent orderings of a partially occupied structure, via two backends:
  - **enumlib** (Hart–Forcade) — fast, for small / medium cells.
  - **SHRY** (Pólya enumeration) — strict, audited, for large partially occupied cells where enumlib exhausts memory.
- **Batch Ewald ranking** — score many CIF/POSCAR structures by Coulombic Ewald energy, then sort and group them for pre-screening.

## Table of Contents

- [Installation](#installation)
- [Quick Start](#quick-start)
- [Command Overview](#command-overview)
- [Core Commands](#core-commands) — `mark`, `skeleton`, `build`, `info`, `fetch`
- [Batch Ewald Ranking](#batch-ewald-ranking)
- [Enumerating Occupational Disorder](#enumerating-occupational-disorder) — enumlib vs SHRY, when to use which
- [TUI (Interactive Mode)](#tui-interactive-mode)
- [Reference](#reference) — dummy atoms, modes, tolerance, space groups
- [Development](#development)

## Installation

Requires Python 3.10+ and [uv](https://docs.astral.sh/uv/). xtalkit is split into **three independent tiers** — install only what you need.

### Tier 1 — Core (`mark` / `skeleton` / `build` / `info` / `fetch` / TUI)

No compilation, no heavy dependencies (gemmi + rich only):

```bash
git clone https://github.com/hydrogen1222/xtalkit.git && cd xtalkit
uv sync
uv pip install -e .
xtalkit --version          # xtalkit 0.1.0
```

### Tier 2 — +pymatgen features (`xtalkit enumerate`, `xtalkit ewald`, SHRY postprocess ranking)

Adds pymatgen + source-compiled enumlib Fortran binaries. Prerequisites: `gfortran make git` (e.g. `sudo apt install gfortran make git` on Debian/Ubuntu; on Windows use [WSL](https://learn.microsoft.com/windows/wsl/)).

```bash
uv sync --extra enumerate          # adds pymatgen (>=2024.5)
bash scripts/build_enumlib.sh      # compiles enum.x + makestr.x (one-time, no root)
```

`scripts/build_enumlib.sh` clones [msg-byu/enumlib](https://github.com/msg-byu/enumlib), builds with `gfortran`, and installs to `~/.local/share/xtalkit/bin/`. xtalkit finds the binaries automatically — **no PATH setup needed**. Verify:

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif --max-cell-size 2
# Expect: 3 ordered CIFs from Au0.5/Cu0.5
```

The same `uv sync --extra enumerate` step also enables `xtalkit ewald`, which uses pymatgen's Ewald summation backend for batch ranking.

### Tier 3 — +SHRY enumeration (`xtalkit shry`)

`xtalkit shry` shells out to the external [SHRY](https://pypi.org/project/shry/) CLI (pure Python, no compilation). **SHRY pins `pymatgen<=2023.10.4`, which conflicts with Tier 2's `pymatgen>=2024.5`** — so install SHRY in an *isolated* environment, not the xtalkit venv:

```bash
uv tool install shry                       # isolated env; `shry` lands on PATH
export XTALKIT_SHRY_CMD="$(which shry)"    # point xtalkit at it (add to shell profile)
shry --version                             # verify
```

> **Why isolated?** If you `uv pip install shry` into the xtalkit venv, it downgrades pymatgen to 2023.10.4 and `uv sync --extra enumerate` will then fight to upgrade it — the two backends break each other. `uv tool install` gives SHRY its own dependency set; xtalkit calls it as a subprocess via `XTALKIT_SHRY_CMD`. This is by design.

> **Full strict audit also needs Tier 2.** `shry prepare`/`verify` run an independent spglib symmetry audit (plan §3.4) and a two-stage `StructureMatcher` dedup (§10) inside the xtalkit process, which need `pymatgen`+`spglib`. So for the full strict workflow also run `uv sync --extra enumerate` (Tier 2's pymatgen brings spglib along). Without it, `shry count`/`enum` still work (they only shell out to SHRY), but `prepare`/`verify` raise a clear "spglib/pymatgen required" error.

### Installation summary

| Tier | Commands enabled | Extra steps | Conflicts |
|------|------------------|-------------|-----------|
| 1 · Core | `mark`, `skeleton`, `build`, `info`, `fetch`, TUI | none | — |
| 2 · +pymatgen | `enumerate`, `ewald`, SHRY postprocess ranking | `uv sync --extra enumerate` + `build_enumlib.sh` (gfortran) | — |
| 3 · +SHRY | + `shry` | `uv tool install shry` + `XTALKIT_SHRY_CMD`; full audit also `uv sync --extra enumerate` | SHRY's `pymatgen<=2023.10.4` pin — isolated, so none |

Tiers 2 and 3 are **independent** — install either or both.

## Quick Start

Mark two Wyckoff positions in a CIF file:

```bash
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f
# Output: Li6PS5Cl_WYCK.cif
```

Open the result in VESTA — dummy atoms (Xe, Kr, etc.) mark the Wyckoff sites with distinct colors.

---

## Command Overview

xtalkit has two interfaces:

| Interface | How to launch | Best for |
|-----------|--------------|----------|
| **CLI** | `xtalkit <command> ...` | Scripting, batch processing, one-liners |
| **TUI** | `xtalkit` (no arguments) | Interactive exploration, guided workflows |

```
xtalkit <command> [options]
```

Eight subcommands, in two groups:

| Command | Purpose | Tier |
|---------|---------|------|
| `mark` | Mark Wyckoff positions in a CIF file | Core |
| `skeleton` | Generate a pure Wyckoff skeleton (no real atoms) | Core |
| `build` | Build a CIF from refinement parameters (SG + cell + Wyckoff sites + occupancy) | Core |
| `info` | Query Wyckoff positions for a space group | Core |
| `fetch` | Verify space group database integrity | Core |
| `enumerate` | Enumerate symmetry-inequivalent ordered configurations (pymatgen + enumlib) | +enumlib |
| `shry` | Strict SHRY workflow for large partial-occupancy enumeration | +SHRY |
| `ewald` | Batch Ewald electrostatic energy ranking and grouping | +pymatgen |

The **Core** commands (Tier 1) are documented next. The two enumeration backends (`enumerate`, `shry`) are documented together under [Enumerating Occupational Disorder](#enumerating-occupational-disorder). The batch Ewald workflow is documented below in [Batch Ewald Ranking](#batch-ewald-ranking).

---

## Batch Ewald Ranking

`xtalkit ewald` scores many structures by Coulombic Ewald energy, sorts them, optionally keeps the top N, and can copy the selected and remaining structures into separate folders.

```
xtalkit ewald <path> [<path> ...] [options]
```

It supports two directory layouts:

1. `flat` - a directory containing many structure files directly, such as `1.cif`, `2.cif`, `3.cif`.
2. `nested` - a directory containing subdirectories, each of which contains structure files one level deep.

**Common options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--layout flat` | flat | Choose `flat` or `nested` directory scanning |
| `--sort asc` | asc | `asc` = lowest energy first, `desc` = highest energy first |
| `--top-n N` | all | Keep only the first N structures after sorting |
| `--charges Li:1 ...` | (none) | Explicit oxidation states for Ewald scoring |
| `--guess` | off | Let pymatgen guess oxidation states |
| `--per-atom` | off | Rank by energy per atom instead of total energy |
| `--out path` | `<input>_ewald/ranking.csv` | Write the full ranking CSV |
| `--group` | off | Copy the selected and remaining structures into folders |
| `--group-dir dir` | `<input>_ewald` | Grouping root directory |
| `--move` | off | Move files instead of copying them |

**Examples:**

```bash
# Flat folder: many CIF files directly in one directory
xtalkit ewald ./cif_pool \
  --layout flat \
  --charges Li:1 Ge:4 P:5 S:-2 \
  --sort asc --top-n 100 --group

# One-level nested folder tree
xtalkit ewald ./structures \
  --layout nested \
  --guess \
  --per-atom \
  --sort desc --top-n 50
```

By default the command writes `ranking.csv` and, if `--group` is set, creates:

```text
<output-dir>/
  ranking.csv
  selected/
  rest/
```

For `nested` mode the original one-level folder structure is preserved under `selected/` and `rest/`.

---

## Core Commands

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
| `--format cif` | cif | Output format(s): `cif`, `xyz`, or comma-separated `cif,xyz` |
| `-o base` | `{name}_WYCK` | Output base path (extension added per format) |

**Examples:**

```bash
# Mark 4a and 24f in a cubic F-43m structure (overlay mode)
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff 4a,24f

# Mark all 8 Wyckoff positions in all output formats
xtalkit mark structure.cif --sg 216 --wyckoff all --format cif,xyz

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
| `--format cif` | cif | `cif`, `xyz`, or `cif,xyz` |
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
    --cell "5.5 6.3 8.1 90 108.5 90" --format cif,xyz
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

### `build` — Build a CIF from Refinement Parameters

```
xtalkit build --sg <N> --cell "<a b c α β γ>" --atom "<spec>" [--atom ...] [options]
```

Assemble a standards-compliant CIF from XRD-refinement results: a space group, a unit cell, and a list of atomic sites. Each site is an element on a Wyckoff position (e.g. `16e`) with its free-coordinate values and an optional occupancy. The crystal system is **derived** from the space group (not asked separately) and used only to sanity-check the cell. Occupancy defaults to 1.0; partial/mixed occupancy is supported and produces a disordered CIF that `xtalkit enumerate` can later order.

All 230 space groups are supported (the bundled Wyckoff dataset, derived from International Tables via pyxtal and verified against gemmi).

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `--sg N` | Space group number, 1–230 |
| `--cell "a b c α β γ"` | Unit-cell parameters |
| `--atom SPEC` | Atom spec (repeatable): `element wyckoff [free...] [occ]` |
| `--spec FILE` | JSON spec file (alternative to `--sg`/`--cell`/`--atom`) |

**Atom spec format:** `element wyckoff [free-coordinate values] [occupancy]`. The free values are entered in the order the variables appear in the position's coordinate template; the parser knows how many each Wyckoff site expects. The last number is taken as occupancy iff one extra is given.

| Spec | Meaning |
|------|---------|
| `Na 4a` | Na on 4a (no free params), occupancy 1.0 |
| `Li 16e 0.3` | Li on 16e `(x,x,x)` with x=0.3, occupancy 1.0 |
| `Li 16e 0.3 0.5` | Li on 16e with x=0.3, occupancy 0.5 |
| `S 48h 0.25 0.3` | S on 48h `(x,x,z)` with x=0.25, z=0.3, occupancy 1.0 |
| `Li 4a 0.5` + `Cu 4a 0.5` | Mixed Li/Cu occupancy on 4a (disorder) |

**Options:** `--format cif[,xyz]` (default `cif`), `-o/--output` (base path; default `SG{N}_built`).

**Fractional-coordinate mode (`--atom-frac`)** — easier when you have a refinement table. Instead of Wyckoff letters + free parameters, give each atom's final fractional coordinates directly: `element x y z [occ]`. The tool auto-detects which Wyckoff orbit each atom sits on (and reports it), so non-canonical representatives (e.g. `(0, 1/2, z)` instead of the canonical `(1/2, 0, 1/4+z)` for SG 137 4d — related by a 4-fold rotation) work as-is, with no offset arithmetic. Atoms with occupancy 0 are skipped.

```bash
xtalkit build --sg 137 --cell "8.694 8.694 12.599 90 90 90" \
    --atom-frac "Li 0.2563 0.2718 0.1832 0.691" \
    --atom-frac "Li 0 0.5 0.9446 1" \
    --atom-frac "S 0 0.1843 0.4103 1" \
    ...
```

**Example — build NaCl (Fm-3m, 225):**

```bash
xtalkit build --sg 225 --cell "5.64 5.64 5.64 90 90 90" \
    --atom "Na 4a" --atom "Cl 4b" -o NaCl
# [OK] SG #225 (Fm-3m), crystal system: cubic, formula: NaCl
#      Saved to: NaCl.cif
```

**Example — build a disordered Li/Cu site, then enumerate its orderings:**

```bash
xtalkit build --sg 225 --cell "4 4 4 90 90 90" \
    --atom "Li 4a 0.5" --atom "Cu 4a 0.5" -o AuCu_disordered
xtalkit enumerate AuCu_disordered.cif --max-cell-size 2
```

The CIF lists only the asymmetric-unit representatives plus the space group's symmetry operations; VESTA, gemmi, and pymatgen expand them to the full unit cell. After writing, xtalkit prints the formula (e.g. `NaCl`, `Li6PS5Cl`) computed from multiplicity × occupancy — compare it against your refinement to confirm the input.

**JSON spec** (`--spec`) for reproducible/scriptable builds:

```json
{"sg": 225,
 "cell": {"a": 5.64, "b": 5.64, "c": 5.64, "alpha": 90, "beta": 90, "gamma": 90},
 "atoms": [
   {"element": "Na", "wyckoff": "4a", "occ": 1.0},
   {"element": "Cl", "wyckoff": "4b", "occ": 1.0}
 ]}
```

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
  4a       4      -43m       0,0,0
  4b       4      -43m       1/2,1/2,1/2
  4c       4      -43m       1/4,1/4,1/4
  4d       4      -43m       3/4,3/4,3/4
  16e      16     .3m        x,x,x
  24f      24     2.mm       x,0,0
  24g      24     2.mm       x,0,1/2
  48h      48     ..m        x,x,z
```

---

### `fetch` — Verify Database

```
xtalkit fetch
```

Verifies that all bundled space group data is intact. Output:

```
[OK] Space group data intact (230/230 space groups supported).
```

---

## Enumerating Occupational Disorder

Both `xtalkit enumerate` (enumlib) and `xtalkit shry` (SHRY) solve the same problem: given a CIF with partial occupancy, generate every **symmetry-inequivalent** ordering of the disordered sites. Both collapse symmetry-equivalent configs into one, so you get only genuinely distinct structures. They differ in backend, scale, and safety:

| | `xtalkit enumerate` (enumlib) | `xtalkit shry` (SHRY) |
|---|---|---|
| **Algorithm** | Hart–Forcade (supercell + enum.x backtracking) | Pólya enumeration theorem (cycle index, exact count) |
| **Best for** | Small/medium primitive cells; binary alloys, simple disorder | Large partially-occupied cells; argyrodite/LGPS-type Li disorder |
| **Count before generating?** | No — runs the full search, then writes | Yes — `shry count` is a cheap safety gate before `shry enum` |
| **Cell** | Auto-converts to primitive internally | Uses the CIF as given (supply a primitive CIF for F/I-centered cells) |
| **Audit trail** | Minimal | Per-stage manifests, formula/dedup/degeneracy verification |
| **Failure mode** | Non-integer stoichiometry can make `enum.x` run away and OOM | Refuses cleanly; `--expect-count` locks generation to the count |
| **Install** | `uv sync --extra enumerate` + `build_enumlib.sh` (gfortran) | `uv tool install shry` + `XTALKIT_SHRY_CMD` |
| **Rigor** | Hart–Forcade is exact; pymatgen-internal dedup | Pólya count is exact; generation is canonical-form (one rep per orbit) |

**Which should I use?**

- **Small cell, want it fast** → `enumerate`. It auto-primitivizes and is well-tested.
- **Large cell, partial occupancy, want a safety gate** → `shry`. Count first; if the number is sane, generate. The strict workflow refuses to write unless count and generation agree.
- **F-/I-centered cell with many disordered sites** → `shry` with a **primitive** CIF (enumlib auto-primitivizes, but `shry` does not — see the SHRY section below).

Both backends have been cross-checked against each other and against brute-force on LPSC (48 inequivalent Li/vacancy orderings; degeneracy sum = 924 = C(12,6)).

---

### `enumerate` — Enumerate Ordered Configurations (enumlib)

Enumerate all **symmetry-inequivalent** ordered configurations of a disordered CIF using [pymatgen](https://pymatgen.org/) + [enumlib](https://github.com/msg-byu/enumlib) (Hart-Forcade algorithm). Structures related by a symmetry operation of the parent are automatically collapsed into one — you get only the genuinely distinct orderings, with no manual deduplication.

This works for any material with **occupational disorder** (partial occupancy on crystallographic sites): binary alloys (Au/Cu), solid solutions, doped systems, and battery electrolytes with mobile-ion disorder (e.g. argyrodite Li₆PS₅Cl, LGPS). It is **opt-in** — it requires the `enumerate` uv extra and source-compiled enumlib binaries (see [Enumeration setup](#enumeration-setup)).

```
xtalkit enumerate <input.cif> [options]
```

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `input.cif` | Path to a CIF with partial/disordered occupancy (e.g. Au0.5/Cu0.5, or Li0.5 + vacancy) |

**Options:**

| Option | Default | Description |
|--------|---------|-------------|
| `--min-cell-size N` | 1 | Minimum supercell size to enumerate |
| `--max-cell-size N` | 2 | Maximum supercell size to enumerate (larger = more orderings + slower; see below) |
| `--symm-prec TOL` | 0.1 | Symmetry tolerance for spacegroup analysis |
| `--vacancy-symbol S` | `X` | DummySpecies symbol used internally for vacancies |
| `--output-dir DIR` | `{name}_enum/` | Output directory for enumerated files |
| `--max-structures N` | unlimited | Cap output and post-enumeration `makestr.x` generation (does not cap `enum.x`) |
| `--timeout MIN` | none | Timeout in minutes for the enumlib subprocess |
| `--format F` | cif | Output format: `cif` or `xyz` |
| `--jobs N` | 1 | Parallel workers for structure generation (`0` = auto/cpu_count) |
| `--batch-size N` | 256 | Structures per batch — caps peak memory |
| `--scratch-dir DIR` | system temp | Directory for enumlib scratch files (try `/dev/shm` for tmpfs) |
| `--skip-preflight` | off | Skip the non-integer-stoichiometry pre-check (dangerous — can OOM) |

**Examples:**

```bash
# Enumerate a 50/50 Au/Cu binary — quick first look with a 1× supercell
xtalkit enumerate AuCu_disordered.cif --max-cell-size 1

# Full 1×–2× enumeration, keep only the first 5 structures
xtalkit enumerate disordered.cif --max-cell-size 2 --max-structures 5

# Enumerate a battery-electrolyte CIF, writing to a custom directory
xtalkit enumerate Li6PS5Cl.cif --max-cell-size 1 --output-dir ./runs/exp01/
```

> **What gets enumerated?** You never tell xtalkit which element to enumerate — there is no flag for it. xtalkit auto-detects from the CIF and enumerates **every** site whose occupancy is < 1 (or that carries a species mix); sites at full occupancy (1.0) are held fixed. If the CIF has several disordered sites (say a Li/vacancy site *and* a Ge/P mixed site), they are enumerated **jointly** — each output structure fixes all of them at once. So `enumerate Li6PS5Cl.cif` enumerates the Li site (Li/vacancy); `enumerate AuCu_disordered.cif` enumerates the Au/Cu site.
>
> **Choosing what to study.** Because enumeration is driven entirely by the CIF, you select the disorder under study by what you put in the CIF — not by a command-line flag. To study only Li/vacancy, leave the other sites fully ordered (occupancy 1); to study only Ge/P mixing, leave Li fully occupied. Conversely, a fully-ordered CIF (every site at occupancy 1) has nothing to enumerate and yields 0 structures — you must edit the CIF to introduce the disorder you want first (see *When enumlib returns 0 structures* below).

**How it works:**

1. Reads the CIF via `pymatgen.core.Structure.from_file` (handles partial occupancy + mmCIF natively).
2. **Converts to the primitive cell** (`Structure.get_primitive_structure()`). This matters: enumlib does not auto-detect the primitive cell, and a conventional F-/I-centered cell can have 2×–4× the candidate sites of the primitive cell — enough to make enumeration impractically slow or to overflow enumlib's internal arrays. Working in the primitive cell keeps the candidate count small.
3. Augments any partial-occupancy site with an explicit vacancy species (`DummySpecies`) so enumlib has a concrete species to place, e.g. `Li0.5` → `Li0.5 + X0.5`; `Au0.3/Cu0.3` → `Au0.3 + Cu0.3 + X0.4`. **Fully-occupied sites (occupancy 1) are left untouched** — they are fixed spectators; only the disordered sites get enumerated.
4. Calls `pymatgen.command_line.enumlib_caller.EnumlibAdaptor`, which shells out to the Fortran `enum.x` and `makestr.x` binaries.
5. Writes each distinct ordered configuration as `<basename>_<NNN>.cif` (in the primitive-cell setting, reported as P1).

#### Choosing `--max-cell-size`

enumlib enumerates ordered configurations in **supercells** of the primitive cell, from `--min-cell-size` up to `--max-cell-size`. A larger max size explores bigger supercells and finds more orderings — but the count (and runtime) can grow explosively.

- **Start small.** For a first look, use `--max-cell-size 1`. If the parent's occupancy is already integerizable in the primitive cell (e.g. 1/2, 1/3, 2/3), the 1× cell already yields the full set of distinct orderings and finishes in seconds.
- **Go bigger only if needed.** A 2× supercell can reveal additional orderings that don't fit in the 1× cell, but it can also produce thousands of structures and run for minutes to hours. Pair it with `--max-structures` and `--timeout` to stay bounded.
- **Occupancy must be integerizable.** The count of each species in the chosen supercell must come out integer. Occupancy 1/2 needs at least 2 sites (a 1× cell of a 2-site primitive works); 1/3 needs a 3× cell; a value like 0.56 (= 14/25) would need a 25× cell, which is impractical. If your CIF has such a value, round it to a nearby simple fraction first (see below).

#### Memory, disk, and parallelism

A large enumeration (thousands of structures) can use a lot of memory and disk. xtalkit processes structures in **streamed batches** rather than holding them all in memory, and exposes three knobs:

- **`--max-structures N`** caps the post-enumeration `makestr.x` + parse + write phase, not the `enum.x` search itself. It avoids generating thousands of intermediate `vasp.*`/CIF files when you only need a sample, but `enum.x` must still finish the full symmetry search first.
- **`--batch-size N`** sets how many structures are generated, parsed, and written per batch. Peak memory scales with the batch size, not the total count. Lower it (e.g. `--batch-size 32`) if you run out of memory; raise it for slightly fewer `makestr.x` invocations.
- **`--jobs N`** runs batches in parallel worker processes, speeding up the structure-generation (makestr + parse + write) phase. `--jobs 0` auto-uses all CPUs. **Each worker loads its own copy of pymatgen (~200 MB), so memory scales with `--jobs`** — pick a value that fits your RAM.
- **`--scratch-dir /dev/shm`** puts enumlib's scratch files (`struct_enum.out`, per-batch `vasp.*`) on tmpfs, avoiding disk I/O. `/dev/shm` is usually capped at ~half of RAM, so only use it if `struct_enum.out` fits.

> **`enum.x` is single-threaded.** The enumeration core (`enum.x`) is a serial Fortran backtracking search — `--jobs` does **not** speed it up. `--jobs` parallelises only the makestr + parse + write phase that follows. If `enum.x` dominates (long enumeration of a huge search space), the wall-clock improvement from `--jobs` will be modest; streaming and `--max-structures` only help after `enum.x` has produced `struct_enum.out`. Parallelising `enum.x` itself would require splitting by cell size, which is not implemented.

#### How CIF partial occupancy works

In a CIF, the same crystallographic coordinate can be listed on **multiple atom_site rows**, each with its own `_atom_site_occupancy` value. A row whose occupancy < 1 means "this atom is here only that fraction of the time". Two cases:

**Case A — mixed occupancy** (sum of rows at one site = 1.0): describes *occupational disorder* — the site is always occupied, but by different species with given probabilities. The diffraction experiment sees a spatial/time average.

```
Au1 Au 0.0 0.0 0.0 0.5      ← 50% Au
Cu1 Cu 0.0 0.0 0.0 0.5      ← 50% Cu  (same coords, sums to 1.0 → site is full)
```

**Case B — single partial occupancy** (one row, occ < 1.0): describes a *true vacancy*. The remaining (1 − occ) is empty space.

```
Li1 Li 0.0 0.0 0.0 0.5      ← 50% Li, 50% vacancy (implicit)
```

xtalkit's `enumerate` handles both:
- **Case A** (multi-species mix, occupancies sum to 1.0): enumlib enumerates which species goes where directly — no augmentation needed. (If a multi-species site sums to < 1, xtalkit adds the shortfall as a vacancy species first.)
- **Case B** (single species + true vacancy): xtalkit auto-augments with an explicit `DummySpecies("X")` at `1 − occ`, turning it into Case A internally (`Li0.5 → Li0.5 + X0.5`), then enumlib enumerates Li-vs-vacancy orderings.

#### Worked example: introducing disorder into an ordered CIF

A CIF from a database (e.g. Materials Project) is usually fully ordered — every site at occupancy 1. `xtalkit enumerate` on such a CIF yields 0 structures, because there is nothing to enumerate. To study disorder you edit the CIF first. Two edit patterns cover every case:

**Edit 1 — single-species partial occupancy (Li/vacancy, Case B).** Change that site's occupancy from `1` to a fraction; xtalkit adds the vacancy automatically.

```
Li  Li0  8  0.229  0.273  0.295  1        →        Li  Li0  8  0.229  0.273  0.295  0.5
```

**Edit 2 — multi-species mixed occupancy (Ge/P, Case A).** Replace the single row with two rows at the *same* fractional coordinates, occupancies summing to 1 (pymatgen merges same-coordinate rows into one mixed site).

```
Ge  Ge4  2  0.5  0.5  0.301  1        →        Ge  Ge4   2  0.5  0.5  0.301  0.5
                                               P   P4a   2  0.5  0.5  0.301  0.5
```

**The occupancy must be integerizable.** For a site of multiplicity M at occupancy p, the count `M × p` must be an integer — that's how many positions enumlib actually fills. On an M=8 site, valid p are 1/8, 1/4, 3/8, 1/2, …; on an M=2 site, only 1/2 (or 1). `0.5` is just the simplest value that works on both — not the only choice. If your desired p doesn't integerize in 1×, raise `--max-cell-size` (p=1/3 needs a 3× cell) or pick a nearby fraction.

**The CIF's occupancies fix the composition of every output.** enumlib shuffles *which* sites hold each species but keeps the counts fixed, so every enumerated structure has exactly the composition defined by your CIF. Two consequences:

- To preserve stoichiometry while introducing **antisite** disorder (e.g. Ge/P), make *both* the Ge site and a P site mixed (each {Ge:0.5, P:0.5}); the Ge and P totals then come out unchanged.
- **Li/vacancy disorder necessarily changes the Li count** (a Li site at 0.5 holds fewer Li). Studying Li disorder *at the stoichiometric composition* requires split Li positions (more Li sites than the ordered model, each partial) — an ordered MP CIF doesn't have these, so you need the disordered site model from the paper you're reproducing.

**Concrete LGPS example.** From an ordered Li10GeP2S12 CIF (P4_2mc, Z=2; Li0/Li1 at M=8, Li2/Li3/Ge4/P5/P6 at M=2):

- Ge/P antisite on the Ge4 and P6 sites (both → {Ge:0.5, P:0.5}) keeps the composition at Li10GeP2S12 and is integerizable in 1×.
- Setting Li0 to 0.5 gives 4 Li + 4 vacancy there (integerizable in 1×) but drops the composition to Li8GeP2S12 — use this only if your target study is Li-deficient.

Run with `xtalkit enumerate <edited>.cif --max-cell-size 1`. After editing, always recompute the composition and confirm it matches the paper you are reproducing before trusting the enumeration.

#### When enumlib returns 0 structures

The most common cause is that the parent CIF's occupancy **can't be integerized** in the supercell size you chose (see above). xtalkit surfaces this with a clear error. Two fixes:

1. **Raise `--max-cell-size`** so the supercell is large enough to hold an integer count of each species (e.g. occupancy 1/3 needs `--max-cell-size 3`).
2. **Prepare a "clean" parent CIF** with a rational occupancy close to the real value. Experimental CIFs often report awkward fractions (e.g. a Li site at occupancy 0.56 = 14/25). Rounding to a nearby simple fraction (0.56 → 0.5) gives a tractable parent that enumerates cleanly; everything else in the CIF (cell, other atoms, space group) is left unchanged. Rounding to a nearby rational stoichiometry is a standard, widely-used approximation when enumerating disordered crystals.

**Non-integer stoichiometry pre-flight check.** Before running `enum.x`, xtalkit checks that every species' count (`multiplicity × occupancy × cell_size`) is an integer for some `cell_size` in your `[--min-cell-size, --max-cell-size]` range. If not, it refuses with a table showing each offending species and the nearest valid fraction — e.g.:

```
[ERR] Cannot enumerate: non-integer stoichiometry (enumlib would likely run away and exhaust memory).
  species  mult  occ      mult*occ   nearest valid (at cell_size 1)
  Li      16    0.6910   11.0560    -> 0.6875
  Li      8     0.6430   5.1440     -> 0.6250
  ...
```

This guard exists because **`enumlib` does not fail cleanly on non-integer stoichiometry** — given the chance it tends to allocate memory until the kernel kills it (which can crash WSL/low-memory systems rather than producing an error). Rebuild the CIF with the rounded occupancies (`xtalkit build --atom-frac ...` using the suggested values) and re-run. `--skip-preflight` bypasses the check (dangerous — only if you are sure enumlib can handle your stoichiometry).

**Known limitations:**

- **Non-integer stoichiometry**: as above, an occupancy that isn't a simple fraction (e.g. 0.56) can't be integerized in a small supercell and yields 0 structures. Round to a nearby fraction or raise `--max-cell-size`.
- **Large multi-site enumerations can exhaust memory even with clean fractions.** A primitive cell with many atoms and several disordered sites (e.g. LGPS — P42/nmc, 50 atoms, Li on 16h+4d+8f plus Ge/P mixing) is a big search space: `enum.x` is single-threaded and can need multiple GB and many minutes. If it returns 0 / crashes: (a) set `--timeout` so runaway searches are killed cleanly, (b) **reduce scope** — fix the disorder you don't need by setting those sites to a single ordering (occupancy 1) in the parent CIF and enumerate one disordered site at a time, (c) keep `--max-cell-size 1`; `--max-structures` limits only the post-`enum.x` output stage.
- **Platform note**: `scripts/build_enumlib.sh` targets Linux/macOS (system `gfortran`). Windows users can compile under [WSL](https://learn.microsoft.com/windows/wsl/), or follow the legacy conda + `m2w64-gcc-fortran` path (place `enum.x`/`makestr.x` in the env's `Library/mingw-w64/bin`).


---

#### Enumeration setup

xtalkit's core (`mark`, `skeleton`, `info`, `fetch`, `build`) does **not** require pymatgen. Only `enumerate` does, and it's an opt-in uv extra — the core stays lightweight (gemmi + rich only). No conda, no root.

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
# Expect: all enumerator tests passed
```

On Windows, `xtalkit/_env.py` additionally applies three runtime workarounds (these do not run on Linux/macOS):

1. Appends `.X` and `.PY` to `PATHEXT` so `shutil.which("enum.x")` finds the binary
2. Calls `os.add_dll_directory(env/Library/bin)` so scipy's native extension loads
3. Monkey-patches `shutil.which` to return absolute paths (Windows otherwise returns `.\makestr.x`, which `subprocess.Popen` cannot launch)

---

### `shry` — Strict SHRY Enumeration Workflow

`xtalkit shry` is a **staged, audited** workflow for large partially occupied structures where enumlib can exhaust memory. It shells out to the external [SHRY](https://pypi.org/project/shry/) CLI (Pólya enumeration) and adds a count-before-generate safety gate plus per-stage manifests. Install it first: see [Tier 3 — +SHRY](#tier-3--shry-enumeration-xtalkit-shry) (`uv tool install shry` + `XTALKIT_SHRY_CMD`).

> **You need a primitive CIF for F-/I-centered cells.** Unlike `enumerate`, `shry` does **not** auto-primitivize — it enumerates the cell you give it. A conventional F-centered argyrodite cell (e.g. LPSC, Li1 on 48h @ 0.5 = 48 Li-sites) gives C(48,24) ≈ 10¹³ candidates and will hang. The primitive cell (12 Li-sites → C(12,6) = 924 → **48 inequivalent** under F-43m) is tractable. Supply a primitive CIF, or convert with spglib (`standardize_cell(to_primitive=True)`) preserving the space-group operations. xtalkit's `prepare` keeps the CIF's declared symmetry, so the primitive CIF must carry F-43m operations (not P1) — see [docs/user/shry-enumeration.md](docs/user/shry-enumeration.md).

The workflow has five stages — `prepare → count → enum → verify → postprocess`:

```bash
# 1) prepare: turn the partially occupied CIF into a SHRY-ready CIF.
#    Fills vacancy rows, checks occupancy sums to 1, checks integerizability,
#    audits parent symmetry, writes manifest + orbit-grouping sidecars.
xtalkit shry prepare LPSC_prim.cif \
  --out LPSC_shry_ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 216 \
  --target-formula Li6PS5Cl \
  --scaling-matrix 1 1 1

# 2) count: cheap safety gate. Counts symmetry-inequivalent configs (Pólya).
#    Inspect LPSC_shry_count.json -> "count_only_result" before generating.
xtalkit shry count LPSC_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 --angle-tolerance 5 --atol 1e-5 \
  --out LPSC_shry_count.json

# 3) enum: generate. Paste the count from step 2 into --expect-count.
#    Strict mode REFUSES to write if generation != count (catches misses).
xtalkit shry enum LPSC_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count 48 \
  --out LPSC_SHRY \
  --remove-vacancy X \
  --target-formula Li6PS5Cl \
  --write-cif --write-poscar --write-degeneracy

# 4) verify: count, formula, residual-vacancy, dedup, symprec-stability,
#    and (optionally) degeneracy-sum == raw combination count.
xtalkit shry verify LPSC_SHRY \
  --check-count --check-formula --check-dedup --check-degeneracy \
  --target-formula Li6PS5Cl \
  --symprec-list 1e-4 1e-3 1e-2

# 5) postprocess (optional): rank or generate DFT job folders.
xtalkit shry postprocess LPSC_SHRY --shortest-distance --pair Li Li
xtalkit shry postprocess LPSC_SHRY --write-tblite --write-slurm
```

**Output layout** (`LPSC_SHRY/`):

```
input/        # copied SHRY-ready CIF, command.txt, mod-only audit, pip freeze
raw_shry/     # SHRY's raw ordered CIFs (still contain X vacancy)
clean_cif/    # cleaned CIFs — X removed, element symbols, P1 supercell
poscar/       # one POSCAR per config (xtalkit-generated; SHRY has no POSCAR writer)
checks/       # verify.json, shortest_distance.json, ewald.json, symprec scans
manifest.json / manifest.jsonl   # full audit trail
```

**Subcommands:**

| Subcommand | Purpose |
|------------|---------|
| `shry prepare` | Build SHRY-ready CIF + audit (occupancy, integerizability, symmetry orbits) |
| `shry count` | Pólya count of inequivalent configs (no generation) |
| `shry enum` | Generate configs; strict mode requires `--expect-count` matching `count` |
| `shry verify` | Check count/formula/vacancy/dedup/symprec-stability/degeneracy |
| `shry postprocess` | Rank by shortest distance or Ewald energy; write tblite/CP2K folders + Slurm array |

**Strict mode.** All stages default to `--strict` (audited): `enum` requires `--expect-count`, runs a `--mod-only` chemistry audit, and refuses to write if generation ≠ count. Pass `--no-strict` to relax for exploratory runs (e.g. when you don't yet have a count).

If you start from refinement data, run `xtalkit build ... -o LPSC` first to make the partially-occupied CIF, then feed it to `shry prepare`. For a beginner-oriented walkthrough of every file and step, see [docs/user/shry-enumeration.md](docs/user/shry-enumeration.md); for the module design, [docs/dev/shry-module-design.md](docs/dev/shry-module-design.md).

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
║  [6] Build       — Build CIF from        ║
║                     refinement params    ║
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

  Output format: [1] cif  [2] xyz  [3] all  > 3

  Tolerance in fractional coords (default 0.5): [Enter for default]

  Element override (e.g. '4a:Xe,16e:Kr') or Enter to skip: [Enter]

  Output base path [default: D:\structures\Li6PS5Cl_WYCK]: [Enter]

  ✓ Saved to: D:\structures\Li6PS5Cl_WYCK.cif,
               D:\structures\Li6PS5Cl_WYCK.xyz

  Press Enter to continue...
```

---

## Reference

### Dummy Atom System

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

### Overlay vs. Replace Modes

#### Overlay (default)

```
Before:  Li at (0,0,0)   P at (0.25,0.25,0.25)
After:   Li at (0,0,0)   P at (0.25,0.25,0.25)
         WYCK_4a (Xe) at (0,0,0)
         WYCK_4c (Kr) at (0.25,0.25,0.25)
         ... more dummy atoms for other requested Wyckoff positions
```

Real atoms are preserved. Dummy atoms are added on top. In VESTA, you see both.

#### Replace

```
Before:  Li at (0,0,0)   P at (0.25,0.25,0.25)
After:   WYCK_4a (Xe) at (0,0,0)    ← Li replaced (matched 4a)
         P at (0.25,0.25,0.25)      ← P stayed (didn't match 4c within tolerance)
```

Atoms that match a requested Wyckoff position (within tolerance) are **replaced** by dummy atoms. Other atoms are left untouched.

---

### Matching Tolerance

The tolerance (`--tol`, default 0.5) controls how closely an atom's coordinates must match a Wyckoff position's theoretical coordinates to be considered occupying that position.

- **Larger values** (0.5–1.0): lenient matching, suitable for experimental structures with slight coordinate deviations
- **Smaller values** (0.01–0.1): strict matching, only atoms very close to ideal Wyckoff positions are considered

The tolerance is applied in **fractional coordinate space** (not angstroms). For cubic cells with roughly equal axes, 0.5 in fractional space ≈ 0.5 × a in real distance along each axis.

---

### Workflow Recipes

#### Recipe 1: Study Li₆PS₅Cl (SG 216, F-43m) Wyckoff occupancy

```bash
# Step 1: Download CIF from Materials Project
# (you already have the file)

# Step 2: See what Wyckoff positions exist in F-43m
xtalkit info --sg 216

# Step 3: Mark ALL Wyckoff positions, overlay mode, all formats
xtalkit mark Li6PS5Cl.cif --sg 216 --wyckoff all --format cif,xyz

# Step 4: Open Li6PS5Cl_WYCK.cif in VESTA
# → All 8 Wyckoff positions are now visible as colored dummy atoms
# → Real atoms (Li, P, S, Cl) are still shown
# → You can toggle atoms in VESTA to compare
```

#### Recipe 2: Create a Wyckoff reference skeleton

```bash
# Generate a skeleton for F-43m with real cell parameters
xtalkit skeleton --sg 216 --wyckoff all \
    --cell "9.85 9.85 9.85 90 90 90" \
    --format cif

# Open SG216_skeleton.cif in VESTA
# → See exactly where each Wyckoff position sits in the unit cell
# → No real atoms — pure reference template
```

#### Recipe 3: Check which atoms occupy specific Wyckoff positions

```bash
# Mark only the Wyckoff positions you care about
xtalkit mark structure.cif --sg 225 --wyckoff 4a,8c --mode replace

# In replace mode, atoms at 4a and 8c are swapped with dummy atoms
# → Instantly see in VESTA: "are there atoms at these positions?"
```

#### Recipe 4: Batch process multiple structures

```bash
# All .cif files in a directory, same space group
for f in *.cif; do
    xtalkit mark "$f" --sg 216 --wyckoff all -o "${f%.cif}_WYCK"
done
```

---

### Supported Space Groups

Wyckoff position data is available for **all 230 space groups**. The bundled dataset (`xtalkit/data/wyckoff.json`) is derived from International Tables for Crystallography Vol. A via [pyxtal](https://pyxtal.readthedocs.io) (MIT) and verified against gemmi's symmetry operations — see `scripts/build_wyckoff_db.py` to regenerate it.

`xtalkit fetch` confirms the dataset is intact (230/230).

---

## Development

```bash
uv sync                     # Install core dependencies
uv sync --extra enumerate   # Also enable `enumerate` (pulls pymatgen)
uv run pytest               # Run all tests
```

To run the SHRY tests against the real backend (not just the fake), install SHRY per [Tier 3](#tier-3--shry-enumeration-xtalkit-shry) and export `XTALKIT_SHRY_CMD`; otherwise the SHRY tests use an in-process fake backend.

### Project structure

```
xtalkit/
├── xtalkit/
│   ├── __init__.py          # Package, version
│   ├── cli.py               # argparse CLI + 7 subcommands
│   ├── tui.py               # rich-based interactive TUI
│   ├── spacegroup.py        # Gemmi space group queries
│   ├── matcher.py           # Atom → Wyckoff position matching
│   ├── marker.py            # Core: mark Wyckoff in CIF
│   ├── skeleton.py          # Pure Wyckoff skeleton generation
│   ├── builder.py           # Build CIF from refinement params
│   ├── exporter.py          # .cif / .xyz writers
│   ├── enumerator.py        # enumlib wrapper (lazy pymatgen import)
│   ├── enumeration/         # SHRY workflow package
│   │   ├── shry_prepare.py / shry_count.py / shry_enum.py / shry_verify.py
│   │   ├── shry_backend.py  # subprocess wrapper for the `shry` CLI
│   │   ├── cif_io.py        # esd-aware CIF reader/writer (gemmi)
│   │   ├── occupancy.py / formula.py / degeneracy.py / fingerprint.py
│   │   └── postprocess.py   # ranking + tblite/CP2K/Slurm job generation
│   ├── _env.py              # enumlib binary discovery + Windows env fixes
│   └── utils.py             # Shared helpers
├── tests/
│   ├── fixtures/            # simple.cif, disordered_binary.cif
│   ├── test_marker.py / test_skeleton.py / test_builder.py / test_matcher.py
│   ├── test_spacegroup.py / test_exporter.py / test_cli.py / test_tui.py
│   ├── test_enumerator.py   # enumlib integration (skips without pymatgen)
│   └── test_shry_module.py  # SHRY workflow (uses a fake backend by default)
├── docs/
│   ├── user/shry-enumeration.md   # beginner SHRY walkthrough
│   └── dev/shry-module-design.md  # SHRY module design & bug history
├── scripts/
│   ├── build_enumlib.sh     # Compile enumlib (enum.x, makestr.x) from source
│   └── build_wyckoff_db.py  # Regenerate the bundled Wyckoff dataset
├── LPSC.cif                 # Example: argyrodite Li6PS5Cl (F-43m, Li partial occ.)
├── pyproject.toml
└── README.md
```

---

## Dependencies

| Package | Purpose | Required? |
|---------|---------|-----------|
| [gemmi](https://gemmi.readthedocs.io/) | Space group data, CIF I/O | Yes (core) |
| [rich](https://rich.readthedocs.io/) | TUI formatting (tables, panels, colors) | Yes (core) |
| [pymatgen](https://pymatgen.org/) ≥2024.5 | enumlib wrapper for `enumerate` | `enumerate` extra only |
| [enumlib](https://github.com/msg-byu/enumlib) | Symmetry-inequivalent enumeration (Fortran, source-compiled) | Only for `enumerate` |
| [SHRY](https://pypi.org/project/shry/) 1.1.x | Pólya enumeration for large partial-occupancy cells (isolated install) | Only for `shry` |
| pytest (dev) | Test framework | Yes (dev) |

`enumerate` lazy-imports pymatgen and `shry` shells out to the external `shry` CLI, so the core toolkit works without either. SHRY is **not** declared as a pip dependency because it pins `pymatgen≤2023.10.4`, which conflicts with the `enumerate` extra — install it isolated via `uv tool install shry` (see [Tier 3](#tier-3--shry-enumeration-xtalkit-shry)).
