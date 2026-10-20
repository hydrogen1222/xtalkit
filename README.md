# xtalkit

[中文说明](README.zh-CN.md) · [SHRY guide](docs/user/shry-enumeration.md)

**In one line: turn a refinement table into a CIF, then turn a partially occupied CIF into a set of ordered structures you can actually compute.**

```text
refinement table ──build──▶ disordered CIF ──enumerate / shry──▶ ordered structures ──ewald──▶ low-energy picks ──▶ DFT
```

There are also a few small helpers for looking at Wyckoff positions. xtalkit does not refine,
relax, or choose oxidation states and supercells for you; when a physical choice is needed it
asks for it or stops with an error instead of guessing.

## Commands at a glance

| I want to… | Command |
|---|---|
| Build a CIF from refinement data | `build` |
| Enumerate ordered structures from partial occupancy (small/medium) | `enumerate` |
| Same, but the count is huge, so count first, then generate | `shry` |
| Rank a set of ordered structures by electrostatic energy | `ewald` |
| List the Wyckoff positions of a space group | `info` |
| Highlight Wyckoff positions in an existing CIF (viewing only) | `mark` |
| Draw a space group's positions without a real structure | `skeleton` |
| Check the bundled 230-space-group data | `fetch` |

Every command has `--help`, for example `uv run xtalkit build --help`.
Running `uv run xtalkit` with no subcommand opens an interactive interface.

## Install

```bash
git clone https://github.com/hydrogen1222/xtalkit.git
cd xtalkit
uv sync
uv run xtalkit --version
```

That is enough for `build`, `info`, `mark`, `skeleton`, and `fetch`.
`enumerate`, `shry`, and `ewald` need extra components; see [Optional components](#optional-components).

---

## Main workflow: model from refinement data

The example is Li₆PS₅Cl (argyrodite, F-43m, No. 216), which has half-occupied Li and mixed Cl/S sites.

### Step 1: copy the refinement table into JSON

A JSON file is easier to edit, review, and keep with a project than a long command line.
A ready-to-copy template is in [examples/argyrodite_template.json](examples/argyrodite_template.json).

```json
{
  "sg": 216,
  "cell": {"a": 9.86, "b": 9.86, "c": 9.86, "alpha": 90, "beta": 90, "gamma": 90},
  "atoms": [
    {"element": "Li", "xyz": [0.182, 0.182, 0.022], "occ": 0.5},
    {"element": "P",  "xyz": [0.5,   0.5,   0.5]},
    {"element": "S",  "xyz": [0.62,  0.62,  0.62]},
    {"element": "Cl", "xyz": [0,     0,     0],    "occ": 0.5},
    {"element": "S",  "xyz": [0,     0,     0],    "occ": 0.5},
    {"element": "S",  "xyz": [0.75,  0.75,  0.75], "occ": 0.5},
    {"element": "Cl", "xyz": [0.75,  0.75,  0.75], "occ": 0.5}
  ]
}
```

> The numbers are illustrative. Replace them with the refined values from the paper you cite.

| Refinement table column | JSON | Notes |
|---|---|---|
| Space group | `"sg"` | Number only. The paper's setting (origin, axes) must be the standard one |
| a, b, c, α, β, γ | `"cell"` | All six, even for cubic |
| Atom / element | `"element"` | Element symbol such as `Li`, not `Li1` or `Li+` |
| x, y, z | `"xyz"` | Copy as printed; any representative works and the Wyckoff orbit is detected |
| Occ / SOF | `"occ"` | 0–1, default 1 |
| Wyckoff | not needed | Detected and printed so you can compare with the paper |
| Uiso / Biso | not needed | Not written to the CIF |

**Mixed sites**: write two rows with the same coordinates (Cl 0.5 + S 0.5 on 4a above).
The occupancies on one site must not sum above 1.

If the paper gives a Wyckoff label and free parameters instead, write
`{"element": "S", "wyckoff": "16e", "free": [0.62]}`. Both forms can be mixed in one file.

### Step 2: run

```bash
uv run xtalkit build --spec examples/argyrodite_template.json -o LPSC
```

```text
[info] SG #216 (F-43m), cubic — detected Wyckoff sites:
        Li  at (0.1820,0.1820,0.0220) occ=0.5   -> 48h (mult 48)
        P   at (0.5000,0.5000,0.5000) occ=1.0   -> 4b (mult 4)
        S   at (0.6200,0.6200,0.6200) occ=1.0   -> 16e (mult 16)
        ...
[OK] SG #216 (F-43m), crystal system: cubic, formula: Li6PS5Cl
     Saved to: LPSC.cif
```

### Step 3: check (do not skip)

1. **Wyckoff labels match the paper.** A different label or `?` usually means a typo or a non-standard setting.
2. **The formula is right.** It is multiplicity × occupancy. Charge neutrality is **not** checked.
   Elements appear in the order of their first row, so list Li first to get `Li6PS5Cl`.
3. Open `LPSC.cif` in VESTA.

For a quick try without a JSON file:

```bash
uv run xtalkit build --sg 216 --cell "9.86 9.86 9.86 90 90 90" \
  --atom-frac "P 0.5 0.5 0.5" \
  --atom-frac "Cl 0 0 0 0.5" --atom-frac "S 0 0 0 0.5" \
  -o test
# --atom-frac "element x y z [occ]"; --atom "element wyckoff [free...] [occ]"; both can be combined
```

### Common build errors

| Message contains | Cause and fix |
|---|---|
| `total occupancy ... exceeds 1.0` | Site occupancies sum above 1, often from rounding (0.385 + 0.616). Renormalize |
| `unit cell is inconsistent` | Cell does not fit the crystal system (e.g. cubic with a ≠ c) |
| `Invalid Wyckoff label` | Wyckoff form only: use the full label such as `16e`, not `e` |
| `Unknown element` | Symbols like `Li1`, `li`, `Li+` |
| A site detected as `?` | The coordinates match no Wyckoff position: check for typos or a wrong space group/setting |

---

## Next: from disorder to ordered structures

`build` gives the average structure with fractional occupancy. DFT needs a supercell in which every site is either occupied or empty.

**Small or medium disorder: `enumerate`** (needs enumlib):

```bash
uv run xtalkit enumerate LPSC.cif \
  --min-cell-size 1 --max-cell-size 1 \
  --timeout 10 --max-structures 20 --output-dir LPSC_enum
```

- `--min/max-cell-size` are volume multiples of the **primitive** cell, not per-axis factors.
- If occupancies cannot become integer atom counts at that size, xtalkit refuses before starting
  enumlib and suggests the nearest feasible occupancy. Avoid `--skip-preflight`: bad input can make enumlib exhaust memory.

**Large counts: `shry`**, in reviewed stages `prepare → count → enum → verify`:

```bash
uv run xtalkit shry prepare LPSC.cif --out ready.cif \
  --parent-spacegroup 216 --target-formula Li6PS5Cl --scaling-matrix 1 1 1
uv run xtalkit shry count ready.cif --scaling-matrix 1 1 1 --out count.json
# read count_only_result in count.json, then:
uv run xtalkit shry enum ready.cif --scaling-matrix 1 1 1 --expect-count <COUNT> \
  --target-formula Li6PS5Cl --write-cif --write-poscar --out SHRY_run
uv run xtalkit shry verify SHRY_run --target-formula Li6PS5Cl
```

The optional postprocess stage and recovery steps are in the [SHRY guide](docs/user/shry-enumeration.md).

> **Argyrodite caveat**: 48h Li sites come in pairs less than 1 Å apart, which cannot both hold Li.
> The enumerators do not know this and will produce such configurations. Filter by the shortest
> Li–Li distance after enumeration.

## Finally: rank by Ewald energy

```bash
uv run xtalkit ewald LPSC_enum/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --top-n 20 --group --group-dir ewald_results
```

- The full ranking goes to a CSV; `--group` **copies** the top N into `selected/` (`--move` moves instead).
- Only rank structures with the same composition and atom count; they must be charge neutral.
- Ewald energy is a screen, not a DFT energy.

---

## Other helpers

```bash
uv run xtalkit info --sg 216                                  # Wyckoff table
uv run xtalkit mark input.cif --sg 216 --wyckoff 4a,48h -o m  # highlight sites in a CIF
uv run xtalkit skeleton --sg 216 --wyckoff all \
  --cell "9.86 9.86 9.86 90 90 90" -o skel                    # draw positions
uv run xtalkit fetch                                          # self-check the database
```

`mark` and `skeleton` outputs carry a `visualization_only` tag and are **not** for calculations.
`mark` checks `--sg` against the space group declared in the CIF; if they differ, find out
which setting the CIF uses instead of changing the number until it passes.

## Optional components

### enumlib (for `enumerate`) and pymatgen (for `ewald`)

```bash
sudo apt install gfortran make git
uv sync --extra enumerate
bash scripts/build_enumlib.sh           # installs to ~/.local/share/xtalkit/bin/, no PATH change
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --max-cell-size 2 --output-dir enum_smoke     # expect three ordered CIFs
```

The script pins the audited enumlib commit (`67217859c1cf…`) and compiles with one job.
Override deliberately with `ENUMLIB_REF=<commit>` and `XTALKIT_BUILD_JOBS=<N>`.

### SHRY (for `shry`)

SHRY 1.1.8 pins an older pymatgen, so install it **separately**, never with `uv pip install shry`:

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
# if not found on a cluster: export XTALKIT_SHRY_CMD="$(command -v shry)"
```

### tblite (optional)

xtalkit writes tblite job folders but does not bundle the executable. Install it in an isolated
conda-forge prefix ([installation](https://tblite.readthedocs.io/en/stable/installation.html)):

```bash
MAMBA_ROOT_PREFIX="$HOME/.local/share/micromamba" \
  micromamba create --yes --prefix "$HOME/.local/share/xtalkit/tblite-0.7.0" \
  --channel conda-forge --strict-channel-priority 'tblite=0.7.0'
ln -sfn "$HOME/.local/share/xtalkit/tblite-0.7.0/bin/tblite" "$HOME/.local/bin/tblite"
```

Check total charge and spin for your system; `--charge 0 --spin 0` is not always right.

## Output files

| Command | Output |
|---|---|
| `build` | Physical CIF with symmetry operations, occupancy and formula; optional XYZ |
| `enumerate` | `<input>_enum/` by default, one ordered configuration per file |
| `shry enum` | `clean_cif/` (vacancies removed), `manifest.json`, `raw_shry/` |
| `ewald` | Ranking CSV; with `--group`, also `selected/` and `rest/` |
| `mark`, `skeleton` | Visualization-only CIF/XYZ |

Files with the same name are overwritten; pass `-o` or a new output directory for important runs.
A `shry enum` output directory must be new or empty.

## Troubleshooting

- **`VIRTUAL_ENV=... does not match`**: another environment is active; run `deactivate`, then `uv run`.
- **`enum.x` not found**: rerun `bash scripts/build_enumlib.sh`.
- **`non-integer stoichiometry`**: occupancies do not give whole atoms at that cell size; raise `--max-cell-size` or check the occupancies.
- **Space group mismatch**: check that the CIF's symbol, number, setting and operations agree. xtalkit will not guess.
- **SHRY count is huge**: go back to `shry count`; check for an F/I-centred conventional cell or a wrong scaling matrix.
- **Ewald says not neutral**: check oxidation states and the formula. No compensating background is added.

## Development

```bash
uv sync --extra enumerate
uv run pytest
uvx ruff check xtalkit tests
```
