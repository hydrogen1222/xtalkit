# xtalkit: crystal Wyckoff and occupational-disorder tools

[中文手册](README.zh-CN.md) · [Detailed SHRY guide](docs/user/shry-enumeration.md)

xtalkit handles a small set of recurring crystallographic tasks:

- inspect space groups and Wyckoff positions;
- mark positions in an existing CIF or draw a Wyckoff skeleton;
- build a CIF from refinement tables;
- enumerate ordered configurations from partial occupancies;
- rank comparable ordered structures with an Ewald electrostatic model.

It is not a refinement program, a structure relaxer, or a substitute for
physical choices such as oxidation states, spin, supercell size, and DFT
settings. When those choices matter, xtalkit asks for them or stops with an
error.

## Start here

Install the core application from the repository:

```bash
git clone https://github.com/hydrogen1222/xtalkit.git
cd xtalkit
uv sync
uv run xtalkit --version
```

Try a complete operation with a bundled test CIF:

```bash
uv run xtalkit info --sg 216
uv run xtalkit mark tests/fixtures/simple.cif \
  --sg 216 --wyckoff 4a,4c -o quickstart_marked
```

Open `quickstart_marked.cif` in VESTA. Marked and skeleton files carry a
`visualization_only` tag: use them to inspect sites, not as physical structures
for energy or composition calculations.

Run without a subcommand for the terminal interface:

```bash
uv run xtalkit
```

In an interactive terminal this opens the guided interface. With redirected
I/O or in CI, it prints help instead of waiting for input. The interface covers
the common space-group, marking, building, skeleton, enumlib, and Ewald tasks.
SHRY remains a command-line workflow because its count and manifest should be
reviewed between stages.

## Choose a command

| Goal | Command | Result |
|---|---|---|
| Inspect a space group | `info` | Wyckoff table |
| Show selected sites in a CIF | `mark` | visualization CIF/XYZ |
| Draw positions without a real structure | `skeleton` | visualization CIF/XYZ |
| Build a structure from refinement data | `build` | physical CIF/XYZ |
| Enumerate small or medium disorder | `enumerate` | ordered structures via enumlib |
| Count and enumerate larger disorder | `shry` | staged SHRY run and manifests |
| Rank ordered structures electrostatically | `ewald` | ranking CSV and optional groups |
| Check the bundled database | `fetch` | 230-space-group integrity result |

The CLI help is the option reference:

```bash
uv run xtalkit --help
uv run xtalkit build --help
uv run xtalkit shry enum --help
```

## Install only what you need

The core environment contains gemmi and Rich. It supports `info`, `mark`,
`skeleton`, `build`, `fetch`, and the terminal interface.

### enumlib and Ewald

Ewald scoring needs pymatgen/spglib. Enumeration additionally needs the
Fortran enumlib executables:

```bash
uv sync --extra enumerate
bash scripts/build_enumlib.sh
```

The build requires `gfortran`, `make`, and `git`. On Debian or Ubuntu:

```bash
sudo apt install gfortran make git
```

The script uses one compiler process by default and installs into
`~/.local/share/xtalkit/bin/`; xtalkit discovers that directory without a PATH
change. Its default source revision is pinned to the audited enumlib commit:

```text
67217859c1cfe6e1d8cfe25b76a1d14a28c6030e
```

Verify the installation with a small real enumeration:

```bash
cat ~/.local/share/xtalkit/bin/enumlib-commit.txt
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --max-cell-size 2 --output-dir enum_smoke
# Expected: three ordered Au/Cu CIFs
```

Advanced users may deliberately override the source or build concurrency with
`ENUMLIB_REF=<commit-or-tag>` and `XTALKIT_BUILD_JOBS=<positive-integer>`.

### SHRY

SHRY 1.1.8 pins an older pymatgen, so keep it outside the xtalkit environment:

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
uv tool list
```

Do not use `uv pip install shry` in this project. If the isolated executable is
not already discoverable, set:

```bash
export XTALKIT_SHRY_CMD="$(command -v shry)"
```

The strict in-process symmetry and duplicate checks still need
`uv sync --extra enumerate`; they do not require compiled enumlib.

### tblite (optional)

xtalkit writes tblite job folders but does not bundle the executable. The
official binary distribution is available through
[conda-forge](https://tblite.readthedocs.io/en/stable/installation.html). An
isolated micromamba prefix avoids mixing it with either Python environment:

```bash
# Install micromamba first if needed:
# https://mamba.readthedocs.io/en/stable/installation/micromamba-installation.html
MAMBA_ROOT_PREFIX="$HOME/.local/share/micromamba" \
  micromamba create --yes \
  --prefix "$HOME/.local/share/xtalkit/tblite-0.7.0" \
  --channel conda-forge --strict-channel-priority 'tblite=0.7.0'
mkdir -p "$HOME/.local/bin"
ln -sfn "$HOME/.local/share/xtalkit/tblite-0.7.0/bin/tblite" \
  "$HOME/.local/bin/tblite"
tblite --version
```

Put `~/.local/bin` on `PATH`. Before running a generated job, set charge and
unpaired-electron count for the actual system; neutral singlet defaults are not
universally valid.

## Common workflows

### Inspect a space group

```bash
uv run xtalkit info --sg 225
```

The default cell in this table is only a crystal-system-compatible drawing
cell. It is not an experimental lattice for a material.

### Mark positions in a CIF

```bash
uv run xtalkit mark input.cif \
  --sg 216 --wyckoff 4a,16e,48h \
  --mode overlay --format cif -o input_marked
```

`overlay` keeps the real atoms and adds visual markers. `replace` changes the
element and label of matched atoms while retaining their refined coordinates
and occupancies. Both modes are for visualization.

Useful options:

- `--wyckoff all` selects every position in the group;
- `--tol 0.01` controls the fractional-coordinate constraint residual;
- `--map 4a:Xe,16e:Kr` chooses marker elements;
- `--format cif,xyz` writes both formats;
- nonzero `--offset` is a visual displacement away from the exact site.

The requested `--sg` must agree with the group declared by the CIF. Resolve a
setting mismatch instead of changing the number merely to bypass the check.

### Draw a Wyckoff skeleton

```bash
uv run xtalkit skeleton \
  --sg 216 --wyckoff all \
  --cell "9.86 9.86 9.86 90 90 90" \
  -o SG216_skeleton
```

A skeleton contains marker species only. If `--cell` is omitted, xtalkit uses
a drawing cell compatible with the crystal system.

### Build a CIF from refinement data

Use Wyckoff labels and free coordinates when the table supplies them:

```bash
uv run xtalkit build \
  --sg 225 \
  --cell "5.64 5.64 5.64 90 90 90" \
  --atom "Na 4a" --atom "Cl 4b" \
  -o NaCl
```

For a `16e` template `(x,x,x)`, `"Li 16e 0.30"` supplies `x`; use
`"Li 16e 0.30 0.6875"` to append an occupancy.

When a refinement table gives fractional coordinates directly, prefer
`--atom-frac`:

```bash
uv run xtalkit build \
  --sg 137 \
  --cell "8.69407 8.69407 12.5994 90 90 90" \
  --atom-frac "Li 0.2563 0.2718 0.1832 0.6875" \
  --atom-frac "P 0 0 0.5 1" \
  -o refinement_model
```

The detected Wyckoff labels are printed so that coordinate or setting errors
are visible. Multiple species may share one disordered site, but their
occupancies must not sum above one.

For repeatable jobs, use a JSON spec:

```json
{
  "sg": 225,
  "cell": {
    "a": 5.64, "b": 5.64, "c": 5.64,
    "alpha": 90, "beta": 90, "gamma": 90
  },
  "atoms": [
    {"element": "Na", "wyckoff": "4a"},
    {"element": "Cl", "wyckoff": "4b"}
  ]
}
```

```bash
uv run xtalkit build --spec nacl.json -o NaCl
```

`--spec` is an alternative to all of `--sg`, `--cell`, `--atom`, and
`--atom-frac`.

### Enumerate disorder with enumlib

```bash
uv run xtalkit enumerate tests/fixtures/disordered_binary.cif \
  --min-cell-size 1 --max-cell-size 2 \
  --output-dir AuCu_enum
```

Cell size is a volume multiplier relative to the primitive cell, not a repeat
count along each axis. xtalkit converts to a primitive structure before calling
enumlib.

Start unfamiliar systems conservatively:

```bash
uv run xtalkit enumerate disordered.cif \
  --min-cell-size 1 --max-cell-size 1 \
  --timeout 10 --max-structures 20 \
  --output-dir trial_enum
```

`--max-structures` limits output generation after `enum.x` has completed; it
does not cap the enumlib search. `--jobs` parallelizes only result generation.
The default of one worker and batched output has the lowest memory risk.

Before launching enumlib, xtalkit checks whether occupancies can become integer
atom counts in the requested cell-size range. Treat `--skip-preflight` as a
dangerous expert escape hatch: invalid inputs can make enumlib consume extreme
amounts of memory.

### Enumerate larger disorder with SHRY

The SHRY workflow deliberately has review points:

```text
prepare → count → enum → verify → optional postprocess
```

```bash
uv run xtalkit shry prepare input.cif \
  --out ready.cif --parent-spacegroup 216 \
  --target-formula Li6PS5Cl --scaling-matrix 1 1 1

uv run xtalkit shry count ready.cif \
  --scaling-matrix 1 1 1 --out count.json

# Read count_only_result from count.json before continuing.
uv run xtalkit shry enum ready.cif \
  --scaling-matrix 1 1 1 --expect-count <COUNT> \
  --target-formula Li6PS5Cl --write-cif --write-poscar \
  --out SHRY_run

uv run xtalkit shry verify SHRY_run --target-formula Li6PS5Cl
```

Large F- or I-centred conventional cells can create astronomical candidate
counts. Read the [detailed SHRY guide](docs/user/shry-enumeration.md) before a
production run; it covers primitive cells, manifests, failure recovery,
degeneracy, and tblite post-processing.

### Rank structures with Ewald energy

Ewald scoring requires fully ordered, charge-neutral periodic structures with
known oxidation states:

```bash
uv run xtalkit ewald ordered_cifs/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --sort asc --out ewald_results/ranking.csv
```

Copy the lowest 20 into a selected group:

```bash
uv run xtalkit ewald ordered_cifs/ \
  --charges Li:1 P:5 S:-2 Cl:-1 \
  --top-n 20 --group --group-dir ewald_results
```

Grouping copies by default; `--move` removes originals after transfer. A
single ranking must contain one composition. Total energies additionally need
the same atom count. Use `--per-atom` only for commensurate cells of that same
composition. Ewald model energies are not relaxed DFT total energies.

## Output and recovery

- `mark` and `skeleton` write visualization-only files;
- `build` writes physical CIFs with symmetry, occupancy, and unit-cell formula;
- `enumerate` writes one ordered structure per file under `<input>_enum/` by
  default;
- `shry enum` keeps copied inputs, raw output, cleaned CIFs, and a run manifest;
- `ewald` writes the complete ranking to CSV and optionally copies/moves files.

Single CIF/XYZ outputs may overwrite a same-named file. Give important runs an
explicit `-o` or a fresh directory. SHRY output directories must be new or
empty so that separate runs cannot be mixed.

## Troubleshooting

`VIRTUAL_ENV=... does not match .../.venv`
: Deactivate the unrelated environment, then use `uv run`. Use uv's `--active`
  only when targeting the already-active environment intentionally.

`enum.x` or `makestr.x` not found
: Re-run `bash scripts/build_enumlib.sh`. xtalkit discovers its user-local
  install automatically.

`non-integer stoichiometry`
: The occupancy, orbit multiplicity, and requested cell size cannot produce
  integer atom counts. Increase `--max-cell-size` or correct the occupancy from
  the source data; do not bypass the guard merely to make the command run.

Space-group mismatch
: Check the CIF symbol, International Tables number, setting, and symmetry
  operations. A wrong group changes orbits, formulas, and enumeration counts.

SHRY count is unexpectedly huge
: Return to the count stage. Check for an F/I conventional cell, the scaling
  matrix, and the target formula before generation.

Ewald reports a net charge
: Correct the oxidation states or structure composition. xtalkit does not add a
  compensating background charge.

The terminal interface re-prompts invalid critical numbers rather than silently
falling back to an unsafe default. Press `Ctrl+C` to leave a workflow; already
written files are not deleted.

## Development

```bash
uv sync --extra enumerate
uv run pytest
uvx ruff check xtalkit tests
uv build
```

Tests cover all 230 space groups and the crystallographic validation paths.
When real SHRY, enumlib, or tblite executables are available, their integration
smoke tests run automatically; otherwise those tests are explicitly skipped.
