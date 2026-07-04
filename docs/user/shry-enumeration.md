# SHRY Enumeration

`xtalkit shry` is a strict workflow for large partially occupied structures
where enumlib can become too memory-intensive. It does not replace the existing
`xtalkit enumerate` command; both backends remain available.

The workflow is intentionally staged:

```bash
xtalkit shry prepare input.cif \
  --out ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 137 \
  --target-formula Li20Ge2P4S24

xtalkit shry count ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 --angle-tolerance 5 --atol 1e-5 \
  --out count.json

xtalkit shry enum ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count <COUNT> \
  --out shry_enum \
  --remove-vacancy X \
  --target-formula Li20Ge2P4S24 \
  --write-cif --write-poscar --write-degeneracy

xtalkit shry verify shry_enum \
  --check-count --check-formula --check-dedup \
  --target-formula Li20Ge2P4S24 \
  --symprec-list 1e-4 1e-3 1e-2
```

`prepare` writes explicit `X` vacancy rows for partially occupied sites, checks
occupancy sums with exact fractions, checks that occupancies are integerizable
for the requested scaling matrix, preserves symmetry operations, and writes
sidecar manifest/orbit-grouping JSON files.

`count` and `enum` call the external `shry` executable through subprocess. If
SHRY is installed in a separate virtual environment, set:

```bash
export XTALKIT_SHRY_CMD=/path/to/shry
```

Strict mode requires `enum --expect-count`; this prevents long generation runs
whose output count cannot be audited. `enum` runs SHRY `--mod-only` first, then
streams generated CIF files from `raw_shry/`, removes vacancy species, checks
the target formula, and writes `manifest.jsonl`.

For LGPS-scale systems, run `count` first. Do not start exhaustive generation
until the count and parent-spacegroup/orbit audit are reasonable.

Optional independent checks and post-processing:

```bash
# Independent backend count, if the external supercell program is configured.
export XTALKIT_SUPERCELL_CMD=/path/to/supercell
xtalkit shry verify shry_enum --check-count --cross-backend supercell

# Geometry ranking and job-folder generation.
xtalkit shry postprocess shry_enum --shortest-distance --pair Li Li
xtalkit shry postprocess shry_enum --write-tblite --write-slurm
xtalkit shry postprocess shry_enum --write-cp2k --cp2k-template cp2k.inp

# Ewald ranking requires explicit oxidation states and pymatgen.
xtalkit shry postprocess shry_enum --ewald \
  --ewald-charges Li:1 Ge:4 P:5 S:-2
```

`--cross-backend supercell` and `--ewald` are explicit opt-in actions. They do
not run during normal verification because they require external executables or
additional scientific assumptions.
