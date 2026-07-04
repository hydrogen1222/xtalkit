# SHRY Enumeration: Beginner Guide

This page explains how to replace the old enumlib command with the safer SHRY
workflow. It is written for the LGPS case in `LGPS.txt`, but the same pattern
works for other partially occupied CIF files.

## 0. What changes compared with enumlib?

Your old workflow had two commands:

```bash
# Build the average/disordered CIF from refinement data.
xtalkit build ... -o LGPS

# Old enumlib route.
xtalkit enumerate LGPS.cif --max-cell-size 1 --output-dir ./LGPS
```

With SHRY, keep the first `xtalkit build ... -o LGPS` command. Replace only the
second command with these stages:

1. `prepare`: make a SHRY-ready CIF and audit the parent symmetry.
2. `count`: count symmetry-inequivalent structures before generating them.
3. `enum`: generate structures only after you have a count.
4. `verify`: check count, formula, vacancies, duplicates, and symmetry sensitivity.
5. `postprocess`: optional ranking or job-file generation.

This is deliberately more verbose than one enumlib command because each step
creates an audit file you can inspect before spending time and memory.

## 1. Install or point to SHRY

`xtalkit shry` calls the external `shry` command. SHRY is a pure-Python
[PyPI package](https://pypi.org/project/shry/) (no compilation), but it pins
`pymatgen<=2023.10.4`, which **conflicts** with xtalkit's `enumerate` extra
(`pymatgen>=2024.5`). So install SHRY in an *isolated* environment, not the
xtalkit venv:

```bash
uv tool install shry                       # isolated env; `shry` goes on PATH
export XTALKIT_SHRY_CMD="$(which shry)"    # point xtalkit at it
shry --version                             # verify (e.g. SHRY 1.1.8)
```

Add the `export` line to your shell profile so it persists. If `shry` is
already on `PATH` (e.g. from a conda env you activate manually), the
`XTALKIT_SHRY_CMD` export is optional — xtalkit falls back to `shutil.which("shry")`.

> Do **not** `uv pip install shry` into the xtalkit venv: it downgrades
> pymatgen to 2023.10.4 and breaks `uv sync --extra enumerate`.

## 1b. Provide a primitive CIF for F-/I-centered cells

Unlike `xtalkit enumerate`, **`shry` does not auto-primitivize** — it
enumerates the cell you give it. A conventional F-centered argyrodite cell
(LPSC, Li1 on 48h @ 0.5 = 48 Li-sites) gives C(48,24) ≈ 10¹³ candidates and
will hang. The primitive cell (12 Li-sites → C(12,6) = 924 → **48
inequivalent** under F-43m) is tractable.

So before `prepare`, supply a **primitive** CIF that carries the F-43m
symmetry operations (not P1 — `prepare` keeps the CIF's declared symmetry,
and strict mode refuses P1). Convert with spglib:

```python
# Save as make_primitive.py
import spglib, numpy as np
from pymatgen.core import Structure
s = Structure.from_file("LPSC.cif")
nums = [list(site.species.elements)[0].Z for site in s]
cell = (s.lattice.matrix, s.frac_coords, np.array(nums))
lat, pos, num = spglib.standardize_cell(cell, to_primitive=True)
sym = spglib.get_symmetry((lat, pos, num))
# ... then write lat/pos/sym to a CIF (see scripts/ in the repo, or use the
# xtalkit SHRY test helper). The key points: list only the asymmetric-unit
# atoms, and include the primitive symops as `_symmetry_equiv_pos_as_xyz`.
```

If your input is already a primitive CIF (or a P1 cell with few sites), skip
this step. For the rest of this guide we assume `LPSC_prim.cif` is such a
primitive F-43m CIF with Li1 @ 0.5.

## 2. Build `LGPS.cif`

Use the first command from `LGPS.txt` unchanged. In short, it should end with:

```bash
xtalkit build ... -o LGPS
```

Expected output:

```text
LGPS.cif
```

This is the average partially occupied structure from your refinement data.

## 3. Prepare a SHRY-ready CIF

```bash
xtalkit shry prepare LGPS.cif \
  --out LGPS_shry_ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 137 \
  --target-formula Li20Ge2P4S24 \
  --scaling-matrix 1 1 1
```

What this does:

- checks that the CIF says parent space group 137;
- fills missing vacancy occupancy with explicit `X`;
- checks every disorder group sums to occupancy 1;
- checks the 1x cell can hold integer counts;
- writes audit files.

Expected output files:

```text
LGPS_shry_ready.cif
LGPS_shry_ready.cif.manifest.json
LGPS_shry_ready.cif.orbit_grouping.json
```

If this step fails, do not continue. Fix the input CIF or the occupancy choices
first.

## 4. Count before generating

```bash
xtalkit shry count LGPS_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 \
  --angle-tolerance 5 \
  --atol 1e-5 \
  --out LGPS_shry_count.json
```

Open `LGPS_shry_count.json` and find:

```json
"count_only_result": 1234
```

The number above is an example. Copy your real number.

Decision point:

- If the count is small or reasonable, continue.
- If the count is enormous, stop and rethink the disorder definition before
  generating structures.

## 5. Generate structures

Paste the count from `LGPS_shry_count.json` into `--expect-count`.

```bash
xtalkit shry enum LGPS_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count 1234 \
  --out LGPS_SHRY \
  --remove-vacancy X \
  --target-formula Li20Ge2P4S24 \
  --write-cif \
  --write-poscar \
  --write-degeneracy
```

Replace `1234` with your real `count_only_result`.

Main output directory:

```text
LGPS_SHRY/
  input/
  raw_shry/
  clean_cif/
  poscar/
  checks/
  manifest.json
  manifest.jsonl
```

The cleaned structures are in:

```text
LGPS_SHRY/clean_cif/
```

The vacancy species `X` is removed from those cleaned CIFs.

## 6. Verify the result

```bash
xtalkit shry verify LGPS_SHRY \
  --check-count \
  --check-formula \
  --check-dedup \
  --target-formula Li20Ge2P4S24 \
  --symprec-list 1e-4 1e-3 1e-2
```

This checks:

- generated count equals count-only result;
- every cleaned structure has formula `Li20Ge2P4S24`;
- no `X` vacancy species remains;
- no obvious duplicate structure buckets exist;
- SHRY counts are stable across a few symmetry tolerances.

If this fails, inspect `LGPS_SHRY/checks/`.

## 7. Optional post-processing

Rank by shortest Li-Li distance:

```bash
xtalkit shry postprocess LGPS_SHRY --shortest-distance --pair Li Li
```

Write tblite folders and a Slurm array script:

```bash
xtalkit shry postprocess LGPS_SHRY --write-tblite --write-slurm
```

Write CP2K folders:

```bash
xtalkit shry postprocess LGPS_SHRY --write-cp2k --cp2k-template cp2k.inp
```

Ewald ranking requires explicit oxidation states:

```bash
xtalkit shry postprocess LGPS_SHRY --ewald \
  --ewald-charges Li:1 Ge:4 P:5 S:-2
```

## 8. Quick command summary for LGPS

```bash
# Keep your existing build command from LGPS.txt first.
xtalkit build ... -o LGPS

xtalkit shry prepare LGPS.cif \
  --out LGPS_shry_ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 137 \
  --target-formula Li20Ge2P4S24 \
  --scaling-matrix 1 1 1

xtalkit shry count LGPS_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 --angle-tolerance 5 --atol 1e-5 \
  --out LGPS_shry_count.json

# Replace <COUNT> with count_only_result from LGPS_shry_count.json.
xtalkit shry enum LGPS_shry_ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count <COUNT> \
  --out LGPS_SHRY \
  --remove-vacancy X \
  --target-formula Li20Ge2P4S24 \
  --write-cif --write-poscar --write-degeneracy

xtalkit shry verify LGPS_SHRY \
  --check-count --check-formula --check-dedup \
  --target-formula Li20Ge2P4S24 \
  --symprec-list 1e-4 1e-3 1e-2
```
