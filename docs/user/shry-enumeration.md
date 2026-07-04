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

`xtalkit shry` calls the external `shry` command. If `shry` is already on
`PATH`, nothing else is needed.

If SHRY is installed in a separate environment, point xtalkit to it:

```bash
export XTALKIT_SHRY_CMD=/path/to/shry
```

Check:

```bash
$XTALKIT_SHRY_CMD --version
```

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
