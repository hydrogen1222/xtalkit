# SHRY Module Design

This module implements Phase 1 of `xtalkit_SHRY模块开发方案_修订版.md`.

Implemented scope:

- `xtalkit.enumeration` package with SHRY prepare/count/enum/verify workflows.
- SHRY is invoked only through CLI subprocesses via `ShryBackend`.
- `XTALKIT_SHRY_CMD` can point to an isolated SHRY environment.
- `prepare` uses exact `Fraction` occupancy checks, fills explicit vacancy
  rows, preserves CIF symmetry operations, records symmetry orbit data, and
  writes manifest sidecars.
- `count` records `symprec`, `angle_tolerance`, and `atol`.
- `enum` requires `expect_count` in strict mode, runs `--mod-only`, streams raw
  CIF outputs, removes vacancy species, checks formula, and writes
  `manifest.jsonl`.
- `verify` uses hash buckets for duplicate detection instead of all-pairs
  `StructureMatcher`, supports SHRY symprec count scans, optional supercell
  cross-backend count checks, and optional degeneracy-sum checks.
- `supercell_backend` wraps the external supercell CLI through
  `XTALKIT_SUPERCELL_CMD`.
- `degeneracy` computes `|G| / |stabilizer|` from parent CIF symmetry
  operations and can write per-configuration values during `shry enum`.
- `postprocess` supports shortest-distance ranking, optional Ewald ranking,
  tblite/CP2K job folders, and Slurm array script generation.

Non-goals:

- No SHRY Python API calls.
- No random sampling or energy filtering.
- No automatic large LGPS enumeration in tests.
- No bundled SHRY or supercell executables; both are configured externally.
- No automatic production LGPS enumeration in tests.

The legacy enumlib module remains in `xtalkit.enumerator` and is not removed.

## Bug history — 2026-07-04 (SHRY backend made runnable)

The Phase 1 code assumed a SHRY CLI that did not match the real SHRY 1.1.x.
The module compiled but could not run end-to-end. Fixed and verified on LPSC
(Li-vacancy, primitive F-43m): prepare → count(48) → enum(48) → verify all
pass, degeneracy_sum = 924 = C(12,6), and a brute-force cross-check confirmed
the 48 generated configs cover all 48 ground-truth orbits (no miss, no dup).

Root causes and fixes:

- **CIF esds** — `cif_io.read_cif` now strips `0.3148(19)`-style standard
  uncertainties from cell and atom-site numeric fields.
- **Occupancy / species notation** — `write_cif` emits decimal occupancies
  (`1/2` → `0.5`) and element-normalized type_symbols (`Li1+` → `Li`); both
  are required for pymatgen (which SHRY uses) to parse the CIF.
- **Non-existent SHRY flags** — dropped `--output-dir`, `--write-cif`,
  `--write-poscar` (none exist in SHRY). SHRY is run with `cwd=raw_dir`; it
  writes `shry-<base>-<ts>/slice*/` there. Clean CIFs and POSCARs are built
  by xtalkit in post-processing (SHRY has no POSCAR writer).
- **Count parsing** — `parse_count_output` now matches `Expected unique
  patterns is N` (the old `inequivalent structures: N` regex also mis-matched
  the raw `Total number of combinations is N` line).
- **Output discovery** — `iter_shry_outputs` only yields `slice*/` CIFs,
  skipping SHRY's disordered modified-parent CIF that crashed
  `counts_from_atom_rows`.
- **mod-only audit** — finds SHRY's modified CIF on disk (was checking stdout
  for `"data_"`, which never appears) and compares reduced, element-normalized
  composition; `Li1+` vs `Li+` no longer false-positives.
- **`--strict` flag** — was `store_true` + `default=True` (a no-op); switched
  to `BooleanOptionalAction` so `--no-strict` works.

### SHRY CLI contract (reference for future maintenance)

SHRY 1.1.x behaviour the integration depends on:

- No `--output-dir` / `--write-cif` / `--write-poscar` flags.
- Writes `shry-<basename>-<unixtime>/sliceN/<formula>_<i>_<weight>.cif`
  (ordered configs) + `<basename>-<scaling>.cif` (disordered modified parent)
  + `sub.log` into its CWD.
- `--count-only` prints `Expected unique patterns is N` (and
  `Total number of combinations is M`); writes nothing.
- `--mod-only` writes only the modified-parent CIF.
- `--scaling-matrix` accepts 1, 3, or 9 integers.
- POSCAR output is not supported — xtalkit generates POSCARs itself.

If SHRY changes any of these, `shry_backend.py`, `iter_shry_outputs`, and
`parse_count_output` are the coupling points.
