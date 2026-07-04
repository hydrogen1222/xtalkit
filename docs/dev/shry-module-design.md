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
