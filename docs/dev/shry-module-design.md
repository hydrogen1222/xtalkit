# SHRY integration: maintainer notes

This document describes the invariants and external interfaces that must remain
true when the SHRY workflow is changed. User instructions live in
`docs/user/shry-enumeration.md` and its Chinese translation.

## Module boundaries

- `shry_backend.py` is the only direct wrapper around the external `shry`
  executable. `XTALKIT_SHRY_CMD` may select an isolated installation.
- `shry_prepare.py` owns CIF normalization, exact occupancy checks, vacancy
  completion, target-formula validation, and independent symmetry auditing.
- `shry_count.py` runs count-only and records the input hash, matrix,
  tolerances, raw combination count, inequivalent count, and backend version.
- `shry_enum.py` owns output-directory isolation, mod-only chemistry checking,
  raw-output discovery, vacancy removal, formula validation, and manifests.
- `shry_verify.py` owns count, formula, vacancy, duplicate, tolerance-scan,
  optional cross-backend, and optional degeneracy checks.
- `postprocess.py` ranks cleaned structures and writes periodic tblite jobs.

SHRY is never imported as a Python library. This prevents its pinned pymatgen
from contaminating xtalkit's current pymatgen environment.

## Scientific invariants

Preparation must:

1. parse occupancies exactly with `Fraction`;
2. normalize real element symbols while keeping an explicit vacancy species;
3. require each disorder group to sum to one after vacancy completion;
4. require integer atom counts for the requested scaling determinant;
5. compare the multiplicity-weighted ordered-cell formula with the exact
   user-supplied target formula;
6. detect the parent group independently with spglib;
7. reject a detected group that differs from `--parent-spacegroup`;
8. reject a spglib orbit split across incompatible CIF coordinate groups.

Generation must use a new or empty directory. A failure after layout creation
must retain diagnostics and write a failed manifest; it must never leave a
success manifest. Strict mode requires an expected count and rejects a
mod-only chemistry change.

Duplicate detection uses exact composition and site count only to choose
candidate pairs. Every candidate pair is decided by pymatgen
`StructureMatcher`; coordinate fingerprints are diagnostic and may not exclude
a pair.

Degeneracy is reported only for determinant-one matrices when the CIF-declared
group equals the independently detected enumeration group. Parent-cell
fractional operations are insufficient for a genuine supercell.

Ewald ranking requires ordered, charge-neutral structures. tblite jobs use
periodic POSCAR input and explicit charge/unpaired-electron settings.

## External SHRY 1.1.x contract

The wrapper currently relies on this observed CLI behaviour:

- there are no `--output-dir`, `--write-cif`, or `--write-poscar` options;
- output is written below `shry-<basename>-<timestamp>/` in the process working
  directory;
- ordered structures are below `slice*/`;
- the modified disordered parent is written outside `slice*/`;
- `--count-only` prints `Expected unique patterns is N` and may also print
  `Total number of combinations is M`;
- `--count-only` writes no structure files;
- `--mod-only` writes only the modified parent;
- a scaling matrix is accepted as one, three, or nine integers;
- SHRY does not write POSCAR files.

If that contract changes, inspect `ShryBackend`, `parse_count_output`,
`iter_shry_outputs`, and the mod-only output discovery together. Do not patch a
single parser without running a real-backend integration test.

## Manifests and evidence

The SHRY interpreter is recovered from the executable shebang when possible.
Backend version and isolated-environment package metadata must describe that
interpreter, not the xtalkit interpreter.

Individual reports belong under `checks/`; `verify.json` is the aggregate
result. The input hash, actual command line, tolerances, matrix, count, formula,
and backend version are required reproducibility fields.

## Testing rules

- Unit tests use a deterministic fake backend and must cover failure manifests,
  count parsing, chemistry changes, vacancy removal, and duplicate logic.
- A real SHRY executable enables a two-site count-and-enumerate smoke test.
- Do not put a production LGPS or argyrodite enumeration in the test suite.
- Keep every real smoke structure small and constrain numerical-library thread
  counts in CI or local validation.
- Changes to CIF parsing require round-trip tests through gemmi and pymatgen.

External SHRY, supercell, enumlib, and tblite executables remain optional and
are never bundled with the Python package.
