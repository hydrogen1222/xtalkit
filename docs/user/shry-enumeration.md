# Running a SHRY enumeration safely

This guide is for a real calculation, not just a CLI demonstration. Read it
once before submitting a large occupational-disorder problem.

SHRY and enumlib solve the same broad problem but have different working
styles. `xtalkit enumerate` is convenient for small and medium primitive
structures. The SHRY path is better when you need an exact count before
generation and a recorded check at every stage.

```text
input CIF
   │
   ├─ prepare  → corrected/checked partial-occupancy CIF
   ├─ count    → number of symmetry-inequivalent structures
   ├─ enum     → raw and cleaned ordered structures
   ├─ verify   → count, formula, vacancy, and duplicate report
   └─ postprocess (optional) → rankings or tblite job folders
```

## 1. Install the two environments

SHRY 1.1.8 requires an older pymatgen than xtalkit. Keep the executable in an
isolated Python 3.11 uv tool environment:

```bash
uv python install 3.11
uv tool install --python 3.11 'shry==1.1.8'
shry --help
uv tool list
```

Do not install SHRY with `uv pip install` inside this repository. If the
command is not on `PATH`, point xtalkit to it explicitly:

```bash
export XTALKIT_SHRY_CMD="$(command -v shry)"
```

The xtalkit-side audits use current pymatgen and spglib:

```bash
uv sync --extra enumerate
```

Compiled enumlib is not required unless you also request the optional
cross-backend check.

## 2. Check the input before running anything

You should know all of the following:

1. Which cell is in the CIF: primitive or conventional.
2. The intended parent space group and setting.
3. Which rows describe the same disordered crystallographic site.
4. The exact ordered-cell formula after applying the scaling matrix.
5. Whether every partial occupancy can become an integer atom count.

Do not infer these from the filename. Inspect the CIF tags, refinement table,
and chemical formula.

### F- and I-centred conventional cells

SHRY enumerates the cell it receives. It does not automatically convert an
F- or I-centred conventional cell to a primitive one. That distinction can be
enormous: choosing 24 atoms among 48 sites has roughly (10^{13}) raw
combinations, while choosing 6 among 12 has only 924.

Whenever possible, export a primitive cell from the refinement or structure
program that produced the data. If you must convert with pymatgen, the
following produces a P1 CIF containing every site of the primitive cell:

```python
import gemmi
from pymatgen.core import Structure
from pymatgen.io.cif import CifWriter
from pymatgen.symmetry.analyzer import SpacegroupAnalyzer

source = Structure.from_file("conventional.cif")
analyzer = SpacegroupAnalyzer(source, symprec=0.01, angle_tolerance=5)
primitive = analyzer.get_primitive_standard_structure()

# P1 output lists the complete primitive cell. The prepare stage below will
# independently detect its actual parent symmetry. Pymatgen gives the species
# rows of one disordered site the same label, so make every CIF row unique.
CifWriter(primitive, symprec=None).write_file("primitive_p1_raw.cif")
document = gemmi.cif.read_file("primitive_p1_raw.cif")
labels = document.sole_block().find_loop("_atom_site_label")
for index in range(len(labels)):
    labels[index] = f"site{index + 1}"
document.write_file("primitive_p1.cif")

print("conventional sites:", len(source))
print("primitive sites:", len(primitive))
print("primitive formula:", primitive.composition)
```

This conversion is not self-validating. Compare the site count, composition,
lattice, and partial occupancies with the source, then remove or archive the
intermediate `primitive_p1_raw.cif`. In the next step, always pass
`--parent-spacegroup`; xtalkit will reject the final file if independent
detection does not recover the intended group.

Do not attach conventional F/I symmetry operations to arbitrary primitive
vectors. A P1 file listing the complete primitive cell is less elegant but
unambiguous.

## 3. Prepare the SHRY input

The target formula is the **exact formula of the cell that SHRY will
enumerate**, after multiplying by the scaling-matrix determinant and removing
vacancies. It is not automatically reduced.

```bash
uv run xtalkit shry prepare primitive_p1.cif \
  --out ready.cif \
  --vacancy-symbol X \
  --parent-spacegroup 216 \
  --target-formula Li6PS5Cl \
  --scaling-matrix 1 1 1
```

Preparation does four important things:

- completes a partially occupied site with explicit `X` vacancy rows;
- rejects site occupancies below zero, above one, or summing above one;
- checks integer atom counts and the exact target formula;
- detects the parent symmetry independently and records site orbits.

It writes:

```text
ready.cif
ready.cif.manifest.json
ready.cif.orbit_grouping.json
```

Read the manifest before continuing. Check at least
`parent_spacegroup_detected`, `supercell_matrix`, `site_orbits`, and
`target_formula_after_removing_vacancy`.

### Deliberate occupancy idealization

Only idealize a refined decimal when the crystallographic model justifies it.
For example:

```bash
uv run xtalkit shry prepare input.cif \
  --out ready.cif \
  --set-occupancy "Li1:11/16 M1:Ge1/2P1/2" \
  --parent-spacegroup 137 \
  --target-formula Li20Ge2P4S24
```

The label before each colon must exist in the input CIF. The exact fractions
and the resulting formula are stored in the manifest; no silent rounding is
performed.

## 4. Count before generating

```bash
uv run xtalkit shry count ready.cif \
  --scaling-matrix 1 1 1 \
  --symprec 0.01 \
  --angle-tolerance 5 \
  --atol 1e-5 \
  --out count.json
```

Inspect the result with any JSON viewer, for example:

```bash
python -m json.tool count.json
```

`count_only_result` is the number of inequivalent output structures.
`total_combinations`, when reported by SHRY, is the unsymmetrized combination
count. Also check the input hash, matrix, tolerances, and backend version.

Do not proceed merely because the command finished. Estimate the storage and
post-processing cost of `count_only_result`. If the result is unexpectedly
large, re-check the cell choice and disorder model rather than starting
generation.

## 5. Generate into a new directory

Copy the exact count into `--expect-count`:

```bash
uv run xtalkit shry enum ready.cif \
  --scaling-matrix 1 1 1 \
  --expect-count <COUNT_FROM_JSON> \
  --remove-vacancy X \
  --target-formula Li6PS5Cl \
  --write-cif --write-poscar \
  --out SHRY_run
```

`SHRY_run` must be new or empty. A successful run contains:

```text
SHRY_run/
├── input/       copied input, command, SHRY environment, mod-only result
├── raw_shry/    original SHRY CIFs, including vacancy pseudo-species
├── clean_cif/   ordered CIFs after vacancy removal
├── poscar/      optional POSCAR files
├── checks/      verification and ranking reports
└── manifest.json
```

Before full generation, xtalkit runs SHRY's `--mod-only` path and compares the
expanded composition with the prepared input. A chemistry change stops strict
mode.

If generation fails, the directory is intentionally retained. Read
`manifest.json` for `failure_stage` and inspect the copied command and raw
output. Do not mix a retry into the same directory.

## 6. Verify the result

```bash
uv run xtalkit shry verify SHRY_run \
  --target-formula Li6PS5Cl
```

The default verification checks:

- generated count against the expected count;
- exact ordered formula and occupancy one;
- removal of the vacancy pseudo-species;
- duplicates confirmed by pymatgen `StructureMatcher`.

For a tolerance-sensitivity check:

```bash
uv run xtalkit shry verify SHRY_run \
  --target-formula Li6PS5Cl \
  --symprec-list 1e-4 1e-3 1e-2
```

All requested tolerances must give a stable count. A change is a reason to
inspect nearly symmetric coordinates, not a value to average away.

### Degeneracy

`--write-degeneracy` and `--check-degeneracy` are valid only when the scaling
matrix determinant is one and the CIF-declared group equals the independently
detected enumeration group. A genuine supercell needs translational operations
that are not represented by the parent-cell action, so xtalkit refuses to
claim a rigorous degeneracy there.

A P1 primitive file whose higher symmetry is detected independently can still
be enumerated, but it cannot satisfy that declared-group requirement. Leave
degeneracy disabled unless you have a correctly declared representation.

## 7. Optional post-processing

Rank by shortest distance for a chosen pair:

```bash
uv run xtalkit shry postprocess SHRY_run \
  --shortest-distance --pair Li Li
```

Rank by charge-neutral ordered Ewald energy:

```bash
uv run xtalkit shry postprocess SHRY_run \
  --ewald --ewald-charges Li:1 P:5 S:-2 Cl:-1
```

Generate periodic tblite POSCAR jobs and a Slurm array:

```bash
uv run xtalkit shry postprocess SHRY_run \
  --write-tblite --charge 0 --uhf 0
uv run xtalkit shry postprocess SHRY_run --write-slurm
```

Install tblite as described in the main README. Review charge and number of
unpaired electrons in every generated `run.sh`; the defaults describe a
neutral singlet and may be wrong for your material.

## 8. Reproducibility checklist

Before reporting a production enumeration, retain and review:

- the original and prepared CIFs;
- exact target formula and scaling matrix;
- parent space group, setting, and detected Hall symbol;
- SHRY version and isolated package snapshot;
- count JSON and successful generation manifest;
- verification report and any tolerance scan;
- the code revision used for the run.

The count and generated count must agree. Every cleaned structure must have the
target formula, no vacancy pseudo-species, and no confirmed duplicate.

## 9. Common failures

**Detected group differs from `--parent-spacegroup`**
: Check the cell, setting, coordinates, and symmetry tolerance. Do not accept
  the detected group until it agrees with the intended parent model.

**Occupancy cannot be integerized**
: The orbit multiplicity and supercell determinant cannot produce integer atom
  counts. Correct the disorder model or choose a justified scaling matrix.

**`--mod-only` changed the chemistry**
: Compare `input/shry_ready.cif` and `input/shry_mod_only.cif`. Do not continue
  strict generation until the change is understood.

**Generated count differs from count-only result**
: Keep the failed directory and compare command lines, tolerances, input hash,
  and SHRY version. Never merge partial output with another run.

**Output directory is not empty**
: Choose a fresh path. This rule prevents stale structures from being counted
  as results of the current run.
