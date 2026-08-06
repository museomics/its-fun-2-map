# its-fun-2-map: Assessment, Improvements & Packaging Recommendation

## Context

The `its-fun-2-map` repository implements a bioinformatics pipeline for extracting and validating fungal ITS barcode sequences from museum specimen genome skims. It covers 12 steps: QC → Reference Retrieval → BWA Mapping → SPAdes Assembly → BLAST (rounds 1 & 2) → ITS Extraction → Summary.

---

## Critical Missing Dependencies

Two structural problems block the pipeline for any new user:

1. **`metahist_tools` package is absent from the repository.** Every script imports from it (`clean_and_tar`, `run_command`, `xlsx2csv`, `setup_logging`, `repair_reads`). The README claims it is "included in the repository" — it is not. `metahist_tools` will be released as its own PyPI package and become a declared dependency of `its-fun-2-map`.

2. **`pipeline.sh` and `its_a_summary_compiler.py` are referenced in the README but do not exist.** There is no way to run the full pipeline end-to-end.

---

## Open Bugs

### `its_primer_binding.py` (line 500)
`logger.info("\n{region_name} sequences %d:", len(samples_with_seqs))` — `{region_name}` is inside a plain string, not an f-string, so it is printed literally rather than interpolated. Fix: add the `f` prefix.

---

## Open Improvements

### `fastp_module.py`
- Q30 quality threshold is hardcoded; expose as `--quality_threshold` argument.
- The `--threads` help text should clarify that it controls the number of parallel fastp jobs, not the number of threads per job.

### `UNITEd.py`
- Samples are processed sequentially. `ThreadPoolExecutor` (already used elsewhere in the pipeline) would reduce wall-clock time for multi-sample runs.
- The UNITE database path is a required argument on every invocation. A default via an environment variable (`UNITE_DB`) would reduce friction.

### `assembly_module.py`
- When `--summary_csv` is omitted, no message is emitted; the skip is silent. Add an info log.
- The 3-stage SPAdes fallback strategy is non-obvious and should be documented in the README or in-script comments.

### `blast_round1.py` / `blast_round2.py` / parsers
- `blast_round1.py` uses `--column_name` (the tracking sheet column); the parser scripts use `--id_column` (the taxonomy CSV column). These refer to different things, but the naming is confusing across consecutive pipeline steps. Standardise the names.
- `blast_task()` in `its_fun_tools.py` accepts an `evalue` parameter but neither `blast_round1.py` nor `blast_round2.py` exposes `--evalue` as a CLI argument. The parsers default to `1e-5` while `blast_task()` defaults to `1e-10`. Wire `--evalue` through the BLAST scripts.
- Both parsers silently continue when no BLAST result file exists for a sample. Add a per-sample warning log.
- `blast_round1_parser.py` has `--evalue_cutoff` and `--allow_all`; `blast_output_parser.py` does not. Since these parse equivalent output from consecutive BLAST rounds, their interfaces should either be unified or the distinction clearly documented.

### `its_primer_binding.py`
- Processing is sequential per sample. `ThreadPoolExecutor` would reduce wall-clock time for large sample sets.

### `its_fun_tools.py`
- `log_and_print()` uses the root logger. It should accept a logger instance as a parameter to match the convention of other functions in the module.

### `its_decision_making.R`
- The file has no `if (!interactive())` block, no `commandArgs()`, and no `optparse`/`argparse` usage. It can only be sourced interactively — it cannot be called by a workflow manager as a subprocess via `Rscript`. A CLI entry point is needed.

### `tutorial.md`
- Contains only empty section headers. Needs worked examples with real arguments for each pipeline step.

---

## Cross-Cutting Issues

| Issue | Affected files |
|---|---|
| `metahist_tools` package absent from repository | All scripts — nothing runs without it |
| `pipeline.sh` and `its_a_summary_compiler.py` missing | README references both; no end-to-end entrypoint exists |
| Argument naming inconsistency (`--column_name` vs `--id_column`) | `blast_round1.py`, parsers |
| No resume/checkpoint — mid-run failure restarts from scratch | All scripts |
| `tutorial.md` empty | Docs |

---

## Snakemake vs PyPI Recommendation

### Snakemake
**Pros:**
- Native checkpointing — failed steps restart from the last completed step, not from scratch
- Explicit dependency graph makes parallelism automatic and safe
- Resource declarations (threads, memory) per rule replace the current manual thread maths
- `snakemake --use-conda` handles environment activation
- Standard in bioinformatics; collaborators and reviewers will recognise it
- Forces documentation of inputs/outputs per step

**Cons:**
- Learning curve for users unfamiliar with Snakemake
- Requires a Snakefile alongside the existing Python scripts (some duplication of CLI parsing)
- Less suitable if users want to run individual steps interactively

### PyPI Package
**Pros:**
- `pip install its-fun-2-map` is maximally accessible
- Enables `python -m its_fun_2_map.fastp_module` style invocation
- Resolves the `metahist_tools` import problem cleanly — declare it as a pip dependency
- Works well for the individual-script use case

**Cons:**
- Does not provide orchestration, dependency tracking, or checkpointing
- The user still needs to call each step manually or maintain their own shell script
- Non-Python tools (fastp, BWA, SPAdes, BLAST) are not installable via pip — conda or system packages are still required

### Recommendation: Both, in combination

The two approaches are not mutually exclusive and together solve different problems:

1. **PyPI package first**: Package the Python modules as `its-fun-2-map`, with `metahist_tools` as a declared PyPI dependency. This fixes the import problem for all scripts immediately and makes installation reproducible. Each script becomes a `console_scripts` entry point (e.g. `itsfun-qc`, `itsfun-map`). The conda YAML continues to manage non-Python tools.

2. **Snakemake wrapper on top**: A `Snakefile` wraps each step as a rule using the installed CLI entry points. Users who want push-button end-to-end execution use `snakemake`; users who want to run individual steps use the CLI entry points directly. This is the pattern used by tools like Prokka, BUSCO, and CheckM2.

This delivers:
- `pip install its-fun-2-map` (or `conda install`) for all Python logic
- `snakemake --cores N` for full pipeline runs with checkpointing
- Individual step CLIs for interactive/exploratory use
- No duplication: Snakemake rules call the same entry points as manual use

---

## Packaging Requirements

### PyPI package

1. **Package directory structure.** All `.py` scripts must move into a `its_fun_2_map/` subdirectory with an `__init__.py`. Currently all scripts are in the repo root, which is not installable as a Python package. All `from metahist_tools import ...` calls must be updated to `from its_fun_2_map.metahist_tools import ...`.

2. **`pyproject.toml`.** Needs: package metadata (name, version, authors, license); Python dependencies (`biopython`, `pandas`, `openpyxl`, `rpy2`); non-Python tools (fastp, BWA-MEM2, SPAdes, BLAST+, seqkit, BBTools) listed under `[project.optional-dependencies]` or documented as conda/system prerequisites; `console_scripts` entry points for each pipeline step (e.g. `itsfun-qc = its_fun_2_map.fastp_module:main`, `itsfun-blast = its_fun_2_map.blast_round1:main`).

3. **Package data for R scripts.** `parse_fastp_json.R` and `its_decision_making.R` must be declared as package data in `pyproject.toml` (e.g. `[tool.setuptools.package-data] its_fun_2_map = ["*.R"]`). Without this, the `Path(__file__).parent` path resolution in `fastp_module.py` will fail after `pip install` because the R files will not be copied into the installed package.

4. **`metahist_tools` as a declared dependency.** `metahist_tools` will be released as its own PyPI package. Once published, add it to `pyproject.toml` under `[project.dependencies]`. Until it is published, the five functions it provides (`xlsx2csv`, `setup_logging`, `run_command`, `clean_and_tar`, `repair_reads`) are the single biggest blocker for any user attempting to run the pipeline.

5. **Library-safe logging.** Some scripts call `setup_logging()` at module level rather than inside `main()`. When the package is imported rather than run as a CLI script, this reconfigures the calling application's root logger. Add `logging.NullHandler()` to `its_fun_2_map/__init__.py` and ensure `setup_logging()` is only called from within `main()` or `if __name__ == "__main__"` blocks.

### Snakemake / Nextflow pipeline

1. **Pipeline config YAML.** Each script currently takes its own independent set of CLI flags. A single `config.yaml` mapping all pipeline-wide parameters (sample sheet path, database paths, top-level output directory, thread counts per step, e-value thresholds) must be defined. Without it, every workflow rule would need to hard-code or duplicate arguments.

2. **Standardised output directory layout.** There is no declared convention for how step N's output directory maps to step N+1's input directory. A fixed subdirectory layout under a single top-level directory (e.g. `01_qc/`, `02_reference/`, `03_mapping/`, `04_assembly/`, `05_blast_round1/`, etc.) would allow a Snakefile to wire steps together automatically and makes output browsable.

3. **`its_decision_making.R` CLI entry point.** Snakemake and Nextflow invoke scripts as subprocesses via `Rscript`. The current file has no argument parsing — it cannot be called from a workflow rule. Add a CLI entry point accepting at minimum `--input_csv` and `--output_csv`.

4. **`its_a_summary_compiler.py` must exist.** The final workflow rule cannot be written until this script has a defined interface. At minimum it needs `--input_dir`, `--output_file`, and `--log_file` arguments. It calls `its_outcome()` from `its_decision_making.R` via rpy2 and writes the consolidated report.

5. **BLAST parser interface parity.** `blast_round1_parser.py` has `--evalue_cutoff` and `--allow_all`; `blast_output_parser.py` does not. A workflow config cannot share a single parameterisation for both parsing steps with diverged interfaces. Unify them or document the distinction explicitly.

---

## Implementation Plan

### Phase 1 — Complete
All per-script bugs listed in the original assessment have been resolved:
- Removed NHM-specific conda `prefix:` from `its-fun-2-map.yaml`
- Fixed missing imports in `its_fun_tools.py` (`time`, `random`, `urllib.error`)
- Fixed `load_name_ids` return-type inconsistency in `its_fun_tools.py`
- Fixed args scoping and undefined variable in `pull_ncbi_lineage.py`
- Fixed `qseqid` header inclusion in `blast_round2.py`
- Fixed scenario-6 copy-paste and `df$df$` double-prefix in `its_decision_making.R`
- Fixed R script relative path in `fastp_module.py`
- Fixed hardcoded CSV fieldnames in `its_primer_binding.py`
- Removed duplicate `get_ncbi_lineage` from `UNITEd.py`
- Exposed BWA thread count as `--bwa_threads` in `mapping_module.py`

Remaining Phase 1 items:
- `metahist_tools` package: still absent from repository (biggest blocker)
- `its_primer_binding.py` line 500: missing f-string prefix (newly identified)

### Phase 2 — Open
- Fix `its_primer_binding.py` line 500 f-string bug
- Standardise argument naming across BLAST scripts and parsers
- Expose `--evalue` through `blast_round1.py` and `blast_round2.py`
- Align `blast_round1_parser.py` and `blast_output_parser.py` argument sets
- Add `ThreadPoolExecutor` to `its_primer_binding.py` and `UNITEd.py`
- Add info log when `--summary_csv` is omitted in `assembly_module.py`
- Add per-sample warning in BLAST parsers when no result file is found
- Standardise `log_and_print()` to accept a logger instance
- Write `tutorial.md` content with worked examples

### Phase 3 — Open
- Publish `metahist_tools` to PyPI; add it as a dependency in `pyproject.toml`
- Create `its_fun_2_map/` package directory with `__init__.py`; move scripts in; update imports
- Write `pyproject.toml` with metadata, dependencies, and `console_scripts` entry points
- Declare `parse_fastp_json.R` and `its_decision_making.R` as package data
- Ensure `setup_logging()` is only called inside `main()`
- Add CLI entry point to `its_decision_making.R`
- Write `its_a_summary_compiler.py` with a defined CLI interface
- Define `config.yaml` schema and output directory layout convention
- Write `Snakefile` with one rule per pipeline step
- Recreate `pipeline.sh` as a thin wrapper over `snakemake`

---

## Verification

- `conda env create -f its-fun-2-map.yaml` on a non-NHM machine succeeds
- `pip install -e .` from repo root succeeds; `itsfun-qc --help` prints usage
- `python -c "import its_fun_2_map"` succeeds with no logging side-effects
- `python blast_round2.py ...` does not pass `qseqid` to seqkit
- `Rscript its_decision_making.R --input_csv test.csv --output_csv out.csv` exits 0
- `snakemake -n --configfile config.yaml` produces the expected DAG without errors
