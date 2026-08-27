# its-fun-2-map: Assessment, Improvements & Packaging Recommendation

## Context

The `its-fun-2-map` repository implements a bioinformatics pipeline for extracting and validating fungal ITS barcode sequences from museum specimen genome skims. It covers 12 steps: QC → Reference Retrieval → BWA Mapping → SPAdes Assembly → BLAST (rounds 1 & 2) → ITS Extraction → Summary.

---

## Critical Missing Dependencies

One structural problem blocks the pipeline for any new user:

1. **`pipeline.sh` is referenced in the README but does not exist.** The README states "All pipeline steps can be run sequentially using the `pipeline.sh` script" (line 67) but the file is absent. There is no way to run the full pipeline end-to-end without manually chaining each step.

`its_a_summary_compiler.py` was previously missing but has now been added to the repository (Step 12).

`seqpy-tools` (previously `metahist_tools`) is now available on PyPI and declared in the conda environment YAML — all import lines across the repository have been updated to `from seqpy_tools import ...`.

---

## Open Improvements

### `UNITEd.py`
- Samples are processed sequentially. `ThreadPoolExecutor` (already used in `blast_round1.py` and `blast_round2.py`) would reduce wall-clock time for multi-sample runs.
- The UNITE database path is a required argument on every invocation. A default via an environment variable (`UNITE_DB`) would reduce friction.

### `assembly_module.py`
- When `--summary_csv` is omitted, no message is emitted; the skip is silent. Add an `else: logger.info("--summary_csv not specified, skipping CSV summary")` branch.
- The 3-stage SPAdes fallback strategy is non-obvious and should be documented in the README or in-script comments.

### `blast_round1.py` / `blast_round2.py` / parsers
- `blast_task()` in `its_fun_tools.py` accepts an `evalue` parameter but neither `blast_round1.py` nor `blast_round2.py` exposes `--evalue` as a CLI argument. The parsers default to `1e-5` while `blast_task()` defaults to `1e-10`. Wire `--evalue` through the BLAST scripts.
- Both parsers silently continue when no BLAST result file exists for a sample. Add a per-sample `logger.warning` (already done in `blast_round1_parser.py` for the assembly file; do the same for the BLAST TSV).
- `blast_round1_parser.py` has `--evalue_cutoff` and `--allow_all`; `blast_output_parser.py` does not. Since these parse equivalent output from consecutive BLAST rounds, their interfaces should either be unified or the distinction clearly documented.

### `its_primer_binding.py`
- Processing is sequential per sample. `ThreadPoolExecutor` would reduce wall-clock time for large sample sets.
- **The sample categorisation block is not region-agnostic.** `failed_samples`, `complete_samples`, `its1_only_samples` and `its2_only_samples` are derived from hardcoded `results[sample]['ITS_complete']`, `['ITS1']` and `['ITS2']` lookups. Because `results` is a nested `defaultdict(int)`, a custom `--regions_tsv` run does not crash — it silently reports every sample as failed. The three closing `logger.info` lines listing output directories are hardcoded to the same three regions. As of v3.0.1 `extraction_summary.csv` is correct for custom regions but `summary_report.txt` is not.
- `complete_samples`, `its1_only_samples` and `its2_only_samples` are populated but never read. Either report them or remove them.
- The per-sample processing loop uses bare `print()` rather than the logger, so progress output does not reach the log file.

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
| `pipeline.sh` missing | README references it; no end-to-end entrypoint exists |
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
- Resolves the `seqpy-tools` import problem cleanly — declare it as a pip dependency
- Works well for the individual-script use case

**Cons:**
- Does not provide orchestration, dependency tracking, or checkpointing
- The user still needs to call each step manually or maintain their own shell script
- Non-Python tools (fastp, BWA, SPAdes, BLAST+) are not installable via pip — conda or system packages are still required

### Recommendation: Both, in combination

The two approaches are not mutually exclusive and together solve different problems:

1. **PyPI package first**: Package the Python modules as `its-fun-2-map`, with `seqpy-tools` as a declared PyPI dependency. Each script becomes a `console_scripts` entry point (e.g. `itsfun-qc`, `itsfun-map`). The conda YAML continues to manage non-Python tools.

2. **Snakemake wrapper on top**: A `Snakefile` wraps each step as a rule using the installed CLI entry points. Users who want push-button end-to-end execution use `snakemake`; users who want to run individual steps use the CLI entry points directly. This is the pattern used by tools like Prokka, BUSCO, and CheckM2.

This delivers:
- `pip install its-fun-2-map` (or `conda install`) for all Python logic
- `snakemake --cores N` for full pipeline runs with checkpointing
- Individual step CLIs for interactive/exploratory use
- No duplication: Snakemake rules call the same entry points as manual use

---

## Packaging Requirements

### PyPI package

1. **Package directory structure.** All `.py` scripts must move into a `its_fun_2_map/` subdirectory with an `__init__.py`. Currently all scripts are in the repo root, which is not installable as a Python package.

2. **`pyproject.toml`.** Needs: package metadata (name, version, authors, license); Python dependencies (`biopython`, `pandas`, `openpyxl`, `rpy2`, `seqpy-tools`); non-Python tools (fastp, BWA-MEM2, SPAdes, BLAST+, seqkit, BBTools) listed under `[project.optional-dependencies]` or documented as conda/system prerequisites; `console_scripts` entry points for each pipeline step (e.g. `itsfun-qc = its_fun_2_map.fastp_module:main`, `itsfun-blast = its_fun_2_map.blast_round1:main`).

3. **Package data for R scripts.** `parse_fastp_json.R` and `its_decision_making.R` must be declared as package data in `pyproject.toml` (e.g. `[tool.setuptools.package-data] its_fun_2_map = ["*.R"]`). Without this, the `Path(__file__).parent` path resolution in `fastp_module.py` will fail after `pip install` because the R files will not be copied into the installed package.

4. **Library-safe logging.** Some scripts call `setup_logging()` at module level rather than inside `main()`. When the package is imported rather than run as a CLI script, this reconfigures the calling application's root logger. Add `logging.NullHandler()` to `its_fun_2_map/__init__.py` and ensure `setup_logging()` is only called from within `main()` or `if __name__ == "__main__"` blocks.

### Snakemake / Nextflow pipeline

1. **Pipeline config YAML.** Each script currently takes its own independent set of CLI flags. A single `config.yaml` mapping all pipeline-wide parameters (sample sheet path, database paths, top-level output directory, thread counts per step, e-value thresholds) must be defined.

2. **Standardised output directory layout.** There is no declared convention for how step N's output directory maps to step N+1's input directory. A fixed subdirectory layout under a single top-level directory (e.g. `01_qc/`, `02_reference/`, `03_mapping/`, `04_assembly/`, `05_blast_round1/`, etc.) would allow a Snakefile to wire steps together automatically.

3. **`its_decision_making.R` CLI entry point.** Snakemake and Nextflow invoke scripts as subprocesses via `Rscript`. The current file has no argument parsing — it cannot be called from a workflow rule. Add a CLI entry point accepting at minimum `--input_csv` and `--output_csv`.

4. **BLAST parser interface parity.** `blast_round1_parser.py` has `--evalue_cutoff` and `--allow_all`; `blast_output_parser.py` does not. A workflow config cannot share a single parameterisation for both parsing steps with diverged interfaces. Unify them or document the distinction explicitly.

---

## Implementation Plan

### Phase 1 — Complete
All per-script bugs have been resolved:
- Removed NHM-specific conda `prefix:` from `its-fun-2-map.yaml`
- Fixed missing imports in `its_fun_tools.py` (`time`, `random`, `urllib.error`)
- Fixed `load_name_ids` return-type inconsistency in `its_fun_tools.py`
- Fixed args scoping and undefined variable in `pull_ncbi_lineage.py`
- Fixed `qseqid` header inclusion in `blast_round2.py`
- Fixed scenario-6 copy-paste and `df$df$` double-prefix in `its_decision_making.R`
- Fixed R script relative path in `fastp_module.py`
- Removed duplicate `get_ncbi_lineage` from `UNITEd.py`
- Exposed BWA thread count as `--bwa_threads` in `mapping_module.py`
- Replaced all `from metahist_tools import ...` with `from seqpy_tools import ...` across all 11 scripts; added `seqpy-tools` to the pip section of `its-fun-2-map.yaml`
- Exposed `--quality_threshold` (default 30) and `--fastp_threads` (default 3) in `fastp_module.py`; clarified `--threads` help text to describe parallel jobs not threads per job
- Added `its_a_summary_compiler.py` (Step 12 — previously missing)

`its_primer_binding.py` v3.0.0:
- Fixed `read_tracking_sheet()` logger default of `None` causing `AttributeError` on all three error paths (missing column, `FileNotFoundError`, generic read failure); `logger` is now a required positional and is passed from `main()`

`its_primer_binding.py` v3.0.1:
- Fixed missing f-string prefix on the per-region `logger.info` call in the console summary block, which printed `{region_name}` literally
- Derived `extraction_summary.csv` `fieldnames` from `REGIONS.keys()` instead of hardcoding `ITS1`/`ITS2`/`ITS_complete`. Custom `--regions_tsv` runs previously raised `ValueError: dict contains fields not in fieldnames` on the first `writer.writerow()`
- Fixed `run_seqkit_amplicon()` logger default of `None`; a `CalledProcessError` from seqkit raised `AttributeError` on `logger.error()` instead of logging and returning `False`. `logger` is now a required positional and is passed from `main()`
- Added per-sample warning in `blast_output_parser.py` when no BLAST result file found for a sample in the taxonomy mapping
- Added info log in `assembly_module.py` when `--summary_csv` is omitted
- Added per-sample warning in `blast_round1_parser.py` when no BLAST result file found for a sample in the taxonomy mapping
- Exposed `--evalue` (default `1e-5`) through `blast_round1.py` and `blast_round2.py`; wired through to `blast_task()` in `its_fun_tools.py`
- Reported `complete_samples`, `its1_only_samples`, `its2_only_samples` in `its_primer_binding.py` summary report and console log (previously populated but silently discarded)
- Standardised `log_and_print()` in `its_fun_tools.py` to accept an optional `logger` instance (falls back to root logger when omitted — fully backward-compatible)
- Confirmed no bare `print()` calls exist in `its_primer_binding.py` sample loop; all output goes through logger
- Made `its_primer_binding.py` categorisation block and closing directory log lines region-agnostic: failed detection now checks all `REGIONS` keys; ITS-specific sub-categories (`complete_samples`, `its1_only_samples`, `its2_only_samples`) are computed and reported only when running with the default ITS regions (`ITS_complete`, `ITS1`, `ITS2` all present); closing output-directory log lines now iterate over `REGIONS` instead of hardcoding ITS paths

### Phase 2 — Open
- Write `tutorial.md` content with worked examples

### Phase 3 — Open
- Create `its_fun_2_map/` package directory with `__init__.py`; move scripts in; update imports
- Write `pyproject.toml` with metadata, dependencies, and `console_scripts` entry points
- Declare `parse_fastp_json.R` and `its_decision_making.R` as package data
- Ensure `setup_logging()` is only called inside `main()`
- Add CLI entry point to `its_decision_making.R`
- Define `config.yaml` schema and output directory layout convention
- Write `Snakefile` with one rule per pipeline step
- Write `pipeline.sh` as a thin wrapper over `snakemake`

### Stretch Goals / Future Development
- Add `ThreadPoolExecutor` to `its_primer_binding.py` — extract a `process_sample()` function covering file lookup + seqkit subprocess calls per region; collect results with `as_completed`; `args.threads` controls worker count
- Add `ThreadPoolExecutor` to `UNITEd.py` — collect all sample rows first, then fan out NCBI Entrez calls; limit `max_workers` to 2 without an API key, ~4 with one to respect rate limits; `unite_index` is read-only so safe to share across threads
- Add `--evalue_cutoff` and `--allow_all` to `blast_output_parser.py` to match `blast_round1_parser.py` interface (unifies parameterisation across both BLAST parsing steps; required for a shared workflow config)
- Add `UNITE_DB` environment variable default to `UNITEd.py` (removes need to pass `--unite_db` on every invocation when the path is stable)

---

## Verification

- `conda env create -f its-fun-2-map.yaml` on a non-NHM machine succeeds
- `pip install -e .` from repo root succeeds; `itsfun-qc --help` prints usage
- `python -c "import its_fun_2_map"` succeeds with no logging side-effects
- `Rscript its_decision_making.R --input_csv test.csv --output_csv out.csv` exits 0
- `snakemake -n --configfile config.yaml` produces the expected DAG without errors
