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

### `assembly_module.py`
- The 3-stage SPAdes fallback strategy (full assembly → contigs only → k55 forced) is non-obvious to anyone reading the script cold. It should be documented with an inline comment explaining what failure condition triggers each fallback and why k55 is the last resort.

### `its_decision_making.R`
- The file defines a single function (`its_outcome()`) that is sourced at runtime by `its_a_summary_compiler.py` via rpy2 — no standalone CLI is needed for the current pipeline. A CLI entry point would only be required if `its_decision_making.R` were invoked as a direct `Rscript` subprocess (e.g. in an alternative Nextflow pipeline), which is a stretch goal.

### `tutorial.md`
- Contains only empty section headers. Needs worked examples with real command-line invocations for each pipeline step (tracked in Phase 2).

All other per-script improvements identified in the original assessment have been resolved or moved to Stretch Goals — see Phase 1 complete list and Stretch Goals section below.

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

4. **Library-safe logging.** All scripts currently call `setup_logging()` inside their entry-point functions rather than at module level, so there is no immediate logging side-effect on import. When packaged, add `logging.NullHandler()` to `its_fun_2_map/__init__.py` as a belt-and-braces measure, and verify that no future script adds a top-level `setup_logging()` call.

### Snakemake / Nextflow pipeline

1. **Pipeline config YAML.** Each script currently takes its own independent set of CLI flags. A single `config.yaml` mapping all pipeline-wide parameters (sample sheet path, database paths, top-level output directory, thread counts per step, e-value thresholds) must be defined.

2. **Standardised output directory layout.** There is no declared convention for how step N's output directory maps to step N+1's input directory. A fixed subdirectory layout under a single top-level directory (e.g. `01_qc/`, `02_reference/`, `03_mapping/`, `04_assembly/`, `05_blast_round1/`, etc.) would allow a Snakefile to wire steps together automatically.

3. **`its_decision_making.R` path resolution.** `its_a_summary_compiler.py` currently sources the R file with a bare filename (`r['source']('its_decision_making.R')`), which relies on the R working directory matching the repo root. Fixed in Phase 1 to use `Path(__file__).parent / 'its_decision_making.R'`, matching the pattern already applied to `fastp_module.py`.

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
- Fixed `its_decision_making.R` path resolution in `its_a_summary_compiler.py`: bare `r['source']('its_decision_making.R')` replaced with `Path(__file__).parent / 'its_decision_making.R'`, matching the pattern already applied to `fastp_module.py`

`its_fun_tools.py`, `pull_ncbi_lineage.py`, `blast_round1.py`, `blast_round2.py`, `assembly_module.py`, `fastp_module.py`, `blast_round1_parser.py`, `its_a_summary_compiler.py` v1.x:
- Fixed `load_name_ids` XLSX branch in `its_fun_tools.py`: `xlsx2csv()` returns a path string — previously assigned directly to `df`, causing `AttributeError` on column access; fixed with `pd.read_csv(xlsx2csv(...))`
- Same XLSX fix applied to `add_ncbi_lineages_to_csv` in `pull_ncbi_lineage.py`
- Fixed `process_renaming_df` in `its_a_summary_compiler.py`: `logger` was referenced as a free variable but only existed in `main()`; `logger` is now a required positional argument, passed from the call site
- Fixed task exceptions silently swallowed in `blast_round1.py`: added `for future in as_completed(tasks): future.result()` inside the `with ThreadPoolExecutor` block so BLAST errors are logged rather than lost
- Fixed `--sheet` not forwarded to `load_name_ids` in `blast_round2.py`: added `sheet=None` parameter to `run_seqkit_and_blast`, wired through to `load_name_ids`, and passed `sheet=args.sheet` at the call site
- Removed `logger` from `assembly_module.py` job tuples and `run_spades`/`run_single_spades` signatures: logger was passed into worker functions (`ProcessPoolExecutor`) but never called there — all logging already happens in the main process; removing it avoids Logger pickling issues
- Fixed `setup_logging` in `pull_ncbi_lineage.py` called with unsupported `log_dir="./"` kwarg; removed to match the interface used by all other scripts
- Fixed `FileNotFoundError` in `fastp_module.py` raised with `%s` placeholder string rather than an f-string: message was never interpolated so the path was invisible in the exception text
- Fixed `blast_round1_parser.py` discarding the `setup_logging` return value: removed module-level `logger = logging.getLogger(__name__)` and replaced with `global logger; logger = setup_logging(...)` inside `main()`, ensuring all module-level functions use the properly configured logger
- Added `debug` level support to `log_and_print()` in `its_fun_tools.py`
- Removed duplicate FAIL-filter assignment in `its_a_summary_compiler.py` (identical `filtered_df = df[...]` appeared twice consecutively)

### Phase 2 — Open
- Write `tutorial.md` content with worked examples

### Phase 3 — Open
- Create `its_fun_2_map/` package directory with `__init__.py`; move scripts in; update imports
- Write `pyproject.toml` with metadata, dependencies, and `console_scripts` entry points
- Declare `parse_fastp_json.R` and `its_decision_making.R` as package data
- Ensure `setup_logging()` is only called inside `main()`
- Define `config.yaml` schema and output directory layout convention
- Write `Snakefile` with one rule per pipeline step
- Write `pipeline.sh` as a thin wrapper over `snakemake`

### Stretch Goals / Future Development
- Add `ThreadPoolExecutor` to `its_primer_binding.py` — extract a `process_sample()` function covering file lookup + seqkit subprocess calls per region; collect results with `as_completed`; `args.threads` controls worker count
- Add `ThreadPoolExecutor` to `UNITEd.py` — collect all sample rows first, then fan out NCBI Entrez calls; limit `max_workers` to 2 without an API key, ~4 with one to respect rate limits; `unite_index` is read-only so safe to share across threads
- Add `--evalue_cutoff` and `--allow_all` to `blast_output_parser.py` to match `blast_round1_parser.py` interface (unifies parameterisation across both BLAST parsing steps; required for a shared workflow config)
- Add `UNITE_DB` environment variable default to `UNITEd.py` (removes need to pass `--unite_db` on every invocation when the path is stable)
- Add CLI entry point to `its_decision_making.R` (only needed if calling it as a direct `Rscript` subprocess outside of `its_a_summary_compiler.py`, e.g. in a Nextflow alternative)

---

## Verification

- `conda env create -f its-fun-2-map.yaml` on a non-NHM machine succeeds
- `pip install -e .` from repo root succeeds; `itsfun-qc --help` prints usage
- `python -c "import its_fun_2_map"` succeeds with no logging side-effects
- `Rscript its_decision_making.R --input_csv test.csv --output_csv out.csv` exits 0
- `snakemake -n --configfile config.yaml` produces the expected DAG without errors
