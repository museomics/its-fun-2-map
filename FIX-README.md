# its-fun-2-map: Assessment, Improvements & Packaging Recommendation

## Context

The `its-fun-2-map` repository implements a bioinformatics pipeline for extracting and validating fungal ITS barcode sequences from museum specimen genome skims. It covers 12 steps: QC → Reference Retrieval → BWA Mapping → SPAdes Assembly → BLAST (rounds 1 & 2) → ITS Extraction → Summary. The user requested: (1) a per-script assessment identifying bugs and improvements, and (2) a recommendation on whether to wrap in Snakemake or distribute via PyPI.

---

## Critical Missing Dependencies

1. **`metahist_tools` package is absent from the repository.** Every script imports from it (`clean_and_tar`, `run_command`, `xlsx2csv`, `setup_logging`, `repair_reads`). README claims it is "included in the repository" — it is not. **this will be released as a PyPi package, and therefore will become a dependency of the pipeline.**

2. **`pipeline.sh` and `its_a_summary_compiler.py` are referenced in the README but do not exist.** There is no way to run the full pipeline end-to-end.

---

## Per-Script Assessment

### `its-fun-2-map.yaml` (environment)
- **Bug**: `prefix: /mnt/apps/users/dparsons/conda/envs/its-fun-2-map` — absolute NHM-specific path breaks conda environment creation for any other user.
- **Fix**: Remove the `prefix:` key; conda will use the default env location.

### `fastp_module.py`
- **Bug**: `parse_fastp_json.R` is sourced via a relative path (`rpy2.robjects.r.source("parse_fastp_json.R")`). This breaks if the script is not run from the repository root.
- **Fix**: Derive the R script path relative to `__file__` using `pathlib.Path(__file__).parent / "parse_fastp_json.R"`.
- **Improvement**: Q30 quality threshold is hardcoded; expose as `--quality_threshold` argument.
- **Improvement**: `--threads` controls number of parallel fastp jobs, but the per-job thread count is implicit. Clarify in help text.

### `UNITEd.py`
- **Bug**: `get_ncbi_lineage()` is duplicated here AND in `its_fun_tools.py`. Duplicate should be removed; all callers should use the shared version.
- **Improvement**: Samples are processed sequentially. Adding `ThreadPoolExecutor` (already used elsewhere) would speed up multi-sample runs.
- **Improvement**: UNITE database path is a required argument every run. A default env variable (`UNITE_DB`) would reduce friction.

### `mapping_module.py`
- **Bug**: Per-job BWA thread count is hardcoded as `8` regardless of the `--threads` argument passed. Fix: derive per-job threads as `max(1, args.threads // n_jobs)` or add a dedicated `--bwa_threads` argument.
- **Note in file**: "metahist_tools has an absolute path at the moment until modules are fully packaged" — confirms the packaging debt is known.

### `assembly_module.py`
- No logic bugs identified.
- **Improvement**: `--summary_csv` is optional but its absence silently skips summary output with no message. Add an info log when it is omitted.
- **Improvement**: The 3-stage SPAdes fallback is well-designed; document it in comments or README since it is non-obvious.

### `blast_round1.py` / `blast_round2.py`
- **Bug (`blast_round2.py` line ~40)**: `ids = {line.split("\t")[0] for line in f if line.strip()}` reads a TSV and includes the header row `qseqid` as a sequence ID passed to `seqkit grep`. Fix: skip the first line, or filter out `qseqid`.
- **Inconsistency**: `blast_round1.py` uses `--column_name`; `blast_round1_parser.py` and `blast_output_parser.py` use `--id_column`. Standardise to one name across all scripts.
- **Improvement (`blast_task()` in `its_fun_tools.py`)**: `-evalue 1e-10` is hardcoded. The parsers default to `1e-5`. Expose `--evalue` as an argument and pass it through `blast_task()`.

### `blast_round1_parser.py` / `blast_output_parser.py`
- No logic bugs identified beyond the `--id_column` vs `--column_name` inconsistency noted above.
- **Improvement**: Both parsers silently continue when no BLAST hits exist for a sample. A warning log per missing sample would help debugging.

### `its_primer_binding.py`
- **Bug**: When custom primers are supplied via TSV, the summary CSV fieldnames are still hardcoded as `ITS1`, `ITS2`, `ITS_complete`. This produces misleadingly named columns for non-ITS amplicons.
- **Fix**: Derive fieldnames dynamically from the primer TSV `name` column.
- **Improvement**: No parallelism. Processing is sequential per sample. `ThreadPoolExecutor` would help for large sample sets.
- **Improvement**: Mixes `print()` and `logger` calls inconsistently. Standardise to logger.

### `its_fun_tools.py` (shared utilities)
- **Bug (`get_ncbi_lineage`)**: Uses `urllib.error`, `random`, and `time` without importing them — would raise `NameError` at runtime. Fix: add missing imports.
- **Bug (`load_name_ids`)**: Returns a DataFrame for XLSX input and a list for CSV input — inconsistent return type. Fix: always return a list (apply `.tolist()` in the XLSX branch too).
- **Improvement (`log_and_print`)**: Uses the root logger. Should accept a logger instance as a parameter to be consistent with other functions.
- **Improvement (`blast_task`)**: Hardcoded `-evalue 1e-10` should be a parameter.

### `pull_ncbi_lineage.py`
- **Bug (`main()`)**: References `args.log_file`, `args.input_csv`, etc. before `args = parser.parse_args()` is called. The `parse_args()` call is inside the `if __name__ == "__main__"` block but `main()` is called from there — works at runtime but `main()` has no `args` in its local scope. Fix: pass args explicitly to `main(args)`.
- **Bug (`add_ncbi_lineages_to_csv`)**: Line 44 raises `ValueError(f"...{tracking_sheet}")` but `tracking_sheet` is undefined in that scope (parameter is `input_csv`). Fix: replace with `input_csv`.

### `its_decision_making.R`
- **Bug (line ~59)**: Inside the `s6` scenario block, `df[s2, "Final_contig_desc"]` references `s2` instead of `s6` — copy-paste error that silently updates the wrong rows.
- **Bug (line ~94)**: `df$df$extraction_ITS_complete_path` — double `df$` prefix is a typo. Fix: `df$extraction_ITS_complete_path`.
- **Improvement**: The file only defines `its_outcome()` with no standalone execution entrypoint (`if (!interactive()) { ... }`). It is called from the missing `its_a_summary_compiler.py`. The R function and its Python caller should both exist in the repository.

### `tutorial.md`
- Contains only empty section headers. Needs real usage examples, argument descriptions, and worked examples.

---

## Cross-Cutting Issues

| Issue | Affected files | Fix |
|---|---|---|
| `metahist_tools` absent | All scripts | Include package or inline its functions |
| Missing `pipeline.sh` | README | Write or regenerate the shell orchestration |
| Missing `its_a_summary_compiler.py` | README, R script | Recreate or document its expected interface |
| Conda prefix hardcoded | `its-fun-2-map.yaml` | Remove `prefix:` key |
| Argument name inconsistency (`--column_name` vs `--id_column`) | Multiple | Standardise across all scripts |
| No resume/checkpoint | All | Snakemake handles this; without it a mid-run failure requires restarting from scratch |
| R script sourced with relative path | `fastp_module.py` | Use `Path(__file__).parent` |
| `tutorial.md` empty | Docs | Write content |

---

## Snakemake vs PyPI Recommendation

### Snakemake
**Pros:**
- Native checkpointing — failed steps restart from the last completed step, not from scratch
- Explicit dependency graph makes parallelism automatic and safe
- Resource declarations (threads, memory) per rule replace the current manual thread math
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
- Solves the `metahist_tools` import problem cleanly — bundle it as a submodule
- Works well for the individual-script use case

**Cons:**
- Does not provide orchestration, dependency tracking, or checkpointing
- The user still needs to call each step manually or maintain their own shell script
- Non-Python tools (fastp, BWA, SPAdes, BLAST) are not installable via pip — conda or system packages are still required

### Recommendation: **Both, in combination**

The two approaches are not mutually exclusive and together solve different problems:

1. **PyPI package first**: Package the Python modules (including `metahist_tools` as a sub-package) as `its-fun-2-map`. This fixes the import problem for all scripts immediately and makes installation reproducible. Each script becomes a `console_scripts` entry point (e.g. `itsfun-qc`, `itsfun-map`). The conda YAML continues to manage non-Python tools.

2. **Snakemake wrapper on top**: A `Snakefile` (≈150 lines) wraps each step as a rule using the installed CLI entry points. Users who want push-button end-to-end execution use `snakemake`; users who want to run individual steps use the CLI entry points directly. This is the pattern used by tools like Prokka, BUSCO, and CheckM2.

This delivers:
- `pip install its-fun-2-map` (or `conda install`) for all Python logic
- `snakemake --cores N` for full pipeline runs with checkpointing
- Individual step CLIs for interactive/exploratory use
- No duplication: Snakemake rules call the same entry points as manual use

---

## Implementation Plan

### Phase 1 – Fix critical bugs (all files, no new features)
1. Add `metahist_tools` content to repo (or inline its used functions)
2. Fix conda prefix in YAML
3. Fix `its_fun_tools.py`: add missing imports, fix `load_name_ids` return type
4. Fix `pull_ncbi_lineage.py`: args scoping bug, undefined `tracking_sheet`
5. Fix `blast_round2.py`: skip TSV header when building ID set
6. Fix `its_decision_making.R`: s2→s6 copy-paste, `df$df$` double prefix
7. Fix `fastp_module.py`: resolve R script path relative to `__file__`
8. Fix `its_primer_binding.py`: dynamic fieldnames for custom primer regions
9. Remove duplicate `get_ncbi_lineage` from `UNITEd.py`

### Phase 2 – Usability improvements
10. Standardise `--column_name` / `--id_column` across all scripts
11. Expose `--evalue` in `blast_task()` and thread count in `mapping_module.py`
12. Add parallelism to `its_primer_binding.py` and `UNITEd.py`
13. Standardise logging (remove bare `print()` calls)
14. Write `tutorial.md` content

### Phase 3 – Packaging
15. Create `pyproject.toml` / `setup.cfg` for `its-fun-2-map` PyPI package
16. Bundle `metahist_tools` as `its_fun_2_map.metahist_tools`
17. Register each script as a `console_scripts` entry point
18. Write `Snakefile` with one rule per pipeline step
19. Recreate `pipeline.sh` as a thin wrapper over `snakemake`
20. Write or stub `its_a_summary_compiler.py`

---

## Verification

- `conda env create -f its-fun-2-map.yaml` on a non-NHM machine succeeds
- `python -c "from its_fun_2_map import its_fun_tools; print('ok')"` passes
- `python pull_ncbi_lineage.py --input_csv test.csv ...` runs without NameError
- `python blast_round2.py ...` does not pass `qseqid` to seqkit
- `snakemake -n` (dry-run) on a sample config file produces the expected DAG
- R: `source("its_decision_making.R"); its_outcome(test_df)` returns correct scenario labels for s2 and s6 cases
