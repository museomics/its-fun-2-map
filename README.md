# its-fun-2-map

A pipeline for processing fungal genome skims from museum specimens — includes quality control, mapping, assembly, BLAST taxonomic assignment, and validation.

## Contents
1. [Dependencies](#dependencies)
2. [Pipeline Overview](#pipeline-overview)
3. [Workflow](#workflow)
4. [Output Directory Structure](#output-directory-structure)
5. [Configuration](#configuration)
6. [Usage Examples](#usage-examples)
7. [Authors](#authors)

## Dependencies
### Languages
- r-base (4.3.0 (2023-04-21 ucrt) -- "Already Tomorrow" +) 
- python 3.11+

### Conda packages
| Package  | Version              |
|----------|----------------------|
| spades   | 4.0.0                |
| seqkit   | 2.2.0                |
| samtools | 1.21                 |
| bwa      | 0.7.18-r1243-dirty   |
| fastp    | 0.24.0               |
| htslib   | 1.21                 |
| blast    | 2.2.31+              |

### Python packages (not built-in)
- pandas (`pip install pandas`)
- Bio (Biopython) (`pip install biopython`)
- rpy2.robjects (`pip install rpy2`)
- pgzip (`pip install pgzip`)

### R packages
- [jsonlite](https://cran.r-project.org/web/packages/jsonlite/index.html)

### Custom tool packages (included in repository)
 - its_fun_tools
 - metahist_tools

### Databases
- BLASTn databases:
  - [UNITE Full "UNITE+INSD" database](https://unite.ut.ee/repository.php#panel6a) — for reference sequence retrieval
  - [UNITE General Release (sh_general_release_dynamic)](https://unite.ut.ee/repository.php#panel5a) — for BLAST Round 1
  - [UCHIME ITS1 reference dataset](https://unite.ut.ee/repository.php#panel7a) — for BLAST Round 2b
  - [UCHIME ITS2 reference dataset](https://unite.ut.ee/repository.php#panel7a) — for BLAST Round 2a

## Installing dependencies
All dependencies are included in  `its-2-map-fun.yaml`. A conda environment can be created from this YAML using the following command after [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions) is installed:
```
conda env create -f its-fun-2-map.yaml
```

## Pipeline Overview
The pipeline processes raw paired-end sequencing reads from fungal museum specimens through quality control, reference-guided read enrichment, assembly, and multi-round BLAST validation to extract and validate ITS barcode sequences.

```
Raw Reads → QC → Reference Retrieval → Read Mapping → Assembly → BLAST Validation → ITS Extraction
```
**All pipeline steps can be run sequentially using the `pipeline.sh` script.**
---


## Step 1: Quality Control (`fastp_module.py`)

Quality control processing of raw paired-end reads using fastp with a two-stage approach.

- **Input:** Raw paired-end sequencing reads
- **Output:**
  -  `{sample}_trimmed_1.fq`, `{sample}_trimmed_2.fq` — trimmed paired reads
  - `{sample}_merged.fq` — merged overlapping reads
  - `{sample}_unmerged_1.fq`, `{sample}_unmerged_2.fq` — reads that couldn't be merged
  - `{sample}_trim.json`, `{sample}_merge.json` — QC metrics
  - `{sample}_overlaps.html` — overlap distribution visualisation
  - `fastp_summary.csv` — aggregated summary statistics
- **Features:**
  - Adapter trimming and filtering (auto-detection for PE reads)
  - Quality filtering (Q30 threshold)
  - Read merging (Minimum merged length: 30 bp)
  - Read statistics generation

---

## Step 2: Pseudo-Reference Sequence Retrieval (`UNITEd.py`)
Retrieves taxonomically relevant reference sequences from the UNITE database using NCBI taxonomy matching.

**Process:**
1. Queries NCBI Entrez API for taxonomic lineage of each specimen
2. Matches NCBI taxonomy against UNITE database taxonomy
3. Extracts sequences based on user-specified parameters
4. Supports fallback to higher taxonomic ranks if matches not found

**Key Features:**
 - **Diversity mode:** Distributes sequences across child taxa, prioritising target species matches (up to 25% of requested sequences)
 - **Traverse mode:** Searches up the taxonomic hierarchy to collect target number of sequences
 - **Flexible matching:** Case-insensitive matching with normalisation for species name variations

- **Input:** Tracking sheet with specimen taxonomy
- **Output:**
  - `{sample_id}_seed.fasta` — reference sequences for each specimen
  - Summary CSV with retrieval statistics
- **Database:** [UNITE public (Full "UNITE+INSD") database](https://unite.ut.ee/repository.php#panel6a)

**Parameters:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--tax_rank` | — | Taxonomic rank to search (species/genus/family/order/class/phylum) |
| `--number` | — | Number of reference sequences to retrieve |
| `--diversity` | False | Distribute sequences across child taxa |
| `--traverse` | — | Maximum rank to traverse up to if insufficient sequences found |
| `--all` | False | Retrieve all matching sequences |

---

## Step 3: Read Mapping & Baiting (`mapping_module.py`)
Maps quality-filtered reads to retrieved reference sequences using BWA to enrich for target ITS regions. Supports both modern DNA (BWA-MEM) and ancient/historical DNA (BWA-ALN) alignment strategies.

**Alignment Algorithms:**
| Algorithm | Use Case | Parameters |
|-----------|----------|------------|
| `bwa-mem` | Modern/historical DNA | Default BWA-MEM settings with `-M` flag |
| `bwa-aln` | Ancient/historical DNA | `-l 16500 -n 0.01 -o 2` (Green et al. 2010) |

**Process:**
1. Automatically indexes reference FASTA files if not already indexed
2. Maps merged reads (single-end mode)
3. Maps unmerged read pairs (paired-end mode)
4. Extracts mapped reads to FASTQ format
5. Repairs paired-end read files to ensure synchronisation
6. Generates mapping statistics

**Input:**
  - Quality-filtered FASTQ files from Step 1 (`*_merged.fq`, `*_unmerged_1.fq`, `*_unmerged_2.fq`)
  - Per-sample reference FASTA files from Step 2 (`*_seed.fasta`)

**Output:**
  - `{sample}_mapped.fastq` — mapped merged reads
  - `{sample}_mapped_unmerged_1.fastq`, `{sample}_mapped_unmerged_2.fastq` — mapped paired reads (repaired)
  - `{sample}_mapped_flagstats.txt` — alignment statistics for merged reads
  - `{sample}_unmerged_1_flagstats.txt`, `{sample}_unmerged_2_flagstats.txt` — alignment statistics for paired reads
  - `mapping_summary.csv` — aggregated mapping metrics for all samples

---

## Step 4: Contig Assembly (`assembly_module.py`)
Assembles mapped reads into contigs using SPAdes with a multi-stage fallback strategy to maximise assembly success.

**Assembly Strategy:**
| Stage | Method | K-mer Sizes | Description |
|-------|--------|-------------|-------------|
| Initial | Single-reads (i.e. merged) | 21,33,55 | Uses mapped merged reads only |
| Fallback 1 | Single-reads (i.e. merged) | 21 | Retries with minimal k-mer if initial fails |
| Fallback 2 | Merged + Unmerged PE | 21 | Uses all available reads if Fallback 1 fails |

**Process:**
1. Attempts assembly with merged reads using default k-mer sizes
2. If assembly fails or produces no valid scaffolds, retries with k=21 only
3. If still unsuccessful, includes unmerged paired reads in the assembly
4. Calculates assembly metrics including N50, scaffold counts, and lengths

**Input:**
  - Mapped merged reads from Step 3 (`*_mapped.fastq`)
  - Original unmerged reads from Step 1 (`*_unmerged_1.fq`, `*_unmerged_2.fq`)

**Output:**
- `{sample}.spades.out/` — SPAdes output directory containing:
  - `scaffolds.fasta` — assembled scaffolds
  - Assembly graphs and logs
- `assembly_summary.csv` — metrics for all samples including:
  - Assembly method used (INITIAL/FALLBACK1/FALLBACK2)
  - Number of scaffolds
  - Scaffold lengths
  - Mean scaffold length
  - N50 value
 
---

## Step 5: BLASTn — Round 1 (General Database) (`blast_pipe.py`)
Searches assembled scaffolds against the general UNITE database for initial taxonomic identification.

**Process:**
1. Optionally creates BLAST database from reference FASTA
2. Locates `scaffolds.fasta` files for each sample in assembly directories
3. Runs BLASTn searches in parallel
4. Outputs results in tabular format (outfmt 6)

**Input:**
- Assembled scaffolds from Step 4 (`*.spades.out/scaffolds.fasta`)
- UNITE general release database
- Tracking sheet with sample IDs

**Output:**
- `{sample}_blast.tsv` — BLAST results in tabular format with columns:
  - qseqid, sseqid, pident, length, mismatch, gapopen, qstart, qend, sstart, send, evalue, bitscore

**Database:** UNITE general release (sh_general_release_dynamic)

---

## Step 6: BLASTn Round 1 Parsing & Taxonomic Validation (`blast_round1_parser.py`)
Parses BLAST Round 1 output and validates taxonomic assignments against expected taxonomy with cascading taxonomic rank matching.

**Taxonomic Validation Strategy:**
The parser attempts to match BLAST hits against expected taxonomy using a cascading approach:
1. **Species-level match** — checks if hits match expected species
2. **Genus-level match** — if no species match, checks genus
3. **Family-level match** — if no genus match, checks family
4. If no matches at any level, the sample is marked as FAIL

**Filtering Criteria:**
| Parameter | Default | Description |
|-----------|---------|-------------|
| `--min_len` | 100 bp | Minimum alignment length |
| `--min_pident` | 85% | Minimum percent identity |
| `--evalue_cutoff` | 1e-5 | Maximum e-value threshold |

**Validation Outcomes:**
- **PASS**: Exactly one contig matches expected taxonomy after filtering
- **FAIL**: Zero contigs match, or multiple contigs remain after filtering

**Input:**
- BLAST TSV files from Step 5
- Taxonomy CSV with expected taxonomic assignments (must contain ID column and taxonomic hierarchy: Kingdom, Phylum, Class, Order, Family, Genus, Species)
- Assembly directory for contig extraction

**Output:**
- `blast_validation_summary.csv` — comprehensive validation results including:
  - `ID` — sample identifier
  - `blast_round1_found_taxon` — taxonomic level where match was found
  - `blast_round1_expected_taxonomy` — full expected taxonomy string
  - `blast_round1_n_contigs_in` — total contigs in assembly
  - `blast_round1_n_contigs_hits` — contigs matching expected taxonomy
  - `blast_round1_contigs` — contig IDs that matched
  - `blast_round1_hit_taxonomy` — taxonomy of best BLAST hit
  - `blast_round1_hit_coords` — alignment coordinates
  - `blast_round1_correct_taxonomy` — PASS/FAIL status
- `filtered_{sample}_blast.tsv` — filtered BLAST results for passing samples
- `{sample}_parsed_contig.fasta` — extracted contig sequences for passing samples
---

## Step 7: BLASTn — Round 2a (ITS2-Specific) (`blast_round2.py`)
Searches validated contigs against the UCHIME ITS2 reference database for region-specific validation. This step uses `seqkit grep` to extract only the contigs that passed Round 1 validation before running BLAST.

**Process:**
1. Reads BLAST Round 1 TSV results to identify validated contig IDs
2. Extracts unique contig headers from each sample's Round 1 results
3. Uses `seqkit grep` to extract matching scaffolds from assembly FASTA files
4. Runs BLASTn against ITS2-specific UCHIME database
5. Processes samples in parallel

**Input:**
- Assembly scaffolds from Step 4 (`*.spades.out/scaffolds.fasta`)
- BLAST Round 1 TSV results from Step 5/6
- UCHIME ITS2 reference database
- Tracking sheet with sample IDs

**Output:**
- `{sample}_blast.tsv` — ITS2-specific BLAST results
- Temporary matched scaffold FASTAs (cleaned up after processing)

**Database:** UCHIME ITS2 reference dataset

---

## Step 8: BLASTn Parsing & Validation (ITS2) (`blast_output_parser.py`)
Parses ITS2-specific BLAST results and applies multi-tier taxonomic validation with cascading fallbacks from family to genus to species level.

**Multi-Tier Taxonomy Matching Logic:**
When expected family is present:
1. **0 family matches** → FAIL
2. **1 family match** → PASS
3. **>1 family matches** → Check genus level:
   - No expected genus → FAIL
   - 1 genus match → PASS
   - 0 or >1 genus matches → Check species level:
     - No expected species → FAIL
     - 1 species match → PASS
     - 0 or >1 species matches → FAIL

When expected family is empty (fallback to genus):
1. **0 genus matches** → FAIL
2. **1 genus match** → PASS
3. **>1 genus matches** → Check species level (same logic as above)

**Input:**
- ITS2 BLAST TSV files from Step 7
- Taxonomy CSV with expected taxonomic assignments
- Assembly directory for contig extraction

**Output:**
- `taxonomy_validation_summary.csv` — validation results including:
  - `ID` — sample identifier
  - `expected_family` — expected taxonomic family (or genus/species if family unavailable)
  - `expected_taxonomy` — full expected taxonomy string
  - `n_contigs_in` — total contigs in assembly
  - `n_contigs_hits` — contigs matching expected taxonomy
  - `contigs` — matching contig IDs
  - `hit_taxonomy` — taxonomy from best BLAST hit
  - `its_coords` — ITS region coordinates on contig
  - `correct_taxonomy` — PASS/FAIL status
  - `contig_path` — path to extracted contig FASTA
- `{sample}_parsed_contig.fasta` — extracted contig sequences for passing samples

---

## Step 9: BLASTn — Round 2b (ITS1-Specific) (`blast_round2.py`)
Parallel analysis searching against UCHIME ITS1 reference database. Uses the same workflow as Step 7 but with the ITS1-specific database.

**Input:**
- Assembly scaffolds from Step 4
- BLAST Round 1 TSV results from Step 5/6
- UCHIME ITS1 reference database
- Tracking sheet with sample IDs

**Output:**
- `{sample}_blast.tsv` — ITS1-specific BLAST results

**Database:** UCHIME ITS1 reference dataset

---

## Step 10: BLASTn Parsing & Validation (ITS1) (`blast_output_parser.py`)
Parses ITS1-specific BLAST results using the same multi-tier taxonomic validation as Step 8.

**Input:**
- ITS1 BLAST TSV files from Step 9
- Taxonomy CSV with expected taxonomic assignments
- Assembly directory for contig extraction

**Output:**
- `taxonomy_validation_summary.csv` — ITS1 validation results (same format as Step 8)
- `{sample}_parsed_contig.fasta` — extracted contig sequences for passing samples

---

### Step 11: ITS Primer Binding & Extraction (`its_primer_binding.py`)
Identifies primer binding sites and extracts ITS sequences from validated contigs using `seqkit amplicon`. Uses standard ITS primers from White et al. (1990) by default, with support for custom primer sets.

**Default Primers (White et al. 1990):**
| Primer | Sequence (5' → 3') |
|--------|-------------------|
| ITS1 | `TCCGTAGGTGAACCTGCGG` |
| ITS2 | `GCTGCGTTCTTCATCGATGC` |
| ITS3 | `GCATCGATGAAGAACGCAGC` |
| ITS4 | `TCCTCCGCTTATTGATATGC` |

**Target Regions:**
| Region | Forward Primer | Reverse Primer | Description |
|--------|---------------|----------------|-------------|
| ITS1 | ITS1 | ITS2 | ITS1 spacer region |
| ITS2 | ITS3 | ITS4 | ITS2 spacer region |
| ITS_complete | ITS1 | ITS4 | Full ITS region (ITS1 + 5.8S + ITS2) |

**Process:**
1. Reads sample IDs from tracking sheet
2. Locates FASTA files for each sample
3. Runs `seqkit amplicon` with each primer pair (max 2 mismatches, strict mode)
4. Generates extraction summary statistics

**Input:**
- Validated contig FASTAs from Steps 8/10 (`*_parsed_contig.fasta`)
- Tracking sheet with sample IDs
- Optional: Custom primers TSV and regions TSV

**Output:**
- `ITS1/{sample}_ITS1.fa` — extracted ITS1 sequences
- `ITS2/{sample}_ITS2.fa` — extracted ITS2 sequences
- `ITS_complete/{sample}_ITS_complete.fa` — extracted complete ITS sequences
- `extraction_summary.csv` — extraction success per sample and region
- `summary_report.txt` — detailed extraction statistics

**Custom Primers TSV Format:**
```
name	sequence
FWD1	ATCGATCGATCG
REV1	GCTAGCTAGCTA
```

**Custom Regions TSV Format:**
```
region	forward	reverse
MyRegion	FWD1	REV1
```

---

## Step 12: Summary Metrics Aggregation (`merge_outputs.py`)
Compiles comprehensive summary statistics across all pipeline steps, performs contig analysis to resolve ITS1/ITS2 results, and prepares final output FASTAs for submission.

**Process:**
1. **CSV Discovery:** Finds and merges all pipeline summary CSVs (fastp, UNITEd, mapping, assembly, BLAST, extraction)
2. **Contig Analysis:** Compares ITS1 and ITS2 results to identify:
   - Contigs with both ITS1 and ITS2 on the same scaffold
   - Overlapping coordinate detection (potential assembly issues)
   - Resolution when ITS1 and ITS2 are on different contigs
3. **Decision Making:** Runs R-based decision logic (`its_decision_making.R`) to determine final outcomes
4. **FASTA Preparation:** Renames and organises final FASTAs for submission with appropriate headers

**Contig Analysis Logic:**
- **Both ITS1 and ITS2 PASS on same contig:** Check for coordinate overlap
  - Non-overlapping → PASS (contiguous)
  - Overlapping → WARNING (potential issue)
- **Both PASS on different contigs:** Default to ITS2 contig
- **Only ITS1 PASS:** Use ITS1 contig
- **Only ITS2 PASS:** Use ITS2 contig
- **Both FAIL:** Mark as failed

**Input:**
- All pipeline summary CSV files from previous steps
- Optional: Naming TSV for custom FASTA header formatting

**Output:**
- `final_results_dir/merged_summary_{date}.csv` — comprehensive merged summary with columns:
  - UNITEd metrics (reference retrieval)
  - Trimming/merging metrics (fastp)
  - Assembly metrics (SPAdes)
  - BLAST Round 1 results
  - ITS1 and ITS2 BLAST results
  - Extraction results
  - Mapping statistics
  - Contig analysis results (overlapping_coords, same_contig)
  - Final decision (Final_outcome, Final_contig_desc, Final_contig)
- `final_results_dir/pass_fastas/` — renamed FASTAs for samples that passed all validation
- `final_results_dir/manual_verification/` — FASTAs requiring manual review

**Naming TSV Format (optional):**
```
SAMPLE001	Genus_species_voucher123
SAMPLE002	Genus_species_voucher456
```

---

## Workflow
<img width="1278" height="1222" alt="image" src="https://github.com/user-attachments/assets/44ada96d-73f3-4545-8290-14239e7d2b3e" />

## Output Directory Structure
```
output_base/
├── fastp_processed/               # Step 1: Quality-filtered reads
│   ├── {sample}_trimmed_1.fq
│   ├── {sample}_trimmed_2.fq
│   ├── {sample}_merged.fq
│   ├── {sample}_unmerged_1.fq
│   ├── {sample}_unmerged_2.fq
│   ├── {sample}_trim.json
│   ├── {sample}_merge.json
│   ├── {sample}_overlaps.html
│   └── fastp_summary.csv
├── UNITEd/                        # Step 2: Retrieved reference sequences
│   ├── {sample}_seed.fasta
│   └── retrieval_summary.csv
├── mapped_reads/                  # Step 3: Reads mapped to references
│   ├── {sample}_mapped.fastq
│   ├── {sample}_mapped_unmerged_1.fastq
│   ├── {sample}_mapped_unmerged_2.fastq
│   ├── {sample}_mapped_flagstats.txt
│   ├── {sample}_unmerged_1_flagstats.txt
│   ├── {sample}_unmerged_2_flagstats.txt
│   └── mapping_summary.csv
├── assemblies/                    # Step 4: Assembled contigs
│   ├── {sample}.spades.out/
│   │   ├── scaffolds.fasta
│   │   ├── contigs.fasta
│   │   └── ...
│   └── assembly_summary.csv
├── blast_round1/                  # Step 5: Initial BLAST results
│   └── {sample}_blast.tsv
├── blast_parsed1/                 # Step 6: Filtered and validated results
│   ├── blast_validation_summary.csv
│   ├── filtered_{sample}_blast.tsv
│   └── {sample}_parsed_contig.fasta
├── blast_round2a-ITS2/            # Step 7: ITS2-specific BLAST results
│   └── {sample}_blast.tsv
├── blast_parsed2a-ITS2/           # Step 8: Validated ITS2 results
│   ├── taxonomy_validation_summary.csv
│   └── {sample}_parsed_contig.fasta
├── blast_round2b-ITS1/            # Step 9: ITS1-specific BLAST results
│   └── {sample}_blast.tsv
├── blast_parsed2b-ITS1/           # Step 10: Validated ITS1 results
│   ├── taxonomy_validation_summary.csv
│   └── {sample}_parsed_contig.fasta
├── its_extraction/                # Step 11: Extracted ITS sequences
│   ├── ITS1/
│   │   └── {sample}_ITS1.fa
│   ├── ITS2/
│   │   └── {sample}_ITS2.fa
│   ├── ITS_complete/
│   │   └── {sample}_ITS_complete.fa
│   ├── extraction_summary.csv
│   └── summary_report.txt
├── logs/                          # Detailed logs for each step
└── final_results_dir/             # Step 12: Final aggregated outputs
    ├── merged_summary_{date}.csv
    ├── pass_fastas/
    │   └── {sample}_{description}.fasta
    └── manual_verification/
        └── {sample}_manual_verification_needed.fasta
```


## Configuration
### Default Parameters

| Parameter | Default | Description |
|-----------|---------|-------------|
| `MIN_LENGTH` | 100 bp | Minimum sequence length for BLAST filtering |
| `MIN_PIDENT` | 85% | Minimum percent identity for BLAST hits |
| `EVALUE_CUTOFF` | 1e-5 | BLAST e-value threshold |
| `NUMBER_REFS` | 40 | Number of reference sequences to retrieve |
| `TAX_RANK` | genus | Taxonomic rank for reference retrieval |

### Tracking Sheet Format
The pipeline accepts CSV or XLSX tracking sheets with the following columns:
**Required columns:**
- Sample ID column (configurable via `--column_name` or `--id_header`)
- Taxid column (for UNITEd step, configurable via `--taxid_header`)
**Optional columns (for fastp_module.py):**
- `forward` / `fwd` — path to forward reads
- `reverse` / `rev` — path to reverse reads
**For taxonomy validation (blast parsers):**
- `Kingdom`, `Phylum`, `Class`, `Order`, `Family`, `Genus`, `Species` — taxonomic hierarchy columns

If forward/reverse columns are present, file paths are read directly from the tracking sheet. Otherwise, files are searched for in the input directory using sample IDs.

## Authors
**M. Kamouyiaros & D. Parsons** @ the Natural History Museum, London (NHMUK)
