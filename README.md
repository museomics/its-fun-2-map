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

### Installing dependencies
All dependencies are included in  `its-2-map-fun.yaml`. A conda environment can be created from this YAML using the following command after [conda](https://www.anaconda.com/docs/getting-started/miniconda/install#quickstart-install-instructions) is installed:
```
conda env create -f its-fun-2-map.yaml
```

## Pipeline Overview
The pipeline processes raw paired-end sequencing reads from fungal museum specimens through quality control, reference-guided read enrichment, assembly, and multi-round BLAST validation to extract and validate ITS barcode sequences.

```
Raw Reads → QC → Reference Retrieval → Read Mapping → Assembly → BLAST Validation → ITS Extraction
```

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


### Step 7: BLASTn - Round 2 (ITS2-Specific)
Searches contigs against [UCHIME ITS2 reference database](https://unite.ut.ee/repository.php#panel7a) for region-specific validation
- **Input:** Assembled contigs from Step 4
- **Output:** ITS2-specific BLAST results
- **Database:** UCHIME ITS2 reference dataset
---

### Step 8: BLASTn Parsing & Validation (ITS2)
Parses ITS2-specific BLAST results and applies quality filters
- **Input:** ITS2 BLAST results from Step 7
- **Output:** Validated ITS2 sequences with taxonomic confirmation
- **Filters:** Same as Step 6
- **Features:**
  - Generates BLAST outcome summary
---

### Step 9: BLASTn - Round 2B (ITS1-Specific)
Parallel analysis searching against [UCHIME ITS1 reference database](https://unite.ut.ee/repository.php#panel7a) for region-specific validation
- **Input:** Assembled contigs from Step 4
- **Output:** ITS1-specific BLAST results
- **Database:** UCHIME ITS1 reference dataset
---

### Step 10: BLASTn Parsing & Validation (ITS1)
Parses ITS1-specific BLAST results and applies quality filters.
- **Input:** ITS1 BLAST results from Step 9
- **Output:** Validated ITS1 sequences with taxonomic confirmation
- **Filters:** Same as Step 6
- **Features:**
  - Generates BLAST outcome summary
---

### Step 11: ITS Primer Binding & Extraction
Identifies primer binding sites and extracts ITS sequences from validated contigs using seqkit amplicon
- **Input:** Validated ITS2 contigs from Step 8
- **Output:** Extracted ITS sequences with primer coordinates
- **Features:**
  - Generates ITS primer-based extraction outcome summary
---

### Step 12: Summary meitrics aggregation
Compiles comprehensive summary statistics across all pipeline steps
- **Input:** Results from all previous steps
- **Output:** Integrated metrics summary
- **Metrics include:**
  - Quality control statistics (fastp)
  - Reference sequence retrieval (UNITEd)
  - Mapping rates
  - Assembly metrics
  - BLAST validation results (ITS1 and ITS2)
  - Extraction success rates

## Workflow
<img width="1278" height="1222" alt="image" src="https://github.com/user-attachments/assets/44ada96d-73f3-4545-8290-14239e7d2b3e" />

## Output Directory Structure
```
output_base/
├── fastp_processed/           # Quality-filtered reads
├── UNITEd/                    # Retrieved reference sequences
├── mapped_reads/              # Reads mapped to references
├── assemblies/                # Assembled contigs
├── blast_round1/              # Initial BLAST results (general database)
├── blast_parsed1/             # Filtered and validated round 1 results
├── blast_round2a-ITS2/        # ITS2-specific BLAST results
├── blast_parsed2a-ITS2/       # Validated ITS2 results
├── blast_round2b-ITS1/        # ITS1-specific BLAST results
├── blast_parsed2b-ITS1/       # Validated ITS1 results
├── its_primer_extraction/     # Extracted ITS sequences
└── logs/                      # Detailed logs for each step
└── final_results_dir          # Contains a final aggregated summary csv and final FASTAs
```

## Configuration
Key parameters can be adjusted in the pipeline script:
- `MIN_LENGTH`: Minimum sequence length for filtering (default: 100 bp)
- `MIN_PIDENT`: Minimum percent identity (default: 85%)
- `EVALUE_CUTOFF`: BLAST e-value threshold (default: 1e-5)
- `NUMBER_REFS`: Number of reference sequences to retrieve (default: 40)
- `TAX_RANK`: Taxonomic rank for reference retrieval (default: genus)


## Authors
**M. Kamouyiaros & D. Parsons** @ the Natural History Museum, London (NHMUK)
