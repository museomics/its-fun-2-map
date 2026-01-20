# its-fun-2-map

A pipeline for processing fungal genome skims from museum specimens — includes quality control, mapping, assembly, BLAST taxonomic assignment, and validation.

## Contents
1. [Dependencies](#dependencies)
2. [Pipeline Overview](#pipeline-overview)
3. [Detailed Process](#detailed-process)
4. [Workflow](#workflow)
5. [Output Directory Structure](#output-directory-structure)
6. [Configuration](#configuration)
7. [Usage Examples](#usage-examples)
8. [Authors](#authors)

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
- (jsonlite)[https://cran.r-project.org/web/packages/jsonlite/index.html] 

### Custom tool packages (included in repository)
 - its_fun_tools
 - metahist_tools

### Databases
- BLASTn databases:
 - [UNITE Full "UNITE+INSD" database](https://unite.ut.ee/repository.php#panel6a) — for reference sequence retrieval
 - [UNITE General Release (sh_general_release_dynamic)](https://unite.ut.ee/repository.php#panel5a) — for BLAST Round 1
 - [UCHIME ITS1 reference dataset](https://unite.ut.ee/repository.php#panel7a) — for BLAST Round 2b
 - [UCHIME ITS2 reference dataset](https://unite.ut.ee/repository.php#panel7a) — for BLAST Round 2a


## Pipeline Overview

The pipeline processes raw paired-end sequencing reads from fungal museum specimens through quality control, reference-guided read enrichment, assembly, and multi-round BLAST validation to extract and validate ITS barcode sequences.

```
Raw Reads → QC → Reference Retrieval → Read Mapping → Assembly → BLAST Validation → ITS Extraction
```

### Step 1: Quality Control (`fastp_module.py`)

Quality control processing of raw paired-end reads using fastp with a two-stage approach.

### Step 1: Quality Control
Quality control processing of raw PE reads using fastp
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

### Step 2: Pseudo-Reference Sequence Retrieval (`UNITEd.py`)

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
 - Optional summary CSV with retrieval statistics
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







### Step 3: Read Mapping & Baiting
Maps quality-filtered reads to retrieved reference sequences to enrich for target regions using bwa-mem
- **Input:** 
  - Quality-filtered reads from Step 1
  - Reference sequences from Step 2
- **Output:** Mapped reads enriched for ITS regions
- **Aligner:** BWA-MEM
- **Features:**
  - Outputs read mapping summary statistics
---

### Step 4: Contig Assembly
Assembles mapped reads into contigs using SPAdes
- **Input:** 
  - Mapped reads from Step 3
  - Unmerged reads from Step 1
- **Output:** Assembled contigs
- **Features:**
  - Generates assembly summary statistics
---

### Step 5: BLASTn - Round 1 (General Database)
Searches assembled contigs against the [general UNITE database](https://unite.ut.ee/repository.php#panel5a) for initial taxonomic identification.
- **Input:** Assembled contigs from Step 4
- **Output:** Initial BLAST hits against comprehensive fungal database
- **Database:** UNITE general release (sh_general_release_dynamic)
---

### Step 6: BLASTn Round 1 Parsing & Taxonomic Validation
Parses BLAST round 1 output and validates taxonomic assignments against expected taxonomy

- **Input:** BLAST results from Step 5
- **Output:** Filtered and validated BLAST hits
- **Filters:**
  - Minimum sequence length: 100 bp
  - Minimum percent identity: 85%
  - E-value cutoff: 1e-5
- **Features:**
  - Generates BLAST outcome summary
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
