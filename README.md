# its-fun-2-map
A pipeline for processing fungal genome skims from museum specimens - includes quality control, mapping, assembly, BLAST taxonomic assignment, and validation.

## Contents
1. [Dependencies](#Dependencies)
2. [Process](#Process)
3. [Output Directory Structure](#Output-Directory-Structure)
4. [Configuration](#Configuration)
5. [Authorship](#Authors)

## Dependencies
- Python 3.11+
- BLAST+
- SeqKit
- Samtools
- SPAdes
- BWA
- fastp
- R (with jsonlite package)


## Process 
### Step 1: Quality Control
Quality control processing of raw PE reads using fastp
- **Input:** Raw paired-end sequencing reads
- **Output:** Quality-filtered reads
- **Features:**
  - Adapter trimming
  - Quality filtering
  - Read statistics generation
---

### Step 2: Pseudo-Reference Sequence Retrieval
Retrieves taxonomically relevant reference sequences from the UNITE database using a bespoke python script (UNITEd.py)
- **Input:** Tracking sheet with specimen taxonomy
- **Output:** Reference sequences for each specimen
- **Database:** [UNITE public (Full "UNITE+INSD") database](https://unite.ut.ee/repository.php#panel6a)
- **Parameters:**
  - Taxonomic rank: genus
  - Number of reference sequences: 40
  - Diversity mode enabled
  - Traverses up to family level if required
- **Features:**
  - Generates reference sequence retrieval summary metrics
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
└── combined_summary.csv       # Final aggregated summary file
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
