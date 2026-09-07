#!/bin/bash
#SBATCH --job-name=ITS-FUN-2-MAP
#SBATCH --mem=150G
#SBATCH --cpus-per-task=16
#SBATCH --partition=[PLACEHOLDER]


#===============================================================================
# ITS-FUN-2-MAP
#===============================================================================
#
# Pipeline for extracting the ITS region from fungal genome skims
# from museum specimens. Includes quality control, mapping, assembly,
# BLAST, and summary.
#
# Author: M. Kamouyiaros & D. Parsons (NHMUK)
# Version: 1.0.0
#
#===============================================================================


#===============================================================================
# CONFIGURATION
#===============================================================================

# Conda environment activation
# Source conda.sh
#source /mnt/apps/users/$USER/conda/etc/profile.d/conda.sh
# Activate env
conda activate its-fun-2-map



#### USER-SPECIFIED VARIABLES/PARAMETERS/PATHS
# Main paths and variables
INPUT_DIR="[PLACEHOLDER]"
TRACKING_SHEET="[PLACEHOLDER]"
OUTPUT_BASE="[PLACEHOLDER]"
COLUMN_NAME="[PLACEHOLDER]" # E.g. "ID"
TAX_COLUMN="[PLACEHOLDER]" # E.g. "taxid"

# UNITEd-specific parameters
need_lineage="NO"
UNITE_DB_PATH="[PLACEHOLDER]"
EMAIL="[PLACEHOLDER]"
API_KEY="[PLACEHOLDER]"

# Fastp parameters
spreadsheet_filepaths="YES"


#### DEFAULTS - Change for advanced use only
## Output directories - CONSTANTS
LOGS="${OUTPUT_BASE}/00_logs"
FASTP_OUTPUT="${OUTPUT_BASE}/01_fastp_processed"
UNITED_OUTPUT="${OUTPUT_BASE}/02_UNITEd"
MAPPING_OUTPUT="${OUTPUT_BASE}/03_mapped_reads"
ASSEMBLY_OUTPUT="${OUTPUT_BASE}/04_assemblies"
BLAST_ROUND1_OUTPUT="${OUTPUT_BASE}/05_blast_round1"
BLAST_ROUND2A_OUTPUT="${OUTPUT_BASE}/06_blast_round2a-ITS2"
BLAST_ROUND2B_OUTPUT="${OUTPUT_BASE}/07_blast_round2b-ITS1"
BLAST_PARSED_OUTPUT_ROUND1="${OUTPUT_BASE}/05a_blast_parsed1"
BLAST_PARSED_OUTPUT_ITS2="${OUTPUT_BASE}/06a_blast_parsed2a-ITS2"
BLAST_PARSED_OUTPUT_ITS1="${OUTPUT_BASE}/07a_blast_parsed2b-ITS1"
ITS_EXTRACTION_OUTPUT="${OUTPUT_BASE}/08_its_primer_extraction"
# Created by its_a_summary_compiler.py, not by create_dir below
FINAL_RESULTS_OUTPUT="${OUTPUT_BASE}/final_results_dir"

## Configurable parameters
# UNITEd.py parameters
NUMBER_REFS="20"
TAX_RANK="genus"
UNITED_SUMMARY="UNITEd_summary.csv"

# Mapping parameters
ALIGNER="bwa-mem"

# Assembly parameters
ASSEMBLY_CSV="${ASSEMBLY_OUTPUT}/assembly_summary.csv"

# BLAST database paths
BLAST_DB_GENERAL="[PLACEHOLDER]"
BLAST_DB_ITS2="[PLACEHOLDER]"
BLAST_DB_ITS1="[PLACEHOLDER]"

# BLAST output parser parameters
MIN_LENGTH="100"
MIN_PIDENT="85"
PARSER_ITS2_SUMMARY="${BLAST_PARSED_OUTPUT_ITS2}/blast_parser_its2_summary.csv"
PARSER_ITS1_SUMMARY="${BLAST_PARSED_OUTPUT_ITS1}/blast_parser_its1_summary.csv"

# Summary parameters (optional)
different_naming="NO"
#NAMING_TSV="./its-fun-2-map/naming_test.tsv"


#===============================================================================
# LOGGING SETUP
#===============================================================================

# Create log file with timestamp
LOG_FILE="${OUTPUT_BASE}_pipeline_$(date +"%Y%m%d_%H%M%S").log"
exec > >(tee -a "${LOG_FILE}")
exec 2>&1

# Logging functions
log_info() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [INFO] $1"
}

log_warn() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [WARN] $1"
}

log_error() {
    echo "[$(date +'%Y-%m-%d %H:%M:%S')] [ERROR] $1"
}

log_section() {
    echo ""
    echo "==============================================================================="
    echo " $1"
    echo "==============================================================================="
}

#===============================================================================
# UTILITY FUNCTIONS
#===============================================================================

# Function to check if command succeeded
check_status() {
    local rc=$?
    if [ $rc -eq 0 ]; then
        log_info "$1 completed successfully"
    else
        log_error "$1 failed with exit code $rc"
        exit $rc
    fi
}

# Function to create directory if it doesn't exist
create_dir() {
    if [ ! -d "$1" ]; then
        mkdir -p "$1"
        log_info "Created directory: $1"
    else
        log_info "Directory already exists: $1"
    fi
}

# Function to check if file exists
check_file() {
    if [ ! -f "$1" ]; then
        log_error "Required file not found: $1"
        exit 1
    else
        log_info "Found required file: $1"
    fi
}

#===============================================================================
# PIPELINE INITIALISATION
#===============================================================================

log_section "PIPELINE SET UP"

log_info "ITS-FUN-2-MAP Pipeline started"
log_info "Job ID: ${SLURM_JOB_ID:-LOCAL}"
log_info "Node: ${SLURMD_NODENAME:-$(hostname)}"
log_info "Working directory: $(pwd)"
log_info "Log file: ${LOG_FILE}"

log_info "Conda env: ${CONDA_PREFIX}"
log_info "which python = $(which python)"
log_info "python version = $(python --version 2>&1)"

# Check for required files
log_info "Checking required input files..."
check_file "${TRACKING_SHEET}"
check_file "${UNITE_DB_PATH}"

# BLAST databases
log_info "Checking BLAST databases..."
check_file "${BLAST_DB_GENERAL}"
check_file "${BLAST_DB_ITS2}"
check_file "${BLAST_DB_ITS1}"

# Create output base directory
log_info "Creating output directories..."
create_dir "${OUTPUT_BASE}"
create_dir "${LOGS}"
create_dir "${FASTP_OUTPUT}"
create_dir "${UNITED_OUTPUT}"
create_dir "${MAPPING_OUTPUT}"
create_dir "${ASSEMBLY_OUTPUT}"
create_dir "${BLAST_ROUND1_OUTPUT}"
create_dir "${BLAST_ROUND2A_OUTPUT}"
create_dir "${BLAST_ROUND2B_OUTPUT}"
create_dir "${BLAST_PARSED_OUTPUT_ROUND1}"
create_dir "${BLAST_PARSED_OUTPUT_ITS2}"
create_dir "${BLAST_PARSED_OUTPUT_ITS1}"
create_dir "${ITS_EXTRACTION_OUTPUT}"

#===============================================================================
# PIPELINE EXECUTION
#===============================================================================
# Step 1: Quality Control with fastp
log_section "STEP 1: QUALITY CONTROL (FASTP)"
log_info "Starting fastp quality control processing..."
log_info "Input directory: ${INPUT_DIR}"
log_info "Output directory: ${FASTP_OUTPUT}"

if [ "$spreadsheet_filepaths" = "YES" ]; then
   log_info "Filepaths specified in ${TRACKING_SHEET}"
   python fastp_module.py \
       --output "${FASTP_OUTPUT}" \
       --tracking_sheet "${TRACKING_SHEET}" \
       --column_name "${COLUMN_NAME}" \
       --fastp_threads 4 \
       --threads 4 \
       --log_file "${LOGS}/fastp.log"
   check_status "fastp quality control"
fi

if [ "$spreadsheet_filepaths" = "NO" ]; then
   log_info "Filepaths not specified; using input directory ${INPUT_DIR}"
   python fastp_module.py \
       --output "${FASTP_OUTPUT}" \
       --tracking_sheet "${TRACKING_SHEET}" \
       --column_name "${COLUMN_NAME}" \
       --fastp_threads 4 \
       --threads 4 \
       --log_file "${LOGS}/fastp.log" \
       --input_dir "${INPUT_DIR}"
   check_status "fastp quality control"
fi


#===============================================================================

# Step 2: UNITEd.py Analysis
log_section "STEP 2: UNITEd PSEUDO-REFERENCE SEQUENCE RETRIEVAL"
log_info "Starting UNITEd.py search..."
log_info "Database: ${UNITE_DB_PATH}"
log_info "Taxonomic rank: ${TAX_RANK}"
log_info "Number of sequences: ${NUMBER_REFS}"

if [ "$need_lineage" = "YES" ]; then
   python pull_ncbi_lineage.py \
       --input_csv "${TRACKING_SHEET}" \
       --output_csv "${UNITED_OUTPUT}/tracking_with_lineage.csv" \
       --email "${EMAIL}" \
       --api_key "${API_KEY}" \
       --taxcolumn "${TAX_COLUMN}" \
       --log_file "${LOGS}/pull_ncbi_lineage.log"
   check_status "pull_ncbi_lineage"
   TRACKING_SHEET="${UNITED_OUTPUT}/tracking_with_lineage.csv"
   log_info "Using lineage-enriched tracking sheet: ${TRACKING_SHEET}"
fi


python UNITEd.py \
   --tracking_sheet "${TRACKING_SHEET}" \
   --unite_db "${UNITE_DB_PATH}" \
   --output "${UNITED_OUTPUT}" \
   --email "${EMAIL}" \
   --api "${API_KEY}" \
   --number "${NUMBER_REFS}" \
   --tax_rank "${TAX_RANK}" \
   --summary_csv "${UNITED_SUMMARY}" \
   --diversity \
   --log_file "${LOGS}/UNITEd.log" \
   --traverse family
check_status "UNITEd.py sequence retrieval"

#===============================================================================

# Step 3: Read Mapping
log_section "STEP 3: READ MAPPING & BAITING"
log_info "Starting read mapping to reference sequences..."
log_info "Input directory: ${FASTP_OUTPUT}"
log_info "Reference directory: ${UNITED_OUTPUT}"
log_info "Output directory: ${MAPPING_OUTPUT}"
log_info "Aligner is: ${ALIGNER}"

python mapping_module.py \
   --input_dir "${FASTP_OUTPUT}" \
   --ref_dir "${UNITED_OUTPUT}" \
   --aligner "${ALIGNER}" \
   --output_dir "${MAPPING_OUTPUT}" \
   --tracking_sheet "${TRACKING_SHEET}" \
   --column_name "${COLUMN_NAME}" \
   --log_file "${LOGS}/mapping.log" \
   --threads 2
check_status "read mapping"

#===============================================================================

# Step 4: Contig Assembly
log_section "STEP 4: CONTIG ASSEMBLY"
log_info "Starting Contig Assembly..."
log_info "Merged reads directory: ${MAPPING_OUTPUT}"
log_info "Unmerged reads directory: ${FASTP_OUTPUT}"
log_info "Assembly output directory: ${ASSEMBLY_OUTPUT}"

python assembly_module.py \
   --merged_dir "${MAPPING_OUTPUT}" \
   --unmerged_dir "${FASTP_OUTPUT}" \
   --output_dir "${ASSEMBLY_OUTPUT}" \
   --summary_csv "${ASSEMBLY_CSV}" \
   --log_file "${LOGS}/assembly.log"
check_status "Contig Assembly"

#===============================================================================

# Step 5: BLAST Round 1 (GENERAL UNITE DATABASE)
log_section "STEP 5: BLAST ANALYSIS - ROUND 1 (GENERAL)"
log_info "Starting BLAST search against general UNITE database..."
log_info "Query directory: ${ASSEMBLY_OUTPUT}"
log_info "Database: ${BLAST_DB_GENERAL}"
log_info "Output directory: ${BLAST_ROUND1_OUTPUT}"

python blast_round1.py \
   --query_dir "${ASSEMBLY_OUTPUT}" \
   --database_file "${BLAST_DB_GENERAL}" \
   --output_dir "${BLAST_ROUND1_OUTPUT}" \
   --tracking_sheet "${TRACKING_SHEET}" \
   --column_name "${COLUMN_NAME}" \
   --log_file "${LOGS}/blast_round1.log"
check_status "BLAST round 1"

#===============================================================================

# Step 5.5: BLAST Round 1 Output Parsing & Taxonomic Validation
log_section "STEP 5.5: BLAST ROUND 1 PARSING & TAXONOMIC VALIDATION"
log_info "Starting BLAST round 1 output parsing and taxonomic validation..."
log_info "Input directory: ${BLAST_ROUND1_OUTPUT}"
log_info "Output directory: ${BLAST_PARSED_OUTPUT_ROUND1}"
log_info "Minimum sequence length: ${MIN_LENGTH}"
log_info "Minimum percent identity: ${MIN_PIDENT}"
log_info "Assembly directory: ${ASSEMBLY_OUTPUT}"


python blast_round1_parser.py \
   --input_dir "${BLAST_ROUND1_OUTPUT}" \
   --output_dir "${BLAST_PARSED_OUTPUT_ROUND1}" \
   --taxonomy_csv "${TRACKING_SHEET}" \
   --id_column "${COLUMN_NAME}" \
   --min_len "${MIN_LENGTH}" \
   --min_pident "${MIN_PIDENT}" \
   --evalue_cutoff "1e-5" \
   --assembly_dir "${ASSEMBLY_OUTPUT}" \
   --log_file "${LOGS}/parse_round1.log" \
   --allow_all
check_status "BLAST round 1 parsing"

#===============================================================================

# Step 6: BLAST Round 2a (UCHIME ITS2 DATABASE)
log_section "STEP 6: BLAST ANALYSIS - ROUND 2A (ITS2)"
log_info "Starting BLAST search against UCHIME reference ITS2-specific database..."
log_info "Query directory: ${ASSEMBLY_OUTPUT}"
log_info "Previous BLAST results: ${BLAST_ROUND1_OUTPUT}"
log_info "Database: ${BLAST_DB_ITS2}"
log_info "Output directory: ${BLAST_ROUND2A_OUTPUT}"

python blast_round2.py \
   --query_dir "${ASSEMBLY_OUTPUT}" \
   --blast_dir "${BLAST_ROUND1_OUTPUT}" \
   --database_file "${BLAST_DB_ITS2}" \
   --output_dir "${BLAST_ROUND2A_OUTPUT}" \
   --tracking_sheet "${TRACKING_SHEET}" \
   --column_name "${COLUMN_NAME}" \
   --log_file "${LOGS}/blast_round2a.log"
check_status "BLAST round 2a (ITS2)"

#===============================================================================

# Step 6.5: BLAST Output Parsing (ITS2)
log_section "STEP 6.5: BLAST OUTPUT PARSING & TAXONOMIC VALIDATION (ITS2)"
log_info "Starting BLAST output parsing and taxonomic validation for ITS2..."
log_info "Input directory: ${BLAST_ROUND2A_OUTPUT}"
log_info "Output directory: ${BLAST_PARSED_OUTPUT_ITS2}"
log_info "Minimum sequence length: ${MIN_LENGTH}"
log_info "Minimum percent identity: ${MIN_PIDENT}"
log_info "Assembly directory: ${ASSEMBLY_OUTPUT}"

python blast_output_parser.py \
   --input_dir "${BLAST_ROUND2A_OUTPUT}" \
   --output_dir "${BLAST_PARSED_OUTPUT_ITS2}" \
   --min_len "${MIN_LENGTH}" \
   --min_pident "${MIN_PIDENT}" \
   --taxonomy_csv "${TRACKING_SHEET}" \
   --assembly_dir "${ASSEMBLY_OUTPUT}" \
   --id_column "${COLUMN_NAME}" \
   --summary_csv "${PARSER_ITS2_SUMMARY}" \
   --log_file "${LOGS}/parse_round2a.log"
check_status "BLAST output parsing (ITS2)"

#===============================================================================

# Step 7: BLAST Round 2b (UCHIME ITS1 DATABASE)
log_section "STEP 7: BLAST ANALYSIS - ROUND 2B (ITS1)"
log_info "Starting BLAST search against UCHIME reference ITS1-specific database..."
log_info "Query directory: ${ASSEMBLY_OUTPUT}"
log_info "Previous BLAST results: ${BLAST_ROUND1_OUTPUT}"
log_info "Database: ${BLAST_DB_ITS1}"
log_info "Output directory: ${BLAST_ROUND2B_OUTPUT}"

python blast_round2.py \
   --query_dir "${ASSEMBLY_OUTPUT}" \
   --blast_dir "${BLAST_ROUND1_OUTPUT}" \
   --database_file "${BLAST_DB_ITS1}" \
   --output_dir "${BLAST_ROUND2B_OUTPUT}" \
   --tracking_sheet "${TRACKING_SHEET}" \
   --column_name "${COLUMN_NAME}" \
   --log_file "${LOGS}/blast_round2b.log"
check_status "BLAST round 2b (ITS1)"

#===============================================================================

# Step 7.5: BLAST Output Parsing (ITS1)
log_section "STEP 7.5: BLAST OUTPUT PARSING & TAXONOMIC VALIDATION (ITS1)"
log_info "Starting BLAST output parsing and taxonomic validation for ITS1..."
log_info "Input directory: ${BLAST_ROUND2B_OUTPUT}"
log_info "Output directory: ${BLAST_PARSED_OUTPUT_ITS1}"
log_info "Minimum sequence length: ${MIN_LENGTH}"
log_info "Minimum percent identity: ${MIN_PIDENT}"
log_info "Assembly directory: ${ASSEMBLY_OUTPUT}"

python blast_output_parser.py \
   --input_dir "${BLAST_ROUND2B_OUTPUT}" \
   --output_dir "${BLAST_PARSED_OUTPUT_ITS1}" \
   --min_len "${MIN_LENGTH}" \
   --min_pident "${MIN_PIDENT}" \
   --taxonomy_csv "${TRACKING_SHEET}" \
   --assembly_dir "${ASSEMBLY_OUTPUT}" \
   --id_column "${COLUMN_NAME}" \
   --summary_csv "${PARSER_ITS1_SUMMARY}" \
   --log_file "${LOGS}/parse_round2b.log"
check_status "BLAST output parsing (ITS1)"

#===============================================================================

# Step 8: ITS Primer Alignment
log_section "STEP 8: ITS PRIMER BINDING & EXTRACTION"
log_info "Starting ITS primer binding and extraction..."
log_info "Input directory: ${BLAST_PARSED_OUTPUT_ITS2}"
log_info "I.e. Using contigs confirmed to contain ITS2 from the correct taxon"
log_info "Output directory: ${ITS_EXTRACTION_OUTPUT}"

python its_primer_binding.py \
   --input "${BLAST_PARSED_OUTPUT_ITS2}" \
   --tracking_sheet "${TRACKING_SHEET}" \
   --column_name "${COLUMN_NAME}" \
   --output "${ITS_EXTRACTION_OUTPUT}" \
   --log_file "${LOGS}/primer_binding.log"
check_status "ITS primer binding"

#===============================================================================
# Step 9: Parse the parsed parsing results
log_section "STEP 9: AGGREGATING METRICS"
log_info "Starting metrics parsing, summarisation, and selection..."
log_info "Summarising across outputs from ${OUTPUT_BASE}"
log_info "Output directory: ${FINAL_RESULTS_OUTPUT}"

if [[ "$different_naming" == "YES" ]]; then
    log_info "Running with renaming option using: ${NAMING_TSV}"
    python its_a_summary_compiler.py \
        "${OUTPUT_BASE}" \
        --naming_tsv "${NAMING_TSV}"
else
    log_info "Running without renaming"
    python its_a_summary_compiler.py \
        "${OUTPUT_BASE}"
fi
check_status "ITS summary compiler"

#===============================================================================
# PIPELINE COMPLETION
#===============================================================================

log_section "PIPELINE COMPLETION"
echo "Start time: $(head -n 50 "${LOG_FILE}" | grep "Pipeline started" | cut -d']' -f1 | tr -d '[')"
echo "End time: $(date +'%Y-%m-%d %H:%M:%S')"
echo "Runtime: $SECONDS seconds"
echo "Input samples: $(wc -l < "${TRACKING_SHEET}") (including header)"
echo "Log file: ${LOG_FILE}"
echo "Output directories created:"
echo "  - Logs: ${LOGS}"
echo "  - Quality control: ${FASTP_OUTPUT}"
echo "  - UNITEd: ${UNITED_OUTPUT}"
echo "  - Mapped reads: ${MAPPING_OUTPUT}"
echo "  - Assemblies: ${ASSEMBLY_OUTPUT}"
echo "  - BLAST round 1: ${BLAST_ROUND1_OUTPUT}"
echo "  - BLAST round 1 parsed: ${BLAST_PARSED_OUTPUT_ROUND1}"
echo "  - BLAST round 2a (ITS2): ${BLAST_ROUND2A_OUTPUT}"
echo "  - ITS2 parsed: ${BLAST_PARSED_OUTPUT_ITS2}"
echo "  - BLAST round 2b (ITS1): ${BLAST_ROUND2B_OUTPUT}"
echo "  - ITS1 parsed: ${BLAST_PARSED_OUTPUT_ITS1}"
echo "  - ITS primer extraction: ${ITS_EXTRACTION_OUTPUT}"
echo "  - Final results: ${FINAL_RESULTS_OUTPUT}"
echo "==============================================================================="

log_info "Pipeline execution complete. Check ${LOG_FILE} for detailed logs."
