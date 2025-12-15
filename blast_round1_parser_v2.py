import os
import sys
import gzip
import argparse
import re
import logging
import pandas as pd
from pathlib import Path
from collections import defaultdict
from datetime import datetime
from metahist_tools import setup_logging

# BLAST output parser v0.0.0
# Processes BLAST TSV files (outfmt 6) to extract best candidate contig based on blast results.
# Contigs are filtered based on minimum length and percent identity, as well as e-value cutoff.
# Outputs a new TSV file with filtered results.
# This includes an option to extract the top N hits per query sequence
# Also the script generates a taxonomy validation summary CSV file
# with contig analysis results based on a user-provided taxonomy CSV file.
# A FASTA file with the contigs that passed the taxonomy validation is also created.
#
# Example usage:
#     python blast_round1_parser.py --input_dir /path/to/blast_files --taxonomy_csv taxonomy.csv \
#       --assembly_dir /path/to/directory/of/assemblies -o output_dir --min_len 100 --allow_all
#
# Taxonomy CSV format: The CSV file must contain 'ID' and taxonomic hierarchy columns.
# If taxonomic heirarchy columns are missing, you can use "pull_ncbi_lineage" to generate an
# updated spreadsheet for this.
#
# Authors: Maria Kamouyiaros & Dan Parsons @ NHMUK
# Date: 2025-11-28
# License: MIT
# Version: 2.0.0

logger = logging.getLogger(__name__)

def detect_blast_format(input_file):
    """
    Detects BLAST outfmt6 format and whether it has a valid header line.
    Returns: (has_header: bool, columns: list[str])
    """
    standard_blast_columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]

    try:
        with open(input_file, 'r', encoding='utf-8') as f:
            # Skip leading comment lines that start with '#'
            first_line = None
            for line in f:
                line = line.strip()
                if not line:
                    continue
                if line.startswith('#'):
                    continue
                first_line = line
                break

            if not first_line:
                # empty file or only comments
                return False, standard_blast_columns

            columns = first_line.split('\t')

            blast_col_lower = [col.lower() for col in standard_blast_columns]
            first_line_lower = [col.lower() for col in columns]

            if len(columns) >= 12 and all(
                col in blast_col_lower for col in first_line_lower[:12]
            ):
                logger.info("  -> Detected header in file: %s", first_line)
                return True, columns[:12]

            return False, standard_blast_columns

    except Exception as e:
        logger.warning("Warning: Could not detect format for %s : %s", input_file, e)
        return False, standard_blast_columns

def extract_full_taxonomy_from_sseqid(sseqid):
    ''' Extracts full taxonomic lineage from sseqid string if "k__" pattern is present.'''

    # Look for the taxonomic part (contains k__ pattern)
    parts = sseqid.split('|')
    for part in parts:
        if 'k__' in part and ';' in part:
            return part.strip()
    return None

def build_taxonomy_string(row, taxonomy_columns):
    ''' Function to build a full taxonomy lineage string from a row of a dataframe. '''

    taxonomy_parts = []
    # Map column names (case-insensitive) to their prefixes
    prefix_map = {
        'kingdom': 'k__',
        'phylum': 'p__',
        'class': 'c__',
        'order': 'o__',
        'family': 'f__',
        'subfamily': 'sf__',
        'tribe': 't__',
        'genus': 'g__',
        'species': 's__',
        'subspecies': 'ss__'
    }

    for col in taxonomy_columns:
        if col in row.index:
            value = row[col]
            # Check if value exists and is not empty/NaN
            if pd.notna(value) and str(value).strip() and str(value).strip().lower() != 'nan':
                # Get prefix based on column name (case-insensitive)
                col_lower = col.lower()
                prefix = prefix_map.get(col_lower, '')
                taxonomy_parts.append(f"{prefix}{str(value).strip()}")

    return ';'.join(taxonomy_parts) if taxonomy_parts else ''

def load_taxonomy_mapping(taxonomy_csv, id_col):
    """
    Loads taxonomy mapping from a CSV file into a dictionary.
    Mapping keys are IDs, values are dicts with all taxonomic ranks + full_taxonomy string.
    """

    try:
        df = pd.read_csv(taxonomy_csv)

        # Define expected taxonomy columns in order
        taxonomy_columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family',
                            'Subfamily', 'Tribe', 'Genus', 'Species', 'Subspecies']

        # Make column matching case-insensitive
        column_mapping = {}
        available_columns = list(df.columns)

        # Find ID column
        for col in df.columns:
            if col.lower() == id_col.lower():
                column_mapping[id_col] = col
                break

        if id_col not in column_mapping:
            logger.error("CSV file must contain %s column (case-insensitive). Found columns: %s",
                         id_col, available_columns)
            raise ValueError(f"ID column '{id_col}' not found in CSV")

        # Find taxonomy columns (case-insensitive)
        for tax_col in taxonomy_columns:
            for col in df.columns:
                if col.lower() == tax_col.lower():
                    column_mapping[tax_col] = col
                    break

        found_tax_columns = [tax_col for tax_col in taxonomy_columns if tax_col in column_mapping]
        if not found_tax_columns:
            logger.error("CSV file must contain at least one taxonomic hierarchy column: %s. Found columns: %s",
                         taxonomy_columns, available_columns)
            raise ValueError("No taxonomy columns found in CSV")

        logger.info("Found taxonomy columns: %s", found_tax_columns)

        # Create mapping dictionary
        id_col_actual = column_mapping[id_col]  # Column name as it appears in CSV
        mapping = {}

        for _, row in df.iterrows():
            sample_id = row[id_col_actual]
            tax_data = {}

            # Store all found taxonomic ranks
            for tax_col in found_tax_columns:
                value = row[column_mapping[tax_col]]
                tax_data[tax_col.lower()] = str(value).lower().strip() if pd.notna(value) else None

            # Build full taxonomy string
            mapped_columns = [column_mapping[tax_col] for tax_col in found_tax_columns]
            tax_data['full_taxonomy'] = build_taxonomy_string(row, mapped_columns)

            mapping[sample_id] = tax_data

        logger.info("Loaded taxonomy mapping for %d IDs from %s", len(mapping), taxonomy_csv)

        # Debug: show first 5 entries
        logger.debug("DEBUG: Taxonomy mapping (first 5 entries):")
        for i, (sample_id, tax_data) in enumerate(mapping.items()):
            if i >= 5:
                break
            logger.info("  %s -> %s", sample_id, tax_data)

        return mapping

    except Exception:
        logger.exception("Error loading taxonomy CSV %s", taxonomy_csv)
        raise

def _grep(x, patterns_list):
    '''takes a string (x) and checks whether any string from a given 
    list of strings (patterns_list) exists in `x`'''

    x_norm = x.lower().replace("_", " ")
    for text in patterns_list:
        text_norm = text.lower().strip()
        if re.search(rf'\b{text_norm}\b', x_norm):
            return True
    return False

def filter_df(df, patterns_list):
    ''' Filters a blast result dataframe to keep only rows where sseqid 
    contains any of the patterns. Assumes headers'''

    mask = df['sseqid'].apply(_grep, patterns_list=patterns_list)
    df2 = df[mask]
    if df2.empty:
        logger.info("    -> No matches found after filtering dataframe.")
        return pd.DataFrame()  # Return empty DataFrame

    return df2

def get_expected_taxonomy_from_filename(filename, taxonomy_mapping, level):
    """
    Grab the expected taxonomy from a file name using the taxonomy mapping dictionary.
    
    Arguments:
        filename: str or Path ; the file name to match
        taxonomy_mapping: dict ; mapping of sample_id -> taxonomic info
        level: str ; taxonomic rank to return ('family', 'genus', 'species', etc.)
    
    Returns:
        (taxon, full_taxonomy, matched_id)
        taxon: str or None ; the requested taxonomic rank for this sample
        full_taxonomy: str or None ; full taxonomy string
        matched_id: str or None ; the matched sample ID
    """

    filename_base = Path(filename).name

    # Ensure level is lowercase to match mapping keys
    level = level.lower()

    for sample_id, tax_data in taxonomy_mapping.items():
        # Exact match of sample ID in filename
        if str(sample_id) == filename_base or filename_base.startswith(f"{sample_id}_"):
            taxon = tax_data.get(level)  # May be None if not available
            full_taxonomy = tax_data.get('full_taxonomy')
            return taxon, full_taxonomy, sample_id

    return None, None, None


def process_blast_results(input_dir, output_dir, taxonomy_mapping,
    min_len=100, min_pident=80, evalue_cutoff=1e-5, allow_all=False):
    '''Creates a taxonomy validation summary CSV file based on parsed BLAST output files.
    Functions by checking contig hits against expected taxonomy with the following logic:
    pull blast tsv
    search species
        if yes pull all species hits, filter to unique contigs, 
        if --allow_all: output FASTA of contigs
        else: filter contigs down by parameters set (min_len, perc_ident)
    if no species, search genus
        if yes pull all genus hits, filter to unique contigs, output FASTA of contigs
        if --allow_all: output FASTA of contigs
        else: filter contigs down by parameters set (min_len, perc_ident)
    if no genus, search family
        if yes pull all family hits, filter to unique contigs, output FASTA of contigs
       if --allow_all: output FASTA of contigs
       else: filter contigs down by parameters set (min_len, perc_ident)
   if no family report none found
    '''

    logger.info("\nCreating taxonomy validation summary...")
    logger.info("Searching for files in %s", input_dir)

    # Find all parsed blast output files (*.tsv) in input_dir
    parsed_files = list(Path(input_dir).glob('*.tsv'))

    if not parsed_files:
        logger.warning("No parsed BLAST output files found. \
            Skipping taxonomy validation summary creation.")
        return

    summary_results = []
    processed_count = 0

    # For each parsed file, extract expected taxonomy and compare with BLAST hits
    for parsed_file in sorted(parsed_files):
        logger.info("  Processing %s...", parsed_file.name)

        # Extract original filename to get expected taxonomy
        original_name = re.sub(r'-(top_\d+_hits|all_hits|blast_results|hits)$', '',
            parsed_file.stem) + '.tsv'
        expected_family, expected_full_taxonomy, matched_id = get_expected_taxonomy_from_filename(
            original_name, taxonomy_mapping, level="family")
        
        # Get all taxonomy levels for this sample
        expected_family = None
        expected_genus = None
        expected_species = None
        if matched_id:
            expected_family = taxonomy_mapping[matched_id].get('family')
            expected_genus = taxonomy_mapping[matched_id].get('genus')
            expected_species = taxonomy_mapping[matched_id].get('species')

        output_file_path = Path(output_dir) / f"filtered_{parsed_file.stem}.tsv"

        # Initialise result dictionary with blast_round1_ prefix column names
        result = {
            'ID': matched_id if matched_id else 'NA',
            'blast_round1_found_taxon': 'NA',
            'blast_round1_expected_taxonomy': expected_full_taxonomy if expected_full_taxonomy else 'NA',
            'blast_round1_n_contigs_hits': 'NA',
            'blast_round1_n_contigs_in': 'NA',
            'blast_round1_contigs': 'NA',
            'blast_round1_hit_taxonomy': 'NA',
            'blast_round1_hit_coords': 'NA',
            'blast_round1_correct_taxonomy': 'NA'
        }

        # Read the parsed BLAST output tsv file
        try:
            headers = ['qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
                       'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore']
            logger.info("    Processing %s...", parsed_file)
            has_header, detected_cols = detect_blast_format(parsed_file)
            if has_header:
                logger.info("  Detected BLAST TSV with header matching: %s", detected_cols)
            else:
                logger.info("  Detected BLAST TSV without header matching: %s", detected_cols)

            # Read the BLAST results
            if has_header:
                df = pd.read_csv(parsed_file, sep='\t', header=0)
                logger.info("  -> Reading file with headers")
            else:
                df = pd.read_csv(parsed_file, sep='\t', header=None, names=headers)
                logger.info("  -> Reading file without headers")

            df.columns = df.columns.astype(str).str.strip().str.replace("\ufeff", "")

            # number of unique contigs parsed == number of unique qseqid in TSV file
            result['blast_round1_n_contigs_in'] = df['qseqid'].nunique()
            if df.empty:
                logger.info("    -> No hits found in file -> FAIL")
                result['blast_round1_correct_taxonomy'] = 'FAIL'
                summary_results.append(result)
                processed_count += 1
                continue

        except Exception:
            logger.exception("    -> Error reading %s", parsed_file)
            summary_results.append(result)
            continue

        # Search taxa in BLAST results with cascading fallback:
        # Try species first, then genus, then family
        df2 = pd.DataFrame()
        unique_contigs = 0

        # Try species-level match
        if expected_species is not None:
            df2 = filter_df(df, [expected_species])
            if not df2.empty:
                logger.info("    -> Found species matches for expected species %s",
                    expected_species)
                unique_contigs = df2['qseqid'].nunique()
                logger.info("    -> Unique contigs matching expected species: %d",
                    unique_contigs)
                result['blast_round1_found_taxon'] = expected_species

        # If no species match, try genus-level match
        if df2.empty and expected_genus is not None:
            logger.info("    -> No species matches found for expected species %s",
                expected_species)
            df2 = filter_df(df, [expected_genus])
            if not df2.empty:
                logger.info("    -> Found genus matches for expected genus %s",
                    expected_genus)
                unique_contigs = df2['qseqid'].nunique()
                logger.info("    -> Unique contigs matching expected genus: %d",
                    unique_contigs)
                result['blast_round1_found_taxon'] = expected_genus

        # If no genus match, try family-level match
        if df2.empty and expected_family is not None:
            logger.info("    -> No genus matches found for expected genus %s",
                expected_genus)
            df2 = filter_df(df, [expected_family])
            if not df2.empty:
                logger.info("    -> Found family matches for expected family %s",
                    expected_family)
                unique_contigs = df2['qseqid'].nunique()
                logger.info("    -> Unique contigs matching expected family: %d",
                    unique_contigs)
                result['blast_round1_found_taxon'] = expected_family

        # If no matches at any level
        if df2.empty:
            logger.info("    -> No species matches found for expected species %s",
                expected_species)
            logger.info("    -> No genus matches found for expected genus %s",
                expected_genus)
            logger.info("    -> No family matches found for expected family %s",
                expected_family)
            result['blast_round1_correct_taxonomy'] = 'FAIL'
            summary_results.append(result)
            processed_count += 1
            continue

        # Process matched contigs
        if unique_contigs >= 1 and not allow_all:
            logger.info("    -> %d contigs match expected taxonomy -> applying filters",
                unique_contigs)
            df_filtered = filter_blast(min_len=min_len, min_pident=min_pident,
                input_df=df2, evalue_cutoff=evalue_cutoff)
            unique_contigs_new = df_filtered['qseqid'].nunique()
            if unique_contigs_new == 0:
                logger.info("    -> No contigs remain after filtering -> FAIL")
                result['blast_round1_correct_taxonomy'] = 'FAIL'
            elif unique_contigs_new > 1:
                logger.info("    -> Multiple contigs (%d) remain after filtering -> FAIL",
                    unique_contigs_new)
                result['blast_round1_n_contigs_hits'] = unique_contigs_new
                result['blast_round1_contigs'] = ';'.join(sorted(df_filtered['qseqid'].unique()))
                result['blast_round1_correct_taxonomy'] = 'FAIL'
            else:
                logger.info("    -> Exactly one contig remains after filtering -> PASS")
                qseqid = df_filtered['qseqid'].iloc[0]
                hit_row = df_filtered[df_filtered['qseqid'] == qseqid].iloc[0]
                result['blast_round1_n_contigs_hits'] = 1
                result['blast_round1_contigs'] = qseqid
                result['blast_round1_correct_taxonomy'] = 'PASS'
                result['blast_round1_hit_taxonomy'] = extract_full_taxonomy_from_sseqid(hit_row['sseqid'])
                result['blast_round1_hit_coords'] = f"{hit_row['qstart']}-{hit_row['qend']}"
                hit_row.to_frame().T.to_csv(output_file_path, sep="\t", index=False)
        elif unique_contigs >= 1 and allow_all:
            logger.info("    -> %d contigs match expected taxonomy -> allowing all",
                unique_contigs)
            result['blast_round1_n_contigs_hits'] = unique_contigs
            result['blast_round1_contigs'] = ';'.join(sorted(df2['qseqid'].unique()))
            result['blast_round1_correct_taxonomy'] = 'PASS'
            # Get hit taxonomy from first matching row
            first_hit = df2.iloc[0]
            result['blast_round1_hit_taxonomy'] = extract_full_taxonomy_from_sseqid(first_hit['sseqid'])
            result['blast_round1_hit_coords'] = f"{first_hit['qstart']}-{first_hit['qend']}"
            df2.to_csv(output_file_path, sep="\t", index=False)

        summary_results.append(result)
        processed_count += 1

    # Create summary CSV
    summary_file = Path(output_dir) / 'blast_validation_summary.csv'
    summary_df = pd.DataFrame(summary_results)

    # Define column order with blast_round1_ prefix names
    column_order = ['ID', 'blast_round1_found_taxon', 'blast_round1_expected_taxonomy',
        'blast_round1_n_contigs_in', 'blast_round1_n_contigs_hits', 'blast_round1_contigs', 
        'blast_round1_hit_taxonomy', 'blast_round1_hit_coords',
        'blast_round1_correct_taxonomy']
    summary_df = summary_df[column_order]

    # Save to CSV
    summary_df.to_csv(summary_file, index=False)

    logger.info("  Processed %d files", processed_count)
    logger.info("  Saved taxonomy validation summary to: %s", summary_file)

    # Summary statistics
    total_pass = len(summary_df[summary_df['blast_round1_correct_taxonomy'] == 'PASS'])
    total_fail = len(summary_df[summary_df['blast_round1_correct_taxonomy'] == 'FAIL'])
    total_with_taxonomy = len(summary_df[summary_df['blast_round1_found_taxon'] != 'NA'])

    logger.info("  Summary statistics:")
    logger.info("    -> Total samples: %d", len(summary_df))
    logger.info("    -> Samples with taxonomy mapping: %d", total_with_taxonomy)
    logger.info("    -> PASS (exactly 1 contig hit): %d", total_pass)
    logger.info("    -> FAIL (0 or >1 contig hits): %d", total_fail)

    # Return the summary dataframe for further processing if needed
    return summary_df

def filter_blast(input_df, top_n=None, min_len=None, min_pident=None, df=None, evalue_cutoff=1e-5):
    '''Filter a single BLAST TSV file based on the following order: length, pident, evalue,
    then (optional) extract top N hits per query sequence.'''

    logger.info("Filtering BLAST results...")
    df = input_df.copy()
    nrows_before = len(df)

    logger.info("  -> Initial number of hits: %d", nrows_before)

    # Apply minimum length filter if specified
    if min_len is not None:
        df_filtered = df[df['length'] >= min_len]
        nrows_after = len(df_filtered)
        filtered_out = nrows_before - nrows_after
        logger.info("  -> Applied minimum length filter %s: %s hits filtered out %d remaining", min_len, filtered_out, nrows_after)
    else:
        df_filtered = df

    # Apply minimum percent identity filter if specified
    if min_pident is not None:
        df_filtered2 = df_filtered[df_filtered['pident'] >= min_pident]
        nrows_after = len(df_filtered2)
        nrows_before = len(df_filtered)
        filtered_out = nrows_before - nrows_after
        logger.info("  -> Applied minimum pident filter %s: %s hits filtered out %d remaining", min_pident, filtered_out, nrows_after)
        df = df_filtered2
    else:
        df_filtered2 = df_filtered

    # Apply minimum cutoff evalue (Default evalue cutoff is 1e-5)
    if evalue_cutoff is not None:
        df_filtered3 = df_filtered2[df_filtered2['evalue'] <= evalue_cutoff]
        nrows_after = len(df_filtered3)
        nrows_before = len(df_filtered2)
        filtered_out = nrows_before - nrows_after
        logger.info("  -> Applied evalue cutoff filter %s: %s hits filtered out %d remaining", evalue_cutoff, filtered_out, nrows_after)
    else:
        df_filtered3 = df_filtered2

    # Sort by qseqid first, then by our ranking criteria
    # pident (desc), length (desc), evalue (asc)
    df_sorted = df_filtered3.sort_values([
        'qseqid',           # Group by query sequence
        'length',           # Longer alignments are better
        'evalue'            # Lower e-value is better
    ], ascending=[True, False, True])

    # Group by qseqid and take top N hits for each (or all hits if top_n is None)
    if top_n is not None:
        top_hits = df_sorted.groupby('qseqid').head(top_n)
    else:
        top_hits = df_sorted

    logger.info("  -> Final number of hits after filtering: %d", len(top_hits))
    return top_hits

def write_blast_out(input_df, output_dir, input_file, top_n=None):
    '''Writes blast filtered data frame to a new TSV file.'''

    # Select only the columns we want in the output
    output_columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'qstart', 'qend']
    result = input_df[output_columns]

    try:
        if input_df.empty:
            logger.warning("No hits to write for %s. \
                Skipping file creation.", input_file)
            return 0, 0
        input_name = Path(input_file).stem
        if top_n is not None:
            output_file = Path(output_dir) / f"{input_name}-top_{top_n}_hits.tsv"
        else:
            output_file = Path(output_dir) / f"{input_name}-all_hits.tsv"

            # Save the results with column headers
        result.to_csv(output_file, sep='\t', index=False)
        logger.info("  -> Saved %d hits from %d queries to %s", \
            len(result), input_df['qseqid'].nunique(), output_file)

    except Exception:
        logger.exception("Error processing %s", input_file)
        return 0, 0

def open_fasta(filename):
    '''Opens a FASTA file, handling gzip if necessary.'''

    filename = str(filename)  # Ensure filename is a string
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r', encoding='utf-8', errors='ignore')

def grep_fasta_header(fasta_file, patterns, output_file):
    '''Search FASTA headers for patterns and write matching records.
    
    Args:
        fasta_file (str): Path to input FASTA file
        patterns (list): List of patterns to search for in headers
        output_file (str): Path to output file
    '''

    logger.info("Searching for %s in %s...", patterns, fasta_file)
    with open_fasta(fasta_file) as fa, open(output_file, 'w', encoding='utf-8') as out:
        write_seq = False
        buffer = []
        for line in fa:
            if line.startswith(">"):
                # flush old buffer if it was a match
                if write_seq and buffer:
                    out.writelines(buffer)
                buffer = [line]
                header = line.strip()
                write_seq = any(pat in header for pat in patterns)
            else:
                buffer.append(line)
        # flush last record
        if write_seq and buffer:
            logger.info("Writing = %s for %s", write_seq, patterns)
            out.writelines(buffer)

def extract_passing_contigs(df, output_dir, assembly_dir):
    '''Extracts FASTA sequences for contigs that passed taxonomy validation.
    This works for real lists in the format of ["NODE_5", "NODE_9"] etc.'''

    logger.info("\nExtracting FASTA sequences for passed contigs...")

    if not assembly_dir or not Path(assembly_dir).is_dir():
        logger.warning("Assembly directory '{assembly_dir}' does not exist \
            or is not a directory. Skipping FASTA extraction.")
        return

    dfkeep = df.loc[df["blast_round1_correct_taxonomy"] == "PASS"]

    extracted_count = 0
    for _, row in dfkeep.iterrows():
        filename = row["ID"]
        if not filename or pd.isna(filename):
            continue

        # wrap in list so grep function works
        contig_patterns = row["blast_round1_contigs"].split(";")
        output_file = Path(output_dir) / f"{filename}_parsed_contig.fasta"

        assembly_file = Path(assembly_dir) / f"{filename}.spades.out" / "scaffolds.fasta"
        if not assembly_file.is_file():
            assembly_file = Path(assembly_dir) / f"{filename}.spades.out" / "contigs.fasta"
            if not assembly_file.is_file():
                logger.warning("  -> Assembly file not found for %s", filename)
                continue

        grep_fasta_header(assembly_file, contig_patterns, output_file)
        if output_file.exists():
            extracted_count += 1
            logger.info("  -> Extracted contig for %s", filename)

    logger.info("  Extracted FASTA sequences for %d passed samples", extracted_count)
    return extracted_count

def update_summary_with_contig_paths(output_dir):
    '''Updates the taxonomy validation summary CSV with absolute contig file paths.'''

    # Convert output_dir to absolute path first
    output_dir = Path(output_dir).resolve()

    summary_file = output_dir / 'blast_validation_summary.csv'
    if not summary_file.exists():
        return

    df = pd.read_csv(summary_file)

    updated_count = 0
    for idx, row in df.iterrows():
        if row['blast_round1_correct_taxonomy'] == 'PASS' and row['ID']:
            contig_fasta = output_dir / f"{row['ID']}_parsed_contig.fasta"
            if contig_fasta.exists():
                df.at[idx, 'contig_path'] = str(contig_fasta)
                updated_count += 1

    # Save updated CSV
    df.to_csv(summary_file, index=False)
    logger.info("  Updated %d contig paths in summary CSV", updated_count)

def main(args):
    '''Main function to parse arguments and process BLAST TSV files.'''

    # Initialise set up:
    # - Create output directory if it doesn't exist
    # - Setup log file and log directory
    os.makedirs(args.output_dir, exist_ok=True)

    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'blast_round1_processing_{timestamp}.log'
    log_dir = os.path.join(args.output_dir, 'logs')
    setup_logging(log_dir=log_dir, log_file=args.log_file)

    # Load taxonomy columns from input CSV to be used in mapping to blast output files
    # output of load_taxonomy_mapping is a dictionary of {ID: {'family': family,
    # 'full_taxonomy': full_taxonomy}}
    try:
        taxonomy_mapping = load_taxonomy_mapping(args.taxonomy_csv, args.id_column)
    except Exception:
        logger.exception("Failed to load taxonomy CSV. Exiting.")
        return

    # Find all BLAST TSV files in input directory
    tsv_files = list(Path(args.input_dir).glob('*.tsv'))

    if not tsv_files:
        logger.error("No TSV files found in %s or its subdirectories", args.input_dir)
        return

    logger.info("BLAST Results parsing started with the following parameters:")
    logger.info("    Input directory: %s", args.input_dir)
    logger.info("    Output directory: %s", args.output_dir)
    logger.info("    Taxonomy CSV: %s", args.taxonomy_csv)
    logger.info("    Log file: %s", args.log_file)
    logger.info("    Found %d TSV files to process", len(tsv_files))
    if args.top_n is not None:
        logger.info("    Keeping top %d hits per query sequence", args.top_n)
    else:
        logger.info("    Keeping all hits per query sequence")
    if args.min_len is not None:
        logger.info("    Minimum alignment length filter: %s", args.min_len)
    if args.min_pident is not None:
        logger.info("    Minimum percentage identity filter: %s", args.min_pident)
    if args.evalue_cutoff is not None:
        logger.info("    E-value cutoff filter: %s", args.evalue_cutoff)
    if args.assembly_dir:
        logger.info("    Assembly directory for FASTA extraction: %s", args.assembly_dir)
    logger.info("-" * 60)
    logger.info("Auto-detecting files with/without headers")
    logger.info("Taxonomy validation summary will be created for all processed files")

    # Process each BLAST TSV file using process_blast_results function
    # the output of this is a csv file with a full summary of contigs based on
    # input parameters and a dataframe of the file results
    processed_blast_df = process_blast_results(args.input_dir, args.output_dir, taxonomy_mapping,
        min_len = args.min_len, min_pident=args.min_pident, evalue_cutoff=args.evalue_cutoff,
        allow_all=args.allow_all)

   # Exit if no results were processed
    if processed_blast_df is None or processed_blast_df.empty:
        logger.error("No valid BLAST results processed. Exiting.")
        sys.exit(1)

    logger.info("\nBLAST Results parsing completed.")
    logger.info("Total number of processed files: %d", len(tsv_files))
    logger.info("=" * 60)

    # Assembly directory argument means FASTAs of passed contigs should be created
    if args.assembly_dir:
        logger.info("\n" + "=" * 60)
        logger.info("ASSEMBLY DIRECTORY PROVIDED, CREATING FASTAS OF ALL 'PASS' CONTIGS:")
        logger.info("=" * 60)

        # Use the existing extract_passing_contigs function
        extracted_count = extract_passing_contigs(processed_blast_df, args.output_dir,
            args.assembly_dir)

        if extracted_count:
            # Update the summary file with contig paths
            update_summary_with_contig_paths(args.output_dir)

if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser = argparse.ArgumentParser()
        parser.print_help()
        sys.exit(1)

    parser = argparse.ArgumentParser(description='Process BLAST TSV files \
        to extract top hits per query')
    parser.add_argument('-i', '--input_dir', required=True, help='Directory containing BLAST TSV \
        files (required)')
    parser.add_argument('--taxonomy_csv', required=True,
        help='CSV file with ID and taxonomic hierarchy columns (required)')
    parser.add_argument('-o', '--output_dir', default='parsed_blast_output',
        help='Output directory (default: parsed_blast_output)')
    parser.add_argument('-n', '--top_n', type=int, default=None,
        help='Number of top hits to keep per query (default: all hits)')
    parser.add_argument('--min_len', type=int, default=None,
        help='Minimum alignment length to consider (optional)')
    parser.add_argument('--min_pident', type=int, default=None,
        help='Minimum percentage identity to consider (optional)')
    parser.add_argument('--evalue_cutoff', type=float, default=1e-5,
        help='E-value cutoff to consider (default: 1e-5)')
    parser.add_argument('--log_file', default=None,
        help='Log file path (default: blast_processor_YYYYMMDD_HHMMSS.log)')
    parser.add_argument('--assembly_dir', default=None, help='Directory containing \
        assemblies for FASTA extraction (default: current directory)')
    parser.add_argument('--id_column', default = "ID", help='Column name in taxonomy CSV that contains \
        the IDs matching filenames')
    parser.add_argument('--allow_all', action='store_true',
        help='If set, allows multiple contig hits to PASS without filtering')

    args = parser.parse_args()
    main(args)
