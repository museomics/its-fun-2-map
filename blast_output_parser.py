#!/usr/bin/env python3
"""
BLAST output parser

Processes BLAST TSV files (outfmt 6) to extract top hits per query sequence.
For each query sequence, sorts hits by pident (desc), length (desc), evalue (asc)
and keeps the top N hits (default: all). Also generates a taxonomy validation summary CSV file 
with contig analysis results based on a user-provided taxonomy CSV file.
A FASTA file with the contigs that passed the taxonomy validation is also created.

Multi-tier taxonomy matching decision process:
    If expected_family is empty:
        - Fall back to genus-level matching
        - 0 genus matches: FAIL
        - 1 genus match: PASS
        - >1 genus matches: Check species level
            - If no expected species data: FAIL
            - If exactly 1 contig matches genus AND species: PASS
            - If 0 or >1 contigs match genus AND species: FAIL
    
    If expected_family is present:
        - 0 family matches: FAIL
        - 1 family match: PASS   
        - >1 family matches: Check genus level
            - If no expected genus data: FAIL
            - If exactly 1 contig matches both family AND genus: PASS
            - If 0 or >1 contigs match both family AND genus: Check species level
                - If no expected species data: FAIL
                - If exactly 1 contig matches family AND genus AND species: PASS
                - If 0 or >1 contigs match family AND genus AND species: FAIL

Example usage:
    python blast_output_parser.py --input_dir /path/to/blast_files --taxonomy_csv taxonomy.csv --assembly_dir /path/to/directory/of/assemblies -o output_dir -n 50 --min_len 100

Taxonomy CSV format: The CSV file must contain 'ID' and taxonomic hierarchy columns:
ID,Accession,taxid,forward,reverse,Kingdom,Phylum,Class,Order,Family,Subfamily,Tribe,Genus,Species,Subspecies

Author: Dan Parsons & Maria Kamouyiaros @ NHMUK
Date: 2024-06-20
License: MIT
"""

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
from its_fun_tools import log_and_print    


def count_contigs_in_assembly(assembly_dir, sample_id):
    """
    Count total contigs in the scaffolds.fasta file for a sample.
    
    Args:
        assembly_dir: Path to directory containing assembly outputs
        sample_id: Sample identifier
    
    Returns:
        int: Number of contigs (sequences) in the scaffolds.fasta file, or 0 if not found
    """
    if not assembly_dir:
        return 0
    
    scaffolds_file = Path(assembly_dir) / f"{sample_id}.spades.out" / "scaffolds.fasta"
    
    if not scaffolds_file.exists():
        log_and_print(f"  -> No scaffolds.fasta found for {sample_id}", level='debug')
        return 0
    
    count = 0
    try:
        with open(scaffolds_file, 'r') as f:
            for line in f:
                if line.startswith('>'):
                    count += 1
    except Exception as e:
        log_and_print(f"  -> Error counting contigs in {scaffolds_file}: {e}", level='warning')
        return 0
    
    return count

        
def detect_blast_format(input_file):
    ''' Detects blast outfmt 6 format and whether it has a valid header line.'''
    # Standard BLAST outfmt 6 column names
    standard_blast_columns = [
        'qseqid', 'sseqid', 'pident', 'length', 'mismatch', 'gapopen',
        'qstart', 'qend', 'sstart', 'send', 'evalue', 'bitscore'
    ]
    
    try:
        with open(input_file, 'r') as f:
            first_line = f.readline().strip()
            
            # Check if first line looks like a header
            if first_line.startswith('#'):
                # Comment line, skip and check next
                first_line = f.readline().strip()
            
            # Split the first line and check if it matches expected column names
            columns = first_line.split('\t')
            
            # Check if first line contains column names (case-insensitive)
            blast_col_lower = [col.lower() for col in standard_blast_columns]
            first_line_lower = [col.lower() for col in columns]
            
            # If the first line matches BLAST column names, it's a header
            if (len(columns) >= 12 and 
                all(col in blast_col_lower for col in first_line_lower[:12])):
                log_and_print(f"  -> Detected header in file: {first_line}")
                return True, columns[:12]  # Use the actual column names from file
            else:
                # No header detected, use standard names
                return False, standard_blast_columns
                
    except Exception as e:
        log_and_print(f"Warning: Could not detect format for {input_file}: {e}", level='warning')
        return False, standard_blast_columns


def extract_species_from_taxonomy(sseqid):
    ''' Extracts species from sseqid string if ";s__" pattern is present.'''

    if ';s__' in sseqid:
        # Find the species part
        species_part = sseqid.split(';s__')[1].split(';')[0]
        species_part = re.sub("_", " ", species_part)
        # Convert to lowercase and strip whitespace for standardisation
        return species_part.lower().strip()
    return None
    
    
def extract_genus_from_taxonomy(sseqid):
    ''' Extracts genus from sseqid string if ";g__" pattern is present.'''

    if ';g__' in sseqid:
        # Find the genus part
        genus_part = sseqid.split(';g__')[1].split(';')[0]
        # Convert to lowercase and strip whitespace for standardisation
        return genus_part.lower().strip()
    return None
  
  
def extract_family_from_taxonomy(sseqid):
    ''' Extracts family from sseqid string if ";f__" pattern is present.'''

    if ';f__' in sseqid:
        # Find the family part
        family_part = sseqid.split(';f__')[1].split(';')[0]
        # Convert to lowercase and strip whitespace for standardisation
        return family_part.lower().strip()
    return None

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
    '''Loads taxonomy mapping from a CSV file into a dictionary. 
    Mapping keys are IDs, values are dicts with 'family' and 'full_taxonomy'.'''

    try:
        df = pd.read_csv(taxonomy_csv)
        
        # Define expected taxonomy columns in order
        taxonomy_columns = ['Kingdom', 'Phylum', 'Class', 'Order', 'Family', 'Subfamily', 'Tribe', 'Genus', 'Species', 'Subspecies']
        
        # Make column matching case-insensitive
        column_mapping = {}
        available_columns = list(df.columns)
        
        # Find ID column
        for col in df.columns:
            if col.lower() == id_col.lower():
                column_mapping[id_col] = col
                break
        
        # Find taxonomy columns (case-insensitive)
        for tax_col in taxonomy_columns:
            for col in df.columns:
                if col.lower() == tax_col.lower():
                    column_mapping[tax_col] = col
                    break
        
        # Check if we found ID column
        if id_col not in column_mapping:
            raise ValueError(f"CSV file must contain {id_col} column (case-insensitive). Found columns: {available_columns}")
        
        # Check which taxonomy columns we found
        found_tax_columns = [tax_col for tax_col in taxonomy_columns if tax_col in column_mapping]
        if not found_tax_columns:
            raise ValueError(f"CSV file must contain at least one taxonomic hierarchy column: {taxonomy_columns}. Found columns: {available_columns}")
        
        log_and_print(f"Found taxonomy columns: {found_tax_columns}")
        
        # Create mapping dictionary
        id_col = column_mapping[id_col] # Column ID from column_mapping
        mapping = {}
        
        for _, row in df.iterrows():
            sample_id = row[id_col]
            
            # Get family (for backward compatibility)
            family = ''
            if 'Family' in column_mapping:
                family_value = row[column_mapping['Family']]
                if pd.notna(family_value):
                    family = str(family_value).lower().strip()
            
            # Build full taxonomy string
            mapped_columns = [column_mapping[tax_col] for tax_col in found_tax_columns if tax_col in column_mapping]
            full_taxonomy = build_taxonomy_string(row, mapped_columns)
            
            mapping[sample_id] = {
                'family': family,
                'full_taxonomy': full_taxonomy
            }
        
        log_and_print(f"Loaded taxonomy mapping for {len(mapping)} IDs from {taxonomy_csv}")
        
        # Debug: show what we loaded
        #log_and_print("DEBUG: Taxonomy mapping (first 5 entries):")
        #for i, (sample_id, tax_data) in enumerate(mapping.items()):
        #    if i >= 5:  # Only show first 5 entries
        #        break
        #    log_and_print(f"  {sample_id} -> family: {tax_data['family']}, full_taxonomy: {tax_data['full_taxonomy']}")
        
        return mapping
        
    except Exception as e:
        log_and_print(f"Error loading taxonomy CSV {taxonomy_csv}: {e}", level='error')
        raise


def get_expected_taxonomy_from_filename(filename, taxonomy_mapping):
    ''' Grabs the expected taxonomy from a file name using the taxonomy mapping dictionary and returns family, full_taxonomy, and matched ID.'''
    filename_base = Path(filename).name
    
    # Try to match each ID in the mapping against the filename
    for sample_id, tax_data in taxonomy_mapping.items():
        if str(sample_id) in filename_base:
            return tax_data['family'], tax_data['full_taxonomy'], sample_id
    
    return None, None, None
    
def create_taxonomy_validation_summary(output_dir, taxonomy_mapping, assembly_dir, summary_csv='taxonomy_validation_summary.csv'):
    '''Creates a taxonomy validation summary CSV file based on parsed BLAST output files.
    Functions by checking contig hits against expected taxonomy starting from family level,
    with fallbacks to genus and species levels as needed.'''
    
    # Find all parsed blast output files (*_hits.tsv) in output_dir
    parsed_files = list(Path(output_dir).glob("*hits.tsv"))
    
    if not parsed_files:
        log_and_print("No parsed BLAST output files found. Skipping taxonomy validation summary creation.", level='warning')
        return None
    
    summary_results = []
    processed_count = 0
    
    # For each parsed file, extract expected taxonomy and compare with BLAST hits
    for parsed_file in sorted(parsed_files):
        # Extract original filename to get expected taxonomy
        original_name = re.sub(r'-(top_\d+_hits|all_hits)$', '', parsed_file.stem) + '.tsv'
        expected_family, expected_full_taxonomy, matched_id = get_expected_taxonomy_from_filename(original_name, taxonomy_mapping)
        
        log_and_print(f"  Processing {parsed_file.name}...")
        
        # Extract expected genus and species from taxonomy mapping
        expected_genus = None
        expected_species = None
        if matched_id and matched_id in taxonomy_mapping:
            # Look for genus in the full taxonomy string
            full_tax = taxonomy_mapping[matched_id].get('full_taxonomy', '')
            if ';g__' in full_tax:
                expected_genus = full_tax.split(';g__')[1].split(';')[0]
                expected_genus = re.sub("_", " ", expected_genus)
                expected_genus = expected_genus.lower().strip()
            if ';s__' in full_tax:
                expected_species = full_tax.split(';s__')[1].split(';')[0]
                expected_species = re.sub("_", " ", expected_species)
                expected_species = expected_species.lower().strip()
        
        # Determine display value for expected_family column
        expected_family_display = expected_family if expected_family else ''
        if not expected_family and expected_genus:
            expected_family_display = f"{expected_genus} (genus)"
        if not expected_family and not expected_genus and expected_species:
            expected_family_display = f"{expected_species} (species)"

        # Count total contigs in the assembly file (input to BLAST)
        n_contigs_in_assembly = count_contigs_in_assembly(assembly_dir, matched_id)

        # Initialise result dictionary
        result = {
            'ID': matched_id if matched_id else '',
            'expected_family': expected_family_display,
            'expected_taxonomy': expected_full_taxonomy if expected_full_taxonomy else '',
            'n_contigs_in': n_contigs_in_assembly if n_contigs_in_assembly > 0 else 'NA',
            'n_contigs_hits': 0,
            'contigs': '',
            'hit_taxonomy': 'NA',
            'its_coords': '',
            'correct_taxonomy': 'FAIL',
            'contig_path': ''
        }
        
        # Read the parsed BLAST output tsv file
        try:
            df = pd.read_csv(parsed_file, sep='\t')
            if df.empty:
                log_and_print(f"    -> No hits found in file -> FAIL")
                summary_results.append(result)
                processed_count += 1
                continue

            # Gather all families found in BLAST results
            families_found = set()

            # First, collect all families present in the BLAST results
            for _, row in df.iterrows():
                sseqid = row['sseqid']
                family = extract_family_from_taxonomy(sseqid)
                if family:
                    families_found.add(family)

            # If expected_family is empty OR not found in BLAST results use genus fallback
            invalid_family = (not expected_family) or (expected_family and expected_family not in families_found)
            use_genus_fallback = invalid_family and expected_genus

            if use_genus_fallback:
                log_and_print(f"    -> Expected family is empty, using genus-level fallback")
                log_and_print(f"    -> Expected genus: {expected_genus}")
                log_and_print(f"    -> Expected species: {expected_species}")
                
                # GENUS-LEVEL FALLBACK LOGIC
                # First pass: find all contigs that match the expected genus
                genus_matching_contigs = {}
                genera_found = set()
                
                for _, row in df.iterrows():
                    qseqid = row['qseqid']
                    sseqid = row['sseqid']
                    genus = extract_genus_from_taxonomy(sseqid)
                    species = extract_species_from_taxonomy(sseqid)
                    
                    if genus:
                        genera_found.add(genus)
                    
                    # If this contig matches the expected genus
                    if expected_genus and genus and genus == expected_genus:
                        if qseqid not in genus_matching_contigs:
                            # Store the first genus match for this contig
                            genus_matching_contigs[qseqid] = {
                                'sseqid': sseqid,
                                'genus': genus,
                                'species': species,
                                'full_taxonomy': extract_full_taxonomy_from_sseqid(sseqid),
                                'qstart': row.get('qstart', ''),
                                'qend': row.get('qend', '')
                            }
                        else:
                            # Check if this is a better match (species also matches expected species)
                            stored_species = genus_matching_contigs[qseqid]['species']
                            
                            # Update if the new hit matches expected species and stored one doesn't
                            should_update = False
                            
                            if expected_species and species and species == expected_species:
                                # New hit matches species
                                if not (stored_species and stored_species == expected_species):
                                    should_update = True
                            
                            if should_update:
                                genus_matching_contigs[qseqid] = {
                                    'sseqid': sseqid,
                                    'genus': genus,
                                    'species': species,
                                    'full_taxonomy': extract_full_taxonomy_from_sseqid(sseqid),
                                    'qstart': row.get('qstart', ''),
                                    'qend': row.get('qend', '')
                                }
                
                log_and_print(f"    -> DEBUG: Genera found in BLAST results: {sorted(genera_found)}")
                log_and_print(f"    -> Contigs matching expected genus: {len(genus_matching_contigs)}")
                
                if len(genus_matching_contigs) == 0:
                    # No contigs match expected genus -> FAIL
                    result['correct_taxonomy'] = 'FAIL'
                    log_and_print(f"    -> No contigs match expected genus -> FAIL")
                    
                elif len(genus_matching_contigs) == 1:
                    # Exactly one contig matches expected genus -> PASS
                    qseqid = list(genus_matching_contigs.keys())[0]
                    hit_data = genus_matching_contigs[qseqid]
                    
                    result['n_contigs_hits'] = 1
                    result['contigs'] = qseqid
                    result['correct_taxonomy'] = 'PASS'
                    
                    if hit_data['full_taxonomy']:
                        result['hit_taxonomy'] = hit_data['full_taxonomy']
                    
                    if hit_data['qstart'] and hit_data['qend']:
                        result['its_coords'] = f"{hit_data['qstart']}-{hit_data['qend']}"
                    
                    log_and_print(f"    -> Exactly 1 contig matches expected genus -> PASS")
                    
                else:
                    # Multiple contigs match expected genus -> check species level
                    log_and_print(f"    -> Multiple contigs ({len(genus_matching_contigs)}) match expected genus")
                    
                    if not expected_species:
                        result['correct_taxonomy'] = 'FAIL'
                        result['n_contigs_hits'] = len(genus_matching_contigs)
                        result['contigs'] = ';'.join(sorted(genus_matching_contigs.keys()))
                        log_and_print(f"    -> No expected species data available -> FAIL")
                    else:
                        # Check which contigs also match the expected species
                        species_matching_contigs = {}
                        
                        for qseqid, hit_data in genus_matching_contigs.items():
                            hit_species = hit_data['species']
                            if hit_species and hit_species == expected_species:
                                species_matching_contigs[qseqid] = hit_data
                        
                        log_and_print(f"    -> Contigs matching expected species '{expected_species}': {len(species_matching_contigs)}")
                        
                        if len(species_matching_contigs) == 1:
                            # Exactly one contig matches genus + species -> PASS
                            qseqid = list(species_matching_contigs.keys())[0]
                            hit_data = species_matching_contigs[qseqid]
                            
                            result['n_contigs_hits'] = 1
                            result['contigs'] = qseqid
                            result['correct_taxonomy'] = 'PASS'
                            
                            if hit_data['full_taxonomy']:
                                result['hit_taxonomy'] = hit_data['full_taxonomy']
                            
                            if hit_data['qstart'] and hit_data['qend']:
                                result['its_coords'] = f"{hit_data['qstart']}-{hit_data['qend']}"
                            
                            log_and_print(f"    -> Exactly 1 contig matches genus + species -> PASS")
                        
                        else:
                            # 0 or >1 contigs match genus + species -> FAIL
                            result['correct_taxonomy'] = 'FAIL'
                            result['n_contigs_hits'] = len(species_matching_contigs)
                            result['contigs'] = ';'.join(sorted(species_matching_contigs.keys()))
                            
                            if species_matching_contigs:
                                best_hit = list(species_matching_contigs.values())[0]
                                if best_hit['full_taxonomy']:
                                    result['hit_taxonomy'] = best_hit['full_taxonomy']
                            
                            log_and_print(f"    -> {len(species_matching_contigs)} contigs match genus + species -> FAIL")
            
            else:
                # STANDARD FAMILY-LEVEL LOGIC
                # First pass: find all contigs that match the expected family
                family_matching_contigs = {}
                families_found = set()
                
                for _, row in df.iterrows():
                    qseqid = row['qseqid']
                    sseqid = row['sseqid']
                    family = extract_family_from_taxonomy(sseqid)
                    genus = extract_genus_from_taxonomy(sseqid)
                    species = extract_species_from_taxonomy(sseqid)
                    
                    if family:
                        families_found.add(family)
                    
                    # If this contig matches the expected family
                    if expected_family and family and family == expected_family:
                        if qseqid not in family_matching_contigs:
                            # Store the first family match for this contig
                            family_matching_contigs[qseqid] = {
                                'sseqid': sseqid,
                                'family': family,
                                'genus': genus,
                                'species': species,
                                'full_taxonomy': extract_full_taxonomy_from_sseqid(sseqid),
                                'qstart': row.get('qstart', ''),
                                'qend': row.get('qend', '')
                            }
                        else:
                            # Check if this is a better match (genus and/or species also matches)
                            stored_genus = family_matching_contigs[qseqid]['genus']
                            stored_species = family_matching_contigs[qseqid]['species']
                            
                            # Update if the new hit is "better":
                            # - New hit matches expected species and stored one doesn't, OR
                            # - New hit matches expected genus and stored one doesn't (and no species match exists)
                            should_update = False
                            
                            if expected_species and species and species == expected_species:
                                # New hit matches species
                                if not (stored_species and stored_species == expected_species):
                                    should_update = True
                            elif expected_genus and genus and genus == expected_genus:
                                # New hit matches genus (but not species)
                                if not (stored_genus and stored_genus == expected_genus):
                                    should_update = True
                            
                            if should_update:
                                family_matching_contigs[qseqid] = {
                                    'sseqid': sseqid,
                                    'family': family,
                                    'genus': genus,
                                    'species': species,
                                    'full_taxonomy': extract_full_taxonomy_from_sseqid(sseqid),
                                    'qstart': row.get('qstart', ''),
                                    'qend': row.get('qend', '')
                                }
                
                # Debug logging for family matches
                log_and_print(f"    -> Expected family: {expected_family}")
                log_and_print(f"    -> Expected genus: {expected_genus}")
                log_and_print(f"    -> Expected species: {expected_species}")
                log_and_print(f"    -> Families found in BLAST results: {sorted(families_found)}")
                log_and_print(f"    -> Contigs matching expected family: {len(family_matching_contigs)}")
                
                # Now apply the multi-level taxonomy logic
                if len(family_matching_contigs) == 0:
                    # No contigs match expected family -> FAIL
                    result['correct_taxonomy'] = 'FAIL'
                    log_and_print(f"    -> No contigs match expected family -> FAIL")
                    
                elif len(family_matching_contigs) == 1:
                    # Exactly one contig matches expected family -> PASS
                    qseqid = list(family_matching_contigs.keys())[0]
                    hit_data = family_matching_contigs[qseqid]
                    log_and_print(f"    == Parsing contig ID: {qseqid}")

                    result['n_contigs_hits'] = 1
                    result['contigs'] = qseqid
                    result['correct_taxonomy'] = 'PASS'
                    
                    if hit_data['full_taxonomy']:
                        result['hit_taxonomy'] = hit_data['full_taxonomy']
                    
                    # Store coordinates
                    if hit_data['qstart'] and hit_data['qend']:
                        result['its_coords'] = f"{hit_data['qstart']}-{hit_data['qend']}"
                    
                    log_and_print(f"    -> Exactly 1 contig matches expected family -> PASS")
                    
                else:
                    # Multiple contigs match expected family -> check genus level
                    log_and_print(f"    -> Multiple contigs ({len(family_matching_contigs)}) match expected family")
                    
                    # If we don't have expected genus data, fall back to family-only logic (FAIL)
                    if not expected_genus:
                        result['correct_taxonomy'] = 'FAIL'
                        result['n_contigs_hits'] = len(family_matching_contigs)
                        result['contigs'] = ';'.join(sorted(family_matching_contigs.keys()))
                        log_and_print(f"    -> No expected genus data available -> FAIL")
                    else:
                        # Check which contigs also match the expected genus
                        genus_matching_contigs = {}
                        
                        for qseqid, hit_data in family_matching_contigs.items():
                            hit_genus = hit_data['genus']
                            if hit_genus and hit_genus == expected_genus:
                                genus_matching_contigs[qseqid] = hit_data
                        
                        log_and_print(f"    -> Contigs matching expected genus '{expected_genus}': {len(genus_matching_contigs)}")
                        
                        if len(genus_matching_contigs) == 1:
                            # Exactly one contig matches both family and genus -> PASS
                            qseqid = list(genus_matching_contigs.keys())[0]
                            hit_data = genus_matching_contigs[qseqid]
                            
                            result['n_contigs_hits'] = 1
                            result['contigs'] = qseqid
                            result['correct_taxonomy'] = 'PASS'
                            
                            if hit_data['full_taxonomy']:
                                result['hit_taxonomy'] = hit_data['full_taxonomy']
                            
                            # Store coordinates
                            if hit_data['qstart'] and hit_data['qend']:
                                result['its_coords'] = f"{hit_data['qstart']}-{hit_data['qend']}"
                            
                            log_and_print(f"    -> Exactly 1 contig matches both family and genus -> PASS")
                        
                        elif len(genus_matching_contigs) == 0:
                            # 0 contigs match both family and genus -> FAIL
                            result['correct_taxonomy'] = 'FAIL'
                            result['n_contigs_hits'] = len(family_matching_contigs)
                            result['contigs'] = ';'.join(sorted(family_matching_contigs.keys()))
                            
                            # Use the best hit from family-matching contigs for hit_taxonomy
                            if family_matching_contigs:
                                best_hit = list(family_matching_contigs.values())[0]
                                if best_hit['full_taxonomy']:
                                    result['hit_taxonomy'] = best_hit['full_taxonomy']
                            
                            log_and_print(f"    -> 0 contigs match both family and genus -> FAIL")
                        
                        else:
                            # >1 contigs match both family and genus -> check species level
                            log_and_print(f"    -> Multiple contigs ({len(genus_matching_contigs)}) match both family and genus")
                            
                            # If we don't have expected species data, fall back (FAIL)
                            if not expected_species:
                                result['correct_taxonomy'] = 'FAIL'
                                result['n_contigs_hits'] = len(genus_matching_contigs)
                                result['contigs'] = ';'.join(sorted(genus_matching_contigs.keys()))
                                log_and_print(f"    -> No expected species data available -> FAIL")
                            else:
                                # Check which contigs also match the expected species
                                species_matching_contigs = {}
                                
                                for qseqid, hit_data in genus_matching_contigs.items():
                                    hit_species = hit_data['species']
                                    if hit_species and hit_species == expected_species:
                                        species_matching_contigs[qseqid] = hit_data
                                
                                log_and_print(f"    -> Contigs matching expected species '{expected_species}': {len(species_matching_contigs)}")
                                
                                if len(species_matching_contigs) == 1:
                                    # Exactly one contig matches family + genus + species -> PASS
                                    qseqid = list(species_matching_contigs.keys())[0]
                                    hit_data = species_matching_contigs[qseqid]
                                    
                                    result['n_contigs_hits'] = 1
                                    result['contigs'] = qseqid
                                    result['correct_taxonomy'] = 'PASS'
                                    
                                    if hit_data['full_taxonomy']:
                                        result['hit_taxonomy'] = hit_data['full_taxonomy']
                                    
                                    # Store coordinates
                                    if hit_data['qstart'] and hit_data['qend']:
                                        result['its_coords'] = f"{hit_data['qstart']}-{hit_data['qend']}"
                                    
                                    log_and_print(f"    -> Exactly 1 contig matches family + genus + species -> PASS")
                                
                                else:
                                    # 0 or >1 contigs match family + genus + species -> FAIL
                                    result['correct_taxonomy'] = 'FAIL'
                                    result['n_contigs_hits'] = len(genus_matching_contigs)
                                    result['contigs'] = ';'.join(sorted(genus_matching_contigs.keys()))
                                    
                                    # Use the best hit from genus-matching contigs for hit_taxonomy
                                    if genus_matching_contigs:
                                        best_hit = list(genus_matching_contigs.values())[0]
                                        if best_hit['full_taxonomy']:
                                            result['hit_taxonomy'] = best_hit['full_taxonomy']
                                    
                                    log_and_print(f"    -> {len(species_matching_contigs)} contigs match family + genus + species -> FAIL")
            
            processed_count += 1
            
        except Exception as e:
            log_and_print(f"    -> Error processing {parsed_file}: {e}", level='error')
        

        summary_results.append(result)
        
    # Create summary CSV
    summary_file = Path(output_dir) / summary_csv  # Use custom filename
    summary_df = pd.DataFrame(summary_results)
    
    # Define column order
    column_order = ['ID', 'expected_family', 'expected_taxonomy', 'n_contigs_in', 'n_contigs_hits', 'contigs', 'hit_taxonomy', 'its_coords', 'correct_taxonomy', 'contig_path']
    summary_df = summary_df[column_order]
    
    # Save to CSV
    summary_df.to_csv(summary_file, index=False)
    
    # Add contig paths if assembly_dir was provided
    if assembly_dir and Path(assembly_dir).is_dir():
        for idx, row in summary_df.iterrows():
            if row['correct_taxonomy'] == 'PASS' and row['ID']:
                contig_fasta = Path(output_dir) / f"{row['ID']}_parsed_contig.fasta"
                if contig_fasta.exists():
                    summary_df.at[idx, 'contig_path'] = str(contig_fasta.resolve())
        summary_df.to_csv(summary_file, index=False)
    
    log_and_print(f"  Processed {processed_count} files")
    log_and_print(f"  Saved taxonomy validation summary to: {summary_file}")
    
    # Summary statistics
    total_pass = len(summary_df[summary_df['correct_taxonomy'] == 'PASS'])
    total_fail = len(summary_df[summary_df['correct_taxonomy'] == 'FAIL'])
    total_with_taxonomy = len(summary_df[summary_df['expected_family'] != ''])
    
    log_and_print(f"  Summary statistics:")
    log_and_print(f"    -> Total samples: {len(summary_df)}")
    log_and_print(f"    -> Samples with taxonomy mapping: {total_with_taxonomy}")
    log_and_print(f"    -> PASS (exactly 1 contig hit): {total_pass}")
    log_and_print(f"    -> FAIL (0 or >1 contig hits): {total_fail}")
    
    return summary_file  # Return the path to summary file
    
def process_blast_file(input_file, output_dir, top_n=None, min_len=None, min_pident=None):
    '''Filter a single BLAST TSV file based on length and pident,
    then extract top N hits per query sequence.'''

    try:
        # Detect file format
        log_and_print(f"Processing {input_file}...")
        has_header, column_names = detect_blast_format(input_file)
        
        # Read the BLAST results
        if has_header:
            df = pd.read_csv(input_file, sep='\t', header=0)
            log_and_print(f"  -> Reading file with headers")
        else:
            df = pd.read_csv(input_file, sep='\t', header=None, names=column_names)
            log_and_print(f"  -> Reading file without headers")
        
        # Ensure we have the minimum required columns
        required_columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue']
        missing_columns = [col for col in required_columns if col not in df.columns]
        
        if missing_columns:
            log_and_print(f"Error: Missing required columns {missing_columns} in {input_file}", level='error')
            return 0, 0
        
        original_hits = len(df)
        
        # Apply minimum length filter if specified
        if min_len is not None:
            df_filtered = df[df['length'] >= min_len]
            filtered_out = len(df) - len(df_filtered)
            log_and_print(f"  -> Applied min_len filter ({min_len}): {filtered_out} hits filtered out ({len(df_filtered)} remaining)")
            df = df_filtered

        # Apply minimum percent identity filter if specified
        if min_pident is not None:
            df_filtered = df[df['pident'] >= min_pident]
            filtered_out = len(df) - len(df_filtered)
            log_and_print(f"  -> Applied min_pident filter ({min_pident}): {filtered_out} hits filtered out ({len(df_filtered)} remaining)")
            df = df_filtered
        
        # Sort by qseqid first, then by our ranking criteria
        # pident (desc), length (desc), evalue (asc)
        df_sorted = df.sort_values([
            'qseqid',           # Group by query sequence
            'length',           # Longer alignments are better
            'evalue'            # Lower e-value is better
        ], ascending=[True, False, True])
        
        # Group by qseqid and take top N hits for each (or all hits if top_n is None)
        if top_n is not None:
            top_hits = df_sorted.groupby('qseqid').head(top_n)
        else:
            top_hits = df_sorted
        
        # Select only the columns we want in the output
        output_columns = ['qseqid', 'sseqid', 'pident', 'length', 'evalue', 'qstart', 'qend']
        result = top_hits[output_columns]
        
        # Generate output filename
        input_name = Path(input_file).stem
        if top_n is not None:
            output_file = Path(output_dir) / f"{input_name}-top_{top_n}_hits.tsv"
        else:
            output_file = Path(output_dir) / f"{input_name}-all_hits.tsv"
        
        # Save the results with column headers
        result.to_csv(output_file, sep='\t', index=False)
        
        log_and_print(f"  -> Saved {len(result)} hits from {df['qseqid'].nunique()} queries to {output_file}")
        
        return len(result), df['qseqid'].nunique()
        
    except Exception as e:
        log_and_print(f"Error processing {input_file}: {e}", level='error')
        return 0, 0
 
def open_fasta(filename):
    filename = str(filename)  # Ensure filename is a string
    if filename.endswith('.gz'):
        return gzip.open(filename, 'rt')
    return open(filename, 'r')

def grep_fasta_header(fasta_file, patterns, output_file):
    """Search FASTA headers for patterns and write matching records.
    
    Args:
        fasta_file (str): Path to input FASTA file
        patterns (list): List of patterns to search for in headers
        output_file (str): Path to output file
    """

    log_and_print(f"Searching for {patterns} in {fasta_file}...")
    with open_fasta(fasta_file) as fa, open(output_file, 'w') as out:
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
            log_and_print(f"Writing = {write_seq} for {patterns}")
            out.writelines(buffer)
         
def extract_passing_contigs(output_dir, assembly_dir, taxonomy_mapping, summary_csv='taxonomy_validation_summary.csv'):
    '''Extracts FASTA sequences for contigs that passed taxonomy validation.'''
    
    log_and_print("\nExtracting FASTA sequences for passed contigs...")
    
    if not assembly_dir or not Path(assembly_dir).is_dir():
        log_and_print(f"Assembly directory '{assembly_dir}' does not exist or is not a directory. Skipping FASTA extraction.", level='warning')
        return
    
    # Ensure output_dir is a Path object
    output_dir = Path(output_dir).resolve()
    
    # Read the summary file
    summary_file = output_dir / summary_csv  # Use the passed filename
    
    if not summary_file.exists():
        log_and_print(f"Taxonomy validation summary file not found at: {summary_file}", level='warning')
        return
    
    df = pd.read_csv(summary_file)
    dfkeep = df.loc[df["correct_taxonomy"] == "PASS"]
    
    extracted_count = 0
    for _, row in dfkeep.iterrows():
        filename = row["ID"]
        if not filename or pd.isna(filename):
            continue
            
        contig_patterns = [row["contigs"]]
        output_file = output_dir / f"{filename}_parsed_contig.fasta"

        assembly_file = Path(assembly_dir) / f"{filename}.spades.out" / "scaffolds.fasta"
        if not assembly_file.is_file():
            log_and_print(f"  -> Assembly file not found for {filename}", level='warning')
            continue
        
        grep_fasta_header(assembly_file, contig_patterns, output_file)
        if output_file.exists():
            extracted_count += 1
            log_and_print(f"  -> Extracted contig for {filename}")
    
    log_and_print(f"  Extracted FASTA sequences for {extracted_count} passed samples")
    return extracted_count

def update_summary_with_contig_paths(output_dir, summary_csv='taxonomy_validation_summary.csv'):
    '''Updates the taxonomy validation summary CSV with absolute contig file paths.'''
    
    # Convert output_dir to absolute path first
    output_dir = Path(output_dir).resolve()
    
    summary_file = output_dir / summary_csv  # Use passed filename
    if not summary_file.exists():
        return
    
    df = pd.read_csv(summary_file)
    
    # Ensure contig_path is string dtype to avoid FutureWarning
    df['contig_path'] = df['contig_path'].astype(str)
    
    updated_count = 0
    for idx, row in df.iterrows():
        if row['correct_taxonomy'] == 'PASS' and row['ID']:
            contig_fasta = output_dir / f"{row['ID']}_parsed_contig.fasta"
            if contig_fasta.exists():
                df.at[idx, 'contig_path'] = str(contig_fasta)
                updated_count += 1
    
    # Save updated CSV
    df.to_csv(summary_file, index=False)
    log_and_print(f"  Updated {updated_count} contig paths in summary CSV")    
    
def main():
    parser = argparse.ArgumentParser(description='Process BLAST TSV files to extract top hits per query')
    parser.add_argument('--input_dir', required=True, help='Directory containing BLAST TSV files (required)')
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
    parser.add_argument('--log_file', default=None,
                       help='Log file path (default: blast_processor_YYYYMMDD_HHMMSS.log)')
    parser.add_argument('--assembly_dir', required=True,
                       help='Directory containing assemblies (*.spades.out) for contig counting and FASTA extraction')
    parser.add_argument('--id_column', default='ID',
                   help='Column name in taxonomy CSV that contains the IDs matching filenames (default: ID)')
    parser.add_argument('--summary_csv', default='taxonomy_validation_summary.csv',
                   help='Name for the taxonomy validation summary CSV file (default: taxonomy_validation_summary.csv)')

    args = parser.parse_args()
    
    # Set up log file with timestamp if not provided
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'blast_processor_{timestamp}.log'

    # Set up logging
    logger = setup_logging(log_file=args.log_file)

    # Load taxonomy mapping
    try:
        taxonomy_mapping = load_taxonomy_mapping(args.taxonomy_csv, args.id_column)
    except Exception:
        log_and_print("Failed to load taxonomy CSV. Exiting.", level='error')
        return
    
    # Create output directory if it doesn't exist
    output_dir = Path(args.output_dir)
    
    # Check if output_dir exists as a file (not directory) and handle it
    if output_dir.exists() and not output_dir.is_dir():
        log_and_print(f"Warning: '{output_dir}' exists as a file, not a directory. Removing file and creating directory.", level='warning')
        output_dir.unlink()  # Remove the file
    
    output_dir.mkdir(exist_ok=True)
    
    # Find all BLAST TSV files recursively
    input_dir = Path(args.input_dir)
    tsv_files = list(input_dir.glob('*.tsv'))
    
    if not tsv_files:
        log_and_print(f"No TSV files found in {input_dir} or its subdirectories", level='error')
        return
        
    log_and_print(f"BLAST Results Processor Started")
    log_and_print(f"Input directory: {input_dir}")
    log_and_print(f"Output directory: {output_dir}")
    log_and_print(f"Taxonomy CSV: {args.taxonomy_csv}")
    log_and_print(f"Assembly directory: {args.assembly_dir}")
    log_and_print(f"Log file: {args.log_file}")
    log_and_print(f"Found {len(tsv_files)} TSV files to process")
    if args.top_n is not None:
        log_and_print(f"Keeping top {args.top_n} hits per query sequence")
    else:
        log_and_print(f"Keeping all hits per query sequence")
    if args.min_len is not None:
        log_and_print(f"Minimum alignment length filter: {args.min_len}")
    if args.min_pident is not None:
        log_and_print(f"Minimum percent identity filter: {args.min_pident}")
    log_and_print("-" * 60)
    
    total_hits = 0
    total_queries = 0
    processed_files = 0
    file_results = []  # Store individual file results
    
    # Process each file
    for tsv_file in sorted(tsv_files):
        # Get taxonomy info for this file
        expected_family, expected_full_taxonomy, matched_id = get_expected_taxonomy_from_filename(tsv_file.name, taxonomy_mapping)
        
        hits, queries = process_blast_file(tsv_file, output_dir, args.top_n, args.min_len, args.min_pident)

        total_hits += hits
        total_queries += queries
        
        if hits > 0:
            processed_files += 1
        
        # Store individual file result
        file_results.append({
            'filename': tsv_file.name,
            'matched_id': matched_id,
            'expected_family': expected_family,
            'expected_full_taxonomy': expected_full_taxonomy,
            'hits': hits,
            'queries': queries,
            'processed': hits > 0
        })
        
    log_and_print("-" * 60)
    log_and_print(f"Processing complete!")
    log_and_print(f"Files processed: {processed_files}/{len(tsv_files)}")
    log_and_print(f"Total queries processed: {total_queries}")
    log_and_print(f"Total hits extracted: {total_hits}")
    
    # Print individual file results
    log_and_print("=" * 60)
    log_and_print("PROCESSING RESULTS BY TAXONOMY ID:")
    log_and_print("-" * 60)
    
    # Group results by matched vs unmatched
    matched_files = [r for r in file_results if r['matched_id'] is not None]
    unmatched_files = [r for r in file_results if r['matched_id'] is None]
    
    # Show matched files
    for result in matched_files:
        status = "Ã¢Å“â€œ" if result['processed'] else "Ã¢Å“â€”"
        log_and_print(f"   {result['matched_id']} {status}")
        log_and_print(f"-> File: {result['filename']}")
        log_and_print(f"   Expected family: {result['expected_family']}")
        log_and_print(f"   Expected full taxonomy: {result['expected_full_taxonomy']}")
        log_and_print(f"   Queries (contigs) processed: {result['queries']}")
        log_and_print(f"   Hits extracted: {result['hits']}")
    
    # Show unmatched files
    if unmatched_files:
        log_and_print(f"\nUnmatched Files:")
        for result in unmatched_files:
            status = "Ã¢Å“â€œ" if result['processed'] else "Ã¢Å“â€”"
            log_and_print(f"  {result['filename']} {status} - No taxonomy mapping found")
            if result['processed']:
                log_and_print(f"    Queries processed: {result['queries']}, Hits extracted: {result['hits']}") 

    # Create taxonomy validation summary CSV
    log_and_print("=" * 60)
    log_and_print("CREATING TAXONOMY VALIDATION SUMMARY:")
    log_and_print("-" * 60)
    summary_file = create_taxonomy_validation_summary(output_dir, taxonomy_mapping, args.assembly_dir, args.summary_csv)
    
    log_and_print(f"\nLog saved to: {args.log_file}")
    
    # Extract FASTAs of passed contigs
    if summary_file:
        log_and_print("\n" + "=" * 60)
        log_and_print("CREATING FASTAS OF ALL 'PASS' CONTIGS:")
        log_and_print("=" * 60)
        
        # Use the existing extract_passing_contigs function
        extracted_count = extract_passing_contigs(output_dir, args.assembly_dir, taxonomy_mapping, args.summary_csv)
        
        if extracted_count:
            # Update the summary file with contig paths
            update_summary_with_contig_paths(output_dir, args.summary_csv)

if __name__ == "__main__":
    main()
