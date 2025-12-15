#!/usr/bin/env python3
"""
UNITE to NCBI taxonomy matcher & sequence retriever

This script replaces go_fetch.py functionality by:
1. Parsing UNITE FASTA files to extract taxonomy information
2. Getting NCBI lineage information for input taxids
3. Matching NCBI taxonomy against UNITE taxonomy with fallback strategy
4. Extracting the longest corresponding sequence for each taxon
5. Creating output files compatible with skim2ribo pipeline

Usage:
    python UNITEd.py --tracking_sheet samples.csv --unite_db unite.fasta --output output_dir --email your@email.com --api your_api_key
    
Required Arguments:
--tracking_sheet: CSV file with ID and taxid columns (see CSV Format section below)
--unite_db: Local UNITE FASTA database file
--output: Output directory for results
--email: Email address (required by NCBI Entrez API)

Optional Arguments:
--api: NCBI API key (increases query limits)
--all: Retrieve ALL sequences for each match instead of just the longest (with fallback, (species->genus->family->order->class->phylum)
--number N: Retrieve exactly N sequences for each match (overrides --all)
--tax_rank RANK: Search only at specific rank (species, genus, family, order, class, phylum) without fallback
--traverse RANK: Traverse up to specified rank if '--number N' cannot be found at the specified '--tax_rank [rank]' (requires both --tax_rank and --number)
--diversity: Distribute sequences across child taxa, prioritising target species matches first (requires --number and --tax_rank)
--summary_csv FILENAME: Generate detailed CSV summary with specified filename (placed in output directory)
--version: Show version information

Input CSV Format: 
ID Column (sample identifiers) accepts any of: id, ID, Id, sample_id, sample_ID, Sample_ID, SampleID, sample, Sample
Taxid Column (NCBI taxonomy IDs) accepts any of: taxid, taxID, TaxID, TAXID, tax_id, Tax_ID, TAX_ID, taxonomy_id, ncbi_taxid, NCBI_TaxID

Output: 
Creates files named {sample_id}_seed.fasta in the specified output directory using a flat directory structure.
"""

import os
import sys
import csv
import argparse
import logging
import re
import time
import urllib
import random
import traceback
from collections import defaultdict
from Bio import Entrez
from datetime import datetime

from metahist_tools import setup_logging

# Increase time between and number of tries used by entrez (from go_fetch.py)
Entrez.sleep_between_tries = 20 
Entrez.max_tries = 20


def setup_argument_parser():
    """
    Set up argument parser with traverse and diversity functionality.
    """
    parser = argparse.ArgumentParser(
        description="Extract sequences from UNITE FASTA based on NCBI taxonomy matching",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Example usage:
    python unite_fetch.py --tracking_sheet samples.csv --unite_db unite_db.fasta --output results --email your@email.com --api your_api_key

    # Use a custom ID column name:
    python unite_fetch.py --tracking_sheet samples.csv --unite_db unite_db.fasta --output results --email your@email.com --id-header "MyCustomID"

    # Get 20 sequences distributed across species within a genus:
    python unite_fetch.py --tracking_sheet samples.csv --unite_db unite_db.fasta --output results --email your@email.com --api your_api_key --tax_rank genus --number 20 --diversity

    # Get 20 sequences with diversity, traverse up to family if insufficient child taxa:
    python unite_fetch.py --tracking_sheet samples.csv --unite_db unite_db.fasta --output results --email your@email.com --api your_api_key --tax_rank genus --number 20 --diversity --traverse family

    # Get 20 sequences starting at genus, traverse up to order if needed (no diversity):
    python unite_fetch.py --tracking_sheet samples.csv --unite_db unite_db.fasta --output results --email your@email.com --api your_api_key --tax_rank genus --number 20 --traverse order

Diversity functionality:
    --diversity: When used with --tax_rank and --number, distributes sequences across child taxa.
                 Prioritises finding target species matches first (up to 25% of --number), then
                 distributes remaining sequences across other child taxa.
                 
    Example: "--tax_rank genus --number 20 --diversity"
    - First checks if target species exists in the genus (takes up to 5 sequences from it)
    - Finds remaining species within the genus
    - Distributes remaining 15 sequences proportionally across those species
    - If 4 other species exist, tries to get ~3-4 sequences from each
    - Prioritises even distribution over total sequence length
        """
    )
    
    parser.add_argument("--tracking_sheet", required=True, 
                        help="CSV file with 'ID' and 'taxid' columns")
    parser.add_argument("--unite_db", required=True, 
                        help="UNITE FASTA file with taxonomy information")
    parser.add_argument("--output", required=True, 
                        help="Output directory for results")
    parser.add_argument("--email", required=True, 
                        help="Email address for NCBI Entrez API")
    parser.add_argument("--api", default="None", 
                        help="NCBI API key (optional)")
    parser.add_argument("--all", action="store_true", 
                        help="Retrieve ALL sequences for each match instead of just the longest (requires more memory)")
    parser.add_argument("--number", type=int, 
                        help="Retrieve exactly N sequences for each match (overrides --all flag)")
    parser.add_argument("--tax_rank", 
                        choices=['species', 'genus', 'family', 'order', 'class', 'phylum'],
                        help="Search only at this specific taxonomic rank (no fallback to higher ranks unless --traverse used)")
    parser.add_argument("--traverse", 
                        choices=['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom'],
                        help="When used with --tax_rank and --number, traverse up taxonomic hierarchy "
                             "to this rank to collect target number of sequences if insufficient sequences found at specified rank "
                             "(e.g., --tax_rank genus --traverse family will search genus then family if needed)")
    parser.add_argument("--diversity", action="store_true",
                        help="When used with --tax_rank and --number, distribute sequences across child taxa "
                             "instead of taking longest sequences from parent taxon. Prioritises target species matches first "
                             "(up to 25%% of --number), then distributes remaining across other child taxa")
    parser.add_argument("--summary_csv", 
                        help="Generate detailed CSV summary with specified filename (will be placed in output directory)")
    parser.add_argument("--id-header", 
                        help="Specify the exact column name to use for sample IDs in the input CSV")
    parser.add_argument("--version", action="version", version="1.2.0")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `UNITEd_log` file will be created")

    return parser


def parse_fasta_header(header):
    """
    Parse a FASTA header to extract taxonomy information regardless of its position.
    Looks for patterns like k__Fungi;p__Ascomycota;... anywhere in the header.
    """
    # Remove the '>' character at the beginning
    header = header.lstrip('>')
    
    # Look taxonomy pattern in the header, e.g. k__Fungi;p__Ascomycota;...
    taxonomy_pattern = r'([a-z]__[^;|]+(?:;[a-z]__[^;|]+)*)'
    taxonomy_match = re.search(taxonomy_pattern, header)
    
    if not taxonomy_match:
        return None
    
    taxonomy_str = taxonomy_match.group(1)
    
    # Extract sequence accession (use the first part before | as default)
    id_match = re.match(r'^([^|]+)', header)
    id_str = id_match.group(1) if id_match else "unknown"
    
    # Parse the taxonomy string
    taxonomy = {}
    for tax_part in taxonomy_str.split(';'):
        if not tax_part:
            continue
        
        # Match patterns like k__Fungi, p__Ascomycota, etc.
        match = re.match(r'([a-z])__(.+)', tax_part)
        if match:
            rank_code, name = match.groups()
            taxonomy[rank_code] = name
    
    return {
        'id': id_str,
        'header': header,
        'taxonomy': taxonomy
    }


def build_unite_taxonomy_index(fasta_file, logger, keep_all_sequences=False):
    """
    Parse UNITE FASTA file and build taxonomy index for fast lookups.
    """
    if keep_all_sequences:
        taxonomy_index = defaultdict(lambda: defaultdict(list))
        logger.info("Building index with ALL sequences per taxon (higher memory usage)")
    else:
        taxonomy_index = defaultdict(dict)
        logger.info("Building index with LONGEST sequence per taxon (memory efficient)")
    
    # Define rank order
    rank_order = ['s', 'g', 'f', 'o', 'c', 'p', 'k']
    rank_names = {
        'k': 'kingdom',
        'p': 'phylum', 
        'c': 'class',
        'o': 'order',
        'f': 'family',
        'g': 'genus',
        's': 'species'
    }
    
    logger.info(f"Parsing UNITE FASTA file: {fasta_file}")
    
    sequences_processed = 0
    sequences_with_taxonomy = 0
    current_header = None
    current_sequence = []
    
    def process_sequence(header, sequence_str):
        """Process a single sequence and update index."""
        nonlocal sequences_processed, sequences_with_taxonomy
        
        sequences_processed += 1
        if sequences_processed % 10000 == 0:
            logger.info(f"Processed {sequences_processed} sequences...")
        
        taxonomy_info = parse_fasta_header(header)
        if not taxonomy_info or not taxonomy_info['taxonomy']:
            return
        
        sequences_with_taxonomy += 1
        sequence_length = len(sequence_str)
        
        # Index this sequence at all taxonomic ranks it has
        for rank_code in rank_order:
            if rank_code in taxonomy_info['taxonomy']:
                taxon_name = taxonomy_info['taxonomy'][rank_code]
                rank_name = rank_names[rank_code]
                
                # Create sequence info
                seq_info = {
                    'header': header,
                    'sequence': sequence_str,
                    'length': sequence_length,
                    'taxonomy_info': taxonomy_info,
                    'rank': rank_name,
                    'taxon_name': taxon_name
                }
                
                if keep_all_sequences:
                    # Store all sequences (higher memory usage)
                    taxonomy_index[rank_name][taxon_name].append(seq_info)
                else:
                    # Only keep the longest sequence (memory efficient)
                    if (taxon_name not in taxonomy_index[rank_name] or 
                        sequence_length > taxonomy_index[rank_name][taxon_name]['length']):
                        taxonomy_index[rank_name][taxon_name] = seq_info
    
    try:
        with open(fasta_file, 'r') as f:
            for line_num, line in enumerate(f, 1):
                line = line.strip()
                
                if line.startswith('>'):
                    # Process previous sequence if exists
                    if current_header is not None and current_sequence:
                        sequence_str = ''.join(current_sequence)
                        process_sequence(current_header, sequence_str)
                    
                    # Start new sequence
                    current_header = line
                    current_sequence = []
                    
                elif line:  # Non-empty sequence line
                    current_sequence.append(line)
            
            # Process the last sequence
            if current_header is not None and current_sequence:
                sequence_str = ''.join(current_sequence)
                process_sequence(current_header, sequence_str)
        
        logger.info(f"Processed {sequences_processed} sequences, {sequences_with_taxonomy} with taxonomy")
        
        # Log statistics
        for rank_name, taxa in taxonomy_index.items():
            if keep_all_sequences:
                total_sequences = sum(len(seq_list) for seq_list in taxa.values())
                logger.info(f"Found {len(taxa)} unique taxa at {rank_name} level ({total_sequences} total sequences)")
            else:
                logger.info(f"Found {len(taxa)} unique taxa at {rank_name} level")
        
        return taxonomy_index
        
    except Exception as e:
        logger.error(f"Error parsing UNITE FASTA file at line {line_num if 'line_num' in locals() else 'unknown'}: {str(e)}")
        import traceback
        logger.error(f"Full traceback: {traceback.format_exc()}")
        raise



def get_ncbi_lineage(taxid, email, logger, api_key=None):
    """
    Get NCBI taxonomic lineage for a given taxid with retry logic.
    """
    Entrez.email = email
    if api_key and api_key != "None":
        Entrez.api_key = api_key
    
    logger.info(f"Getting NCBI lineage for taxid: {taxid}")
    
    max_retries = 5
    base_delay = 10
    
    for attempt in range(max_retries):
        try:
            # Get taxonomy information
            handle = Entrez.efetch(db="Taxonomy", id=taxid, retmode="xml")
            record = Entrez.read(handle)
            
            if not record:
                raise ValueError(f"No taxonomy record found for taxid {taxid}")
            
            # Extract lineage information
            tax_record = record[0]
            lineage_list = tax_record.get("LineageEx", [])
            
            # Add the current taxon to the lineage
            lineage_list.append({
                'TaxId': tax_record['TaxId'],
                'ScientificName': tax_record['ScientificName'],
                'Rank': tax_record['Rank']
            })
            
            # Convert to our format
            lineage = {}
            for item in lineage_list:
                rank = item['Rank']
                name = item['ScientificName']
                
                # Map NCBI ranks
                if rank in ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']:
                    lineage[rank] = name
            
            logger.info(f"Retrieved lineage for {taxid}: {', '.join([f'{k}:{v}' for k,v in lineage.items()])}")
            return lineage
            
        except (ValueError, urllib.error.HTTPError, urllib.error.URLError) as e:
            delay = base_delay * (2 ** attempt) + random.uniform(0, 1)
            logger.warning(f"Attempt {attempt+1}/{max_retries}: Error getting lineage for {taxid}: {str(e)}")
            
            if attempt < max_retries - 1:
                logger.info(f"Retrying in {delay:.2f} seconds...")
                time.sleep(delay)
            else:
                logger.error(f"Failed to get lineage for {taxid} after {max_retries} attempts.")
                raise
    
    raise Exception(f"Failed to get lineage for {taxid}")


def get_most_specific_taxon(ncbi_lineage):
    """
    Get the most specific taxon name from NCBI lineage.
    Returns the name at the most specific rank available.
    """
    rank_order = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    
    for rank in rank_order:
        if rank in ncbi_lineage:
            return ncbi_lineage[rank]
    
    return "Unknown"


def search_rank_for_taxon(unite_index, rank, taxon_name):
    """
    Helper function to search for a taxon at a specific rank with various matching strategies.
    """
    # Look for exact matches first
    if taxon_name in unite_index[rank]:
        return unite_index[rank][taxon_name]
    
    # Look for case-insensitive matches
    taxon_lower = taxon_name.lower()
    for unite_taxon, sequences in unite_index[rank].items():
        if unite_taxon.lower() == taxon_lower:
            return sequences
    
    # Look for partial matches (species name without underscore/space differences)
    if rank == 'species':
        # Handle variations in species names (spaces vs underscores)
        taxon_normalised = re.sub(r'[_\s]+', ' ', taxon_name.lower().strip())
        for unite_taxon, sequences in unite_index[rank].items():
            unite_normalised = re.sub(r'[_\s]+', ' ', unite_taxon.lower().strip())
            if taxon_normalised == unite_normalised:
                return sequences
    
    return None


def get_rank_hierarchy():
    """Get the taxonomic rank hierarchy in order from most specific to most general."""
    return ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']


def get_child_rank(parent_rank):
    """Get the child rank for a given parent rank."""
    hierarchy = get_rank_hierarchy()
    try:
        parent_index = hierarchy.index(parent_rank)
        if parent_index > 0:
            return hierarchy[parent_index - 1]
        else:
            return None  # species has no child rank
    except ValueError:
        return None


def find_child_taxa(unite_index, parent_rank, parent_taxon_name, logger):
    """
    Find all child taxa under a given parent taxon.
    """
    child_rank = get_child_rank(parent_rank)
    if not child_rank:
        logger.warning(f"No child rank exists for {parent_rank}")
        return {}
    
    child_taxa = {}
    parent_rank_code = parent_rank[0]  # Get first letter (genus -> 'g')
    
    # Look through all sequences at the child rank
    if child_rank in unite_index:
        for child_taxon_name, child_sequences in unite_index[child_rank].items():
            # For each child taxon, check if any of its sequences also belong to our parent taxon
            if isinstance(child_sequences, list):
                # Multiple sequences case
                matching_sequences = []
                for seq_info in child_sequences:
                    if (parent_rank_code in seq_info['taxonomy_info']['taxonomy'] and 
                        seq_info['taxonomy_info']['taxonomy'][parent_rank_code] == parent_taxon_name):
                        matching_sequences.append(seq_info)
                if matching_sequences:
                    child_taxa[child_taxon_name] = matching_sequences
            else:
                # Single sequence case
                seq_info = child_sequences
                if (parent_rank_code in seq_info['taxonomy_info']['taxonomy'] and 
                    seq_info['taxonomy_info']['taxonomy'][parent_rank_code] == parent_taxon_name):
                    child_taxa[child_taxon_name] = [seq_info]
    
    logger.info(f"Found {len(child_taxa)} child taxa under {parent_rank} '{parent_taxon_name}': {list(child_taxa.keys())}")
    return child_taxa


def distribute_sequences_across_taxa(child_taxa, target_count, logger):
    """
    Distribute a target number of sequences across child taxa as evenly as possible.
    Redistributes unused allocations when taxa have insufficient sequences.
    """
    if not child_taxa:
        return []
    
    # First pass: calculate ideal distribution
    num_taxa = len(child_taxa)
    base_per_taxon = target_count // num_taxa
    extra_sequences = target_count % num_taxa
    
    logger.info(f"Distributing {target_count} sequences across {num_taxa} taxa: "
               f"{base_per_taxon} base sequences per taxon, {extra_sequences} taxa get +1 extra")
    
    # Create list of taxa with their allocations and available sequences
    taxa_info = []
    for i, (taxon_name, sequences) in enumerate(child_taxa.items()):
        allocation = base_per_taxon + (1 if i < extra_sequences else 0)
        available = len(sequences)
        taxa_info.append({
            'name': taxon_name,
            'sequences': sequences,
            'allocation': allocation,
            'available': available,
            'can_take': min(allocation, available)
        })
    
    # Second pass: redistribute unused allocations
    total_can_take = sum(info['can_take'] for info in taxa_info)
    if total_can_take < target_count:
        # We can't reach target even with perfect distribution
        logger.warning(f"Cannot reach target {target_count} sequences - only {total_can_take} available across all taxa")
    else:
        # Redistribute unused allocations to taxa that can take more
        unused = sum(info['allocation'] - info['can_take'] for info in taxa_info if info['can_take'] < info['allocation'])
        if unused > 0:
            logger.info(f"Redistributing {unused} unused allocations from taxa with insufficient sequences")
            
            # Find taxa that can take more than their current allocation
            expandable_taxa = [info for info in taxa_info if info['available'] > info['can_take']]
            if expandable_taxa:
                # Distribute unused sequences among expandable taxa
                while unused > 0 and expandable_taxa:
                    for info in expandable_taxa[:]:
                        if unused <= 0:
                            break
                        if info['can_take'] < info['available']:
                            info['can_take'] += 1
                            unused -= 1
                        if info['can_take'] >= info['available']:
                            expandable_taxa.remove(info)
    
    # Third pass: collect sequences
    selected_sequences = []
    for info in taxa_info:
        taxon_name = info['name']
        sequences = info['sequences']
        take_count = info['can_take']
        
        # Sort sequences by length (longest first) and take the requested number
        sorted_sequences = sorted(sequences, key=lambda x: x['length'], reverse=True)
        
        for seq_info in sorted_sequences[:take_count]:
            # Create a copy and add source information
            seq_copy = seq_info.copy()
            seq_copy['diversity_source_taxon'] = taxon_name
            # The diversity_source_rank should be the rank of the sequences we're collecting,
            # not the child rank of that rank (which would be one level down)
            seq_copy['diversity_source_rank'] = seq_info['rank']
            selected_sequences.append(seq_copy)
        
        logger.info(f"  {taxon_name}: took {take_count}/{len(sequences)} sequences")
    
    logger.info(f"Total sequences collected with diversity: {len(selected_sequences)}")
    return selected_sequences


def find_sequences_with_diversity(ncbi_lineage, unite_index, target_rank, target_count, logger, traverse_limit=None):
    """
    Find sequences using diversity strategy - distribute across child taxa.
    Prioritises target species matches first (up to 25% of target_count).
    When traverse_limit is set, collects cumulatively across ranks until target is reached or limit hit.
    """
    hierarchy = get_rank_hierarchy()
    
    # Validate target rank exists in lineage
    if target_rank not in ncbi_lineage:
        logger.error(f"Target rank '{target_rank}' not found in NCBI lineage")
        return None, None
    
    # Determine search ranks based on traverse_limit
    start_index = hierarchy.index(target_rank)
    if traverse_limit:
        try:
            limit_index = hierarchy.index(traverse_limit)
            if limit_index < start_index:
                logger.error(f"Traverse limit '{traverse_limit}' is more specific than start rank '{target_rank}'")
                return None, None
            search_ranks = hierarchy[start_index:limit_index + 1]
        except ValueError:
            logger.error(f"Invalid traverse limit rank: {traverse_limit}")
            return None, None
    else:
        search_ranks = [target_rank]
    
    logger.info(f"Starting diversity search at '{target_rank}' rank, target: {target_count} sequences")
    if traverse_limit:
        logger.info(f"Will traverse up hierarchy collecting cumulatively up to '{traverse_limit}': {' -> '.join(search_ranks)}")
    else:
        logger.info("No traverse mode - will only search at specified rank")
    
    # For cumulative collection when traversing
    all_collected_sequences = []
    collection_info = []
    target_species_sequences = []
    
    # First, try to find target species matches (up to 25% of target)
    target_species_name = ncbi_lineage.get('species')
    if target_species_name:
        max_species_sequences = max(1, int(target_count * 0.25))
        logger.info(f"Checking for target species '{target_species_name}' matches (max {max_species_sequences} sequences, 25% of {target_count})")
        
        species_matches = search_rank_for_taxon(unite_index, 'species', target_species_name)
        if species_matches:
            # Convert to list if needed
            if isinstance(species_matches, dict):
                species_seq_list = [species_matches]
            else:
                species_seq_list = species_matches
            
            # Sort by length and take up to 25%
            species_seq_list = sorted(species_seq_list, key=lambda x: x['length'], reverse=True)
            target_species_sequences = species_seq_list[:max_species_sequences]
            
            # Tag these sequences
            for seq in target_species_sequences:
                seq_copy = seq.copy()
                seq_copy['diversity_source_taxon'] = target_species_name
                seq_copy['diversity_source_rank'] = 'species'
                seq_copy['target_species_match'] = True
                all_collected_sequences.append(seq_copy)
            
            logger.info(f"Found target species match! Collected {len(target_species_sequences)} sequences from '{target_species_name}'")
        else:
            logger.info(f"No sequences found for target species '{target_species_name}'")
    
    # Calculate remaining sequences needed
    remaining_needed = target_count - len(all_collected_sequences)
    logger.info(f"Need {remaining_needed} more sequences after target species collection")
    
    if remaining_needed <= 0:
        logger.info(f"Target reached with species-level matches alone")
        diversity_info = {
            'strategy': 'diversity_species_only',
            'sequences_collected': len(all_collected_sequences),
            'target_species_sequences': len(target_species_sequences),
            'target_count': target_count
        }
        return all_collected_sequences, diversity_info
    
    # Now proceed with diversity distribution for remaining sequences
    for current_rank in search_ranks:
        if current_rank not in ncbi_lineage:
            continue
        
        # Calculate how many more sequences we need
        current_remaining = target_count - len(all_collected_sequences)
        if current_remaining <= 0:
            logger.info(f"Target {target_count} sequences already reached, stopping traverse")
            break
        
        current_taxon_name = ncbi_lineage[current_rank]
        logger.info(f"Searching for child taxa under {current_rank}: {current_taxon_name} "
                   f"(need {current_remaining} more sequences)")
        
        # Find child taxa under current rank
        child_taxa = find_child_taxa(unite_index, current_rank, current_taxon_name, logger)
        
        # Remove target species from child taxa if it exists (already processed)
        if target_species_name and target_species_name in child_taxa:
            logger.info(f"Removing target species '{target_species_name}' from child taxa (already processed)")
            del child_taxa[target_species_name]
        
        if not child_taxa:
            logger.warning(f"No child taxa found under {current_rank}: {current_taxon_name}")
            if not traverse_limit:
                break
            continue
        
        # Apply diversity distribution for this rank
        rank_sequences = distribute_sequences_across_taxa(child_taxa, current_remaining, logger)
        
        if rank_sequences:
            # Add to cumulative collection
            all_collected_sequences.extend(rank_sequences)
            collection_info.append({
                'rank': current_rank,
                'taxon_name': current_taxon_name,
                'child_rank': get_child_rank(current_rank),
                'num_child_taxa': len(child_taxa),
                'sequences_added': len(rank_sequences)
            })
            logger.info(f"Added {len(rank_sequences)} sequences from {current_rank} level "
                       f"(cumulative: {len(all_collected_sequences)}/{target_count})")
            
            # Check if we have enough sequences now
            if len(all_collected_sequences) >= target_count:
                logger.info(f"Target {target_count} sequences reached through cumulative collection")
                break
        else:
            logger.warning(f"No sequences collected at {current_rank} level")
            if not traverse_limit:
                break
    
    if all_collected_sequences:
        # Check if we reached target
        if len(all_collected_sequences) < target_count:
            if traverse_limit:
                logger.warning(f"Traverse limit '{traverse_limit}' reached but only collected "
                             f"{len(all_collected_sequences)}/{target_count} sequences")
            else:
                logger.warning(f"Could not reach target {target_count}, collected {len(all_collected_sequences)} sequences")
        
        # Return cumulative results
        diversity_info = {
            'strategy': 'diversity_traverse' if traverse_limit else 'diversity',
            'sequences_collected': len(all_collected_sequences),
            'target_species_sequences': len(target_species_sequences),
            'target_count': target_count,
            'collection_details': collection_info,
            'ranks_searched': [info['rank'] for info in collection_info],
            'traverse_limit': traverse_limit,
            'target_reached': len(all_collected_sequences) >= target_count
        }
        
        logger.info(f"Diversity {'traverse ' if traverse_limit else ''}search completed: {len(all_collected_sequences)} sequences")
        if target_species_sequences:
            logger.info(f"  Target species: {len(target_species_sequences)} sequences")
        for info in collection_info:
            logger.info(f"  {info['rank']}: {info['sequences_added']} sequences from "
                       f"{info['num_child_taxa']} {info['child_rank']} taxa")
        
        return all_collected_sequences, diversity_info
    
    logger.warning("Diversity search failed - no suitable child taxa found")
    return None, None


def find_best_unite_match(ncbi_lineage, unite_index, logger, return_all_sequences=False, max_sequences=None, target_rank=None, traverse_limit=None, diversity=False):
    """
    Find the best matching UNITE sequence(s) for an NCBI lineage with diversity support.
    """
    # Define valid ranks in hierarchical order
    valid_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    
    # Validate diversity usage
    if diversity:
        if not target_rank:
            logger.error("--diversity requires --tax_rank to be specified")
            return None, None
        if max_sequences is None:
            logger.error("--diversity requires --number to be specified")
            return None, None
        
        # Use diversity strategy
        return find_sequences_with_diversity(ncbi_lineage, unite_index, target_rank, max_sequences, logger, traverse_limit)
    
    # Original logic for non-diversity cases
    if target_rank is not None:
        if target_rank not in valid_ranks:
            logger.error(f"Invalid target rank '{target_rank}'. Valid ranks: {', '.join(valid_ranks)}")
            return None, None
        
        # Check if the requested rank exists in the NCBI lineage
        if target_rank not in ncbi_lineage:
            logger.error(f"Requested rank '{target_rank}' not found in NCBI lineage for this taxid")
            logger.info(f"Available ranks: {', '.join(ncbi_lineage.keys())}")
            return None, None
        
        if traverse_limit and max_sequences is not None:
            # Traverse up taxonomy to collect enough sequences
            return find_sequences_with_traverse(ncbi_lineage, unite_index, target_rank, max_sequences, logger, traverse_limit)
        else:
            # Original behavior: search only at specified rank
            rank_priority = [target_rank]
            logger.info(f"Searching only at rank '{target_rank}' (no fallback)")
    else:
        # Default behavior: search all ranks with fallback
        rank_priority = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
        logger.info("Using default search strategy with fallback (species -> genus -> family -> ...)")
    
    # Original search logic for non-traverse cases
    for rank in rank_priority:
        if rank in ncbi_lineage and rank in unite_index:
            taxon_name = ncbi_lineage[rank]
            sequences = search_rank_for_taxon(unite_index, rank, taxon_name)
            
            if sequences:
                return process_sequences_for_output(sequences, rank, taxon_name, return_all_sequences, max_sequences, logger)
    
    logger.warning(f"No match found for lineage: {', '.join([f'{k}:{v}' for k,v in ncbi_lineage.items()])}")
    return None, None


def find_sequences_with_traverse(ncbi_lineage, unite_index, start_rank, target_count, logger, traverse_limit):
    """
    Traverse up taxonomic hierarchy to collect target number of sequences.
    Stops at traverse_limit rank even if target not reached.
    """
    valid_ranks = ['species', 'genus', 'family', 'order', 'class', 'phylum', 'kingdom']
    
    # Validate traverse_limit
    start_index = valid_ranks.index(start_rank)
    try:
        limit_index = valid_ranks.index(traverse_limit)
    except ValueError:
        logger.error(f"Invalid traverse limit rank: {traverse_limit}")
        return None, None
    
    if limit_index < start_index:
        logger.error(f"Traverse limit '{traverse_limit}' is more specific than start rank '{start_rank}'")
        return None, None
    
    # Create rank list from start to limit
    traverse_ranks = valid_ranks[start_index:limit_index + 1]
    
    collected_sequences = []
    matched_ranks_info = []
    
    logger.info(f"Starting traverse search at '{start_rank}' rank, target: {target_count} sequences, limit: '{traverse_limit}'")
    
    for rank in traverse_ranks:
        if len(collected_sequences) >= target_count:
            break
            
        if rank not in ncbi_lineage or rank not in unite_index:
            continue
            
        taxon_name = ncbi_lineage[rank]
        sequences = search_rank_for_taxon(unite_index, rank, taxon_name)
        
        if sequences:
            # Convert single sequence to list for uniform processing
            if isinstance(sequences, dict):
                rank_sequences = [sequences]
            else:
                rank_sequences = sequences
            
            # Sort by length (longest first)
            rank_sequences = sorted(rank_sequences, key=lambda x: x['length'], reverse=True)
            
            # Calculate how many sequences we need from this rank
            remaining_needed = target_count - len(collected_sequences)
            sequences_to_take = min(len(rank_sequences), remaining_needed)
            
            # Add sequences from this rank and TAG them with source information
            taken_sequences = rank_sequences[:sequences_to_take]
            for seq in taken_sequences:
                # Create a copy to avoid modifying the original
                seq_copy = seq.copy()
                seq_copy['source_rank'] = rank
                seq_copy['source_taxon'] = taxon_name
                collected_sequences.append(seq_copy)
            
            # Record which rank these came from
            matched_ranks_info.append({
                'rank': rank,
                'taxon_name': taxon_name,
                'sequences_taken': len(taken_sequences),
                'total_available': len(rank_sequences)
            })
            
            logger.info(f"Found {len(rank_sequences)} sequences at {rank} level: {taxon_name}")
            logger.info(f"Taking {sequences_to_take} sequences (collected {len(collected_sequences)}/{target_count})")
    
    if not collected_sequences:
        logger.warning(f"No sequences found during traverse search from {start_rank} to {traverse_limit}")
        return None, None
    
    # Check if we reached target
    if len(collected_sequences) < target_count:
        logger.warning(f"Traverse limit '{traverse_limit}' reached but only collected {len(collected_sequences)}/{target_count} sequences")
    
    # Log traverse summary
    logger.info(f"Traverse search completed: collected {len(collected_sequences)} sequences from {len(matched_ranks_info)} taxonomic levels")
    for rank_info in matched_ranks_info:
        logger.info(f"  {rank_info['rank']}: {rank_info['taxon_name']} "
                   f"({rank_info['sequences_taken']}/{rank_info['total_available']} sequences)")
    
    # Create summary of matched ranks for return value
    matched_ranks_summary = " + ".join([f"{info['rank']}({info['sequences_taken']})" for info in matched_ranks_info])
    
    return collected_sequences, matched_ranks_summary


def process_sequences_for_output(sequences, rank, taxon_name, return_all_sequences, max_sequences, logger):
    """
    Process sequences according to output requirements (all, longest, or specific number).
    """
    # Convert single sequence to list for uniform processing
    if isinstance(sequences, dict):
        sequence_list = [sequences]
    else:
        sequence_list = sequences
    
    # Sort by length (longest first) to ensure we get the best sequences
    sequence_list = sorted(sequence_list, key=lambda x: x['length'], reverse=True)
    
    # Determine how many sequences to return
    if max_sequences is not None:
        # Specific number requested
        selected_sequences = sequence_list[:max_sequences]
        count_msg = f"{len(selected_sequences)} sequences (requested {max_sequences})"
    elif return_all_sequences:
        # All sequences requested
        selected_sequences = sequence_list
        count_msg = f"{len(selected_sequences)} sequences (all)"
    else:
        # Only longest sequence (default behavior)
        selected_sequences = [sequence_list[0]]
        count_msg = f"1 sequence (longest, length: {selected_sequences[0]['length']})"
    
    logger.info(f"Found match at {rank} level: {taxon_name} ({count_msg})")
    
    # Return in the expected format
    if not return_all_sequences and max_sequences is None:
        # Original behavior: return single sequence as dict
        return selected_sequences[0], rank
    else:
        # Return as list
        return selected_sequences, rank


def extract_unite_sequence_id(header):
    """
    Extract the UNITE sequence ID from the original header.
    Expected format: '>MK092089|k__Fungi;p__...|SH1320883.10FU'
    """
    # Remove the '>' character
    header = header.lstrip('>')
    
    # Extract the first part before the first '|'
    parts = header.split('|')
    if parts:
        return parts[0]
    else:
        return "unknown"


def write_output_fasta(sequence_info, matched_rank, sample_id, taxid, output_file, logger, summary_data=None):
    """
    Write the selected sequence(s) to output FASTA file with custom header.
    Handles both single sequences, traverse results, and diversity results.
    Aggregates data for summary CSV instead of writing directly.
    """
    # Create output directory if it doesn't exist
    os.makedirs(os.path.dirname(output_file), exist_ok=True)
    
    # Handle both single sequence and multiple sequences
    if isinstance(sequence_info, list):
        sequences = sequence_info
    else:
        sequences = [sequence_info]
    
    total_length = 0
    
    # Lists to collect data for summary
    reference_names = []
    taxa_list = []
    accessions = []
    lengths = []
    
    with open(output_file, 'w') as f:
        for seq_info in sequences:
            # Extract UNITE sequence ID from original header
            unite_seq_id = extract_unite_sequence_id(seq_info['header'])
            
            # Determine which rank and taxon to use in the header
            if 'target_species_match' in seq_info and seq_info['target_species_match']:
                # This is a target species match from diversity mode
                header_rank = seq_info['diversity_source_rank']
                header_taxon = seq_info['diversity_source_taxon']
                header_suffix = f"div-target-{header_rank}-{header_taxon}"
            elif 'diversity_source_taxon' in seq_info:
                # This sequence came from diversity mode
                header_rank = seq_info['diversity_source_rank']
                header_taxon = seq_info['diversity_source_taxon']
                header_suffix = f"div-{header_rank}-{header_taxon}"
            elif 'source_rank' in seq_info:
                # This sequence came from traverse mode
                header_rank = seq_info['source_rank']
                header_taxon = seq_info['source_taxon']
                header_suffix = f"trav-{header_rank}-{header_taxon}"
            else:
                # This sequence came from regular mode
                if isinstance(matched_rank, dict) and 'search_rank' in matched_rank:
                    # Diversity mode matched_rank info
                    header_rank = matched_rank['child_rank']
                    header_taxon = seq_info.get('taxon_name', 'unknown')
                    header_suffix = f"div-{header_rank}-{header_taxon}"
                else:
                    # Regular mode
                    header_rank = matched_rank if isinstance(matched_rank, str) else 'unknown'
                    header_taxon = seq_info.get('taxon_name', 'unknown')
                    header_suffix = f"{header_rank}-{header_taxon}"
            
            # Create new header: {sample_id}-{taxid}-{strategy_info}-{unite_seq_id}
            new_header = f"{sample_id}-{taxid}-{header_suffix}-{unite_seq_id}"
            # Replace spaces with underscores to ensure BWA compatibility
            new_header = new_header.replace(' ', '_')

            header_line = f">{new_header}\n"
            f.write(header_line)
            
            # Write sequence with line breaks every 80 characters
            sequence = seq_info['sequence']
            for i in range(0, len(sequence), 80):
                seq_line = sequence[i:i+80] + '\n'
                f.write(seq_line)
            
            total_length += len(sequence)
            
            # Collect data for summary CSV
            if summary_data is not None:
                reference_names.append(seq_info['header'])
                taxa_list.append(header_taxon)
                accessions.append(unite_seq_id)
                lengths.append(seq_info['length'])
    
    logger.info(f"Wrote {len(sequences)} sequence(s) to {output_file} (total length: {total_length} bp)")
    
    # Add aggregated data to summary_data dictionary
    if summary_data is not None and sample_id not in summary_data:
        summary_data[sample_id] = {
            'reference_names': reference_names,
            'taxa': taxa_list,
            'accessions': accessions,
            'lengths': lengths,
            'taxid': taxid
        }
        logger.info(f"Added aggregated data to summary for sample {sample_id}")


def normalise_csv_headers(fieldnames, logger, id_header=None):
    """
    Create a mapping from normalised headers to actual headers in the CSV.
    Returns a dictionary mapping expected column names to actual column names.
    """
    
    # Define acceptable variations for each required column
    id_variations = ['id', 'ID', 'Id', 'sample_id', 'sample_ID', 'Sample_ID', 'SampleID', 'sample', 'Sample']
    taxid_variations = ['taxid', 'taxID', 'TaxID', 'TAXID', 'tax_id', 'Tax_ID', 'TAX_ID', 'taxonomy_id', 'ncbi_taxid', 'NCBI_TaxID']
    
    # Find ID column - prioritise user-specified column
    id_column = None
    if id_header is not None:
        if id_header in fieldnames:
            id_column = id_header
        else:
            logger.warning(f"Specified ID column '{id_header}' not found in CSV. Will try standard variations.")
    
    # If no user-specified column or it wasn't found, try standard variations
    if id_column is None:
        for variation in id_variations:
            if variation in fieldnames:
                id_column = variation
                break
    
    # Find taxid column
    taxid_column = None
    for variation in taxid_variations:
        if variation in fieldnames:
            taxid_column = variation
            break
    
    return {
        'ID': id_column,
        'taxid': taxid_column
    }


def validate_input_files(csv_file, fasta_file, logger, id_header=None):
    """
    Validate that input files exist and have required format.
    """
    # Check CSV file
    if not os.path.exists(csv_file):
        logger.error(f"Input CSV file does not exist: {csv_file}")
        sys.exit(1)
    
    try:
        with open(csv_file, 'r') as f:
            reader = csv.DictReader(f)
            header_mapping = normalise_csv_headers(reader.fieldnames,logger, id_header, )
            
            if header_mapping['ID'] is None:
                logger.error(f"CSV file must contain an ID column. Found columns: {list(reader.fieldnames)}")
                if id_header:
                    logger.error(f"Specified ID column '{id_header}' was not found in the CSV.")
                logger.error(f"Acceptable ID column names: id, ID, Id, sample_id, sample_ID, Sample_ID, SampleID, sample, Sample")
                sys.exit(1)
            
            if header_mapping['taxid'] is None:
                logger.error(f"CSV file must contain a taxid column. Found columns: {list(reader.fieldnames)}")
                logger.error(f"Acceptable taxid column names: taxid, taxID, TaxID, TAXID, tax_id, Tax_ID, TAX_ID, taxonomy_id, ncbi_taxid, NCBI_TaxID")
                sys.exit(1)
            
            logger.info(f"Found ID column: '{header_mapping['ID']}', taxid column: '{header_mapping['taxid']}'")
            
    except Exception as e:
        logger.error(f"Error reading CSV file: {str(e)}")
        sys.exit(1)
    
    # Check FASTA file
    if not os.path.exists(fasta_file):
        logger.error(f"UNITE FASTA file does not exist: {fasta_file}")
        sys.exit(1)
    
    logger.info("Input file validation passed")
    return header_mapping


def write_summary_csv(summary_data, output_path, logger):
    """
    Write the aggregated summary data to CSV file.
    Each sample appears once with semicolon-delimited lists and length statistics.
    """
    fieldnames = ['ID', 'searched_taxid', 'searched_taxa', 'number', 'reference_name', 
                  'taxa', 'accession', 'min_length', 'max_length', 'avg_length']
    
    try:
        with open(output_path, 'w', newline='') as csvfile:
            writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
            writer.writeheader()
            
            for sample_id, data in summary_data.items():
                # Calculate length statistics
                lengths = data['lengths']
                min_length = min(lengths) if lengths else 0
                max_length = max(lengths) if lengths else 0
                avg_length = sum(lengths) / len(lengths) if lengths else 0
                
                # Create semicolon-delimited strings
                reference_name_str = ';'.join(data['reference_names'])
                taxa_str = ';'.join(data['taxa'])
                accession_str = ';'.join(data['accessions'])
                
                row = {
                    'ID': sample_id,
                    'searched_taxid': data['taxid'],
                    'searched_taxa': data.get('searched_taxa', ''),
                    'number': len(data['reference_names']),
                    'reference_name': reference_name_str,
                    'taxa': taxa_str,
                    'accession': accession_str,
                    'min_length': min_length,
                    'max_length': max_length,
                    'avg_length': round(avg_length, 2)
                }
                
                writer.writerow(row)
        
        logger.info(f"Successfully wrote summary CSV to {output_path}")
        logger.info(f"Summary contains {len(summary_data)} samples")
        
    except Exception as e:
        logger.error(f"Failed to write summary CSV: {str(e)}")
        raise


def main():
    """Main function to orchestrate the UNITE sequence extraction process."""

    # Set up argument parser
    parser = setup_argument_parser()
    args = parser.parse_args()

    # Set up logger
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'UNITEd_log_{timestamp}.log'
    
    # Set up logging    
    log_dir_path = os.path.join(args.output, "logs")
    logger = setup_logging(log_dir=log_dir_path, log_file=args.log_file)
  
    # Validate argument combinations
    if args.all and args.number is not None:
        logger.warning("Both --all and --number specified. --number takes precedence.")
        args.all = False
    
    # Validate diversity argument usage
    if args.diversity:
        if not args.tax_rank:
            logger.error("--diversity requires --tax_rank to be specified")
            sys.exit(1)
        if args.number is None:
            logger.error("--diversity requires --number to be specified")
            sys.exit(1)
        if args.all:
            logger.warning("--diversity overrides --all flag")
            args.all = False
        logger.info(f"Diversity mode enabled: will distribute {args.number} sequences across child taxa of {args.tax_rank}")
        logger.info(f"  -> Will prioritise target species matches first (up to 25% of {args.number} = {max(1, int(args.number * 0.25))} sequences)")
        if args.traverse:
            logger.info(f"  -> Traverse mode also enabled: will traverse up to '{args.traverse}' if insufficient child taxa found")
    
    # Validate traverse argument usage
    if args.traverse and not args.diversity:
        if not args.tax_rank:
            logger.error("--traverse requires --tax_rank to be specified")
            sys.exit(1)
        if args.number is None:
            logger.error("--traverse requires --number to be specified (what target number of sequences?)")
            sys.exit(1)
        logger.info(f"Traverse mode enabled: will search from {args.tax_rank} up to {args.traverse} to collect {args.number} sequences")
    
    logger.info("Starting UNITE sequence extraction")
    logger.info(f"Input CSV: {args.tracking_sheet}")
    logger.info(f"UNITE FASTA: {args.unite_db}")
    logger.info(f"Output directory: {args.output}")
    logger.info(f"Email: {args.email}")
    
    if args.number is not None:
        logger.info(f"Retrieve specific number of sequences: {args.number}")
        keep_all_sequences = True  # Need to keep all sequences to select top N or apply diversity
    elif args.all:
        logger.info("Retrieve all sequences: True")
        keep_all_sequences = True
    else:
        logger.info("Retrieve all sequences: False (longest only)")
        keep_all_sequences = False
    
    if args.tax_rank:
        strategy_msg = f"Searching at taxonomic rank: {args.tax_rank}"
        if args.diversity:
            strategy_msg += " (with diversity across child taxa, prioritising target species)"
        if args.traverse:
            strategy_msg += f" (with traverse up to {args.traverse} if needed)" if args.diversity else f" (with traverse up to {args.traverse} if needed)"
        if not args.diversity and not args.traverse:
            strategy_msg += " (no fallback)"
        logger.info(strategy_msg)
    else:
        logger.info("Using default search strategy with fallback (species -> genus -> family -> ...)")
    
    if keep_all_sequences:
        logger.warning("WARNING: Sequence collection mode enabled. This will use more memory for large UNITE databases.")
    
    # Validate input files and get header mapping
    header_mapping = validate_input_files(args.tracking_sheet, args.unite_db, logger, args.id_header)
    
    # Parse UNITE FASTA and build taxonomy index
    try:
        unite_index = build_unite_taxonomy_index(args.unite_db, logger, keep_all_sequences=keep_all_sequences)
        if not unite_index:
            logger.error("No taxonomy information found in UNITE FASTA file")
            sys.exit(1)
    except Exception as e:
        logger.error(f"Failed to parse UNITE FASTA file: {str(e)}")
        sys.exit(1)
    
    # Process each sample in the CSV
    logger.info(f"Reading input CSV: {args.tracking_sheet}")
    
    total_samples = 0
    successful_matches = 0
    failed_matches = 0
    
    # Initialise summary CSV writer if requested
    summary_data = {}  # Dictionary to aggregate data by sample ID
    if args.summary_csv:
        try:
            # Place summary CSV in the output directory
            summary_csv_path = os.path.join(args.output, args.summary_csv)
            
            # Create output directory if it doesn't exist
            os.makedirs(args.output, exist_ok=True)
            
            # We'll write the summary at the end after aggregating all data
            logger.info(f"Will create summary CSV: {summary_csv_path}")
        except Exception as e:
            logger.error(f"Failed to set up summary CSV file: {str(e)}")
            sys.exit(1)
    
    try:
        with open(args.tracking_sheet, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            
            for row in reader:
                total_samples += 1
                sample_id = row[header_mapping['ID']].strip()
                taxid = row[header_mapping['taxid']].strip()
                
                logger.info(f"Processing sample {sample_id} with taxid {taxid}")
                
                try:
                    # Get NCBI lineage
                    ncbi_lineage = get_ncbi_lineage(taxid, args.email, logger, args.api)
                    
                    # Get most specific taxon name for summary
                    searched_taxa = get_most_specific_taxon(ncbi_lineage)
                    
                    # Find best UNITE match (now with diversity support)
                    sequences, matched_rank = find_best_unite_match(
                        ncbi_lineage, 
                        unite_index, 
                        logger,
                        return_all_sequences=args.all, 
                        max_sequences=args.number,
                        target_rank=args.tax_rank,
                        traverse_limit=args.traverse,
                        diversity=args.diversity
                    )
                    
                    # Define output file path - CHANGED: flat structure with sample ID
                    output_file = os.path.join(args.output, f"{sample_id}_seed.fasta")
                    
                    if sequences is None:
                        logger.warning(f"No UNITE sequence found for sample {sample_id} (taxid: {taxid})")
                        # Create empty output file with comment
                        os.makedirs(os.path.dirname(output_file), exist_ok=True)
                        with open(output_file, 'w') as f:
                            f.write(f"# No UNITE sequence found for taxid {taxid}\n")
                            f.write(f"# NCBI lineage: {', '.join([f'{k}:{v}' for k,v in ncbi_lineage.items()])}\n")
                        
                        failed_matches += 1
                    else:
                        # Write successful match(es) and collect summary data
                        write_output_fasta(sequences, matched_rank, sample_id, taxid, output_file, logger, 
                                         summary_data if args.summary_csv else None)
                        
                        # Add searched_taxa to summary data
                        if args.summary_csv and sample_id in summary_data:
                            summary_data[sample_id]['searched_taxa'] = searched_taxa
                        
                        logger.info(f"Successfully wrote seed.fasta for sample {sample_id} (taxid: {taxid})")
                        
                        # Log diversity/traverse strategy results
                        if isinstance(matched_rank, dict) and 'strategy' in matched_rank:
                            if matched_rank.get('target_species_sequences', 0) > 0:
                                logger.info(f"  -> Target species priority: {matched_rank['target_species_sequences']} sequences")
                            if not matched_rank.get('target_reached', True):
                                logger.warning(f"  -> Target not fully reached: {matched_rank['sequences_collected']}/{matched_rank['target_count']}")
                        
                        successful_matches += 1
                    
                except Exception as e:
                    logger.error(f"Error processing sample {sample_id} (taxid: {taxid}): {str(e)}")
                    # Create empty output file for failed cases
                    output_file = os.path.join(args.output, f"{sample_id}_seed.fasta")
                    os.makedirs(os.path.dirname(output_file), exist_ok=True)
                    with open(output_file, 'w') as f:
                        f.write(f"# Error processing taxid {taxid}: {str(e)}\n")
                    
                    failed_matches += 1

    except Exception as e:
        logger.error(f"Error reading CSV file: {str(e)}")
        sys.exit(1)
    
    # Write summary CSV if requested
    if args.summary_csv and summary_data:
        try:
            write_summary_csv(summary_data, summary_csv_path, logger)
        except Exception as e:
            logger.error(f"Failed to write summary CSV: {str(e)}")
    
    # Log final statistics
    logger.info("="*50)
    logger.info("UNITE sequence extraction completed")
    logger.info(f"Total samples processed: {total_samples}")
    logger.info(f"Successful matches: {successful_matches}")
    logger.info(f"Failed matches: {failed_matches}")
    logger.info(f"Success rate: {(successful_matches/total_samples)*100:.1f}%")
    logger.info("="*50)


if __name__ == "__main__":
    main()
