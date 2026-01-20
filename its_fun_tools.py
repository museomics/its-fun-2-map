import subprocess
import glob
import os
import argparse
import sys
import pandas as pd
import pathlib
import shutil
import logging
from Bio import Entrez

from metahist_tools import xlsx2csv

#### Helper functions used across scripts

def blast_task(scaff_path, database_file, output_dir, name_id, prefix=None):
    parent_dir = os.path.basename(os.path.dirname(scaff_path))
    # Use name_id as prefix if no prefix provided
    file_prefix = f"{prefix}_" if prefix else ""
    output_file = f"{file_prefix}{name_id}_blast_results.tsv"
    output_path = os.path.join(output_dir, output_file)

    print(f"[INFO] Running BLAST: query={scaff_path} -> {output_path}")
    subprocess.run([
        "blastn",
        "-query", scaff_path,
        "-db", database_file,
        "-outfmt", "6 qseqid sseqid pident length mismatch gapopen qstart qend sstart send evalue bitscore",
        "-evalue", "1e-10",
        "-out", output_path
    ], check=True)




def log_and_print(message, level='info'):
    """
    Log a message and print it to console.
    """
    logger = logging.getLogger()
    if level == 'info':
        logger.info(message)
    elif level == 'error':
        logger.error(message)
    elif level == 'warning':
        logger.warning(message)


def cleanup_temp_dir(temp_dir):
    if os.path.exists(temp_dir):
        try:
            shutil.rmtree(temp_dir)
        except Exception as e:
            raise ValueError(f"Failed to clean up temporary directory {temp_dir}: {e}")


def load_name_ids(tracking_sheet, column_name, sheet=None):
    filetype = pathlib.Path(tracking_sheet).suffix.lower()
    if filetype == ".xlsx":
        df = xlsx2csv(tracking_sheet, sheet=sheet)
    elif filetype == ".csv":
        df = pd.read_csv(tracking_sheet)
    else:
        raise ValueError(f"Unsupported tracking sheet format: {tracking_sheet}")
    name_ids = df[column_name].dropna().astype(str).tolist()
    return name_ids


def get_ncbi_lineage(taxid, email, logger, api_key=None):
    """
    Get NCBI taxonomic lineage for a given taxid with retry logic.
    Parameters:
    taxid (str): NCBI taxonomy ID.
    email (str): Email address for NCBI Entrez.
    logger (logging.Logger): Logger for logging messages.
    api_key (str, optional): NCBI API key for increased rate limits.

    Returns:
    dict: A dictionary with taxonomic ranks as keys and names as values.
    
    """
    Entrez.email = email
    if api_key and api_key != "None":
        Entrez.api_key = api_key
    
    log_and_print(f"Getting NCBI lineage for taxid: {taxid}")
    
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
            log_and_print(f"Attempt {attempt+1}/{max_retries}: Error getting lineage for {taxid}: {str(e)}", level="warning")
            
            if attempt < max_retries - 1:
                logger.info(f"Retrying in {delay:.2f} seconds...")
                time.sleep(delay)
            else:
                log_and_print(f"Failed to get lineage for {taxid} after {max_retries} attempts.", level="error")
                raise
    
    raise Exception(f"Failed to get lineage for {taxid}")

