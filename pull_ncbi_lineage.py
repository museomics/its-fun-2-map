import os
import pathlib
import pandas as pd
import urllib.error
import time
import random
import logging
import argparse
from Bio import Entrez
from datetime import datetime

from its_fun_tools import get_ncbi_lineage, log_and_print
from metahist_tools import xlsx2csv, setup_logging


# Increase time between and number of tries used by entrez (from go_fetch.py)
Entrez.sleep_between_tries = 20
Entrez.max_tries = 20


def add_ncbi_lineages_to_csv(input_csv, output_csv, taxcolumn, email, logger, api_key=None, sheet=1):
    """
    Add NCBI taxonomic lineages to a CSV file based on taxids.

    Parameters:
    input_csv (str): Path to the input CSV file. XLSX is also an accepted format.
    output_csv (str): Path to the output CSV file.
    taxcolumn (str): Column name in CSV that contains taxids.
    email (str): Email address for NCBI Entrez.
    logger (logging.Logger): Logger for logging messages.
    api_key (str, optional): NCBI API key for increased rate limits.
    sheet (int, optional): Sheet number to read if input is XLSX. Default is 1.
    """

    Entrez.email = email
    if api_key and api_key != "None":
        Entrez.api_key = api_key

    filetype = pathlib.Path(input_csv).suffix.lower()
    if filetype == ".xlsx":
        df = xlsx2csv(input_csv, sheet=sheet)
    elif filetype == ".csv":
        df = pd.read_csv(input_csv)
    else:
        raise ValueError(f"Unsupported tracking sheet format: {tracking_sheet}")

    if taxcolumn not in df.columns:
        logger.error(f"Taxid column '{taxcolumn}' not found in input CSV.")
        return
        # Prepare columns for lineage data
    ranks = ['kingdom', 'phylum', 'class', 'order', 'family', 'genus', 'species']
    for rank in ranks:
        df[rank] = None

    # Process each row
    for index, row in df.iterrows():
        taxid = str(row[taxcolumn])
        try:
            lineage = get_ncbi_lineage(taxid, Entrez.email, logger, Entrez.api_key)
            for rank in ranks:
                if rank in lineage:
                    df.at[index, rank] = lineage[rank]
        except Exception as e:
            logger.error(f"Could not retrieve lineage for taxid {taxid}: {str(e)}")

    # Write output CSV
    df.to_csv(output_csv, index=False)
    logger.info(f"Lineage data added and saved to {output_csv}")

def main():

    # Set up logging
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'ncbi_pull_lineage_{timestamp}.log'

    logger = setup_logging(log_dir="./", log_file=args.log_file)

    # Add lineages to CSV
    add_ncbi_lineages_to_csv(args.input_csv, args.output_csv, args.taxcolumn, args.email, logger, args.api_key, args.sheet)

if __name__ == "__main__":

    # Set up argument parser
    parser = argparse.ArgumentParser(description="Add NCBI taxonomic lineages to a CSV file based on taxids.")
    parser.add_argument("--input_csv", help="Path to the input CSV file.")
    parser.add_argument("--output_csv", help="Path to the output CSV file.")
    parser.add_argument("--email", help="Email address for NCBI Entrez.")
    parser.add_argument("--api_key", help="NCBI API key for increased rate limits.", default=None)
    parser.add_argument("--taxcolumn", help="Column name in CSV that contains taxids.", default="taxid")
    parser.add_argument("--log_file", help="Path to log file.", default="ncbi_lineage.log")
    parser.add_argument("--sheet", help="Sheet number to read if input is XLSX. Default is 1.", type=int, default=1)    
    args = parser.parse_args()

    main()
