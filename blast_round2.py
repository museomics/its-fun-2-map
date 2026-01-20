import os
import glob
import subprocess
import shutil
import pathlib
from pathlib import Path
from concurrent.futures import ThreadPoolExecutor, as_completed
import pandas as pd
import argparse
import sys
from datetime import datetime

from its_fun_tools import load_name_ids, cleanup_temp_dir, blast_task
from metahist_tools import xlsx2csv, setup_logging


### Author: Maria Kamouyiaros
### 2025-08-28 
### Loads sample names from spreadsheet.
### Skips hard-coded scaffolds (misc/broken_scaffolds.fasta, K55/scaffolds.fasta).
### Matches scaffolds to TSV, extracts unique headers and creates a FASTA in tmp dir.
### Runs blast_task for each name_id on a specified database.


def find_tsv_for_id(blast_dir, sample_id):
    # Match any file that starts with the ID and ends with .tsv
    pattern = os.path.join(blast_dir, f"{sample_id}*.tsv")
    matches = glob.glob(pattern)
    return matches  # list of matching TSVs


def run_seqkit_and_blast(tracking_sheet, column_name,query_dir, blast_dir, database_file, output_dir, log_file="blast_round2.log", prefix=None, max_workers=4, makeblastdb=False):

    os.makedirs(output_dir, exist_ok=True)
    tmp_fasta_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_fasta_dir, exist_ok=True)

    # Set up log file with timestamp if not provided
    if log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = f'blast2_log_{timestamp}.log'
    
    logger = setup_logging(log_file=log_file)

    if makeblastdb:
        logger.info(f"[INFO] Creating BLAST DB for {database_file}")
        subprocess.run(["makeblastdb", "-in", database_file, "-dbtype", "nucl"], check=True)

    # Find only scaffolds.fasta in the immediate *.spades.out directories
    scaffold_paths = glob.glob(os.path.join(query_dir, "*.spades.out", "scaffolds.fasta"))
    name_ids = load_name_ids(tracking_sheet, column_name)


    tasks = []
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for name_id in name_ids:
            # Find scaffolds for this sample
            matched_scaffolds = [p for p in scaffold_paths if name_id in p]
            if not matched_scaffolds:
                logger.info(f"[SKIP] No scaffolds.fasta found for {name_id}")
                continue

            # Find TSVs for this sample
            tsv_files = find_tsv_for_id(blast_dir, name_id)
            if not tsv_files:
                logger.info(f"[SKIP] No TSV found for {name_id}")
                continue

            tsv_file = tsv_files[0]
            logger.info(f"[FOUND] TSV {tsv_file} for {name_id}")

            # Extract unique headers from TSV
            unique_ids_file = os.path.join(tmp_fasta_dir, f"{name_id}_unique_headers.txt")
            with open(tsv_file) as f, open(unique_ids_file, "w") as out_f:
                ids = {line.split("\t")[0] for line in f if line.strip()}
                out_f.write("\n".join(sorted(ids)))

            # Run extraction + BLAST for each scaffold
            for scaff_path in matched_scaffolds:
                parent_dir = Path(scaff_path).parent.name

                output_fasta = os.path.join(
                    tmp_fasta_dir, f"{name_id}_matched_scaffolds.fasta"
                )
                cmd = [
                    "seqkit",
                    "grep",
                    "-f",
                    unique_ids_file,
                    scaff_path,
                    "-o",
                    output_fasta,
                ]
                logger.info(f"[INFO] Running seqkit: {' '.join(cmd)}")
                subprocess.run(cmd, check=True)

                # Submit BLAST task
                tasks.append(
                    executor.submit(
                        blast_task,
                        output_fasta,
                        database_file,
                        output_dir,
                        name_id,
                        prefix,
                    )
                )

        for future in as_completed(tasks):
            try:
                future.result()
            except Exception as e:
                logger.info(f"[ERROR] Task failed: {e}")

    # Clean up tmp dir
    cleanup_temp_dir(tmp_fasta_dir)
    logger.info("[DONE] All tasks complete.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run seqkit+BLAST pipeline on scaffolds.")
    parser.add_argument("--tracking_sheet", required=True, help="CSV/XLSX tracking sheet")
    parser.add_argument("--column_name", required=True, help="Column name in tracking sheet with sample IDs")
    parser.add_argument("--sheet", type=int, required=False, help="(For XLSX input only) Sheet index (0-based)")
    parser.add_argument("--query_dir", required=True, help="Directory containing scaffolds.fasta files")
    parser.add_argument("--blast_dir", required=True, help="Directory containing BLAST TSV results")
    parser.add_argument("--database_file", required=True, help="Database file for BLAST")
    parser.add_argument("--makeblastdb", action="store_true", help="Create BLAST DB from --database_file")
    parser.add_argument("--output_dir", default="blast_results", help="Directory to store BLAST output files")
    parser.add_argument("--prefix", required=False, help="Optional prefix to add to output files")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum parallel workers")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `blast2_log` file will be created")

    args = parser.parse_args()

    run_seqkit_and_blast(
        tracking_sheet=args.tracking_sheet,
        column_name=args.column_name,
        query_dir=args.query_dir,
        blast_dir=args.blast_dir,
        database_file=args.database_file,
        makeblastdb=args.makeblastdb,
        output_dir=args.output_dir,
        prefix=args.prefix,
        max_workers=args.max_workers,
        log_file=args.log_file
    )
