import subprocess
import glob
import os
import argparse
import sys
import pandas as pd
import pathlib
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from its_fun_tools import load_name_ids, blast_task
from metahist_tools import xlsx2csv, setup_logging

### Author: Maria Kamouyiaros & Dan Parsons (NHMUK)

def run_blast_pipeline(database_file, makeblastdb, query_dir, output_dir, prefix, name_ids, max_workers=4, log_file=None):
    os.makedirs(output_dir, exist_ok=True)

    if log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        log_file = f'Blast1_log_{timestamp}.log'

    logger = setup_logging(log_file=log_file)

    if makeblastdb:
        logger.info(f"[INFO] Creating BLAST DB for {database_file}")
        subprocess.run(["makeblastdb", "-in", database_file, "-dbtype", "nucl"], check=True)

    # Find only scaffolds.fasta in the immediate *.spades.out directories
    # NOTE: We intentionally do NOT fall back to contigs.fasta, as failed assemblies
    # may have contigs.fasta but no valid scaffolds.fasta
    scaffold_paths = glob.glob(os.path.join(query_dir, "*.spades.out", "scaffolds.fasta"))

    tasks = []
    skipped_samples = []

    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        for name_id in name_ids:
            # Only look for scaffolds.fasta
            matched_scaffolds = [p for p in scaffold_paths if name_id in p]

            if matched_scaffolds:
                logger.info(f"Found scaffolds.fasta for {name_id}")
                for scaff_path in matched_scaffolds:
                    tasks.append(executor.submit(blast_task, scaff_path, database_file, output_dir, name_id, prefix))
            else:
                logger.warning(f"[SKIP] No scaffolds.fasta found for {name_id} - assembly may have failed")
                skipped_samples.append(name_id)

    # Log summary of skipped samples
    if skipped_samples:
        logger.info(f"\n[SUMMARY] Skipped {len(skipped_samples)} samples with no scaffolds.fasta:")
        for sample in skipped_samples:
            logger.info(f"  - {sample}")

    logger.info(f"[DONE] All tasks complete. Processed {len(tasks)} samples, skipped {len(skipped_samples)} samples.")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run BLAST for all scaffolds.fasta matches per genome ID.")
    parser.add_argument("--database_file", required=True, help="Database file")
    parser.add_argument("--makeblastdb", action="store_true", help="Create BLAST DB from --database_file")
    parser.add_argument("--query_dir", required=True, help="Directory containing query FASTA files")
    parser.add_argument("--output_dir", required=True, help="Directory to store BLAST output files")
    parser.add_argument("--tracking_sheet", required=True, help="Tracking spreadsheet with relevant metadata")
    parser.add_argument("--sheet", type=int, required=False, help="(For XLSX input only) Sheet index (0-based)")
    parser.add_argument("--column_name", required=True, help="Column name in tracking sheet with library names")
    parser.add_argument("--prefix", required=False, help="(Optional) Prefix to add to output files")
    parser.add_argument("--max_workers", type=int, default=4, help="Maximum parallel workers")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `blast1_log` file will be created")

    args = parser.parse_args()

    name_ids = load_name_ids(args.tracking_sheet, args.column_name, sheet=args.sheet)
    if not name_ids:
        raise ValueError(f"No valid mappings found in {args.tracking_sheet}")

    run_blast_pipeline(
        database_file=args.database_file,
        makeblastdb=args.makeblastdb,
        query_dir=args.query_dir,
        output_dir=args.output_dir,
        prefix=args.prefix,
        name_ids=name_ids,
        max_workers=args.max_workers,
        log_file=args.log_file
    )
