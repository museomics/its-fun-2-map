import os
import glob
import subprocess
import argparse
import shutil
import csv
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from metahist_tools import setup_logging

#### Author: Maria Kamouyiaros & Daniel Parsons @ NHM
#### Date: 2025-08-04
#### Assumes that naming is consistent so it is hard-coded into searching for files with the terms merged and unmerged. 
#### Update 2025-08-05 - added in -k 21,33,55 because of getting invalid k-mer error for 1,000,000 read ts
#### Update 2025-08-27 - added in concurrent.futures parallelisation 
#### Update 2025-10-10 - Edited version of assembly_wrapper.py
####                   - Cleaned up and added in consistent logging (to match mapping and blast_parser) 
####                   - Added in assembly check and re-run with other k-mers if needed
#### Update 2025-10-31 - fixed path errors 
####		       - fixed fall-back logic to remove directory of failed assembly
#### Update 2025-11-14 - fixed return values in check functions
####                   - improved logging to show file existence, size, and sequence counts
####                   - fixed check functions to properly validate assembly output
#### Update 2025-11-17 - improved logging for better progress tracking per sample
####                   - fixed fallback 1 to use single-reads k=21 instead of paired-end
####                   - added clear stage labels and summary reporting
####                   - added optional CSV summary output with assembly metrics
#### Update 2025-11-18 - added N50 calculation to assembly metrics
#### VERSION: 7.1.0

def get_scaffold_metrics(scaffolds_file):
    """
    Parse a scaffolds.fasta file and return metrics including N50.
    
    Returns:
        tuple: (n_scaffolds, scaffold_lengths_list, mean_length, n50)
               Returns (0, [], 0, 0) if file doesn't exist or is invalid
    """
    if not os.path.isfile(scaffolds_file) or os.path.getsize(scaffolds_file) == 0:
        return 0, [], 0, 0
    
    scaffold_lengths = []
    current_length = 0
    
    try:
        with open(scaffolds_file, 'r') as f:
            for line in f:
                line = line.strip()
                if line.startswith('>'):
                    # New sequence - save previous length if any
                    if current_length > 0:
                        scaffold_lengths.append(current_length)
                    current_length = 0
                else:
                    # Accumulate sequence length
                    current_length += len(line)
            
            # Don't forget the last sequence
            if current_length > 0:
                scaffold_lengths.append(current_length)
        
        n_scaffolds = len(scaffold_lengths)
        mean_length = sum(scaffold_lengths) / n_scaffolds if n_scaffolds > 0 else 0
        
        # Calculate N50
        if n_scaffolds > 0:
            # Sort scaffolds by length (descending)
            sorted_lengths = sorted(scaffold_lengths, reverse=True)
            total_length = sum(sorted_lengths)
            target_length = total_length / 2
            
            cumulative_length = 0
            n50 = 0
            for length in sorted_lengths:
                cumulative_length += length
                if cumulative_length >= target_length:
                    n50 = length
                    break
        else:
            n50 = 0
        
        return n_scaffolds, scaffold_lengths, mean_length, n50
    
    except Exception as e:
        return 0, [], 0, 0


def run_spades(merged_path, unmerged_1, unmerged_2, output_path, logger, k="21,33,55"):
    """Run one SPAdes job with paired-end mode."""
    cmd = [
        "spades.py",
        "--merged", merged_path,
        "-1", unmerged_1,
        "-2", unmerged_2,
        "-k", k,
        "-o", output_path
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return output_path


def run_single_spades(mapped_path, output_path, logger, k="21,33,55"):
    """Run one SPAdes job with single reads."""
    cmd = [
        "spades.py",
        "-s", mapped_path,
        "-k", k,
        "-o", output_path
    ]
    logger.info(f"Running: {' '.join(cmd)}")
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return output_path


def check_spades_out1(merged_path, unmerged_1, unmerged_2, output_path, logger, k):
    """Check if SPAdes output is valid; if not, rerun with single-reads k=21."""
    
    scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
    
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    logger.info(f"[{sample_name}] Checking assembly output")
    
    # Check scaffolds file
    if os.path.isfile(scaffolds_file):
        file_size = os.path.getsize(scaffolds_file)
        
        if file_size > 0:
            # Count number of sequences and get metrics
            with open(scaffolds_file, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
            
            if seq_count > 0:
                n_scaffolds, scaffold_lengths, mean_length, n50 = get_scaffold_metrics(scaffolds_file)
                logger.info(f"[{sample_name}] ✓ Initial assembly successful - {seq_count} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp ({file_size} bytes)")
                return True
            else:
                logger.warning(f"[{sample_name}] Scaffolds file is empty (no sequences found)")
        else:
            logger.warning(f"[{sample_name}] Scaffolds file exists but is empty (0 bytes)")
    else:
        logger.warning(f"[{sample_name}] Scaffolds file does not exist")
    
    # If we get here, scaffolds file has no valid content - try fallback 1
    logger.warning(f"[{sample_name}] Initial assembly failed. Initiating Fallback 1: single-reads k=21 assembly")
    logger.info(f"[{sample_name}] Removing failed assembly directory: {output_path}")
    shutil.rmtree(output_path)
    logger.info(f"[{sample_name}] Starting Fallback 1: single-reads k=21 assembly")
    run_single_spades(merged_path, output_path, logger, k="21")
    logger.info(f"[{sample_name}] Fallback 1 assembly finished")
    return False


def check_spades_out2(merged_path, unmerged_1, unmerged_2, output_path, logger, k):
    """Check if SPAdes fallback output is valid; if not, rerun with paired-end k=21."""

    scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
    
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    logger.info(f"[{sample_name}] Checking Fallback 1 assembly output")
    
    # Check scaffolds file
    if os.path.isfile(scaffolds_file):
        file_size = os.path.getsize(scaffolds_file)
        
        if file_size > 0:
            # Count number of sequences and get metrics
            with open(scaffolds_file, 'r') as f:
                seq_count = sum(1 for line in f if line.startswith('>'))
            
            if seq_count > 0:
                n_scaffolds, scaffold_lengths, mean_length, n50 = get_scaffold_metrics(scaffolds_file)
                logger.info(f"[{sample_name}] ✓ Fallback 1 assembly successful - {seq_count} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp ({file_size} bytes)")
                return True
            else:
                logger.warning(f"[{sample_name}] Scaffolds file is empty (no sequences found)")
        else:
            logger.warning(f"[{sample_name}] Scaffolds file exists but is empty (0 bytes)")
    else:
        logger.warning(f"[{sample_name}] Scaffolds file does not exist")
    
    # If we get here, fallback 1 failed - try final fallback with paired-end k=21
    logger.error(f"[{sample_name}] Fallback 1 assembly failed. Initiating Fallback 2: merged + unmerged pairs k=21 assembly")
    logger.info(f"[{sample_name}] Removing failed assembly directory: {output_path}")
    shutil.rmtree(output_path)
    logger.info(f"[{sample_name}] Starting Fallback 2: merged + unmerged pairs k=21 assembly")
    run_spades(merged_path, unmerged_1, unmerged_2, output_path, logger, k="21")
    logger.info(f"[{sample_name}] Fallback 2 assembly finished")
    return False


def spades_wrapper(job):
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    try:
        logger.info(f"[{sample_name}] Starting initial assembly: single-reads k={k} assembly")
        run_single_spades(merged_path, output_path, logger, k)
        logger.info(f"[{sample_name}] Initial assembly finished")
        return (merged_path, True, sample_name)
    except subprocess.CalledProcessError as e:
        logger.error(f"[{sample_name}] Initial assembly failed with error: {str(e)}")
        return (merged_path, False, sample_name)


def check_wrapper(job):
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    success = check_spades_out1(merged_path, unmerged_1, unmerged_2, output_path, logger, k)
    return (merged_path, success, sample_name)


def check2_wrapper(job):
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    success = check_spades_out2(merged_path, unmerged_1, unmerged_2, output_path, logger, k)
    return (merged_path, success, sample_name)


def write_summary_csv(jobs, initial_pass, fallback1_pass, fallback2_pass, fallback2_fail, 
                      output_dir, csv_path, logger):
    """Write a summary CSV with assembly metrics for all samples."""
    
    logger.info(f"\nWriting summary CSV to: {csv_path}")
    
    csv_data = []
    
    for job in jobs:
        merged_path, unmerged_1, unmerged_2, output_path, _, k = job
        sample_name = os.path.basename(output_path).replace(".spades.out", "")
        scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
        
        # Determine assembly method and status
        if sample_name in initial_pass:
            assembly_method = "INITIAL"
            assembly_status = "PASS"
        elif sample_name in fallback1_pass:
            assembly_method = "FALLBACK1"
            assembly_status = "PASS"
        elif sample_name in fallback2_pass:
            assembly_method = "FALLBACK2"
            assembly_status = "PASS"
        elif sample_name in fallback2_fail:
            assembly_method = "FAIL"
            assembly_status = "FAIL"
        else:
            # Should not happen, but handle gracefully
            assembly_method = "FAIL"
            assembly_status = "UNKNOWN"
        
        # Get scaffold metrics if assembly passed
        if assembly_status == "PASS":
            n_scaffolds, scaffold_lengths, mean_length, n50 = get_scaffold_metrics(scaffolds_file)
            scaffold_lengths_str = ";".join(map(str, scaffold_lengths))
            mean_scaffold_length = f"{mean_length:.2f}"
            n50_str = str(n50)
        else:
            n_scaffolds = ""
            scaffold_lengths_str = ""
            mean_scaffold_length = ""
            n50_str = ""
        
        csv_data.append({
            "ID": sample_name,
            "assembly_method": assembly_method,
            "assembly_status": assembly_status,
            "n_scaffolds": n_scaffolds,
            "scaffold_lengths": scaffold_lengths_str,
            "mean_scaffold_length": mean_scaffold_length,
            "n50": n50_str
        })
    
    # Write CSV
    fieldnames = ["ID", "assembly_method", "assembly_status", "n_scaffolds", 
                  "scaffold_lengths", "mean_scaffold_length", "n50"]
    
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_data)
    
    logger.info(f"Summary CSV written successfully with {len(csv_data)} samples")


def main(args):    
    # Setup
    os.makedirs(args.output_dir, exist_ok=True)
    log_dir = os.path.join(args.output_dir, "logs")

    # Handle log file naming
    if args.log_file:
        log_file_name = os.path.basename(args.log_file)  # Extract just filename
    else:
        log_file_name = f'assembly_log_{datetime.now().strftime("%Y%m%d_%H%M%S")}.log'

    logger = setup_logging(log_dir=log_dir, log_file=log_file_name)

    # Get all mapped files recursively
    merged_files = glob.glob(os.path.join(args.merged_dir, "**/*_mapped.fastq"), recursive=True)

    jobs = []
    for merged_path in merged_files:
        filename = os.path.basename(merged_path).replace("_mapped.fastq", "")

        # Search for matching unmerged files
        unmerged_1_matches = glob.glob(
            os.path.join(args.unmerged_dir, "**", f"{filename}*_unmerged_1.f*q"),
            recursive=True
        )
        unmerged_2_matches = glob.glob(
            os.path.join(args.unmerged_dir, "**", f"{filename}*_unmerged_2.f*q"),
            recursive=True
        )

        if not unmerged_1_matches or not unmerged_2_matches:
            logger.error(f"Missing unmerged files for: {filename}")
            continue

        unmerged_1 = unmerged_1_matches[0]
        unmerged_2 = unmerged_2_matches[0]

        output_path = os.path.join(args.output_dir, f"{filename}.spades.out")

        # Store only needed args for multiprocessing
        jobs.append((merged_path, unmerged_1, unmerged_2, output_path, logger, args.k))

    logger.info(f"-" * 40)
    logger.info(f"Starting SPAdes assembly for {len(jobs)} samples")
    logger.info(f"Strategy: single-reads k={args.k} → single-reads k=21 → merged+unmerged k=21")

    # Track results at each stage
    # Stage tracking: sample_name -> status
    initial_pass = set()  # Samples that passed initial assembly validation
    initial_fail = set()  # Samples that failed initial assembly (needed Fallback 1)
    fallback1_pass = set()  # Samples that passed Fallback 1 validation
    fallback1_fail = set()  # Samples that failed Fallback 1 (needed Fallback 2)
    fallback2_pass = set()  # Samples that passed Fallback 2 validation
    fallback2_fail = set()  # Samples that failed Fallback 2 (ultimate failure)
    
    # Initial SPAdes run
    logger.info(f"{'-' * 80}")
    logger.info(f"STAGE 1: Initial assembly (single-reads k={args.k})")
    logger.info(f"{'-' * 80}")
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(spades_wrapper, job) for job in jobs]
        for future in as_completed(futures):
            try:
                merged_path, success, sample_name = future.result()
                # Note: success here just means SPAdes didn't crash, not that output is valid
            except Exception as e:
                logger.error(f"Job failed with error: {str(e)}")

    # First check to validate and fall-back if needed (single-reads k=21)
    logger.info(f"{'-' * 80}")
    logger.info(f"STAGE 2: Validation and Fallback 1 (single-reads k=21)")
    logger.info(f"{'-' * 80}")
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(check_wrapper, job) for job in jobs]
        for future in as_completed(futures):
            try:
                merged_path, success, sample_name = future.result()
                if success:
                    # Passed initial validation - no fallback needed
                    initial_pass.add(sample_name)
                else:
                    # Failed initial validation - Fallback 1 was triggered
                    initial_fail.add(sample_name)
            except Exception as e:
                logger.error(f"Check 1 failed with error: {str(e)}")

    # Second check to validate Fallback 1 and trigger Fallback 2 if needed
    logger.info(f"{'-' * 80}")
    logger.info(f"STAGE 3: Validation and Fallback 2 (merged+unmerged k=21)")
    logger.info(f"{'-' * 80}")
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = [executor.submit(check2_wrapper, job) for job in jobs]
        for future in as_completed(futures):
            try:
                merged_path, success, sample_name = future.result()
                # Only samples that needed Fallback 1 reach this check
                if sample_name in initial_fail:
                    if success:
                        # Passed Fallback 1 validation
                        fallback1_pass.add(sample_name)
                    else:
                        # Failed Fallback 1 validation - Fallback 2 was triggered
                        fallback1_fail.add(sample_name)
            except Exception as e:
                logger.error(f"Check 2 failed with error: {str(e)}")

    # Final validation to check Fallback 2 results
    for job in jobs:
        merged_path, unmerged_1, unmerged_2, output_path, _, k = job
        sample_name = os.path.basename(output_path).replace(".spades.out", "")
        
        # Only check samples that went through Fallback 2
        if sample_name in fallback1_fail:
            scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
            
            has_valid_output = False
            if os.path.isfile(scaffolds_file) and os.path.getsize(scaffolds_file) > 0:
                with open(scaffolds_file, 'r') as f:
                    if sum(1 for line in f if line.startswith('>')) > 0:
                        has_valid_output = True
            
            if has_valid_output:
                fallback2_pass.add(sample_name)
                n_scaffolds, scaffold_lengths, mean_length, n50 = get_scaffold_metrics(scaffolds_file)
                logger.info(f"[{sample_name}] ✓ Fallback 2 assembly successful - {n_scaffolds} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp")
            else:
                fallback2_fail.add(sample_name)
                logger.error(f"[{sample_name}] ✗ Fallback 2 assembly failed - no valid output")

    # Calculate totals
    total_samples = len(jobs)
    total_successful = len(initial_pass) + len(fallback1_pass) + len(fallback2_pass)
    total_failed = len(fallback2_fail)

    # Print summary
    logger.info(f"\n{'-' * 80}")
    logger.info(f"ASSEMBLY SUMMARY")
    logger.info(f"{'-' * 80}")
    logger.info(f"Total samples processed: {total_samples}")
    logger.info(f"")
    logger.info(f"Initial assembly (single-reads k={args.k}):")
    logger.info(f"  ✓ Passed: {len(initial_pass)}")
    logger.info(f"  ✗ Failed (needed Fallback 1): {len(initial_fail)}")
    logger.info(f"")
    
    if initial_fail:
        logger.info(f"Fallback 1 (single-reads k=21):")
        logger.info(f"  ✓ Passed: {len(fallback1_pass)}")
        logger.info(f"  ✗ Failed (needed Fallback 2): {len(fallback1_fail)}")
        logger.info(f"")
    
    if fallback1_fail:
        logger.info(f"Fallback 2 (merged+unmerged k=21):")
        logger.info(f"  ✓ Passed: {len(fallback2_pass)}")
        logger.info(f"  ✗ Failed: {len(fallback2_fail)}")
        logger.info(f"")
    
    logger.info(f"FINAL RESULTS:")
    logger.info(f"  ✓ Total successful assemblies: {total_successful}")
    logger.info(f"  ✗ Total failed assemblies: {total_failed}")
    
    if initial_fail:
        logger.info(f"")
        logger.info(f"Samples that needed Fallback 1:")
        for sample in sorted(initial_fail):
            logger.info(f"  - {sample}")
    
    if fallback1_fail:
        logger.info(f"")
        logger.info(f"Samples that needed Fallback 2:")
        for sample in sorted(fallback1_fail):
            logger.info(f"  - {sample}")
    
    if fallback2_fail:
        logger.error(f"")
        logger.error(f"Samples that failed all assembly attempts:")
        for sample in sorted(fallback2_fail):
            logger.error(f"  - {sample}")
    
    logger.info(f"{'-' * 80}")
    
    # Write CSV summary if requested
    if args.summary_csv:
        write_summary_csv(jobs, initial_pass, fallback1_pass, fallback2_pass, 
                         fallback2_fail, args.output_dir, args.summary_csv, logger)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Run SPAdes on merged + unmerged reads in parallel.")
    parser.add_argument("--merged_dir", required=True, help="Path to merged mapped reads (e.g., *_mapped.fastq)")
    parser.add_argument("--unmerged_dir", required=True, help="Path to unmerged read pairs (e.g., *_unmerged_1.fq, *_unmerged_2.fq)")
    parser.add_argument("--output_dir", required=True, help="Directory to store SPAdes output.")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel SPAdes jobs.")
    parser.add_argument("--k", type=str, default="21,33,55", help="Comma-separated k-mer sizes for SPAdes.")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `assembly_log` file will be created.")
    parser.add_argument("--summary_csv", help="Optional path to write a CSV summary of assembly results.")

    args = parser.parse_args()
    main(args)

## Example usage: 
#python assembly_module_v7.py \
#  --merged_dir ./TrueITSref_mapped_reads \
#  --unmerged_dir ./fastp_processed_ts \
#  --output_dir ./TrueITSref_assemblies \
#  --summary_csv ./assembly_summary.csv
