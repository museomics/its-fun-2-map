import os
import glob
import subprocess
import argparse
import shutil
import csv
from datetime import datetime
from concurrent.futures import ProcessPoolExecutor, as_completed
from metahist_tools import setup_logging

#### Author: Maria Kamouyiaros & Daniel Parsons @ NHMUK
#### Date: 2025-08-04
#### VERSION: 7.3.0

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


def check_assembly_output(scaffolds_file):
    """
    Check if a scaffolds.fasta file exists and contains valid sequences.
    
    Returns:
        tuple: (is_valid, n_scaffolds, mean_length, n50)
    """
    if not os.path.isfile(scaffolds_file):
        return False, 0, 0, 0
    
    if os.path.getsize(scaffolds_file) == 0:
        return False, 0, 0, 0
    
    n_scaffolds, scaffold_lengths, mean_length, n50 = get_scaffold_metrics(scaffolds_file)
    
    if n_scaffolds > 0:
        return True, n_scaffolds, mean_length, n50
    else:
        return False, 0, 0, 0


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
    subprocess.run(cmd, check=True, stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    return output_path


def run_initial_assembly(job):
    """Run initial assembly (Stage 1) for a sample."""
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    
    try:
        run_single_spades(merged_path, output_path, logger, k)
        
        # Check if assembly produced valid output
        scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
        is_valid, n_scaffolds, mean_length, n50 = check_assembly_output(scaffolds_file)
        
        if is_valid:
            return (sample_name, "PASS", n_scaffolds, mean_length, n50)
        else:
            return (sample_name, "FAIL", 0, 0, 0)
            
    except subprocess.CalledProcessError as e:
        return (sample_name, "ERROR", 0, 0, 0)


def run_fallback1_assembly(job):
    """Run Fallback 1 assembly (single-reads k=21) for a sample."""
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    
    try:
        # Remove failed assembly directory
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        
        run_single_spades(merged_path, output_path, logger, k="21")
        
        # Check if assembly produced valid output
        scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
        is_valid, n_scaffolds, mean_length, n50 = check_assembly_output(scaffolds_file)
        
        if is_valid:
            return (sample_name, "PASS", n_scaffolds, mean_length, n50)
        else:
            return (sample_name, "FAIL", 0, 0, 0)
            
    except subprocess.CalledProcessError as e:
        return (sample_name, "ERROR", 0, 0, 0)


def run_fallback2_assembly(job):
    """Run Fallback 2 assembly (merged + unmerged pairs k=21) for a sample."""
    merged_path, unmerged_1, unmerged_2, output_path, logger, k = job
    sample_name = os.path.basename(output_path).replace(".spades.out", "")
    
    try:
        # Remove failed assembly directory
        if os.path.exists(output_path):
            shutil.rmtree(output_path)
        
        run_spades(merged_path, unmerged_1, unmerged_2, output_path, logger, k="21")
        
        # Check if assembly produced valid output
        scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
        is_valid, n_scaffolds, mean_length, n50 = check_assembly_output(scaffolds_file)
        
        if is_valid:
            return (sample_name, "PASS", n_scaffolds, mean_length, n50)
        else:
            return (sample_name, "FAIL", 0, 0, 0)
            
    except subprocess.CalledProcessError as e:
        return (sample_name, "ERROR", 0, 0, 0)


def write_summary_csv(jobs, sample_paths, sample_metrics, output_dir, csv_path, logger):
    """Write a summary CSV with assembly metrics for all samples."""
    
    logger.info(f"Writing summary CSV to: {csv_path}")
    
    csv_data = []
    
    for job in jobs:
        merged_path, unmerged_1, unmerged_2, output_path, _, k = job
        sample_name = os.path.basename(output_path).replace(".spades.out", "")
        
        path = sample_paths.get(sample_name, [])
        metrics = sample_metrics.get(sample_name, (0, 0, 0))
        
        # Determine final status and method
        if path and path[-1].endswith("PASS"):
            assembly_status = "PASS"
            # Get the method from the last successful stage
            if "INITIAL: PASS" in path:
                assembly_method = "INITIAL"
            elif "FALLBACK1: PASS" in path:
                assembly_method = "FALLBACK1"
            elif "FALLBACK2: PASS" in path:
                assembly_method = "FALLBACK2"
            else:
                assembly_method = "UNKNOWN"
            
            n_scaffolds, mean_length, n50 = metrics
            scaffold_lengths_str = ""
            
            # Get scaffold lengths from file
            scaffolds_file = os.path.join(output_path, "scaffolds.fasta")
            if os.path.isfile(scaffolds_file):
                _, scaffold_lengths, _, _ = get_scaffold_metrics(scaffolds_file)
                scaffold_lengths_str = ";".join(map(str, scaffold_lengths))
            
            mean_scaffold_length = f"{mean_length:.2f}"
            n50_str = str(n50)
        else:
            assembly_status = "FAIL"
            assembly_method = "FAIL"
            n_scaffolds = ""
            scaffold_lengths_str = ""
            mean_scaffold_length = ""
            n50_str = ""
        
        csv_data.append({
            "ID": sample_name,
            "method": assembly_method,
            "status": assembly_status,
            "n_scaffolds": n_scaffolds,
            "scaffold_lengths": scaffold_lengths_str,
            "mean_scaffold_length": mean_scaffold_length,
            "n50": n50_str
        })
    
    # Write CSV
    fieldnames = ["ID", "method", "status", "n_scaffolds", 
                  "scaffold_lengths", "mean_scaffold_length", "n50"]
    
    with open(csv_path, 'w', newline='') as csvfile:
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(csv_data)
    
    logger.info(f"Summary CSV written successfully with {len(csv_data)} samples")


def main(args):    
    # Setup
    os.makedirs(args.output_dir, exist_ok=True)

    # Handle log file naming
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'assembly_log_{timestamp}.log'

    logger = setup_logging(log_file=args.log_file)

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
            logger.error(f"[{filename}] Missing unmerged files - skipping sample")
            continue

        unmerged_1 = unmerged_1_matches[0]
        unmerged_2 = unmerged_2_matches[0]

        output_path = os.path.join(args.output_dir, f"{filename}.spades.out")
        jobs.append((merged_path, unmerged_1, unmerged_2, output_path, logger, args.k))

    # Track sample paths and metrics
    sample_paths = {os.path.basename(job[3]).replace(".spades.out", ""): [] for job in jobs}
    sample_metrics = {}  # sample_name -> (n_scaffolds, mean_length, n50)
    
    # Track which samples need each stage
    samples_for_fallback1 = set()
    samples_for_fallback2 = set()
    
    logger.info("=" * 80)
    logger.info("SPADES ASSEMBLY PIPELINE")
    logger.info("=" * 80)
    logger.info(f"Total samples to process: {len(jobs)}")
    logger.info(f"Strategy: INITIAL (single-reads k={args.k})")
    logger.info(f"          â†’ FALLBACK1 (single-reads k=21)")
    logger.info(f"          â†’ FALLBACK2 (merged+unmerged k=21)")
    logger.info("=" * 80)

    # =========================================================================
    # STAGE 1: Initial Assembly
    # =========================================================================
    logger.info("")
    logger.info("=" * 80)
    logger.info("STAGE 1: INITIAL ASSEMBLY (single-reads k={})".format(args.k))
    logger.info("=" * 80)
    
    with ProcessPoolExecutor(max_workers=args.threads) as executor:
        futures = {executor.submit(run_initial_assembly, job): job for job in jobs}
        
        for future in as_completed(futures):
            try:
                sample_name, status, n_scaffolds, mean_length, n50 = future.result()
                
                if status == "PASS":
                    logger.info(f"[{sample_name}] INITIAL: PASS âœ“ ({n_scaffolds} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp)")
                    sample_paths[sample_name].append("INITIAL: PASS")
                    sample_metrics[sample_name] = (n_scaffolds, mean_length, n50)
                elif status == "FAIL":
                    logger.warning(f"[{sample_name}] INITIAL: FAIL âœ— â†’ triggering FALLBACK1")
                    sample_paths[sample_name].append("INITIAL: FAIL")
                    samples_for_fallback1.add(sample_name)
                else:  # ERROR
                    logger.error(f"[{sample_name}] INITIAL: ERROR âœ— (SPAdes crashed) â†’ triggering FALLBACK1")
                    sample_paths[sample_name].append("INITIAL: ERROR")
                    samples_for_fallback1.add(sample_name)
                    
            except Exception as e:
                logger.error(f"Job failed with unexpected error: {str(e)}")

    # Stage 1 summary
    stage1_pass = len(jobs) - len(samples_for_fallback1)
    logger.info("")
    logger.info(f"STAGE 1 COMPLETE: {stage1_pass}/{len(jobs)} passed, {len(samples_for_fallback1)} need FALLBACK1")

    # =========================================================================
    # STAGE 2: Fallback 1 (if needed)
    # =========================================================================
    if samples_for_fallback1:
        logger.info("")
        logger.info("=" * 80)
        logger.info("STAGE 2: FALLBACK1 ASSEMBLY (single-reads k=21)")
        logger.info("=" * 80)
        logger.info(f"Processing {len(samples_for_fallback1)} samples that failed INITIAL assembly")
        logger.info("")
        
        # Filter jobs for samples that need fallback1
        fallback1_jobs = [job for job in jobs 
                         if os.path.basename(job[3]).replace(".spades.out", "") in samples_for_fallback1]
        
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = {executor.submit(run_fallback1_assembly, job): job for job in fallback1_jobs}
            
            for future in as_completed(futures):
                try:
                    sample_name, status, n_scaffolds, mean_length, n50 = future.result()
                    
                    if status == "PASS":
                        logger.info(f"[{sample_name}] FALLBACK1: PASS âœ“ ({n_scaffolds} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp)")
                        sample_paths[sample_name].append("FALLBACK1: PASS")
                        sample_metrics[sample_name] = (n_scaffolds, mean_length, n50)
                    elif status == "FAIL":
                        logger.warning(f"[{sample_name}] FALLBACK1: FAIL âœ— â†’ triggering FALLBACK2")
                        sample_paths[sample_name].append("FALLBACK1: FAIL")
                        samples_for_fallback2.add(sample_name)
                    else:  # ERROR
                        logger.error(f"[{sample_name}] FALLBACK1: ERROR âœ— (SPAdes crashed) â†’ triggering FALLBACK2")
                        sample_paths[sample_name].append("FALLBACK1: ERROR")
                        samples_for_fallback2.add(sample_name)
                        
                except Exception as e:
                    logger.error(f"Fallback 1 job failed with unexpected error: {str(e)}")

        # Stage 2 summary
        stage2_pass = len(samples_for_fallback1) - len(samples_for_fallback2)
        logger.info("")
        logger.info(f"STAGE 2 COMPLETE: {stage2_pass}/{len(samples_for_fallback1)} passed, {len(samples_for_fallback2)} need FALLBACK2")

    # =========================================================================
    # STAGE 3: Fallback 2 (if needed)
    # =========================================================================
    if samples_for_fallback2:
        logger.info("")
        logger.info("=" * 80)
        logger.info("STAGE 3: FALLBACK2 ASSEMBLY (merged+unmerged k=21)")
        logger.info("=" * 80)
        logger.info(f"Processing {len(samples_for_fallback2)} samples that failed FALLBACK1 assembly")
        logger.info("")
        
        # Filter jobs for samples that need fallback2
        fallback2_jobs = [job for job in jobs 
                         if os.path.basename(job[3]).replace(".spades.out", "") in samples_for_fallback2]
        
        with ProcessPoolExecutor(max_workers=args.threads) as executor:
            futures = {executor.submit(run_fallback2_assembly, job): job for job in fallback2_jobs}
            
            for future in as_completed(futures):
                try:
                    sample_name, status, n_scaffolds, mean_length, n50 = future.result()
                    
                    if status == "PASS":
                        logger.info(f"[{sample_name}] FALLBACK2: PASS âœ“ ({n_scaffolds} scaffolds, N50={n50} bp, mean={mean_length:.2f} bp)")
                        sample_paths[sample_name].append("FALLBACK2: PASS")
                        sample_metrics[sample_name] = (n_scaffolds, mean_length, n50)
                    elif status == "FAIL":
                        logger.error(f"[{sample_name}] FALLBACK2: FAIL âœ— (all assembly attempts exhausted)")
                        sample_paths[sample_name].append("FALLBACK2: FAIL")
                    else:  # ERROR
                        logger.error(f"[{sample_name}] FALLBACK2: ERROR âœ— (SPAdes crashed, all attempts exhausted)")
                        sample_paths[sample_name].append("FALLBACK2: ERROR")
                        
                except Exception as e:
                    logger.error(f"Fallback 2 job failed with unexpected error: {str(e)}")

        # Stage 3 summary
        stage3_pass = sum(1 for s in samples_for_fallback2 
                        if sample_paths[s][-1] in ["FALLBACK2: PASS"])
        stage3_fail = len(samples_for_fallback2) - stage3_pass
        logger.info("")
        logger.info(f"STAGE 3 COMPLETE: {stage3_pass}/{len(samples_for_fallback2)} passed, {stage3_fail} failed all attempts")

    # =========================================================================
    # FINAL SUMMARY
    # =========================================================================
    logger.info("")
    logger.info("=" * 80)
    logger.info("ASSEMBLY PIPELINE COMPLETE")
    logger.info("=" * 80)
    
    # Count final results
    total_pass = sum(1 for paths in sample_paths.values() if paths and paths[-1].endswith("PASS"))
    total_fail = len(jobs) - total_pass
    
    initial_pass = sum(1 for paths in sample_paths.values() if "INITIAL: PASS" in paths)
    fallback1_pass = sum(1 for paths in sample_paths.values() if "FALLBACK1: PASS" in paths)
    fallback2_pass = sum(1 for paths in sample_paths.values() if "FALLBACK2: PASS" in paths)
    
    logger.info("")
    logger.info("SUMMARY BY STAGE:")
    logger.info(f"  INITIAL:   {initial_pass} passed")
    if samples_for_fallback1:
        logger.info(f"  FALLBACK1: {fallback1_pass} passed")
    if samples_for_fallback2:
        logger.info(f"  FALLBACK2: {fallback2_pass} passed")
    logger.info("")
    logger.info("FINAL RESULTS:")
    logger.info(f"  âœ“ Total successful: {total_pass}/{len(jobs)}")
    logger.info(f"  âœ— Total failed:     {total_fail}/{len(jobs)}")
    
    # Per-sample path summary
    logger.info("")
    logger.info("-" * 80)
    logger.info("PER-SAMPLE ASSEMBLY PATH:")
    logger.info("-" * 80)
    
    for sample_name in sorted(sample_paths.keys()):
        path = sample_paths[sample_name]
        path_str = " â†’ ".join(path)
        
        # Determine final status symbol
        if path and path[-1].endswith("PASS"):
            final_status = "âœ“"
        else:
            final_status = "âœ—"
        
        logger.info(f"  {sample_name}: {path_str} {final_status}")
    
    logger.info("-" * 80)
    
    # Write CSV summary if requested
    if args.summary_csv:
        logger.info("")
        write_summary_csv(jobs, sample_paths, sample_metrics, 
                         args.output_dir, args.summary_csv, logger)


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
#python assembly_module.py \
#  --merged_dir ./mapped_reads \
#  --unmerged_dir ./fastp_processed \
#  --output_dir ./assemblies \
#  --summary_csv ./assembly_summary.csv
