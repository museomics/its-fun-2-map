import os 
import sys
import subprocess
import logging
import pathlib
import shutil
import glob
import argparse
import pandas as pd
import re
from datetime import datetime
from concurrent.futures import ThreadPoolExecutor, as_completed
from its_fun_tools import load_name_ids, cleanup_temp_dir, log_and_print

sys.path.append('/home/mkamouyi/scratch/private/defra-fungi/PRJEB81712')

from metahist_tools import repair_reads, setup_logging

#### Date: 2025-08-04
#### Purpose: tweaked version of decontam module to keep mapped reads instead of bin them (changed samtools flag and removed previous phix filtering step)
####        : possible update to decontam module: make phix optional and add --mapped versus --unmapped reads
#### Author: Maria Kamouyiaros
#### Important note: metahist_tools has an absolute path at the moment until modules are fully packaged - will need to change this after python package is finished
#### Update 2025-08-28: Doesn't use --suffix, --prefix, --paired anymore; main edited so that reference (seed) is paired with each sample FASTQ based on ID
####              : Includes step to index reference files if they're not already indexed
#### Update 2025-10-08: Tweaked it so that it maps merged and unmerged reads to references and checks that there are no more than 3 found files (1 merged, 2 unmerged) 
####              : Added logging to a file with timestamp if not provided by user and kept consistent with blast_parser logging
#### Update 2025-11-17: Added --aligner argument to support both bwa-mem and bwa-aln algorithms
#### Update 2025-11-18: Modified bwa-aln to output separate flagstats files and added mapping_summary.csv output
#### Update 2025-11-26: Added empty file checks, improved logging, and ancient DNA-specific BWA ALN parameters (-l 16500 -n 0.01 -o 2)


def parse_flagstats(flagstats_file):
    """Parse a samtools flagstats file to extract mapped read count and percentage.
    Returns tuple: (n_mapped, mapped_percentage)
    """
    if not os.path.exists(flagstats_file):
        return (None, None)
    
    with open(flagstats_file, 'r') as f:
        lines = f.readlines()
        n_mapped = None
        mapped_perc = None
        
        # Look for the "mapped" line (not "primary mapped")
        for i, line in enumerate(lines):
            if '\tmapped\n' in line:
                parts = line.strip().split('\t')
                n_mapped = int(parts[0])
            elif 'mapped %' in line and 'primary mapped %' not in line:
                parts = line.strip().split('\t')
                perc_str = parts[0].replace('%', '').strip()
                if perc_str == 'N/A':
                    mapped_perc = None
                else:
                    mapped_perc = float(perc_str)
    
    return (n_mapped, mapped_perc)


def run_bwa_mem_and_samtools(sample_id, input_file, output_dir, temp_dir, ref, threads, logger):

    ''' Function that runs BWA MEM and SAMtools to map reads to a reference and parse the mapped reads to an output FASTQ file.'''

    if input_file.endswith("merged.fq"):
        bam_file = os.path.join(temp_dir, f"{sample_id}_mapped.bam")
        sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_sorted.bam")
        mapped_fastq = os.path.join(output_dir, f"{sample_id}_mapped.fastq")
        stats_file = os.path.join(output_dir, f"{sample_id}_mapped_flagstats.txt")
    elif input_file.endswith("unmerged_1.fq"):
        bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_unmerged_1.bam")
        sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_sorted_unmerged_1.bam")
        mapped_fastq = os.path.join(output_dir, f"{sample_id}_mapped_unmerged_1.fastq")
        stats_file = os.path.join(output_dir, f"{sample_id}_unmerged_1_flagstats.txt")
    elif input_file.endswith("unmerged_2.fq"):
        bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_unmerged_2.bam")
        sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_sorted_unmerged_2.bam")
        mapped_fastq = os.path.join(output_dir, f"{sample_id}_mapped_unmerged_2.fastq")
        stats_file = os.path.join(output_dir, f"{sample_id}_unmerged_2_flagstats.txt")

    try:
        logger.info(f"Running BWA MEM and SAMtools for {input_file} with {ref}")
        bwa_cmd = ["bwa", "mem", "-M", "-t", str(threads), ref, input_file]
        with open(bam_file, "wb") as bam_out:
            bwa_proc = subprocess.Popen(bwa_cmd, stdout=subprocess.PIPE)
            view_proc = subprocess.Popen(["samtools", "view", "-b", "-"], stdin=bwa_proc.stdout, stdout=bam_out)
            bwa_proc.wait()
            view_proc.communicate()

        subprocess.run(["samtools", "sort", "-o", sorted_bam_file, bam_file], check=True)
        logger.info(f"Saved {sorted_bam_file}")

        with open(mapped_fastq, "w") as fastq_out:
            subprocess.run(
                ["samtools", "fastq", "-F", "4", sorted_bam_file],
                stdout=fastq_out,
                check=True
            )

        with open(stats_file, "w") as out:
            subprocess.run(["samtools", "flagstat", "-O", "tsv", sorted_bam_file], stdout=out, check=True)
            logger.info(f"Finished BWA MEM/SAMtools for {sample_id}")
    except Exception as e:
        logger.error(f"Error processing {input_file}: {e}")
        raise


def run_bwa_aln_and_samtools(sample_id, input_file, output_dir, temp_dir, ref, threads, logger, paired_file=None):
    
    """Function that runs BWA ALN and SAMtools to map reads to a reference and parse the mapped reads to an output FASTQ file.
    For paired-end reads, paired_file should be provided.
    
    Uses ancient DNA-specific parameters (Green et al. 2010):
    - -l 16500: Disable seeding (avoids rejecting reads with damage in seed region)
    - -n 0.01: Stringent edit distance (reduces false positives with seeding disabled)
    - -o 2: Allow 2 gap opens (accommodates small indels from degradation)
    """

    if input_file.endswith("merged.fq"):
        # Single-end merged reads
        sai_file = os.path.join(temp_dir, f"{sample_id}_mapped.sai")
        bam_file = os.path.join(temp_dir, f"{sample_id}_mapped.bam")
        sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_sorted.bam")
        mapped_fastq = os.path.join(output_dir, f"{sample_id}_mapped.fastq")
        stats_file = os.path.join(output_dir, f"{sample_id}_mapped_flagstats.txt")
        
        try:
            logger.info(f"[{sample_id}] Starting BWA ALN (single-end) for {os.path.basename(input_file)}")
            
            # Get input file size for logging
            input_size_mb = os.path.getsize(input_file) / (1024 * 1024)
            logger.info(f"[{sample_id}] Input file size: {input_size_mb:.1f} MB")
            
            # Run bwa aln with ancient DNA-specific parameters (Green et al. 2010)
            aln_cmd = [
                "bwa", "aln",
                "-t", str(threads),
                "-l", "16500",  # Disable seeding - critical for ancient DNA with damage patterns
                "-n", "0.01",   # Stringent edit distance (works because seeding is disabled)
                "-o", "2",      # Allow 2 gap opens for degradation-related indels
                ref, input_file
            ]
            logger.info(f"[{sample_id}] BWA ALN command: {' '.join(aln_cmd)}")
            
            with open(sai_file, "w") as sai_out:
                subprocess.run(aln_cmd, stdout=sai_out, check=True)
            
            logger.info(f"[{sample_id}] BWA ALN complete, running samse...")
            
            # Run bwa samse and pipe to samtools view
            with open(bam_file, "wb") as bam_out:
                samse_proc = subprocess.Popen(["bwa", "samse", ref, sai_file, input_file], stdout=subprocess.PIPE)
                view_proc = subprocess.Popen(["samtools", "view", "-b", "-"], stdin=samse_proc.stdout, stdout=bam_out)
                samse_proc.wait()
                view_proc.communicate()
            
            logger.info(f"[{sample_id}] Sorting BAM...")
            
            # Sort BAM
            subprocess.run(["samtools", "sort", "-o", sorted_bam_file, bam_file], check=True)
            logger.info(f"[{sample_id}] Saved {sorted_bam_file}")
            
            # Extract mapped reads to FASTQ
            with open(mapped_fastq, "w") as fastq_out:
                subprocess.run(
                    ["samtools", "fastq", "-F", "4", sorted_bam_file],
                    stdout=fastq_out,
                    check=True
                )
            
            # Generate flagstats
            with open(stats_file, "w") as out:
                subprocess.run(["samtools", "flagstat", "-O", "tsv", sorted_bam_file], stdout=out, check=True)
            
            # Clean up .sai file
            if os.path.exists(sai_file):
                os.remove(sai_file)
            
            logger.info(f"[{sample_id}] Finished BWA ALN/SAMtools for merged reads")
            
        except Exception as e:
            logger.error(f"[{sample_id}] Error processing {input_file}: {e}")
            raise
            
    elif input_file.endswith("unmerged_1.fq") and paired_file and paired_file.endswith("unmerged_2.fq"):
        # Paired-end unmerged reads - process both together
        sai_file_1 = os.path.join(temp_dir, f"{sample_id}_mapped_unmerged_1.sai")
        sai_file_2 = os.path.join(temp_dir, f"{sample_id}_mapped_unmerged_2.sai")
        bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_unmerged.bam")
        sorted_bam_file = os.path.join(temp_dir, f"{sample_id}_mapped_sorted_unmerged.bam")
        mapped_fastq_1 = os.path.join(output_dir, f"{sample_id}_mapped_unmerged_1.fastq")
        mapped_fastq_2 = os.path.join(output_dir, f"{sample_id}_mapped_unmerged_2.fastq")
        stats_file_1 = os.path.join(output_dir, f"{sample_id}_unmerged_1_flagstats.txt")
        stats_file_2 = os.path.join(output_dir, f"{sample_id}_unmerged_2_flagstats.txt")
        
        try:
            logger.info(f"[{sample_id}] Starting BWA ALN (paired-end) for {os.path.basename(input_file)} and {os.path.basename(paired_file)}")
            
            # Get input file sizes for logging
            input1_size_mb = os.path.getsize(input_file) / (1024 * 1024)
            input2_size_mb = os.path.getsize(paired_file) / (1024 * 1024)
            logger.info(f"[{sample_id}] Input file sizes: R1={input1_size_mb:.1f} MB, R2={input2_size_mb:.1f} MB")
            
            # BWA ALN command with ancient DNA-specific parameters (Green et al. 2010)
            aln_cmd_base = [
                "bwa", "aln",
                "-t", str(threads),
                "-l", "16500",  # Disable seeding
                "-n", "0.01",   # Stringent edit distance
                "-o", "2",      # Allow 2 gap opens
                ref
            ]
            
            # Run bwa aln for read 1
            logger.info(f"[{sample_id}] Running BWA ALN for R1...")
            with open(sai_file_1, "w") as sai_out:
                subprocess.run(aln_cmd_base + [input_file], stdout=sai_out, check=True)
            
            # Run bwa aln for read 2
            logger.info(f"[{sample_id}] Running BWA ALN for R2...")
            with open(sai_file_2, "w") as sai_out:
                subprocess.run(aln_cmd_base + [paired_file], stdout=sai_out, check=True)
            
            logger.info(f"[{sample_id}] BWA ALN complete, running sampe...")
            
            # Run bwa sampe and pipe to samtools view
            with open(bam_file, "wb") as bam_out:
                sampe_proc = subprocess.Popen(
                    ["bwa", "sampe", ref, sai_file_1, sai_file_2, input_file, paired_file], 
                    stdout=subprocess.PIPE
                )
                view_proc = subprocess.Popen(["samtools", "view", "-b", "-"], stdin=sampe_proc.stdout, stdout=bam_out)
                sampe_proc.wait()
                view_proc.communicate()
            
            logger.info(f"[{sample_id}] Sorting BAM...")
            
            # Sort BAM
            subprocess.run(["samtools", "sort", "-o", sorted_bam_file, bam_file], check=True)
            logger.info(f"[{sample_id}] Saved {sorted_bam_file}")
            
            # Extract mapped reads to FASTQ (paired-end output)
            # samtools fastq writes directly to files specified by -1 and -2
            subprocess.run(
                ["samtools", "fastq", "-F", "4", "-1", mapped_fastq_1, "-2", mapped_fastq_2, "-0", "/dev/null", "-s", "/dev/null", sorted_bam_file],
                check=True
            )
            
            # Generate separate flagstats for read 1 and read 2
            # For read 1: filter for first in pair (flag 64)
            with open(stats_file_1, "w") as out:
                r1_proc = subprocess.Popen(["samtools", "view", "-b", "-f", "64", sorted_bam_file], 
                                          stdout=subprocess.PIPE)
                subprocess.run(["samtools", "flagstat", "-O", "tsv"], 
                             stdin=r1_proc.stdout, stdout=out, check=True)
                r1_proc.wait()
            
            # For read 2: filter for second in pair (flag 128)
            with open(stats_file_2, "w") as out:
                r2_proc = subprocess.Popen(["samtools", "view", "-b", "-f", "128", sorted_bam_file], 
                                          stdout=subprocess.PIPE)
                subprocess.run(["samtools", "flagstat", "-O", "tsv"], 
                             stdin=r2_proc.stdout, stdout=out, check=True)
                r2_proc.wait()
            
            # Clean up .sai files
            for sai in [sai_file_1, sai_file_2]:
                if os.path.exists(sai):
                    os.remove(sai)
            
            logger.info(f"[{sample_id}] Finished BWA ALN/SAMtools for unmerged reads")
            
        except Exception as e:
            logger.error(f"[{sample_id}] Error processing {input_file} and {paired_file}: {e}")
            raise
    else:
        logger.error(f"[{sample_id}] Invalid file configuration for bwa aln: {input_file}, paired_file={paired_file}")
        raise ValueError(f"Invalid file configuration for bwa aln")


def generate_mapping_summary(sample_ids, output_dir, aligner, logger):
    """Generate a summary CSV with mapping statistics for all samples."""
    
    summary_data = []
    
    for sid in sample_ids:
        merged_flagstats = os.path.join(output_dir, f"{sid}_mapped_flagstats.txt")
        unmerged_1_flagstats = os.path.join(output_dir, f"{sid}_unmerged_1_flagstats.txt")
        unmerged_2_flagstats = os.path.join(output_dir, f"{sid}_unmerged_2_flagstats.txt")
        
        # Parse merged stats
        n_merged, perc_merged = parse_flagstats(merged_flagstats)
        
        # Parse unmerged stats
        n_unmerged_1, _ = parse_flagstats(unmerged_1_flagstats)
        n_unmerged_2, _ = parse_flagstats(unmerged_2_flagstats)
        
        summary_data.append({
            'ID': sid,
            'alignment_method': aligner,
            'mapped_n_reads': n_merged if n_merged is not None else 'NA',
            'mapped_perc': perc_merged if perc_merged is not None else 'NA',
            'mapped_n_reads_unmerged_1': n_unmerged_1 if n_unmerged_1 is not None else 'NA',
            'mapped_n_reads_unmerged_2': n_unmerged_2 if n_unmerged_2 is not None else 'NA'
        })
    
    # Create DataFrame and save to CSV
    df = pd.DataFrame(summary_data)
    summary_file = os.path.join(output_dir, "mapping_summary.csv")
    df.to_csv(summary_file, index=False)
    logger.info(f"Mapping summary saved to {summary_file}")


def main(args):
    
    # Setup
    os.makedirs(args.output_dir, exist_ok=True)
    temp_dir = os.path.join(args.output_dir, "tmp")
    os.makedirs(temp_dir, exist_ok=True)

    # Set up log file with timestamp if not provided
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'mapping_log_{timestamp}.log'
    
    # Set up logging    
    log_dir_path = os.path.join(args.output_dir, "logs")
    logger = setup_logging(log_dir=log_dir_path, log_file=args.log_file)

    logger.info(f"Using aligner: {args.aligner}")
    if args.aligner == "bwa-aln":
        logger.info("BWA ALN parameters: -l 16500 -n 0.01 -o 2 (ancient DNA-specific, Green et al. 2010)")
    else:
        logger.info("Using default BWA MEM to map reads to references.")

    # Step 1: Get sample IDs from tracking sheet
    if not args.tracking_sheet:
        logger.error("You must provide a tracking sheet with sample IDs.")
        sys.exit(1)

    logger.info(f"Reading tracking sheet: {args.tracking_sheet}")
    try:
        sample_ids = load_name_ids(args.tracking_sheet, args.column_name, sheet=args.sheet)
        logger.info(f"Found {len(sample_ids)} sample IDs")
    except Exception as e:
        logger.error(f"Failed to load sample IDs: {e}")
        sys.exit(1)
        

    # Step 2: Gather input files
    jobs = []
    for sid in sample_ids:
        input_matches = glob.glob(os.path.join(args.input_dir, f"{sid}*merged*.fq"))
        ref_matches = glob.glob(os.path.join(args.ref_dir, f"{sid}*_seed.fasta"))
        logger.info(f"[{sid}] Found {len(input_matches)} FASTQ files and {len(ref_matches)} reference files")

        # Check reference FASTAs
        if len(ref_matches) == 0:
            logger.error(f"[{sid}] No reference FASTA found in {args.ref_dir}; skipping")
            continue
        elif len(ref_matches) > 1:
            logger.error(f"[{sid}] Multiple references found: {ref_matches}; skipping")
            continue

        ref_file = ref_matches[0]

    # Check find all input FASTQs
        if len(input_matches) == 0:
            logger.error(f"[{sid}] No FASTQ found in {args.input_dir}; skipping")
            continue
        elif len(input_matches) > 3:
            logger.error(f"[{sid}] Too many FASTQs found: {input_matches}; skipping")
            continue

        for input_file in input_matches:
            jobs.append((sid, input_file, ref_file))

    if not jobs:
        logger.error("No valid input/ref pairs found. Exiting.")
        sys.exit(1)

    logger.info(f"Total jobs to process: {len(jobs)}")

    #Step 3: Index seed files
    logger.info("Indexing reference FASTAs with BWA")
    indexed_refs = set()
    for _, _, ref_file in jobs:
        if ref_file in indexed_refs:
            continue
        # BWA creates files with extensions: .amb, .ann, .bwt, .pac, .sa
        index_prefix = ref_file
        if not all(os.path.exists(f"{index_prefix}.{ext}") for ext in ["amb","ann","bwt","pac","sa"]):
            logger.info(f"Indexing {ref_file}")
            try:
                subprocess.run(["bwa", "index", ref_file], check=True)
            except Exception as e:
                logger.error(f"Failed to index {ref_file}: {e}")
                sys.exit(1)
        else:
            logger.info(f"Index already exists for {ref_file}")
        indexed_refs.add(ref_file)


    # Step 4: Run BWA/SAMtools based on aligner choice
    if args.aligner == "bwa-mem":
        logger.info("Mapping reads to references with BWA MEM and SAMtools")
        max_workers = args.threads or os.cpu_count()
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = [
                executor.submit(run_bwa_mem_and_samtools, sid, infile, args.output_dir, temp_dir, reffile, 8, logger)
                for sid, infile, reffile in jobs
            ]

            for future in as_completed(futures):
                try:
                    future.result()
                except Exception as e:
                    logger.error(f"BWA MEM/SAMtools future failed: {e}")
                    
    elif args.aligner == "bwa-aln":
        logger.info("Mapping reads to references with BWA ALN and SAMtools")
        
        # Group jobs by sample_id to handle paired-end reads together
        jobs_by_sample = {}
        for sid, infile, reffile in jobs:
            if sid not in jobs_by_sample:
                jobs_by_sample[sid] = {"ref": reffile, "files": []}
            jobs_by_sample[sid]["files"].append(infile)
        
        logger.info(f"Processing {len(jobs_by_sample)} samples")
        
        max_workers = args.threads or os.cpu_count()
        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {}
            for sid, data in jobs_by_sample.items():
                reffile = data["ref"]
                files = data["files"]
                
                # Separate merged and unmerged files
                merged_files = [f for f in files if "merged.fq" in f and "unmerged" not in f]
                unmerged_1_files = [f for f in files if "unmerged_1.fq" in f]
                unmerged_2_files = [f for f in files if "unmerged_2.fq" in f]
                
                # Process merged file (single-end)
                for merged_file in merged_files:
                    future = executor.submit(
                        run_bwa_aln_and_samtools, sid, merged_file, args.output_dir, temp_dir, reffile, 8, logger, None
                    )
                    futures[future] = f"{sid}_merged"
                
                # Process paired-end unmerged files together
                if unmerged_1_files and unmerged_2_files:
                    if len(unmerged_1_files) == 1 and len(unmerged_2_files) == 1:
                        future = executor.submit(
                            run_bwa_aln_and_samtools, sid, unmerged_1_files[0], args.output_dir, temp_dir, reffile, 8, logger, unmerged_2_files[0]
                        )
                        futures[future] = f"{sid}_unmerged"
                    else:
                        logger.warning(f"[{sid}] Unexpected number of unmerged files: {len(unmerged_1_files)} _1 files and {len(unmerged_2_files)} _2 files")

            completed = 0
            failed = 0
            for future in as_completed(futures):
                job_name = futures[future]
                try:
                    future.result()
                    completed += 1
                    logger.info(f"Completed {completed}/{len(futures)} jobs ({job_name})")
                except Exception as e:
                    failed += 1
                    logger.error(f"[{job_name}] Job failed: {e}")
            
            logger.info(f"BWA ALN Step 4 complete: {completed} succeeded, {failed} failed")


    # Step 5: Repair pairs (with empty file check)
    logger.info("Re-pairing unmerged read pairs")

    max_workers = args.threads or os.cpu_count()
    repair_count = 0
    skip_count = 0
    
    with ThreadPoolExecutor(max_workers=max_workers) as executor:
        futures = {}
        for sid in sample_ids:
            unmerged_1 = os.path.join(args.output_dir, f"{sid}_mapped_unmerged_1.fastq")
            unmerged_2 = os.path.join(args.output_dir, f"{sid}_mapped_unmerged_2.fastq")
            
            # Check if files exist
            if not os.path.exists(unmerged_1) or not os.path.exists(unmerged_2):
                logger.warning(f"[{sid}] Unmerged pair files not found, skipping repair.")
                skip_count += 1
                continue
            
            # Check if files are empty or too small
            size_1 = os.path.getsize(unmerged_1)
            size_2 = os.path.getsize(unmerged_2)
            
            if size_1 == 0 or size_2 == 0:
                logger.warning(f"[{sid}] Empty unmerged files (R1={size_1} bytes, R2={size_2} bytes), skipping repair.")
                skip_count += 1
                continue
            
            logger.info(f"[{sid}] Submitting repair job (R1={size_1} bytes, R2={size_2} bytes)")
            future = executor.submit(repair_reads, unmerged_1, unmerged_2, args.output_dir, sid, logger)
            futures[future] = sid
            repair_count += 1

        logger.info(f"Submitted {repair_count} repair jobs, skipped {skip_count} samples")

        for future in as_completed(futures):
            sid = futures[future]
            try:
                future.result()
                logger.info(f"[{sid}] Repair completed")
            except Exception as e:
                logger.error(f"[{sid}] Repair failed: {e}")

    # Step 6: Generate mapping summary
    logger.info("Generating mapping summary CSV")
    generate_mapping_summary(sample_ids, args.output_dir, args.aligner, logger)

    # Step 7: Clean up
    pattern = os.path.join(args.output_dir, f"*_mapped_unmerged_*.fastq")

    for f in glob.glob(pattern):
        try:
            os.remove(f)
            logger.info(f"Removed: {f}")
        except OSError as e:
            logger.error(f"Error removing {f}: {e}")

    cleanup_temp_dir(temp_dir)
    logger.info("Compressing and cleaning up data directories")
    logger.info("All samples processed!")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Map merged reads to per-sample references",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input_dir", required=True, help="Directory with FASTQ files")
    parser.add_argument("--ref_dir", required=True, help="Directory with per-sample reference FASTAs")
    parser.add_argument("--output_dir", required=True, help="Directory for outputs")
    parser.add_argument("--tracking_sheet", required=True, help="CSV or XLSX tracking sheet")
    parser.add_argument("--column_name", required=True, help="Column with sample IDs in tracking sheet")
    parser.add_argument("--aligner", required=True, choices=["bwa-mem", "bwa-aln"], help="BWA algorithm to use: bwa-mem or bwa-aln")
    parser.add_argument("--sheet", type=int, default=0, help="Sheet index if XLSX file is used")
    parser.add_argument("--threads", type=int, default=4, help="Number of parallel jobs (not threads per job)")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `mapping_log` file will be created")

    args = parser.parse_args()
    main(args)
