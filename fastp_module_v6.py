import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import os
import sys
import glob
import subprocess
import argparse
import shutil
from pathlib import Path
import pandas as pd
from datetime import datetime
import pathlib

from metahist_tools import clean_and_tar, run_command, xlsx2csv, setup_logging
# add get_read_ids2 to metahist_tools later

def run_fastp_trim(sample_ids, r1_path, r2_path, output_dir, extra_args=None):
    """Run initial fastp trimming and filtering on paired-end reads."""

    tmp_dir = os.path.join(output_dir, "tmp")
    os.makedirs(tmp_dir, exist_ok=True)

    fastp_command = [
        "fastp", "-i", r1_path, "-I", r2_path,
        "--unpaired1", f"{tmp_dir}/{sample_ids}_unpaired_1.fq",
        "--unpaired2", f"{tmp_dir}/{sample_ids}_unpaired_2.fq",
        "--out1", f"{output_dir}/{sample_ids}_trimmed_1.fq",
        "--out2", f"{output_dir}/{sample_ids}_trimmed_2.fq",
        "--detect_adapter_for_pe",
        "--qualified_quality_phred=30", "--trim_poly_g",
        "--correction", "--dedup",
        "--html", f"{output_dir}/{sample_ids}_trim.html",
        "--json", f"{output_dir}/{sample_ids}_trim.json"
    ]

    if extra_args:
        fastp_command += extra_args

    # Redirect output to sample-specific log file in output directory
    log_file = os.path.join(output_dir, f"{sample_ids}_trim.log")
    with open(log_file, 'w', encoding='utf-8') as log:
        subprocess.run(fastp_command, stdout=log, stderr=log, check=True)

def run_fastp_merge(sample_ids, output_dir):
    """Run fastp merge on trimmed paired-end reads
    This is designed to rely on the output from run_fastp_trim.
    Assumes run_fastp_trim naming convention for input pairs."""

    fastp_command = [
        "fastp", "--in1", f"{output_dir}/{sample_ids}_trimmed_1.fq",
        "--in2", f"{output_dir}/{sample_ids}_trimmed_2.fq",
        "--merge", "--merged_out", f"{output_dir}/{sample_ids}_merged.fq",
        "--out1", f"{output_dir}/{sample_ids}_unmerged_1.fq",
        "--out2", f"{output_dir}/{sample_ids}_unmerged_2.fq",
        "--length_required", "30",
        "--html", f"{output_dir}/{sample_ids}_merge.html",
        "--json", f"{output_dir}/{sample_ids}_merge.json"
    ]

    # Redirect output to sample-specific log file in output directory
    log_file = os.path.join(output_dir, f"{sample_ids}_merge.log")
    with open(log_file, 'w', encoding='utf-8') as log:
        subprocess.run(fastp_command, stdout=log, stderr=log, check=True)

def generate_seqkit_stats(output_dir, stats_output):
    """Generate seqkit stats summary of all output files in output_dir"""

    fastq_files = glob.glob(os.path.join(output_dir, "*.f*q"))
    seqkit_command = ["seqkit", "stats", *fastq_files]
    with open(stats_output, "w", encoding="utf-8") as stats_file:
        subprocess.run(seqkit_command, stdout=stats_file, check=True)

def run_fastp_overlap_plot(sample_ids, r1_path, r2_path, output_dir):
    """Uses fastp merge without merge output to generate overlap plots as an html file."""

    command = [
        "fastp", "--in1", r1_path, "--in2", r2_path,
        "--stdout", "--merge", "-A", "-G", "-Q", "-L",
        "--json", "/dev/null",
        "--html", os.path.join(output_dir, f"{sample_ids}_overlaps.html")
    ]
    
    # Redirect output to sample-specific log file in output directory
    log_file = os.path.join(output_dir, f"{sample_ids}_overlaps.log")
    with open(log_file, 'w', encoding="utf-8") as log:
        subprocess.run(command, stdout=log, stderr=log)

def run_fastp_json_summary(
    json_dir,
    r_script,
    logger,
    output_file="combined_fastp_json_out.csv"
):
    """Run the R fastp JSON parser on all JSON files in a directory."""

    json_dir = Path(json_dir)
    r_script = Path(r_script)

    # Confirm script exists
    if not r_script.exists():
        logger.error(f"R script not found: {r_script}")
        raise FileNotFoundError(f"R script not found: {r_script}")

    # Collect input JSONs
    json_files = sorted(json_dir.glob("*_trim.json"))
    if not json_files:
        logger.error(f"No JSON files found in: {json_dir}")
        raise RuntimeError(f"No JSON files found in: {json_dir}")

    cmd = ["Rscript", str(r_script)] + [str(f) for f in json_files]

    logger.info(f"Running fastp JSON parser: {' '.join(cmd)}")

    subprocess.run(cmd, check=True)

    logger.info(f"Combined JSON output written to: {output_file}")



def process_sample(sample_id, r1, r2, output_dir, logger, extra_args=None):
    """The function for order of fastp operations per sample."""

    logger.info(f"Processing {sample_id}")
    run_fastp_trim(sample_id, r1, r2, output_dir, extra_args)
    run_fastp_merge(sample_id, output_dir)
    run_fastp_overlap_plot(sample_id, r1, r2, output_dir)
    logger.info(f"Finished processing {sample_id}")

def find_paths_case_insensitive(df, target_names):
    """
    Find read path columns in the DataFrame by checking multiple possible names (case-insensitve)
    """
    columns_lower = {col.lower(): col for col in df.columns}
    for target in target_names:
        if target.lower() in columns_lower:
            return columns_lower[target.lower()]
    return None

def get_read_ids2(input_source, logger, paired=False, prefix=None, suffix=None, column_name=None, sheet=0):
    """
    Function to extract read IDs from a directory or spreadsheet.

    input_source : str
        Path to a directory or CSV/XLSX file.
    paired : bool, default=True
        If True, return only paired IDs (occurring twice); otherwise, return singletons.
    prefix : str, optional
        Optional filename prefix to match.
    suffix : str, optional
        Optional filename suffix to match before file extension.
    column_name : str, optional
        Required for spreadsheet input to specify the column of IDs.
    sheet : int or str, default=0
        Sheet index for Excel files.
    """

    input_path = pathlib.Path(input_source)

    # Spreadsheet input
    if input_path.is_file():
        if input_source.endswith(".xlsx"):
            xlsx_read = xlsx2csv(input_source, sheet=sheet)
            df = pd.read_csv(xlsx_read)
            logger.info(f"Using sheet: {xlsx_read}")
        else:
            df = pd.read_csv(input_source)

        if prefix:
            logger.error("Ignoring prefix due to spreadsheet input.")

        if column_name is None or column_name not in df.columns:
            logger.error("You must provide a valid column_name for spreadsheet input.")
        ids_to_include = sorted(df[column_name].dropna().astype(str).tolist())
        return ids_to_include

    # Directory input
    if input_path.is_dir():
        all_files = []
        for root, _, files in os.walk(input_source):
            for fname in files:
                all_files.append(os.path.join(root, fname))

        logger.info(f"Total files found: {len(all_files)}")

        # Filter all_files to keep FASTQ and based on prefix/suffix
        matched_files = []
        matched_ids = []
        names = [os.path.basename(f) for f in all_files]
        all_files_ids = [re.sub(r"(?:_[Rr]?[12])?\.f(?:ast)?q(?:\.gz)?$", "", n) for n in names]
        logger.info(f"File names before applying prefix/suffix filter: {len(names)}")
        if prefix or suffix:
            # Escape prefix/suffix for regex if present
            pre_pattern = prefix if prefix else ""
            suf_pattern = suffix if suffix else ""
            pattern = rf"({pre_pattern}.+?{suf_pattern})(?:_[Rr]?[12])?\.f(?:ast)?q(?:\.gz)?$"
            my_list = [re.match(pattern, n) is not None for n in names]
            logger.info(f"Files matching prefix/suffix pattern: {sum(my_list)}")
            if sum(my_list) == 0:
                logger.error("No files matched the given prefix/suffix pattern.")

            for i, match in enumerate(my_list):
                if match:
                    matched_files.append(names[i])
                    matched_ids.append(all_files_ids[i])
        else:
            matched_files = names
            matched_ids = all_files_ids
        logger.info(f"Files after applying prefix/suffix filter: {len(matched_files)}")
        if len(matched_files) == 0:
            raise ValueError("No files matched the given prefix/suffix criteria.")
        logger.info(f"Sample matched files: {matched_files} ...")
        logger.info(f"Sample matched IDs: {matched_ids} ...")
        logger.info(f"Total matched files: {len(matched_files)}")
        logger.info(f"Total matched IDs: {len(matched_ids)}")
        if len(matched_files) != len(matched_ids):
            raise ValueError("Mismatch between matched files and IDs.")
        logger.info("Counting occurrences of matched filenames...")

        # Count occurrences of full matched filenames (including extension)
        id_counts = Counter(matched_ids)
        # Build sets of paired and single files
        paired_ids = set(id_ for id_, c in id_counts.items() if c == 2)
        single_ids = set(id_ for id_, c in id_counts.items() if c == 1)
        for id_, c in id_counts.items():
            if c > 2:
                logger.info(f"Warning: File '{id_}' appears {c} times. Consider adding a suffix to filter duplicates.")

        # Return match_ids for paired or single files
        if paired:
            ids_to_include = paired_ids
        else:
            ids_to_include = single_ids
        logger.info(f"Final IDs: {ids_to_include} ...")
        return ids_to_include

    raise FileNotFoundError(f"'{input_source}' is neither a valid file nor directory")

def main(args):
    """The main function to handle argument parsing and workflow of module."""

    # Set up
    os.makedirs(args.output_dir, exist_ok=True)
    
    # Create logs directory inside output_dir instead of current directory
    log_dir = os.path.join(args.output_dir, "logs")
    os.makedirs(log_dir, exist_ok=True)

    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'fastp_log_{timestamp}.log'
    
    logger = setup_logging(log_dir=log_dir, log_file=args.log_file)

    jobs = []
    paired_files = {}

    # Begin processing
    if args.tracking_sheet:
        logger.info(f"Using tracking sheet: {args.tracking_sheet}")
        
        # Read the CSV to get file paths directly
        if args.tracking_sheet.endswith(".xlsx"):
            xlsx_read = xlsx2csv(args.tracking_sheet, sheet=args.sheet)
            df = pd.read_csv(xlsx_read)
            logger.info(f"Using sheet: {xlsx_read}")
        else:
            df = pd.read_csv(args.tracking_sheet)
        
        # Check if column_name is valid
        if args.column_name is None or args.column_name not in df.columns:
            logger.error(f"You must provide a valid column_name. Available columns: {df.columns.tolist()}")
            sys.exit(1)
        
        # Filter by specific IDs if provided
        if args.ids:
            ids_list = args.ids.split(",")
            df = df[df[args.column_name].astype(str).isin(ids_list)]
            logger.info(f"Filtered to {len(df)} samples based on --ids argument")
        
        # Check for forward/fwd and reverse/rev columns (case-insensitive)
        forward_col = find_paths_case_insensitive(df, ['forward', 'fwd'])
        reverse_col = find_paths_case_insensitive(df, ['reverse', 'rev'])
        
        if forward_col and reverse_col:
            # Use the file paths directly from CSV
            logger.info(f"Using '{forward_col}' and '{reverse_col}' columns from tracking sheet for file paths")
            
            for _, row in df.iterrows():
                sample_id = str(row[args.column_name])
                r1_path = str(row[forward_col])
                r2_path = str(row[reverse_col])
                
                # Check if files exist
                if os.path.exists(r1_path) and os.path.exists(r2_path):
                    logger.info(f"Found files for {sample_id}: {r1_path} and {r2_path}")
                    paired_files[sample_id] = [r1_path, r2_path]
                    jobs.append((sample_id, r1_path, r2_path))
                else:
                    if not os.path.exists(r1_path):
                        logger.error(f"Error: Forward file does not exist for ID '{sample_id}': {r1_path}")
                    if not os.path.exists(r2_path):
                        logger.error(f"Error: Reverse file does not exist for ID '{sample_id}': {r2_path}")
        else:
            # Fall back to the existing glob search method
            logger.info(f"Forward/reverse columns not found in tracking sheet (checked: forward, fwd, reverse, rev)")
            logger.info("Falling back to searching for files in input directory using sample IDs")
            
            if not args.input_dir:
                logger.error("--input_dir is required when tracking sheet doesn't have forward/reverse columns")
                sys.exit(1)
            
            # Get IDs from tracking sheet
            ids = df[args.column_name].dropna().astype(str).tolist()
            logger.info(f"Extracted {len(ids)} sample IDs from tracking sheet")
            
            # Search for files in input directory using the existing matching logic
            for sample_id in ids:
                id_pair = glob.glob(os.path.join(args.input_dir, "**", f"*{sample_id}*"), recursive=True)
                # Filter to only fastq files
                id_pair = [f for f in id_pair if re.search(r'\.f(ast)?q(\.gz)?$', f, re.IGNORECASE)]
                
                if len(id_pair) == 2:
                    r1_path, r2_path = sorted(id_pair)
                    logger.info(f"Found files for {sample_id}: {r1_path} and {r2_path}")
                    paired_files[sample_id] = [r1_path, r2_path]
                    jobs.append((sample_id, r1_path, r2_path))
                else:
                    logger.error(f"Error: ID '{sample_id}' has {len(id_pair)} matching files in {args.input_dir}")
    
    else:
        # Original behavior: use input directory with get_read_ids2
        logger.info(f"Using input directory: {args.input_dir}")
        
        if args.ids is None:
            ids = get_read_ids2(
                args.input_dir,
                logger=logger,
                paired=True,
                prefix=args.prefix,
                suffix=args.suffix,
                column_name=args.column_name,
                sheet=args.sheet
            )
        else:
            ids = args.ids.split(",")

        # Search for files in input directory
        for sample_id in ids:
            id_pair = glob.glob(os.path.join(args.input_dir, "**", f"{sample_id}*"), recursive=True)
            if len(id_pair) == 2:
                r1_path, r2_path = sorted(id_pair)
                logger.info(f"Running files: {r1_path} and {r2_path}")
                paired_files[sample_id] = [r1_path, r2_path]
                jobs.append((sample_id, r1_path, r2_path))
            else:
                logger.error(f"Error: ID '{sample_id}' has {len(id_pair)} matching files.")

    # Check if we have any jobs to run
    if not jobs:
        logger.error("No valid file pairs found. Exiting.")
        sys.exit(1)
    
    logger.info(f"Total samples to process: {len(jobs)}")

    # Run in parallel after jobs are built
    with ThreadPoolExecutor(max_workers=args.threads) as executor:
        futures = [
            executor.submit(
                process_sample,
                sample_id,
                r1,
                r2,
                args.output_dir,
                logger,
                args.fastp_extra_args
            )
            for sample_id, r1, r2 in jobs
        ]
        for future in futures:
            future.result()

    # Generate summary stats and convert to CSV
    logger.info("Generating summary statistics...")
    generate_seqkit_stats(args.output_dir, f"{args.output_dir}/fastp_summary.out")
    generate_seqkit_stats(f"{args.output_dir}/tmp", f"{args.output_dir}/tmp/fastp_summary_intermediate.out")
    df = pd.read_csv(f'{args.output_dir}/fastp_summary.out', sep=r'\s+', header=None)
    df.to_csv(f'{args.output_dir}/fastp_summary.csv', header=None, index=False)
    
    # Generate summary stats from JSON files
    pattern = os.path.join(args.output_dir, "*_trim.json")
    json_paths = glob.glob(pattern, recursive=False)
    if not json_paths:
        raise FileNotFoundError(f"No *_trim.json files found under {args.output_dir}")
    run_fastp_json_summary(
        os.path.dirname(json_paths[0]),
        r_script="/mnt/shared/scratch/museomix/its-fun-2-map/parse_fastp_json.R",
        logger=logger,
        output_file=os.path.join(args.output_dir, "combined_fastp_json_out.csv")
    )
    
    # Remove tmp directory if it exists and finish
    tmp_dir = os.path.join(args.output_dir, "tmp")
    if os.path.exists(tmp_dir):
        shutil.rmtree(tmp_dir)
        logger.info(f"Removed directory: {tmp_dir}")
    logger.info("All samples processed!")


if __name__ == "__main__":
    if len(sys.argv) == 1:
        parser = argparse.ArgumentParser()
        parser.print_help()
        sys.exit(1)

    parser = argparse.ArgumentParser(description="Pre-process and process raw read data using fastp.")
    parser.add_argument("--input_dir", type=str, required=False, help="Path to input directory containing FASTQ files. Required when tracking sheet doesn't have forward/reverse columns, or when not using a tracking sheet.")
    parser.add_argument("--output_dir", type=str, required=True, help="Path to output directory.")
    parser.add_argument("--prefix", required=False, help="Prefix to use to find specific files in the input directory.")
    parser.add_argument("--ids", required=False, help="List of specific file IDs to use. Comma-separated. If not provided, all IDs will be used.")
    parser.add_argument("--tracking_sheet", required=False, help="CSV/XLSX file with tracking metadata. If it has 'forward'/'fwd' and 'reverse'/'rev' columns (case-insensitive), file paths will be read directly. Otherwise, IDs will be extracted and matched against input_dir.")
    parser.add_argument("--column_name", required=False, help="Column name in tracking sheet with library names/IDs. Required when using --tracking_sheet.")
    parser.add_argument("--sheet", type=int, default=0, help="(For XLSX input only) Optional sheet index to be used as a tracking sheet. Default is the first sheet (--sheet 0).")
    parser.add_argument("--suffix", required=False, help="Suffix to use to find specific files in the input directory.")
    parser.add_argument("--threads", type=int, default=4, help="Number of threads for parallel processing.")
    parser.add_argument("--fastp_extra_args", nargs=argparse.REMAINDER, help="Extra arguments for fastp to happen before merging.")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `fastp_log` file will be created")

    args = parser.parse_args()
    
    # Validate argument combinations
    if args.tracking_sheet and not args.column_name:
        parser.error("--column_name is required when using --tracking_sheet")
    
    main(args)

##Example usage:
#python fastp_module.py --input_dir raw_reads --output_dir fastp_processed --fastp_extra_args --trim_front1 10 --trim_front2 10
#python fastp_module.py --tracking_sheet samples.csv --column_name ID --output_dir fastp_processed
#python fastp_module.py --tracking_sheet samples.csv --column_name ID --input_dir raw_reads --output_dir fastp_processed
