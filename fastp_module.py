import re
from collections import Counter
from concurrent.futures import ThreadPoolExecutor
import os
import sys
import glob
import shlex
import subprocess
import argparse
import shutil
from pathlib import Path
import pandas as pd
from datetime import datetime
import pathlib
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, default_converter
from seqpy_tools import clean_and_tar, run_command, xlsx2csv, setup_logging

# fastp flags that this module sets itself. Passing any of these via
# --fastp_extra_args would place a duplicate on the fastp command line, so they
# are rejected at startup in favour of the dedicated CLI arguments.
MANAGED_FASTP_FLAGS = {
    '-q': '--qualified_quality_phred',
    '--qualified_quality_phred': '--qualified_quality_phred',
    '-w': '--fastp_threads',
    '--thread': '--fastp_threads',
}


def parse_fastp_extra_args(extra_args_string):
    """Split a quoted --fastp_extra_args string into a list of fastp arguments.

    Returns an empty list when no extras were supplied.
    """

    if not extra_args_string:
        return []

    return shlex.split(extra_args_string)


def check_fastp_extra_args(extra_args, parser=None, logger=None):
    """Reject extra fastp arguments that duplicate flags this module manages."""

    conflicts = []

    for token in extra_args:
        # Handle both "--flag value" and "--flag=value" forms
        flag = token.split('=', 1)[0]
        if flag in MANAGED_FASTP_FLAGS:
            conflicts.append((flag, MANAGED_FASTP_FLAGS[flag]))

    if conflicts:
        details = "; ".join(
            f"'{flag}' is set by this script - use {replacement} instead"
            for flag, replacement in conflicts
        )
        message = f"Conflicting fastp arguments in --fastp_extra_args: {details}"

        if parser is not None:
            parser.error(message)
        if logger is not None:
            logger.error(message)
        raise ValueError(message)

    return True


def run_fastp_trim(
    sample_ids,
    r1_path,
    r2_path,
    output_dir,
    qualified_quality_phred=30,
    fastp_threads=3,
    extra_args=None
):
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
        f"--qualified_quality_phred={qualified_quality_phred}", "--trim_poly_g",
        "--correction", "--dedup",
        "--thread", str(fastp_threads),
        "--html", f"{output_dir}/{sample_ids}_trim.html",
        "--json", f"{output_dir}/{sample_ids}_trim.json"
    ]

    if extra_args:
        fastp_command += extra_args

    # Redirect output to sample-specific log file in output directory
    log_file = os.path.join(output_dir, f"{sample_ids}_trim.log")
    with open(log_file, 'w', encoding='utf-8') as log:
        subprocess.run(fastp_command, stdout=log, stderr=log, check=True)

def run_fastp_merge(sample_ids, output_dir, fastp_threads=3):
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
        "--thread", str(fastp_threads),
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

def run_fastp_overlap_plot(sample_ids, r1_path, r2_path, output_dir, fastp_threads=3):
    """Uses fastp merge without merge output to generate overlap plots as an html file."""

    command = [
        "fastp", "--in1", r1_path, "--in2", r2_path,
        "--stdout", "--merge", "-A", "-G", "-Q", "-L",
        "--thread", str(fastp_threads),
        "--json", "/dev/null",
        "--html", os.path.join(output_dir, f"{sample_ids}_overlaps.html")
    ]
    
    # Redirect output to sample-specific log file in output directory
    log_file = os.path.join(output_dir, f"{sample_ids}_overlaps.log")
    with open(log_file, 'w', encoding="utf-8") as log:
        subprocess.run(command, stdout=log, stderr=log)

def run_fastp_json_summary(
    json_dir,
    outdir,
    logger
):
    """Run the R fastp JSON parser on all JSON files in a directory."""

    json_dir = Path(json_dir)
    r_script = Path(__file__).parent / 'parse_fastp_json.R'

    # Confirm script exists
    if not r_script.exists():
        logger.error("R script not found: %s", r_script)
        raise FileNotFoundError(f"R script not found: {r_script}")

    logger.info(f"Running fastp JSON parser")

    # Assess for which sequences/contigs can be taken forward for submission
    r['source'](str(r_script))
    json2csv = r['json2csv']

    with localconverter(default_converter + pandas2ri.converter):
        result = json2csv(str(json_dir), str(outdir))

    logger.info("Combined JSON output written to: %s", outdir)



def process_sample(
    sample_id,
    r1,
    r2,
    output_dir,
    logger,
    qualified_quality_phred=30,
    fastp_threads=3,
    extra_args=None
):
    """The function for order of fastp operations per sample."""

    logger.info(f"Processing {sample_id}")
    run_fastp_trim(
        sample_id,
        r1,
        r2,
        output_dir,
        qualified_quality_phred=qualified_quality_phred,
        fastp_threads=fastp_threads,
        extra_args=extra_args
    )
    run_fastp_merge(sample_id, output_dir, fastp_threads=fastp_threads)
    run_fastp_overlap_plot(sample_id, r1, r2, output_dir, fastp_threads=fastp_threads)
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
    
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'fastp_log_{timestamp}.log'
    
    logger = setup_logging(log_file=args.log_file)

    # Parse and validate any extra fastp arguments
    fastp_extra_args = parse_fastp_extra_args(args.fastp_extra_args)
    check_fastp_extra_args(fastp_extra_args, logger=logger)
    if fastp_extra_args:
        logger.info(f"Extra fastp arguments: {' '.join(fastp_extra_args)}")

    # Report effective parallelism: these multiply, so make the product explicit
    logger.info(
        "Parallelism: %d concurrent samples x %d fastp threads = up to %d worker threads",
        args.threads,
        args.fastp_threads,
        args.threads * args.fastp_threads
    )
    logger.info("Quality filtering: --qualified_quality_phred %d", args.qualified_quality_phred)

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
                args.qualified_quality_phred,
                args.fastp_threads,
                fastp_extra_args
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
    else: 
        run_fastp_json_summary(
            json_dir=args.output_dir,
            outdir=args.output_dir,
            logger=logger
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
    parser.add_argument("--threads", type=int, default=4, help="Number of samples to process in parallel (i.e. concurrent fastp processes). This is NOT the number of threads used by each fastp job - see --fastp_threads. Total worker threads is approximately --threads x --fastp_threads. Default: 4.")
    parser.add_argument("--fastp_threads", type=int, default=3, help="Number of worker threads used within each fastp process (fastp's -w/--thread). Default: 3, which matches fastp's own default. Set to 1 to make --threads an accurate cap on total CPU use.")
    parser.add_argument("--qualified_quality_phred", type=int, default=30, help="Phred quality score below which a base is counted as unqualified during trimming (fastp's -q/--qualified_quality_phred). Default: 30. Note fastp's own default is 15.")
    parser.add_argument("--fastp_extra_args", type=str, required=False, help="Extra arguments passed to the trimming fastp call, given as a single quoted string, e.g. --fastp_extra_args \"--trim_front1 10 --trim_front2 10\". Cannot be used to set -q/--qualified_quality_phred or -w/--thread; use the dedicated arguments instead.")
    parser.add_argument("--log_file", help="Log file path; if not provided, a timestamped `fastp_log` file will be created")

    args = parser.parse_args()
    
    # Validate argument combinations
    if args.tracking_sheet and not args.column_name:
        parser.error("--column_name is required when using --tracking_sheet")

    if args.threads < 1:
        parser.error("--threads must be at least 1")

    if args.fastp_threads < 1:
        parser.error("--fastp_threads must be at least 1")

    # Fail fast on conflicting fastp extras, before any logging or processing
    check_fastp_extra_args(parse_fastp_extra_args(args.fastp_extra_args), parser=parser)

    main(args)

##Example usage:
#python fastp_module.py --input_dir raw_reads --output_dir fastp_processed --fastp_extra_args "--trim_front1 10 --trim_front2 10"
#python fastp_module.py --tracking_sheet samples.csv --column_name ID --output_dir fastp_processed
#python fastp_module.py --tracking_sheet samples.csv --column_name ID --input_dir raw_reads --output_dir fastp_processed
#python fastp_module.py --input_dir raw_reads --output_dir fastp_processed --threads 8 --fastp_threads 1 --qualified_quality_phred 20
