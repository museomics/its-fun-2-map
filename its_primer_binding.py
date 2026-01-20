#!/usr/bin/env python3

"""
Primer binding script, employing seqkit amplicon to extract multiple target regions input FASTA. 
Default usage is for ITS region extraction using 4 primers from White et al (1990),
able to extract the full ITS region (ITS1->ITS4), the ITS1 region (ITS1->ITS2), and the
ITS2 region (ITS3->ITS4). ITS3 primer is a reverse compliment of ITS2 (see White et al., 1990). 

Custom primer pairs and target regions can be specified via TSV files.

Author: D. Parsons & M. KAMOUYIAROS (@ NHMUK)
"""

import os
import sys
import subprocess
import argparse
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

from metahist_tools import setup_logging

# ITS primer sequences from White et al. 1990
DEFAULT_PRIMERS = {
    'ITS1': 'TCCGTAGGTGAACCTGCGG',
    'ITS2': 'GCTGCGTTCTTCATCGATGC',
    'ITS3': 'GCATCGATGAAGAACGCAGC',
    'ITS4': 'TCCTCCGCTTATTGATATGC'
}

# Primer combinations for each region
DEFAULT_REGIONS = {
    'ITS1': ('ITS1', 'ITS2'),
    'ITS2': ('ITS3', 'ITS4'),
    'ITS_complete': ('ITS1', 'ITS4')
}


def load_primers_tsv(tsv_path):
    """ Load custom primer pairs from a TSV file (Defaults are currently ITS1, ITS2, 
    ITS3, ITS4 from White et al., 1990).
    TSV should have 2 columns: name and sequence """

    primers = {}

    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"name", "sequence"}
        if not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Primer TSV must contain columns: {required}"
            )

        for row in reader:
            name = row["name"].strip()
            seq = row["sequence"].strip().upper()

            if not name or not seq:
                raise ValueError(f"Invalid primer row: {row}")

            primers[name] = seq

    return primers


def load_regions_tsv(tsv_path, primers):
    """ Load region names that are covered by the different primer pairs (Defaults are currently 
    'ITS1': ('ITS1', 'ITS2') 'ITS2': ('ITS3', 'ITS4') 'ITS_complete': ('ITS1', 'ITS4')) 
    TSV should have 3 columns: region, forward, reverse """

    regions = {}

    with open(tsv_path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")

        required = {"region", "forward", "reverse"}
        if not required.issubset(reader.fieldnames):
            raise ValueError(
                f"Regions TSV must contain columns: {required}"
            )

        for row in reader:
            region = row["region"].strip()
            fwd = row["forward"].strip()
            rev = row["reverse"].strip()

            if fwd not in primers or rev not in primers:
                raise ValueError(
                    f"Region '{region}' references undefined primer(s): "
                    f"{fwd}, {rev}"
                )

            regions[region] = (fwd, rev)

    return regions


def load_primers_and_regions(
    primers_tsv=None,
    regions_tsv=None,
    logger=None
):
    """ Load primers and regions from TSV files or use defaults """    

    # Enforce both-or-neither rule
    if (primers_tsv and not regions_tsv) or (regions_tsv and not primers_tsv):
        raise ValueError(
            "Both --primers_tsv and --regions_tsv must be provided together, "
            "or neither (to use built-in defaults)."
        )

    primers = DEFAULT_PRIMERS.copy()
    regions = DEFAULT_REGIONS.copy()

    # Override defaults only if both TSVs are provided
    if primers_tsv and regions_tsv:
        if logger:
            logger.info(
                "Loading primers from %s and regions from %s",
                primers_tsv,
                regions_tsv
            )

        primers = load_primers_tsv(primers_tsv)
        regions = load_regions_tsv(regions_tsv, primers)

    return primers, regions


def read_tracking_sheet(tracking_sheet_path, column_name='ID', logger=None):
    """Read tracking sheet and return set of sample IDs"""
    sample_ids = set()

    try:
        with open(tracking_sheet_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)

            if column_name not in reader.fieldnames:
                logger.error("ERROR: Tracking sheet must contain %s column", column_name)
                sys.exit(1)

            for row in reader:
                sample_id = row[column_name].strip()
                if sample_id:
                    sample_ids.add(sample_id)

        return sample_ids

    except FileNotFoundError:
        logger.error("ERROR: Tracking sheet not found: %s", tracking_sheet_path)
        sys.exit(1)
    except Exception as e:
        logger.error("ERROR: Failed to read tracking sheet: %s", e)
        sys.exit(1)


def extract_sample_name(filename):
    """Extract sample name from filename (everything before first underscore)"""
    base = Path(filename).stem
    return base.split('_')[0]


def find_fasta_for_sample(sample_id, input_dir):
    """Find FASTA file matching the sample ID"""
    # Look for files starting with the sample ID
    fasta_files = list(input_dir.rglob(f'{sample_id}*.fa')) + list(input_dir.rglob(f'{sample_id}*.fasta'))

    if not fasta_files:
        return None

    # If multiple files found, prefer exact match or first file
    for fasta_file in fasta_files:
        if extract_sample_name(fasta_file.name) == sample_id:
            return fasta_file

    return fasta_files[0]


def count_sequences(fasta_file):
    """Count number of sequences in a FASTA file"""
    if not os.path.exists(fasta_file):
        return 0

    # Check if file is empty
    if os.path.getsize(fasta_file) == 0:
        return 0

    try:
        with open(fasta_file, 'r') as f:
            return sum(1 for line in f if line.startswith('>'))
    except:
        return 0


def run_seqkit_amplicon(input_file, output_file, forward_primer, reverse_primer, logger = None):
    """Run seqkit amplicon for a specific primer pair"""
    cmd = [
        'seqkit', 'amplicon',
        '--forward', forward_primer,
        '--reverse', reverse_primer,
        '--max-mismatch', '2',
        '--strict-mode',
        '--threads', '4',
        input_file
    ]

    try:
        with open(output_file, 'w') as out:
            subprocess.run(cmd, stdout=out, stderr=subprocess.PIPE, check=True)
        return True
    except subprocess.CalledProcessError as e:
        logger.error("Error running seqkit amplicon: %s", e)
        return False


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='ITS primer binding script using seqkit amplicon',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""

Basic usage:
  python its_primer_binding.py --input /path/to/fastas --output /path/to/results --tracking_sheet samples.csv
  python its_primer_binding.py -i input_dir/ -o output_dir/ -t tracking.csv

Advanced usage with custom primers and regions:
    python its_primer_binding.py --input /path/to/fastas --output /path/to/results --primers_tsv primers.tsv --regions_tsv regions.tsv --tracking_sheet samples.csv
        """
    )

    parser.add_argument(
        '-i', '--input',
        required=True,
        type=Path,
        help='Input directory containing FASTA files (.fa or .fasta)'
    )

    parser.add_argument(
        '-o', '--output',
        required=True,
        type=Path,
        help='Output directory for results'
    )

    parser.add_argument(
        '-t', '--tracking_sheet',
        required=True,
        type=Path,
        help='CSV file containing sample IDs in "ID" column'
    )

    parser.add_argument(
        '-c', '--column_name',
        default='ID',
        help='Column name in tracking sheet for sample IDs (default: ID)'
    )

    parser.add_argument(
        '-p', '--primers_tsv',
        help='TSV file containing primer information (default: ITS White et al., 1990 primers). \
            TSV should have columns: name, sequence (must be provided if --regions_tsv is used)'
    )

    parser.add_argument(
        '-r', '--regions_tsv',
        help='TSV file containing region information (default: use built-in ITS regions) \
            TSV should have columns: name, forward_primer, reverse_primer \
            (must be provided if --primers_tsv is used)'
    )

    parser.add_argument(
        '--log_file',
        help='Log file name (default: its_primer_binding_<TIMESTAMP>.log)'
    )
    args = parser.parse_args()

    # Validate input directory exists
    if not args.input.exists():
        parser.error(f"Input directory does not exist: {args.input}")

    if not args.input.is_dir():
        parser.error(f"Input path is not a directory: {args.input}")

    # Validate tracking sheet exists
    if not args.tracking_sheet.exists():
        parser.error(f"Tracking sheet does not exist: {args.tracking_sheet}")

    # Validate primers and regions TSVs
    if (args.primers_tsv and not args.regions_tsv) or (args.regions_tsv and not args.primers_tsv):
        parser.error("Both --primers_tsv and --regions_tsv must be provided together, or neither.")
    if args.primers_tsv is not None and not args.primers_tsv.exists():
        parser.error(f"Primers TSV file does not exist: {args.primers_tsv}")
    if args.regions_tsv is not None and not args.regions_tsv.exists():
        parser.error(f"Regions TSV file does not exist: {args.regions_tsv}")


    return args


def main():
    args = parse_args()

    input_dir = args.input
    output_dir = args.output
    tracking_sheet = args.tracking_sheet
    column_name = args.column_name if hasattr(args, 'column_name') else 'ID'

    # Create logs directory inside output_dir and set up logging
    os.makedirs(output_dir, exist_ok=True)

    # Handle log file naming
    if args.log_file is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log_file = f'its_primer_binding_log_{timestamp}.log'

    logger = setup_logging(log_file=args.log_file)

    # Read sample IDs from tracking sheet
    logger.info("Reading sample IDs from tracking sheet: %s", tracking_sheet)
    all_samples = read_tracking_sheet(tracking_sheet, column_name=column_name)

    if not all_samples:
        logger.error("No sample IDs found in tracking sheet")
        sys.exit(1)

    logger.info("Found %d samples in tracking sheet", len(all_samples))
    logger.info("Samples: %s", sorted(all_samples))

    # Create output directories

    PRIMERS, REGIONS  = load_primers_and_regions(primers_tsv=args.primers_tsv,
        regions_tsv=args.regions_tsv, logger=logger)

    for region in REGIONS.keys():
        (output_dir / region).mkdir(parents=True, exist_ok=True)

    logger.info("Starting ITS primer alignment at %s", datetime.now())
    logger.info("Processing files from directory: %s", input_dir)

    # Track which samples have FASTA files
    samples_with_fasta = set()
    samples_without_fasta = set()

    # Process each sample from tracking sheet
    for sample_id in sorted(all_samples):
        # Find FASTA file for this sample
        fasta_file = find_fasta_for_sample(sample_id, input_dir)

        if fasta_file is None:
            print(f"WARNING: No FASTA file found for sample: {sample_id}")
            samples_without_fasta.add(sample_id)
            continue

        samples_with_fasta.add(sample_id)
        print(f"\nProcessing sample {sample_id}: {fasta_file}")

        # Run amplicon extraction for each region
        for region_name, (forward_key, reverse_key) in REGIONS.items():
            output_file = output_dir / region_name / f"{sample_id}_{region_name}.fa"

            forward_primer = PRIMERS[forward_key]
            reverse_primer = PRIMERS[reverse_key]

            run_seqkit_amplicon(
                str(fasta_file),
                str(output_file),
                forward_primer,
                reverse_primer
            )

            count = count_sequences(output_file)
            print(f"  {region_name}: {count} sequences")

    # Report samples without FASTA files
    if samples_without_fasta:
        logger.warning("%d samples from tracking sheet have no FASTA files:", len(samples_without_fasta))
        for sample in sorted(samples_without_fasta):
            logger.warning("  %s", sample)

    # Generate summary statistics
    logger.info("\nGenerating summary report...")

    results = defaultdict(lambda: defaultdict(int))

    # Count sequences for each sample and region
    for sample_id in all_samples:
        for region_name in REGIONS.keys():
            output_file = output_dir / region_name / f"{sample_id}_{region_name}.fa"
            count = count_sequences(output_file)
            results[sample_id][region_name] = count

    # Write CSV summary
    csv_file = output_dir / "extraction_summary.csv"
    with open(csv_file, 'w', newline='') as csvfile:
        fieldnames = [
            'ID',
            'ITS1',
            'ITS1_path',
            'ITS2',
            'ITS2_path',
            'ITS_complete',
            'ITS_complete_path'
        ]
        writer = csv.DictWriter(csvfile, fieldnames=fieldnames)
        writer.writeheader()

        for sample_id in sorted(all_samples):
            row = {'ID': sample_id}

            for region_name in REGIONS.keys():
                output_file = output_dir / region_name / f"{sample_id}_{region_name}.fa"
                count = results[sample_id][region_name]

                # Determine if extracted (YES if count > 0, NO otherwise)
                extracted = 'YES' if count > 0 else 'NO'

                # Get absolute path only if file has sequences
                file_path = str(output_file.resolve()) if count > 0 else ''

                # Add to row based on region
                row[f'{region_name}'] = extracted  # Instead of f'{region_name}_extracted'
                row[f'{region_name}_path'] = file_path

            writer.writerow(row)

    logger.info("CSV summary written to: %s", csv_file)

    # Write summary report
    report_file = output_dir / "summary_report.txt"
    with open(report_file, 'w') as report:
        report.write("ITS Primer Alignment Summary\n")
        report.write(f"Generated on: {datetime.now()}\n")
        report.write(f"Total samples in tracking sheet: {len(all_samples)}\n")
        report.write(f"Samples with FASTA files: {len(samples_with_fasta)}\n")
        report.write(f"Samples without FASTA files: {len(samples_without_fasta)}\n")
        report.write("=" * 70 + "\n\n")

        # Total sequences per region
        for region_name in REGIONS.keys():
            total = sum(results[s][region_name] for s in all_samples)
            report.write(f"{region_name} region: Total sequences: {total}\n")

        report.write("\n" + "=" * 70 + "\n")
        report.write("SAMPLE CATEGORISATION\n")
        report.write("=" * 70 + "\n\n")

        # Per-sample results for each region
        for region_name in REGIONS.keys():
            samples_with_seqs = [s for s in sorted(all_samples) if results[s][region_name] > 0]
            report.write(f"{region_name} sequences ({len(samples_with_seqs)}):\n")
            for sample in samples_with_seqs:
                count = results[sample][region_name]
                if region_name == 'ITS_complete':
                    report.write(f"  {sample}\n")
                else:
                    report.write(f"  {sample}: {count}\n")
            report.write("\n")

        # Categorize samples
        failed_samples = []
        complete_samples = []
        its1_only_samples = []
        its2_only_samples = []

        for sample in sorted(all_samples):
            complete_count = results[sample]['ITS_complete']
            its1_count = results[sample]['ITS1']
            its2_count = results[sample]['ITS2']

            if complete_count == 0 and its1_count == 0 and its2_count == 0:
                failed_samples.append(sample)
            elif complete_count > 0:
                complete_samples.append(sample)
            elif its1_count > 0 and its2_count == 0:
                its1_only_samples.append(sample)
            elif its2_count > 0 and its1_count == 0:
                its2_only_samples.append(sample)

        # Report failed samples
        report.write(f"Samples with no ITS regions extracted ({len(failed_samples)}):\n")
        if not failed_samples:
            report.write("  (none)\n")
        else:
            for sample in failed_samples:
                report.write(f"  {sample}\n")

        # Report samples without FASTA files
        if samples_without_fasta:
            report.write(f"\nSamples without FASTA files ({len(samples_without_fasta)}):\n")
            for sample in sorted(samples_without_fasta):
                report.write(f"  {sample}\n")

    # Print console summary
    logger.info("\n" + "=" * 70)
    logger.info("SAMPLE CATEGORISATION")
    logger.info("=" * 70)

    for region_name in REGIONS.keys():
        samples_with_seqs = [s for s in sorted(all_samples) if results[s][region_name] > 0]
        logger.info("\n{region_name} sequences %d:", len(samples_with_seqs))
        for sample in samples_with_seqs:
            count = results[sample][region_name]
            if region_name == 'ITS_complete':
                logger.info("  %s", sample)
            else:
                logger.info("  %s: %d", sample, count)

    logger.info("\nSamples with no ITS regions extracted (%d):", len(failed_samples))
    if not failed_samples:
        logger.info("  (none)")
    else:
        for sample in failed_samples:
            logger.info("  %s", sample)

    logger.info(f"\nJob completed at {datetime.now()}")
    logger.info(f"\nSummary report is available at: {report_file}")
    logger.info(f"CSV summary is available at: {csv_file}")
    logger.info("\nTo examine specific results:")
    logger.info(f"- ITS1 region: {output_dir / 'ITS1'}/")
    logger.info(f"- ITS2 region: {output_dir / 'ITS2'}/")
    logger.info(f"- Complete ITS region: {output_dir / 'ITS_complete'}/")

if __name__ == "__main__":
    main()
