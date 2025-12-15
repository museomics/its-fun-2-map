#!/usr/bin/env python3

"""
ITS primer binding script, employing seqkit amplicon and 4 primers from White et al (1990),
able to extract the full ITS region (ITS1->ITS4), the ITS1 region (ITS1->ITS2), and the
ITS2 region (ITS3->ITS4). ITS3 primer is a reverse compliment of ITS2 (see White et al., 1990)

Author: D. Parsons & M. Kamouyiaros (NHMUK)
Version: 2.2.0
"""

import os
import sys
import subprocess
import argparse
import csv
from pathlib import Path
from datetime import datetime
from collections import defaultdict

# ITS primer sequences from White et al. 1990
PRIMERS = {
    'ITS1': 'TCCGTAGGTGAACCTGCGG',
    'ITS2': 'GCTGCGTTCTTCATCGATGC',
    'ITS3': 'GCATCGATGAAGAACGCAGC',
    'ITS4': 'TCCTCCGCTTATTGATATGC'
}

# Primer combinations for each region
REGIONS = {
    'ITS1': ('ITS1', 'ITS2'),
    'ITS2': ('ITS3', 'ITS4'),
    'ITS_complete': ('ITS1', 'ITS4')
}


def read_tracking_sheet(tracking_sheet_path, column_name='ID'):
    """Read tracking sheet and return set of sample IDs"""
    sample_ids = set()
    
    try:
        with open(tracking_sheet_path, 'r') as csvfile:
            reader = csv.DictReader(csvfile)
            
            if column_name not in reader.fieldnames:
                print(f"ERROR: Tracking sheet must contain {column_name} column", file=sys.stderr)
                sys.exit(1)
            
            for row in reader:
                sample_id = row[column_name].strip()
                if sample_id:
                    sample_ids.add(sample_id)
        
        return sample_ids
    
    except FileNotFoundError:
        print(f"ERROR: Tracking sheet not found: {tracking_sheet_path}", file=sys.stderr)
        sys.exit(1)
    except Exception as e:
        print(f"ERROR: Failed to read tracking sheet: {e}", file=sys.stderr)
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


def run_seqkit_amplicon(input_file, output_file, forward_primer, reverse_primer):
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
        print(f"Error running seqkit amplicon: {e}", file=sys.stderr)
        return False


def parse_args():
    """Parse command line arguments"""
    parser = argparse.ArgumentParser(
        description='ITS primer binding script using seqkit amplicon and primers from White et al (1990)',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  %(prog)s --input /path/to/fastas --output /path/to/results --tracking_sheet samples.csv
  %(prog)s -i input_dir/ -o output_dir/ -t tracking.csv

Regions extracted:
  - ITS1 region (ITS1-F + ITS2-R)
  - ITS2 region (ITS3-F + ITS4-R)
  - Complete ITS region (ITS1-F + ITS4-R)
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
    
    args = parser.parse_args()
    
    # Validate input directory exists
    if not args.input.exists():
        parser.error(f"Input directory does not exist: {args.input}")
    
    if not args.input.is_dir():
        parser.error(f"Input path is not a directory: {args.input}")
    
    # Validate tracking sheet exists
    if not args.tracking_sheet.exists():
        parser.error(f"Tracking sheet does not exist: {args.tracking_sheet}")
    
    return args


def main():
    args = parse_args()
    
    input_dir = args.input
    output_dir = args.output
    tracking_sheet = args.tracking_sheet
    column_name = args.column_name if hasattr(args, 'column_name') else 'ID'
    
    # Read sample IDs from tracking sheet
    print(f"Reading sample IDs from tracking sheet: {tracking_sheet}")
    all_samples = read_tracking_sheet(tracking_sheet, column_name=column_name)
    
    if not all_samples:
        print("ERROR: No sample IDs found in tracking sheet")
        sys.exit(1)
    
    print(f"Found {len(all_samples)} samples in tracking sheet")
    print(f"Samples: {sorted(all_samples)}")
    
    # Create output directories
    for region in REGIONS.keys():
        (output_dir / region).mkdir(parents=True, exist_ok=True)
    
    print(f"\nStarting ITS primer alignment at {datetime.now()}")
    print(f"Processing files from directory: {input_dir}")
    
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
        print(f"\nWARNING: {len(samples_without_fasta)} samples from tracking sheet have no FASTA files:")
        for sample in sorted(samples_without_fasta):
            print(f"  {sample}")
    
    # Generate summary statistics
    print("\nGenerating summary report...")
    
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
            'extraction_ITS1',
            'ITS1_path',
            'extraction_ITS2',
            'ITS2_path',
            'extraction_ITS_complete',
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
                row[f'extraction_{region_name}'] = extracted  # Instead of f'{region_name}_extracted'
                row[f'{region_name}_path'] = file_path
            
            writer.writerow(row)
    
    print(f"CSV summary written to: {csv_file}")
    
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
    print("\n" + "=" * 70)
    print("SAMPLE CATEGORISATION")
    print("=" * 70)
    
    for region_name in REGIONS.keys():
        samples_with_seqs = [s for s in sorted(all_samples) if results[s][region_name] > 0]
        print(f"\n{region_name} sequences ({len(samples_with_seqs)}):")
        for sample in samples_with_seqs:
            count = results[sample][region_name]
            if region_name == 'ITS_complete':
                print(f"  {sample}")
            else:
                print(f"  {sample}: {count}")
    
    print(f"\nSamples with no ITS regions extracted ({len(failed_samples)}):")
    if not failed_samples:
        print("  (none)")
    else:
        for sample in failed_samples:
            print(f"  {sample}")
    
    print(f"\nJob completed at {datetime.now()}")
    print(f"\nSummary report is available at: {report_file}")
    print(f"CSV summary is available at: {csv_file}")
    print("\nTo examine specific results:")
    print(f"- ITS1 region: {output_dir / 'ITS1'}/")
    print(f"- ITS2 region: {output_dir / 'ITS2'}/")
    print(f"- Complete ITS region: {output_dir / 'ITS_complete'}/")

if __name__ == "__main__":
    main()
