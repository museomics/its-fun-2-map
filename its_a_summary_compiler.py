#!/usr/bin/env python3

import os
import glob
import argparse
from pathlib import Path
from datetime import datetime
import pandas as pd
from rpy2.robjects import r
from rpy2.robjects.conversion import localconverter
from rpy2.robjects import pandas2ri, default_converter
import logging
from Bio import SeqIO
from seqpy_tools import setup_logging

def find_csv_in_dir(project_dir, logger, custom):
    """Find CSV files in a directory (recursive)."""

    project_dir = str(project_dir)

    if not custom:
        patterns = [
            "**/UNITEd_*summary.csv",
            "**/*_json_*summary.csv",
            "**/mapping_*summary.csv",
            "**/assembly_*summary.csv",
            "**/blast_*summary.csv",
            "**/extraction_*summary.csv"
        ]

        csv_files = []
        for pattern in patterns:
            csv_files.extend(
                glob.glob(os.path.join(project_dir, pattern), recursive=True)
            )

        if not csv_files:
            logger.error("No default CSV files found in directory: %s", project_dir)

        return sorted(set(csv_files))

    if custom:
        csv_files = glob.glob(
            os.path.join(project_dir, "**", "*.csv"), recursive=True
        )
        logger.info(
            "Custom CSV search enabled. "
            "Consistent column renaming and filtering will not happen."
        )
        logger.info("Custom CSV files found: %s", csv_files)

        if not csv_files:
            logger.error("No CSV files found in directory: %s", project_dir)

        return csv_files


def find_id_column(df, csv_path):
    """Find the ID column in a dataframe, accepting various naming conventions."""
    valid_id_names = [
        "ID", "id", "Id",
        "sample_id", "sample_ID", "Sample_ID", "Sample_id", "SAMPLE_ID",
        "sampleid", "sampleID", "SampleID", "Sampleid", "SAMPLEID",
        "sample", "Sample", "SAMPLE",
    ]

    for col_name in valid_id_names:
        if col_name in df.columns:
            return col_name

    raise ValueError(
        f"No ID column found in {csv_path}. "
        f"Expected one of: {', '.join(valid_id_names)}"
    )


def load_and_prefix_df(csv_path, logger):
    """Load a CSV and normalise ID column name when using default file outputs.

    Special handling for fastp files: extracts sample IDs from filenames
    containing '_merged' and keeps only num_seqs and sum_len columns.
    """
    df = pd.read_csv(csv_path)
    prefix = os.path.basename(csv_path)
    prefix = prefix.removesuffix("_summary.csv")

    logger.info("Processing file: %s", csv_path)

    # Specific fastp handling (deprecated)
    if prefix.startswith("fastp"):
        if "file" not in df.columns:
            raise ValueError(f"'file' column not found in fastp CSV: {csv_path}")

        df = df[df["file"].str.contains("_merged", na=False)].copy()
        logger.info("  fastp: %d merged records retained", len(df))

        df["ID"] = df["file"].apply(
            lambda x: os.path.basename(x).split("_merged")[0]
        )

        required_cols = ["num_seqs", "sum_len"]
        missing = [c for c in required_cols if c not in df.columns]
        if missing:
            raise ValueError(
                f"Missing required columns in fastp CSV {csv_path}: {missing}"
            )

        df = df[["ID", "num_seqs", "sum_len"]].rename(
            columns={
                "num_seqs": f"{prefix}_n_sequences",
                "sum_len": f"{prefix}_sum_length",
            }
        )

        return df

    # Specific fastp JSON handling
    elif "json" in prefix:
        if "ID" not in df.columns:
            raise ValueError(f"'ID' column not found in fastp json CSV: {csv_path}")

        # Rename all other columns with prefix
        json_prefix = prefix.split('_')[0]
        new_columns = {col: f"{json_prefix}_{col}" if col != "ID" else col for col in df.columns}
        df = df.rename(columns=new_columns)
        return df


    # Specific blast handling
    elif prefix.startswith("blast"):
        if "round1" in prefix or "validation" in prefix:
            col_prefix = "blast_round1"
            rename_map = {}
            if "contig_path" in df.columns:
                rename_map["contig_path"] = f"{col_prefix}_contig_path"
                df = df.rename(columns=rename_map)
            return df

        else:
            if "its2" in prefix:
                prefix = "blast_its2"
            elif "its1" in prefix:
                prefix = "blast_its1"

            id_col = find_id_column(df, csv_path)
            if id_col != "ID":
                df = df.rename(columns={id_col: "ID"})
             # Rename all other columns with prefix
            new_columns = {col: f"{prefix}_{col}" if col != "ID" else col for col in df.columns}
            df = df.rename(columns=new_columns)
            return df

    else:
        # Rename ID column to "ID" to standardise across all input csvs
        id_col = find_id_column(df, csv_path)
        if id_col != "ID":
            df = df.rename(columns={id_col: "ID"})

        # Rename all other columns with prefix
        new_columns = {col: f"{prefix}_{col}" if col != "ID" else col for col in df.columns}
        df = df.rename(columns=new_columns)

        return df

def merge_outputs(project_dir, logger, custom):
    """Parse and combine csv outputs from multiple pipeline stages."""

    logger.info("Searching for CSV files...")
    files_list = find_csv_in_dir(project_dir, logger, custom)

    if not files_list:
        raise RuntimeError("No CSV files found - aborting merge.")

    logger.info("Found %d CSV files", len(files_list))

    dataframes = {}
    for csv_path in files_list:
        df = load_and_prefix_df(csv_path, logger)
        dataframes[csv_path] = df

    if not custom:
        logger.info("Validating ID columns...")
        if all("ID" in df.columns for df in dataframes.values()):
            dfs = list(dataframes.values())
            summary_df = dfs[0]
            for df in dfs[1:]:
                summary_df = summary_df.merge(df, on="ID", how="left")
        else:
            logger.warning(
                "ID column missing in at least one dataframe - column-binding instead"
            )
            summary_df = pd.concat(dataframes.values(), axis=1)
    else:
        summary_df = pd.concat(dataframes.values(), axis=1)

    return summary_df


def tidy_summary_df(summary_df):
    """
    Keep, order, and subset summary dataframe columns
    according to the defined pipeline order.
    Missing columns are ignored.
    """

    col_order = [
        "ID",
        "trimmed_before_total_reads",
        "trimmed_before_total_bases",
        "trimmed_before_q30_rate",
        "trimmed_after_total_reads",
        "trimmed_after_total_bases",
        "trimmed_after_q30_rate",
        "trimmed_after_read1_mean_length",
        "trimmed_after_read2_mean_length",
        "trimmed_after_gc_content",
        "trimmed_duplication_rate",
	"merged_after_total_reads",
        "merged_after_total_bases",
        "merged_after_q30_rate",
        "merged_after_read1_mean_length",
        "merged_after_gc_content",
        "merged_passed_filter_reads",
        "merged_duplication_rate",
        "UNITEd_searched_taxid",
        "UNITEd_searched_taxa",
        "UNITEd_number",
        "UNITEd_n_unique_taxa",
        "UNITEd_taxa",
        "UNITEd_min_length",
        "UNITEd_max_length",
        "UNITEd_avg_length",
        "mapping_n_reads",
        "mapping_perc",
        "assembly_method",
        "assembly_status",
        "assembly_n_scaffolds",
        "assembly_mean_scaffold_length",
        "assembly_n50",
        "blast_round1_found_taxon",
        "blast_round1_n_contigs_in",
        "blast_round1_n_contigs_hits",
        "blast_round1_hit_taxonomy",
        "blast_round1_correct_taxonomy",
        "blast_its2_n_contigs_in",
        "blast_its2_n_contigs_hits",
        "blast_its2_hit_taxonomy",
        "blast_its2_correct_taxonomy",
        "blast_its1_n_contigs_in",
        "blast_its1_n_contigs_hits",
        "blast_its1_hit_taxonomy",
        "blast_its1_correct_taxonomy",
        "extraction_ITS1",
        "extraction_ITS2",
        "extraction_ITS_complete",
        "overlapping_coords",
        "same_contig",
        "decision_description",
        "Final_outcome",
        "Final_contig_desc",
        "Final_contig"
    ]

    # Keep only columns that exist, in the desired order
    existing_cols = [c for c in col_order if c in summary_df.columns]

    return summary_df.loc[:, existing_cols]


def contig_analysis(summary_df):

    # Add contig analysis columns
    summary_df["Contig"] = None
    summary_df["Contig_desc"] = None
    summary_df["Contig_filepath"] = None


    summary_df["overlapping_coords"] = None
    summary_df["same_contig"] = None

    # ITS1/ITS2 coordinate analysis
    for i, row in summary_df.iterrows():
        its1_pass = row.get("blast_its1_correct_taxonomy") == "PASS"
        its2_pass = row.get("blast_its2_correct_taxonomy") == "PASS"

        # Both PASS
        if its1_pass and its2_pass:
            its1_contig = row.get("blast_its1_contigs")
            its2_contig = row.get("blast_its2_contigs")

            if its1_contig == its2_contig:
                # Compare coordinates
                try:
                    its1_coords = row.get("blast_its1_its_coords", "")
                    its2_coords = row.get("blast_its2_its_coords", "")

                    start1, end1 = map(int, its1_coords.split("-"))
                    start2, end2 = map(int, its2_coords.split("-"))

                    # Check if ranges overlap
                    if not (end1 < start2 or end2 < start1):
                        summary_df.at[i, "Contig_desc"] = "WARNING: Overlapping coordinates"
                        summary_df.at[i, "overlapping_coords"] = "YES"
                        summary_df.at[i, "same_contig"] = "YES"
                    else:
                        summary_df.at[i, "Contig_desc"] = "Contiguous"
                        summary_df.at[i, "overlapping_coords"] = "NO"
                        summary_df.at[i, "Contig"] = its1_contig
                        summary_df.at[i, "same_contig"] = "YES"
                        summary_df.at[i, "Contig_filepath"] = row.get("blast_its1_contig_path")
                except (ValueError, AttributeError) as e:
                    summary_df.at[i, "Contig_desc"] = "WARNING: Could not parse coordinates"
            else:
                summary_df.at[i, "Contig"] = its2_contig
                summary_df.at[i, "Contig_desc"] = "Not on same contig - defaulting to ITS2"
                summary_df.at[i, "same_contig"] = "NO"
                summary_df.at[i, "Contig_filepath"] = row.get("blast_its2_contig_path")

        # ITS1 PASS only
        elif its1_pass and not its2_pass:
            summary_df.at[i, "Contig"] = row.get("blast_its1_contigs")
            summary_df.at[i, "Contig_desc"] = "ITS1 only"
            summary_df.at[i, "Contig_filepath"] = row.get("blast_its1_contig_path")

        # ITS2 PASS only
        elif its2_pass and not its1_pass:
            summary_df.at[i, "Contig"] = row.get("blast_its2_contigs")
            summary_df.at[i, "Contig_desc"] = "ITS2 only"
            summary_df.at[i, "Contig_filepath"] = row.get("blast_its2_contig_path")

        # Both FAIL
        else:
            summary_df.at[i, "Contig"] = None
            summary_df.at[i, "Contig_desc"] = "FAILED CONTIG"
            summary_df.at[i, "Contig_filepath"] = None

    return summary_df

def process_renaming_df(original_df, new_df):
    # Check number of columns
    if new_df.shape[1] == 2:
        # Check required column exists in summary_df
        if "ID" in original_df.columns:
            updated_df = original_df.merge(new_df, on="ID", how="left")
        else:
            logger.warning(
                "ID column missing in summary dataframe; cannot merge naming table."
            )
            updated_df = original_df
    else:
        logger.warning(
            "Naming TSV must have exactly 2 columns. "
            f"Found {new_df.shape[1]} columns skipping renaming."
        )
        updated_df = original_df

    return updated_df

def rename_fastas(df, name_column, out_dir, logger, delim="_"):
    """Rename FASTA headers based on dataframe metadata."""

    pass_dir = Path(out_dir) / "pass_fastas"
    pass_dir.mkdir(parents=True, exist_ok=True)

    man_verification_dir = Path(out_dir) / "manual_verification"
    man_verification_dir.mkdir(parents=True, exist_ok=True)

    # Normalise NAs
    NA_STRINGS = {
        "NA",
        "NA_character_",
        "NA_real_",
        "NA_integer_",
        "NA_logical_"
    }

    string_cols = df.select_dtypes(include=["object", "string"]).columns

    df[string_cols] = (
        df[string_cols]
        .apply(lambda col: col.astype("string"))
        .replace(NA_STRINGS, pd.NA)
   )

    # Filter rows: non-NA Final_contig and PASS outcome only
    filtered_df = df[
        df["Final_contig"].notna() &
        (df["Final_outcome"] == "PASS")
    ]

    for _, row in filtered_df.iterrows():

        fasta_path = Path(str(row["Final_contig"]).strip())

        if not fasta_path.exists():
            logger.warning(f"FASTA not found: {fasta_path}")
            continue

        if fasta_path in (Path("."), Path("")):
            logger.error(f"Final_contig is invalid placeholder: {fasta_path}")
            continue

        if not fasta_path.is_file():
            logger.error(f"Final_contig is not a FASTA file: {fasta_path}")
            continue

        # Read FASTA
        records = list(SeqIO.parse(fasta_path, "fasta"))

        # Enforce single sequence per FASTA
        if len(records) != 1:
            logger.warning(
                f"Multiple sequences found in {fasta_path}. Skipping..."
            )
            continue

        record = records[0]

        # Construct new header
        new_header = f"{row[name_column]}{delim}{row['Final_contig_desc']}"

        # Update FASTA header
        record.id = new_header
        record.name = new_header
        record.description = ""

        # Output FASTA path
        new_fasta_name = (
            pass_dir /
            f"{row[name_column]}_{row['Final_contig_desc']}.fasta"
        )

        # Write FASTA
        SeqIO.write(record, new_fasta_name, "fasta")



    filtered_df = df[
        df["Final_contig"].notna() &
        (df["Final_outcome"] == "FAIL")
    ]
    filtered_df = df[
    df["Final_contig"].notna() &
    (df["Final_outcome"] == "FAIL")
    ]

    logger.info("FAIL rows to process: %d", len(filtered_df))
    for _, row in filtered_df.iterrows():
        logger.info("Processing ID=%s, Final_contig='%s'", row.get("ID"), row.get("Final_contig"))

        fasta_path = Path(row["Final_contig"])

        if not fasta_path.exists():
            logger.warning(f"FASTA not found: {fasta_path}")
            continue

        # Read FASTA
        records = SeqIO.parse(fasta_path, "fasta")

        # Output FASTA path
        manual_verification_name = (
            man_verification_dir /
            f"{row[name_column]}_manual_verification_needed.fasta"
        )

        # Write FASTA
        SeqIO.write(records, manual_verification_name, "fasta")

def n_taxa(y):
    if pd.isna(y) or y == "":
        return 0
    v = y.split(";")
    return len(set(v))


def main():
    parser = argparse.ArgumentParser(
        description="Merge pipeline CSV outputs into a single summary table"
    )
    parser.add_argument(
        "project_dir",
        help="Project directory containing pipeline output CSV files",
    )
    parser.add_argument(
        "-o",
        "--output",
        default=None,
        help="Output CSV file (default: merged_summary_<date>.csv)",
    )
    parser.add_argument(
        "--custom",
        action="store_true",
        help="Use custom list of CSV files found (disabled for using all files)",
    )
    parser.add_argument(
        "--log",
        default="merge_outputs.log",
        help="Log file name",
    )

    parser.add_argument(
        "--naming_tsv",
        help="(Optional) Tab-delimited text file with ID and new_name (no headers) for renaming FASTA outputs for submission.",
    )

    parser.add_argument(
        "--delim",
        default = "|",
        help = "Optional delimiter naming output FASTA files depending on submission location."
        "See BOlD for futher information: https://v3.boldsystems.org/index.php/resources/handbook?chapter=3_submissions.html")


    args = parser.parse_args()

    out_dir = Path(args.project_dir) / "final_results_dir"
    os.makedirs(out_dir, exist_ok=True)
    delim = args.delim

    # Set up log file with timestamp if not provided
    if args.log is None:
        timestamp = datetime.now().strftime('%Y%m%d_%H%M%S')
        args.log = f'merged_outputs_{timestamp}.log'
    
    logger = setup_logging(log_file=args.log)

    logger.info("Starting merge at %s", datetime.now().isoformat())

    # Merge all summary files from pipeline into single df
    summary_df = merge_outputs(
        args.project_dir,
        logger,
        custom=args.custom
    )

    # Check for contiguity and coordinate overlap from blastn results
    summary_df = contig_analysis(summary_df)

    # Calculate the n. unique taxa calculation for UNITEd
    summary_df["UNITEd_n_unique_taxa"] = summary_df["UNITEd_taxa"].apply(n_taxa)

    paths_cols = [c for c in summary_df.columns if 'path' in c.lower()]
    logger.info("Path columns available: %s", paths_cols)
 
    # Assess for which sequences/contigs can be taken forward for submission
    r_script = Path(__file__).parent / 'its_decision_making.R'
    r['source'](str(r_script))
    decision_making = r['its_outcome']

    # normalise column types (expected mixture of float and str for str columns with NaN)
    summary_df = summary_df.fillna("").astype("string")

    with localconverter(default_converter + pandas2ri.converter):
        result = decision_making(summary_df)

    summary_df = pd.DataFrame(result)

    # Write new FASTAs out with edited headers and file naming
    if args.naming_tsv:
        logger.info("Found naming file %s.", args.naming_tsv)
        naming_tsv = Path(args.naming_tsv)
        new_names_index = pd.read_csv(naming_tsv, header = None, sep = "\t", names = ["ID", "New_Name"])

        # Clean naming and ID column in summary df of whitespaces
        new_names_index["ID"] = new_names_index["ID"].str.strip()
        new_names_index["New_Name"] = new_names_index["New_Name"].str.strip()
        summary_df["ID"] = summary_df["ID"].astype(str).str.strip()

        naming_df = process_renaming_df(summary_df, new_names_index)
        n_matched = naming_df["New_Name"].notna().sum()
        logger.info("Renaming table matched %d / %d samples", n_matched, len(naming_df))

        rename_fastas(naming_df, 'New_Name', out_dir, logger, delim=delim)
    else:
        naming_df = summary_df
        rename_fastas(naming_df, 'ID', out_dir, logger, delim=delim)


    # Clean up summary dataframe
    summary_df = tidy_summary_df(naming_df)

    if args.output is None:
        date = datetime.now().strftime("%Y%m%d")
        args.output = f"merged_summary_{date}.csv"

    outfile = out_dir / args.output
    summary_df.to_csv(outfile, index=False)
    logger.info("Merged summary written to %s", outfile)


if __name__ == "__main__":
    main()
