"""Tests for the bundled R scripts, run through rpy2 exactly as the pipeline does.

Skipped when rpy2 cannot start an embedded R (for example, no R installed).
CI installs R, so these always run there.
"""

import json
from importlib.resources import files

import pandas as pd
import pytest

try:
    from rpy2.robjects import default_converter, pandas2ri, r
    from rpy2.robjects.conversion import localconverter
    r("1")  # fail here, not mid-test, if the embedded R cannot start
except Exception as exc:  # ImportError, or rpy2 failing to load libR
    pytest.skip(f"rpy2/R not available: {exc}", allow_module_level=True)


def fastp_report(total_reads, read2_mean_length=None):
    """Minimal fastp JSON report with the fields parse_fastp_json.R reads."""
    after = {
        "total_reads": total_reads // 2,
        "total_bases": total_reads * 75,
        "q30_bases": total_reads * 70,
        "q30_rate": 0.93,
        "read1_mean_length": 150,
        "gc_content": 0.48,
    }
    if read2_mean_length is not None:
        after["read2_mean_length"] = read2_mean_length
    return {
        "summary": {
            "before_filtering": {
                "total_reads": total_reads,
                "total_bases": total_reads * 150,
                "q30_bases": total_reads * 140,
                "q30_rate": 0.91,
                "gc_content": 0.47,
            },
            "after_filtering": after,
        },
        "filtering_result": {
            "passed_filter_reads": total_reads // 2,
            "low_quality_reads": 10,
            "too_short_reads": 5,
        },
        "duplication": {"rate": 0.12},
    }


def test_fastp_json_summary(tmp_path, logger):
    from its_fun_2_map.fastp_module import run_fastp_json_summary

    reports = {
        "NHM001_trim.json": fastp_report(1000, read2_mean_length=148),
        "NHM002_trim.json": fastp_report(2000, read2_mean_length=147),
        # Merged reads are single-end, so fastp reports no read2 length
        "NHM001_merge.json": fastp_report(400),
    }
    for name, report in reports.items():
        (tmp_path / name).write_text(json.dumps(report))

    run_fastp_json_summary(json_dir=tmp_path, outdir=tmp_path, logger=logger)

    trimmed = pd.read_csv(tmp_path / "trimmed_json_summary.csv")
    merged = pd.read_csv(tmp_path / "merged_json_summary.csv")

    assert sorted(trimmed["ID"]) == ["NHM001", "NHM002"]
    assert list(merged["ID"]) == ["NHM001"]

    row = trimmed.set_index("ID").loc["NHM002"]
    assert row["before_total_reads"] == 2000
    assert row["after_read2_mean_length"] == 147
    assert row["duplication_rate"] == pytest.approx(0.12)
    assert pd.isna(merged["after_read2_mean_length"].iloc[0])


# Columns its_outcome reads. Values not set per scenario default to "".
DECISION_COLUMNS = [
    "ID", "Contig_desc",
    "blast_round1_correct_taxonomy", "blast_round1_contig_path",
    "blast_its2_correct_taxonomy", "blast_its2_contig_path",
    "blast_its1_correct_taxonomy", "blast_its1_contig_path",
    "overlapping_coords", "same_contig",
    "extraction_ITS_complete", "extraction_ITS_complete_path",
    "extraction_ITS2", "extraction_ITS2_path",
]

BOTH_PASS = {
    "blast_round1_correct_taxonomy": "PASS",
    "blast_its2_correct_taxonomy": "PASS",
    "blast_its1_correct_taxonomy": "PASS",
}

SCENARIOS = {
    # id: (row values, expected decision_description, Final_outcome, Final_contig)
    "s1_putative": (
        {**BOTH_PASS, "overlapping_coords": "NO", "same_contig": "YES", "extraction_ITS_complete": "NO"},
        "Putative ITS region", "PASS", "its2.fa",
    ),
    "s2_complete": (
        {**BOTH_PASS, "overlapping_coords": "NO", "same_contig": "YES", "extraction_ITS_complete": "YES"},
        "Complete ITS region", "PASS", "complete.fa",
    ),
    "s10_its2_only": (
        {**BOTH_PASS, "blast_its1_correct_taxonomy": "FAIL",
         "extraction_ITS_complete": "NO", "extraction_ITS2": "YES"},
        "ITS2 found only", "PASS", "its2_extract.fa",
    ),
    "s13_not_contiguous": (
        {**BOTH_PASS, "same_contig": "NO", "extraction_ITS_complete": "NO"},
        "WARNING: Blast results not contiguous - defaulting to ITS2 contig", "PASS", "its2.fa",
    ),
    "s18_manual_curation": (
        {"blast_round1_correct_taxonomy": "PASS", "blast_its2_correct_taxonomy": "FAIL",
         "blast_its1_correct_taxonomy": "FAIL"},
        "Not enough support - manual curation required", "FAIL", "round1.fa",
    ),
    "s19_round1_fail": (
        {"blast_round1_correct_taxonomy": "FAIL"},
        "Failed all checks", "FAIL", None,
    ),
    "s20_failed_contig": (
        {"Contig_desc": "FAILED CONTIG"},
        "Failed all checks", "FAIL", None,
    ),
}


def run_its_outcome(rows):
    """Call its_decision_making.R the way its_a_summary_compiler does."""
    df = pd.DataFrame(rows, columns=DECISION_COLUMNS).fillna("").astype("string")
    r["source"](str(files("its_fun_2_map") / "its_decision_making.R"))
    with localconverter(default_converter + pandas2ri.converter):
        result = r["its_outcome"](df)
    return pd.DataFrame(result)


def test_its_outcome_decision_table():
    paths = {
        "blast_round1_contig_path": "round1.fa",
        "blast_its2_contig_path": "its2.fa",
        "blast_its1_contig_path": "its1.fa",
        "extraction_ITS_complete_path": "complete.fa",
        "extraction_ITS2_path": "its2_extract.fa",
    }
    rows = [{"ID": sid, **paths, **values} for sid, (values, *_) in SCENARIOS.items()]

    result = run_its_outcome(rows).set_index("ID")

    for sid, (_, description, outcome, contig) in SCENARIOS.items():
        row = result.loc[sid]
        assert row["decision_description"] == description, sid
        assert row["Final_outcome"] == outcome, sid
        if contig is None:
            # R's NA comes back as an rpy2 NA object, not a pandas NA (see below)
            assert row["Final_contig"] not in set(paths.values()), sid
        else:
            assert row["Final_contig"] == contig, sid


def test_rename_fastas_with_passing_and_failed_samples(tmp_path, logger):
    from its_fun_2_map.its_a_summary_compiler import rename_fastas

    contig = tmp_path / "NHM001_its2.fa"
    contig.write_text(">contig1\nACGTACGTAC\n")
    passing, _, _, _ = SCENARIOS["s1_putative"]
    failing, _, _, _ = SCENARIOS["s19_round1_fail"]
    rows = [
        {"ID": "NHM001", "blast_its2_contig_path": str(contig), **passing},
        {"ID": "NHM002", **failing},
    ]
    result = run_its_outcome(rows)

    # rename_fastas normalises rpy2's NA_character_ strings, so the failed
    # sample is skipped and only the passing contig is written
    rename_fastas(result, "ID", tmp_path / "out", logger)

    assert [p.name for p in (tmp_path / "out" / "pass_fastas").iterdir()] == [
        "NHM001_putative_ITS_region.fasta"
    ]
    assert not any((tmp_path / "out" / "manual_verification").iterdir())


@pytest.mark.xfail(
    raises=TypeError,
    strict=False,
    reason="Known issue: when every sample fails, R returns Final_contig as an "
           "all-NA column, which rpy2 converts to integers (-2147483648) rather than "
           "strings. rename_fastas only normalises NA in string columns, so it "
           "passes the integer to Path() and raises TypeError.",
)
def test_rename_fastas_when_every_sample_fails(tmp_path, logger):
    from its_fun_2_map.its_a_summary_compiler import rename_fastas

    rows = [
        {"ID": "NHM001", **SCENARIOS["s19_round1_fail"][0]},
        {"ID": "NHM002", **SCENARIOS["s20_failed_contig"][0]},
    ]
    result = run_its_outcome(rows)

    rename_fastas(result, "ID", tmp_path / "out", logger)

    assert not any((tmp_path / "out" / "pass_fastas").iterdir())
