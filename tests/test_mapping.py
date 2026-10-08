"""Tests for samtools flagstat parsing and read-file classification in mapping_module."""

import pytest

from its_fun_2_map.mapping_module import is_merged, is_unmerged_1, is_unmerged_2, parse_flagstats

# Excerpt of `samtools flagstat -O tsv` output (passed QC, failed QC, category)
FLAGSTAT_TSV = (
    "1000\t0\tprimary\n"
    "850\t0\tmapped\n"
    "85.00%\tN/A\tmapped %\n"
    "840\t0\tprimary mapped\n"
    "84.00%\tN/A\tprimary mapped %\n"
)


def test_parse_flagstats(write_file):
    flagstats = write_file("sample.flagstat", FLAGSTAT_TSV)

    # "primary mapped" lines must not overwrite the overall mapped values
    assert parse_flagstats(str(flagstats)) == (850, 85.0)


def test_parse_flagstats_na_percentage(write_file):
    flagstats = write_file("sample.flagstat", "0\t0\tmapped\nN/A\tN/A\tmapped %\n")

    assert parse_flagstats(str(flagstats)) == (0, None)


def test_parse_flagstats_missing_file(tmp_path):
    assert parse_flagstats(str(tmp_path / "missing.flagstat")) == (None, None)


@pytest.mark.parametrize("name, merged, r1, r2", [
    ("S1_merged.fq", True, False, False),
    ("S1_merged.fastq.gz", True, False, False),
    ("S1_unmerged_1.fq", False, True, False),
    ("S1_unmerged_2.fastq.gz", False, False, True),
    ("S1_trimmed_1.fq", False, False, False),
    ("S1_merged.txt", False, False, False),
])
def test_read_file_classification(name, merged, r1, r2):
    assert is_merged(name) is merged
    assert is_unmerged_1(name) is r1
    assert is_unmerged_2(name) is r2
