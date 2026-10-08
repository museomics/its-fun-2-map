"""Tests for primer and region loading in its_primer_binding."""

import pytest

from its_fun_2_map.its_primer_binding import (
    DEFAULT_PRIMERS,
    DEFAULT_REGIONS,
    extract_sample_name,
    load_primers_and_regions,
    load_primers_tsv,
    load_regions_tsv,
)

PRIMERS_TSV = "name\tsequence\nfwd1\tacgtacgt\nrev1\t TTGGCCAA \n"
REGIONS_TSV = "region\tforward\treverse\nmy_region\tfwd1\trev1\n"


def test_defaults_are_white_et_al_its_primers():
    primers, regions = load_primers_and_regions()

    assert primers == DEFAULT_PRIMERS
    assert regions == DEFAULT_REGIONS
    assert set(regions) == {"ITS1", "ITS2", "ITS_complete"}
    # Every default region must reference defined primers
    for fwd, rev in regions.values():
        assert fwd in primers and rev in primers


def test_load_primers_strips_and_uppercases(write_file):
    primers = load_primers_tsv(write_file("primers.tsv", PRIMERS_TSV))

    assert primers == {"fwd1": "ACGTACGT", "rev1": "TTGGCCAA"}


def test_load_primers_requires_columns(write_file):
    with pytest.raises(ValueError, match="Primer TSV must contain columns"):
        load_primers_tsv(write_file("primers.tsv", "id\tseq\nfwd1\tACGT\n"))


def test_load_regions(write_file):
    primers = {"fwd1": "ACGT", "rev1": "TTGG"}

    assert load_regions_tsv(write_file("regions.tsv", REGIONS_TSV), primers) == {"my_region": ("fwd1", "rev1")}


def test_region_with_undefined_primer_is_rejected(write_file):
    with pytest.raises(ValueError, match="undefined primer"):
        load_regions_tsv(write_file("regions.tsv", REGIONS_TSV), {"fwd1": "ACGT"})


def test_custom_primers_and_regions_replace_defaults(write_file):
    primers, regions = load_primers_and_regions(
        write_file("primers.tsv", PRIMERS_TSV), write_file("regions.tsv", REGIONS_TSV)
    )

    assert set(primers) == {"fwd1", "rev1"}
    assert regions == {"my_region": ("fwd1", "rev1")}


@pytest.mark.parametrize("which", ["primers", "regions"])
def test_primers_and_regions_must_be_given_together(which, write_file):
    path = write_file("only.tsv", PRIMERS_TSV)
    kwargs = {"primers_tsv": path} if which == "primers" else {"regions_tsv": path}

    with pytest.raises(ValueError, match="must be provided together"):
        load_primers_and_regions(**kwargs)


@pytest.mark.parametrize("filename, expected", [
    ("NHM001_scaffolds.fasta", "NHM001"),
    ("dir/NHM001_ITS2_contig.fa", "NHM001"),
    ("NHM001.fasta", "NHM001"),
])
def test_extract_sample_name(filename, expected):
    assert extract_sample_name(filename) == expected
