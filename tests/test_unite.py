"""Tests for UNITE header parsing and taxon lookup helpers in UNITEd."""

import pytest

from its_fun_2_map.UNITEd import (
    get_child_rank,
    get_most_specific_taxon,
    parse_fasta_header,
    search_rank_for_taxon,
)

UNITE_HEADER = ">Fusarium_oxysporum|KX123456|SH1234567.10FU|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;o__Hypocreales;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum"


def test_parse_unite_header():
    parsed = parse_fasta_header(UNITE_HEADER)

    assert parsed["id"] == "Fusarium_oxysporum"
    assert parsed["header"] == UNITE_HEADER.lstrip(">")
    assert parsed["taxonomy"] == {
        "k": "Fungi",
        "p": "Ascomycota",
        "c": "Sordariomycetes",
        "o": "Hypocreales",
        "f": "Nectriaceae",
        "g": "Fusarium",
        "s": "Fusarium_oxysporum",
    }


def test_header_without_taxonomy_returns_none():
    assert parse_fasta_header(">KX123456 some description") is None


@pytest.mark.parametrize("lineage, expected", [
    ({"genus": "Fusarium", "species": "Fusarium oxysporum"}, "Fusarium oxysporum"),
    ({"kingdom": "Fungi", "family": "Nectriaceae"}, "Nectriaceae"),
    ({}, "Unknown"),
])
def test_most_specific_taxon(lineage, expected):
    assert get_most_specific_taxon(lineage) == expected


@pytest.mark.parametrize("parent, child", [
    ("genus", "species"),
    ("kingdom", "phylum"),
    ("species", None),
    ("not_a_rank", None),
])
def test_child_rank(parent, child):
    assert get_child_rank(parent) == child


@pytest.fixture
def unite_index():
    return {
        "genus": {"Fusarium": ["seq_g1"]},
        "species": {"Fusarium_oxysporum": ["seq_s1"]},
    }


def test_search_exact_and_case_insensitive(unite_index):
    assert search_rank_for_taxon(unite_index, "genus", "Fusarium") == ["seq_g1"]
    assert search_rank_for_taxon(unite_index, "genus", "fusarium") == ["seq_g1"]
    assert search_rank_for_taxon(unite_index, "genus", "Aspergillus") is None


def test_search_species_ignores_space_underscore_differences(unite_index):
    assert search_rank_for_taxon(unite_index, "species", "Fusarium oxysporum") == ["seq_s1"]
