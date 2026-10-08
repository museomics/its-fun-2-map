"""Tests for the BLAST round 1 and round 2 parser helpers."""

import logging

import pandas as pd
import pytest

from its_fun_2_map import blast_output_parser, blast_round1_parser

STANDARD_COLUMNS = [
    "qseqid", "sseqid", "pident", "length", "mismatch", "gapopen",
    "qstart", "qend", "sstart", "send", "evalue", "bitscore",
]
# UNITE-style subject ID: the taxonomy string is the last |-separated field
SSEQID = "Fusarium_oxysporum|KX123|SH1|reps|k__Fungi;p__Ascomycota;c__Sordariomycetes;f__Nectriaceae;g__Fusarium;s__Fusarium_oxysporum"
HIT_ROW = "\t".join(["contig1", SSEQID, "99.5", "500", "2", "0", "1", "500", "1", "500", "1e-50", "900"])


@pytest.fixture(autouse=True)
def round1_logger(monkeypatch):
    # blast_round1_parser's helpers log through a module-level `logger` that
    # only run() creates (via `global logger`), so set one for direct calls.
    monkeypatch.setattr(
        blast_round1_parser, "logger", logging.getLogger("its_fun_2_map.tests"), raising=False
    )


@pytest.mark.parametrize("parser", [blast_round1_parser, blast_output_parser])
def test_detect_header(parser, write_file):
    tsv = write_file("with_header.tsv", "\t".join(STANDARD_COLUMNS) + "\n" + HIT_ROW + "\n")

    has_header, columns = parser.detect_blast_format(str(tsv))

    assert has_header is True
    assert columns == STANDARD_COLUMNS


@pytest.mark.parametrize("parser", [blast_round1_parser, blast_output_parser])
def test_detect_headerless(parser, write_file):
    tsv = write_file("no_header.tsv", HIT_ROW + "\n")

    has_header, columns = parser.detect_blast_format(str(tsv))

    assert has_header is False
    assert columns == STANDARD_COLUMNS


def test_round1_detect_skips_comment_lines(write_file):
    tsv = write_file("comments.tsv", "# BLASTN 2.17.0\n# Query: x\n" + "\t".join(STANDARD_COLUMNS) + "\n")

    assert blast_round1_parser.detect_blast_format(str(tsv))[0] is True


def test_taxonomy_extractors():
    assert blast_output_parser.extract_species_from_taxonomy(SSEQID) == "fusarium oxysporum"
    assert blast_output_parser.extract_genus_from_taxonomy(SSEQID) == "fusarium"
    assert blast_output_parser.extract_family_from_taxonomy(SSEQID) == "nectriaceae"
    assert blast_output_parser.extract_full_taxonomy_from_sseqid(SSEQID).startswith("k__Fungi;")

    assert blast_output_parser.extract_species_from_taxonomy("no_taxonomy_here") is None
    assert blast_output_parser.extract_full_taxonomy_from_sseqid("no_taxonomy_here") is None


def test_grep_matches_whole_words_and_ignores_underscores():
    assert blast_round1_parser._grep("s__Fusarium_oxysporum", ["fusarium oxysporum"])
    assert blast_round1_parser._grep("g__Fusarium", ["Fusarium"])
    # "Fusar" is only part of a word, so it must not match
    assert not blast_round1_parser._grep("g__Fusarium", ["Fusar"])


def test_filter_df_keeps_matching_rows():
    df = pd.DataFrame({"sseqid": ["g__Fusarium;s__x", "g__Aspergillus;s__y"]})

    assert list(blast_round1_parser.filter_df(df, ["Fusarium"])["sseqid"]) == ["g__Fusarium;s__x"]
    assert blast_round1_parser.filter_df(df, ["Penicillium"]).empty


def blast_df():
    return pd.DataFrame({
        "qseqid": ["q1", "q1", "q1", "q2"],
        "sseqid": ["a", "b", "c", "d"],
        "pident": [99.0, 95.0, 80.0, 99.0],
        "length": [300, 500, 600, 100],
        "evalue": [1e-50, 1e-60, 1e-10, 1e-3],
    })


def test_filter_blast_applies_thresholds():
    result = blast_round1_parser.filter_blast(blast_df(), min_len=200, min_pident=90.0)

    # c fails pident, d fails length (and the default 1e-5 evalue cutoff)
    assert set(result["sseqid"]) == {"a", "b"}


def test_filter_blast_top_n_prefers_longest_alignment():
    result = blast_round1_parser.filter_blast(blast_df(), top_n=1, evalue_cutoff=None)

    assert dict(zip(result["qseqid"], result["sseqid"])) == {"q1": "c", "q2": "d"}


def test_round1_expected_taxonomy_from_filename():
    mapping = {
        "NHM001": {"family": "Nectriaceae", "genus": "Fusarium", "full_taxonomy": "f__Nectriaceae;g__Fusarium"},
    }

    assert blast_round1_parser.get_expected_taxonomy_from_filename(
        "results/NHM001_blast_results.tsv", mapping, "Genus"
    ) == ("Fusarium", "f__Nectriaceae;g__Fusarium", "NHM001")
    # A sample ID that is only a prefix of another ID must not match
    assert blast_round1_parser.get_expected_taxonomy_from_filename(
        "NHM0012_blast_results.tsv", mapping, "genus"
    ) == (None, None, None)


def test_round2_expected_taxonomy_from_filename():
    mapping = {"NHM001": {"family": "Nectriaceae", "full_taxonomy": "f__Nectriaceae"}}

    assert blast_output_parser.get_expected_taxonomy_from_filename(
        "ITS2_NHM001_blast_results.tsv", mapping
    ) == ("Nectriaceae", "f__Nectriaceae", "NHM001")
    assert blast_output_parser.get_expected_taxonomy_from_filename("other.tsv", mapping) == (None, None, None)
