"""Tests for the shared helpers in its_fun_tools.

get_ncbi_lineage is tested with Bio.Entrez faked out, so no request ever
reaches NCBI and the retry back-off does not actually sleep.
"""

import urllib.error

import pytest

from its_fun_2_map import its_fun_tools
from its_fun_2_map.its_fun_tools import cleanup_temp_dir, get_ncbi_lineage, load_name_ids

TAXONOMY_RECORD = [{
    "TaxId": "5507",
    "ScientificName": "Fusarium oxysporum",
    "Rank": "species",
    "LineageEx": [
        {"TaxId": "2759", "ScientificName": "Eukaryota", "Rank": "superkingdom"},
        {"TaxId": "4751", "ScientificName": "Fungi", "Rank": "kingdom"},
        {"TaxId": "4890", "ScientificName": "Ascomycota", "Rank": "phylum"},
        {"TaxId": "147550", "ScientificName": "Sordariomycetes", "Rank": "class"},
        {"TaxId": "5125", "ScientificName": "Hypocreales", "Rank": "order"},
        {"TaxId": "110618", "ScientificName": "Nectriaceae", "Rank": "family"},
        {"TaxId": "5506", "ScientificName": "Fusarium", "Rank": "genus"},
        {"TaxId": "0", "ScientificName": "unranked thing", "Rank": "no rank"},
    ],
}]


def test_load_name_ids_from_csv(write_file):
    sheet = write_file("tracking.csv", "ID,taxid\nNHM001,5507\n,1\n1002,5506\n")

    # Blank IDs are dropped and numeric IDs come back as strings
    assert load_name_ids(str(sheet), "ID") == ["NHM001", "1002"]


def test_load_name_ids_rejects_unknown_format(write_file):
    sheet = write_file("tracking.tsv", "ID\nNHM001\n")

    with pytest.raises(ValueError, match="Unsupported tracking sheet format"):
        load_name_ids(str(sheet), "ID")


def test_cleanup_temp_dir(tmp_path):
    temp = tmp_path / "tmp"
    (temp / "nested").mkdir(parents=True)

    cleanup_temp_dir(str(temp))
    assert not temp.exists()
    # Cleaning up a directory that is already gone is not an error
    cleanup_temp_dir(str(temp))


@pytest.fixture
def fake_entrez(monkeypatch):
    """Replace Entrez.efetch/read with a scripted sequence of responses.

    Each item in `responses` is either a record to return or an exception to
    raise. Returns a dict recording how many fetches and sleeps happened.
    """
    calls = {"efetch": 0, "sleep": 0}
    responses = []

    def efetch(**kwargs):
        calls["efetch"] += 1
        assert kwargs["db"] == "Taxonomy"
        return object()

    def read(handle):
        response = responses.pop(0)
        if isinstance(response, Exception):
            raise response
        return response

    def sleep(seconds):
        calls["sleep"] += 1

    monkeypatch.setattr(its_fun_tools.Entrez, "efetch", efetch)
    monkeypatch.setattr(its_fun_tools.Entrez, "read", read)
    monkeypatch.setattr(its_fun_tools.time, "sleep", sleep)
    return responses, calls


def test_lineage_keeps_only_standard_ranks(fake_entrez, logger):
    responses, calls = fake_entrez
    responses.append(TAXONOMY_RECORD)

    lineage = get_ncbi_lineage("5507", "test@example.org", logger)

    assert lineage == {
        "kingdom": "Fungi",
        "phylum": "Ascomycota",
        "class": "Sordariomycetes",
        "order": "Hypocreales",
        "family": "Nectriaceae",
        "genus": "Fusarium",
        "species": "Fusarium oxysporum",
    }
    assert calls == {"efetch": 1, "sleep": 0}


def test_lineage_retries_after_transient_error(fake_entrez, logger):
    responses, calls = fake_entrez
    responses.extend([urllib.error.URLError("timed out"), TAXONOMY_RECORD])

    lineage = get_ncbi_lineage("5507", "test@example.org", logger)

    assert lineage["genus"] == "Fusarium"
    assert calls == {"efetch": 2, "sleep": 1}


def test_lineage_gives_up_after_five_attempts(fake_entrez, logger):
    responses, calls = fake_entrez
    # An empty record counts as a failure, like a network error
    responses.extend([[]] * 5)

    with pytest.raises(ValueError, match="No taxonomy record found"):
        get_ncbi_lineage("999999999", "test@example.org", logger)
    assert calls == {"efetch": 5, "sleep": 4}
