"""Tests for small helpers in its_a_summary_compiler."""

import pandas as pd
import pytest

from its_fun_2_map.its_a_summary_compiler import find_id_column, n_taxa


@pytest.mark.parametrize("value, expected", [
    ("Fusarium;Fusarium;Aspergillus", 2),
    ("Fusarium", 1),
    ("", 0),
    (float("nan"), 0),
])
def test_n_taxa_counts_unique(value, expected):
    assert n_taxa(value) == expected


@pytest.mark.parametrize("column", ["ID", "sample_id", "SampleID", "SAMPLE"])
def test_find_id_column_accepts_variants(column):
    df = pd.DataFrame({column: ["NHM001"], "other": [1]})

    assert find_id_column(df, "file.csv") == column


def test_find_id_column_raises_when_missing():
    with pytest.raises(ValueError, match="No ID column found in file.csv"):
        find_id_column(pd.DataFrame({"name": ["x"]}), "file.csv")
