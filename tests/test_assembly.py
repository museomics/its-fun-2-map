"""Tests for scaffold metrics in assembly_module."""

import pytest

from its_fun_2_map.assembly_module import check_assembly_output, get_scaffold_metrics


def test_metrics_for_three_scaffolds(write_file):
    # Lengths 100, 60 (split over two lines) and 40: total 200, half is 100,
    # so the longest scaffold alone reaches it and N50 = 100
    fasta = write_file(
        "scaffolds.fasta",
        ">s1\n" + "A" * 100 + "\n>s2\n" + "C" * 30 + "\n" + "G" * 30 + "\n>s3\n" + "T" * 40 + "\n",
    )

    n, lengths, mean, n50 = get_scaffold_metrics(str(fasta))

    assert n == 3
    assert lengths == [100, 60, 40]
    assert mean == pytest.approx(200 / 3)
    assert n50 == 100


def test_n50_when_longest_is_below_half(write_file):
    # Lengths 50, 40, 30, 30, 30: total 180, half is 90; 50 + 40 = 90 so N50 = 40
    seqs = "".join(f">s{i}\n{'A' * n}\n" for i, n in enumerate([50, 40, 30, 30, 30]))
    fasta = write_file("scaffolds.fasta", seqs)

    assert get_scaffold_metrics(str(fasta))[3] == 40


def test_missing_and_empty_files(tmp_path, write_file):
    empty = write_file("empty.fasta", "")

    assert get_scaffold_metrics(str(tmp_path / "missing.fasta")) == (0, [], 0, 0)
    assert get_scaffold_metrics(str(empty)) == (0, [], 0, 0)


def test_check_assembly_output(tmp_path, write_file):
    fasta = write_file("scaffolds.fasta", ">s1\nACGT\n>s2\nAC\n")
    headers_only = write_file("headers.fasta", ">s1\n>s2\n")

    assert check_assembly_output(str(fasta)) == (True, 2, 3.0, 4)
    assert check_assembly_output(str(headers_only)) == (False, 0, 0, 0)
    assert check_assembly_output(str(tmp_path / "missing.fasta")) == (False, 0, 0, 0)
