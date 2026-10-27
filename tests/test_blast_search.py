"""Tests for the BLAST search commands, itsfun-blast1 and itsfun-blast2.

blast_task and the database check are faked out, so no BLAST+ or seqkit
binaries are needed. The tests check that a failed BLAST run makes the
command exit non-zero, so run_pipeline.sh stops instead of reporting a FAIL.
"""

import subprocess

import pytest

from its_fun_2_map import blast_round1, blast_round2


def failing_blast(*args):
    raise subprocess.CalledProcessError(2, ["blastn"])


def passing_blast(*args):
    return None


@pytest.fixture
def blast_inputs(tmp_path, write_file, monkeypatch):
    """One sample with assembled scaffolds and a round 1 hit."""
    write_file("tracking.csv", "ID\ntest1\n")
    write_file("assemblies/test1.spades.out/scaffolds.fasta", ">NODE_1\nACGT\n")
    write_file("round1/test1_blast_results.tsv",
               "NODE_1\thit\t99\t100\t0\t0\t1\t100\t1\t100\t1e-50\t200\n")
    for module in (blast_round1, blast_round2):
        monkeypatch.setattr(module, "ensure_blast_db", lambda *args, **kwargs: None)
    return tmp_path


def common_args(path, n):
    return [
        "--database_file", str(path / "db.fasta"),
        "--query_dir", str(path / "assemblies"),
        "--output_dir", str(path / f"out{n}"),
        "--tracking_sheet", str(path / "tracking.csv"),
        "--column_name", "ID",
        "--log_file", str(path / f"blast{n}.log"),
    ]


@pytest.mark.parametrize("blast, exit_code", [(failing_blast, 1), (passing_blast, 0)])
def test_round1_exit_code(blast_inputs, monkeypatch, blast, exit_code):
    monkeypatch.setattr(blast_round1, "blast_task", blast)

    assert blast_round1.main(common_args(blast_inputs, 1)) == exit_code


@pytest.mark.parametrize("blast, exit_code", [(failing_blast, 1), (passing_blast, 0)])
def test_round2_exit_code(blast_inputs, monkeypatch, blast, exit_code):
    monkeypatch.setattr(blast_round2, "blast_task", blast)
    # The seqkit extraction step is not under test
    monkeypatch.setattr(blast_round2.subprocess, "run", lambda *args, **kwargs: None)
    args = common_args(blast_inputs, 2) + ["--blast_dir", str(blast_inputs / "round1")]

    assert blast_round2.main(args) == exit_code