"""Tests for --fastp_extra_args parsing and validation in fastp_module."""

import argparse

import pytest

from its_fun_2_map.fastp_module import check_fastp_extra_args, parse_fastp_extra_args


@pytest.mark.parametrize("value", [None, ""])
def test_no_extra_args(value):
    assert parse_fastp_extra_args(value) == []


def test_quoted_string_is_split_like_a_shell():
    assert parse_fastp_extra_args('--trim_front1 10 --adapter_sequence "AGAT CGG"') == [
        "--trim_front1", "10", "--adapter_sequence", "AGAT CGG",
    ]


def test_unmanaged_flags_are_accepted():
    assert check_fastp_extra_args(["--trim_front1", "10", "--length_required=50"]) is True


@pytest.mark.parametrize("flag", ["-q", "--qualified_quality_phred", "-w", "--thread"])
def test_managed_flags_are_rejected(flag):
    with pytest.raises(ValueError, match="Conflicting fastp arguments"):
        check_fastp_extra_args([flag, "20"])


def test_managed_flag_in_equals_form_is_rejected():
    with pytest.raises(ValueError, match="--fastp_threads"):
        check_fastp_extra_args(["--thread=8"])


def test_rejection_goes_through_parser_when_given():
    parser = argparse.ArgumentParser(prog="itsfun-qc")

    # parser.error() prints usage and exits with status 2
    with pytest.raises(SystemExit) as excinfo:
        check_fastp_extra_args(["-q", "20"], parser=parser)
    assert excinfo.value.code == 2
