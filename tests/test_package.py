"""Package-level checks: metadata, importability and bundled R scripts."""

import importlib
from importlib.metadata import version
from importlib.resources import files

import pytest

import its_fun_2_map

MODULES = [
    "UNITEd",
    "assembly_module",
    "blast_output_parser",
    "blast_round1",
    "blast_round1_parser",
    "blast_round2",
    "fastp_module",
    "its_a_summary_compiler",
    "its_fun_tools",
    "its_primer_binding",
    "mapping_module",
    "pull_ncbi_lineage",
]


def test_installed_version_matches_dunder_version():
    assert version("its-fun-2-map") == its_fun_2_map.__version__


@pytest.mark.parametrize("module", MODULES)
def test_module_imports(module):
    importlib.import_module(f"its_fun_2_map.{module}")


@pytest.mark.parametrize("script", ["parse_fastp_json.R", "its_decision_making.R"])
def test_r_scripts_are_bundled(script):
    resource = files("its_fun_2_map") / script
    assert resource.is_file()
    assert resource.read_text().strip()


def test_parse_fastp_json_does_not_install_packages():
    # R dependencies come from the conda environment, never installed mid-run
    text = (files("its_fun_2_map") / "parse_fastp_json.R").read_text()
    assert "install.packages" not in text
