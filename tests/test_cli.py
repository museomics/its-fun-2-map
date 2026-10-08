"""Smoke tests for the itsfun-* console scripts.

The command list is read from the installed package metadata, so a command
added to pyproject.toml is tested automatically.
"""

import subprocess
import sys
from importlib.metadata import entry_points
from pathlib import Path

import pytest

import its_fun_2_map

COMMANDS = sorted(
    ep.name
    for ep in entry_points(group="console_scripts")
    if ep.value.startswith("its_fun_2_map.")
)


def run_command(name, *args):
    # Console scripts are installed next to the running interpreter
    executable = Path(sys.executable).parent / name
    return subprocess.run(
        [str(executable), *args], capture_output=True, text=True, timeout=120
    )


def test_all_commands_are_registered():
    assert COMMANDS == sorted([
        "itsfun-assemble",
        "itsfun-blast-parse",
        "itsfun-blast1",
        "itsfun-blast1-parse",
        "itsfun-blast2",
        "itsfun-extract",
        "itsfun-lineage",
        "itsfun-map",
        "itsfun-qc",
        "itsfun-refs",
        "itsfun-summary",
    ])


@pytest.mark.parametrize("command", COMMANDS)
def test_help_exits_cleanly(command):
    result = run_command(command, "--help")
    assert result.returncode == 0, result.stderr
    assert f"usage: {command}" in result.stdout


@pytest.mark.parametrize("command", COMMANDS)
def test_version_matches_package(command):
    result = run_command(command, "--version")
    assert result.returncode == 0, result.stderr
    assert result.stdout.strip() == f"{command} {its_fun_2_map.__version__}"
