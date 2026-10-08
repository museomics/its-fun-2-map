"""Shared fixtures for the its-fun-2-map test suite.

All test data is generated on the fly in pytest's tmp_path, so the suite needs
no external data files, network access or bioinformatics tools.
"""

import logging

import pytest


@pytest.fixture
def logger():
    """A named logger for functions that require one. pytest's caplog still
    captures its records."""
    return logging.getLogger("its_fun_2_map.tests")


@pytest.fixture
def write_file(tmp_path):
    """Write text to a file under tmp_path and return its Path."""

    def _write(name, content):
        path = tmp_path / name
        path.parent.mkdir(parents=True, exist_ok=True)
        path.write_text(content)
        return path

    return _write
