"""Shared fixtures for the its-fun-2-map test suite.

All test data is generated on the fly in pytest's tmp_path, so the suite needs
no external data files, network access or bioinformatics tools.
"""

import logging
import os

import pytest

# Environment as it was before any test module imports rpy2. Starting the
# embedded R rewrites LD_LIBRARY_PATH to R's library dirs (which include
# /usr/lib/x86_64-linux-gnu), so a Python launched afterwards can load the
# system libpython instead of its own and fail to find installed packages.
# conftest.py is imported before test modules are collected, so this copy is
# untouched; subprocesses in the tests should use it.
CLEAN_ENV = os.environ.copy()


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
