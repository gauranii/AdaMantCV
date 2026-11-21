# -----------------------------------------------------------------------------
# File:    python/tests/conftest.py
# Paper:   Gauran, Xue, Chen, Ombao, Yu (2022). Front. Neurosci. 16:836100.
# Purpose: Settings of the pytest run of python/tests/: makes the package
#          `adamant` importable, names the two markers and skips the tests
#          they mark when their condition does not hold.
#            against_r  the test runs R (Rscript) and compares with it, bit
#                       for bit; skipped when Rscript is not on the PATH.
#            full       the long comparisons of test_glmnet.py (the full
#                       penalty grids of Figure 1); run only with --full or
#                       ADAMANT_FULL_TESTS=1.
# Inputs:  None.
# Outputs: None.
# Usage:   python -m pytest python/tests            # inside `nix develop`
#          python -m pytest python/tests --full
#          python -m pytest python/tests -m "not against_r"
# -----------------------------------------------------------------------------
import os
import shutil
import sys

import pytest

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))


def pytest_addoption(parser):
    """Add the option --full."""
    parser.addoption("--full", action="store_true", default=False,
                     help="also run the long comparisons with R (marker full)")


def pytest_configure(config):
    """Name the markers."""
    config.addinivalue_line("markers", "against_r: runs R and compares with it bit for bit")
    config.addinivalue_line("markers", "full: long comparison, run with --full")


def pytest_collection_modifyitems(config, items):
    """Skip the tests whose marker asks for something this run does not have."""
    full = config.getoption("--full") or os.environ.get("ADAMANT_FULL_TESTS") == "1"
    has_r = shutil.which("Rscript") is not None
    for item in items:
        if "full" in item.keywords and not full:
            item.add_marker(pytest.mark.skip(reason="long comparison: run with --full"))
        if "against_r" in item.keywords and not has_r:
            item.add_marker(pytest.mark.skip(reason="Rscript is not on the PATH"))


@pytest.fixture(autouse=True)
def close_figures():
    """Close the figures a test opened."""
    yield
    if "matplotlib.pyplot" in sys.modules:
        sys.modules["matplotlib.pyplot"].close("all")
