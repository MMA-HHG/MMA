"""Pytest fixtures and configuration for Hankel test suite."""

import os
import sys
import pytest

# Ensure local modules are accessible
current_dir = os.path.dirname(os.path.abspath(__file__))
mma_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
for p in ["shared_python", "Hankel", "CUPRAD/python", "1DTDSE"]:
    full_path = os.path.join(mma_root, p)
    if full_path not in sys.path:
        sys.path.insert(0, full_path)

HANKEL_DATA_DIR = os.path.join(current_dir, "data")
HANKEL_INPUT_PATH = os.path.join(HANKEL_DATA_DIR, "results.h5")
HANKEL_REF_PATH = os.path.join(HANKEL_DATA_DIR, "results_Hankel.h5")


@pytest.fixture(scope="session")
def hankel_data_dir():
    """Path to Hankel permanent reference test data directory."""
    return HANKEL_DATA_DIR


@pytest.fixture(scope="session")
def hankel_input_h5(ensure_dataset_path):
    """Path to Hankel input dataset (from teach_me_mma.ipynb pipeline)."""
    return ensure_dataset_path(HANKEL_INPUT_PATH)


@pytest.fixture(scope="session")
def hankel_ref_h5(ensure_dataset_path):
    """Path to Hankel reference output file."""
    return ensure_dataset_path(HANKEL_REF_PATH)

