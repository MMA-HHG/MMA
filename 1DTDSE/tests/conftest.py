"""Pytest configuration and fixtures for the 1DTDSE test suite."""

import os
import sys
import pytest

# Ensure local modules are accessible
current_dir = os.path.dirname(os.path.abspath(__file__))
mma_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
for p in ["shared_python", "1DTDSE", "CUPRAD/python", "Hankel"]:
    full_path = os.path.join(mma_root, p)
    if full_path not in sys.path:
        sys.path.insert(0, full_path)

TDSE_DATA_DIR = os.path.join(current_dir, "data")
IONIZATION_H5 = os.path.join(TDSE_DATA_DIR, "ionization.h5")


@pytest.fixture(scope="session")
def tdse_data_dir():
    """Return path to 1DTDSE test data directory."""
    return TDSE_DATA_DIR


@pytest.fixture(scope="session")
def ionization_h5():
    """Return path to ionization.h5 reference dataset."""
    if os.path.exists(IONIZATION_H5):
        return IONIZATION_H5
    fallback = os.path.join(current_dir, "ionization.h5")
    if os.path.exists(fallback):
        return fallback
    raise FileNotFoundError(f"Reference ionization.h5 not found at {IONIZATION_H5} or {fallback}")

