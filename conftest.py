"""Root pytest configuration and CLI build-path fixtures for MMA test suites."""

import os
import pytest


def pytest_addoption(parser):
    """Register custom CLI options to specify build and installation directories."""
    parser.addoption(
        "--cuprad-build",
        action="store",
        default=os.environ.get("CUPRAD_BUILD", "CUPRAD/build"),
        help="Path to CUPRAD build directory containing cuprad.e and make_start.e (default: $CUPRAD_BUILD or 'CUPRAD/build')",
    )
    parser.addoption(
        "--tdse-build",
        action="store",
        default=os.environ.get("TDSE_1D_BUILD", "1DTDSE/build"),
        help="Path to 1DTDSE build directory containing TDSE.e or libTDSE.so (default: $TDSE_1D_BUILD or '1DTDSE/build')",
    )


@pytest.fixture(scope="session")
def cuprad_build_dir(request):
    """Return the absolute path to the CUPRAD build directory."""
    path = request.config.getoption("--cuprad-build")
    return os.path.abspath(path)


@pytest.fixture(scope="session")
def tdse_build_dir(request):
    """Return the absolute path to the TDSE build directory."""
    path = request.config.getoption("--tdse-build")
    return os.path.abspath(path)

