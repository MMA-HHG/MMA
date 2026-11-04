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
        help="Path to 1DTDSE build directory or DLL (default: $TDSE_1D_BUILD or '1DTDSE/build')",
    )
    parser.addoption(
        "--tdse-dll",
        action="store",
        default=os.environ.get("TDSE_DLL", None),
        help="Explicit path to CTDSE dynamic library (libsingleTDSE.so / libsingleTDSE.dylib). Overrides --tdse-build if set.",
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
    abs_path = os.path.abspath(path)
    if os.path.isfile(abs_path):
        return os.path.dirname(abs_path)
    return abs_path


@pytest.fixture(scope="session")
def tdse_dll_path(request):
    """Return the absolute path to the CTDSE dynamic library."""
    dll_opt = request.config.getoption("--tdse-dll")
    if dll_opt:
        return os.path.abspath(dll_opt)

    build_path = request.config.getoption("--tdse-build")
    abs_path = os.path.abspath(build_path)
    if os.path.isfile(abs_path):
        return abs_path

    # Try common library filenames in the build directory
    for candidate in ["libsingleTDSE.dylib", "libsingleTDSE.so"]:
        full_candidate = os.path.join(abs_path, candidate)
        if os.path.isfile(full_candidate):
            return full_candidate

    return os.path.join(abs_path, "libsingleTDSE.so")

