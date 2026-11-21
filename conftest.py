"""Root pytest configuration, CLI build-path fixtures, and git-annex dataset handling."""

import os
import shutil
import subprocess
import sys
import pytest

_REPO_ROOT = os.path.dirname(os.path.abspath(__file__))


def pytest_addoption(parser):
    """Register custom CLI options to specify build directories and dataset handling."""
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
    parser.addoption(
        "--datasets",
        action="store",
        choices=("auto", "local", "skip"),
        default="auto",
        help=(
            "Control git-annex test datasets: "
            "'auto' (default) downloads missing datasets required by selected tests via git annex get; "
            "'local' runs only with locally present datasets and skips missing ones; "
            "'skip' deselects all @pytest.mark.dataset tests."
        ),
    )


def is_unfetched_annex_file(path: str) -> bool:
    """Return True if path is a broken git-annex symlink or an unfetched annex pointer file."""
    if os.path.islink(path):
        return not os.path.exists(path)
    if os.path.isfile(path) and os.path.getsize(path) < 1024:
        try:
            with open(path, "rb") as f:
                snippet = f.read(512).decode("utf-8", errors="ignore").strip()
        except OSError:
            return False
        return snippet.startswith("/annex/objects/") or "/.git/annex/objects/" in snippet
    return False


def find_unfetched_annex_files(paths) -> list[str]:
    """Return sorted repo-relative paths for unfetched git-annex files under paths."""
    missing = []
    for raw_path in paths:
        path = os.path.abspath(raw_path)
        if os.path.isdir(path) and not os.path.islink(path):
            for dirpath, _, filenames in os.walk(path):
                for fname in filenames:
                    candidate = os.path.join(dirpath, fname)
                    if is_unfetched_annex_file(candidate):
                        missing.append(os.path.relpath(candidate, _REPO_ROOT).replace(os.sep, "/"))
        elif is_unfetched_annex_file(path):
            missing.append(os.path.relpath(path, _REPO_ROOT).replace(os.sep, "/"))
    return sorted(set(missing))


def _git_env_with_fallback_identity() -> dict[str, str]:
    """Return process env with fallback Git author/committer identity if unset."""
    env = os.environ.copy()
    has_name = bool(env.get("GIT_AUTHOR_NAME") or env.get("GIT_COMMITTER_NAME"))
    has_email = bool(env.get("GIT_AUTHOR_EMAIL") or env.get("GIT_COMMITTER_EMAIL"))

    if not has_name:
        res = subprocess.run(
            ["git", "-C", _REPO_ROOT, "-c", f"safe.directory={_REPO_ROOT}", "config", "--get", "user.name"],
            capture_output=True,
            text=True,
        )
        has_name = res.returncode == 0 and bool(res.stdout.strip())
    if not has_email:
        res = subprocess.run(
            ["git", "-C", _REPO_ROOT, "-c", f"safe.directory={_REPO_ROOT}", "config", "--get", "user.email"],
            capture_output=True,
            text=True,
        )
        has_email = res.returncode == 0 and bool(res.stdout.strip())

    if not has_name:
        env["GIT_AUTHOR_NAME"] = "MMA Test Runner"
        env["GIT_COMMITTER_NAME"] = "MMA Test Runner"
    if not has_email:
        env["GIT_AUTHOR_EMAIL"] = "mma-tests@localhost"
        env["GIT_COMMITTER_EMAIL"] = "mma-tests@localhost"
    return env


def fetch_annex_files(rel_paths: list[str]) -> None:
    """Download missing git-annex files from configured remotes."""
    if not rel_paths:
        return

    has_git_annex = shutil.which("git-annex") is not None
    if not has_git_annex:
        probe = subprocess.run(["git", "annex", "version"], capture_output=True, text=True)
        has_git_annex = probe.returncode == 0
    if not has_git_annex:
        raise pytest.UsageError(
            "git-annex is required to download missing test datasets "
            f"({', '.join(rel_paths)}). Install it via `pip install git-annex`, "
            "or run pytest with `--datasets=local`, `--datasets=skip`, or `-m \"not dataset\"`."
        )

    sys.stderr.write(
        f"[git-annex] Downloading {len(rel_paths)} missing test dataset file(s): "
        f"{', '.join(rel_paths)}\n"
    )
    sys.stderr.flush()

    res = subprocess.run(
        ["git", "-C", _REPO_ROOT, "-c", f"safe.directory={_REPO_ROOT}", "annex", "get", "--", *rel_paths],
        capture_output=True,
        text=True,
        env=_git_env_with_fallback_identity(),
    )
    if res.returncode != 0:
        raise pytest.UsageError(
            "Failed to download missing test datasets with `git annex get`.\n"
            f"Files: {', '.join(rel_paths)}\n"
            f"Stdout:\n{res.stdout}\n"
            f"Stderr:\n{res.stderr}\n"
            "Check network access or ensure the `git-annex` branch is fetched "
            "(`git fetch origin git-annex`), or run with `--datasets=local`, "
            "`--datasets=skip`, or `-m \"not dataset\"`."
        )


def pytest_collection_modifyitems(config, items):
    """Deselect @pytest.mark.dataset tests when --datasets=skip is active."""
    if config.getoption("--datasets") != "skip":
        return

    kept = []
    deselected = []
    for item in items:
        if item.get_closest_marker("dataset") is not None:
            deselected.append(item)
        else:
            kept.append(item)

    if deselected:
        config.hook.pytest_deselected(items=deselected)
        items[:] = kept


def pytest_collection_finish(session):
    """Fetch missing git-annex datasets required by the selected tests in auto mode."""
    if session.config.option.collectonly:
        return
    if session.config.getoption("--datasets") != "auto":
        return

    data_dirs = sorted(
        {
            data_dir
            for item in session.items
            if item.get_closest_marker("dataset") is not None
            for data_dir in [os.path.join(os.path.dirname(str(item.path)), "data")]
            if os.path.isdir(data_dir)
        }
    )
    missing = find_unfetched_annex_files(data_dirs)
    if missing:
        fetch_annex_files(missing)


@pytest.fixture(scope="session")
def ensure_dataset_path(request):
    """Return a helper ensuring a dataset file or directory is fetched per --datasets mode."""
    mode = request.config.getoption("--datasets")

    def _ensure(path: str) -> str:
        missing = find_unfetched_annex_files([path])
        if missing:
            if mode == "local":
                pytest.skip(
                    f"Annexed test dataset not present locally: {', '.join(missing)}. "
                    f"Run without --datasets=local or run `git annex get {' '.join(missing)}`."
                )
            if mode == "auto":
                fetch_annex_files(missing)
        if not os.path.exists(path):
            pytest.skip(f"Test dataset not found at {path}.")
        return path

    return _ensure


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
