"""Validation suite and CLI tool for CUPRAD (cuprad.e, make_start.e) simulation outputs.
"""

import os
import sys
import shutil
import subprocess
import h5py
import numpy as np
import pytest
import MMA_administration as MMA

# Ensure local modules are accessible
current_dir = os.path.dirname(os.path.abspath(__file__))
mma_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
for p in ["shared_python", "CUPRAD/python", "Hankel", "1DTDSE"]:
    full_path = os.path.join(mma_root, p)
    if full_path not in sys.path:
        sys.path.insert(0, full_path)


def validate_cuprad_output(h5_source, edge_decay_ratio=0.05, verbose=True):
    """Validate structural, numerical, and physical sanity of CUPRAD output."""

    with h5py.File(h5_source, "r") as f:
        report = {}

        # Structural Check: Verify required groups and datasets exist
        required_paths = [
            MMA.paths["CUPRAD_outputs"] + "/zgrid",
            MMA.paths["CUPRAD_outputs"] + "/rgrid",
            MMA.paths["CUPRAD_outputs"] + "/tgrid",
            MMA.paths["CUPRAD_outputs"] + "/output_field",
            MMA.paths["CUPRAD_outputs"] + "/output_plasma",
        ]
        for path in required_paths:
            assert path in f, f"Missing required CUPRAD dataset: {path}"

        zgrid = f[MMA.paths["CUPRAD_outputs"] + "/zgrid"][:]
        rgrid = f[MMA.paths["CUPRAD_outputs"] + "/rgrid"][:]
        tgrid = f[MMA.paths["CUPRAD_outputs"] + "/tgrid"][:]
        field = f[MMA.paths["CUPRAD_outputs"] + "/output_field"][:]
        plasma = f[MMA.paths["CUPRAD_outputs"] + "/output_plasma"][:]

        Nz, Nr, Nt = len(zgrid), len(rgrid), len(tgrid)
        report["grid_sizes"] = {"Nz": Nz, "Nr": Nr, "Nt": Nt}

        # Shape and Alignment Checks
        # Note: In raw CUPRAD HDF5, output_field is stored as (Nz, Nt, Nr) or (Nz, Nr, Nt).
        # Arrays may be over-allocated along the z-axis by CUPRAD (e.g. allocated for full run).
        if field.shape[0] > Nz:
            field = field[:Nz]
            plasma = plasma[:Nz]

        if field.shape == (Nz, Nt, Nr):
            field_zrt = np.transpose(field, (0, 2, 1))
            plasma_zrt = np.transpose(plasma, (0, 2, 1))
        elif field.shape == (Nz, Nr, Nt):
            field_zrt = field
            plasma_zrt = plasma
        else:
            raise AssertionError(
                f"Field array shape {field.shape} does not match grids ({Nz}, {Nr}, {Nt}) or ({Nz}, {Nt}, {Nr})"
            )

        assert field.shape == plasma.shape, (
            f"Shape mismatch between output_field {field.shape} and output_plasma {plasma.shape}"
        )

        # Numerical Integrity: No NaNs or Infs
        assert np.isfinite(field_zrt).all(), "output_field contains NaN or Inf values!"
        assert np.isfinite(plasma_zrt).all(), "output_plasma contains NaN or Inf values!"

        # Physical Sanity: Non-trivial field and non-negative plasma
        peak_field = float(np.max(np.abs(field_zrt)))
        assert peak_field > 0.0, "output_field is completely zero (trivial/empty simulation)!"
        report["peak_field_V_m"] = peak_field

        # Plasma density cannot be negative
        min_plasma = float(np.min(plasma_zrt))
        max_plasma = float(np.max(plasma_zrt))
        assert min_plasma >= -1e-6, f"Unphysical negative plasma density detected: {min_plasma}"
        report["plasma_range"] = (min_plasma, max_plasma)

        # Boundary Decay Check (avoid strong box reflections)
        # Check field amplitude at the largest radial coordinate compared to on-axis peak
        edge_field_peak = float(np.max(np.abs(field_zrt[:, -1, :])))
        decay_ratio = edge_field_peak / peak_field
        report["boundary_edge_decay_ratio"] = decay_ratio
        if decay_ratio > edge_decay_ratio:
            msg = (
                f"Warning: Significant field amplitude at radial boundary ({decay_ratio:.3%} of peak). "
                "Consider increasing radial box window."
            )
            if verbose:
                print(msg)

        # Monotonic or sensible propagation grids
        assert np.all(np.diff(rgrid) > 0), "rgrid must be strictly monotonically increasing"
        assert np.all(np.diff(tgrid) > 0), "tgrid must be strictly monotonically increasing"
        if Nz > 1:
            assert np.all(np.diff(zgrid) >= 0), "zgrid must be non-decreasing"

        report["validation_status"] = "PASSED"
        if verbose:
            print(f"[CUPRAD VALIDATION SUCCESSFUL] Grids: Nz={Nz}, Nr={Nr}, Nt={Nt}, Peak E={peak_field:.2e} V/m")

        return report


def test_make_start_execution(tmp_path, cuprad_build_dir, input_h5, make_start_data_dir):
    """Test running make_start.e with a copy of INPUT.h5 in an isolated temporary directory,
    and verify that all expected output files match the contents of make_start reference directory.
    """
    make_start_bin = os.path.join(cuprad_build_dir, "make_start.e")
    if not os.path.exists(make_start_bin):
        pytest.skip(f"make_start.e binary not found in {cuprad_build_dir}")

    # Copy INPUT.h5 into the temporary working directory (leaves original data/INPUT.h5 untouched)
    target_input = tmp_path / "INPUT.h5"
    shutil.copy2(input_h5, target_input)

    # Run make_start.e passing INPUT.h5 via standard input
    try:
        res = subprocess.run(
            [make_start_bin],
            input="INPUT.h5\n",
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
    except OSError as e:
        pytest.skip(f"Cannot execute make_start.e on this platform: {e}")

    assert res.returncode == 0, f"make_start.e failed:\nStdout: {res.stdout}\nStderr: {res.stderr}"

    # Compare with make_start reference directory: verify all expected files exist
    assert os.path.exists(make_start_data_dir), f"Reference directory {make_start_data_dir} does not exist"
    expected_files = sorted(os.listdir(make_start_data_dir))
    assert len(expected_files) > 0, f"No reference files found in {make_start_data_dir}"

    missing_files = []
    for fname in expected_files:
        produced_file = tmp_path / fname
        if not produced_file.exists():
            missing_files.append(fname)

    assert not missing_files, (
        f"make_start did not produce the following expected files in {tmp_path}: {missing_files}"
    )

    # Verify msg.tmp contains INPUT.h5
    msg_tmp = tmp_path / "msg.tmp"
    with open(msg_tmp, "r") as f:
        assert "INPUT.h5" in f.read().strip()

    # Verify the temporary INPUT.h5 was updated with /CUPRAD/inputs/calculated group
    with h5py.File(target_input, "r") as f:
        assert MMA.paths["CUPRAD_inputs"] + "/calculated" in f
        # Original inputs file in data/ remains untouched
        with h5py.File(input_h5, "r") as f_orig:
            assert MMA.paths["CUPRAD_inputs"] + "/calculated" not in f_orig

    # Verify that tmp_path/INPUT.h5 and reference make_start/INPUT.h5 match closely
    ref_input_h5 = os.path.join(make_start_data_dir, "INPUT.h5")
    assert os.path.exists(ref_input_h5), f"Reference INPUT.h5 not found in {make_start_data_dir}"

    with h5py.File(ref_input_h5, "r") as f_ref, h5py.File(target_input, "r") as f_tmp:
        def compare_datasets(name):
            ref_obj = f_ref[name]
            assert name in f_tmp, f"Dataset or group '{name}' missing from produced INPUT.h5"
            tmp_obj = f_tmp[name]

            if isinstance(ref_obj, h5py.Dataset):
                assert isinstance(tmp_obj, h5py.Dataset), f"'{name}' expected to be dataset in produced INPUT.h5"
                d_ref = ref_obj[()]
                d_tmp = tmp_obj[()]

                if np.issubdtype(d_ref.dtype, np.number):
                    # For field arrays (peak ~ 2.4e10 V/m), double-precision rounding across runs
                    # has absolute noise ~ 1e-4 - 1e-3 (relative diff < 1e-13).
                    peak = float(np.max(np.abs(d_ref))) if d_ref.size > 0 else 0.0
                    atol = max(1e-6 * peak, 1e-3)
                    assert np.allclose(
                        d_ref, d_tmp, rtol=1e-7, atol=atol, equal_nan=True
                    ), f"Numerical mismatch in dataset '{name}': max diff = {np.max(np.abs(d_ref - d_tmp))}"
                else:
                    assert np.array_equal(
                        d_ref, d_tmp
                    ), f"Non-numerical dataset mismatch in '{name}'"

        f_ref.visit(compare_datasets)


@pytest.mark.mpi
def test_cuprad_mpi_execution(tmp_path, cuprad_build_dir, make_start_data_dir, cuprad_data_dir):
    """Test running cuprad.e via MPI with a copy of make_start reference directory,
    and validate the resulting output file using validate_cuprad_output and reference comparison.
    """
    cuprad_bin = os.path.join(cuprad_build_dir, "cuprad.e")
    if not os.path.exists(cuprad_bin):
        pytest.skip(f"cuprad.e binary not found in {cuprad_build_dir}")

    # Determine MPI launcher (mpirun or mpiexec)
    mpirun_bin = shutil.which("mpirun") or shutil.which("mpiexec")
    if not mpirun_bin:
        pytest.skip("MPI runner (mpirun/mpiexec) not available on this system")

    # Determine number of MPI processes
    # Check CUPRAD/MPI related env variables, defaulting to 4 (power of 2 required by CUPRAD)
    env_cores = (
        os.environ.get("NUM_PROC_DEFAULT_CUPRAD")
        or os.environ.get("CUPRAD_NUM_PROC")
        or os.environ.get("CUPRAD_CORES")
        or os.environ.get("NUM_PROC")
        or "4"
    )
    try:
        num_cores = int(env_cores)
    except ValueError:
        num_cores = 4

    # Copy entire make_start directory contents into isolated tmp_path
    assert os.path.exists(make_start_data_dir), f"make_start data directory {make_start_data_dir} does not exist"
    for fname in os.listdir(make_start_data_dir):
        src_path = os.path.join(make_start_data_dir, fname)
        dst_path = tmp_path / fname
        if os.path.isdir(src_path):
            shutil.copytree(src_path, dst_path)
        else:
            shutil.copy2(src_path, dst_path)

    # Run cuprad.e under MPI
    cmd = [mpirun_bin, "-n", str(num_cores), cuprad_bin]
    try:
        res = subprocess.run(
            cmd,
            cwd=str(tmp_path),
            capture_output=True,
            text=True,
        )
    except OSError as e:
        pytest.skip(f"Cannot execute MPI command {cmd}: {e}")

    assert res.returncode == 0, f"cuprad.e MPI execution failed:\nStdout: {res.stdout}\nStderr: {res.stderr}"

    # Output file results
    output_h5 = tmp_path / "INPUT.h5"

    assert output_h5.exists(), f"Expected CUPRAD output file (INPUT.h5) not found in {tmp_path}"

    # Validate the output file using the adapted validation logic
    report = validate_cuprad_output(str(output_h5), verbose=False)
    assert report["validation_status"] == "PASSED"
    assert report["peak_field_V_m"] > 0

    # If permanent reference data exists in data/cuprad, compare outputs
    ref_cuprad_h5 = os.path.join(cuprad_data_dir, "INPUT.h5")
    if os.path.exists(ref_cuprad_h5):
        with h5py.File(ref_cuprad_h5, "r") as f_ref, h5py.File(str(output_h5), "r") as f_out:
            # Check field output
            if (
                MMA.paths["CUPRAD_outputs"] + "/output_field" in f_ref
                and MMA.paths["CUPRAD_outputs"] + "/output_field" in f_out
            ):
                f_ref_field = f_ref[MMA.paths["CUPRAD_outputs"] + "/output_field"][:]
                f_out_field = f_out[MMA.paths["CUPRAD_outputs"] + "/output_field"][:]
                Nz_check = min(f_ref_field.shape[0], f_out_field.shape[0])
                peak = float(np.max(np.abs(f_ref_field[:Nz_check])))
                # Scale absolute tolerance with peak amplitude (e.g. 1e-6 of peak).
                # For fields with peak ~ 2.5e10 V/m, numerical noise floor at the boundary
                # is ~ 100-200 V/m, which is rel_diff ~ 8e-9.
                atol = max(1e-6 * peak, 1e-4)
                assert np.allclose(
                    f_ref_field[:Nz_check], f_out_field[:Nz_check], rtol=1e-5, atol=atol
                ), "Simulated output_field differs from reference data/cuprad/INPUT.h5"
