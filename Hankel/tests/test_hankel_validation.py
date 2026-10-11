"""Validation suite and test tools for the Hankel transform module.
"""

import os
import sys
import h5py
import numpy as np
import pytest
import MMA_administration as MMA
import mynumerics as mn
import units

# Ensure local modules are accessible
current_dir = os.path.dirname(os.path.abspath(__file__))
mma_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
for p in ["shared_python", "Hankel", "CUPRAD/python", "1DTDSE"]:
    full_path = os.path.join(mma_root, p)
    if full_path not in sys.path:
        sys.path.insert(0, full_path)

import Hankel_transform as HT


def validate_hankel_output(h5_source, verbose=True):
    """Validate structural, numerical, and physical sanity of Hankel transform output."""
    report = {}
    with h5py.File(h5_source, "r") as f:
        hankel_out_grp = MMA.paths["Hankel_outputs"]
        assert hankel_out_grp in f, f"Missing Hankel outputs group: {hankel_out_grp}"

        required_datasets = [
            hankel_out_grp + "/FF_integrated",
            hankel_out_grp + "/ogrid",
            hankel_out_grp + "/rgrid",
            hankel_out_grp + "/entry_plane_transform",
            hankel_out_grp + "/exit_plane_transform",
        ]
        for dpath in required_datasets:
            assert dpath in f, f"Missing required Hankel dataset: {dpath}"

        ogrid = f[hankel_out_grp + "/ogrid"][:]
        rgrid = f[hankel_out_grp + "/rgrid"][:]
        ff_int = f[hankel_out_grp + "/FF_integrated"][:]
        entry_t = f[hankel_out_grp + "/entry_plane_transform"][:]
        exit_t = f[hankel_out_grp + "/exit_plane_transform"][:]

        Nr, No = len(rgrid), len(ogrid)
        report["grid_sizes"] = {"Nr": Nr, "No": No}

        # Grid Monotonicity Checks
        assert Nr > 1, f"Expected Nr > 1, got {Nr}"
        assert No > 1, f"Expected No > 1, got {No}"
        assert np.all(np.diff(rgrid) > 0), "rgrid must be strictly monotonically increasing"
        assert np.all(np.diff(ogrid) > 0), "ogrid must be strictly monotonically increasing"
        assert rgrid[0] >= 0.0, "rgrid must start at or above 0"

        # Shape Checks: complex field stored as (Nr, No, 2)
        assert ff_int.shape == (Nr, No, 2), (
            f"Shape mismatch for FF_integrated: expected ({Nr}, {No}, 2), got {ff_int.shape}"
        )
        assert entry_t.shape == (Nr, No, 2), (
            f"Shape mismatch for entry_plane_transform: expected ({Nr}, {No}, 2), got {entry_t.shape}"
        )
        assert exit_t.shape == (Nr, No, 2), (
            f"Shape mismatch for exit_plane_transform: expected ({Nr}, {No}, 2), got {exit_t.shape}"
        )

        # Numerical Integrity: No NaNs or Infs
        assert np.isfinite(ff_int).all(), "FF_integrated contains NaN or Inf values!"
        assert np.isfinite(entry_t).all(), "entry_plane_transform contains NaN or Inf values!"
        assert np.isfinite(exit_t).all(), "exit_plane_transform contains NaN or Inf values!"

        # Physical Sanity: Non-trivial field amplitude
        ff_complex = ff_int[..., 0] + 1j * ff_int[..., 1]
        peak_intensity = float(np.max(np.abs(ff_complex) ** 2))
        assert peak_intensity > 0.0, "FF_integrated is completely zero (trivial simulation)!"
        report["peak_intensity"] = peak_intensity

        # Optional cumulative fields
        if hankel_out_grp + "/cumulative_field" in f:
            cum_field = f[hankel_out_grp + "/cumulative_field"][:]
            assert cum_field.ndim == 4, f"cumulative_field expected 4D array, got {cum_field.ndim}D"
            assert cum_field.shape[1:] == (Nr, No, 2), (
                f"cumulative_field spatial-spectral shape mismatch: {cum_field.shape}"
            )
            assert np.isfinite(cum_field).all(), "cumulative_field contains NaN or Inf values!"
            report["cumulative_planes"] = cum_field.shape[0]

        if hankel_out_grp + "/zgrid" in f:
            zgrid = f[hankel_out_grp + "/zgrid"][:]
            assert np.all(np.diff(zgrid) > 0), "Hankel output zgrid must be strictly monotonically increasing"
            report["zgrid_len"] = len(zgrid)

        report["validation_status"] = "PASSED"
        if verbose:
            print(f"[HANKEL VALIDATION PASSED] Screen: Nr={Nr}, No={No}, Peak={peak_intensity:.3e}")

        return report


def test_hankel_reference_data_integrity(hankel_ref_h5):
    """Validate that the permanent reference file results_Hankel.h5 passes all structural and physical checks."""
    report = validate_hankel_output(hankel_ref_h5, verbose=False)
    assert report["validation_status"] == "PASSED"
    assert report["grid_sizes"]["Nr"] == 150
    assert report["grid_sizes"]["No"] == 209
    assert report["peak_intensity"] > 0.0


def test_hankel_input_integrity(hankel_input_h5):
    """Verify that results.h5 contains all required inputs and TDSE dipole outputs for Hankel transform."""
    with h5py.File(hankel_input_h5, "r") as f:
        # Check Hankel inputs
        h_inp = MMA.paths["Hankel_inputs"]
        assert h_inp in f, f"Missing group {h_inp}"
        for key in ["distance_FF", "rmax_FF", "Nr_FF", "Harmonic_range", "Nr_max"]:
            assert key in f[h_inp], f"Missing required Hankel parameter: {key}"

        # Check CTDSE outputs (macroscopic dipoles and grids)
        c_out = MMA.paths["CTDSE_outputs"]
        assert c_out in f, f"Missing CTDSE outputs group: {c_out}"
        for key in ["omegagrid", "rgrid_coarse", "zgrid_coarse", "FSourceTerm"]:
            assert key in f[c_out], f"Missing CTDSE output dataset: {key}"

        zgrid = f[c_out + "/zgrid_coarse"][:]
        rgrid = f[c_out + "/rgrid_coarse"][:]
        ogrid = f[c_out + "/omegagrid"][:]
        fsource = f[c_out + "/FSourceTerm"]

        assert len(zgrid) >= 2, f"Hankel requires at least 2 z-planes, found {len(zgrid)}"
        assert fsource.shape == (len(zgrid), len(rgrid), len(ogrid), 2), (
            f"FSourceTerm shape mismatch: {fsource.shape}"
        )


@pytest.mark.slow
def test_hankel_computation_matches_reference(hankel_input_h5, hankel_ref_h5):
    """Compute Hankel transform directly from input file and verify numerical equivalence with reference output."""
    omega_au2SI = mn.ConvertPhoton(1.0, "omegaau", "omegaSI")

    with h5py.File(hankel_input_h5, "r") as in_f, h5py.File(hankel_ref_h5, "r") as ref_f:
        zgrid = in_f[MMA.paths["CTDSE_outputs"] + "/zgrid_coarse"][:]
        rgrid = in_f[MMA.paths["CTDSE_outputs"] + "/rgrid_coarse"][:]
        ogrid = in_f[MMA.paths["CTDSE_outputs"] + "/omegagrid"][:]

        inp_group = in_f[MMA.paths["Hankel_inputs"]]
        hrange = inp_group["Harmonic_range"][:]
        wavelength = mn.readscalardataset(in_f, MMA.paths["CUPRAD_inputs"] + "/laser_wavelength", "N")
        omega0 = mn.ConvertPhoton(1e-2 * wavelength, "lambdaSI", "omegaau")
        ko_min = mn.FindInterval(ogrid / omega0, hrange[0])
        ko_max = mn.FindInterval(ogrid / omega0, hrange[-1])
        ko_step = mn.readscalardataset(inp_group, "ko_step", "N")
        kr_step = mn.readscalardataset(inp_group, "kr_step", "N")
        kr_max = mn.readscalardataset(inp_group, "Nr_max", "N")
        distance_ff = mn.readscalardataset(inp_group, "distance_FF", "N")
        rmax_ff = mn.readscalardataset(inp_group, "rmax_FF", "N")
        nr_ff = mn.readscalardataset(inp_group, "Nr_FF", "N")

        rgrid_ff = np.linspace(0.0, rmax_ff, nr_ff)

        target = HT.FSources_provider(
            zgrid,
            rgrid,
            omega_au2SI * ogrid,
            h5_handle=in_f,
            h5_path=MMA.paths["CTDSE_outputs"] + "/FSourceTerm",
            data_source="dynamic",
            ko_min=ko_min,
            ko_max=ko_max,
            ko_step=ko_step,
            kr_max=kr_max,
            kr_step=kr_step,
        )

        inv_gv = in_f[MMA.paths["CUPRAD_logs"] + "/inverse_group_velocity_SI"][()]
        eff_ir_index = inv_gv * units.c_light
        pressure = MMA.pressure_constructor(in_f)
        preset_gas = mn.readscalardataset(in_f, MMA.paths["global_inputs"] + "/gas_preset", "S")

        res = HT.Hankel_long(
            target,
            distance_ff,
            rgrid_ff,
            preset_gas=preset_gas,
            pressure=pressure,
            absorption_tables="Henke",
            dispersion_tables="NIST",
            effective_IR_refrective_index=eff_ir_index,
            store_cumulative_result=True,
        )

        # Compare Far-Field Integrated signal
        computed_ff = res.FF_integrated
        ref_ff_raw = ref_f[MMA.paths["Hankel_outputs"] + "/FF_integrated"][:]
        ref_ff_complex = ref_ff_raw[..., 0] + 1j * ref_ff_raw[..., 1]

        assert np.allclose(computed_ff, ref_ff_complex, rtol=1e-5, atol=1e-8), (
            f"Mismatch between computed FF_integrated and reference: max diff = {np.max(np.abs(computed_ff - ref_ff_complex))}"
        )

