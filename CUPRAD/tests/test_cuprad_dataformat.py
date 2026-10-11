import h5py
import numpy as np
import pytest
import dataformat_CUPRAD as dfC
import mynumerics as mn
import HHG
import units
import MMA_administration as MMA


def test_load_cuprad_data(synthetic_cuprad_h5):
    """Test standard loading of CUPRAD simulation data via dataformat_CUPRAD."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        data = dfC.get_data(f)

        assert data.Nz == 3
        assert data.Nr == 32
        assert data.Nt == 64
        assert data.E_zrt.shape == (3, 32, 64)
        assert len(data.zgrid) == 3
        assert len(data.rgrid) == 32
        assert len(data.tgrid) == 64

        # Validate physics metadata
        expected_omega0 = mn.ConvertPhoton(800e-9, "lambdaSI", "omegaSI")
        assert np.isclose(data.omega0, expected_omega0, rtol=1e-4)
        assert np.isclose(data.pressure_mbar, 100.0)
        assert data.gas_type == "Ar"
        assert data.effective_neutral_particle_density > 0


def test_load_plasma(synthetic_cuprad_h5):
    """Test plasma density extraction and physical constraints."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        data = dfC.get_data(f)
        data.get_plasma(f)

        assert hasattr(data, "plasma")
        assert data.plasma.value_zrt.shape == (3, 32, 64)

        # Physical invariants: plasma density cannot be negative or exceed neutral density
        assert np.all(data.plasma.value_zrt >= 0.0)
        assert np.max(data.plasma.value_zrt) <= data.effective_neutral_particle_density


def test_computed_fluence(synthetic_cuprad_h5):
    """Test computed fluence integration over time."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        data = dfC.get_data(f)
        data.get_Fluence(f, fluence_source="computed")

        assert hasattr(data, "Fluence")
        assert data.Fluence.value.shape == (data.Nr, data.Nz)
        assert data.Fluence.units == "J/m2"
        assert np.all(data.Fluence.value >= 0.0)
        # On-axis fluence at focus must be greater than edge fluence
        assert data.Fluence.value[0, 0] > data.Fluence.value[-1, 0]


def test_compute_spectrum(synthetic_cuprad_h5):
    """Test spectrum computation from electric field traces."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        data = dfC.get_data(f)
        FE_zrt, dE_domega = data.compute_spectrum(output="return", compute_dE_domega=True)

        assert FE_zrt.shape[0] == data.Nz
        assert FE_zrt.shape[1] == data.Nr
        assert dE_domega.shape[0] == data.Nz
        assert np.all(dE_domega >= 0.0)

        # Peak frequency should be near omega0
        peak_omega_idx = np.argmax(dE_domega[0, :])
        peak_omega = data.ogrid[peak_omega_idx]
        assert np.isclose(peak_omega, data.omega0, rtol=0.1)


def test_subsampled_radial_resolution(synthetic_cuprad_h5):
    """Test loading CUPRAD data with custom radial downsampling."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        dr_requested = 3e-6
        rmax_requested = 30e-6
        data_downsampled = dfC.get_data(f, r_resolution=[False, dr_requested, rmax_requested])

        assert data_downsampled.rgrid[-1] <= rmax_requested + 1e-12
        assert data_downsampled.Nr < 32
        assert data_downsampled.E_zrt.shape == (3, data_downsampled.Nr, 64)


def test_notebook_analysis_helpers(synthetic_cuprad_h5):
    """Test the cutoff and symmetrization helpers as used in teach_me_mma.ipynb (Cell 36)."""
    with h5py.File(synthetic_cuprad_h5, "r") as f:
        data = dfC.get_data(f)

        # Test radial symmetrization matching teach_me_mma.ipynb:
        # array_xy is transposed to have r along axis 1
        field_tr = data.E_zrt[0, :, :].T
        r_grid_sym, sym_field = mn.symmetrize_y(data.rgrid, field_tr)
        assert sym_field.shape[1] == 2 * data.Nr - 1
        assert sym_field.shape[0] == data.Nt
        assert np.isclose(r_grid_sym[0], -data.rgrid[-1])
        assert np.isclose(r_grid_sym[-1], data.rgrid[-1])

        # Test HHG cutoff mapping
        intensity_si = mn.FieldToIntensitySI(data.E_zrt[0, 0, :])
        intensity_au = intensity_si / units.INTENSITYau
        omega_au = mn.ConvertPhoton(data.omega0, "omegaSI", "omegaau")
        ip_au = mn.ConvertPhoton(data.Ip_eV, "eV", "omegaau")

        cutoffs = HHG.ComputeCutoff(intensity_au, omega_au, ip_au)[1]
        assert np.all(cutoffs >= 1.0)


def test_notebook_inputs_creation(input_h5):
    """Test that input_h5 fixture accurately prepares teach_me_mma.ipynb input schema."""
    with h5py.File(input_h5, "r") as f:
        # Check required groups from inputs_transformer
        assert MMA.paths["global_inputs"] in f
        assert MMA.paths["CUPRAD_inputs"] in f
        assert MMA.paths["CTDSE_inputs"] in f
        assert MMA.paths["Hankel_inputs"] in f

        # Verify key parameters from teach_me_mma.ipynb
        assert f[MMA.paths["global_inputs"] + "/gas_preset"][()].decode() == "Ar"
        assert np.isclose(f[MMA.paths["global_inputs"] + "/medium_pressure_in_bar"][()], 0.1)
        assert f[MMA.paths["CUPRAD_inputs"] + "/numerics_number_of_points_in_r"][()] == 256
        assert f[MMA.paths["CUPRAD_inputs"] + "/numerics_number_of_points_in_t"][()] == 256
        assert f[MMA.paths["Hankel_inputs"] + "/distance_FF"][()] == 1.0


def test_notebook_reference_results_loading(reference_results_h5):
    """Test loading the user's reference CUPRAD simulation results when placed in CUPRAD/tests/data/."""
    with h5py.File(reference_results_h5, "r") as f:
        data = dfC.get_data(f)
        assert data.gas_type == "Ar"
        assert data.Nr == 256
        assert data.Nt == 256
        assert data.E_zrt.shape[1] == 256
        assert data.E_zrt.shape[2] == 256
