import h5py
import numpy as np
import pytest
import os
import sys
import units
import mynumerics as mn
import MMA_administration as MMA
import HHG
from inputs_transformer import add_variables2hdf5

# Ensure local modules are accessible
current_dir = os.path.dirname(os.path.abspath(__file__))
mma_root = os.path.abspath(os.path.join(current_dir, "..", ".."))
for p in ["shared_python", "CUPRAD/python", "Hankel", "1DTDSE"]:
    full_path = os.path.join(mma_root, p)
    if full_path not in sys.path:
        sys.path.insert(0, full_path)


def generate_synthetic_cuprad_data(h5_path, Nz=3, Nr=32, Nt=64):
    """Generate a self-contained, valid synthetic CUPRAD simulation archive.

    The archive contains an analytical Gaussian laser pulse and a synthetic
    plasma density profile adhering to the MMA CUPRAD HDF5 specification.
    """
    wavelength = 800e-9  # m
    omega0 = mn.ConvertPhoton(wavelength, "lambdaSI", "omegaSI")
    w0 = 20e-6  # m
    tau = 10e-15  # s
    I0 = 1e18  # W/m2 (~1e14 W/cm2)
    E0 = np.sqrt(2.0 * I0 / (units.c_light * units.eps0))

    zgrid = np.linspace(0.0, 40e-6, Nz)
    rgrid = np.linspace(0.0, 45e-6, Nr)
    tgrid = np.linspace(-25e-15, 25e-15, Nt)

    # Coordinates for meshgrid: (Nz, Nt, Nr) to match CUPRAD's raw output layout
    Z, T, R = np.meshgrid(zgrid, tgrid, rgrid, indexing="ij")

    # Analytical Gaussian beam: E(z, t, r)
    zR = np.pi * w0**2 / wavelength
    w_z = w0 * np.sqrt(1.0 + (Z / zR)**2)
    envelope = (w0 / w_z) * np.exp(-(R / w_z)**2) * np.exp(-(T / tau)**2)
    output_field = E0 * envelope * np.cos(omega0 * T)

    # Synthetic plasma density (non-negative, delayed buildup following intensity profile)
    rho_neutral = 2.504e25  # m^-3 (effective neutral density)
    ionization_fraction = 0.05 * np.exp(-2.0 * (R / w0)**2) * (1.0 + np.tanh(T / (0.5 * tau)))
    output_plasma = rho_neutral * ionization_fraction

    with h5py.File(h5_path, "w") as f:
        # Global inputs
        f.create_dataset(MMA.paths["global_inputs"] + "/gas_preset", data=np.bytes_("Ar"))
        f.create_dataset(MMA.paths["global_inputs"] + "/medium_pressure_in_bar", data=0.1)

        # CUPRAD inputs
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_wavelength", data=1e2 * wavelength)  # cm
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/medium_pressure_in_bar", data=0.1)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_beamwaist_entry", data=w0)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_pulse_duration_in_1_e_Efield", data=tau)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_focus_position_Gaussian", data=0.0)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_focus_beamwaist_Gaussian", data=w0)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/laser_focus_intensity_Gaussian", data=I0)
        f.create_dataset(MMA.paths["CUPRAD_inputs"] + "/calculated/laser_intensity_entry", data=I0)
        f.create_dataset(
            MMA.paths["CUPRAD_inputs"] + "/ionization_ionization_potential_of_neutral_molecules",
            data=15.76,  # Argon Ip in eV
        )
        f.create_dataset(
            MMA.paths["CUPRAD_inputs"] + "/calculated/medium_effective_density_of_neutral_molecules",
            data=1e-6 * rho_neutral,  # cm^-3 in HDF5
        )

        # CUPRAD longstep diagnostics
        f.create_dataset(MMA.paths["CUPRAD"] + "/longstep/energy", data=np.full(Nz, 1e-4))
        f.create_dataset(MMA.paths["CUPRAD"] + "/longstep/z_buff", data=zgrid)

        # CUPRAD logs
        f.create_dataset(MMA.paths["CUPRAD_logs"] + "/inverse_group_velocity_SI", data=1.0 / units.c_light)

        # CUPRAD outputs
        f.create_dataset(MMA.paths["CUPRAD_outputs"] + "/zgrid", data=zgrid)
        f.create_dataset(MMA.paths["CUPRAD_outputs"] + "/rgrid", data=rgrid)
        f.create_dataset(MMA.paths["CUPRAD_outputs"] + "/tgrid", data=tgrid)
        f.create_dataset(MMA.paths["CUPRAD_outputs"] + "/output_field", data=output_field)
        f.create_dataset(MMA.paths["CUPRAD_outputs"] + "/output_plasma", data=output_plasma)

    return h5_path


@pytest.fixture
def synthetic_cuprad_h5(tmp_path):
    """Fixture providing a temporary synthetic CUPRAD HDF5 archive."""
    file_path = str(tmp_path / "synthetic_cuprad.h5")
    generate_synthetic_cuprad_data(file_path)
    return file_path

def prepare_input(h5_path):
    """Prepare results.h5 inputs using the exact parameters from teach_me_mma.ipynb (Cells 6-24)."""
    gas = "Ar"
    medium_length = 40e-6  # [m]
    medium_pressure = 100e-3  # [mbar]

    ionisation_model = "PPT"
    XUV_dispersion_tables = "NIST"
    XUV_absorption_tables = "Henke"

    laser_wavelength = 800e-9  # [m]
    reference_Gaussian_focus = 0.0
    reference_Gaussian_waist = 20e-6  # [m]
    laser_pulse_duration = 10e-15  # [s]

    peak_intensity = (
        HHG.ComputeInvCutoff_gas(
            20.0, mn.ConvertPhoton(laser_wavelength, "lambdaSI", "omegaau"), gas=gas
        )
        * units.INTENSITYau
    )
    reference_Gaussian_focus_intensity = peak_intensity

    XUV_camera_distance = 1.0  # [m]
    XUV_camera_harmonic_range = np.asarray([14.0, 28.0])
    XUV_camera_radial_range = 0.007  # [m]

    number_of_points_in_r = 256
    number_of_points_in_t = 256
    operators_t = 2
    first_delta_z = 0.01  # [mm]
    phase_threshold_for_decreasing_delta_z = 0.002  # [rad]
    length_of_window_for_r_normalized_to_beamwaist = 4.0
    length_of_window_for_t_normalized_to_pulse_duration = 4.0
    number_of_absorber_points_in_time = 16
    physical_output_distance_for_plasma_and_Efield = 10e-6  # [m]
    output_distance_in_z_steps_for_fluence_and_power = 100
    radius_for_diagnostics = 0.1  # [cm]
    run_time_in_hours = 5.0  # [h]

    kz_step = 1
    kr_step_CTDSE = 4
    Nr_max = 110
    dt_TDSE = 0.25  # [a.u.]
    dx = 0.4  # [a.u.]
    Nx_max = 7000
    x_int = 2.0  # [a.u.]
    absorber_type = 0
    CV_criterion_of_ground_state = 1e-12

    list_of_CTDSE_outputs = [
        "electric field",
        "<dj/dt> (Fourier)",
        "ground-state population (projected)",
        "<x>",
    ]

    store_cumulative_field = True
    kr_step_Hankel = 1
    ko_step = 1
    Nr_max_Hankel_integration = 256
    XUV_camera_number_of_r_points = 150
    Nthreads = int(os.environ.get("NUM_PROC_DOCKER", "1"))

    global_inps = {
        "gas_preset": (np.bytes_(gas), "[-]"),
        "medium_pressure_in_bar": (medium_pressure, "[bar]"),
    }

    cuprad_inps = {
        "laser_wavelength": (1e2 * laser_wavelength, "[cm]"),
        "laser_pulse_duration_in_1_e_Efield": (1e15 * laser_pulse_duration, "[fs]"),
        "laser_focus_intensity_Gaussian": (reference_Gaussian_focus_intensity, "[W/m2]"),
        "laser_focus_beamwaist_Gaussian": (reference_Gaussian_waist, "[m]"),
        "laser_focus_position_Gaussian": (reference_Gaussian_focus, "[m]"),
        "medium_physical_distance_of_propagation": (medium_length, "[m]"),
        "ionization_model": (np.bytes_(ionisation_model), "[-]"),
        "numerics_number_of_points_in_r": (number_of_points_in_r, "[-]"),
        "numerics_number_of_points_in_t": (number_of_points_in_t, "[-]"),
        "numerics_operators_t_t-1": (operators_t, "[-]"),
        "numerics_physical_first_stepwidth": (first_delta_z, "[mm]"),
        "numerics_phase_threshold_for_decreasing_delta_z": (phase_threshold_for_decreasing_delta_z, "[rad]"),
        "numerics_length_of_window_for_r_normalized_to_beamwaist": (
            length_of_window_for_r_normalized_to_beamwaist,
            "[-]",
        ),
        "numerics_length_of_window_for_t_normalized_to_pulse_duration": (
            length_of_window_for_t_normalized_to_pulse_duration,
            "[-]",
        ),
        "numerics_number_of_absorber_points_in_time": (number_of_absorber_points_in_time, "[-]"),
        "numerics_physical_output_distance_for_plasma_and_Efield": (
            physical_output_distance_for_plasma_and_Efield,
            "[m]",
        ),
        "numerics_output_distance_in_z-steps_for_fluence_and_power": (
            output_distance_in_z_steps_for_fluence_and_power,
            "[-]",
        ),
        "numerics_radius_for_diagnostics": (radius_for_diagnostics, "[cm]"),
        "numerics_run_time_in_hours": (run_time_in_hours, "[s]"),
    }

    ctdse_inps = {
        "x_int": (x_int, "[a.u.]"),
        "Nr_max": (Nr_max, "[-]"),
        "kr_step": (kr_step_CTDSE, "[-]"),
        "kz_step": (kz_step, "[-]"),
        "dx": (dx, "[a.u.]"),
        "Nx_max": (Nx_max, "[a.u.]"),
        "dt": (dt_TDSE, "[a.u.]"),
        "CV_criterion_of_GS": (CV_criterion_of_ground_state, "[a.u.]"),
        "absorber_type": (absorber_type, "[-]"),
    }

    ctdse_outputs = {
        "print_Efield": "electric field",
        "print_F_Efield": "electric field (Fourier)",
        "print_Source_Term": "<dj/dt>",
        "print_F_Source_Term": "<dj/dt> (Fourier)",
        "print_GS_population": "ground-state population (projected)",
        "print_integrated_population": "ground-state population (integrated)",
        "print_x_expectation_value": "<x>",
    }

    hankel_inps = {
        "distance_FF": (XUV_camera_distance, "[m]"),
        "rmax_FF": (XUV_camera_radial_range, "[m]"),
        "Nr_FF": (XUV_camera_number_of_r_points, "[-]"),
        "XUV_table_type_dispersion": (np.bytes_(XUV_dispersion_tables), "[-]"),
        "XUV_table_type_absorption": (np.bytes_(XUV_absorption_tables), "[-]"),
        "kr_step": (kr_step_Hankel, "[-]"),
        "ko_step": (ko_step, "[-]"),
        "Nr_max": (Nr_max_Hankel_integration, "[-]"),
        "Harmonic_range": (XUV_camera_harmonic_range, "[harmonic_order]"),
        "store_cumulative_result": (int(store_cumulative_field), "[-]"),
        "Nthreads": (Nthreads, "[-]"),
    }

    with h5py.File(h5_path, "w") as f:
        add_variables2hdf5(
            f,
            global_input_names_to_jupyter_variables=global_inps,
            CUPRAD_names_to_jupyter_variables=cuprad_inps,
            CTDSE_names_to_jupyter_variables=ctdse_inps,
            CTDSE_outputs_to_jupyter_names=ctdse_outputs,
            list_of_CTDSE_outputs=list_of_CTDSE_outputs,
            Hankel_names_to_jupyter_variables=hankel_inps,
        )

    return h5_path


DATA_DIR = os.path.join(current_dir, "data")
MAKE_START_DATA_DIR = os.path.join(DATA_DIR, "make_start")
CUPRAD_DATA_DIR = os.path.join(DATA_DIR, "cuprad")
INPUT_PATH = os.path.join(DATA_DIR, "INPUT.h5")
CUPRAD_REF_INPUT_PATH = os.path.join(CUPRAD_DATA_DIR, "INPUT.h5")


@pytest.fixture(scope="session")
def make_start_data_dir():
    """Fixture providing path to the permanent make_start pre-processed data directory."""
    return MAKE_START_DATA_DIR


@pytest.fixture(scope="session")
def cuprad_data_dir():
    """Fixture providing path to the permanent CUPRAD reference data directory."""
    return CUPRAD_DATA_DIR


@pytest.fixture(scope="session")
def input_h5():
    """Fixture ensuring the notebook inputs file is generated in CUPRAD/tests/data/INPUT.h5."""
    os.makedirs(DATA_DIR, exist_ok=True)
    if not os.path.exists(INPUT_PATH):
        prepare_input(INPUT_PATH)
    return INPUT_PATH


@pytest.fixture
def reference_results_h5():
    """Fixture pointing to the designated reference results file"""
    if os.path.exists(CUPRAD_REF_INPUT_PATH):
        return CUPRAD_REF_INPUT_PATH

    pytest.skip(
        f"Reference CUPRAD file not found at {CUPRAD_REF_INPUT_PATH}."
    )
