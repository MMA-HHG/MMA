"""Pytest suite for CTDSE / 1DTDSE dynamic library and Python interface."""

from collections import namedtuple
import multiprocessing
import pickle
from typing import Tuple

import numpy as np
import pytest

from PythonTDSE import (
    TDSE_DLL,
    InitializedStructureError,
    PythonOwnedError,
    byref,
    c_double,
    c_int,
    inputs_def,
    outputs_def,
)

TDSESimulation = namedtuple(
    "TDSESimulation", ["dll", "inputs", "output", "wavefunction"]
)


def set_time_and_field() -> Tuple[np.ndarray, np.ndarray]:
    """Generate analytical laser field temporal grid and waveform."""
    # Field amplitude
    E_0 = 0.14
    # Fundamental frequency
    omega_0 = 0.07
    # Number of cycles
    Nc = 4
    # Period
    T = 2 * np.pi / omega_0
    # Pulse length
    T_max = Nc * T
    # Timestep
    dt = 0.25
    # Number of time points
    N_t = int(T_max / dt) + 1
    # Temporal grid
    t = np.linspace(0, T_max, N_t)
    # Sine squared envelope
    sin_2 = lambda t: np.sin(np.pi * t / T_max) ** 2

    # Field
    Efield = E_0 * sin_2(t) * np.cos(omega_0 * t)

    return t, Efield


def _multiprocessing_error_callback(e):
    """Callback for error logging during multiprocessing pool execution."""
    print(e)


@pytest.fixture(scope="module")
def tdse_sim(tdse_dll_path):
    """Initialize TDSE inputs, ground state, and run baseline 1DTDSE simulation."""
    dll = TDSE_DLL(tdse_dll_path)

    inputs = inputs_def()
    inputs.init_default_inputs(trg_a=1.0, CV=1e-15, num_r=8000)

    t, efield = set_time_and_field()

    ### Init variables
    inputs.init_time_and_field(t=t, E=efield)
    ### Set writing true
    inputs.analy.writewft = c_int(1)
    ### Set wavefunction writing each 10 au in time
    inputs.analy.tprint = c_double(10.0)

    ### Init ground state
    dll.init_GS(inputs)

    ### Define outputs
    output = outputs_def()

    ### Compute TDSE
    dll.call1DTDSE(inputs, output)

    wavefunction = output.get_wavefunction(inputs, grids=False)
    dll.free_mtrx(byref(output.psi), len(wavefunction))

    return TDSESimulation(
        dll=dll, inputs=inputs, output=output, wavefunction=wavefunction
    )


def test_inputs(tdse_sim, ionization_h5):
    """Verify input grid, ground state wavefunction, and initial energy against reference."""
    ### Load data from the HDF5 archive
    inputs = inputs_def()
    inputs.load_from_hdf5(ionization_h5)

    assert inputs._python_owned

    ### Check x_grid
    assert np.allclose(inputs.get_xgrid(), tdse_sim.inputs.get_xgrid())
    ### Check GS
    assert np.allclose(inputs.get_GS(), tdse_sim.inputs.get_GS())
    ### Check Energy
    assert np.isclose(inputs.Einit, tdse_sim.inputs.Einit)


def test_outputs(tdse_sim, ionization_h5):
    """Verify TDSE computed outputs against reference archive."""
    # Wavefunction assertions from baseline simulation
    assert tdse_sim.output._has_wavefunction
    assert len(tdse_sim.wavefunction) == tdse_sim.output._len_wavefunction
    assert not bool(tdse_sim.output.psi)

    ### Load data from the HDF5 archive
    inputs = inputs_def()
    output = outputs_def()
    inputs.load_from_hdf5(ionization_h5)
    output.load_from_hdf5(ionization_h5)

    assert output._python_owned

    ### Check output length
    assert len(output.get_sourceterm()) == len(tdse_sim.output.get_sourceterm())
    ### Check fields
    assert np.allclose(output.get_Efield(), tdse_sim.output.get_Efield())
    ### Check source term
    assert np.allclose(output.get_sourceterm(), tdse_sim.output.get_sourceterm())
    assert np.allclose(output.get_Fsourceterm(), tdse_sim.output.get_Fsourceterm())
    ### Check expectation value of x
    assert np.allclose(output.get_expval(), tdse_sim.output.get_expval())
    ### Check population of GS
    assert np.allclose(output.get_PopTot(), tdse_sim.output.get_PopTot())
    assert np.allclose(output.get_PopInt(), tdse_sim.output.get_PopInt())
    ### Check wavefunction
    wavefunction = output.get_wavefunction(inputs, grids=False)
    assert np.allclose(wavefunction, tdse_sim.wavefunction)


def test_delete(tdse_sim):
    """Verify C pointer management and Python memory deallocation on outputs."""
    output = outputs_def()

    i, o = tdse_sim.dll.call1DTDSE(tdse_sim.inputs, output)

    assert all([i == tdse_sim.inputs, o == output])

    assert not output._python_owned

    try:
        del output
        assert True
    except Exception:
        assert False


def test_list_comprehension(tdse_sim):
    """Verify batched initialization and execution using list comprehensions."""
    N = 2

    t, Efield = set_time_and_field()

    ### Init variables
    inputs = [inputs_def() for _ in range(N)]
    for inp in inputs:
        inp.init_default_inputs(trg_a=1.0, CV=1e-15, num_r=8000)
        inp.init_time_and_field(t=t, E=Efield)
        tdse_sim.dll.init_GS(inp)

    outputs = [outputs_def() for _ in range(N)]

    _ = [tdse_sim.dll.call1DTDSE(inp, out) for (inp, out) in zip(inputs, outputs)]


def test_pickling(tdse_sim):
    """Verify serialization and deserialization of inputs and outputs via pickle."""
    # Pickle and unpickle inputs
    pkl_inputs = pickle.loads(pickle.dumps(tdse_sim.inputs))
    assert pkl_inputs._python_owned
    assert np.allclose(pkl_inputs.get_xgrid(), tdse_sim.inputs.get_xgrid())
    assert np.allclose(pkl_inputs.get_GS(), tdse_sim.inputs.get_GS())
    assert np.allclose(pkl_inputs.get_tgrid(), tdse_sim.inputs.get_tgrid())
    assert np.allclose(pkl_inputs.get_Efield(), tdse_sim.inputs.get_Efield())

    # Pickle and unpickle output
    pkl_output = pickle.loads(pickle.dumps(tdse_sim.output))
    assert pkl_output._python_owned
    assert np.allclose(pkl_output.get_tgrid(), tdse_sim.output.get_tgrid())
    assert np.allclose(pkl_output.get_Efield(), tdse_sim.output.get_Efield())
    assert np.allclose(pkl_output.get_sourceterm(), tdse_sim.output.get_sourceterm())
    assert np.allclose(pkl_output.get_PopTot(), tdse_sim.output.get_PopTot())
    assert np.allclose(pkl_output.get_PopInt(), tdse_sim.output.get_PopInt())


def test_multicore_TDSE(tdse_sim):
    """Verify concurrent TDSE propagation using multiprocessing.Pool."""
    N = 2

    pool = multiprocessing.Pool(N)
    outputs = [outputs_def() for _ in range(N)]
    res = pool.starmap_async(
        tdse_sim.dll.call1DTDSE,
        [(tdse_sim.inputs, o) for o in outputs],
        error_callback=_multiprocessing_error_callback,
    )

    pool.close()
    pool.join()

    for r in res.get():
        for s in r:
            assert s._python_owned

    results = res.get()
    assert np.allclose(tdse_sim.output.get_Efield(), results[0][1].get_Efield())
    assert np.allclose(tdse_sim.output.get_expval(), results[0][1].get_expval())
    assert np.allclose(tdse_sim.output.get_Fsourceterm(), results[0][1].get_Fsourceterm())

    assert res.successful()


def test_python_errors(tdse_sim):
    """Verify error checking for double-initialization and ownership flags."""
    inputs = inputs_def()
    inputs.init_default_inputs(trg_a=1.0, CV=1e-15, num_r=8000)

    inputs = tdse_sim.dll.init_GS(inputs)

    with pytest.raises(InitializedStructureError):
        tdse_sim.dll.init_GS(inputs)

    inputs._python_owned = True

    with pytest.raises(PythonOwnedError):
        tdse_sim.dll.init_GS(inputs)

    inputs._python_owned = False

    t, Efield = set_time_and_field()

    # Init variables
    inputs.init_time_and_field(t=t, E=Efield)

    with pytest.raises(InitializedStructureError):
        inputs.init_time_and_field(t=t, E=Efield)

    inputs, outputs = tdse_sim.dll.call1DTDSE(inputs, outputs_def())

    with pytest.raises(InitializedStructureError):
        inputs, outputs = tdse_sim.dll.call1DTDSE(inputs, outputs)
