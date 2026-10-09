import unittest
import numpy as np
from PythonTDSE import *
import argparse
import sys
import multiprocessing
import pickle
import logging

def configure_logging(verbose: bool) -> None:
    level = logging.DEBUG if verbose else logging.INFO
    logging.basicConfig(level=level, format="[%(levelname)s] %(message)s")

class TestTDSE(unittest.TestCase):
    ### test inputs
    @staticmethod
    def set_time_and_field() -> Tuple[np.ndarray]:
        # Field amplitude
        E_0 = 0.14
        # Fundamental frequency
        omega_0 = 0.07
        # Number of cycles
        Nc = 4
        # Period
        T = 2*np.pi/omega_0
        # Pulse length
        T_max = Nc*T
        # Timestep
        dt = 0.25
        # Number of time points
        N_t = int(T_max/dt) + 1
        # Temporal grid
        t = np.linspace(0, T_max, N_t)
        # Sine squared envelope
        sin_2 = lambda t: np.sin(np.pi*t/T_max)**2

        # Field
        Efield = E_0*sin_2(t)*np.cos(omega_0*t)

        return t, Efield

    @classmethod
    def setUpClass(cls):
        cls.DLL = TDSE_DLL(DLL_path)

        cls.inputs = inputs_def()
        cls.inputs.init_default_inputs(trg_a=1., CV = 1e-15, num_r=8000)

        t, Efield = cls.set_time_and_field()

        ### Init variables
        cls.inputs.init_time_and_field(t = t, E = Efield)
        ### Set writing true
        cls.inputs.analy.writewft = c_int(1)
        ### Set wavefunction writing each 10 au in time
        cls.inputs.analy.tprint = c_double(10.)

        ### Init ground state
        cls.DLL.init_GS(cls.inputs)

        ### Define outputs
        cls.output = outputs_def()

        ### Compute TDSE
        cls.DLL.call1DTDSE(cls.inputs, cls.output)

        cls.wavefunction = cls.output.get_wavefunction(cls.inputs, grids=False)

        assert cls.output._has_wavefunction
        assert len(cls.wavefunction) == cls.output._len_wavefunction

        cls.DLL.free_mtrx(byref(cls.output.psi), len(cls.wavefunction))

        assert not bool(cls.output.psi)

    def test_inputs(self):
        ### Load data from the HDF5 archive
        inputs = inputs_def()
        inputs.load_from_hdf5("tests/ionization.h5")

        self.assertTrue(inputs._python_owned)

        ### Check x_grid
        self.assertTrue(np.allclose(inputs.get_xgrid(), self.inputs.get_xgrid()))
        ### Check GS
        self.assertTrue(np.allclose(inputs.get_GS(), self.inputs.get_GS()))
        ### Check Energy
        self.assertAlmostEqual(inputs.Einit, self.inputs.Einit)

    def test_outputs(self):
        ### Load data from the HDF5 archive
        inputs = inputs_def()
        output = outputs_def()
        inputs.load_from_hdf5("tests/ionization.h5")
        output.load_from_hdf5("tests/ionization.h5")

        self.assertTrue(output._python_owned)

        ### Check output length
        self.assertEqual(len(output.get_sourceterm()), len(self.output.get_sourceterm()))
        ### Check fields
        self.assertTrue(np.allclose(output.get_Efield(), self.output.get_Efield()))
        ### Check source term
        self.assertTrue(np.allclose(output.get_sourceterm(), self.output.get_sourceterm()))
        self.assertTrue(np.allclose(output.get_Fsourceterm(), self.output.get_Fsourceterm()))
        ### Check expectation value of x
        self.assertTrue(np.allclose(output.get_expval(), self.output.get_expval()))
        ### Check population of GS
        self.assertTrue(np.allclose(output.get_PopTot(), self.output.get_PopTot()))
        self.assertTrue(np.allclose(output.get_PopInt(), self.output.get_PopInt()))
        ### Check wavefunction
        wavefunction = output.get_wavefunction(inputs, grids = False)
        self.assertTrue(np.allclose(wavefunction, self.wavefunction))

    def test_delete(self):
        output = outputs_def()

        i, o = self.DLL.call1DTDSE(self.inputs, output)

        self.assertTrue(all([i == self.inputs, o == output]))

        self.assertFalse(output._python_owned)

        try:
            del output
            assert True
        except:
            assert False

    def test_list_comprehension(self):
        N = 2

        inputs = inputs_def()
        inputs.init_default_inputs(trg_a=1., CV = 1e-15, num_r=8000)

        t, Efield = self.set_time_and_field()

        ### Init variables
        inputs = [inputs_def() for i in range(N)]
        for input in inputs:
            input.init_default_inputs(trg_a=1., CV = 1e-15, num_r=8000)
            input.init_time_and_field(t = t, E = Efield)
            self.DLL.init_GS(input)

        outputs = [outputs_def() for i in range(N)]

        _ = [self.DLL.call1DTDSE(input, output) for (input, output) in zip(inputs, outputs)]

    def test_pickling(self):
        # Pickle and unpickle self.inputs
        pkl_inputs = pickle.loads(pickle.dumps(self.inputs))
        self.assertTrue(pkl_inputs._python_owned)
        self.assertTrue(np.allclose(pkl_inputs.get_xgrid(), self.inputs.get_xgrid()))
        self.assertTrue(np.allclose(pkl_inputs.get_GS(), self.inputs.get_GS()))
        self.assertTrue(np.allclose(pkl_inputs.get_tgrid(), self.inputs.get_tgrid()))
        self.assertTrue(np.allclose(pkl_inputs.get_Efield(), self.inputs.get_Efield()))

        # Pickle and unpickle self.output
        pkl_output = pickle.loads(pickle.dumps(self.output))
        self.assertTrue(pkl_output._python_owned)
        self.assertTrue(np.allclose(pkl_output.get_tgrid(), self.output.get_tgrid()))
        self.assertTrue(np.allclose(pkl_output.get_Efield(), self.output.get_Efield()))
        self.assertTrue(np.allclose(pkl_output.get_sourceterm(), self.output.get_sourceterm()))
        self.assertTrue(np.allclose(pkl_output.get_PopTot(), self.output.get_PopTot()))
        self.assertTrue(np.allclose(pkl_output.get_PopInt(), self.output.get_PopInt()))

    def test_multicore_TDSE(self):
        N = 2

        pool = multiprocessing.Pool(N)
        outputs = [outputs_def() for i in range(N)]
        res = pool.starmap_async(
            self.DLL.call1DTDSE,
            [(self.inputs, o) for o in outputs],
            error_callback=self.callback
        )

        pool.close()
        pool.join()

        for r in res.get():
            for s in r:
                self.assertTrue(s._python_owned)

        assert res.successful()

    def test_python_errors(self):
        inputs = inputs_def()
        inputs.init_default_inputs(trg_a=1., CV = 1e-15, num_r=8000)

        inputs = self.DLL.init_GS(inputs)

        with self.assertRaises(InitializedStructureError):
            self.DLL.init_GS(inputs)

        inputs._python_owned = True

        with self.assertRaises(PythonOwnedError):
            self.DLL.init_GS(inputs)

        inputs._python_owned = False

        t, Efield = self.set_time_and_field()

        # Init variables
        inputs.init_time_and_field(t = t, E = Efield)

        with self.assertRaises(InitializedStructureError):
            inputs.init_time_and_field(t = t, E = Efield)

        inputs, outputs = self.DLL.call1DTDSE(inputs, outputs_def())

        with self.assertRaises(InitializedStructureError):
            inputs, outputs = self.DLL.call1DTDSE(inputs, outputs)

    @staticmethod
    def callback(e):
        print(e)
        return

    @classmethod
    def tearDownClass(cls):
        pass

if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument(
        "-d", "--debug",
        action="store_true",
        help="Enable debug logging output"
    )
    ap.add_argument("-l", "--dll", required=True, help="Path do the C-TDSE DLL.")
    ap.add_argument('unittest_args', nargs='*')
    args = vars(ap.parse_args())

    configure_logging(args["debug"])

    DLL_path = args["dll"]

    args = ap.parse_args()
    sys.argv[1:] = args.unittest_args
    unittest.main()