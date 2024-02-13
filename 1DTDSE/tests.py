import unittest
import numpy as np
from PythonTDSE import *
import argparse
import sys

class TestTDSE(unittest.TestCase):
    ### test inputs
    @classmethod
    def setUpClass(cls):
        cls.DLL = TDSE_DLL(DLL_path)
        cls.inputs = inputs_def()
        cls.inputs.init_default_inputs(trg_a=1., CV = 1e-15, num_r=8000)
        ### Field amplitude
        E_0 = 0.14
        ### Fundamental frequency
        omega_0 = 0.07
        ### Number of cycles
        Nc = 4
        ### Period
        T = 2*np.pi/omega_0
        ### Pulse length
        T_max = Nc*T
        ### Number of time points
        N_t = int(T_max/cls.inputs.dt) + 1
        ### Temporal grid
        t = np.linspace(0, T_max, N_t)
        ### Sine squared envelope
        sin_2 = lambda t: np.sin(np.pi*t/T_max)**2

        ### Field
        Efield = E_0*sin_2(t)*np.cos(omega_0*t)

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

        cls.DLL.free_mtrx(cls.output.psi, len(cls.wavefunction))
        

    def test_inputs(self):
        ### Load data from the HDF5 archive
        inputs = inputs_def()
        inputs.init_inputs("ionization.h5")
        
        ### Check x_grid
        self.assertTrue(np.allclose(inputs.get_xgrid(), self.inputs.get_xgrid()))
        ### Check GS
        self.assertTrue(np.allclose(inputs.get_GS(), self.inputs.get_GS()))
        ### Check Energy
        self.assertAlmostEqual(inputs.Einit, self.inputs.Einit)

        inputs.delete(self.DLL)

    def test_outputs(self):
        ### Load data from the HDF5 archive
        inputs = inputs_def()
        output = outputs_def()
        inputs.init_inputs("ionization.h5")
        output.load_from_hdf5("ionization.h5")

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
        
        #self.DLL.free_mtrx(output.psi, len(wavefunction))
        inputs.delete(self.DLL)
        output.delete(self.DLL)

    @classmethod
    def tearDownClass(cls):
        cls.inputs.delete(cls.DLL)
        cls.output.delete(cls.DLL)

        
if __name__ == '__main__':
    ap = argparse.ArgumentParser()
    ap.add_argument("-d", "--dll", required=True, help="Path do the C-TDSE DLL.")
    ap.add_argument('unittest_args', nargs='*')
    args = vars(ap.parse_args())

    DLL_path = args['dll']

    args = ap.parse_args()
    sys.argv[1:] = args.unittest_args
    unittest.main()