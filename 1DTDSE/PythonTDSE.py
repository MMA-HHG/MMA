"""
Python TDSE
===========
"""

from ctypes import *
import h5py


### Define structures
class sin2_definition(Structure):
    _fields_ = [
        ("oc", c_double),
        ("o", c_double),
        ("A0", c_double),
        ("nc", c_double),
        ("phi", c_double),
        ("phi0", c_double),
        ("ti", c_double),
        ("E0", c_double)
    ]

class trap_definition(Structure):
    _fields_ = [
        ("omega", c_double),
        ("E0", c_double),
        ("phi", c_double),
        ("ton", c_double),
        ("toff", c_double),
        ("nc", c_int)
    ]


class Esin2_definition(Structure):
    _fields_ = [
        ("oc", c_double),
        ("o", c_double),
        ("A0", c_double),
        ("nc", c_double),
        ("phi", c_double),
        ("phi0", c_double),
        ("ti", c_double),
        ("E0", c_double)
    ]

class flattop1_definition(Structure):
    _fields_ = [
        ("ton", c_double),
        ("toff", c_double),
        ("T", c_double),
        ("o", c_double),
        ("phi", c_double),
        ("A", c_double),
        ("ti", c_double)
    ]

class flattop1chirp_definition(Structure):
    _fields_ = [
        ("ton", c_double),
        ("toff", c_double),
        ("T", c_double),
        ("o", c_double),
        ("phi", c_double),
        ("A", c_double),
        ("ti", c_double),
        ("b", c_double),
        ("c", c_double)
    ]
    
### Field structure
class Efield_var(Structure):
    _fields_ = [
        ("fieldtype", c_int),
        ("fieldtype2", c_int),
        ("tgrid", POINTER(c_double)),
        ("Field", POINTER(c_double)),
        ("dt", c_double),
        ("omega", c_double),
        ("E0", c_double),
        ("phi", c_double),
        ("ton", c_double),
        ("toff", c_double),
        ("Nt", c_int),
        ("Nftlt1", c_int),
        ("Nftlt1ch", c_int),
        ("Nsin2", c_int),
        ("NEsin2", c_int),
        ("nc", c_int),
        ("trap", trap_definition),
        ("sin2", POINTER(sin2_definition)),
        ("Esin2", POINTER(Esin2_definition)),
        ("flt1", POINTER(flattop1_definition)),
        ("flt1ch", POINTER(flattop1chirp_definition)),
    ]

class trg_def(Structure):
    _fields_ = [
        ("a", c_double)
    ]

class analy_def(Structure):
    _fields_ = [
        ("tprint", c_double),
        ("writewft", c_int)
    ]

class output_print_def(Structure):
    _fields_ = [
        ("Efield", c_int),
        ("FEfield", c_int),
        ("sourceterm", c_int),
        ("Fsourceterm", c_int),
        ("FEfieldM2", c_int),
        ("FsourceTermM2", c_int),
        ("PopTot", c_int),
        ("tgrid", c_int),
        ("omegagrid", c_int),
        ("PopInt", c_int),
        ("expval_x", c_int)
    ]

class inputs_def(Structure):
    _fields_ = [
        ("trg", trg_def),
        ("Efield", Efield_var),
        ("Eguess", c_double),
        ("Einit", c_double),
        ("tmin", c_double),
        ("Nt", c_int),
        ("num_t", c_int),
        ("dt", c_double),
        ("num_r", c_int),
        ("num_exp", c_int),
        ("dx", c_double),
        ("psi0", POINTER(c_double)),
        ("psi", POINTER(c_double)),
        ("x", POINTER(c_double)),
        ("ton", c_double),
        ("toff", c_double),
        ("timet", POINTER(c_double)),
        ("dipole", POINTER(c_double)),
        ("gauge", c_int),
        ("transformgauge", c_int),
        ("x_int", c_double),
        ("analy", analy_def),
        ("InterpByDTorNT", c_int),
        ("Ntinterp", c_int),
        ("PrintGaborAndSpectrum", c_int),
        ("PrintOutputMethod", c_int),
        ("textend", c_double),
        ("dtGabor", c_double),
        ("tmin1window", c_double),
        ("tmin2window", c_double),
        ("tmax1window", c_double),
        ("tmax2window", c_double),
        ("a_Gabor", c_double),
        ("omegaMaxGabor", c_double),
        ("Print", output_print_def),
        ("CV", c_double),
        ("precision", c_char * 2)
    ]

    def init_inputs(self, filename):
        with h5py.File(filename, "r") as f:
            self.E_guess = c_double(f["TDSE_inputs/Eguess"][()])
            self.num_r = c_int(f["TDSE_inputs/N_r_grid"][()])
            self.num_exp = c_int(f["TDSE_inputs/N_r_grid_exp"][()])
            self.dx = c_double(f["TDSE_inputs/dx"][()])
            self.InterpByDTorNT = c_int(f["TDSE_inputs/InterpByDTorNT"][()])
            self.dt = c_double(f["TDSE_inputs/dt"][()])
            self.Ntinterp = c_int(f["TDSE_inputs/Ntinterp"][()])
            self.textend = c_double(f["TDSE_inputs/textend"][()])
            self.analy.writewft = c_int(f["TDSE_inputs/analy_writewft"][()])
            self.analy.tprint = c_double(f["TDSE_inputs/analy_tprint"][()])
            self.x_int = c_double(f["TDSE_inputs/x_int"][()])
            self.PrintGaborAndSpectrum = c_int(f["TDSE_inputs/PrintGaborAndSpectrum"][()])
            self.a_Gabor = c_double(f["TDSE_inputs/a_Gabor"][()])
            self.omegaMaxGabor = c_double(f["TDSE_inputs/omegaMaxGabor"][()])
            self.dtGabor = c_double(f["TDSE_inputs/dtGabor"][()])
            self.tmin1window = c_double(f["TDSE_inputs/tmin1window"][()])
            self.tmax1window = c_double(f["TDSE_inputs/tmax1window"][()])
            self.tmin2window = c_double(f["TDSE_inputs/tmin2window"][()])
            self.tmax2window = c_double(f["TDSE_inputs/tmax2window"][()])
            self.PrintOutputMethod = c_int(f["TDSE_inputs/PrintOutputMethod"][()])
            self.trg.a = c_double(f["TDSE_inputs/trg_a"][()])
            self.CV = c_double(f["TDSE_inputs/CV_criterion_of_GS"][()])
            self.gauge = c_int(f["TDSE_inputs/gauge_type"][()])

class outputs_def(Structure):
    _fields_ = [
        ("tgrid", POINTER(c_double)),
        ("tgrid_fftw", POINTER(c_double)),
        ("Efield", POINTER(c_double)),
        ("sourceterm", POINTER(c_double)),
        ("omegagrid", POINTER(c_double)),
        ("FEfield", POINTER(POINTER(c_double))),
        ("FEfield_data", POINTER(POINTER(c_double))),
        ("Fsourceterm", POINTER(POINTER(c_double))),
        ("Fsourceterm_data", POINTER(c_double)),
        ("FEfieldM2", POINTER(c_double)),
        ("FsourcetermM2", POINTER(c_double)),
        ("PopTot", POINTER(c_double)),
        ("sourcetermfiltered", POINTER(c_double)),
        ("PopInt", POINTER(c_double)),
        ("expval", POINTER(c_double)),
        ("Nt", c_int),
        ("Nomega", c_int)
    ]

### Structures
inputs = inputs_def()
#outputs = outputs_def()


with h5py.File("1DTDSE/results.h5", "r") as f:
    inputs.E_guess = c_double(f["TDSE_inputs/Eguess"][()])
    inputs.num_r = c_int(f["TDSE_inputs/N_r_grid"][()])
    inputs.num_exp = c_int(f["TDSE_inputs/N_r_grid_exp"][()])
    inputs.dx = c_double(f["TDSE_inputs/dx"][()])
    inputs.InterpByDTorNT = c_int(f["TDSE_inputs/InterpByDTorNT"][()])
    inputs.dt = c_double(f["TDSE_inputs/dt"][()])
    inputs.Ntinterp = c_int(f["TDSE_inputs/Ntinterp"][()])
    inputs.textend = c_double(f["TDSE_inputs/textend"][()])
    inputs.analy.writewft = c_int(f["TDSE_inputs/analy_writewft"][()])
    inputs.analy.tprint = c_double(f["TDSE_inputs/analy_tprint"][()])
    inputs.x_int = c_double(f["TDSE_inputs/x_int"][()])
    inputs.PrintGaborAndSpectrum = c_int(f["TDSE_inputs/PrintGaborAndSpectrum"][()])
    inputs.a_Gabor = c_double(f["TDSE_inputs/a_Gabor"][()])
    inputs.omegaMaxGabor = c_double(f["TDSE_inputs/omegaMaxGabor"][()])
    inputs.dtGabor = c_double(f["TDSE_inputs/dtGabor"][()])
    inputs.tmin1window = c_double(f["TDSE_inputs/tmin1window"][()])
    inputs.tmax1window = c_double(f["TDSE_inputs/tmax1window"][()])
    inputs.tmin2window = c_double(f["TDSE_inputs/tmin2window"][()])
    inputs.tmax2window = c_double(f["TDSE_inputs/tmax2window"][()])
    inputs.PrintOutputMethod = c_int(f["TDSE_inputs/PrintOutputMethod"][()])
    inputs.trg.a = c_double(f["TDSE_inputs/trg_a"][()])
    inputs.CV = c_double(f["TDSE_inputs/CV_criterion_of_GS"][()])
    inputs.gauge = c_int(f["TDSE_inputs/gauge_type"][()])


### Init input structure from the HDF5 file
inputs.init_inputs("1DTDSE/results.h5")

### Load compiled dynamic library
path = "/Users/tadeasnemec/Programming/Git/CUPRAD_TDSE_Hankel/1DTDSE/singleTDSE.so"
DLL_Func = CDLL(path)

### Init ground state
init_grid = DLL_Func.Initialise_grid_and_ground_state
init_grid.restype = None
init_grid.argtypes = [POINTER(inputs_def)]
init_grid(byref(inputs))


### Run TDSE from Python
call1DTDSE = DLL_Func.call1DTDSE
call1DTDSE.restype = outputs_def
call1DTDSE.argtypes = [POINTER(inputs_def)]
outputs = call1DTDSE(byref(inputs))

#DLL_Func = CDLL(DLL)
#DLL = "/Users/tadeasnemec/Programming/Git/CUPRAD_TDSE_Hankel/1DTDSE/tools.so"
#print(DLL_Func.test(pointer(inputs)))

#inputs.num_r


#test = DLL_Func.test
#test.restype = c_double
#print(test(byref(inputs)))

