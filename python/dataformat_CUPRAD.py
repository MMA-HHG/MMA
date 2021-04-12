import numpy as np
import math
import sys
import units
import mynumerics as mn

# def load_data

class get_data:
    def __init__(self,InputArchive,r_resolution=[True]):
        full_resolution = r_resolution[0]
        self.omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
        self.k0_wave = 2.0*np.pi/mn.ConvertPhoton(self.omega0,'omegaSI','lambdaSI')
        self.tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(self.tgrid)
        rgrid = InputArchive['/outputs/rgrid'][:]; Nr = len(rgrid)            
        self.zgrid = InputArchive['/outputs/zgrid'][:]; Nz = len(self.zgrid)            
        if full_resolution:
            kr_step = 1; Nr_max = Nr
        else:
            dr = r_resolution[1]; rmax = r_resolution[2]
            dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
            rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid) 
            
        self.E_trz = InputArchive['/outputs/output_field'][:,0:Nr_max:kr_step,:Nz] # Arrays may be over-allocated by CUPRAD

        self.inverse_GV = InputArchive['/logs/inverse_group_velocity_SI'][()]
        self.VG_IR = 1.0/self.inverse_GV               
        self.rho0_init = 1e6 * mn.readscalardataset(InputArchive, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
        self.Ip_eV = InputArchive['/inputs/ionization_ionization_potential_of_neutral_molecules'][()]
        self.pressure_mbar = 1e3*InputArchive['/inputs/medium_pressure_in_bar'][()]; self.pressure_string = "{:.1f}".format(self.pressure_mbar)+' mbar'
        self.preionisation_ratio = InputArchive['/pre_ionised/initial_electrons_ratio'][()]; self.preionisation_string = "{:.1f}".format(100*self.preionisation_ratio) + '%'
            
        self.rgrid = rgrid
        self.Nr = Nr; self.Nt = Nt; self.Nz = Nz