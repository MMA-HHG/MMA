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
            
        self.rgrid = rgrid
        self.Nr = Nr; self.Nt = Nt; self.Nz = Nz