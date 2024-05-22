import numpy as np
import math
import sys
import units
import mynumerics as mn
import MMA_administration as MMA

# def load_data

class empty_class:
    pass


class get_data:
    def __init__(self,InputArchive,r_resolution=[True]):
        full_resolution = (r_resolution[0] is True)
        self.omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,
                      MMA.paths['CUPRAD_inputs'] +'/laser_wavelength','N'),'lambdaSI','omegaSI')
        self.k0_wave = 2.0*np.pi/mn.ConvertPhoton(self.omega0,'omegaSI','lambdaSI')
        self.tgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/tgrid'][:]; Nt = len(self.tgrid)
        rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]; Nr = len(rgrid)            
        self.zgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/zgrid'][:]; Nz = len(self.zgrid)            
        if full_resolution:
            kr_step = 1; Nr_max = Nr
        else:
            dr = r_resolution[1]; rmax = r_resolution[2]
            dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
            rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid) 
            
        # self.E_trz = InputArchive['/outputs/output_field'][:,0:Nr_max:kr_step,:Nz] # Arrays may be over-allocated by CUPRAD
        
        self.E_trz = InputArchive[MMA.paths['CUPRAD_outputs'] +'/output_field'][:Nz,:,0:Nr_max:kr_step] # Arrays may be over-allocated by CUPRAD
        
        # hot-fix case of underallocated array, happens rarely
        if (self.E_trz.shape[0] < Nz):
            Nz = self.E_trz.shape[0]
            self.zgrid = self.zgrid[:Nz]
        
        print('oldshape', self.E_trz.shape)
        self.E_trz = self.E_trz.transpose(1,2,0) # hot-fix rearangement due to CUPRAD output
        print('newshape', self.E_trz.shape)

        self.inverse_GV = InputArchive[MMA.paths['CUPRAD_logs'] +'/inverse_group_velocity_SI'][()]
        self.VG_IR = 1.0/self.inverse_GV               
        self.rho0_init = 1e6 * mn.readscalardataset(InputArchive,
                         MMA.paths['CUPRAD_inputs'] +'/calculated/medium_effective_density_of_neutral_molecules','N')
        self.Ip_eV = InputArchive[MMA.paths['CUPRAD_inputs'] +'/ionization_ionization_potential_of_neutral_molecules'][()]
        self.pressure_mbar = 1e3*InputArchive[MMA.paths['CUPRAD_inputs'] +'/medium_pressure_in_bar'][()]; self.pressure_string = "{:.1f}".format(self.pressure_mbar)+' mbar'
        try:
            self.preionisation_ratio = InputArchive['/pre_ionised/initial_electrons_ratio'][()]
        except:
            self.preionisation_ratio = 0
        self.preionisation_string = "{:.1f}".format(100*self.preionisation_ratio) + ' %'
        
        if 'density_mod' in InputArchive[MMA.paths['global_inputs']].keys():
            self.density_mod_profile_mbar = self.pressure_mbar*\
                                            InputArchive[MMA.paths['global_inputs']+'/density_mod/table'][:]
            try: self.density_mod_zgrid = InputArchive[MMA.paths['global_inputs']+'/density_mod/zgrid'][:]
            except: pass
        
            try: self.density_mod_rgrid = InputArchive[MMA.paths['global_inputs']+'/density_mod/rgrid'][:]
            except: pass
        
        self.w0_entry = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_beamwaist_entry',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_beamwaist_entry')
            
        self.rgrid = rgrid
        self.Nr = Nr; self.Nt = Nt; self.Nz = Nz
        
        # Further analyses that may be stored in various directions
        self.Gaussian_focus = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_focus_position_Gaussian',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_focus_position_Gaussian')
        self.Gaussian_focus_string = "{:.1f}".format(1e3*self.Gaussian_focus) + ' mm' # 'z='+"{:.1f}".format(1e3*res.zgrid[0])
        
        self.Intensity_entry = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_intensity_entry',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_intensity_entry')
        self.Intensity_entry_string = "{:.1f}".format(1e-18*self.Intensity_entry) + ' 1e18 W/m2'
        
        # Idealised beam Rayleigh range (defined only if Gaussian waist is available)
        
        self.energy = InputArchive[MMA.paths['CUPRAD'] +'/longstep/energy'][:]
        self.energy_zgrid = InputArchive[MMA.paths['CUPRAD'] +'/longstep/z_buff'][:]
            
        try:
            Gaussian_w0 = mn.readscalardataset(InputArchive, MMA.paths['CUPRAD_inputs'] +'/laser_focus_beamwaist_Gaussian','N')
            self.Gaussian_zR = np.pi*(Gaussian_w0**2)/(1e-2*mn.readscalardataset(InputArchive,MMA.paths['CUPRAD_inputs'] +'/laser_wavelength','N'))
        except:
            self.Gaussian_zR = np.NaN

        try:
            self.Intensity_Gaussian_focus = mn.readscalardataset(InputArchive, MMA.paths['CUPRAD_inputs'] +'/laser_focus_intensity_Gaussian','N')
            self.Intensity_Gaussian_focus_string = "{:.1f}".format(1e-18*self.Intensity_Gaussian_focus) + ' 1e18 W/m2'
        except:
            self.Intensity_Gaussian_focus = np.NaN   
            self.Intensity_Gaussian_focus_string = "xxx"
        
    def vacuum_shift(self,output='replace'):
        E_vac = np.zeros(self.E_trz.shape)   
        for k1 in range(self.Nz):
            delta_z = self.zgrid[k1] # local shift
            delta_t_lab = self.inverse_GV*delta_z # shift to the laboratory frame
            delta_t_vac = delta_t_lab - delta_z/units.c_light # shift to the coordinates moving by c.
            for k2 in range(self.Nr):
                ogrid_nn, FE_s, NF = mn.fft_t_nonorm(self.tgrid, self.E_trz[:,k2,k1]) # transform to omega space        
                FE_s = np.exp(1j*ogrid_nn*delta_t_vac) * FE_s # phase factor        
                tnew, E_s = mn.ifft_t_nonorm(ogrid_nn,FE_s,NF)
                E_vac[:,k2,k1] = E_s.real
        
        if (output == 'replace'):      self.E_trz = E_vac 
        elif (output == 'return'):     return E_vac 
        elif (output == 'add'):        self.E_trz_vac = E_vac 
        else: raise ValueError('wrongly specified output for the vacuum shift.')
        

    def complexify_envel(self,output='return'):
        E_trz_cmplx_envel = np.zeros(self.E_trz.shape,dtype=complex)
        rem_fast_oscillations = np.exp(-1j*self.omega0*self.tgrid)
            
        for k1 in range(self.Nz):
            for k2 in range(self.Nr):
                E_trz_cmplx_envel[:,k2,k1] = rem_fast_oscillations*mn.complexify_fft(self.E_trz[:,k2,k1])
        
        if (output == 'return'):     return E_trz_cmplx_envel
        elif (output == 'add'):      self.E_trz_cmplx_envel = E_trz_cmplx_envel
        else: raise ValueError('wrongly specified output for the vacuum shift.') 
        
    def get_Fluence(self, InputArchive, fluence_source='file'):
        self.Fluence = empty_class()
        if (fluence_source == 'file'):
            self.Fluence.value = InputArchive[MMA.paths['CUPRAD'] +'/longstep/fluence'][:,:]
            self.Fluence.zgrid = InputArchive[MMA.paths['CUPRAD'] +'/longstep/zgrid_analyses2'][:]
            self.Fluence.rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]
            self.Fluence.units = 'C.U.'                
        
        elif (fluence_source == 'computed'):                
            self.Fluence.zgrid = self.zgrid
            self.Fluence.rgrid = self.rgrid
            self.Fluence.value = np.zeros((self.Nr, self.Nz))
            for k1 in range(self.Nz):
                for k2 in range(self.Nr):
                    # self.Fluence.value[k2, k1] = sum(abs(self.E_trz[:, k2, k1])**2)
                    self.Fluence.value[k2, k1] = units.c_light*units.eps0 * np.trapz(abs(self.E_trz[:, k2, k1])**2,self.tgrid)
            self.Fluence.units = 'J/m2'

    def get_plasma(self, InputArchive, r_resolution=[True]): # analogy to the fields
        full_resolution = r_resolution[0]
        self.plasma = empty_class()
        
        self.plasma.tgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/tgrid'][:]; Nt = len(self.tgrid)
        rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]; Nr = len(rgrid)            
        self.plasma.zgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/zgrid'][:]; Nz = len(self.zgrid)            
        if full_resolution:
            kr_step = 1; Nr_max = Nr
        else:
            dr = r_resolution[1]; rmax = r_resolution[2]
            dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
            rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid) 
            
        self.plasma.value_trz = InputArchive[MMA.paths['CUPRAD_outputs'] +'/output_plasma'][:Nz,:,0:Nr_max:kr_step] # [:,0:Nr_max:kr_step,:Nz] # Arrays may be over-allocated by CUPRAD
        self.plasma.value_trz = self.plasma.value_trz.transpose(1,2,0) # (1,2,0) # hot-fix to reshape
        
        self.plasma.rgrid = rgrid
        self.plasma.Nr = Nr; self.plasma.Nt = Nt; self.plasma.Nz = Nz

    def compute_spectrum(self,output='add',compute_dE_domega = False):
        self.ogrid = mn.fft_t(self.tgrid, self.E_trz[:,0,0])[0]
        
        No = len(self.ogrid); Nr = len(self.rgrid); Nz = len(self.zgrid)
                
        FE_trz = np.zeros((No,Nr,Nz),dtype=complex)
        
        for k1 in range(Nz):
            for k2 in range(Nr):
                FE_trz[:,k2,k1] = mn.fft_t(self.tgrid, self.E_trz[:,k2,k1])[1]
  
        
        if compute_dE_domega:
            dE_domega = np.empty((No,Nz))
            
            for k1 in range(Nz):
                for k2 in range(No):
                    dE_domega[k2,k1] = np.trapz(np.abs(FE_trz[k2,:,k1])**2,self.rgrid)
            
            if (output == 'return'):     return FE_trz, dE_domega
            elif (output == 'add'):      self.FE_trz = FE_trz; self.dE_domega = dE_domega
            else: raise ValueError('wrongly specified output for the vacuum shift.') 
        else:
            if (output == 'return'):     return FE_trz
            elif (output == 'add'):      self.FE_trz = FE_trz
            else: raise ValueError('wrongly specified output for the vacuum shift.') 

        
        
def add_print_parameter(parameter,data):
    if (parameter=='pressure'): return data.pressure_string
    elif (parameter=='preionisation'): return data.preionisation_string
    elif (parameter=='focus_in_medium'): return data.Gaussian_focus_string
    elif (parameter=='intensity_entry'): return data.Intensity_entry_string
    elif (parameter=='intensity_Gaussian_focus'): return data.Intensity_Gaussian_focus_string
    else: return ''


def create_param_string(params,data):
    res = ''
    for param in params:
        curr = add_print_parameter(param,data)
        if not(len(curr)==0):
            res = res + ', ' + curr     
    return res   


def measure_beam(grid, beam, measure, *args, measured_axis = 0):
    N0, N1 = beam.shape
    if (measured_axis == 0):
        radius = np.zeros(N1)
        for k1 in range(N1):
            radius[k1] = measure(grid,beam[:,k1],*args)
    if (measured_axis == 1):
        radius = np.zeros(N0)
        for k1 in range(N0):
            radius[k1] = measure(grid,beam[k1,:],*args)
    
    return radius


