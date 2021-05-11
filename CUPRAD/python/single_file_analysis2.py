import numpy as np
import os
import time
# import multiprocessing as mp
import shutil
import h5py
import sys
sys.path.append('D:\git\python_modules')
import units
import mynumerics as mn
import HHG
# import mynumerics as mn
import matplotlib.pyplot as plt
import re
import glob
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
from contextlib import ExitStack
import dataformat_CUPRAD as dfC

# There may be a large memory overload in this implementation, it's for testing putposes

# =============================================================================
# Inputs of the script
arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    results_path = os.path.join("D:\data", "Discharges")
    results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "foci")
    results_path = os.path.join("D:\data", "Discharges", "f_scan")
    results_path = os.path.join("D:\data", "Discharges", "TDSE", "t6")

    
cwd = os.getcwd()

vacuum_frame = True

Horders = [15, 17, 19, 21, 23]

gas_type = 'Kr'
XUV_table_type = 'NIST' # {Henke, NIST}


tlim = [-60.0,60.0]
t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase

t_probe1 = 1e-15*np.asanyarray([-10.0, 0.0, 10.0])

OutPath = 'outputs'

full_resolution = False
rmax = 200e-6 # only for analyses
dr = rmax/200.0

Lcoh_saturation = 0.01
Lcoh_zero = 0.0

file = 'results_1.h5' # 'results_Ar_vac.h5', 'Ar_vac_long.h5' 'results_3.h5' 'results_1.h5'

# file = 'results_f_half.h5' # 'results_Ar_vac.h5', 'Ar_vac_long.h5' 'results_3.h5'


# =============================================================================
# The body of the script
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)


file_path = os.path.join(results_path,file)
print('processing:', file_path)             
with h5py.File(file_path, 'r') as InputArchive:
    # load data
    res = dfC.get_data(InputArchive, r_resolution = [full_resolution, dr, rmax])
            
    
    # shift to the vacuum frame
    if vacuum_frame:
        res.vacuum_shift()
    
    
    ## Beamshape
    # Get the fluence
    
    
    fluence_source = 'computed' #'computed'
    
    # Fluence is either computed from Efield or taken from the file
    if (fluence_source == 'file'):
        Fluence = InputArchive['/longstep/fluence'][:,:]
        zgrid_Fluence = InputArchive['/longstep/zgrid_analyses2'][:]
        rgrid_fluence = InputArchive['/outputs/rgrid'][:]
        Fluence_units = 'SI'                
        
    elif (fluence_source == 'computed'):                
        zgrid_Fluence = res.zgrid
        rgrid_fluence = res.rgrid
        Fluence = np.zeros((res.Nr, res.Nz))
        for k1 in range(res.Nz):
            for k2 in range(res.Nr):
                Fluence[k2, k1] = sum(abs(res.E_trz[:, k2, k1])**2)
        Fluence_units = 'arb.u.'
 
 
# get radii
radius_FWHMs = dfC.measure_beam(rgrid_fluence, Fluence, mn.measure_beam_FWHM_zeromax) 
radius_inv_e2 = dfC.measure_beam(rgrid_fluence, Fluence, mn.measure_beam_max_ratio_zeromax, np.exp(-2.0) )    
radius_E_alpha = dfC.measure_beam(rgrid_fluence, Fluence, mn.measure_beam_E_alpha_zeromax, 0.76 )    
radius_RMS = dfC.measure_beam(rgrid_fluence, Fluence, mn.measure_beam_RMS )   
 
fig = plt.figure()
plt.pcolor(1e3*zgrid_Fluence, 1e6*rgrid_fluence, Fluence, cmap='plasma',shading='nearest') # vmin=0.75*np.max(Fluence)
plt.plot(1e3*zgrid_Fluence,1e6*radius_FWHMs)
           
plt.xlabel('z [mm]')
plt.ylabel('r [mum]')
plt.title('Fluence ['+Fluence_units+']')
plt.colorbar()
plt.savefig('Fluence.png', dpi = 600)
if showplots: plt.show()
plt.close(fig)

fig = plt.figure()
plt.plot(1e3*zgrid_Fluence,1e6*radius_FWHMs)
plt.show()
plt.close(fig)

fig = plt.figure()
plt.plot(1e3*zgrid_Fluence,1e6*radius_inv_e2)
plt.show()
plt.close(fig)

fig = plt.figure()
plt.plot(1e3*zgrid_Fluence,1e6*radius_FWHMs)
plt.plot(1e3*zgrid_Fluence,1e6*radius_RMS)
plt.plot(1e3*zgrid_Fluence,1e6*radius_E_alpha)
# plt.plot(1e3*zgrid_Fluence,1e6*radius_inv_e2)
plt.show()
plt.close(fig)

# # Get intensity, envelope & phase
# E_trz_cmplx_envel = np.zeros(E_trz.shape,dtype=complex)
# rem_fast_oscillations = np.exp(-1j*omega0*tgrid)

# for k2 in range(Nz):
#     for k3 in range(Nr):
#         E_trz_cmplx_envel[:,k3,k2] = rem_fast_oscillations*mn.complexify_fft(E_trz[:,k3,k2])

# k_t = mn.FindInterval(tgrid, t_fix)

# phase = np.angle(E_trz_cmplx_envel[k_t,:,:]) # ordering (r,z)
# Intens = mn.FieldToIntensitySI(abs(E_trz_cmplx_envel))

# grad_z_I = np.gradient(Intens[k_t,:,:],zgrid,axis=1,edge_order=2)

# grad_z_phase = np.zeros(phase.shape)
# for k1 in range(Nr):
#     phase_rfix = np.unwrap(phase[k1,:])
#     grad_z_phase[k1,:] = np.gradient(phase_rfix,zgrid,edge_order=2)
        
# # Get XUV
# # the additional term for the phase derivative is given by k0*(nXUV-1), k0 is the vacuum wavenumber
# nXUV = [];
# NH = len(Horders)
# for k1 in range(NH):
#     q = Horders[k1]
#     f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*omega0, 'omegaSI', 'eV'))
#     nXUV.append(1.0 - rho0_init*units.r_electron_classical*(mn.ConvertPhoton(q*omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))



# # Get FSPA 
# with h5py.File(os.path.join(cwd,'FSPA_tables_Krypton_test.h5'),'r') as h5_FSPA_tables:
#     interp_FSPA_short = HHG.FSPA.get_dphase(h5_FSPA_tables,'Igrid','Hgrid','short/dphi')


# FSPA_alphas = []; grad_z_phase_FSPA = []
# for k1 in range(NH):
#     q = Horders[k1]
#     FSPA_alphas.append(-interp_FSPA_short[q](Intens[k_t,:,:]/units.INTENSITYau))
#     grad_z_phase_FSPA.append(FSPA_alphas[k1]*grad_z_I /units.INTENSITYau)


# # make the maps
# k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')


# fig1, ax1 = plt.subplots()
# fig2, ax2 = plt.subplots()
# fig3, ax3 = plt.subplots()

# for k1 in range(NH):
#     q = Horders[k1]   
    
#     # onax
#     ax1.plot(1e3*zgrid,q*(grad_z_phase[0,:] + k0_wave*(nXUV[k1]-1)),label='H'+str(q))
#     ax1.set_xlabel('z [mm]'); ax1.set_ylabel('dPhi/dz [1/m]'); ax1.set_title('beam phase')      
    
#     ax2.plot(1e3*zgrid,grad_z_phase_FSPA[k1][0,:],label='H'+str(q))
#     ax2.set_xlabel('z [mm]'); ax2.set_ylabel('dPhi/dz [1/m]'); ax2.set_title('FSPA (atom)') 
    
#     ax3.plot(1e3*zgrid,
#              q*(grad_z_phase[0,:] + k0_wave*(nXUV[k1]-1)) + grad_z_phase_FSPA[k1][0,:],
#              label='H'+str(q))
#     ax3.set_xlabel('z [mm]'); ax3.set_ylabel('dPhi/dz [1/m]'); ax3.set_title('full phase') 
    
#     fig4, ax4 = plt.subplots()
    
#     map1 = ax4.pcolor(1e3*zgrid, 1e6*rgrid,
#                       q*(grad_z_phase + k0_wave*(nXUV[k1]-1))+grad_z_phase_FSPA[k1],
#                          shading='auto')
#     ax4.set_xlabel('z [mm]'); ax4.set_ylabel('r [mum]'); ax4.set_title('dPhi/dz, full, H'+str(q)) 
#     fig4.colorbar(map1)
#     fig4.savefig('dPhi_dz_H'+str(q)+'.png', dpi = 600)
    
#     fig5, ax5 = plt.subplots()
#     dum = abs( 1.0/ (q*(grad_z_phase + k0_wave*(nXUV[k1]-1))+grad_z_phase_FSPA[k1] ) )
#     if (np.max(dum) > Lcoh_saturation):
#         map1 = ax5.pcolor(1e3*zgrid, 1e6*rgrid, dum, shading='auto', vmax=Lcoh_saturation)
#     else:
#         map1 = ax5.pcolor(1e3*zgrid, 1e6*rgrid, dum, shading='auto')
    
#     ax5.set_xlabel('z [mm]'); ax5.set_ylabel('r [mum]'); ax5.set_title('Lcoh, H'+str(q))         
#     fig5.colorbar(map1)
#     fig5.savefig('Lcoh_H'+str(q)+'.png', dpi = 600)
    

# ax1.legend(loc='best')
# ax2.legend(loc='best')
# ax3.legend(loc='best')

# fig1.savefig('phase_onax_beam.png', dpi = 600)
# fig2.savefig('phase_onax_FSPA.png', dpi = 600)
# fig3.savefig('phase_onax_full.png', dpi = 600)

# fig6, ax6 = plt.subplots()
# map1 = ax6.pcolor(1e3*zgrid, 1e6*rgrid, Intens[k_t,:,:], shading='auto')
# ax6.set_xlabel('z [mm]'); ax6.set_ylabel('r [mum]'); ax6.set_title('Intensity')
# fig6.colorbar(map1) 
# fig6.savefig('Intensity.png', dpi = 600)

# # Cut-off map
# Cutoff = HHG.ComputeCutoff(Intens[k_t,:,:]/units.INTENSITYau,
#                            mn.ConvertPhoton(omega0,'omegaSI','omegaau'),
#                            mn.ConvertPhoton(Ip_eV,'eV','omegaau')
#                            )[1]
# fig7, ax7 = plt.subplots()
# map1 = ax7.pcolor(1e3*zgrid, 1e6*rgrid, Cutoff, shading='auto')
# ax7.set_xlabel('z [mm]'); ax7.set_ylabel('r [mum]'); ax7.set_title('Cutoff')
# fig7.colorbar(map1) 
# fig7.savefig('Cutoff.png', dpi = 600)
    
os.chdir(cwd)


