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
# import mynumerics as mn
import matplotlib.pyplot as plt
import re
import glob
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
from contextlib import ExitStack

# There may be a large memory overload in this implementation, it's for testing putposes

# =============================================================================
# Inputs of the script

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")
cwd = os.getcwd()
# os.chdir(results_path)
files = ['results_1.h5', 'results_4.h5', 'results_7.h5', 'results_10.h5', 'results_13.h5']

# files = ['results_2.h5', 'results_17.h5']
# files = ['results_1.h5', 'results_2.h5', 'results_3.h5']
files = ['results_1.h5', 'results_16.h5']


# labels = ['p=15 mbar', 'p=35 mbar'], ['Pi=0 %', 'Pi=4 %','Pi=8 %', 'Pi=12 %', 'Pi=16 %'], ['I0=1e14 W/cm2', 'I0=1.75e14 W/cm2', 'I0=2.5e14 W/cm2']
labels = ['a','b','c','d','e','f']
# labels = ['I0=1e14 W/cm2', 'I0=1.75e14 W/cm2', 'I0=2.5e14 W/cm2']
linestyles = ['-','--','-.',':','-','--']
outgraph_name = 'Field_compare'
# os.chdir(cwd)

# files = ['results_1.h5']

out_h5name = 'analyses.h5'

Horders = [19, 21, 23, 25, 27]
q = Horders[2]

tlim = [-60.0,60.0]

OutPath = 'outputs'


# =============================================================================
# The body of the script
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)

Nfiles = len(files)
with h5py.File(out_h5name,'w') as OutFile:
    with ExitStack() as stack:
        
        InArch = [stack.enter_context(h5py.File(os.path.join(results_path,fname), 'r')) for fname in files]
        
        # shifts:
        Efield_onaxis_s=[]
        Efield_onaxis_s_XUV=[]
        
        Efield_onaxis_cmplx_envel = []
        Efield_onaxis_cmplx_envel_XUV = []
        
        phase_onaxis_map = []
        phase_onaxis_map_XUV = []
        dPhi_dz_map = []
        dPhi_dz_map_XUV = []
        
        Efield_s=[]
        Efield_s_XUV=[]
        
        zgrid = []
        tgrid = []
        inverse_GV = []
        pressure = []
        VF_IR = []
        nIR = []
        nXUV = []
        
        dens = []
        
        for k1 in range(Nfiles):
            # dum = InArch[k1]['/outputs/output_field'][:,0,:]
            Efield_onaxis_s.append(InArch[k1]['/outputs/output_field'][:,0,:])
            Efield_onaxis_s_XUV.append(np.zeros(Efield_onaxis_s[k1].shape))
            
            # Efield_s.append(InArch[k1]['/outputs/output_field'][:])
            # Efield_s_XUV.append(np.zeros(Efield_s[k1].shape))
            
            
            zgrid.append(InArch[k1]['/outputs/zgrid'][:])
            tgrid.append(InArch[k1]['/outputs/tgrid'][:])
            omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InArch[0],'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
            inverse_GV.append(InArch[k1]['/logs/inverse_group_velocity_SI'][()])
            pressure.append(InArch[k1]['/inputs/medium_pressure_in_bar'][()])
            
            rho0_atm = 1e6 * mn.readscalardataset(InArch[k1], '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
            
            dens.append(rho0_atm)
            
            nIR.append(IR_index.getsusc('Ar', mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')))
            nIR[k1] = np.sqrt(1.0 + pressure[k1]*nIR[k1])
            VF_IR.append(units.c_light/nIR[k1]) # phase velocity of IR
            
            f1 = XUV_index.getf('Ar', mn.ConvertPhoton(q*omega0, 'omegaSI', 'eV'))[0]
            nXUV.append(1.0 - rho0_atm*units.r_electron_classical*(mn.ConvertPhoton(q*omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))
            VF_XUV = units.c_light/nXUV[k1] # phase velocity of XUV
            
            for k2 in range(len(zgrid[k1])):
                ogrid_nn, FE_s, NF = mn.fft_t_nonorm(tgrid[k1], Efield_onaxis_s[k1][:,k2]) # transform to omega space
                
                delta_z = zgrid[k1][k2] # local shift
                delta_t_lab = inverse_GV[k1]*delta_z # shift to the laboratory frame
                delta_t = delta_t_lab - delta_z/units.c_light # shift to the coordinates moving by c.
                delta_t_XUV = delta_t_lab - delta_z/VF_XUV # shift to the coordinates moving by XUV phase.
                
                FE_s_XUV = np.exp(1j*ogrid_nn*delta_t_XUV) * FE_s # phase factor
                FE_s = np.exp(1j*ogrid_nn*delta_t) * FE_s # phase factor
                
                tnew, E_s = mn.ifft_t_nonorm(ogrid_nn,FE_s,NF)
                tnew, E_s_XUV = mn.ifft_t_nonorm(ogrid_nn,FE_s_XUV,NF)
                
                Efield_onaxis_s[k1][:,k2] = E_s.real
                Efield_onaxis_s_XUV[k1][:,k2] = E_s_XUV.real
                
            ## Complexify the field using fft -> ifft & remove fast oscillations
            Efield_onaxis_cmplx_envel.append(np.zeros(Efield_onaxis_s[k1].shape, dtype=complex))
            Efield_onaxis_cmplx_envel_XUV.append(np.zeros(Efield_onaxis_s[k1].shape, dtype=complex))
            
            rem_fast_oscillations = np.exp(-1j*omega0*tgrid[k1])
            for k2 in range(len(zgrid[k1])):
                Efield_onaxis_cmplx_envel[k1][:,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_onaxis_s[k1][:,k2])
                Efield_onaxis_cmplx_envel_XUV[k1][:,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_onaxis_s_XUV[k1][:,k2])
                
            phase_onaxis_map.append(np.squeeze(np.angle(Efield_onaxis_cmplx_envel[k1])))
            phase_onaxis_map_XUV.append(np.squeeze(np.angle(Efield_onaxis_cmplx_envel_XUV[k1])))
            
            dPhi_dz_map.append(np.zeros(phase_onaxis_map[k1].shape))
            dPhi_dz_map_XUV.append(np.zeros(phase_onaxis_map_XUV[k1].shape))
            
            for k2 in range(len(tgrid[k1])):
                phase_t = np.unwrap(phase_onaxis_map[k1][k2,:])
                phase_t_XUV = np.unwrap(phase_onaxis_map_XUV[k1][k2,:])
                
                dPhi_dz_map[k1][k2, 0] = (phase_t[1] - phase_t[0]) / (zgrid[k1][1] - zgrid[k1][0])
                dPhi_dz_map[k1][k2, -1] = (phase_t[-1] - phase_t[-2]) / (zgrid[k1][-1] - zgrid[k1][-2])
                
                dPhi_dz_map_XUV[k1][k2, 0] = (phase_t_XUV[1] - phase_t_XUV[0]) / (zgrid[k1][1] - zgrid[k1][0])
                dPhi_dz_map_XUV[k1][k2, -1] = (phase_t_XUV[-1] - phase_t_XUV[-2]) / (zgrid[k1][-1] - zgrid[k1][-2])
                for k3 in range(1,len(zgrid[k1])-1):
                    dPhi_dz_map[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],phase_t)
                    dPhi_dz_map_XUV[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],phase_t_XUV)
                    # Efield_onaxis_cmplx_envel[k1][:,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_s[:,k2])
                
    # Field in the vacuum frame
    # plt.plot(1e15*tgrid, Efield_onaxis_s[:,0], linewidth=0.2, label='z=0')
    # plt.plot(1e15*tgrid, Efield_onaxis_s[:,Nz//2], linewidth=0.2, label='z=0.5zmax')
    
    for k1 in range(Nfiles):
        plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, label=labels[k1])
        plt.plot(1e15*tgrid[k1], Efield_onaxis_s_XUV[k1][:,-1], linewidth=0.2, label=labels[k1])
    
    # for k1 in range(Nfiles):
    #     plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, linestyle=linestyles[k1], label=labels[k1])
    # plt.plot(1e15*tgrid[1], Efield_onaxis_s[1][:,-1], linewidth=0.2, linestyle='--', label=labels[1])
    plt.xlim(tlim)
    plt.legend(loc='best')
    plt.xlabel('t [fs]')
    plt.ylabel('E [V/m]')
    plt.savefig(outgraph_name+'.png', dpi = 600)
    plt.show()
    
    print(zgrid[0][-1])
    print(zgrid[1][-1])
    print(1e15*(zgrid[0][-1]-zgrid[1][-1])*(1/units.c_light-1/VF_IR[0]))
    
    dPhi_dz_map_XUV2 = [] 
    for k1 in range(Nfiles):
        plt.plot(zgrid[k1], dPhi_dz_map[k1][len(tgrid[k1])//2,:], linewidth=0.2, label=labels[k1])
        plt.plot(zgrid[k1], dPhi_dz_map_XUV[k1][len(tgrid[k1])//2,:], linewidth=0.2, label=labels[k1])
        
        k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
        dPhi_dz_map_XUV2.append(dPhi_dz_map[k1] + k0_wave*(nXUV[k1]-1))
        plt.plot(zgrid[k1], dPhi_dz_map_XUV2[k1][len(tgrid[k1])//2,:], linestyle = '--', linewidth=0.2, label=labels[k1])
    
    # for k1 in range(Nfiles):
    #     plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, linestyle=linestyles[k1], label=labels[k1])
    # plt.plot(1e15*tgrid[1], Efield_onaxis_s[1][:,-1], linewidth=0.2, linestyle='--', label=labels[1])
    plt.legend(loc='best')
    plt.xlabel('z [m]')
    plt.ylabel('Lcoh [m]')
    plt.savefig('Lcoh.png', dpi = 600)
    plt.show()
    
    # print(zgrid[0][-1])
    # print(zgrid[1][-1])
    # print(1e15*(zgrid[0][-1]-zgrid[1][-1])*(1/units.c_light-1/VF_IR[0]))
            
        
        

os.chdir(cwd)
# print('done')
