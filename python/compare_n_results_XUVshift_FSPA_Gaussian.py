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

with h5py.File('FSPA_tables_Krypton.h5','r') as h5_FSPA_tables:
    interp_FSPA_short = HHG.FSPA.get_dphase(h5_FSPA_tables,'Igrid','Hgrid','short/dphi')

print(interp_FSPA_short[17]([8.47268e-5,9.56486e-5]))

# sys.exit(0)    
    
cwd = os.getcwd()
# os.chdir(results_path)
files = ['results_1.h5', 'results_4.h5', 'results_7.h5', 'results_10.h5', 'results_13.h5']

# files = ['results_2.h5', 'results_17.h5']
# files = ['results_1.h5', 'results_2.h5', 'results_3.h5']
files = ['results_1.h5', 'results_16.h5']
files = ['results_1.h5', 'results_4.h5', 'results_7.h5']

files = ['results_1.h5', 'results_10.h5']
files = ['results_1.h5', 'results_2.h5', 'results_3.h5']

files = ['results_1.h5', 'results_4.h5', 'results_7.h5', 'results_10.h5', 'results_13.h5', 'results_16.h5', 'results_19.h5']

files = ['results_Ar_vac.h5']

files = ['Ar_vac_long.h5']

# files = ['results_12.h5', 'results_13.h5']
# files = ['results_1.h5', 'results_2.h5']

# labels = ['p=15 mbar', 'p=35 mbar'], ['Pi=0 %', 'Pi=4 %','Pi=8 %', 'Pi=12 %', 'Pi=16 %'], ['I0=1e14 W/cm2', 'I0=1.75e14 W/cm2', 'I0=2.5e14 W/cm2']
labels = ['a','b','c','d','e','f']
# labels = ['I0=1e14 W/cm2', 'I0=1.75e14 W/cm2', 'I0=2.5e14 W/cm2']
linestyles = ['-','--','-.',':','-','--']
outgraph_name = 'Field_compare'
# os.chdir(cwd)

# files = ['results_1.h5']

out_h5name = 'analyses.h5'

Horders = [15, 17, 19, 21, 23]

gas_type = 'Kr'
XUV_table_type = 'NIST' # {Henke, NIST}

# q = Horders[2]

tlim = [-60.0,60.0]
t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase

t_probe1 = 1e-15*np.asanyarray([-10.0, 0.0, 10.0])

OutPath = 'outputs'


# =============================================================================
# The body of the script

NH = len(Horders)

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
        
        # test intensities
        Intens_onaxis_envel = []
        Intens_onaxis_envel_XUV = []
        
        dI_dz_map = []
        dI_dz_map_XUV = []
        
        # Copmare to Gaussian
        w0 = []
        zR = []
        zgrid_probe1 = []
        tFWHM = []
        
        
        # Phase maps from the propagation
        phase_onaxis_map = []
        phase_onaxis_map_XUV = []
        dPhi_dz_map = []
        dPhi_dz_map_XUV = []
        
        Efield_s=[]
        Efield_s_XUV=[]
        Efield_probe = []
        
        zgrid = []
        tgrid = []
        rgrid = []
        inverse_GV = []
        pressure = []
        VF_IR = []
        nIR = []
        nXUV = []
        VF_XUV = []
        
        dens = []
        Intensity_entry = []
        Ip_eV = []
        
        for k1 in range(Nfiles):
            # dum = InArch[k1]['/outputs/output_field'][:,0,:]
            
            # Efield_s.append(InArch[k1]['/outputs/output_field'][:])
            # Efield_s_XUV.append(np.zeros(Efield_s[k1].shape))
            
            
            zgrid.append(InArch[k1]['/outputs/zgrid'][:]); Nz_loc=len(zgrid[k1])
            tgrid.append(InArch[k1]['/outputs/tgrid'][:]);  # t_probe1_ind = mn.FindInterval(tgrid[k1], t_probe1)
            rgrid.append(InArch[k1]['/outputs/rgrid'][:]);
            z_probe1_ind = np.asarray([0,Nz_loc//2,Nz_loc-1])
            zgrid_probe1.append(zgrid[k1][z_probe1_ind])
            # zgrid_probe1.append(zgrid[k1])
            
            
            # We arrived to the problem with possible over-allocation of Efield in CUPRAD
            Efield_onaxis_s.append(InArch[k1]['/outputs/output_field'][:,0,:Nz_loc])
            Efield_onaxis_s_XUV.append(np.zeros(Efield_onaxis_s[k1].shape))
            Efield_probe.append(InArch[k1]['/outputs/output_field'][:,:,z_probe1_ind])  
            # Efield_probe.append(InArch[k1]['/outputs/output_field'][:,:,:Nz_loc])
            
            
            omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InArch[0],'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
            inverse_GV.append(InArch[k1]['/logs/inverse_group_velocity_SI'][()])
            pressure.append(InArch[k1]['/inputs/medium_pressure_in_bar'][()])
            
            Intensity_entry.append(InArch[k1]['/inputs/laser_intensity_entry'][()])
            Ip_eV.append(InArch[k1]['/inputs/ionization_ionization_potential_of_neutral_molecules'][()])
            
            rho0_init = 1e6 * mn.readscalardataset(InArch[k1], '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
            
            dens.append(rho0_init)
            
            # Gaussian
            w0.append(1e-2*mn.readscalardataset(InArch[k1],'/inputs/laser_beamwaist','N'))
            zR.append(np.pi * (w0[k1])**2 / mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))
            tFWHM.append(1e-15*mn.readscalardataset(InArch[k1],'/inputs/calculated/laser_pulse_duration_in_FWHM_Intensity','N'))
            # Efocus = np.sqrt(2.0*Intensity_entry[k1]/(units.c_light*units.eps0))
            
            # Efield_onaxis_s[k1] = np.squeeze(
            #                         mn.GaussianBeamEfieldConstruct(np.asarray([0.0]),zgrid[k1],tgrid[k1],Efocus,w0[k1],1e-15*tFWHM,mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))
            #                                 )
            nIR.append(IR_index.getsusc(gas_type, mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')))
            nIR[k1] = np.sqrt(1.0 + pressure[k1]*nIR[k1])
            VF_IR.append(units.c_light/nIR[k1]) # phase velocity of IR
            
            nXUV.append([]); VF_XUV.append([])
            for k2 in range(NH):
                q = Horders[k2]
                f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*omega0, 'omegaSI', 'eV'))
                nXUV[k1].append(1.0 - rho0_init*units.r_electron_classical*(mn.ConvertPhoton(q*omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))
                VF_XUV[k1].append(units.c_light/nXUV[k1][k2]) # phase velocity of XUV
            
            for k2 in range(len(zgrid[k1])):
                ogrid_nn, FE_s, NF = mn.fft_t_nonorm(tgrid[k1], Efield_onaxis_s[k1][:,k2]) # transform to omega space
                
                delta_z = zgrid[k1][k2] # local shift
                delta_t_lab = inverse_GV[k1]*delta_z # shift to the laboratory frame
                delta_t = delta_t_lab - delta_z/units.c_light # shift to the coordinates moving by c.
                delta_t_XUV = delta_t_lab - delta_z/VF_XUV[k1][0] # shift to the coordinates moving by XUV phase.
                
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
                
                dPhi_dz_map[k1][k2, :] = mn.ddx_vec_arb(zgrid[k1],phase_t)
                dPhi_dz_map_XUV[k1][k2, :] = mn.ddx_vec_arb(zgrid[k1],phase_t_XUV)
                
                
                # dPhi_dz_map[k1][k2, 0] = (phase_t[1] - phase_t[0]) / (zgrid[k1][1] - zgrid[k1][0])
                # dPhi_dz_map[k1][k2, -1] = (phase_t[-1] - phase_t[-2]) / (zgrid[k1][-1] - zgrid[k1][-2])
                
                # dPhi_dz_map_XUV[k1][k2, 0] = (phase_t_XUV[1] - phase_t_XUV[0]) / (zgrid[k1][1] - zgrid[k1][0])
                # dPhi_dz_map_XUV[k1][k2, -1] = (phase_t_XUV[-1] - phase_t_XUV[-2]) / (zgrid[k1][-1] - zgrid[k1][-2])
                # for k3 in range(1,len(zgrid[k1])-1):
                #     dPhi_dz_map[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],phase_t)
                #     dPhi_dz_map_XUV[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],phase_t_XUV)
                    # Efield_onaxis_cmplx_envel[k1][:,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_s[:,k2])
                    
            ## Compute the intensity profile to get the atomic phase
            # We have to point out that our actual FSPA gives only dphi/dI! We cannot obtain full phase at the instant
            Efield_onaxis_envel = abs(Efield_onaxis_cmplx_envel[k1])
            Efield_onaxis_envel_XUV = abs(Efield_onaxis_cmplx_envel_XUV[k1])

            Intens_onaxis_envel.append( mn.FieldToIntensitySI(Efield_onaxis_envel) )
            Intens_onaxis_envel_XUV.append( mn.FieldToIntensitySI(Efield_onaxis_envel_XUV) ) 
            
            dI_dz_map.append(np.zeros(Intens_onaxis_envel[k1].shape))
            dI_dz_map_XUV.append(np.zeros(Intens_onaxis_envel_XUV[k1].shape))
            
            for k2 in range(len(tgrid[k1])):
                
                dI_dz_map[k1][k2, :] = mn.ddx_vec_arb(zgrid[k1], Intens_onaxis_envel[k1][k2,:])
                dI_dz_map_XUV[k1][k2, :] = mn.ddx_vec_arb(zgrid[k1], Intens_onaxis_envel_XUV[k1][k2,:])
                
                # dI_dz_map[k1][k2, 0] = (Intens_onaxis_envel[k1][k2,1] - Intens_onaxis_envel[k1][k2,0]) / (zgrid[k1][1] - zgrid[k1][0])
                # dI_dz_map[k1][k2, -1] = (Intens_onaxis_envel[k1][k2,-1] - Intens_onaxis_envel[k1][k2,-2]) / (zgrid[k1][-1] - zgrid[k1][-2])
                
                # dI_dz_map_XUV[k1][k2, 0] = (Intens_onaxis_envel_XUV[k1][k2,1] - Intens_onaxis_envel_XUV[k1][k2,0]) / (zgrid[k1][1] - zgrid[k1][0])
                # dI_dz_map_XUV[k1][k2, -1] = (Intens_onaxis_envel_XUV[k1][k2,-1] - Intens_onaxis_envel_XUV[k1][k2,-2]) / (zgrid[k1][-1] - zgrid[k1][-2])

                
                # for k3 in range(1,len(zgrid[k1])-1):
                #     dI_dz_map[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],Intens_onaxis_envel[k1][k2,:])
                #     dI_dz_map_XUV[k1][k2, k3] = mn.ddx_arb(k3,zgrid[k1],Intens_onaxis_envel_XUV[k1][k2,:])
                #     # Efield_onaxis_cmplx_envel[k1][:,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_s[:,k2])
            
            

    # Intensity compare
    
    fig = plt.figure()
    for k1 in range(Nfiles):
        # plt.plot(1e15*tgrid[k1], Intens_onaxis_envel[k1][:,-1], linewidth=0.2, label=labels[k1])
        # plt.plot(1e15*tgrid[k1], Intens_onaxis_envel_XUV[k1][:,-1], linewidth=0.2, label=labels[k1])
        plt.plot(1e15*tgrid[k1], Intens_onaxis_envel[k1][:,-1], linewidth=0.2)
        plt.plot(1e15*tgrid[k1], Intens_onaxis_envel_XUV[k1][:,-1], linewidth=0.2)
    
    # for k1 in range(Nfiles):
    #     plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, linestyle=linestyles[k1], label=labels[k1])
    # plt.plot(1e15*tgrid[1], Efield_onaxis_s[1][:,-1], linewidth=0.2, linestyle='--', label=labels[1])
    plt.xlim(tlim)
    plt.legend(loc='best')
    plt.xlabel('t [fs]')
    plt.ylabel('I [SI]')
    plt.savefig('Intens_compare'+'.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)

    # FSPA phase
    FSPA_alpha_map = interp_FSPA_short[15](Intens_onaxis_envel_XUV[0]/units.INTENSITYau)
    FSPA_phase_map = np.multiply(Intens_onaxis_envel_XUV[0]/units.INTENSITYau, FSPA_alpha_map)
    dFSPA_phase_map = np.multiply(dI_dz_map_XUV[0]/units.INTENSITYau, FSPA_alpha_map)
 
    
    kt1 = mn.FindInterval(tgrid[0], 1e-15*tlim[0])
    kt2 = mn.FindInterval(tgrid[0], 1e-15*tlim[1])
    
    # plot
    fig = plt.figure()
    plt.pcolor(zgrid[0], tgrid[0][kt1:kt2], dFSPA_phase_map[kt1:kt2,:],shading='auto')
    # plt.xlabel('z [mm]')
    # plt.ylabel('r [mum]')
    plt.title('dphiFSPA/dz')
    # plt.ylim(1e-15*np.asarray(tlim))
    plt.colorbar()
    plt.savefig('dphiFSPA_dz.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)    

    fig = plt.figure()
    plt.pcolor(zgrid[0], tgrid[0][kt1:kt2], dPhi_dz_map[0][kt1:kt2,:],shading='auto')
    # plt.xlabel('z [mm]')
    # plt.ylabel('r [mum]')
    plt.title('dPhi/dz')
    # plt.ylim(1e-15*np.asarray(tlim))
    plt.colorbar()
    plt.savefig('dPhi_dz.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)    
    
    # phase plot
    fig = plt.figure()
    k_t = mn.FindInterval(tgrid[k1], t_fix) 
    plt.plot(zgrid[0], phase_onaxis_map[0][k_t,:], label='IR-phase')
    plt.xlabel('z [m]')
    plt.ylabel('Phi [rad]')
    plt.legend(loc='best')
    plt.title('onaxis')
    plt.savefig('phase_onaxis_test.png', dpi = 600)
    plt.close(fig) 
    
    # linear plots of dPhi/dz
    k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
    for k1 in range(Nfiles):
        k_t = mn.FindInterval(tgrid[k1], t_fix) 
        fig = plt.figure()
        plt.plot(zgrid[k1], dPhi_dz_map[k1][k_t,:], label='only IR')
        for k2 in range(NH):
            dPhi_dz_map_Harm = dPhi_dz_map[k1] + k0_wave*(nXUV[k1][k2]-1)
            plt.plot(zgrid[0], dPhi_dz_map_Harm[k_t,:], label='H'+str(Horders[k2]))
        
        plt.xlabel('z [m]')
        plt.ylabel('dPhi/dz [rad/m]')
        plt.legend(loc='best')
        plt.title('onaxis')
        plt.savefig('dPhidz_onaxis_test_'+str(k1)+'.png', dpi = 600)
        if showplots: plt.show()
        plt.close(fig)
        
        
        fig = plt.figure()
        for k2 in range(NH):
            FSPA_alpha_map = interp_FSPA_short[Horders[k2]](Intens_onaxis_envel_XUV[k1]/units.INTENSITYau)
            FSPA_phase_map = np.multiply(Intens_onaxis_envel_XUV[k1]/units.INTENSITYau, FSPA_alpha_map)
            dFSPA_phase_map = np.multiply(dI_dz_map_XUV[k1]/units.INTENSITYau, FSPA_alpha_map)
            plt.plot(zgrid[0], dFSPA_phase_map[k_t,:], label='H'+str(Horders[k2]))
            
        plt.xlabel('z [m]')
        plt.ylabel('dPhiFSPA/dz [rad/m]')
        plt.legend(loc='best')
        plt.title('onaxis')
        plt.savefig('dPhiFSPA_dz_onaxis_XUV_test_'+str(k1)+'.png', dpi = 600)
        if showplots: plt.show()
        plt.close(fig)
 
    
        fig = plt.figure()
        for k2 in range(NH):
            FSPA_alpha_map = interp_FSPA_short[Horders[k2]](Intens_onaxis_envel[k1]/units.INTENSITYau)
            FSPA_phase_map = np.multiply(Intens_onaxis_envel[k1]/units.INTENSITYau, FSPA_alpha_map)
            dFSPA_phase_map = np.multiply(dI_dz_map[k1]/units.INTENSITYau, FSPA_alpha_map)
            plt.plot(zgrid[0], dFSPA_phase_map[k_t,:], label='H'+str(Horders[k2]))
            
        plt.xlabel('z [m]')
        plt.ylabel('dPhiFSPA/dz [rad/m]')
        plt.legend(loc='best')
        plt.title('onaxis')
        plt.savefig('dPhiFSPA_dz_onaxis_test_'+str(k1)+'.png', dpi = 600)
        if showplots: plt.show()
        plt.close(fig)
        
        

    fig = plt.figure()
    plt.pcolor(zgrid[0], tgrid[0][kt1:kt2], dPhi_dz_map_XUV[0][kt1:kt2,:],shading='auto')
    # plt.xlabel('z [mm]')
    # plt.ylabel('r [mum]')
    plt.title('dPhi/dz')
    plt.ylim(1e-15*np.asarray(tlim))
    plt.colorbar()
    plt.savefig('dPhi_dz_XUV.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
                
    # Field in the vacuum frame
    # plt.plot(1e15*tgrid, Efield_onaxis_s[:,0], linewidth=0.2, label='z=0')
    # plt.plot(1e15*tgrid, Efield_onaxis_s[:,Nz//2], linewidth=0.2, label='z=0.5zmax')
    
    fig = plt.figure()
    for k1 in range(Nfiles):
        # plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, label=labels[k1])
        # plt.plot(1e15*tgrid[k1], Efield_onaxis_s_XUV[k1][:,-1], linewidth=0.2, label=labels[k1])

        plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2)
        plt.plot(1e15*tgrid[k1], Efield_onaxis_s_XUV[k1][:,-1], linewidth=0.2)
        
    # for k1 in range(Nfiles):
    #     plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, linestyle=linestyles[k1], label=labels[k1])
    # plt.plot(1e15*tgrid[1], Efield_onaxis_s[1][:,-1], linewidth=0.2, linestyle='--', label=labels[1])
    plt.xlim(tlim)
    plt.legend(loc='best')
    plt.xlabel('t [fs]')
    plt.ylabel('E [V/m]')
    plt.savefig(outgraph_name+'.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
    
    # print(zgrid[0][-1])
    # print(zgrid[1][-1])
    # print(1e15*(zgrid[0][-1]-zgrid[1][-1])*(1/units.c_light-1/VF_IR[0]))
    
    dPhi_dz_map_XUV2 = [] 
    
    fig = plt.figure()
    for k1 in range(Nfiles):
        # plt.plot(zgrid[k1], dPhi_dz_map[k1][len(tgrid[k1])//2,:], linewidth=0.2, label=labels[k1])
        # plt.plot(zgrid[k1], dPhi_dz_map_XUV[k1][len(tgrid[k1])//2,:], linewidth=0.2, label=labels[k1])
        
        plt.plot(zgrid[k1], dPhi_dz_map[k1][len(tgrid[k1])//2,:], linewidth=0.2)
        plt.plot(zgrid[k1], dPhi_dz_map_XUV[k1][len(tgrid[k1])//2,:], linewidth=0.2)
        
        k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
        dPhi_dz_map_XUV2.append(dPhi_dz_map[k1] + k0_wave*(nXUV[k1][0]-1))
        # plt.plot(zgrid[k1], dPhi_dz_map_XUV2[k1][len(tgrid[k1])//2,:], linestyle = '--', linewidth=0.2, label=labels[k1])
        plt.plot(zgrid[k1], dPhi_dz_map_XUV2[k1][len(tgrid[k1])//2,:], linestyle = '--', linewidth=0.2)
    
    # for k1 in range(Nfiles):
    #     plt.plot(1e15*tgrid[k1], Efield_onaxis_s[k1][:,-1], linewidth=0.2, linestyle=linestyles[k1], label=labels[k1])
    # plt.plot(1e15*tgrid[1], Efield_onaxis_s[1][:,-1], linewidth=0.2, linestyle='--', label=labels[1])
    plt.legend(loc='best')
    plt.xlabel('z [m]')
    plt.ylabel('Lcoh [m]')
    plt.savefig('Lcoh.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
    

    # Plot all orders for the first file
    fig = plt.figure()
    for k2 in range(NH):
        k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
        Harm_map = dPhi_dz_map[0] + k0_wave*(nXUV[0][k2]-1)
        plt.plot(zgrid[0], Harm_map[len(tgrid[0])//2,:], linestyle = '-', linewidth=0.2, label=str(Horders[k2]))
      
    plt.legend(loc='best')
    plt.xlabel('z [m]')
    plt.ylabel('Lcoh [m]')
    plt.savefig('Lcoh_Harm.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
    
    Cutoff_entry = []
    Plateau_threshold = []
    for k1 in range(Nfiles):
        Cutoff_entry.append(HHG.ComputeCutoff(Intensity_entry[k1]/units.INTENSITYau,
                                              mn.ConvertPhoton(omega0,'omegaSI','omegaau'),
                                              mn.ConvertPhoton(Ip_eV[k1],'eV','omegaau')
                                              )[1])
        Plateau_threshold.append(mn.ConvertPhoton(Ip_eV[k1],'eV','omegaau')/mn.ConvertPhoton(omega0,'omegaSI','omegaau'))
    
    
    # test f1 values
    f1_val = []
    f1_grid = np.linspace(15, 23, 1000)
    for k1 in range(len(f1_grid)):
        f1_val.append(XUV_index.getf1(gas_type+'_'+'NIST', mn.ConvertPhoton(f1_grid[k1]*omega0, 'omegaSI', 'eV')))
    f1_val = np.asarray(f1_val)
    fig = plt.figure()
    plt.plot(f1_grid,f1_val)
    # plt.legend(loc='best')
    plt.xlabel('H [-]')
    plt.ylabel('f1 [SI]')
    plt.savefig('f1_vals_NIST.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
    
    f1_val = []
    f1_grid = np.linspace(15, 23, 1000)
    for k1 in range(len(f1_grid)):
        f1_val.append(XUV_index.getf1(gas_type+'_'+'Henke', mn.ConvertPhoton(f1_grid[k1]*omega0, 'omegaSI', 'eV')))
    f1_val = np.asarray(f1_val)
    fig = plt.figure()
    plt.plot(f1_grid,f1_val)
    # plt.legend(loc='best')
    plt.xlabel('H [-]')
    plt.ylabel('f1 [SI]')
    plt.savefig('f1_vals_Henke.png', dpi = 600)
    if showplots: plt.show()
    plt.close(fig)
    
    
    # radial profile of the phase in the comoving frame
            
            
    Efield_probe_cmplx_envel = []; phase_map_probe = []
    
    
    for k1 in range(Nfiles):
        
        Efocus = np.sqrt(2.0*Intensity_entry[0]/(units.c_light*units.eps0))
        Efield_Gauss = np.moveaxis(mn.GaussianBeamEfieldConstruct(
            rgrid[0],zgrid_probe1[0],tgrid[0],Efocus,w0[0],tFWHM[0],mn.ConvertPhoton(omega0, 'omegaSI', 'lambdaSI')
            ),[0,1],[1,2])
    
        Nt, Nr, Nz = Efield_probe[k1].shape
        Efield_probe_cmplx_envel.append(np.zeros(Efield_probe[k1].shape,dtype=complex))
        rem_fast_oscillations = np.exp(-1j*omega0*tgrid[k1])
        
        Efield_Gauss_cmplx_envel = np.zeros(Efield_Gauss.shape,dtype=complex)
        
        for k2 in range(Nz):
            for k3 in range(Nr):
                Efield_probe_cmplx_envel[k1][:,k3,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_probe[k1][:,k3,k2])

                Efield_Gauss_cmplx_envel[:,k3,k2] = rem_fast_oscillations*mn.complexify_fft(Efield_Gauss[:,k3,k2])


        k_t = mn.FindInterval(tgrid[k1], t_fix) # find time index
        
        phase_map_probe.append(np.angle(Efield_probe_cmplx_envel[k1][k_t,:,:])) # ordering (r,z)
        Intens_loc = mn.FieldToIntensitySI(abs(Efield_probe_cmplx_envel[k1]))
         # mn.FieldToIntensitySI(Efield_onaxis_envel)
        
        phase_map_Gauss = np.angle(Efield_Gauss_cmplx_envel[k_t,:,:]) # ordering (r,z)
        Intens_loc_Gauss = mn.FieldToIntensitySI(abs(Efield_Gauss_cmplx_envel))
        
        

        
        fig = plt.figure()    
        plt.plot(rgrid[k1],Intens_loc[k_t,:,0]) 
        plt.plot(rgrid[k1],Intens_loc_Gauss[k_t,:,0],linestyle='--') 
        plt.plot(rgrid[k1],Intens_loc[k_t,:,1])  
        plt.plot(rgrid[k1],Intens_loc_Gauss[k_t,:,1],linestyle='--') 
        plt.plot(rgrid[k1],Intens_loc[k_t,:,2])   
        plt.plot(rgrid[k1],Intens_loc_Gauss[k_t,:,2],linestyle='--') 
        if showplots: plt.show()
        plt.close(fig)
    
    
    fig = plt.figure()  
    plt.title('z-phase')
    plt.plot(zgrid_probe1[0],np.unwrap(phase_map_probe[0][0,:]))    
    
    plt.plot(zgrid_probe1[0],np.unwrap(phase_map_Gauss[0,:]))    
    
    if showplots: plt.show()
    plt.close(fig)
    
    fig = plt.figure()   
    plt.title('r-IR-phase')
    plt.plot(rgrid[0],np.unwrap(phase_map_probe[0][:,1]))    
    
    plt.plot(rgrid[0],np.unwrap(phase_map_Gauss[:,1]))    
    
    if showplots: plt.show()
    plt.close(fig)
    
    # FSPA phase
    FSPA_alpha_map2 = interp_FSPA_short[15](Intens_loc/units.INTENSITYau)
    FSPA_phase_map2 = np.multiply(Intens_loc/units.INTENSITYau, FSPA_alpha_map2)

    fig = plt.figure()    
    plt.title('r-FSPA-phase')
    dum = FSPA_phase_map2[k_t,:,1]
    plt.plot(1e6*rgrid[0],dum-dum[0])
    dum = -(FSPA_alpha_map2[k_t,0,1]/units.INTENSITYau) * Intens_loc[k_t,:,1]
    plt.plot(1e6*rgrid[0],dum-dum[0])
    
    dum = np.unwrap(phase_map_Gauss[:,1])
    plt.plot(1e6*rgrid[0],dum-dum[0])
    
    if showplots: plt.show()
    # plt.close(fig)

os.chdir(cwd)
# print('done')
