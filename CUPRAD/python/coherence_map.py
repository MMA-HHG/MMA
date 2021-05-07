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
import matplotlib.pyplot as plt
import re
# import copy
import glob
import pandas as pd
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index
import subprocess
import dataformat_CUPRAD as dfC

import gc

try:
    univ_input_path = os.environ['UNIV_INPUT_PATH']
except:
    univ_input_path = os.path.join("D:\git", "universal_input")
        

print('path:', univ_input_path)
print('test1')
# sys.exit(0)
# There may be a large memory overload in this implementation, it's for testing putposes

# =============================================================================
# Inputs of the script

arguments = sys.argv

showplots = not('-nodisplay' in arguments)

if ('-here' in arguments):
    results_path = os.getcwd()
else:
    # results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
    # results_path = os.path.join("D:\data", "Discharges")
    results_path = os.path.join("D:\data", "Discharges", "f_scan_Kerr")
    results_path = os.path.join("D:\data", "Discharges","TDSE", "t2")
    # results_path = os.path.join("D:\data", "Discharges", "preion2")

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'

# results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "tests", "sim")

cwd = os.getcwd()
os.chdir(results_path)
files = glob.glob('results_*.h5')
os.chdir(cwd)

vacuum_frame = True

# files = ['results_1.h5','results_25.h5','results_2.h5']
# files = ['results_1.h5','results_2.h5', 'results_3.h5']
# files = ['results.h5']
# files = ['results_19.h5']

# files = ['results_10.h5']

# files = ['results_1.h5','results_10.h5', 'results_15.h5']

# files = ['results_10.h5']

out_h5name = 'analyses.h5'

# create inputs from the free-form inp
# run_args = ['python3', univ_input_path +'/create_universal_HDF5.py',
#             '-i', 'map3.inp', '-ohdf5', 'inputs_tmp.h5',
#             '-g', 'inputs']

# subprocess.run(run_args) # problematic in Spyder
# # runfile(run_args)

print('archive created')
print('my dir is:',os.getcwd())




try:
    with h5py.File('inputs_tmp.h5', 'r') as Parameters:
        gas_type = mn.readscalardataset(Parameters, 'inputs/gas_type', 'S')
        XUV_table_type = mn.readscalardataset(Parameters, 'inputs/XUV_table_type', 'S') # 'NIST' # {Henke, NIST}
        Horders = Parameters['inputs/Horders'][:].tolist() # mn.readscalardataset(Parameters, 'inputs/gas_type', 'S') # [15, 17, 19, 21, 23] # [19, 21, 23, 25, 27]
        Lcoh_saturation = mn.readscalardataset(Parameters, 'inputs/Lcoh_saturation', 'N') #0.02
        Lcoh_zero = mn.readscalardataset(Parameters, 'inputs/Lcoh_zero', 'N') # 0.0
        H_shift_mask = mn.readscalardataset(Parameters, 'inputs/H_shift_mask', 'N') # 4.0
        tlim = Parameters['inputs/tlim'][:].tolist() # [-60.0,60.0]
        t_probe = Parameters['inputs/t_probe'][:]
        rmax = mn.readscalardataset(Parameters, 'inputs/rmax', 'N')
        dr = mn.readscalardataset(Parameters, 'inputs/dr', 'N')
        
        print_parameters = Parameters['inputs/print_parameters'][:].tolist()
        print_parameters = [k1.decode() for k1 in print_parameters]
        
except:
    print('error reading hdf5 file, using defaults') 
    
    gas_type = 'Kr'
    XUV_table_type = 'NIST' # {Henke, NIST}
    Horders = [15, 17, 19, 21, 23] # [19, 21, 23, 25, 27]
    Lcoh_saturation = 0.02
    Lcoh_zero = 0.0
    H_shift_mask = 4.0
    tlim = [-60.0,60.0]
    t_probe = 1e-15*np.asanyarray([-5.0, 0.0, 5.0])
    rmax = 130e-6 # only for analyses
    dr = rmax/40.0

fix_saturation = True

full_resolution = False

invoke_garbage_collector = True

OutPath = 'outputs'

# t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase
# q = 23 # harmonic of our interest
# Coherence_length = True
# Beam_analysis = False
# Efield_analysis = False
# Gaussian_curvature = True # print Gaussian curvature, it is applied only in the first run
# fluence_source = 'computed' # options: 'file', 'computed'
# sys.exit(0)

# Get FSPA 
with h5py.File(os.path.join(cwd,'FSPA_tables_Krypton_test.h5'),'r') as h5_FSPA_tables:
    interp_FSPA_short = HHG.FSPA.get_dphase(h5_FSPA_tables,'Igrid','Hgrid','short/dphi')

# =============================================================================
# The body of the script
if os.path.exists(OutPath) and os.path.isdir(OutPath):
  shutil.rmtree(OutPath)
  print('deleted previous results')
os.mkdir(OutPath)

os.chdir(OutPath)

prop_cycle = plt.rcParams['axes.prop_cycle']
colors_plt = prop_cycle.by_key()['color']
linestyles_plt = ['--','-','-.'] 


NH = len(Horders)
Nt_probe = len(t_probe)

# df_delta_t_columns = pd.MultiIndex.from_arrays([['pressure', 'ionisation', 'Intensity'        ,'IR phase','IR group','XUV phase'],
#                                                 ['{[bar]}',     '[\%]',      '$10^{14}$ [W/cm2]','[fs]',    '[fs]',    '[fs]']])


#  ['pressure \\\ [mbar]', 'ionisation \\\ [\%]', 'Intensity \\\ $10^{14}$ [W/cm2]', 'IR phase \\\ [fs]', 'IR group \\\ [fs]', 'XUV phase \\\ [fs]']
# df_delta_t = pd.DataFrame(columns = df_delta_t_columns)
# df_delta_t.loc[len(df_delta_t)+1] = [10, 15, 20]

with h5py.File(out_h5name,'w') as OutFile: # this file contains numerical analyses
# This is the archive, where we store numerical results at the instant. It may
# be further extended by storing also figure sources etc.
    for fname in files:
    # Here we loop over all result files in the destiantion folder.
        filename = fname
        
        # Using regexp to obtain the rusult number as a literal.
        k_sim = re.findall(r'\d+', filename)
        k_sim = k_sim[0] 
        
        # Prepare group in OutFile and open the repsective result file
        # grp = OutFile.create_group('results_' + str(k_sim))   
        file_path = os.path.join(results_path,filename)
        print('processing:', file_path)             
        with h5py.File(file_path, 'r') as InputArchive:
            
            # load data
            res = dfC.get_data(InputArchive, r_resolution = [full_resolution, dr, rmax])
            
            ### Shift to the vacuum frame
            if vacuum_frame:
                res.vacuum_shift()
         

            # ===============================================
            # Complexify the fields: E(r,z,t) = Re(E_cmplx(r,z,t))
            E_trz_cmplx_envel = np.zeros(res.E_trz.shape,dtype=complex)
            rem_fast_oscillations = np.exp(-1j*res.omega0*res.tgrid)
            
            for k1 in range(res.Nz):
                for k2 in range(res.Nr):
                    E_trz_cmplx_envel[:,k2,k1] = rem_fast_oscillations*mn.complexify_fft(res.E_trz[:,k2,k1])


            # ===============================================
            # Retrieve the phase at the time of our interest
            t_probe_ind = mn.FindInterval(res.tgrid, t_probe) # find time indices
            
            phase = np.angle(E_trz_cmplx_envel[t_probe_ind,:,:]) # ordering (r,z)
            Intens = mn.FieldToIntensitySI(abs(E_trz_cmplx_envel))
            grad_z_I = np.gradient(Intens[t_probe_ind,:,:],res.zgrid,axis=2,edge_order=2)
            
            grad_z_phase = np.zeros(phase.shape)
            for k1 in range(res.Nr):
                for k2 in range(Nt_probe):
                    phase_rfix = np.unwrap(phase[k2,k1,:])
                    grad_z_phase[k2,k1,:] = np.gradient(phase_rfix,res.zgrid,edge_order=2)
                    
            nXUV = [];
            for k1 in range(NH):
                q = Horders[k1]
                f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*res.omega0, 'omegaSI', 'eV'))
                nXUV.append(1.0 - res.rho0_init*units.r_electron_classical*(mn.ConvertPhoton(q*res.omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))
            
            FSPA_alphas = []; grad_z_phase_FSPA = []
            for k1 in range(NH):
                q = Horders[k1]
                FSPA_alphas.append(-interp_FSPA_short[q](Intens[t_probe_ind,:,:]/units.INTENSITYau))
                grad_z_phase_FSPA.append(FSPA_alphas[k1]*grad_z_I /units.INTENSITYau)
                
            # Cut-off map
            Cutoff = HHG.ComputeCutoff(Intens/units.INTENSITYau,
                                       mn.ConvertPhoton(res.omega0,'omegaSI','omegaau'),
                                       mn.ConvertPhoton(res.Ip_eV,'eV','omegaau')
                                       )[1]
              
            # loop outputs over times and harmonic orders
            fig1, ax1 = plt.subplots()
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            title_string = dfC.create_param_string(print_parameters, res)
            # title_string = ', ' + res.preionisation_string + ', ' + res.pressure_string
            for k1 in range(Nt_probe):
              t_string = ', '+"{:.1f}".format(1e15*res.tgrid[t_probe_ind[k1]])+' fs'
            
              Cutoff_loc = Cutoff[t_probe_ind[k1],:,:]

              fig7, ax7 = plt.subplots()
              map1 = ax7.pcolor(1e3*res.zgrid, 1e6*res.rgrid, Cutoff[t_probe_ind[k1],:,:], shading='auto')
              map2 = ax7.contour(1e3*res.zgrid, 1e6*res.rgrid, Cutoff[t_probe_ind[k1],:,:], Horders, colors = "black")
              ax7.set_xlabel('z [mm]'); ax7.set_ylabel('r [mum]'); ax7.set_title('Cutoff'+title_string+t_string )
              fig7.colorbar(map1) 
              fig7.savefig('Cutoff_t'+str(k1)+'_sim'+str(k_sim)+'.png', dpi = 600)
              # if showplots: plt.show(fig7)
              # plt.close(fig7)
              
              for k2 in range(NH):
                q = Horders[k2]
                
                # Create mask
                H_mask = 1.0*Cutoff_loc[:] # instead of copy.deepcopy(Cutoff_loc[:])
                H_mask[H_mask < q-H_shift_mask] = np.nan # https://stackoverflow.com/questions/38800532/set-color-for-nan-values-in-matplotlib/38800580
                H_mask = np.ma.masked_where(np.isnan(H_mask),H_mask) 
                
                # onax
                if (linestyles_plt[k1] == '-'):              
                  ax1.plot(1e3*res.zgrid,q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)),label='H'+str(q),
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])               
                  ax2.plot(1e3*res.zgrid,grad_z_phase_FSPA[k2][k1,0,:],label='H'+str(q),
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])                
                  ax3.plot(1e3*res.zgrid,
                         q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)) + grad_z_phase_FSPA[k2][k1,0,:],
                         label='H'+str(q), color=colors_plt[k2], linestyle=linestyles_plt[k1])
                else:
                  ax1.plot(1e3*res.zgrid,q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)),
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])               
                  ax2.plot(1e3*res.zgrid,grad_z_phase_FSPA[k2][k1,0,:],
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])                
                  ax3.plot(1e3*res.zgrid,
                         q*(grad_z_phase[k1,0,:] + res.k0_wave*(nXUV[k2]-1)) + grad_z_phase_FSPA[k2][k1,0,:],
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])
                  
                  
                fig4, ax4 = plt.subplots()                
                masked = np.ma.masked_where(np.isnan(H_mask), q*(grad_z_phase[k1,:,:] + res.k0_wave*(nXUV[k2]-1))+grad_z_phase_FSPA[k2][k1,:,:] )                
                map1 = ax4.pcolor(1e3*res.zgrid, 1e6*res.rgrid,
                                  masked, # q*(grad_z_phase[k1,:,:] + k0_wave*(nXUV[k2]-1))+grad_z_phase_FSPA[k2][k1,:,:],
                                     shading='auto')
                ax4.set_xlabel('z [mm]'); ax4.set_ylabel('r [mum]'); ax4.set_title('dPhi/dz, full, H'+str(q)+title_string+t_string ) 
                fig4.colorbar(map1)
                fig4.savefig('dPhi_dz_t'+str(k1)+'_H'+str(q)+'_sim'+str(k_sim)+'.png', dpi = 600)
                
                fig8, ax8 = plt.subplots()
                
                masked = np.ma.masked_where(np.isnan(H_mask), abs(q*(grad_z_phase[k1,:,:] + res.k0_wave*(nXUV[k2]-1))+grad_z_phase_FSPA[k2][k1,:,:]))
                map1 = ax8.pcolor(1e3*res.zgrid, 1e6*res.rgrid,
                                  masked,
                                     shading='auto',vmin=0.0)
                ax8.set_xlabel('z [mm]'); ax8.set_ylabel('r [mum]'); ax8.set_title('|dPhi/dz|, full, H'+str(q)+title_string+t_string ) 
                fig8.colorbar(map1)
                fig8.savefig('abs_dPhi_dz_t'+str(k1)+'_H'+str(q)+'_sim'+str(k_sim)+'.png', dpi = 600)
                
                fig5, ax5 = plt.subplots()
                dum = abs( 1.0/ (q*(grad_z_phase[k1,:,:] + res.k0_wave*(nXUV[k2]-1))+grad_z_phase_FSPA[k2][k1,:,:] ) )
                masked = np.ma.masked_where(np.isnan(H_mask), abs( 1.0/ (q*(grad_z_phase[k1,:,:] + res.k0_wave*(nXUV[k2]-1))+grad_z_phase_FSPA[k2][k1,:,:] ) ))
                if ((np.max(dum) > Lcoh_saturation) or fix_saturation):
                    map1 = ax5.pcolor(1e3*res.zgrid, 1e6*res.rgrid, masked, shading='auto', vmax=Lcoh_saturation)
                else:
                    map1 = ax5.pcolor(1e3*res.zgrid, 1e6*res.rgrid, masked, shading='auto')
                
                ax5.set_xlabel('z [mm]'); ax5.set_ylabel('r [mum]'); ax5.set_title('Lcoh, H'+str(q)+title_string+t_string )         
                fig5.colorbar(map1)
                fig5.savefig('Lcoh_t'+str(k1)+'_H'+str(q)+'_sim'+str(k_sim)+'.png', dpi = 600)
                
             
            ax1.set_xlabel('z [mm]'); ax1.set_ylabel('dPhi/dz [1/m]'); ax1.set_title('beam phase'+title_string)   
            ax2.set_xlabel('z [mm]'); ax2.set_ylabel('dPhi/dz [1/m]'); ax2.set_title('FSPA (atom)'+title_string)
            ax3.set_xlabel('z [mm]'); ax3.set_ylabel('dPhi/dz [1/m]'); ax3.set_title('full phase'+title_string) 
            ax1.legend(loc='best')
            ax2.legend(loc='best')
            ax3.legend(loc='best')
            
            fig1.savefig('phase_onax_beam_sim'+str(k_sim)+'.png', dpi = 600)
            fig2.savefig('phase_onax_FSPA_sim'+str(k_sim)+'.png', dpi = 600)
            fig3.savefig('phase_onax_full_sim'+str(k_sim)+'.png', dpi = 600)
            
            
            # Intnesity shape
            fig9, ax9 = plt.subplots()
            ax9.plot(1e15*res.tgrid,Cutoff[:,0,0], color=colors_plt[0],label='z='+"{:.1f}".format(1e3*res.zgrid[0])) 
            ax9.plot(1e15*res.tgrid,Cutoff[:,0,(res.Nz-1)//2], color=colors_plt[1],label='z='+"{:.1f}".format(1e3*res.zgrid[(res.Nz-1)//2])) 
            ax9.plot(1e15*res.tgrid,Cutoff[:,0,-1], color=colors_plt[2],label='z='+"{:.1f}".format(1e3*res.zgrid[-1]))
            ax9.set_xlim(tlim)
            ax9.legend(loc='best')
            ax9.set_xlabel('t [fs]'); ax9.set_ylabel('I [cutoff]'); ax9.set_title('onaxis intensity'+title_string)   
            fig9.savefig('Intens_onax_sim'+str(k_sim)+'.png', dpi = 600)
            
            k_t = mn.FindInterval(res.tgrid, 0.0)
            fig10, ax10 = plt.subplots()
            ax10.plot(1e6*res.rgrid,Cutoff[k_t,:,0], color=colors_plt[0],label='z='+"{:.1f}".format(1e3*res.zgrid[0])) 
            ax10.plot(1e6*res.rgrid,Cutoff[k_t,:,(res.Nz-1)//2], color=colors_plt[1],label='z='+"{:.1f}".format(1e3*res.zgrid[(res.Nz-1)//2])) 
            ax10.plot(1e6*res.rgrid,Cutoff[k_t,:,-1], color=colors_plt[2],label='z='+"{:.1f}".format(1e3*res.zgrid[-1]))
            ax10.legend(loc='best')
            ax10.set_xlabel('r [mum]'); ax10.set_ylabel('I [cutoff]'); ax10.set_title('t=0 fs, intensity'+title_string)   
            fig10.savefig('Intens_tfix_sim'+str(k_sim)+'.png', dpi = 600)           
            
            if showplots: plt.show()
            plt.close()
            
            # Curvature map + plasma map ## IMPLEMENT GAUSSIAN SEPARATELY
            res.get_plasma(InputArchive, r_resolution = [full_resolution, dr, rmax])
            beam_curvature = np.zeros(phase.shape)
            ddr_beam_curvature = np.zeros(phase.shape)
            for k1 in range(Nt_probe):
                for k2 in range(res.Nz):                
                    phase_tz_fix = np.unwrap(phase[k1,:,k2])
                    beam_curvature[k1,:,k2] = phase_tz_fix - phase_tz_fix[0]
                    ddr_beam_curvature[k1,:,k2] = np.gradient(phase_tz_fix,res.rgrid,edge_order=2)


            for k1 in range(Nt_probe):
                t_string = ', '+"{:.1f}".format(1e15*res.tgrid[t_probe_ind[k1]])+' fs'
                
                # curvature
                fig11, ax11 = plt.subplots()
                
                map11 = ax11.pcolor(1e3*res.zgrid,
                                  1e6*res.rgrid, beam_curvature[k1,:,:],
                                  shading='auto', cmap='plasma')
                
                ax11.set_xlabel('z [mm]'); ax11.set_ylabel('r [mum]')
                ax11.set_title('Curvature [rad]'+title_string+t_string ) 
                fig11.colorbar(map11)
                
                fig11.savefig('Curvature_t'+str(k1)+'_sim'+str(k_sim)+'.png', dpi = 600)
                if showplots: plt.show()
                plt.close()
                
                # curvature gradient
                fig13, ax13 = plt.subplots()
                
                map13 = ax13.pcolor(1e3*res.zgrid,
                                  1e6*res.rgrid, ddr_beam_curvature[k1,:,:],
                                  shading='auto', cmap='plasma')
                
                ax13.set_xlabel('z [mm]'); ax11.set_ylabel('r [mum]')
                ax13.set_title('d/dr Curv. [rad/m]'+title_string+t_string ) 
                fig13.colorbar(map13)
                
                fig13.savefig('ddr_Curvature_t'+str(k1)+'_sim'+str(k_sim)+'.png', dpi = 600)
                if showplots: plt.show()
                plt.close()
                
                # plasma
                fig12, ax12 = plt.subplots()
                
                map12 = ax12.pcolor(1e3*res.plasma.zgrid, 1e6*res.plasma.rgrid,
                                  100*res.plasma.value_trz[t_probe_ind[k1],:,:]/res.rho0_init,
                                  shading='auto', cmap='plasma')
                
                ax12.set_xlabel('z [mm]'); ax12.set_ylabel('r [mum]')
                ax12.set_title('el. density [%]'+ title_string + t_string) 
                fig12.colorbar(map12)
                
                fig12.savefig('Plasma_t'+str(k1)+'_sim'+str(k_sim)+'.png', dpi = 600)
                if showplots: plt.show()
                plt.close()

            
            
            fig12, ax12 = plt.subplots()
            
            map12 = ax12.pcolor(1e3*res.plasma.zgrid, 1e6*res.plasma.rgrid,
                              100*res.plasma.value_trz[-1,:,:]/res.rho0_init,
                              shading='auto', cmap='plasma')
            
            ax12.set_xlabel('z [mm]'); ax12.set_ylabel('r [mum]')
            ax12.set_title('el. density [%], end'+ title_string) 
            fig12.colorbar(map12)
            
            fig12.savefig('Plasma_end_sim'+str(k_sim)+'.png', dpi = 600)
            if showplots: plt.show()
            plt.close()
            
            
                
        
            # treat fluence with full precision
            res2 = dfC.get_data(InputArchive)
            res2.get_Fluence(InputArchive, fluence_source='computed')
            
            k_r = mn.FindInterval(res2.rgrid, rmax)
            
            radius_inv_e2 = dfC.measure_beam(
                res2.Fluence.rgrid, res2.Fluence.value, mn.measure_beam_max_ratio_zeromax, np.exp(-2.0) ) 
            radius_RMS = dfC.measure_beam(
                res2.Fluence.rgrid, res2.Fluence.value, mn.measure_beam_RMS )         
            
            fig1, ax1 = plt.subplots()
            
            map1 = ax1.pcolor(1e3*res2.Fluence.zgrid,
                              1e6*res2.Fluence.rgrid[:k_r], res2.Fluence.value[:k_r,:],
                              shading='auto', cmap='plasma')
            ax1.plot(1e3*res2.Fluence.zgrid, 1e6*radius_RMS, '--', linewidth=1, color = 'k')
            ax1.plot(1e3*res2.Fluence.zgrid, 1e6*radius_inv_e2, '-', linewidth=1, color = 'k')
            
            ax1.set_ylim([0,1e6*rmax])
            ax1.set_xlabel('z [mm]'); ax1.set_ylabel('r [mum]'); ax1.set_title('Fluence ['+res2.Fluence.units+']'+title_string)
            fig1.colorbar(map1)
            fig1.savefig('Fluence_sim'+str(k_sim)+'.png', dpi = 600)
            if showplots: plt.show()
            plt.close()
            
        if invoke_garbage_collector:
            del res
            del res2
            gc.collect()
            plt.close('all')
        

os.chdir(cwd)
print('done')