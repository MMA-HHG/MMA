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
import pandas as pd
import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index

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

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'

# results_path = os.path.join("D:\TEMP", "OCCIGEN_CUPRAD", "tests", "sim")

cwd = os.getcwd()
os.chdir(results_path)
files = glob.glob('results_*.h5')
os.chdir(cwd)

vacuum_frame = True

# files = ['results_1.h5','results_25.h5','results_2.h5']
files = ['results_1.h5','results_2.h5', 'results_3.h5']
# files = ['results.h5']
# files = ['results_19.h5']

out_h5name = 'analyses.h5'

# q = 23 # harmonic of our interest
t_fix = 0.0e-15 # the time of our interest to inspect e.g. phase
fluence_source = 'computed' # options: 'file', 'computed'


gas_type = 'Kr'
XUV_table_type = 'NIST' # {Henke, NIST}

Horders = [15, 17, 19, 21, 23]# [19, 21, 23, 25, 27]

Lcoh_saturation = 0.05
Lcoh_zero = 0.0

tlim = [-60.0,60.0]


full_resolution = False
Coherence_length = True
Beam_analysis = False
Efield_analysis = False

rmax = 200e-6 # only for analyses
dr = rmax/40.0

t_probe = 1e-15*np.asanyarray([-5.0, 0.0, 5.0])

Gaussian_curvature = True # print Gaussian curvature, it is applied only in the first run

OutPath = 'outputs'


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

df_delta_t_columns = pd.MultiIndex.from_arrays([['pressure', 'ionisation', 'Intensity'        ,'IR phase','IR group','XUV phase'],
                                                ['{[bar]}',     '[\%]',      '$10^{14}$ [W/cm2]','[fs]',    '[fs]',    '[fs]']])


#  ['pressure \\\ [mbar]', 'ionisation \\\ [\%]', 'Intensity \\\ $10^{14}$ [W/cm2]', 'IR phase \\\ [fs]', 'IR group \\\ [fs]', 'XUV phase \\\ [fs]']
df_delta_t = pd.DataFrame(columns = df_delta_t_columns)
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
        grp = OutFile.create_group('results_' + str(k_sim))   
        file_path = os.path.join(results_path,filename)
        print('processing:', file_path)             
        with h5py.File(file_path, 'r') as InputArchive:
            omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
            k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
            tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(tgrid)
            rgrid = InputArchive['/outputs/rgrid'][:]; Nr = len(rgrid)            
            zgrid = InputArchive['/outputs/zgrid'][:]; Nz = len(zgrid)            
            if full_resolution:
                kr_step = 1; Nr_max = Nr
            else:
                dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
                rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid)                
            E_trz = InputArchive['/outputs/output_field'][:,0:Nr_max:kr_step,:Nz] # Arrays may be over-allocated by CUPRAD
            inverse_GV = InputArchive['/logs/inverse_group_velocity_SI'][()]
            VG_IR = 1.0/inverse_GV               
            rho0_init = 1e6 * mn.readscalardataset(InputArchive, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
            
            ### Shift to the vacuum frame
            if vacuum_frame:
                E_vac = np.zeros(E_trz.shape)   
                for k1 in range(Nz):
                    delta_z = zgrid[k1] # local shift
                    delta_t_lab = inverse_GV*delta_z # shift to the laboratory frame
                    delta_t_vac = delta_t_lab - delta_z/units.c_light # shift to the coordinates moving by c.
                    for k2 in range(Nr):
                        ogrid_nn, FE_s, NF = mn.fft_t_nonorm(tgrid, E_trz[:,k2,k1]) # transform to omega space        
                        FE_s = np.exp(1j*ogrid_nn*delta_t_vac) * FE_s # phase factor        
                        tnew, E_s = mn.ifft_t_nonorm(ogrid_nn,FE_s,NF)
                        E_vac[:,k2,k1] = E_s.real
                    
                E_trz = E_vac          

            # ===============================================
            # Complexify the fields: E(r,z,t) = Re(E_cmplx(r,z,t))
            E_trz_cmplx_envel = np.zeros(E_trz.shape,dtype=complex)
            rem_fast_oscillations = np.exp(-1j*omega0*tgrid)
            
            for k1 in range(Nz):
                for k2 in range(Nr):
                    E_trz_cmplx_envel[:,k2,k1] = rem_fast_oscillations*mn.complexify_fft(E_trz[:,k2,k1])


            # ===============================================
            # Retrieve the phase at the time of our interest
            t_probe_ind = mn.FindInterval(tgrid, t_probe) # find time indices
            
            phase = np.angle(E_trz_cmplx_envel[t_probe_ind,:,:]) # ordering (r,z)
            Intens = mn.FieldToIntensitySI(abs(E_trz_cmplx_envel))
            grad_z_I = np.gradient(Intens[t_probe_ind,:,:],zgrid,axis=2,edge_order=2)
            
            grad_z_phase = np.zeros(phase.shape)
            for k1 in range(Nr):
                for k2 in range(Nt_probe):
                    phase_rfix = np.unwrap(phase[k2,k1,:])
                    grad_z_phase[k2,k1,:] = np.gradient(phase_rfix,zgrid,edge_order=2)
                    
            nXUV = [];
            for k1 in range(NH):
                q = Horders[k1]
                f1 = XUV_index.getf1(gas_type+'_'+XUV_table_type, mn.ConvertPhoton(q*omega0, 'omegaSI', 'eV'))
                nXUV.append(1.0 - rho0_init*units.r_electron_classical*(mn.ConvertPhoton(q*omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi))
            
            FSPA_alphas = []; grad_z_phase_FSPA = []
            for k1 in range(NH):
                q = Horders[k1]
                FSPA_alphas.append(-interp_FSPA_short[q](Intens[t_probe_ind,:,:]/units.INTENSITYau))
                grad_z_phase_FSPA.append(FSPA_alphas[k1]*grad_z_I /units.INTENSITYau)
                
                
            # loop outputs over times and harmonic orders
            fig1, ax1 = plt.subplots()
            fig2, ax2 = plt.subplots()
            fig3, ax3 = plt.subplots()
            for k1 in range(Nt_probe):
              for k2 in range(NH):
                q = Horders[k2]
                
                # onax
                ax1.plot(1e3*zgrid,q*(grad_z_phase[k1,0,:] + k0_wave*(nXUV[k2]-1)),label='H'+str(q),
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])
                ax1.set_xlabel('z [mm]'); ax1.set_ylabel('dPhi/dz [1/m]'); ax1.set_title('beam phase')      
                
                ax2.plot(1e3*zgrid,grad_z_phase_FSPA[k2][k1,0,:],label='H'+str(q),
                         color=colors_plt[k2], linestyle=linestyles_plt[k1])
                ax2.set_xlabel('z [mm]'); ax2.set_ylabel('dPhi/dz [1/m]'); ax2.set_title('FSPA (atom)') 
                
                ax3.plot(1e3*zgrid,
                         q*(grad_z_phase[k1,0,:] + k0_wave*(nXUV[k2]-1)) + grad_z_phase_FSPA[k2][k1,0,:],
                         label='H'+str(q), color=colors_plt[k2], linestyle=linestyles_plt[k1])
                ax3.set_xlabel('z [mm]'); ax3.set_ylabel('dPhi/dz [1/m]'); ax3.set_title('full phase') 
                
            ax1.legend(loc='best')
            ax2.legend(loc='best')
            ax3.legend(loc='best')
            
            fig1.savefig('phase_onax_beam_'+str(k_sim)+'.png', dpi = 600)
            fig2.savefig('phase_onax_FSPA_'+str(k_sim)+'.png', dpi = 600)
            fig3.savefig('phase_onax_full_'+str(k_sim)+'.png', dpi = 600)
                
                
            
            
            
            
            

os.chdir(cwd)
print('done') 

           
            
    #         if Coherence_length:
    #             # Compute dPhi/dz at t = t_fix
    #             dPhi_dz_map = np.zeros((Nr,Nz))
    #             for k1 in range(Nr):
    #                 phase_r = np.unwrap(phase_map[k1,:])
    #                 dPhi_dz_map[k1, 0] = (phase_r[1] - phase_r[0]) / (zgrid[1] - zgrid[0])
    #                 dPhi_dz_map[k1, -1] = (phase_r[-1] - phase_r[-2]) / (zgrid[-1] - zgrid[-2])
    #                 for k2 in range(1, Nz - 1):
    #                     dPhi_dz_map[k1, k2] = mn.ddx_arb(k2,zgrid,phase_r)
                    

    #         # ===============================================
    #         # Coherence length
            
    #         # The coherence lenght taking only the dephasing from CUPRAD in the group-velocity frame
    #             Lcoh_map = np.abs(np.pi/dPhi_dz_map)
                
    #             k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
    #             Lcoh_map_XUV = np.abs(np.pi/(dPhi_dz_map + k0_wave*0.0))

    #         # ===============================================
    #         # Local curvature
            
    #         if Beam_analysis:
            
    #             # Local curvature: phase evolution for a fixed z compared to the on-axis phase
    #             Curvature_Gaussian_map = np.zeros((Nr,Nz)) # Gaussian beam for a comparison
    #             Curvature_map = np.zeros((Nr, Nz)) # Phas from CUPRAD
    #             for k1 in range(Nz):
    #                 if (k1 == 0):
    #                     Curv_coeff = 0
    #                 else:
    #                     Rz = zgrid[k1] + zR ** 2 / zgrid[k1]
    #                     Curv_coeff = np.pi / (Rz * mn.ConvertPhoton(omega0, 'omegaSI', 'lambdaSI'))
    #                 Curvature_Gaussian_map[:, k1] = -Curv_coeff*(rgrid)**2
    #                 Curvature_map[:, k1] = np.unwrap(phase_map[:, k1])
    #                 Curvature_map[:, k1] = Curvature_map[:, k1] - Curvature_map[0, k1] # start always at 0
            
    #             Curvature_map = Curvature_map - np.max(Curvature_map)
            
  
    #         # ===============================================
    #         # Get the fluence
            
    #         # Fluence is either computed from Efield or taken from the file
    #             if (fluence_source == 'file'):
    #                 Fluence = InputArchive['/longstep/fluence'][:,:]
    #                 zgrid_Fluence = InputArchive['/longstep/zgrid_analyses2'][:]
    #                 Fluence_units = 'SI'                
                    
    #             elif (fluence_source == 'computed'):                
    #                 zgrid_Fluence = zgrid
    #                 Fluence = np.zeros((Nr, Nz))
    #                 for k1 in range(Nz):
    #                     for k2 in range(Nr):
    #                         Fluence[k2, k1] = sum(abs(Efield[:, k2, k1])**2)
    #                 Fluence_units = 'arb.u.'
                
            
    #         # ===============================================
    #         # The ionisation map at t = t_fix
            
    #             ionisation_tfix_map = np.zeros((Nr, Nz))
    #             for k1 in range(Nz):
    #                 for k2 in range(Nr):
    #                     ionisation_tfix_map[k2, k1] = electron_density_map[k_t,k2,k1]



    #         # =================================================================
    #         # Print outputs
            
    #         if Coherence_length:
                
    #             # dPhi/dz
    #             for k1 in range(NH):
    #                 k0_wave = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
                    
    #                 fig = plt.figure()
    #                 plt.pcolor(zgrid, rgrid, dPhi_dz_map + k0_wave*(nXUV[k1]-1.0))
    #                 plt.colorbar()
    #                 plt.savefig('dPhidz_map_'+str(k_sim)+'_'+str(Horders[k1])+'.png', dpi = 600)
    #                 if showplots: plt.show()
    #                 plt.close(fig)
                
    #             # Coherence length
    #                 fig = plt.figure()
    #                 Lcoh_map_XUV = np.abs(np.pi/(dPhi_dz_map + k0_wave*(nXUV[k1]-1.0)))
    #                 plt.pcolor(1e3*zgrid, 1e6*rgrid, Lcoh_map_XUV, vmin=Lcoh_zero, vmax=Lcoh_saturation)
    #                 plt.xlabel('z [mm]')
    #                 plt.ylabel('r [mum]')
    #                 plt.title('Lcoh [m]')
    #                 plt.colorbar()
    #                 plt.savefig('Lcoh_map_'+str(k_sim)+'_'+str(Horders[k1])+'.png', dpi = 600)
    #                 plt.show()
    #                 if showplots: plt.show()
    #                 plt.close(fig)
                    
    #             fig = plt.figure()
    #             plt.plot(zgrid, dPhi_dz_map[0,:])
    #             plt.xlabel('z [m]')
    #             plt.ylabel('dPhi/dz')
    #             plt.title('onaxis')
    #             plt.savefig('dPhidz_onaxis_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
                
    #         if Coherence_length or Beam_analysis:
    #             # Phase(r,z,t=t_fix) # not unwrapped, should be not difficult in a smooth case, or use some 2D-unwprapping
    #             fig = plt.figure()
    #             plt.pcolor(zgrid,rgrid,phase_map)
    #             plt.colorbar()
    #             plt.title('Phi [rad]')
    #             plt.savefig('Phase_z_map_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
        
    #         if Beam_analysis:
    #             # Curvature of the beam 
    #             fig = plt.figure()
    #             plt.pcolor(1e3*zgrid, 1e6*rgrid, Curvature_map, vmax = 0)
    #             plt.xlabel('z [mm]')
    #             plt.ylabel('r [mum]')
    #             plt.title('phi_Curv [rad]')
    #             plt.colorbar()
    #             plt.savefig('Curvature_map_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
            
    #             # reference Gaussian curvature 
    #             if Gaussian_curvature:
    #                 Gaussian_curvature = False
    #                 fig = plt.figure()
    #                 plt.pcolor(1e3*zgrid, 1e6*rgrid, Curvature_Gaussian_map, vmax = 0)
    #                 plt.xlabel('z [mm]')
    #                 plt.ylabel('r [mum]')
    #                 plt.title('phi_Curv [rad]')
    #                 plt.colorbar()
    #                 plt.savefig('Curvature_map_Gauss.png', dpi = 600)
    #                 if showplots: plt.show()
    #                 plt.close(fig)
            
    #             # Fluence
    #             fig = plt.figure()
    #             plt.pcolor(1e3*zgrid_Fluence, 1e6*rgrid, Fluence)
    #             plt.xlabel('z [mm]')
    #             plt.ylabel('r [mum]')
    #             plt.title('Fluence ['+Fluence_units+']')
    #             plt.colorbar()
    #             plt.savefig('Fluence_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
            
    #             # ionisation(r,z,t=t_fix)
    #             fig = plt.figure()
    #             plt.pcolor(1e3*zgrid, 1e6*rgrid, 100.0*ionisation_tfix_map/rho0_init)
    #             plt.xlabel('z [mm]')
    #             plt.ylabel('r [mum]')
    #             plt.title('electron density [%]')
    #             plt.colorbar()
    #             plt.savefig('ionisation_thalf_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
            
    #             # ionisation(r=0,z=zmax/2,t)
    #             fig = plt.figure()
    #             plt.plot(1e15*tgrid, 100.0*electron_density_map[:,0,Nz//2]/rho0_init)
    #             plt.xlabel('t [fs]')
    #             plt.ylabel('electron density [%]')
    #             plt.title('z = ' + str(1e3*zgrid[Nz//2]) + ' mm')
    #             plt.savefig('ionisation_middle_'+str(k_sim)+'.png', dpi = 600)
    #             if showplots: plt.show()
    #             plt.close(fig)
            
    #         # if Efield_analysis:
    #         #     # Field in the vacuum frame
    #         #     plt.plot(1e15*tgrid, Efield_onaxis_s[:,0], linewidth=0.2, label='z=0')
    #         #     plt.plot(1e15*tgrid, Efield_onaxis_s[:,Nz//2], linewidth=0.2, label='z=0.5zmax')
    #         #     plt.plot(1e15*tgrid, Efield_onaxis_s[:,-1], linewidth=0.2, label='z=zmax')
    #         #     plt.legend(loc='best')
    #         #     plt.xlabel('t [fs]')
    #         #     plt.ylabel('E [V/m]')
    #         #     plt.savefig('Field_shift_lines_'+str(k_sim)+'.png', dpi = 600)
    #         #     plt.show()
                
    #         #     # E(r=0,z,t) in the vacuum frame
    #         #     plt.pcolor(1e3*zgrid, 1e15*tgrid, Efield_onaxis_s)
    #         #     plt.xlabel('z [mm]')
    #         #     plt.ylabel('t [fs]')
    #         #     plt.title('shifted field [V/m]')
    #         #     plt.colorbar()
    #         #     plt.savefig('Field_shift_'+str(k_sim)+'.png', dpi = 600)
    #         #     plt.show()


    # df_delta_t_columns_tup = df_delta_t_columns.tolist()
    # df_delta_t = df_delta_t.sort_values(by=df_delta_t_columns_tup[0:3]) 
    
    # format_mapping={df_delta_t_columns_tup[0]: '{:,.0f}',
    #                 df_delta_t_columns_tup[1]: '{:,.0f}',
    #                 df_delta_t_columns_tup[2]: '{:,.2f}',
    #                 df_delta_t_columns_tup[3]: '{:,.3f}',
    #                 df_delta_t_columns_tup[4]: '{:,.3f}',
    #                 df_delta_t_columns_tup[5]: '{:,.3f}'}
    
    # for key, value in format_mapping.items():
    #     df_delta_t[key] = df_delta_t[key].apply(value.format)

    # # https://stackoverflow.com/questions/32744997/python-pandas-apply-formatting-to-each-column-in-dataframe-using-a-dict-mapping
    # print(df_delta_t)
    # # print(df_delta_t.to_latex(float_format="%.3f",escape=False,index=False)) 
    # print(df_delta_t.to_latex(escape=False,index=False)) 
    # df_delta_t.to_latex('ltxtable.out',float_format="%.3f",escape=False,index=False)
    # # with open('lout.out','w') as lout:   