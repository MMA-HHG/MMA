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

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")
cwd = os.getcwd()
os.chdir(results_path)
files = glob.glob('results_*.h5')
os.chdir(cwd)

files = ['results_1.h5']

out_h5name = 'analyses.h5'

q = 23 # harmonic of our interest
t_fix = 0.0 # the time of our interest to inspect e.g. phase



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

            # =================================================================
            # Load data

        
            print(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'))
            omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
            
            print(omega0)
            tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(tgrid)
            rgrid = InputArchive['/outputs/rgrid'][:]; Nr = len(rgrid)
            zgrid = InputArchive['/outputs/zgrid'][:]; Nz = len(zgrid)
            electron_density_map = InputArchive['/outputs/output_plasma'][:]
            Efield = InputArchive['/outputs/output_field'][:]
            
            inverse_GV = InputArchive['/logs/inverse_group_velocity_SI'][()]
            VG_IR = 1.0/inverse_GV
            w0 = 1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_beamwaist','N')
            zR = np.pi * w0**2 / mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
            zgridzR = np.linspace(0,zR,100)
        
            rho0_atm = 1e6 * mn.readscalardataset(InputArchive, '/inputs/calculated/medium_effective_density_of_neutral_molecules','N')
            pressure = InputArchive['/inputs/medium_pressure_in_bar'][()]
            
            
            # =================================================================
            # Process the data
            # plasma_frequency_map = np.sqrt((units.elcharge ** 2) / (units.eps0 * units.elmass)) * np.sqrt(electron_density_map)
    
            # ===============================================
            # Refractive indexes of the IR and XUV
            # XUV is retrieved from the tables https://henke.lbl.gov/optical_constants/asf.html
            # See also (2.80) in Attwood: X-Rays and Extreme Ultraviolet Radiation Principles and Applications; 2016
            # IR is given by the Dalgarno, Kingston; 1960
            
            f1 = XUV_index.getf('Ar', mn.ConvertPhoton(q*omega0, 'omegaSI', 'eV'))[0]
            nXUV = 1 - pressure*rho0_atm*units.r_electron_classical*(mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')**2)*f1/(2.0*np.pi)
            VF_XUV = units.c_light/nXUV # phase velocity of XUV
            
            nIR = IR_index.getsusc('Ar', mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))
            nIR = np.sqrt(1.0 + pressure*nIR)
            VF_IR = units.c_light/nIR # phase velocity of IR
            
            # ===============================================
            # Compute delays and store them

            dset_id = grp.create_dataset('delta_t_IR_phase', data=-1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VF_IR))
            dset_id.attrs['units'] = np.string_('[fs]')            
 
            dset_id = grp.create_dataset('delta_t_IR_group', data=-1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VG_IR))
            dset_id.attrs['units'] = np.string_('[fs]')

            dset_id = grp.create_dataset('delta_t_XUV_phase', data=-1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VF_XUV))
            dset_id.attrs['units'] = np.string_('[fs]')

            dset_id = grp.create_dataset('T0_IR', data=1e15*mn.ConvertPhoton(omega0,'omegaSI','T0SI'))
            dset_id.attrs['units'] = np.string_('[fs]')            

            dset_id = grp.create_dataset('T0_XUV', data=1e15*mn.ConvertPhoton(q*omega0,'omegaSI','T0SI'))
            dset_id.attrs['units'] = np.string_('[fs]')

            print('characteristic velocities')
            print(units.c_light,VF_IR,1.0/inverse_GV,VF_XUV, sep='\n')
            print('delta_t_IR_phase =', -1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VF_IR),'fs')      
            print('delta_t_IR_group =', -1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VG_IR),'fs')            
            print('delta_t_XUV_phase =', -1e15*zgrid[-1]*(1.0/units.c_light - 1.0/VF_XUV),'fs')          
            print('T0_IR =', 1e15*mn.ConvertPhoton(omega0,'omegaSI','T0SI'),'fs')           
            print('T0_XUV =', 1e15*mn.ConvertPhoton(q*omega0,'omegaSI','T0SI'),'fs')


            # ===============================================
            # Complexify the fields: E(r,z,t) = Re(E_cmplx(r,z,t))
            
            ## Complexify the field using fft -> ifft & remove fast oscillations
            Efield_cmplx_envel = np.zeros((Nt,Nr,Nz),dtype=complex)
            rem_fast_oscillations = np.exp(-1j*omega0*tgrid)
            for k1 in range(Nz):
                for k2 in range(Nr):
                    Efield_cmplx_envel[:,k2,k1] = rem_fast_oscillations*mn.complexify_fft(Efield[:,k2,k1])


            # ===============================================
            # Retrieve the phase at the time of our interest
            
            k_t = mn.FindInterval(tgrid, t_fix) # find time index
            phase_map = np.angle(Efield_cmplx_envel[k_t,:,:]) # ordering (r,z)
        
            Lcoh_map = np.zeros((Nr//2,Nz-1))
            Lcoh_map_dd = np.zeros((Nr//2,Nz-1))
            for k1 in range(Nr//2):
                phase_r = np.unwrap(phase_map[k1,:])
                Lcoh_map[k1, 0] = (phase_r[1] - phase_r[0]) / (zgrid[1] - zgrid[0])
                Lcoh_map[k1, -1] = (phase_r[-1] - phase_r[-2]) / (zgrid[-1] - zgrid[-2])
                Lcoh_map_dd[k1, 0] = Lcoh_map[k1, 0]
                Lcoh_map_dd[k1, -1] = Lcoh_map[k1, -1]
                for k2 in range(1, Nz - 1):
                    Lcoh_map[k1, k2] = (phase_r[k2+1] - phase_r[k2]) / (zgrid[k2+1] - zgrid[k2])
                    Lcoh_map_dd[k1, k2] = mn.ddx_arb(k2,zgrid,phase_r)
            # Lcoh_map = np.abs(np.pi / Lcoh_map)
        
                # Lcoh_map[k1, 0] = (phase_r[1]-phase_r[0])/(zgrid[1]-zgrid[0])
                # Lcoh_map[k1, -1] = (phase_r[-1]-phase_r[-2])/(zgrid[-1]-zgrid[-2])
                # for k2 in range(1,Nz-1):
                #     Lcoh_map[k1,k2] = mn.ddx_arb(k2,zgrid,phase_r)
        
            Lcoh_map2 = np.abs(np.pi/Lcoh_map)
        
            # compute also local curvature
            Curvature_Gaussian_map = np.zeros((Nr//2,Nz))
            Curvature_map = np.zeros((Nr // 2, Nz))
            for k1 in range(Nz):
                if (k1 == 0): Curv_coeff = 0
                else:
                    Rz = zgrid[k1] + zR ** 2 / zgrid[k1]
                    Curv_coeff = np.pi / (Rz * mn.ConvertPhoton(omega0, 'omegaSI', 'lambdaSI'))
                Curvature_Gaussian_map[:(Nr//2), k1] = -Curv_coeff*(rgrid[:(Nr//2)])**2
                Curvature_map[:(Nr//2), k1] = np.unwrap(phase_map[:(Nr//2), k1])
                Curvature_map[:(Nr // 2), k1] = Curvature_map[:(Nr//2), k1] - Curvature_map[0, k1] # start always at 0
        
            Curvature_map = Curvature_map - np.max(Curvature_map)
        
            ## Fluence estimation
            Fluence = np.zeros((Nr//2, Nz))
            for k1 in range(Nz):
                for k2 in range(Nr//2):
                    Fluence[k2, k1] = sum(abs(Efield[:, k2, k1])**2)
        
            ## thalf_ionisation map
            ionisation_tfix_map = np.zeros((Nr//2, Nz))
            for k1 in range(Nz):
                for k2 in range(Nr//2):
                    ionisation_tfix_map[k2, k1] = electron_density_map[k_t,k2,k1]
        
            # plt.pcolor(zgrid,rgrid,Lcoh_map,vmax=0.1)
            plt.pcolor(zgrid[:-1], rgrid[:(Nr//2)], Lcoh_map)
            plt.colorbar()
            plt.savefig('dPhidz_map_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
            
            plt.pcolor(zgrid[:-1], rgrid[:(Nr//2)], Lcoh_map_dd)
            plt.colorbar()
            plt.savefig('dPhidz_map_dd_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
            
            plt.pcolor(zgrid[:-1], rgrid[:(Nr//2)], Lcoh_map-Lcoh_map_dd)
            plt.colorbar()
            plt.savefig('dPhidz_map_error_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            plt.pcolor(1e3*zgrid[:-1], 1e6*rgrid[:(Nr//2)], Lcoh_map2, vmax=0.3)
            plt.xlabel('z [mm]')
            plt.ylabel('r [mum]')
            plt.title('Lcoh [m]')
            plt.colorbar()
            plt.savefig('Lcoh_map_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            plt.pcolor(zgrid,rgrid,phase_map)
            plt.colorbar()
            plt.savefig('Phase_z_unwrp_map_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            plt.pcolor(1e3*zgrid, 1e6*rgrid[:(Nr//2)], Curvature_map, vmax = 0)
            plt.xlabel('z [mm]')
            plt.ylabel('r [mum]')
            plt.title('phi_Curv [rad]')
            plt.colorbar()
            plt.savefig('Curvature_map_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            # plt.plot(Curvature_map[:,0])
            # plt.savefig('0curv.png', dpi=600)
            # plt.show()
        
            plt.pcolor(1e3*zgrid, 1e6*rgrid[:(Nr//2)], Curvature_Gaussian_map, vmax = 0)
            plt.xlabel('z [mm]')
            plt.ylabel('r [mum]')
            plt.title('phi_Curv [rad]')
            plt.colorbar()
            plt.savefig('Curvature_map_Gauss_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            # Fluence
            plt.pcolor(1e3*zgrid, 1e6*rgrid[:(Nr//2)], Fluence)
            plt.xlabel('z [mm]')
            plt.ylabel('r [mum]')
            plt.title('Fluence [arb.u.]')
            plt.colorbar()
            plt.savefig('Fluence_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            # thalf ionisation
            plt.pcolor(1e3*zgrid, 1e6*rgrid[:(Nr//2)], 100.0*ionisation_tfix_map/rho0_atm)
            plt.xlabel('z [mm]')
            plt.ylabel('r [mum]')
            plt.title('electron density [%]')
            plt.colorbar()
            plt.savefig('ionisation_thalf_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            # thalf ionisation
            plt.plot(1e15*tgrid, 100.0*electron_density_map[:,0,15]/rho0_atm)
            plt.xlabel('t [fs]')
            plt.ylabel('electron density [%]')
            plt.title('z = ' + str(1e3*zgrid[15]) + ' mm')
            plt.savefig('ionisation_middle_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
            
            
            ## Get the field shift in the vacuum frame
            
            Efield_onaxis_s = np.squeeze(Efield[:,0,:])
            for k1 in range(Nz):
                ogrid_nn, FE_s, NF = mn.fft_t_nonorm(tgrid, Efield_onaxis_s[:,k1])
                
                delta_z = zgrid[k1]
                delta_t = inverse_GV*delta_z
                delta_t = delta_t - delta_z/units.c_light
                
                FE_s = np.exp(1j*ogrid_nn*delta_t) * FE_s
                
                tnew, E_s = mn.ifft_t_nonorm(ogrid_nn,FE_s,NF)
                Efield_onaxis_s[:,k1] = E_s.real
                
             # shifted
            plt.plot(1e15*tgrid, Efield_onaxis_s[:,0], linewidth=0.2)
            plt.plot(1e15*tgrid, Efield_onaxis_s[:,Nz//2], linewidth=0.2)
            plt.plot(1e15*tgrid, Efield_onaxis_s[:,-1], linewidth=0.2)
            plt.xlabel('t [fs]')
            plt.ylabel('E [V/m]')
            # plt.title('z = ' + str(1e3*zgrid[15]) + ' mm')
            plt.savefig('Field_shift_lines_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
            
            
            plt.pcolor(1e3*zgrid, 1e15*tgrid, Efield_onaxis_s)
            plt.xlabel('z [mm]')
            plt.ylabel('t [fs]')
            plt.title('shifted field [V/m]')
            plt.colorbar()
            plt.savefig('Field_shift_'+str(k_sim)+'.png', dpi = 600)
            plt.show()
        
            # plt.plot(zgrid, phase_map[0,:],linewidth=0.2)
            # plt.plot(zgrid, phase_map[1, :], linewidth=0.2)
            # plt.plot(zgrid, phase_map[2, :], linewidth=0.2)
            # plt.plot(zgrid, phase_map[3, :], linewidth=0.2)
            # plt.plot(zgrid, phase_map[4, :], linewidth=0.2)
            # plt.savefig('Phases.png', dpi = 600)
            # plt.show()
        
        
        
            ## Get plasma contribution
        
            ## Plasma and geometry cannot be separated with the numerical solution


















    # plasma_frequency_map[512, :, :]

    # onax-thalf phase

    # n_plasma_map = np.sqrt(1-(plasma_frequency_map/omega0)**2)
    #
    # k_plas_map = 25*(2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))*(n_plasma_map-1)
    #
    # plt.pcolor(tgrid,rgrid,k_plas_map[:,:,0].T)
    # plt.colorbar()
    # plt.savefig('kmap.png', dpi = 600)
    # plt.show()
    #
    # plt.pcolor(tgrid,rgrid,abs(np.pi/k_plas_map[:,:,15].T))
    # plt.colorbar()
    # plt.savefig('Lcohmap.png', dpi = 600)
    # plt.show()
    #
    # xxx = electron_density_map[0:3,0,0]
    # yyy = electron_density_map[0, 0:3, 0]


    ## Get geometrical phase contribution



#     E_cmplx_onax = np.zeros((Nz,Nt),dtype=complex)
#     for k1 in range(Nz):
#         E_cmplx_onax[k1,:] = mn.complexify_fft(Efield[:,0,k1])
#
#     E_cmplx_zfirst = np.zeros((Nr, Nt), dtype=complex)
#     E_cmplx_zlast = np.zeros((Nr, Nt), dtype=complex)
#
#     for k1 in range(Nr):
#         E_cmplx_zfirst[k1,:] = mn.complexify_fft(Efield[:,k1,0])
#         E_cmplx_zlast[k1,:] = mn.complexify_fft(Efield[:, k1, 30])
#
#
#     plt.pcolor(tgrid,rgrid,E_cmplx_zfirst.real)
#     plt.colorbar()
#     plt.savefig('map_zfix1cmplx.png', dpi = 600)
#     plt.show()
#
#     plt.pcolor(tgrid,rgrid,E_cmplx_zlast.real)
#     plt.colorbar()
#     plt.savefig('map_zfix2cmplx.png', dpi = 600)
#     plt.show()
#
#     Rz = zgrid[30] + zR**2/zgrid[30]
#     Curv_coeff = np.pi / (Rz*mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))
#
#     E_cmplx_zlast_envel = np.exp(-1j*omega0*tgrid)*E_cmplx_zlast
#     E_cmplx_zlast_angle = np.angle(E_cmplx_zlast_envel)
#
#     plt.plot(rgrid, np.unwrap(E_cmplx_zlast_angle[:,512]),linewidth=0.2)
#     plt.plot(rgrid, -Curv_coeff*rgrid**2, linewidth=0.2)
#     plt.savefig('Phase_zfix2.png', dpi = 600)
#     plt.show()
#
#
#
#
#     E_cmplx_onax_envel = np.exp(-1j*omega0*tgrid)*E_cmplx_onax
#
#     # E_cmplx_onax_angle = np.arctan2(E_cmplx_onax_envel.imag,E_cmplx_onax_envel.real)
#     E_cmplx_onax_angle = np.angle(E_cmplx_onax_envel)
#
#     plt.plot(tgrid, Efield[:,0,0],linewidth=0.2)
#     plt.plot(tgrid, E_cmplx_onax[0,:].real,linewidth=0.2)
#     plt.plot(tgrid, E_cmplx_onax_envel[0,:].real, linewidth=0.2)
#     plt.savefig('Efieldfirstenvel.png', dpi = 600)
#     plt.show()
#
#     plt.plot(zgrid, E_cmplx_onax_angle[:,512],linewidth=0.2)
#     plt.plot(zgrid, np.arctan(zgrid/zR), linewidth=0.2)
#     plt.plot(zgridzR, np.arctan(zgridzR / zR), linewidth=0.2)
#     plt.savefig('Phase_thalf.png', dpi = 600)
#     plt.show()
#
#
#
#     plt.plot(tgrid, Efield[:,0,1],linewidth=0.2)
#     plt.plot(tgrid, Efield[:, 0, 2],linewidth=0.2)
#     plt.plot(tgrid, Efield[:, 0, 30], linewidth=0.2)
#     plt.savefig('Efield.png', dpi = 600)
#     plt.show()
#
#     plt.pcolor(tgrid,rgrid,Efield[:,:,0].T)
#     plt.colorbar()
#     plt.savefig('Efield_zfix.png', dpi = 600)
#     plt.show()
#
#
#
#     plt.pcolor(tgrid,rgrid,plasma_frequency_map[:,:,0].T)
#     plt.colorbar()
#     plt.savefig('map_zfix.png', dpi = 600)
#     plt.show()
#
#
#     plt.pcolor(zgrid,tgrid,np.squeeze(plasma_frequency_map[:,0,:]))
#     plt.colorbar()
#     plt.savefig('map_onax.png', dpi = 600)
#     plt.show()
#
#     plt.pcolor(zgrid,rgrid,np.squeeze(plasma_frequency_map[512,:,:]))
#     plt.colorbar()
#     plt.savefig('map_thalf.png', dpi = 600)
#     plt.show()
#
#
#
# ## refractive index contribs
#
#
#
