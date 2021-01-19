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

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")

filename = "results_1.h5"

file_path = os.path.join(results_path,filename)

with h5py.File(file_path, 'r') as InputArchive:

    print(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'))
    omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_wavelength','N'),'lambdaSI','omegaSI')
    print(omega0)
    tgrid = InputArchive['/outputs/tgrid'][:]; Nt = len(tgrid)
    rgrid = InputArchive['/outputs/rgrid'][:]; Nr = len(rgrid)
    zgrid = InputArchive['/outputs/zgrid'][:]; Nz = len(zgrid)
    electron_density_map = InputArchive['/outputs/output_plasma'][:]
    Efield = InputArchive['/outputs/output_field'][:]

    w0 = 1e-2*mn.readscalardataset(InputArchive,'/inputs/laser_beamwaist','N')
    zR = np.pi * w0**2 / mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')
    zgridzR = np.linspace(0,zR,100)

    plasma_frequency_map = np.sqrt((units.elcharge ** 2) / (units.eps0 * units.elmass)) * np.sqrt(electron_density_map)


    E_cmplx_onax = np.zeros((Nz,Nt),dtype=complex)
    for k1 in range(Nz):
        E_cmplx_onax[k1,:] = mn.complexify_fft(Efield[:,0,k1])

    E_cmplx_zfirst = np.zeros((Nr, Nt), dtype=complex)
    E_cmplx_zlast = np.zeros((Nr, Nt), dtype=complex)

    for k1 in range(Nr):
        E_cmplx_zfirst[k1,:] = mn.complexify_fft(Efield[:,k1,0])
        E_cmplx_zlast[k1,:] = mn.complexify_fft(Efield[:, k1, 30])


    plt.pcolor(tgrid,rgrid,E_cmplx_zfirst.real)
    plt.colorbar()
    plt.savefig('map_zfix1cmplx.png', dpi = 600)
    plt.show()

    plt.pcolor(tgrid,rgrid,E_cmplx_zlast.real)
    plt.colorbar()
    plt.savefig('map_zfix2cmplx.png', dpi = 600)
    plt.show()

    Rz = zgrid[30] + zR**2/zgrid[30]
    Curv_coeff = np.pi / (Rz*mn.ConvertPhoton(omega0,'omegaSI','lambdaSI'))

    E_cmplx_zlast_envel = np.exp(-1j*omega0*tgrid)*E_cmplx_zlast
    E_cmplx_zlast_angle = np.angle(E_cmplx_zlast_envel)

    plt.plot(rgrid, np.unwrap(E_cmplx_zlast_angle[:,512]),linewidth=0.2)
    plt.plot(rgrid, -Curv_coeff*rgrid**2, linewidth=0.2)
    plt.savefig('Phase_zfix2.png', dpi = 600)
    plt.show()




    E_cmplx_onax_envel = np.exp(-1j*omega0*tgrid)*E_cmplx_onax

    # E_cmplx_onax_angle = np.arctan2(E_cmplx_onax_envel.imag,E_cmplx_onax_envel.real)
    E_cmplx_onax_angle = np.angle(E_cmplx_onax_envel)

    plt.plot(tgrid, Efield[:,0,0],linewidth=0.2)
    plt.plot(tgrid, E_cmplx_onax[0,:].real,linewidth=0.2)
    plt.plot(tgrid, E_cmplx_onax_envel[0,:].real, linewidth=0.2)
    plt.savefig('Efieldfirstenvel.png', dpi = 600)
    plt.show()

    plt.plot(zgrid, E_cmplx_onax_angle[:,512],linewidth=0.2)
    plt.plot(zgrid, np.arctan(zgrid/zR), linewidth=0.2)
    plt.plot(zgridzR, np.arctan(zgridzR / zR), linewidth=0.2)
    plt.savefig('Phase_thalf.png', dpi = 600)
    plt.show()

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

    plt.plot(tgrid, Efield[:,0,1],linewidth=0.2)
    plt.plot(tgrid, Efield[:, 0, 2],linewidth=0.2)
    plt.plot(tgrid, Efield[:, 0, 30], linewidth=0.2)
    plt.savefig('Efield.png', dpi = 600)
    plt.show()

    plt.pcolor(tgrid,rgrid,Efield[:,:,0].T)
    plt.colorbar()
    plt.savefig('Efield_zfix.png', dpi = 600)
    plt.show()



    plt.pcolor(tgrid,rgrid,plasma_frequency_map[:,:,0].T)
    plt.colorbar()
    plt.savefig('map_zfix.png', dpi = 600)
    plt.show()


    plt.pcolor(zgrid,tgrid,np.squeeze(plasma_frequency_map[:,0,:]))
    plt.colorbar()
    plt.savefig('map_onax.png', dpi = 600)
    plt.show()

    plt.pcolor(zgrid,rgrid,np.squeeze(plasma_frequency_map[512,:,:]))
    plt.colorbar()
    plt.savefig('map_thalf.png', dpi = 600)
    plt.show()

print('test1')


## refractive index contribs


#### THE MAIN PROGRAM #####


# class NumericalParamsClass:
#   def __init__(self):
#     pass
#
# NumericalParams.z_medium = ParamFile['inputs/'+'jetpositions'][()]

# test matplot

plt.plot([1, 2, 3, 4])
plt.ylabel('some numbers')
plt.show()

