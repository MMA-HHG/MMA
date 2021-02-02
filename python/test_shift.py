import numpy as np
from scipy import integrate
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

# from skimage.restoration import unwrap_phase

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")

filename = "results_31.h5"

file_path = os.path.join(results_path,filename)


def myfft(x,fx):
    if (len(x)%2==1): x = x[0:-1]
    N = len(x)
    dx = x[1]-x[0]
    Fx = (dx/np.sqrt(np.pi))*np.conj(np.fft.fft(fx)[0:((N//2)+1)])
    xi = (np.pi/dx)*np.linspace(0,1,(N//2)+1)
    return xi, Fx

def imyfft(xi,Fx):
    dxi = xi[1]-xi[0]
    Fx = np.append(Fx, np.flip(np.conj(Fx[1:])))
    fx = ((len(Fx)*dxi)/(2.0*np.sqrt(np.pi)))*np.flip(np.fft.ifft(Fx))
    x = np.linspace(0,2.0*np.pi/dxi,len(fx))
    return x, fx


def complexify_fft(fx):
    N = len(fx)
    fx = np.fft.fft(fx)
    for k1 in range((N//2)+1,N):
        fx[k1] = 0.0
    fx = 2.0*np.fft.ifft(fx)
    return fx

with h5py.File(file_path, 'r') as InputArchive:

    print(1e-2*InputArchive['/inputs/laser_wavelength'][()])
    omega0 = mn.ConvertPhoton(1e-2*InputArchive['/inputs/laser_wavelength'][()],'lambdaSI','omegaSI')
    # omega0 = 2379924242424242.5
    print(omega0)
    tgrid = InputArchive['/outputs/tgrid'][:]
    rgrid = InputArchive['/outputs/rgrid'][:]
    zgrid = InputArchive['/outputs/zgrid'][:]
    electron_density_map = InputArchive['/outputs/output_plasma'][:]
    Efield = InputArchive['/outputs/output_field'][:]
    
    inverse_GV = InputArchive['/logs/inverse_group_velocity_SI'][()]
    # inverse_GV  = 1/units.c_light



    Etest = Efield[:,0,10]

    omegauppe = 10.1 * omega0

    ogrid, FEtest = myfft(tgrid,Etest)
    
    FEtest_norm = np.trapz(np.abs(FEtest)**2,ogrid)
    Etest_norm = np.trapz(np.abs(Etest)**2,tgrid)
    

    plt.plot(ogrid,np.abs(FEtest), linewidth=0.2)
    plt.savefig('FEfield.png', dpi = 600)
    plt.show()

    # plt.plot(np.abs(FEtest[75:100]), linewidth=0.2)
    plt.plot(np.abs(FEtest[125:145]), linewidth=0.2)
    plt.savefig('FEfield_pts.png', dpi = 600)
    plt.show()
    
    plt.plot(ogrid[125:145],np.abs(FEtest[125:145]), linewidth=0.2)
    plt.savefig('FEfield_zoom.png', dpi = 600)
    plt.show()


    plt.plot(tgrid, Etest, linewidth=0.2)
    plt.xlim([-30.e-15,30.e-15])
    plt.savefig('Efield.png', dpi = 600)
    plt.show()

    # test 1
    # delta_t = 50.0e-15  # 0.01
    wavenumber = 2.0*np.pi/mn.ConvertPhoton(omega0,'omegaSI','lambdaSI')

    delta_z = 0.001
    delta_n = 0.0004

    delta_phi = np.pi# wavenumber*delta_n*delta_z
    
    delta_t = delta_phi / omega0

    
    # test 2 - move in the laboratory frame adjusted by 'c' (a vacuum case)
    delta_z = 0.00001
    
    delta_t = delta_z/units.c_light
    
    # test 3 - move in the laboratory frame adjusted by 'vg' (a vacuum case)
    delta_z = 0.00001
    
    delta_t = inverse_GV*delta_z
    
    # test 4 - two shifts
    delta_z = 0.01
    delta_t = inverse_GV*delta_z
    delta_t = delta_t - delta_z/units.c_light
    
    # testing shift using nonorm
    ogrid_nn, FEtest_nn, NF = mn.fft_t_nonorm(tgrid,Etest)
    FEtest_nns = np.exp(1j*ogrid_nn*delta_t) * FEtest_nn
    tnew, FFEtest_nns = mn.ifft_t_nonorm(ogrid_nn,FEtest_nns,NF)
    
    plt.plot(tgrid, Etest, linewidth=0.2)
    plt.plot(tgrid, FFEtest_nns.real, linewidth=0.2)
    plt.title('FFE + shift, 2')
    plt.savefig('FFEfield_s2.png', dpi=600)
    plt.show()    
    
    
    # FEtest_shift = np.exp(1j*np.linspace(0,len(ogrid)-1,len(ogrid))*delta_t)*FEtest
    FEtest_shift = np.exp(1j*ogrid*delta_t)*FEtest
    

    tt, FFEtest = imyfft(ogrid,FEtest)
    FFEtest_shift = imyfft(ogrid, FEtest_shift)[1]

    tt_center = tt - (tt[-1]-tt[0])/2.0

    plt.plot(FFEtest, linewidth=0.2)
    plt.title('FFE')
    plt.savefig('FFEfield.png', dpi = 600)
    plt.show()

    plt.plot(tt_center,FFEtest.real, linewidth=0.2)
    plt.plot(tt_center,FFEtest_shift.real, linewidth=0.2)
    plt.title('FFE + shift')
    plt.savefig('FFEfield_s.png', dpi=600)
    plt.show()

    print('test1')

    # show envelopes
    # complexify fields
    rem_fast_oscillations = np.exp(-1j * omega0 * tt_center)
    Ecmplxe = rem_fast_oscillations*mn.complexify_fft(FFEtest.real)
    Ecmplxes = rem_fast_oscillations*mn.complexify_fft(FFEtest_shift.real)

    Ecmplxea = np.abs(Ecmplxe)
    Ecmplxesa = np.abs(Ecmplxes)

    # the amplitude seems to be wrong
    plt.plot(tt_center,Ecmplxea, linewidth=0.2)
    plt.plot(tt_center,Ecmplxesa, linewidth=0.2)
    plt.title('FFE + shift + envel')
    plt.savefig('FFEfield_cmplxe.png', dpi=600)
    plt.show()



    # Ecmplx = complexify_fft(Etest)
    #
    # Eenvelope = np.exp(-1j*omega0*tgrid)*Ecmplx
    #
    # plt.plot(tgrid, Ecmplx.real, linewidth=0.2)
    # plt.savefig('EcmplxR.png', dpi = 600)
    # plt.show()
    #
    # plt.plot(tgrid, Eenvelope.real, linewidth=0.2)
    # plt.savefig('EenvelR.png', dpi = 600)
    # plt.show()
    #
    # plt.plot(tgrid, np.angle(Eenvelope), linewidth=0.2)
    # plt.savefig('EenvelAngle.png', dpi = 600)
    # plt.show()
    #
    # plt.plot(tgrid, np.unwrap(np.angle(Eenvelope)), linewidth=0.2)
    # plt.savefig('EenvelAngleUnwrp.png', dpi = 600)
    # plt.show()
    #
    # plt.plot(tgrid, abs(Ecmplx.real-Etest), linewidth=0.2)
    # plt.savefig('Error.png', dpi = 600)
    # plt.show()


print('done')
