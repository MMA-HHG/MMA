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

filename = "results_5.h5"

file_path = os.path.join(results_path,filename)


def myfft(x,fx):
    if (len(x)%2==1): x = x[0:-1]
    N = len(x)
    dx = x[1]-x[0]
    Fx = (dx/np.sqrt(2.0*np.pi))*np.conj(np.fft.fft(fx)[0:((N//2)+1)])
    xi = (np.pi/dx)*np.linspace(0,1,(N//2)+1)
    return xi, Fx

def imyfft(xi,Fx):
    dxi = xi[1]-xi[0]
    Fx = np.append(Fx, np.flip(np.conj(Fx[1:])))
    fx = ((len(Fx)*dxi)/(np.sqrt(2.0*np.pi)))*np.flip(np.fft.ifft(Fx))
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
    print(omega0)
    tgrid = InputArchive['/outputs/tgrid'][:]
    rgrid = InputArchive['/outputs/rgrid'][:]
    zgrid = InputArchive['/outputs/zgrid'][:]
    electron_density_map = InputArchive['/outputs/output_plasma'][:]
    Efield = InputArchive['/outputs/output_field'][:]

    Etest = Efield[:, 0, 30]

    ogrid, FEtest = myfft(tgrid,Etest)

    plt.plot(ogrid,np.abs(FEtest), linewidth=0.2)
    plt.savefig('FEfield.png', dpi = 600)
    plt.show()

    plt.plot(tgrid, Etest, linewidth=0.2)
    plt.savefig('Efield.png', dpi = 600)
    plt.show()

    tt, FFEtest = imyfft(ogrid,FEtest)

    plt.plot(FFEtest, linewidth=0.2)
    plt.savefig('FFEfield.png', dpi = 600)
    plt.show()

    print('test1')

    Ecmplx = complexify_fft(Etest)

    plt.plot(tgrid, Ecmplx.real, linewidth=0.2)
    plt.savefig('EcmplxR.png', dpi = 600)
    plt.show()

    plt.plot(tgrid, abs(Ecmplx.real-Etest), linewidth=0.2)
    plt.savefig('Error.png', dpi = 600)
    plt.show()


print('test2')


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

