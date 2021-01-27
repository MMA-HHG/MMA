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

import XUV_refractive_index as XUV_index
import IR_refractive_index as IR_index

f1h = XUV_index.index_funct['Ar']['f1'](29.5)

f1, f2 = XUV_index.getf('Ar',29.5)

susc = IR_index.getsusc('Ar',600e-9)

# Energy = XUV_index.Energy

# from skimage.restoration import unwrap_phase

# results_path = os.path.join("/mnt", "d", "data", "Discharges") # 'D:\data\Discharges'
results_path = os.path.join("D:\data", "Discharges")



# filename = "results_1.h5"
#
# file_path = os.path.join(results_path,filename)
#
#
# def myfft(x,fx):
#     if (len(x)%2==1): x = x[0:-1]
#     N = len(x)
#     dx = x[1]-x[0]
#     Fx = (dx/np.sqrt(2.0*np.pi))*np.conj(np.fft.fft(fx)[0:((N//2)+1)])
#     xi = (np.pi/dx)*np.linspace(0,1,(N//2)+1)
#     return xi, Fx
#
# def imyfft(xi,Fx):
#     dxi = xi[1]-xi[0]
#     Fx = np.append(Fx, np.flip(np.conj(Fx[1:])))
#     fx = ((len(Fx)*dxi)/(np.sqrt(2.0*np.pi)))*np.flip(np.fft.ifft(Fx))
#     x = np.linspace(0,2.0*np.pi/dxi,len(fx))
#     return x, fx
#
#
# with h5py.File(file_path, 'r') as InputArchive:
#
#     print(1e-2*InputArchive['/inputs/laser_wavelength'][()])
#     omega0 = mn.ConvertPhoton(1e-2*InputArchive['/inputs/laser_wavelength'][()],'lambdaSI','omegaSI')
#     # omega0 = 2379924242424242.5
#     print(omega0)
#     tgrid = InputArchive['/outputs/tgrid'][:]
#     rgrid = InputArchive['/outputs/rgrid'][:]
#     zgrid = InputArchive['/outputs/zgrid'][:]
#     electron_density_map = InputArchive['/outputs/output_plasma'][:]
#     Efield = InputArchive['/outputs/output_field'][:]
#
#     Etest = Efield[:,0,10]
#
#     omegauppe = 10.1 * omega0
#     # Etest = np.exp(-(tgrid/29.72626301008067e-15)**2)*np.cos(omega0*tgrid+1.0) #Efield[:, 0, 30]
#     # Etest = np.exp(-(tgrid / 29.72626301008067e-15) ** 2 + 1j*(omegauppe - omega0)*tgrid)
#     # Etest = np.exp(-1j*omegauppe*tgrid)*Etest
#     # Etest = Etest.real
#
#     # Etest = Etest.astype(np.float32)
#
#     ogrid, FEtest = myfft(tgrid,Etest)
#
#     plt.plot(ogrid,np.abs(FEtest), linewidth=0.2)
#     plt.savefig('FEfield.png', dpi = 600)
#     plt.show()
#
#     plt.plot(np.abs(FEtest[75:100]), linewidth=0.2)
#     plt.savefig('FEfield_pts.png', dpi = 600)
#     plt.show()
#
#
#     plt.plot(tgrid, Etest, linewidth=0.2)
#     plt.xlim([-30.e-15,30.e-15])
#     plt.savefig('Efield.png', dpi = 600)
#     plt.show()
#
#     tt, FFEtest = imyfft(ogrid,FEtest)
#
#     plt.plot(FFEtest, linewidth=0.2)
#     plt.savefig('FFEfield.png', dpi = 600)
#     plt.show()
#
#     print('test1')
#
#     Ecmplx = complexify_fft(Etest)
#
#     Eenvelope = np.exp(-1j*omega0*tgrid)*Ecmplx
#
#     plt.plot(tgrid, Ecmplx.real, linewidth=0.2)
#     plt.savefig('EcmplxR.png', dpi = 600)
#     plt.show()
#
#     plt.plot(tgrid, Eenvelope.real, linewidth=0.2)
#     plt.savefig('EenvelR.png', dpi = 600)
#     plt.show()
#
#     plt.plot(tgrid, np.angle(Eenvelope), linewidth=0.2)
#     plt.savefig('EenvelAngle.png', dpi = 600)
#     plt.show()
#
#     plt.plot(tgrid, np.unwrap(np.angle(Eenvelope)), linewidth=0.2)
#     plt.savefig('EenvelAngleUnwrp.png', dpi = 600)
#     plt.show()
#
#     plt.plot(tgrid, abs(Ecmplx.real-Etest), linewidth=0.2)
#     plt.savefig('Error.png', dpi = 600)
#     plt.show()
#
#
# print('test2')
#
#
# ## refractive index contribs
#
#
# #### THE MAIN PROGRAM #####
#
#
# # class NumericalParamsClass:
# #   def __init__(self):
# #     pass
# #
# # NumericalParams.z_medium = ParamFile['inputs/'+'jetpositions'][()]
#
# # test matplot
#
# plt.plot([1, 2, 3, 4])
# plt.ylabel('some numbers')
# plt.show()

