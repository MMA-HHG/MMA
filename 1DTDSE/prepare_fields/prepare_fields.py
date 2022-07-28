import numpy as np
import os
import time
import shutil
import h5py
import sys
import units
import mynumerics as mn
import re
# import copy
import glob
import plot_presets as pp

I0_max = 1e18 # SI
lambdaSI = 800e-9
tFWHMSI = 30e-15

E0_max = np.sqrt(I0_max/units.INTENSITYau)
NE = 100

Nt = 1000;
precision = 'd'


param_list = ['phi0','Intensity']
parameters = {param_list[0]: [0,0.5*np.pi, -1, '[rad]'],
              param_list[1]: [0,E0_max,NE, '[a.u.]']}

# grid of input values
N = 1; param_grids = []; param_dims = []
for k1 in range(len(param_list)): # ensure ordering according to 
    param = param_list[k1]
    N_points = parameters[param][2] + 2
    if (parameters[param][2] == -1): param_grids.append(
                                                np.asarray([parameters[param][1]])
                                                )
    else: param_grids.append( np.linspace( *parameters[param][:2], N_points ) )
    N *= N_points; param_dims.append(N_points)    
param_dims = np.asarray(param_dims) 



# sin^2 - envelope
def sin2pulse(t,omega0,omegac,E0,phi0):
    phi_central = - (np.pi * omega0) / (2.0*omegac) # use the peak of the pulse as the reference for the cosine pulse
    return (t>=0) * (t<=np.pi/omegac) * E0*((np.sin(omegac*t))**2) *  np.cos(omega0*t + phi_central + phi0)




## testplot
omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
tgrid = np.linspace(0,np.pi/omegac,Nt)

image = pp.figure_driver()    
image.sf = [pp.plotter() for k1 in range(1)]
image.sf[0].args = [tgrid, sin2pulse(tgrid,mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau'),omegac,E0_max,0.0)]
pp.plot_preset(image)   



out_h5name = 'fields_table.h5'
omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

with h5py.File(out_h5name,'w') as OutFile:
    shape = [N,Nt]
    dset = OutFile.create_dataset('Fields_list', shape, precision)
    dset.attrs['units'] = np.string_('[a.u.]')
    
    for k1 in range(N):
        MultInd = np.unravel_index(k1,param_dims)
        dset[k1,:] = sin2pulse(tgrid,omega0,omegac,param_grids[1][MultInd[1]],param_grids[0][MultInd[0]])

    OutFile.create_dataset('param_list',data=np.string_(param_list))  
    
    # store grids
    for k1 in range(len(param_list)):
        # OutFile.create_dataset('param_'+str(k1), data=param_grids[k1] )
        mn.adddataset(OutFile, 'param_'+str(k1), param_grids[k1] , parameters[param_list[k1]][3] )
        # .create_dataset('param_'+str(k1), data=param_grids[k1] )

print('done')


### Initial point for simulations

# Conditions are:
# No envelop temporally
# Integration on 1 optical cycle
# Trajectory selection (Lewenstein model)
# I0=1e14 W/cm2 // future 7e13
# w0=25 microns // future 10
# Lambda0=800 nm // future 1030 nm
# Spatial mode= TEM00 gaussian
# Pressure=0.1 bar
# Atom: argon// future xenon
# Absorption=1cm response
# Length of medium= 50 micron
# Position of target= @ focus and before z=-100 micron and after z=+100
# Calculation plan: @ farfield this 50 cm from medium output

# + times

# 30 fs - 800 nm
# 200 fs - 1030 nm