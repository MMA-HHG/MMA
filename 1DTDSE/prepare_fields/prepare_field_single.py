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

arguments = sys.argv
showplots = not('-nodisplay' in arguments)

I0_max = 1e18 # SI
lambdaSI = 800e-9
tFWHMSI = 10e-15

E0_max = np.sqrt(I0_max/units.INTENSITYau)
NE = 100

Nt = 1000;
precision = 'd'

# read parameters
# inputfilename = 'TDSE_create_fields.inp'
# with open(inputfilename, 'r') as InputMP:
#     lines = InputMP.readlines()
    
    
fixed_list = []
varying_list = []
parameters = {}
fixed = {}



# processing_fixed_inputs = False;  processing_varying_inputs = False
# for line in lines:
#     sep_line = line.split()  # separate the line
#     if ((len(sep_line) == 0) or (sep_line[0] == '#') or (sep_line[0] == '##')):
#         pass  # print('empty or commented line')
#     else:
#       if (sep_line[0] == '$fixed'):
#           processing_fixed_inputs = True; processing_varying_inputs = False
#       elif (sep_line[0] == '$varying'):
#           processing_fixed_inputs = False; processing_varying_inputs = True
#       elif processing_fixed_inputs: 
#           fixed_list.append(sep_line[0])
#           if (sep_line[2] == 'R'): fixed[sep_line[0]] = float(sep_line[1])
#           elif (sep_line[2] == 'I'): fixed[sep_line[0]] = int(sep_line[1])
#           elif (sep_line[2] == 'S'): fixed[sep_line[0]] = sep_line[1]
#           else: raise TypeError('line:' + line + '\n Specify datatype by R or I.')
#       elif processing_varying_inputs: 
#           varying_list.append(sep_line[0]) 
#           parameters[sep_line[0]] = [float(sep_line[2]), float(sep_line[3]),
#                                      int(sep_line[4]), sep_line[1]]
#       else: raise NotImplementedError('The input fiel must follow the $fixed - $varying structure now')
#         # n_params = n_params + 1
#         # names = names + sep_line[0] + '\t'
#         # if grouped: groups = groups + sep_line[1] + '\t'
#         # dtypes = dtypes + sep_line[shift+1] + '\t'
#         # units = units + sep_line[shift+2] + '\t'
#         # if firstline:
#         #     content = sep_line[shift+3] + '\t' + sep_line[shift+4] + '\t' + sep_line[shift+5]
#         #     firstline = False
#         # else:
#         #     content = content + '\n' + sep_line[shift+3] + '\t' + sep_line[shift+4] + '\t' + sep_line[shift+5]


# # varying_list = ['phi0','Intensity']
# # parameters = {varying_list[0]: [0,0.5*np.pi, -1, '[rad]'],
# #               varying_list[1]: [0,E0_max,NE, '[a.u.]']}




# # grid of input values
# N = 1; param_grids = []; param_dims = []
# for k1 in range(len(varying_list)): # ensure ordering according to 
#     param = varying_list[k1]
#     N_points = parameters[param][2] + 2
#     if (parameters[param][2] == -1): param_grids.append(
#                                                 np.asarray([parameters[param][1]])
#                                                 )
#     else: param_grids.append( np.linspace( *parameters[param][:2], N_points ) )
#     N *= N_points; param_dims.append(N_points)    
# param_dims = np.asarray(param_dims) 



class pulse_types:
    def __init__(self, pulse_type):
        if (pulse_type == 'sin2'): # sin^2 - envelope
            def sin2pulse(t,omega0,omegac,E0,phi0):
                phi_central = - (np.pi * omega0) / (2.0*omegac) # use the peak of the pulse as the reference for the cosine pulse
                return (t>=0) * (t<=np.pi/omegac) * E0*((np.sin(omegac*t))**2) *  np.cos(omega0*t + phi_central + phi0)
            self.pulse = sin2pulse
            self.inputs_list_direct = ['omega0', 'omegac', 'E0', 'phi0']
            self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
            def inputs_converter(lambdaSI,tFWHMSI,E0,phi0):
                omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
                omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
                return omega0, omegac, E0, phi0
            self.inputs_converter = inputs_converter
            
            def construct_tgrid(Nt, **kwargs):
                omegac = 2.0*np.arccos(2**(-0.25))/(kwargs['tFWHMSI']/units.TIMEau)
                tgrid = np.linspace(0,np.pi/omegac,Nt)
                return tgrid
            self.construct_tgrid = construct_tgrid
            
            def make_field(Nt,**kwargs):
                tgrid = construct_tgrid(Nt, **kwargs)
                Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
                return tgrid, Efield
            self.make_field = make_field

        elif (pulse_type == 'Gaussian'): # sin^2 - envelope
            def Gaussian_pulse(t,omega0,tFWHM,E0,phi0):
                return E0* np.exp(-(2.0*np.log(2.0)*t/tFWHM)**2)*  np.cos(omega0*t + phi0)
            self.pulse = Gaussian_pulse
            self.inputs_list_direct = ['omega0', 'T_FWHM', 'E0', 'phi0']
            self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
            def inputs_converter(lambdaSI,tFWHMSI,E0,phi0):
                omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
                return omega0, tFWHMSI/units.TIMEau, E0, phi0
            self.inputs_converter = inputs_converter
            
            def construct_tgrid(Nt, **kwargs):
                tmax = kwargs['t_expand'] * kwargs['tFWHMSI']/units.TIMEau
                tgrid = np.linspace(-0.5*tmax, 0.5*tmax,Nt)
                return tgrid
            self.construct_tgrid = construct_tgrid
            
            def make_field(Nt,**kwargs):
                tgrid = construct_tgrid(Nt, **kwargs)
                Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
                return tgrid, Efield
            self.make_field = make_field
                
        else: raise NotImplementedError('The input fiel must follow the $fixed - $varying structure now')


# dp = pulse_types(fixed['pulse_type'])

# def inputs_wrapper(k, inputs_list):    
#     MultInd = np.unravel_index(k, param_dims)
#     inputs = []
#     for inp in inputs_list:
#         if (inp in fixed_list): inputs.append(fixed[inp])
#         elif (inp in varying_list):
#             k1 = varying_list.index(inp)
#             inputs.append(param_grids[k1][MultInd[k1]])
#     return inputs

mypulse = pulse_types('sin2')        
    
tgrid = mypulse.construct_tgrid(1000, tFWHMSI=tFWHMSI)
Efield = mypulse.pulse(tgrid,*mypulse.inputs_converter(lambdaSI, tFWHMSI, E0_max, 0.0))


## testplot
if showplots:
    image = pp.figure_driver()    
    image.sf = [pp.plotter() for k1 in range(1)]
    image.sf[0].args = [tgrid, Efield]
    pp.plot_preset(image)   


## store


out_h5name = 'field_input.h5'
omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

with h5py.File(out_h5name,'w') as OutFile:
    mn.adddataset(OutFile, 'IRField/tgrid', tgrid , '[a.u.]' )
    mn.adddataset(OutFile, 'IRField/Field', Efield , '[a.u.]' )


# ## testplot
# omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
# tgrid = np.linspace(0,np.pi/omegac,Nt)

# image = pp.figure_driver()    
# image.sf = [pp.plotter() for k1 in range(1)]
# image.sf[0].args = [tgrid, sin2pulse(tgrid,mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau'),omegac,E0_max,0.0)]
# pp.plot_preset(image)   



# out_h5name = 'fields_table.h5'
# omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')

# with h5py.File(out_h5name,'w') as OutFile:
#     shape = [N,Nt]
#     dset = OutFile.create_dataset('fields_list', shape, precision)
#     dset.attrs['units'] = np.string_('[a.u.]')
    
#     for k1 in range(N):
#         MultInd = np.unravel_index(k1,param_dims)
#         dset[k1,:] = sin2pulse(tgrid,omega0,omegac,param_grids[1][MultInd[1]],param_grids[0][MultInd[0]])

#     OutFile.create_dataset('param_list',data=np.string_(varying_list)) 
    
#     mn.adddataset(OutFile, 'tgrid', tgrid , '[a.u.]' )
    
#     # store grids
#     for k1 in range(len(varying_list)):
#         # OutFile.create_dataset('param_'+str(k1), data=param_grids[k1] )
#         mn.adddataset(OutFile, 'param_'+str(k1), param_grids[k1] , parameters[varying_list[k1]][3] )
#         # .create_dataset('param_'+str(k1), data=param_grids[k1] )

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