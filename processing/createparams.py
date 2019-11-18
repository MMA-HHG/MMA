###################
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE)

# this code defines parameters for 1DTDSE calculation for a vacuum Gaussian beam,

###################

from scipy import special
from scipy import integrate
import numpy as np
import struct
import array
import os
import time
#import ray
#import matlab.engine
#import string
import multiprocessing as mp
import math
#import joblib
import mpi4py
from mpi4py import MPI
#import oct2py
import shutil
import h5py



### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factors to atomic units
TIME = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);

### list of all required parameters - allocate in a separate script?
params = [
'Efield.fieldtype', # 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical
'Eguess', # // Energy of the initial state
'num_r', # // Number of points of the initial spatial grid 16000
'num_exp', # // Number of points of the spatial grid for the expansion
'dx', # // resolution for the grid
'InterpByDTorNT', # // refine resolution only for numerical fields (0 - by dt, 1 - by number of points)
'dt', # // resolution in time 0.0625
'Ntinterp', # // number of intermediate points
'textend', # // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
'analy.writewft', # // writewavefunction (1-writting every tprint)
'analy.tprint', # // time spacing for writing the wavefunction
'x_int', # // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
'PrintGaborAndSpectrum', # // print Gabor and partial spectra (1-yes)
'a_Gabor', # // the parameter of the gabor window [a.u.]
'omegaMaxGabor', # // maximal frequency (source term) in Gabor [a.u.]
'dtGabor', # // spacing in Gabor (source term) [a.u.]
'tmin1window', # // analyse 1st part of the dipole
'tmax1window', # // analyse 1st part of the dipole
'tmin2window', # // analyse 2nd part of the dipole
'tmax2window', # // analyse 2nd part of the dipole
'PrintOutputMethod', # // (0 - only text, 1 - only binaries, 2 - both)
'IonisationFilterForTheSourceTerm', # // filter source term by high-ionisation components (1-yes)
'IonFilterThreshold', # // threshold for the ionisation [-]
#target
'trg.a',
# units
'fieldinau', # // 0-input inatomic units, 1 - in femto and GV/m
'input0', # // 0 field starts with 0, 1 in the middle of the file
# definitions of sin^2 for an analytic field
'sin2.E0', # amplitude
'sin2.omega0',   #  frequency [a.u.]\n"
'sin2.ti',   #  initial time\n"
'sin2.Ncyc',   #  # of cycles\n"
'sin2.phase',  # CEP [in radians, reference is a cosine pulse in A]\n"
# GAUGES
'gauge', #; // 0-length, otherwise velocity, velocity available only for analytic field (A needed)
'transformgauge' #		dumint=fscanf(param,"%i %*[^\n]\n",&transformgauge); // 1 - transform also to another gauge during the calculation, (A needed)
]

print(params)
Nparams = len(params)

# create dictionaries with params
### are these parameters constants (0) or variables (1)? # names must match
ParamsType = dict([ (params[k1] , 0 )  for k1 in range(Nparams)]) # all parametrs are preallocated as constants
ParamsType['sin2.phase'] = 1; ParamsType['sin2.E0'] = 1;
print(ParamsType)

### dictionary containing constant parameters; variable parameters are defined in forcomming loops
ParamsValue = dict([ (params[k1] , None )  for k1 in range(Nparams)]) # all parametrs are preallocated as empty


### the definition of the beam
BeamSpecification={
'Nsim': 4096, # simulations in r #4096
'w0z' : 96.0e-6,
'r_extend' : 4.0,
'E0' : 0.075, 
'z' : 0.05,
'lambda' : 810.0e-9,
'phase0' : 0.0, # initial CEP
'N1cycles' : 20.0
}

### constant parameters list
## an example for an analytic field
ParamsValue['Efield.fieldtype'] = 2; # 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical
ParamsValue['Eguess'] = -1.0 ; # // Energy of the initial state
ParamsValue['num_r'] = 128000; # // Number of points of the initial spatial grid 16000
ParamsValue['num_exp'] = 0; # // Number of points of the spatial grid for the expansion
ParamsValue['dx'] = 0.2;  # // resolution for the grid
ParamsValue['InterpByDTorNT'] = 0; # // refine resolution only for numerical fields (0 - by dt, 1 - by number of points)
ParamsValue['dt'] = 0.0625; # // resolution in time 0.0625
ParamsValue['Ntinterp'] = 1; # // number of intermediate points
ParamsValue['textend'] = 200; # // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
ParamsValue['analy.writewft'] = 0;  # // writewavefunction (1-writting every tprint)
ParamsValue['analy.tprint'] = 10.0; # // time spacing for writing the wavefunction
ParamsValue['x_int'] = 2.0; # // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
ParamsValue['PrintGaborAndSpectrum'] = 0; # // print Gabor and partial spectra (1-yes)
ParamsValue['a_Gabor'] = 8.0; # // the parameter of the gabor window [a.u.]
ParamsValue['omegaMaxGabor'] = 15.0; # // maximal frequency (source term) in Gabor [a.u.]
ParamsValue['dtGabor'] = 10.; #// spacing in Gabor (source term) [a.u.]
ParamsValue['tmin1window'] = 2000.; # // analyse 1st part of the dipole
ParamsValue['tmax1window'] = 5000.; # // analyse 1st part of the dipole
ParamsValue['tmin2window'] = 5250.; # // analyse 2nd part of the dipole
ParamsValue['tmax2window'] = 10000.; #// analyse 2nd part of the dipole
ParamsValue['PrintOutputMethod'] = 2; # // (0 - only text, 1 - only binaries, 2 - both)
ParamsValue['IonisationFilterForTheSourceTerm'] = 1; # // filter source term by high-ionisation components (1-yes)
ParamsValue['IonFilterThreshold'] = 0.1; # // threshold for the ionisation [-]
#target
ParamsValue['trg.a'] = 1.189;
# units
ParamsValue['fieldinau'] = 1; # // 0-input inatomic units, 1 - in femto and GV/m
ParamsValue['input0'] = 1; # // 0 field starts with 0, 1 in the middle of the file
# definitions of sin^2 for an analytic field
ParamsValue['sin2.omega0'] = (2.0*np.pi*hbar*inverse_alpha_fine**2)/(BeamSpecification['lambda']*elmass*c_light) # find frequency in atomic units   #  frequency [a.u.]\n"
ParamsValue['sin2.ti'] = 0.0;    #  initial time\n"
ParamsValue['sin2.Ncyc'] = 10.0;   #  # of cycles\n"
# GAUGES
ParamsValue['gauge'] = 0;  #; // 0-length, otherwise velocity, velocity available only for analytic field (A needed)
ParamsValue['transformgauge']= 0; #		dumint=fscanf(param,"%i %*[^\n]\n",&transformgauge); // 1 - transform also to another gauge during the calculation, (A needed)

print(ParamsValue)

f = h5py.File('data.h5','w')
### save constant parameters to scalar datasets in hdf5 file




#### THE MAIN PROGRAM #####
#h5py.run_tests()


## create input file for 1DTDSE
#print('-----------------Creating param file-----------------')


#f = h5py.File('data.h5','w')

## put listing file in the root folder
#listingfile='This is the content of the listing file \nsecond row'
#f.create_dataset('listing', data=[np.string_(listingfile)])

#grp = f.create_group('micro')
## we put in shared parameters for microscpic solvers, they're either floats or integers
#params_real=[np.pi,np.e,5.0]; params_int=[0, 1, 2];
#grp.create_dataset('SharedParamsReal',data=params_real); grp.create_dataset('SharedParamsInt',data=params_int)



## We have 3 dimensions: z, r, t. We know numbers of points in r and t, but we will be appending in z. Based on the microscopic model, last dimension may contain variable parameters for an analytic model
#Nr=2; Nt=3
#dset = grp.create_dataset('FieldsForTDSE', (1,Nr,Nt), dtype='f', maxshape=(None,Nr,Nt)) # we start with a proper number of points
#dset.dims[0].label = 'z [SI]'; dset.dims[1].label = 'r [SI]'; dset.dims[2].label = 't';
#dset.attrs['stored quantity']=np.string_('Electric field [a.u.]') # for avoiding confusion while reading data, we can label each dataset by the used units


##We also provide grids for computations. I'd suggest to use atomic units for time and SI for the rest
#tgrid = (0.0, 2.0, 4.0)
#tgrid = grp.create_dataset('tgrid', data=tgrid)
#tgrid.attrs['units']=np.string_('[SI]')

#rgrid = (0.0, 5.0)
#rgrid = grp.create_dataset('rgrid', data=rgrid)
#rgrid.attrs['units']=np.string_('[SI]')

#zgrid=grp.create_dataset('zgrid', (1,), dtype='f', maxshape=(None,))
#zgrid.attrs['units']=np.string_('[SI]')

## now, I introduce a serial solution of the problem, MPI will be added later
## we are in the first plane and fill some data
#for k1 in range(2): dset[0,k1,:] = np.array([float(k1)+1.0, float(k1)+2.0, float(k1)+3.0]) # the fields are writen now
#zgrid[0]=-5.5;


#print('shape of the array is:', dset.shape) 
#f.close # we close the file to ensure saving

## reopen in the next plane
#f = h5py.File('results.h5','r+') # know file exists and we will write to it
#dset = f['micro/FieldsForTDSE'] # open our dataset
#newshape = np.asarray(dset.shape); newshape[0]=newshape[0]+1 # find new shape of the array and extend the first axis
#dset.resize(newshape) # resize the dataset
#for k1 in range(2): dset[1,k1,:] = np.array([float(k1)+5.0, float(k1)+6.0, float(k1)+7.0]) # fill some data
#print('shape of the array is:', dset.shape)

#zgrid = f['micro/zgrid']; zgrid.resize( ((zgrid.len()+1),) ); zgrid[1]=-4.5; # one-line for extending zgrid
#f.close()


#print('-----------------------------------------------------')
