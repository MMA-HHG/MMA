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


####
# Jan Vabek - ELI-Beamlines, CELIA, CTU in Prague (FNSPE)

# the purpose of that code is provide basic introduction to mpi h5py; next, I present the hierarchy of the data storage

### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factors to atomic units
TIME = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);



#### THE MAIN PROGRAM #####
#h5py.run_tests()


## first, we test MPI
comm = MPI.COMM_WORLD
rank = MPI.COMM_WORLD.rank
Nworkers = comm.Get_size()
print('My rank is: ', rank)
comm.Barrier() # synchro of workers
if (rank == 0): print('Number of workers: ', Nworkers)



## a basic example of HDF5: store just predefined arrays

if (rank == 0): # now, play with dimensions
  print('-----------------Example 1-----------------')
  d1 = np.random.random(size=(10,5))
  d2 = np.random.random(size=(5,2))
  print('Shape of the first array is: ', d1.shape)
  print('Shape of the second array is: ', d2.shape)
  print('Content of the second array is:\n', d2)

  hf = h5py.File('example1-basic.h5','w') # store data in an h5 file
  hf.create_dataset('data1',data=d1)
  hf.create_dataset('data2',data=d2)
  hf.close()

  print('-------------------------------------------')





## now a more advanced example showing the functionality in a combination with mpi
if (rank == 0): print('-----------------Example 2-----------------')
f = h5py.File('example2-parallel.h5', 'w', driver='mpio', comm=MPI.COMM_WORLD) # we need to provide the flags to say that it's a mpi-hdf5 file
## ACCESS VARIABLES
#r 	Readonly, file must exist
#r+ 	Read/write, file must exist
#w 	Create file, truncate if exists
#w- or x 	Create file, fail if exists
#a 	Read/write if exists, create otherwise (default)


# sorting data in groups
grp = f.create_group('group') 


# each worker writes its ID in the dataset
dset = grp.create_dataset('dataset1', (Nworkers,), dtype='i') 
dset[rank] = rank
f.close()


# now we reopen the file and only read data
if (rank == 0):
  f2 = h5py.File('example2-parallel.h5', 'r')
  content=list(f2['group'].keys()) # this lists all the datasets contained in a group
  print('datasets in the group are:', content)
  stored_data = f2['group']['dataset1'][:] # list-like access to data
  print('contained data are: ', stored_data)
  f2.close()
  
  # reopen and append some data dirextly in the root group
  f3 = h5py.File('example2-parallel.h5', 'a')
  dset = f3.create_dataset('dataset2', (3,), dtype='float') #float64, float, float128
  dset[0]=0.0
  dset[1]=np.pi
  dset[2]=np.e

  f3.close()
  

comm.Barrier() # synchro of workers
# now test multiple readers 

#f4 = h5py.File('example2-parallel.h5', 'r')
f4 = h5py.File('example2-parallel.h5', 'r', libver='latest', swmr=True) # now works with both, test some big data. The problem is probably with a writer and the readers. It then has to be noted to the writer that there could be readers accessing the file... See documentation.
stored_data = f4['group/dataset1'][:] # path-like access to data

mydata=stored_data[rank] # each worker reads data corresponding to its index
print('worker: ',rank, 'my data: ', mydata)
f4.close()

comm.Barrier()
if (rank == 0): print('-------------------------------------------')
comm.Barrier() # synchro of workers





 
# now, the implementation in our case
if (rank == 0): print('-----------------Example 3-----------------')
if (rank == 0): 
  f = h5py.File('results.h5','w')

  # put listing file in the root folder
  listingfile='This is the content of the listing file \nsecond row'
  f.create_dataset('listing', data=[np.string_(listingfile)])

  grp = f.create_group('micro')
  # we put in shared parameters for microscpic solvers, they're either floats or integers
  params_real=[np.pi,np.e,5.0]; params_int=[0, 1, 2];
  grp.create_dataset('SharedParamsReal',data=params_real); grp.create_dataset('SharedParamsInt',data=params_int)
  
  

  # We have 3 dimensions: z, r, t. We know numbers of points in r and t, but we will be appending in z. Based on the microscopic model, last dimension may contain variable parameters for an analytic model
  Nr=2; Nt=3
  dset = grp.create_dataset('FieldsForTDSE', (1,Nr,Nt), dtype='f', maxshape=(None,Nr,Nt)) # we start with a proper number of points
  dset.dims[0].label = 'z [SI]'; dset.dims[1].label = 'r [SI]'; dset.dims[2].label = 't';
  dset.attrs['stored quantity']=np.string_('Electric field [a.u.]') # for avoiding confusion while reading data, we can label each dataset by the used units
  

  #We also provide grids for computations. I'd suggest to use atomic units for time and SI for the rest
  tgrid = (0.0, 2.0, 4.0)
  tgrid = grp.create_dataset('tgrid', data=tgrid)
  tgrid.attrs['units']=np.string_('[SI]')

  rgrid = (0.0, 5.0)
  rgrid = grp.create_dataset('rgrid', data=rgrid)
  rgrid.attrs['units']=np.string_('[SI]')

  zgrid=grp.create_dataset('zgrid', (1,), dtype='f', maxshape=(None,))
  zgrid.attrs['units']=np.string_('[SI]')

  # now, I introduce a serial solution of the problem, MPI will be added later
  # we are in the first plane and fill some data
  for k1 in range(2): dset[0,k1,:] = np.array([float(k1)+1.0, float(k1)+2.0, float(k1)+3.0]) # the fields are writen now
  zgrid[0]=-5.5;
  

  print('shape of the array is:', dset.shape) 
  f.close # we close the file to ensure saving

  # reopen in the next plane
  f = h5py.File('results.h5','r+') # know file exists and we will write to it
  dset = f['micro/FieldsForTDSE'] # open our dataset
  newshape = np.asarray(dset.shape); newshape[0]=newshape[0]+1 # find new shape of the array and extend the first axis
  dset.resize(newshape) # resize the dataset
  for k1 in range(2): dset[1,k1,:] = np.array([float(k1)+5.0, float(k1)+6.0, float(k1)+7.0]) # fill some data
  print('shape of the array is:', dset.shape)
  
  zgrid = f['micro/zgrid']; zgrid.resize( ((zgrid.len()+1),) ); zgrid[1]=-4.5; # one-line for extending zgrid
  f.close()


