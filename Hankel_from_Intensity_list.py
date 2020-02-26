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
#from mpi4py import MPI
#import oct2py
import shutil
import h5py
import sys
import units
import mynumerics as mn
import Hfn

#the plan is to
# - load all data from the HDF5 intensity list
# - specify the rest in the code now
# - generate the Gaussian profile at any point
# - add the extra phase
# - use old procedure for the integration in r now
# - we use precomputed dipoles now

#ray.init()


#### THE MAIN PROGRAM #####

#print("Number of procs: ", mp.cpu_count()) # Actually, it should be adjusted by slurm-scheduler, it probably sees "physical HW" and not the allocated resources


###################### THE PARAMETERS OF SIMULATION
#inpath = os.path.join('sims11','z_000002') # path for TDSEs

# IntensityListFile = os.path.join("C:\data","ThinTargets_collab","DipoleIntensityTable_1k.h5")
IntensityListFile = os.path.join("/mnt","c","data","ThinTargets_collab","DipoleIntensityTable_1k.h5")
# IntensityListFile = 'ThinDipoleIntensityTable_5k.h5' # path for fields

# loading
file1 = h5py.File(IntensityListFile, 'r')
Igrid = file1['/Igrid'][:]
omegagrid = file1['/omegagrid'][:]
FSourceterm = file1['/FDipoleAccelerations'][:]
FSourceterm = np.squeeze(FSourceterm[:,:,0] + 1j*FSourceterm[:,:,1]) # convert to complex numbers


## params
LaserParams={ ## define macroscopic gaussian beam # try also fancy reading directly here
'I0' : 4.0e18,
'w0' : 96.0e-6,
'r_extend' : 4.0,
'E0' : 0.075, 
'z' : 0.05,
'lambda' : 800.0e-9, # must correspond
'phase0' : 0.0, # initial CEP
'TFWHM' : 40e-15 # [SI]
}
LaserParams['omega0'] = mn.ConvertPhoton(LaserParams['lambda'] ,'lambdaSI','omegaau') # find frequency in atomic units
LaserParams['zR'] = np.pi*((LaserParams['w0'])**2)/LaserParams['lambda']
omega0 = LaserParams['omega0']; zR = LaserParams['zR'];

# anlyses params # at the moment optimised for t he intensity list, change later

z_medium = -0.003;

rmax = 2.0*LaserParams['w0'];
Nr = 100;

rmax_anal = 0.3*0.008 # [SI] on screen # 0.0001
Nr_anal=100 #750

zmin_anal = 0.001 # !!!!!! now in the reference of the jet
zmax_anal = 0.2
Nz_anal = 200

Hmin_anal = 28.5
Hmax_anal = 29.5
omega_step = 1

## other parameters
integrator = 'Trapezoidal'; # 'Trapezoidal', Simpson
W = mp.cpu_count() # this is the number of workers
W = 8;


#print(PhenomParams)
print(omega0,'omega0 in a.u.')
#quit()


## create grids
Nomega = len(omegagrid)
Hgrid = np.empty([Nomega], dtype=np.double)
Hgrid[:] = omegagrid[:]/LaserParams['omega0']
# Hgrid = np.empty([Nomega], dtype=np.double)
# for k1 in range(Nomega): Hgrid[k1] = omegagrid

rgrid = np.linspace(0.0,rmax,Nr)
zgrid_anal = np.linspace(z_medium+zmin_anal,zmax_anal,Nz_anal)
rgrid_anal = np.linspace(0,rmax_anal,Nr_anal)
omegamin_anal = omega0*Hmin_anal
omegamax_anal = omega0*Hmax_anal
Nomega_anal = mn.FindInterval(omegagrid,omegamax_anal)
Nomega_anal_start = mn.FindInterval(omegagrid,omegamin_anal)

Nomega_points = mn.NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform



## compute fields in our rgrid


print('Computing fields in the grid');
FField_r=np.empty([Nomega,Nr], dtype=np.cdouble)

for k1 in range(Nr): # We use linear interpolation using the intensity-grid at the instant
  I_r, phase_r = mn.GaussianBeam(rgrid[k1],z_medium,0,LaserParams['I0']/units.INTENSITYau,LaserParams['w0'],1,LaserParams['lambda']) # a.u.
  phase_XUV = phase_r*Hgrid

  # find a proper interval in the Igrid, we use linear interp written by hand now
  k2 = mn.FindInterval(Igrid,I_r)
  weight1 = (Igrid[k2+1]-I_r)/(Igrid[k2+1]-Igrid[k2]); weight2 = (I_r-Igrid[k2])/(Igrid[k2+1]-Igrid[k2]);
  FField_r[:, k1] = np.exp(-1j*phase_XUV)*(weight1*FSourceterm[k2,:]+weight2*FSourceterm[k2,:]); # free-form works?


## print some analyses outputs
print('om_max', Nomega_anal)
print('om_min', Nomega_anal_start)
# print('tmax',tgrid[len(tgrid)-1])
print('omax',omegagrid[len(omegagrid)-1])

  
dr = rgrid[1]-rgrid[0]



omegagrid_anal=[]

print('Nomega_points = ', Nomega_points);

  

#########################################################################################
#### MAIN INTEGRATION ####
tic1 = time.process_time()
ttic1 = time.time()

# define output queue
output = mp.Queue()

# passing by reference is unPythonic, we define the extra function though
def FieldOnScreen_handle(z_medium, omegagrid, omega_step, rgrid, FField_r, rgrid_anal, zgrid_anal, k_start, k_num, integrator):
  res = Hfn.FieldOnScreen(z_medium, omegagrid, omega_step, rgrid, FField_r, rgrid_anal, zgrid_anal, k_start, k_num, integrator)
  output.put(res)


# split the workload evenly
N_PointsGrid, N_PointsForProcess = Hfn.ObtainWorkload(Nomega_points,W)

print('----')
print('process grids for workers: starting point + loads')
print(N_PointsGrid)
print(N_PointsForProcess)
print('----')

# now we need to retrieve the number of loops for each worker in general cases...


### we use multiprocessing by assigning each part of the load as one process
# define processes
processes = [mp.Process(target=FieldOnScreen_handle, args=(z_medium, omegagrid, omega_step, rgrid, FField_r, rgrid_anal, zgrid_anal, Nomega_anal_start+N_PointsGrid[k1], N_PointsForProcess[k1],integrator)) for k1 in range(W)]

# run processes
for p in processes: p.start();

# exit the completed processes
#for p in processes: p.join(); # officially recommended but causes some "dead-lock"-like stuff


toc1 = time.process_time()
ttoc1 = time.time()
print('duration before merging 1',toc1-tic1)
print('duration before merging 1',ttoc1-ttic1)


# append results, eventually some more complicated code for directly put in the final matrix
results = [output.get() for p in processes]


toc2 = time.process_time()
ttoc2 = time.time()
print('duration after merging 2',toc2-tic1)
print('duration after merging 2',ttoc2-ttic1)


# create omega grid foe analyses
for k1 in range(Nomega_anal_start,Nomega_anal,omega_step): omegagrid_anal.append(omegagrid[k1])
omegagrid_anal=np.asarray(omegagrid_anal);


# conceneate the results
FHHGOnScreen = np.empty([Nz_anal,Nomega_points,Nr_anal], dtype=np.cdouble)
for k1 in range(W): # loop over unsorted results
  for k2 in range(results[k1][1]): #  # of omegas computed by this worker
    for k3 in range(Nr_anal): # results in the radial grid
      for k4 in range(Nz_anal): # results in the z grid
        FHHGOnScreen[k4, results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k4][k2][k3] # we adjust the field index properly to the field it just sorts matrices the following way [A[1], A[2], ...], the indices are retrieved by the append mapping
#      FHHGOnScreen[ results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k2][k3]/(r_Bohr**2) # eventually fully in atomic units for the integral, but there is still a prefactor of the integral!!!
  



#################################################################################
#### PROCESS OUTPUTS ####

    
#file1=open("Spectrum.dat","w")
np.savetxt(os.path.join(outpath,"Spectrumreal.dat"),FHHGOnScreen[1,:,:].real,fmt="%e")
np.savetxt(os.path.join(outpath,"Spectrumimag.dat"),FHHGOnScreen[1,:,:].imag,fmt="%e")
np.savetxt(os.path.join(outpath,"omegagrid_anal.dat"),omegagrid_anal,fmt="%e")
np.savetxt(os.path.join(outpath,"rgrid_anal.dat"),rgrid_anal,fmt="%e")
np.savetxt(os.path.join(outpath,"zgrid_anal.dat"),zgrid_anal,fmt="%e")

#HDF5 results
f = h5py.File(os.path.join(outpath,"results.h5"),'a')
grp = f.create_group('XUV')

h5spectrum_r = grp.create_dataset('Spectrum_r', data=FHHGOnScreen.real)
#h5spectrum_r = grp.create_dataset('Spectrum_r', data=FHHGOnScreen.real,compression='gzip',compression_opts=9)
h5spectrum_r.attrs['units']=np.string_('[arb.u.]')
h5spectrum_r.dims[0].label = 'z [SI]'; h5spectrum_r.dims[1].label = 'omega [a.u.]'; h5spectrum_r.dims[2].label = 'r [SI]';

h5spectrum_i = grp.create_dataset('Spectrum_i', data=FHHGOnScreen.imag)
#h5spectrum_i = grp.create_dataset('Spectrum_i', data=FHHGOnScreen.imag,compression='gzip',compression_opts=9)
h5spectrum_i.attrs['units']=np.string_('[arb.u.]')

h5_ogrid_anal = grp.create_dataset('omegagrid_screen', data=omegagrid_anal)
h5_ogrid_anal.attrs['units']=np.string_('[a.u.]')

h5_rgrid_anal = grp.create_dataset('rgrid_screen', data=rgrid_anal)
h5_rgrid_anal.attrs['units']=np.string_('[SI]')

h5_zgrid_anal = grp.create_dataset('zgrid_screen', data=zgrid_anal)
h5_zgrid_anal.attrs['units']=np.string_('[SI]')

grp = f.create_group('params')
dum = grp.create_dataset('omega0', data=omega0)
dum.attrs['units']=np.string_('[a.u.]')

dum = grp.create_dataset('zR', data=zR)
dum.attrs['units']=np.string_('[SI]')

f.close()


## write params for reference
file1=open( os.path.join(outpath,'paramHankel.txt') ,"a")
content = "// parameters of Hankel transform\n"
content = content + str(rgrid[Nr-1]) + " : rmax for the integral [SI]\n"
content = content + str(rgrid[0]) + " : rmin for the integral [SI]\n"
content = content + str(dr) + " : dr for the integral [SI]\n"
content = content + str(Nr) + " : # of points in r\n"
file1.write(content)
file1.close()

# some graphical outputs directly?



















############# examples of some pieces of code
#file1=open("FSourceTerm.bin","rb")
#xxx = array.array('d');
#xxx.fromfile(file1,2*Nomega);
##xxx = struct.unpack('d',file1.read(8))
#file1.close();
##print(xxx)
##for k1 in range(Nomega):
##  print(z[k1])

#y=np.empty([Nomega], dtype=np.cdouble)
#for k1 in range(Nomega): y[k1] = xxx[2*k1]+1j*xxx[2*k1+1];


#print(y[Nomega-1])


#folders = []
# r=root, d=directories, f = files
#for r, d, f in os.walk("."):
#    for folder in d:
#        folders.append(os.path.join(r, folder))

#for f in folders:
#    print(f)





#dirs=os.listdir("35c/TDSEs2")
##print(dirs)
#print(dirs[0])


##tstr='xx'
#tstr=str(1)
#tstr='35c/z_000501_'+tstr.zfill(6)

#print(tstr)



#print(results[0][1])
#print(FHHGOnScreen[0,0])
#print(FHHGOnScreen[0,1])
#print(FHHGOnScreen[0,2])
#print(FHHGOnScreen[0,0])
#print(FHHGOnScreen[1,0])
#print(FHHGOnScreen[2,0])
#print(FHHGOnScreen[3,0])

