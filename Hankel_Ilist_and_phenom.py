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
# - we add one loop over different medium positions

#ray.init()


#### THE MAIN PROGRAM #####

#print("Number of procs: ", mp.cpu_count()) # Actually, it should be adjusted by slurm-scheduler, it probably sees "physical HW" and not the allocated resources


###################### THE PARAMETERS OF SIMULATION

LaserParams={ ## define macroscopic gaussian beam # try also fancy reading directly here
'I0' : 4.0e18,
'w0' : 85.0e-6,
'r_extend' : 4.0,
'z' : 0.05,
'lambda' : 800.0e-9, # must correspond
'phase0' : 0.0, # initial CEP
'TFWHM' : 30e-15 # [SI]
}
LaserParams['omega0'] = mn.ConvertPhoton(LaserParams['lambda'] ,'lambdaSI','omegaau') # find frequency in atomic units
LaserParams['zR'] = np.pi*((LaserParams['w0'])**2)/LaserParams['lambda']
omega0 = LaserParams['omega0']; zR = LaserParams['zR'];

# anlyses params # at the moment optimised for t he intensity list, change later

z_medium = np.asarray([-0.025, -0.02, -0.015, -0.01, -0.005, 0.0, 0.01])  # np.array([-0.003, 0.0, 0.003]);

rmax = 3.0*LaserParams['w0'];
Nr = 300;

rmax_anal = 5*1e-3 # [SI] on screen # 0.0001
Nr_anal = 50 #750

zmin_anal = 0.05 # !!!!!! in the reference of the jet, the grid is then reshaped correctly
zmax_anal = 0.5
only_one_plane_bool = False # only zmax used # Nz_anal is overrun
fix_zmax_bool = False
Nz_anal = 3 #200

Hmin_anal = 28.5 # 0.0 #28.5 np.nan
Hmax_anal = 29.5 # 2.5 #29.5 71.5
omega_step = 1

# used only for phenomenological dipoles
tcoeff = 6.0; # extension of tgrid in the units of TFWHM
Nt = 1000;


## other parameters
integrator = 'Trapezoidal'; # 'Trapezoidal', Simpson
dipole_model = 'IntensityList' # 'IntensityList', Phenomenological
W = mp.cpu_count() # this is the number of workers
W = 10;



# outpath = os.path.join("/mnt", "jvabek", "ThinTargets_collab")
# IntensityListFile = os.path.join("/mnt", "jvabek", "ThinTargets_collab", "Ilists", "DipoleIntensityTable_1k.h5") # used only for the list

# local
# outpath = os.path.join("/mnt","c","data","ThinTargets_collab")
# IntensityListFile = os.path.join("/mnt","c","data","ThinTargets_collab","DipoleIntensityTable_1k.h5")# used only for the list

# curta
# outpath = os.path.join("/scratch","jvabek","optics-less-focusing","beams")
# IntensityListFile = os.path.join("/scratch","jvabek","optics-less-focusing","DipoleIntensityTable_1k.h5")# used only for the list

#occigen
outpath = os.path.join("/scratch","cnt0025","cli7594","vabekjan","ThinTargets_collab","beams")
IntensityListFile = os.path.join("/scratch","cnt0025","cli7594","vabekjan","ThinTargets_collab","DipoleIntensityTable_1k.h5")# used only for the list


OutputFileName = "results1.h5" # "results_phenom8.h5"

# IntensityListFile = 'ThinDipoleIntensityTable_5k.h5' # path for fields
# IntensityListFile = os.path.join("C:\data","ThinTargets_collab","DipoleIntensityTable_1k.h5")

#print(PhenomParams)
print(omega0,'omega0 in a.u.')
#quit()



########################################################## THE BODY OF THE PROGRAM


if (dipole_model == 'IntensityList'):
  # loading procedure
  file1 = h5py.File(IntensityListFile, 'r')
  Igrid = file1['/Igrid'][:]
  omegagrid = file1['/omegagrid'][:]
  FSourceterm = file1['/FDipoleAccelerations'][:]
  FSourceterm = np.squeeze(FSourceterm[:,:,0] + 1j*FSourceterm[:,:,1]) # convert to complex numbers

elif (dipole_model == 'Phenomenological'):
  tgrid = np.linspace(-tcoeff * 0.5 * LaserParams['TFWHM'] / units.TIMEau, tcoeff * 0.5 * LaserParams['TFWHM'] / units.TIMEau, Nt)
  Nomega = len(tgrid)//2 + 1
  omegagrid = np.linspace(0,2.0*np.pi*Nomega/(tcoeff*LaserParams['TFWHM']/units.TIMEau),Nomega)





## create grids
Nomega = len(omegagrid)
Hgrid = np.empty([Nomega], dtype=np.double)
Hgrid[:] = omegagrid[:]/LaserParams['omega0']
# Hgrid = np.empty([Nomega], dtype=np.double)
# for k1 in range(Nomega): Hgrid[k1] = omegagrid

Nz_medium = len(z_medium)

rgrid = np.linspace(0.0,rmax,Nr)
if only_one_plane_bool:
  Nz_anal = 1;
  zgrid_anal = np.empty([Nz_medium, 1], dtype=np.double)
  if fix_zmax_bool:
    for k1 in range(Nz_medium): zgrid_anal[k1, :] = np.asarray([zmax_anal])  # zgrid_anal is thus in the IR-focus reference frame
  else:
    for k1 in range(Nz_medium): zgrid_anal[k1, :] = np.asarray([z_medium[k1]+zmax_anal])  # zgrid_anal is thus in the IR-focus reference frame
else:
  zgrid_anal = np.empty([Nz_medium, Nz_anal],dtype=np.double)
  if fix_zmax_bool:
    for k1 in range(Nz_medium): zgrid_anal[k1, :] = np.linspace(z_medium[k1] + zmin_anal, zmax_anal, Nz_anal)  # zgrid_anal is thus in the IR-focus reference frame
  else:
    for k1 in range(Nz_medium): zgrid_anal[k1, :] = np.linspace(z_medium[k1] + zmin_anal, z_medium[k1] + zmax_anal, Nz_anal)  # zgrid_anal is thus in the IR-focus reference frame


rgrid_anal = np.linspace(0,rmax_anal,Nr_anal)

if (np.isnan(Hmin_anal)):
  omegamin_anal = omegagrid[0]
  Nomega_anal_start = 0
else:
  omegamin_anal = omega0*Hmin_anal
  Nomega_anal_start = mn.FindInterval(omegagrid,omegamin_anal)

if (np.isnan(Hmax_anal)):
  omegamax_anal = omegagrid[len(omegagrid)-1]
  Nomega_anal = len(omegagrid)-1
else:
  omegamax_anal = omega0 * Hmax_anal
  Nomega_anal = mn.FindInterval(omegagrid, omegamax_anal)

Nomega_points = mn.NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform



## compute fields in our rgrid
if (dipole_model == 'IntensityList'): FField_r = Hfn.ComputeFieldsInRFromIntensityList(z_medium, rgrid, Hgrid, Nomega, LaserParams, Igrid, FSourceterm)
elif (dipole_model == 'Phenomenological'): FField_r = Hfn.ComputeFieldsPhenomenologicalDipoles(LaserParams['I0'],omega0,LaserParams['TFWHM'],LaserParams['w0'],tgrid,omegagrid,rgrid,z_medium)

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


# create omega grid for analyses
for k1 in range(Nomega_anal_start,Nomega_anal,omega_step): omegagrid_anal.append(omegagrid[k1])
omegagrid_anal=np.asarray(omegagrid_anal);


# coalesce the results
FHHGOnScreen = Hfn.CoalesceResults(results,Nz_medium,Nz_anal,Nomega_anal_start,Nomega_points,Nr_anal,W)



#################################################################################
#### PROCESS OUTPUTS ####

    
#file1=open("Spectrum.dat","w")
np.savetxt(os.path.join(outpath,"Spectrumreal.dat"),FHHGOnScreen[0,0,:,:,0].real,fmt="%e")
np.savetxt(os.path.join(outpath,"Spectrumimag.dat"),FHHGOnScreen[0,0,:,:,1].imag,fmt="%e")
np.savetxt(os.path.join(outpath,"omegagrid_anal.dat"),omegagrid_anal,fmt="%e")
np.savetxt(os.path.join(outpath,"rgrid_anal.dat"),rgrid_anal,fmt="%e")
np.savetxt(os.path.join(outpath,"zgrid_anal.dat"),zgrid_anal[0,:],fmt="%e")

#HDF5 results
f = h5py.File(os.path.join(outpath,OutputFileName),'w')
grp = f.create_group('XUV')

h5spectrum_r = grp.create_dataset('Spectrum', data=FHHGOnScreen)
#h5spectrum_r = grp.create_dataset('Spectrum_r', data=FHHGOnScreen.real,compression='gzip',compression_opts=9)
h5spectrum_r.attrs['units']=np.string_('[arb.u.]')
h5spectrum_r.dims[0].label = 'z [SI]'; h5spectrum_r.dims[1].label = 'omega [a.u.]'; h5spectrum_r.dims[2].label = 'r [SI], real [-]'; h5spectrum_r.dims[3].label = 'r [SI], imag [-]';

## grids
mn.adddataset(grp,'omegagrid_screen',omegagrid_anal,'[a.u.]')
mn.adddataset(grp,'rgrid_screen',rgrid_anal,'[SI]')
mn.adddataset(grp,'zgrid_screen',zgrid_anal,'[SI]')
mn.adddataset(grp,'ThinTargetPositions',z_medium,'[SI]')

## parameters
grp = f.create_group('params')
mn.adddataset(grp,'omega0',omega0,'[a.u.]')
mn.adddataset(grp,'zR',zR,'[SI]')
mn.adddataset(grp,'w0',LaserParams['w0'],'[SI]')
mn.adddataset(grp,'I0',LaserParams['I0'],'[SI]')



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


print('done')


















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

