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


#ray.init()

## functions to work with indices
def n1n2mapping(k1,k2,N1): ## a mapping from (0,...,N1-1) x (0,...,N2-1) -> (0,...,N1*N2)
  return k1+k2*N1

def n1n2mapping_inv(n,N1): ## inverse to the previous one
  k1 = n % N1
  n = n - k1
  k2 = n // N1
  return k1,k2


def NumOfPointsInRange(N1,N2,k): #number of points between two integers following the Python's range-function logic (0,...,N-1), assumes N1<N2
  if N1 != 0:
    return NumOfPointsInRange(N1-N1,N2-N1,k);
  else:
    if (N2 % k) == 0:
      return N2 // k;
    else:
      return (N2 // k) + 1;


### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factor to atomic units
TIME = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);



#### THE MAIN PROGRAM #####

#print("Number of procs: ", mp.cpu_count()) # Actually, it should be adjusted by slurm-scheduler, it probably sees "physical HW" and not the allocated resources


## parameters

inpath = '35c/TDSEs2/' # path for TDSEs
inpath2 = '35c/fields/' # path for fields
outpath = 'res1/' # path for results


## parameters of the screen, etc.
rmax_anal = 0.00002; # [SI] on screen
Nr_anal=5;
D = 1.0 # [SI], screen distance

Nomega_anal = 3000 #3000
Nomega_anal_start = 2980
omega_step = 1

## other parameters
integrator = 0; #(0 - trapezoidal, 1 - Simpson) 
W = mp.cpu_count() # this is the number of workers
W = 4;


### LOAD GRIDS AND FIELDS

## load radial grid
file1=open(inpath2+"rgrid.dat","r")
#if file1.mode == "r":
lines = file1.readlines()
k1=0; rgrid=[]
for line in lines: dum = line.split(); rgrid.append(float(dum[1])); k1=k1+1
Nr = k1; rgrid=np.asarray(rgrid);  
file1.close()  


## retrieve the dimension of TDSEs
file1=open( (inpath+"z_000501_r_000001/GridDimensionsForBinaries.dat") ,"r")
lines = file1.readlines()
file1.close();
Nomega = int(lines[1]);
#print(Nomega)


## omega grid loaded directly in binary form
file1=open( (inpath+"z_000501_r_000001/omegagrid.bin") ,"rb")
omegagrid = array.array('d'); #[a.u.]
omegagrid.fromfile(file1,Nomega);
file1.close();


## read all fields in the binary form, it follows the padding andf naming of the files
# binary fields are stored in the omega domain (real(1),imaginary(1),real(2),imaginary(2),...)
Nfiles=Nr;
FField_r=np.empty([Nomega,Nfiles], dtype=np.cdouble)
for k1 in range(Nfiles):
  fold=str(k1+1); fold = inpath + 'z_000501_r_' + fold.zfill(6) + '/'
  FSourceTermPath = fold+'FSourceTerm.bin'
  file1=open(FSourceTermPath,"rb")
  dum = array.array('d');
  dum.fromfile(file1,2*Nomega);
  for k2 in range(Nomega): FField_r[k2,k1] = dum[2*k2]+1j*dum[2*k2+1];
#  print(k1)



## grids for analyses:
rgrid_anal = np.linspace(0.0,rmax_anal,Nr_anal)
Nomega_points = NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform
omegagrid_anal=[]

print('Nomega_points = ', Nomega_points);




#### MAIN INTEGRATION ####
tic1 = time.clock()
ttic1 = time.time()

# define output queue
output = mp.Queue()


# define function to integrate, there are some global variables! ## THE OUTPUT IS IN THE MIX OF ATOMIC UNITS (field) and SI UNITS (radial coordinate + dr in the integral)
def FieldOnScreen(omegagrid, omega_step, rgrid, Nr, FField_r, rgrid_anal, k_start, k_num):
  FHHGOnScreen = np.empty([k_num,Nr_anal], dtype=np.cdouble)
  k4=0 # # of loops in omega 
  for k1 in range(k_num): #Nomega
    k5 = k_start + k1*omega_step # accesing the grid
    tic = time.clock()
    for k2 in range(Nr_anal): #Nomega
      k_omega =  omegagrid[k5]/(TIME*c_light); # omega divided by time: a.u. -> SI
      integrand = np.empty([Nr], dtype=np.cdouble)
      for k3 in range(Nr): integrand[k3] = rgrid[k3]*FField_r[k5,k3]*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/D); # rescale r to atomic units for spectrum in atomic units! (only scaling)
      if (integrator == 0):
        FHHGOnScreen[k4,k2] = integrate.trapz(integrand,rgrid);
      elif (integrator == 1):
        FHHGOnScreen[k4,k2] = integrate.simps(integrand,rgrid);
    toc = time.clock()
    print('cycle',k1,'duration',toc-tic)
    k4=k4+1
  res = (k_start,k_num,FHHGOnScreen)
  output.put(res)



# Optimal workload is obtained by the same amount of load for each woker if possible, eventually one extra task for the last worker. Otherwise (NOT OPTIMAL!!!), every worker takes more load and some workers may be eventually not used. An optimal routine would be either balance the load among the workers or employ some sophisticated  parallel scheduler.
if ( ( (Nomega_points % W)==0 ) or ( (Nomega_points % W)==1 ) ):
  Nint = Nomega_points//W; # the number of points in an interval (beside the last interval...); NumOfPointsInRange(0,Nomega_points,W);
  N_PointsGrid=[]; N_PointsForProcess=[];
  for k1 in range(W): N_PointsGrid.append(k1*Nint);
  N_PointsGrid.append(Nomega_points);
  for k1 in range(W): N_PointsForProcess.append(N_PointsGrid[k1+1]-N_PointsGrid[k1])
else:
  Nint = (Nomega_points//W) + 1; 
  N_PointsGrid=[]; N_PointsForProcess=[];
  for k1 in range(W+1):
    dum = k1*Nint
    if dum >= Nomega_points:
      print('dum', dum)
      N_PointsGrid.append(Nomega_points);
      W = k1;
      break;
    else:
      N_PointsGrid.append(dum);
  print(N_PointsGrid)
  for k1 in range(W): N_PointsForProcess.append(N_PointsGrid[k1+1]-N_PointsGrid[k1])

# optimal workload is now given by the number of processes

print('----')
print('process grids for workers: starting point + loads')
print(N_PointsGrid)
print(N_PointsForProcess)
print('----')

# now we need to retrieve the number of loops for each worker in general cases...


### we use multiprocessing by assigning each part of the load as one process
# define processes
processes = [mp.Process(target=FieldOnScreen, args=(omegagrid, omega_step, rgrid, Nr, FField_r, rgrid_anal, Nomega_anal_start+N_PointsGrid[k1], N_PointsForProcess[k1])) for k1 in range(W)]

# run processes
for p in processes: p.start();

# exit the completed processes
#for p in processes: p.join(); # officially recommended but causes some "dead-lock"-like stuff


toc1 = time.clock()
ttoc1 = time.time()
print('duration before merging 1',toc1-tic1)
print('duration before merging 1',ttoc1-ttic1)


# append results, eventually some more complicated code for directly put in the final matrix
results = [output.get() for p in processes]


toc2 = time.clock()
ttoc2 = time.time()
print('duration after merging 2',toc2-tic1)
print('duration after merging 2',ttoc2-ttic1)


# create omega grid foe analyses
for k1 in range(Nomega_anal_start,Nomega_anal,omega_step): omegagrid_anal.append(omegagrid[k1])
omegagrid_anal=np.asarray(omegagrid_anal);


# conceneate the results
FHHGOnScreen = np.empty([Nomega_points,Nr_anal], dtype=np.cdouble)
for k1 in range(W): # loop over unsorted results
  for k2 in range(results[k1][1]): #  # of omegas computed by this worker
    for k3 in range(Nr_anal): # results in the radial grid
      FHHGOnScreen[ results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k2][k3] # we adjust the field index properly to the field it just sorts matrices the following way [A[1], A[2], ...], the indices are retrieved by the append mapping
#      FHHGOnScreen[ results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k2][k3]/(r_Bohr**2) # eventually fully in atomic units for the integral, but there is still a prefactor of the integral!!!
  



#### PROCESS OUTPUTS ####

    
#file1=open("Spectrum.dat","w")
np.savetxt(outpath+"Spectrumreal.dat",FHHGOnScreen.real,fmt="%e")
np.savetxt(outpath+"Spectrumimag.dat",FHHGOnScreen.imag,fmt="%e")
np.savetxt(outpath+"omegagrid_anal.dat",omegagrid_anal,fmt="%e")
np.savetxt(outpath+"rgrid_anal.dat",rgrid_anal,fmt="%e")




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

