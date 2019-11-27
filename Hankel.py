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
#import h5py


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

def FindInterval(x,x0): # find an index corresponding to given x0 value interval. ordering <  ), < ),..., < >; throws error otherwise
  N = len(x)
  for k1 in range(N-2):
    if ( (x[k1]<= x0) and (x0 < x[k1+1]) ): return k1
  if ( (x[N-2]<= x0) and (x0 <= x[N-1]) ): return N-2
  sys.exit('out of range in FindInterval')



### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factor to atomic units
TIMEau = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);
INTENSITYau = (inverse_alpha_fine/(8.0*np.pi))*(hbar**3)/((elmass**2)*(r_Bohr**6));



#### THE MAIN PROGRAM #####

#print("Number of procs: ", mp.cpu_count()) # Actually, it should be adjusted by slurm-scheduler, it probably sees "physical HW" and not the allocated resources


## parameters

inpath = os.path.join('sims11','z_000002') # path for TDSEs
inpath2 = 'sims11' # path for fields
outpath = 'res1' #'res8-2-dr16rm4' # path for results -dr2dr2rm2
if os.path.exists(outpath) and os.path.isdir(outpath):
  shutil.rmtree(outpath)
  print('deleted previous results')
os.mkdir(outpath)

## the type of the field
MicroscopicModelType = 1; # 0-TDSE, 1 -phenomenological 

omega0 = 0.057; # [a.u.]
TFWHM = 50e-15; # [SI]
tcoeff = 2.0; # extension of tgrid 
omegawidth = 4.0/np.sqrt(4000.0**2); # roughly corresponds to 100 fs
I0 = 2.5e18;
w0 = 96e-6;
PhenomParams = np.array([
[29, 35, 39], # harmonics
[500., 1775., 3600.], # alphas
[omegawidth, omegawidth, omegawidth]
])
NumHarm = 3; # numbero of harmonics

#print(PhenomParams)
#quit()

## parameters of the screen, etc.
rmax_anal = 0.02; # [SI] on screen # 0.0001
Nr_anal=200; #750
D = 3.0 # [SI], screen distance 1

omegamin_anal = 0.057*34.0 ;
omegamax_anal = 0.057*36.0 # 0.057*40.0 # 0.057*55.0
omega_step = 1

## numerical params
Nr_step = 1 # reshape the grid for the integration in r usw every Nr_step point
rIntegrationFactor = 1.0/2.0;
#rIntegrationFactorMin = 1.0/16.0;
rIntegrationFactorMin = 0; # not implemented yet, need to redefine the integration function

## other parameters
integrator = 0; #(0 - trapezoidal, 1 - Simpson) 
W = mp.cpu_count() # this is the number of workers
W = 1;
  
  
if (MicroscopicModelType == 0):
  print('numerical model')
  ## PUT LOADING the dipoles FROM MODULE HERE
elif (MicroscopicModelType == 1):
#  omegagrid = np.linspace(0.0,0.057*100.0,10000)
  print('analytical model')

  rgrid = np.linspace(0.0,2*w0,1000)
  tgrid = np.linspace(-tcoeff*0.5*TFWHM/TIMEau,tcoeff*0.5*TFWHM/TIMEau,10000)
  Nomega = len(tgrid)//2 + 1
  omegagrid = np.linspace(0,2.0*np.pi*Nomega/(tcoeff*TFWHM/TIMEau),Nomega)
  Nr = len(rgrid)
  FField_r = []; # computed on the fly in this case

  Nomega_anal = FindInterval(omegagrid,omegamax_anal) #3000 3000 1000
  Nomega_anal_start = FindInterval(omegagrid,omegamin_anal)
  print('om_max', Nomega_anal)
  print('om_min', Nomega_anal_start)


  



# reshape rgrid if all points are not used
if ( Nr_step != 1):
  rgridnew=[]
  for k1 in range(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step): rgridnew.append(rgrid[k1])
  rgrid=np.asarray(rgridnew);
  Nr = NumOfPointsInRange(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step)
  
dr = rgrid[1]-rgrid[0]


## define dipole function
def dipoleTimeDomainApp(tgrid,r,I0,PhenomParams,tcoeff,rcoeff,omega0): # some global variables involved
#  tcoeff = 4.0*np.log(2.0)*TIMEau**2 / ( TFWHMSI**2 )
#  rcoeff = 2.0/(w0r**2)
  res = []
  for k1 in range(len(tgrid)):
    res1 = 0.0*1j;
    intens = I0*np.exp(-tcoeff*(tgrid[k1])**2 - rcoeff*r**2)
    for k2 in range(NumHarm): res1 = res1 + intens*np.exp(1j*(omega0*PhenomParams[0,k2]-PhenomParams[1,k2]*intens)) 
    res.append(res1); ## various points in time
  return np.asarray(res)


  


## grids for analyses:
rgrid_anal = np.linspace(0.0,rmax_anal,Nr_anal)
Nomega_points = NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform
omegagrid_anal=[]

print('Nomega_points = ', Nomega_points);


# create dipole grid in r and omega
if (MicroscopicModelType == 1):
  print('Computing FFTs');
  FField_r=np.empty([Nomega,Nr], dtype=np.cdouble)

  tcoeff = 4.0*np.log(2.0)*TIMEau**2 / ( TFWHM**2 )
  rcoeff = 2.0/(w0**2)

  for k1 in range(Nr):    
    dum = dipoleTimeDomainApp(tgrid,rgrid[k1],I0,PhenomParams,tcoeff,rcoeff,omega0)
    dum = np.fft.fft(dum)
    for k2 in range(Nomega):
      FField_r[k2,k1] = dum[k2]

  print('FFT computed');
  


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
      k_omega =  omegagrid[k5]/(TIMEau*c_light); # omega divided by time: a.u. -> SI
      integrand = np.empty([Nr], dtype=np.cdouble)
      if ( (MicroscopicModelType == 0) or (MicroscopicModelType == 1) ):
        for k3 in range(Nr): integrand[k3] = np.exp(-(rgrid[k3]**2)/(2.0*D))*rgrid[k3]*FField_r[k5,k3]*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/D); # rescale r to atomic units for spectrum in atomic units! (only scaling)
      elif (MicroscopicModelType == 2): # we compute now fftw for evwery harmonic independently and then shift it... we shoul use linearity instead
        print('pure gaussian maybe? Now do nothing.')
#        for k3 in range(NumHarm): dip
#        for k3 in range(Nr): integrand[k3] = np.exp(-(rgrid[k3]**2)/(2.0*D))*rgrid[k3]*dipoleph(omegagrid[k5],omega0,I0rprofile(rgrid[k3],I0,w0),PhenomParams)*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/D); # rescale r to atomic units for spectrum in atomic units! (only scaling)
      if (integrator == 0):
        FHHGOnScreen[k4,k2] = integrate.trapz(integrand,rgrid);
      elif (integrator == 1):
        FHHGOnScreen[k4,k2] = integrate.simps(integrand,rgrid);
    toc = time.clock()
#    print('cycle',k1,'duration',toc-tic)
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
np.savetxt(os.path.join(outpath,"Spectrumreal.dat"),FHHGOnScreen.real,fmt="%e")
np.savetxt(os.path.join(outpath,"Spectrumimag.dat"),FHHGOnScreen.imag,fmt="%e")
np.savetxt(os.path.join(outpath,"omegagrid_anal.dat"),omegagrid_anal,fmt="%e")
np.savetxt(os.path.join(outpath,"rgrid_anal.dat"),rgrid_anal,fmt="%e")





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

