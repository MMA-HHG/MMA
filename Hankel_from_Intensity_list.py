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

#the plan is to
# - load all data from the HDF5 intensity list
# - specify the rest in the code now
# - generate the Gaussian profile at any point
# - add the extra phase
# - use old procedure for the integration in r now

#ray.init()


#### THE MAIN PROGRAM #####

#print("Number of procs: ", mp.cpu_count()) # Actually, it should be adjusted by slurm-scheduler, it probably sees "physical HW" and not the allocated resources


###################### THE PARAMETERS OF SIMULATION

inpath = os.path.join('sims11','z_000002') # path for TDSEs
inpath2 = 'sims11' # path for fields
outpath = 'res' #'res8-2-dr16rm4' # path for results -dr2dr2rm2
if os.path.exists(outpath) and os.path.isdir(outpath):
  shutil.rmtree(outpath)
  print('deleted previous results')
os.mkdir(outpath)

## the type of the field
MicroscopicModelType = 1; # 0-TDSE, 1 -phenomenological


 
LaserParams={ ## define macroscopic gaussian beam
'w0' : 96.0e-6,
'r_extend' : 4.0,
'E0' : 0.075, 
'z' : 0.05,
'lambda' : 800.0e-9,
'phase0' : 0.0, # initial CEP
'TFWHM' : 50e-15 # [SI]
}

LaserParams['omega0'] = (2.0*np.pi*units.hbar*units.inverse_alpha_fine**2)/(LaserParams['lambda']*units.elmass*units.c_light) # find frequency in atomic units
LaserParams['zR'] = np.pi*((LaserParams['w0'])**2)/LaserParams['lambda']

#omega0 = 0.057; # [a.u.]
omega0 = LaserParams['omega0']
TFWHM = LaserParams['TFWHM'] #50e-15; # [SI]
tcoeff = 4.5; # extension of tgrid 
#omegawidth = 4.0/np.sqrt(4000.0**2); # roughly corresponds to 100 fs

zR = LaserParams['zR']

I0 = 2.5e18;

z_medium = -0.005;

#w0 = 96e-6;
PhenomParams = np.array([
[1, 29, 35, 39], # harmonics
[0.0, 500., 1775., 3600.], # alphas
[3.0,3.0, 3.0, 3.0]
])
NumHarm = 4; # numbero of harmonics

#print(PhenomParams)
print(omega0,'omega0 in a.u.')
#quit()

## parameters of the screen, etc.
rmax_anal = 0.3*0.008; # [SI] on screen # 0.0001
Nr_anal=100; #750
#zgrid_anal = np.array([0.5, 1.0, 2.0, 3.0]);  #D=3.0 # [SI], screen distance 1
zgrid_anal = np.linspace(z_medium+0.001,0.2,200)
Nz_anal = len(zgrid_anal);

omegamin_anal = omega0*28.5 ;
omegamax_anal = omega0*29.5 # 0.057*40.0 # 0.057*55.0
omega_step = 1

## numerical params
Nr_step = 1 # reshape the grid for the integration in r usw every Nr_step point
rIntegrationFactor = 1.0; #1.0/2.0;
#rIntegrationFactorMin = 1.0/16.0;
rIntegrationFactorMin = 0; # not implemented yet, need to redefine the integration function

## other parameters
integrator = 0; #(0 - trapezoidal, 1 - Simpson) 
W = mp.cpu_count() # this is the number of workers
W = 8;
  
  


#################################### PREPARATORY CALCULATIONS

if (MicroscopicModelType == 0):
  print('numerical model')
  ## PUT LOADING the dipoles FROM MODULE HERE
elif (MicroscopicModelType == 1):
#  omegagrid = np.linspace(0.0,0.057*100.0,10000)
  print('analytical model')

  rgrid = np.linspace(0.0,LaserParams['w0'],100)
  tgrid = np.linspace(-tcoeff*0.5*TFWHM/units.TIMEau,tcoeff*0.5*TFWHM/units.TIMEau,10000)
  Nomega = len(tgrid)//2 + 1
  omegagrid = np.linspace(0,2.0*np.pi*Nomega/(tcoeff*TFWHM/units.TIMEau),Nomega)
  Nr = len(rgrid)
  FField_r = []; # computed on the fly in this case

  Nomega_anal = mn.FindInterval(omegagrid,omegamax_anal) #3000 3000 1000
  Nomega_anal_start = mn.FindInterval(omegagrid,omegamin_anal)
  print('om_max', Nomega_anal)
  print('om_min', Nomega_anal_start)
  print('tmax',tgrid[len(tgrid)-1])
  print('omax',omegagrid[len(omegagrid)-1])


# reshape rgrid if all points are not used
if ( Nr_step != 1):
  rgridnew=[]
  for k1 in range(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step): rgridnew.append(rgrid[k1])
  rgrid=np.asarray(rgridnew);
  Nr = mn.NumOfPointsInRange(int(round(rIntegrationFactorMin*Nfiles)),int(round(rIntegrationFactor*Nfiles)),Nr_step)
  
dr = rgrid[1]-rgrid[0]



## grids for analyses:
rgrid_anal = np.linspace(0.0,rmax_anal,Nr_anal)
#zgrid_anal = np.linspace(zmin_anal,zmax_anal,Nz_anal)
Nomega_points = mn.NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform
omegagrid_anal=[]

print('Nomega_points = ', Nomega_points);


# create dipole grid in r and omega # it's already loaded in the numeric case
if (MicroscopicModelType == 1):
  print('Computing FFTs');
  FField_r=np.empty([Nomega,Nr], dtype=np.cdouble)

  tcoeff = 4.0*np.log(2.0)*units.TIMEau**2 / ( TFWHM**2 )
  rcoeff = 2.0/((LaserParams['w0'])**2)
  dt = tgrid[1]-tgrid[0]

  for k1 in range(Nr):    
    dum = mn.dipoleTimeDomainApp(z_medium,tgrid,rgrid[k1],I0/units.INTENSITYau,PhenomParams,tcoeff,rcoeff,LaserParams)
    if (k1 == 0):
      np.savetxt(os.path.join(outpath,"tgrid.dat"),tgrid,fmt="%e")
      np.savetxt(os.path.join(outpath,"dipoler.dat"),dum.real,fmt="%e")
      np.savetxt(os.path.join(outpath,"dipolei.dat"),dum.real,fmt="%e")
    dum = (dt/np.sqrt(2.0*np.pi))*np.fft.fft(dum)
    if (k1 == 0):
      np.savetxt(os.path.join(outpath,"ogrid_full.dat"),omegagrid,fmt="%e")
      np.savetxt(os.path.join(outpath,"Fdipoler.dat"),dum.real,fmt="%e")
      np.savetxt(os.path.join(outpath,"Fdipolei.dat"),dum.imag,fmt="%e")
    for k2 in range(Nomega):
      FField_r[k2,k1] = dum[k2]

  print('FFT computed');
  

#########################################################################################
#### MAIN INTEGRATION ####
tic1 = time.process_time()
ttic1 = time.time()

# define output queue
output = mp.Queue()


# define function to integrate, there are some global variables! ## THE OUTPUT IS IN THE MIX OF ATOMIC UNITS (field) and SI UNITS (radial coordinate + dr in the integral)
def FieldOnScreen(z_medium, omegagrid, omega_step, rgrid, Nr, FField_r, rgrid_anal, zgrid_anal, k_start, k_num):
  FHHGOnScreen = np.empty([Nz_anal,k_num,Nr_anal], dtype=np.cdouble)
  k4=0 # # of loops in omega 
  for k1 in range(k_num): #Nomega
    k5 = k_start + k1*omega_step # accesing the grid
    tic = time.process_time()
    for k6 in range(Nz_anal):
      for k2 in range(Nr_anal): #Nr
        k_omega =  omegagrid[k5]/(units.TIMEau*units.c_light); # omega divided by time: a.u. -> SI
        integrand = np.empty([Nr], dtype=np.cdouble)
        if ( (MicroscopicModelType == 0) or (MicroscopicModelType == 1) ):
          for k3 in range(Nr): integrand[k3] = np.exp(-(rgrid[k3]**2)/(2.0*(zgrid_anal[k6]-z_medium)))*rgrid[k3]*FField_r[k5,k3]*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/(zgrid_anal[k6]-z_medium)); # rescale r to atomic units for spectrum in atomic units! (only scaling)
        elif (MicroscopicModelType == 2): # we compute now fftw for evwery harmonic independently and then shift it... we shoul use linearity instead
          print('pure gaussian maybe? Now do nothing.')
  #        for k3 in range(NumHarm): dip
  #        for k3 in range(Nr): integrand[k3] = np.exp(-(rgrid[k3]**2)/(2.0*D))*rgrid[k3]*dipoleph(omegagrid[k5],omega0,I0rprofile(rgrid[k3],I0,w0),PhenomParams)*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/D); # rescale r to atomic units for spectrum in atomic units! (only scaling)
        if (integrator == 0):
          FHHGOnScreen[k6,k4,k2] = (1.0/(zgrid_anal[k6]-z_medium))*integrate.trapz(integrand,rgrid);
        elif (integrator == 1):
          FHHGOnScreen[k6,k4,k2] = (1.0/(zgrid_anal[k6]-z_medium))*integrate.simps(integrand,rgrid);
  
      #k2 loop end
    #k6 loop end
    toc = time.process_time()
#    print('cycle',k1,'duration',toc-tic)
    k4=k4+1
  res = (k_start,k_num,FHHGOnScreen)
  output.put(res)



# Optimal workload is obtained by the same amount of load for each woker if possible, eventually one extra task for the last worker. Otherwise (NOT OPTIMAL!!!), every worker takes more load and some workers may be eventually not used. An optimal routine would be either balance the load among the workers or employ some sophisticated  parallel scheduler.
if ( ( (Nomega_points % W)==0 ) or ( (Nomega_points % W)==1 ) ):
  Nint = Nomega_points//W; # the number of points in an interval (beside the last interval...); mn.NumOfPointsInRange(0,Nomega_points,W);
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
processes = [mp.Process(target=FieldOnScreen, args=(z_medium, omegagrid, omega_step, rgrid, Nr, FField_r, rgrid_anal, zgrid_anal, Nomega_anal_start+N_PointsGrid[k1], N_PointsForProcess[k1])) for k1 in range(W)]

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

