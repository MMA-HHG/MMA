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

#ray.init()


def n1n2mapping(k1,k2,N1): 
  return k1+k2*N1

def NumOfPointsInRange(N1,N2,k): #number of points between two integers following the Python's range logic (0,...,N-1), assumes N1<N2
  if N1 != 0:
    return NumOfPointsInRange(N1-N1,N2-N1,k);
  else:
    if (N2 % k) == 0:
      return N2 // k;
    else:
      return (N2 // k) + 1;

print("Number of procs: ", mp.cpu_count())

def myfun():
  print("abc")

### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;

# conversion factor to atomic units
TIME = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);
print(TIME)

inpath = '35c/TDSEs2/'


print("Hello world")

myfun()

a = special.jn(0,1)

print(a)

print("J0=",a)


##file1=open("Spectrum.dat","w")
#xxx = np.array([[1.0e25, 2.0], [3.0, 5.0]])
#print(xxx)
#print(xxx[0,1])
##file1.fprint(xxx)
#np.savetxt("Spectrum.dat",xxx,fmt="%e")
##np.savetxt("Spectrum.dat",xxx)
##file1.close();

file1=open("35c/fields/rgrid.dat","r")
#if file1.mode == "r":
lines = file1.readlines()
k1=0; rgrid=[]
for line in lines: dum = line.split(); rgrid.append(float(dum[1])); k1=k1+1
Nr = k1; rgrid=np.asarray(rgrid);  
file1.close()  


print(rgrid)
#for k1 in range(Nr): print(z[k1])


#file1=open("GridDimensionsForBinaries.dat","r")
file1=open( (inpath+"z_000501_r_000001/GridDimensionsForBinaries.dat") ,"r")
lines = file1.readlines()
file1.close();
Nomega = int(lines[1]);
print(Nomega)



file1=open( (inpath+"z_000501_r_000001/omegagrid.bin") ,"rb")
omegagrid = array.array('d'); #[a.u.]
omegagrid.fromfile(file1,Nomega);
file1.close();
#xxx = struct.unpack('d',file1.read(8))
#print(xxx)

file1=open("FSourceTerm.bin","rb")
xxx = array.array('d');
xxx.fromfile(file1,2*Nomega);
#xxx = struct.unpack('d',file1.read(8))
file1.close();
#print(xxx)
#for k1 in range(Nomega):
#  print(z[k1])

y=np.empty([Nomega], dtype=np.cdouble)
for k1 in range(Nomega): y[k1] = xxx[2*k1]+1j*xxx[2*k1+1];


print(y[Nomega-1])


folders = []

# r=root, d=directories, f = files
#for r, d, f in os.walk("."):
#    for folder in d:
#        folders.append(os.path.join(r, folder))

#for f in folders:
#    print(f)


dirs=os.listdir("35c/TDSEs2")
#print(dirs)
print(dirs[0])


#tstr='xx'
tstr=str(1)
tstr='35c/z_000501_'+tstr.zfill(6)

print(tstr)

Nfiles=Nr;
FField_r=np.empty([Nomega,Nfiles], dtype=np.cdouble)
for k1 in range(Nfiles):
  fold=str(k1+1); fold = inpath + 'z_000501_r_' + fold.zfill(6) + '/'
  FSourceTermPath = fold+'FSourceTerm.bin'
  file1=open(FSourceTermPath,"rb")
  dum = array.array('d');
  dum.fromfile(file1,2*Nomega);
  for k2 in range(Nomega): FField_r[k2,k1] = dum[2*k2]+1j*dum[2*k2+1];
  print(k1)



print(FField_r[0,0])

#print(complex(xxx[1],z[1]))


rmax_anal = 0.00002; # [SI]
Nr_anal=100;
D = 1.0 # [SI], screen distance
rgrid_anal = np.linspace(0.0,rmax_anal,Nr_anal)
Nomega_anal = 3000 #3000
Nomega_anal_start = 2700
omega_step = 2

Nomega_points = NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step);

omegagrid_anal=[]

FHHGOnScreen = np.empty([Nomega_points,Nr_anal], dtype=np.cdouble)
integrand = np.empty([Nr], dtype=np.cdouble)


pool = mp.Pool(mp.cpu_count())

## main integration, first list omegas
k4=0 # # of loops in omega 
for k1 in range(Nomega_anal_start,Nomega_anal,omega_step): #Nomega
  tic = time.clock()
  for k2 in range(Nr_anal): #Nomega
    k_omega =  omegagrid[k1]/(TIME*c_light); # omega divided by time: a.u. -> SI
    for k3 in range(Nr): integrand[k3] = rgrid[k3]*FField_r[k1,k3]*special.jn(0,k_omega*rgrid[k3]*rgrid_anal[k2]/D); # rescale r to atomic units for spectrum in atomic units! (only scaling)
#    integrand = 
    FHHGOnScreen[k4,k2] = integrate.trapz(integrand,rgrid);
#    FHHGOnScreen[k4,k2] = integrate.simps(integrand,rgrid);
  toc = time.clock()
  print('cycle',k1,'duration',toc-tic)
  omegagrid_anal.append(omegagrid[k1]);
  k4=k4+1

omegagrid_anal=np.asarray(omegagrid_anal);

    
#file1=open("Spectrum.dat","w")
np.savetxt("Spectrumreal.dat",FHHGOnScreen.real,fmt="%e")
np.savetxt("Spectrumimag.dat",FHHGOnScreen.imag,fmt="%e")
np.savetxt("omegagrid_anal.dat",omegagrid_anal,fmt="%e")
np.savetxt("rgrid_anal.dat",rgrid_anal,fmt="%e")











