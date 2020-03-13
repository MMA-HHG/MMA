from scipy import special
from scipy import integrate
import numpy as np
import math
import sys
import units
import h5py


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
  # if ( (x0 < x[0]) or (x0 > x[N-1])): sys.exit('out of range in FindInterval')
  # if (N == 2): return k1; # we finished looking
  # else:
  #   if (x0 > x[N//2] ): return FindInterval(x[(N//2):(N-1)],x0,N//2+?); # bookkeeping needed here... best will be additions and subtractions to be in-place
  #   else : return FindInterval(x[0:(N//2)],x0,?);

### low-level routines
def IsPowerOf2(n):
  if ( (n & (n-1) == 0) and (n != 0) ): return True
  else: return False


def romberg(x_length,fx,eps,n0):
  N = len(fx)
  if ( not IsPowerOf2(N-1) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( not IsPowerOf2(n0) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( n0 > (N-1) ): sys.exit("romberg: initial number of points is larger than provided grid")
  dx = x_length/(N-1)
  # Npow = np.rint(np.log2(N-1)) # should work up to 2**62, and above with rounding, be careful
  step = (N-1)//n0 # adjust to n0 points, divisibility already checked
  k1 = 0
  I = [] # empty list of lists to store the integrals
  while (step >= 1):
    I.append([])
    indices = [k2 for k2 in range(0,N,step)]
    for k2 in range(k1+1):
      if (k2 == 0): value = integrate.trapz(fx[indices],dx=step*dx)
      else: value = (4.0**k2 * I[k1][k2-1] - I[k1-1][k2-1]) / (4.0**k2-1.0)
      I[k1].append(value)

    if (k1>0):# convergence test
      Res = abs(I[k1][k1]-I[k1-1][k1-1])/abs(I[k1][k1])
      if (Res <= eps): return k1, I[k1][k1], I, Res

    step = step // 2
    k1 = k1+1

  return -1, I[-1][-1], I, Res # didn't converged in requested precision, returns the last value
  #  return [-1, I[-1][-1]] # didn't converged in requested precision, returns the last value

# xgrid = np.linspace(1.0,2.0,2049)
# fx = 1/(xgrid**2)
# nint, Int, full, err = romberg(1.0,fx,1e-15,4)

############# gaussian beam
def GaussianBeamRayleighRange(w0,lambd):
  return np.pi*w0**2/lambd


def invRadius(z,zR):
  return z/(zR**2+z**2)


def GaussianBeamCurvaturePhase(r,z,k0,zR):
  return 0.5*k0*invRadius(z,zR)*r**2


def waist(z,w0,zR):
  return w0*np.sqrt(1.0+(z/zR)**2);


def GaussianBeam(r,z,t,I0,w0,tFWHM,lambd):
  zR = np.pi*w0**2/lambd;
  w=w0*np.sqrt(1.0+(z/zR)**2);
  I=I0*((w0/w)**2)*np.exp(-2.0*(r/w)**2)*np.exp(-(2.0*np.sqrt(np.log(2.0))*t/tFWHM)**2);
  k0=2.0*np.pi/lambd;
  phase = GaussianBeamCurvaturePhase(r,z,k0,zR);
  return I, phase


# conversion of photons
def ConvertPhoton(x,inp,outp):
  # convert to omega in a.u.
  if (inp == 'omegaau'): omega = x
  elif (inp == 'lambdaSI'): omega = 2.0 * np.pi* units.hbar / (x * units.elmass * units.c_light * units.alpha_fine**2);
  elif (inp == 'lambdaau'): omega = 2.0 * np.pi/(units.alpha_fine*x);
  elif (inp == 'omegaSI'): omega = x * units.TIMEau;
  elif (inp == 'eV'): omega = x * np.elcharge/(units.elmass*units.alpha_fine**2*units.c_light**2);
  elif (inp == 'T0SI'): omega = units.TIMEau*2.0*np.pi/x;
  elif (inp == 'T0au'): omega = 2.0*np.pi/x;
  elif (inp == 'Joule'): omega = x / (units.elmass*units.alpha_fine**2 * units.c_light**2);
  else: sys.exit('Wrong input unit')

  # convert to output
  if (outp == 'omegaau'): return omega;
  elif (outp == 'lambdaSI'): return 2.0*np.pi*units.hbar/(omega*units.elmass*units.c_light*units.alpha_fine**2);
  elif (outp == 'lambdaau'): return 2.0*np.pi/(units.alpha_fine*omega);
  elif (outp == 'omegaSI'): return omega/units.TIMEau;
  elif (outp == 'eV'): return omega/(units.elcharge/(units.elmass*units.alpha_fine**2 * units.c_light**2));
  elif (outp == 'T0SI'): return units.TIMEau*2.0*np.pi/omega;
  elif (outp == 'T0au'): return 2.0*np.pi/omega;
  elif (outp == 'Joule'): return omega*(units.elmass*units.alpha_fine**2 * units.c_light**2);
  else: sys.exit('Wrong output unit')

## define dipole function
# def dipoleTimeDomainApp(z_medium,tgrid,r,I0,PhenomParams,tcoeff,rcoeff,LaserParams): # some global variables involved
# #  tcoeff = 4.0*np.log(2.0)*units.TIMEau**2 / ( TFWHMSI**2 )
# #  rcoeff = 2.0/(w0r**2)
#   omega0 = LaserParams['omega0']
#   kw0 = 2.0*np.pi/LaserParams['lambda']
#   phiIR = IRphase(r,z_medium,kw0,LaserParams['zR'])
#   res = []
#   NumHarm = PhenomParams.shape[1]
#   for k1 in range(len(tgrid)):
#     res1 = 0.0*1j;
#     intens = I0*np.exp(-tcoeff*(tgrid[k1])**2 - rcoeff*r**2)
#     for k2 in range(NumHarm): res1 = res1 + intens*np.exp(1j*(tgrid[k1]*omega0*PhenomParams[0,k2]-PhenomParams[1,k2]*intens + PhenomParams[0,k2]*phiIR))
#     res.append(res1); ## various points in time
#   return np.asarray(res)



# part to compute the field from the intensity list




## handling HDF5 files
def adddataset(h_path,dset_name,dset_data,unit):
  dset_id = h_path.create_dataset(dset_name,data=dset_data)
  dset_id.attrs['units']=np.string_(unit)


def addrealdataset_setprec(h_path, dset_name, dset_data, unit, precision):
  if ( precision == 'f'):
    dset_id = h_path.create_dataset(dset_name, data=dset_data, dtype='f')
    dset_id.attrs['units'] = np.string_(unit)
  elif (precision == 'd'):
    dset_id = h_path.create_dataset(dset_name, data=dset_data, dtype='d')
    dset_id.attrs['units'] = np.string_(unit)












