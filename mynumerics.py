from scipy import special
import numpy as np
import math
import sys
import units


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


############# gaussian beam
def invRadius(z,zR):
  return z/(zR**2+z**2)

def IRphase(r,z,k0,zR):
  return 0.5*k0*invRadius(z,zR)*r**2



## define dipole function
def dipoleTimeDomainApp(z_medium,tgrid,r,I0,PhenomParams,tcoeff,rcoeff,LaserParams): # some global variables involved
#  tcoeff = 4.0*np.log(2.0)*units.TIMEau**2 / ( TFWHMSI**2 )
#  rcoeff = 2.0/(w0r**2)
  omega0 = LaserParams['omega0']
  kw0 = 2.0*np.pi/LaserParams['lambda']
  phiIR = IRphase(r,z_medium,kw0,LaserParams['zR'])
  res = []
  NumHarm = PhenomParams.shape[1]
  for k1 in range(len(tgrid)):
    res1 = 0.0*1j;
    intens = I0*np.exp(-tcoeff*(tgrid[k1])**2 - rcoeff*r**2)
    for k2 in range(NumHarm): res1 = res1 + intens*np.exp(1j*(tgrid[k1]*omega0*PhenomParams[0,k2]-PhenomParams[1,k2]*intens + PhenomParams[0,k2]*phiIR)) 
    res.append(res1); ## various points in time
  return np.asarray(res)





















