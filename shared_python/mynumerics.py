from scipy import special
from scipy import integrate
from scipy import interpolate
import numpy as np
import math
import sys
import units
import h5py
import warnings


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


def FindInterval(x,x0): 
  '''
   find an index corresponding to given x0 value interval.
   ordering <  ), < ),..., < >; throws error otherwise
  '''
  if hasattr(x0, "__len__"):
    indices = []
    for x_loc in x0:
      indices.append(FindInterval(x,x_loc))
    return np.asarray(indices)

  else:    
    N = len(x)
    if ( (x0 > x[-1]) or (x0 < x[0]) ): raise LookupError('out of range in FindInterval')
    k1 = 1; k2 = N; length = N;
    while (length > 2):
      if ( x0 < x[k1 - 1 + length//2]): k2 = k1  + length//2
      else: k1 = k1 + length//2
      length = k2-k1

    if ((length == 2) and (x0 >= x[k1])): return k1
    else: return k1-1
    
  # for k1 in range(N-2):
  #   if ( (x[k1]<= x0) and (x0 < x[k1+1]) ): return k1
  # if ( (x[N-2]<= x0) and (x0 <= x[N-1]) ): return N-2
  # sys.exit('out of range in FindInterval')
  # if ( (x0 < x[0]) or (x0 > x[N-1])): sys.exit('out of range in FindInterval')
  # if (N == 2): return k1; # we finished looking
  # else:
  #   if (x0 > x[N//2] ): return FindInterval(x[(N//2):(N-1)],x0,N//2+?); # bookkeeping needed here... best will be additions and subtractions to be in-place
  #   else : return FindInterval(x[0:(N//2)],x0,?);

# a = []
# xxx = np.asarray([1.0, 2, 3, 4, 5, 6, 7, 8, 9])
# for k1 in range(len(xxx)-1): a.append(FindInterval(xxx,0.1+k1+1))

def find_index_across_arrays(vals,lists): ## inverse to the previous one
  if not(len(vals) == len(lists)): raise ValueError('mismatched lengths of values and lists')
  n_list = len(lists[0])
  n_vals = len(vals)
  candidates = []
  for k1 in range(n_list): # create candidates
    if (vals[0] == lists[0][k1]): candidates.append(k1)
  for k1 in range(1,n_vals): # remove candidates not matching conditions
    if not(len(lists[k1]) == n_list): raise ValueError('an array of different size from the first one present')
    old_candidates = candidates.copy()
    for k2 in range(len(old_candidates)):
        if not(vals[k1] == lists[k1][old_candidates[k2]]): del candidates[k2]
  
  if not(len(candidates) == 1): raise ValueError('there is no or multiple candidates')
  return candidates[0]



### low-level routines

def get_odd_interior_points(interval):
    res = np.asarray(
          list(range(int(np.floor(interval[0])),int(np.ceil(interval[1]))))
          )
    return res[(res%2==1)*(res>=interval[0])*(res<=interval[1])]

def get_divisible_interior_points(interval,divisor):
    res = np.asarray(
          list(range(int(np.floor(interval[0])),int(np.ceil(interval[1])+1)))
          )
    return res[(res%divisor==0)*(res>=interval[0])*(res<=interval[1])]

def IsPowerOf2(n):
  if ( (n & (n-1) == 0) and (n != 0) ): return True
  else: return False


def contains_substrings(string,substrings):
    for substring in substrings:
        if substring in string:
            return True
    return False


def tensor_constructor(arrays, fun,
                           init_args = [],
                           fixed_kwargs = {},
                           outputs_selector=None):
    """
    This constructor returns a recursive list of 'fun' evaluated for all
    combinations of inputs in 'arrays' inputted after 'init_args' i.e.
    result[i][j]... = fun(*init_args,array[0][i],array[1][j],...,**fixed_kwargs).    

    Parameters
    ----------
    arrays : list of arrays
        
    fun : function
        arguments passed to this function: fun(*init_args, 'args from arrays',
        **fixed_kwargs), the number of inputs and arrays has to match.
    init_args : args, optional
        The fixed arguments passed to 'fun'. The default is [].
    fixed_kwargs : TYPE, optional
        The default is {}.
    outputs_selector : TYPE, optional
        Enables to select only desired outputs from fun (fun(.)[outputs_selector]
        is returned). The default is None.

    Returns
    -------
    list
        recursive list containg fun at the final nodes (see the initial description).
        
    Example
    -------
    tensor_constructor([[1,2],[3,4,5]],fun) =
    [[1+3, 1+4, 1+5], [2+3, 3+4, 2+5]]
    for 'fun' being the sum of inputs

    """
    if len(arrays) == 1 :
        if (outputs_selector == None):
            return [fun(
                        *[*init_args, arrays[0][k1]],
                        **fixed_kwargs)
                    for k1 in range(len(arrays[0]))]
        else:
            return [fun(
                        *[*init_args, arrays[0][k1]],
                        **fixed_kwargs)[outputs_selector]
                    for k1 in range(len(arrays[0]))]
    else:
        return [tensor_constructor(arrays[1:], fun,
                init_args=init_args + [arrays[0][k1]],
                fixed_kwargs=fixed_kwargs,
                outputs_selector=outputs_selector)
                for k1 in range(len(arrays[0]))]  


## CALCULUS
def ddx_arb(k,x,fx):
  h2 = x[k+1]-x[k]
  h1 = x[k]-x[k-1]
  ratio = h2/h1
  return (fx[k+1] - fx[k-1]*ratio**2 - (1.0-ratio**2)*fx[k])/(h2*(1.0+ratio))

def ddx_vec_arb(x,fx):
  dfx = np.zeros(len(fx),dtype=fx.dtype)
  dfx[0] = (fx[1]-fx[0])/(x[1]-x[0])
  dfx[-1] = (fx[-1]-fx[-2])/(x[-1]-x[-2])
  for k1 in range(1,len(fx)-1):
    dfx[k1] = ddx_arb(k1,x,fx)
  return dfx

def complexify_fft(fx,convention='+'):
  N = len(fx)
  fx = np.fft.fft(fx)
  for k1 in range((N // 2) + 1, N):
    fx[k1] = 0.0
  if     (convention == '+'): fx = 2.0 * np.fft.ifft(fx)
  elif   (convention == '-'): fx = 2.0 * np.conj(np.fft.ifft(fx))
  else:  raise ValueError('Convention must be "+" or "-"') 
  return fx

def fft_t_nonorm(t, ft):
  Nt = len(t)
  Ft = np.conj(np.fft.fft(ft)[0:((Nt // 2) + 1)])
  omega = np.linspace(0, (np.pi * (Nt - 1) / (t[-1] - t[0])), (Nt // 2) + 1)
  return omega, Ft, Nt

def fft_t(t, ft):
  Nt = len(t)
  t0_ind = FindInterval(t,0.0)
  dt = t[t0_ind+1] - t[t0_ind]
  Ft = (dt/(np.sqrt(2.0*np.pi)))*np.conj(np.fft.fft(ft)[0:((Nt // 2) + 1)])
  omega = np.linspace(0, (np.pi * (Nt - 1) / (t[-1] - t[0])), (Nt // 2) + 1)
  return omega, Ft, Nt

def ifft_t_nonorm(omega, Ft, Nt):
  Ft = np.append(Ft, np.flip(np.conj(Ft[1:(Nt - len(Ft) + 1)])))
  ft = np.flip(np.fft.ifft(Ft))
  t = np.linspace(0, 2.0 * np.pi * (len(omega)) / (omega[-1] - omega[0]), len(ft))
  return t, ft

def integrate_subinterval(fx,x,xlim):
    if ( (xlim[0]<x[0]) or (xlim[-1]>x[-1]) ):
        raise ValueError('integration out of bounds')
    
    k_low = FindInterval(x, xlim[0]) + 1
    k_up = FindInterval(x, xlim[-1])
    
    x_low = x[k_low]; x_up = x[k_up] 
    fx_low = fx[k_low]; fx_up = fx[k_up] 
    
    delta_x = x[k_low] - x[k_low-1]
    fx_min = (fx[k_low-1] * (xlim[0] - x[k_low-1])/delta_x + fx[k_low] * (x[k_low]-xlim[0])/delta_x )
    delta_x = x[k_up] - x[k_up-1]
    fx_max = (fx[k_up] * (xlim[-1] - x[k_up])/delta_x + fx[k_up+1] * (x[k_up+1]-xlim[-1])/delta_x )
    
    # integrate in the intervale and add trpaezoidal parts of outer intervals
    Ifx = np.trapz(fx[k_low:(k_up+1)],x[k_low:(k_up+1)]) + \
          0.5 * (fx_low + fx_min) * (x_low-xlim[0]) + \
          0.5 * (fx_up + fx_max) * (xlim[-1]-x_up)    
    
    return Ifx

# fx = list(range(10))
# x = fx
# a = integrate_subinterval(fx,x,[1.5,5.5])

def romberg(x_length,fx,eps,n0):
  N = len(fx)
  if ( not IsPowerOf2(N-1) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( not IsPowerOf2(n0) ): sys.exit("romberg: initial stepsize isn't 2**k")
  elif ( n0 > (N-1) ): sys.exit("romberg: initial number of points is larger than provided grid")
  dx = x_length/(N-1)
  step = (N-1)//n0 # adjust to n0 points, divisibility already checked
  k1 = 0
  I = [] # empty list of lists to store the integrals
  while (step >= 1):
    I.append([])
    indices = [k2 for k2 in range(0,N,step)]
    for k2 in range(k1+1):
      if (k2 == 0): value = integrate.trapz(fx[indices],dx=step*dx) # this is inefficient, we already keep the previous results, just refine them
      else: value = (4.0**k2 * I[k1][k2-1] - I[k1-1][k2-1]) / (4.0**k2-1.0)
      I[k1].append(value)

    if (k1>0):# convergence test
      Res = abs(I[k1][k1]-I[k1-1][k1-1])/abs(I[k1][k1])
      if (Res <= eps): return k1, I[k1][k1], Res

    step = step // 2
    k1 = k1+1

  return -1, I[-1][-1], Res # didn't converged in requested precision, returns the last value
  #  return [-1, I[-1][-1]] # didn't converged in requested precision, returns the last value


def romberg(x_length,fx,eps,n0):
  N = len(fx)
  if ( not IsPowerOf2(N-1) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( not IsPowerOf2(n0) ): sys.exit("romberg: initial stepsize isn't 2**k")
  elif ( n0 > (N-1) ): sys.exit("romberg: initial number of points is larger than provided grid")
  dx = x_length/(N-1)
  step = (N-1)//n0 # adjust to n0 points, divisibility already checked
  k1 = 0
  I = [] # empty list of lists to store the integrals
  while (step >= 1):
    I.append([])
    indices = [k2 for k2 in range(0,N,step)]
    for k2 in range(k1+1):
      if (k2 == 0): value = integrate.trapz(fx[indices],dx=step*dx) # this is inefficient, we already keep the previous results, just refine them
      else: value = (4.0**k2 * I[k1][k2-1] - I[k1-1][k2-1]) / (4.0**k2-1.0)
      I[k1].append(value)

    if (k1>0):# convergence test
      Res = abs(I[k1][k1]-I[k1-1][k1-1])/abs(I[k1][k1])
      if (Res <= eps): return k1, I[k1][k1], Res

    step = step // 2
    k1 = k1+1

  return -1, I[-1][-1], Res # didn't converged in requested precision, returns the last value
  #  return [-1, I[-1][-1]] # didn't converged in requested precision, returns the last value

# xgrid = np.linspace(1.0,2.0,2049)
# fx = 1/(xgrid**2)
# nint, Int, full, err = romberg(1.0,fx,1e-15,4)


## Working with field
# conversion of photons
def ConvertPhoton(x,inp,outp):
  """
  available I/O: 'omegaau', 'omegaSI', 'lambdaau', 'lambdaSI', 'eV', 'T0au',
                 'T0SI', 'Joule'     
  """
  # convert to omega in a.u.
  if (inp == 'omegaau'): omega = x
  elif (inp == 'lambdaSI'): omega = 2.0 * np.pi* units.hbar / (x * units.elmass * units.c_light * units.alpha_fine**2);
  elif (inp == 'lambdaau'): omega = 2.0 * np.pi/(units.alpha_fine*x);
  elif (inp == 'omegaSI'): omega = x * units.TIMEau;
  elif (inp == 'eV'): omega = x * units.elcharge/(units.elmass*units.alpha_fine**2*units.c_light**2);
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
  
  
def FieldToIntensitySI(Efield):
    return 0.5*units.c_light*units.eps0*Efield**2

def IntensityToFieldSI(Intensity):
    return np.sqrt((2.0/(units.c_light*units.eps0))*Intensity)

def Spectrum_lambda2omega(lambdagrid,Spectrum_lambda, include_Jacobian = True):
    ogrid = 2.*np.pi*units.c_light/lambdagrid
    if include_Jacobian: Spectrum_omega = ((ogrid**2)*Spectrum_lambda)/(2.*np.pi*units.c_light)
    else: Spectrum_omega = Spectrum_lambda    
    return ogrid, Spectrum_omega

## Gaussian beam
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

def GaussianBeamEfield(r,z,t,E0,w0,tFWHM,lambd, comoving = True):
  # Gaussian beam in the comoving frame with c
  if (not(comoving)): t = t - z/units.c_light
  omega0 = ConvertPhoton(lambd, 'lambdaSI', 'omegaSI')
  k0 = 2.0*np.pi/lambd
  zR = np.pi*w0**2/lambd
  w=w0*np.sqrt(1.0+(z/zR)**2)
  phase_Gouy = np.arctan(z/zR)
  phase_curv = GaussianBeamCurvaturePhase(r,z,k0,zR)
  return E0*(w0/w)*np.exp(-(r/w)**2)*np.exp(-(2.0*np.log(2.0)*t/tFWHM)**2)*np.cos(omega0*t -phase_curv + phase_Gouy)
 
 

# def GaussianBeamEfieldConstruct(r,z,t,E0,w0,tFWHM,lambd): # (t,r,z)
#   Nz = len(z); Nr = len(r); Nt = len(t) 
#   Efield = np.zeros((Nt,Nr,Nz))
#   for k1 in range(Nt):
#     for k2 in range(Nr):
#       for k3 in range(Nz):
#         Efield[k1,k2,k3] = GaussianBeamEfield(r[k2],z[k3],t[k1],E0,w0,tFWHM,lambd)
#   # retarded time
#   return Efield

def GaussianBeamEfieldConstruct(rgrid,zgrid,tgrid,E0,w0,tFWHM,lambd, comoving = True):
  r, z, t = np.meshgrid(rgrid,zgrid,tgrid, indexing='ij')
  if (not(comoving)): t = t - z/units.c_light
  omega0 = ConvertPhoton(lambd, 'lambdaSI', 'omegaSI')
  k0 = 2.0*np.pi/lambd
  zR = np.pi*w0**2/lambd
  w=w0*np.sqrt(1.0+(z/zR)**2)
  phase_Gouy = np.arctan(z/zR)
  invRadius = z/(zR**2+z**2)
  phase_curv = 0.5*k0*invRadius*r**2
  return E0*(w0/w)*np.exp(-(r/w)**2 - (2.0*np.log(2.0)*t/tFWHM)**2)*np.cos(omega0*t -phase_curv + phase_Gouy)


def Gaussian_phase_map(z,r,w0,lambd,n=1.0,vacuum_frame=True,
                       incl_curv = True, incl_Gouy = True,
                       incl_lin = True):
    k = 2.0*np.pi*n/(lambd)
    zR = np.pi*(w0**2)*n/lambd
    inv_curv_radius = z/(zR**2+z**2)
    phi_G = np.arctan(z/zR)
    if vacuum_frame: k_corr = 2.0*np.pi*(n-1.0)/(lambd)
    else: k_corr = k
        
    phase = 0.
    if incl_lin: phase += k_corr*z
    if incl_curv: phase += 0.5*k*(r**2)*inv_curv_radius 
    if incl_Gouy: phase += -phi_G
    
    return phase


def Gaussian_E0_map(z,r,w0,E0,lambd,n=1.0,
                        incl_z_profile = True,
                        incl_radial_wz_profile = True,
                        create_mesh = True):
    
    # raise wrning if potentially wrong shapes are passed
    if ( (np.shape(z) != np.shape(r))
                    and
         ((np.shape(r) != (1,)) and (np.shape(r) != ()))
                    and
         ((np.shape(z) != (1,)) and (np.shape(z) != ()))
        ):
        warnings.warn("The 'z'- and 'r'-shapes doesn't match. It might produce undesired result. Consider using 'np.meshgrid'.")
    
    zR = np.pi*(w0**2)*n/lambd
    w_z = lambda z: w0*np.sqrt(1+(z/zR)**2)
        
    E0_rz = E0
    if incl_z_profile: E0_rz *= (w0/w_z(z))
    
    if incl_radial_wz_profile: E0_rz *= np.exp(-(r/w_z(z))**2)
    else: E0_rz *= np.exp(-(r/w0)**2)
    
    # resize to get z-dimension if no modifiers are applied
    if (not(incl_z_profile) and not(incl_radial_wz_profile) and (np.shape(r)==(1,))):
        E0_rz = E0_rz[0]*np.ones(np.shape(z))
    if (not(incl_z_profile) and not(incl_radial_wz_profile) and (np.shape(r)==())):
        E0_rz = E0_rz*np.ones(np.shape(z))
    
    return E0_rz

class pulse_types:
    def __init__(self, pulse_type):
        if (pulse_type == 'sin2'): # sin^2 - envelope
            def sin2pulse(t,omega0,omegac,E0,phi0):
                phi_central = - (np.pi * omega0) / (2.0*omegac) # use the peak of the pulse as the reference for the cosine pulse
                return (t>=0) * (t<=np.pi/omegac) * E0*((np.sin(omegac*t))**2) *  np.cos(omega0*t + phi_central + phi0)
            self.pulse = sin2pulse
            # self.inputs_list_direct = ['omega0', 'omegac', 'E0', 'phi0']
            # self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
            def inputs_converter1(lambdaSI,tFWHMSI,E0,phi0):
                omega0 = ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
                omegac = 2.0*np.arccos(2**(-0.25))/(tFWHMSI/units.TIMEau)
                return omega0, omegac, E0, phi0
            self.inputs_converter1 = inputs_converter1
            
            # def construct_tgrid(Nt, **kwargs):
            #     omegac = 2.0*np.arccos(2**(-0.25))/(kwargs['tFWHMSI']/units.TIMEau)
            #     tgrid = np.linspace(0,np.pi/omegac,Nt)
            #     return tgrid
            # self.construct_tgrid = construct_tgrid
            
            def make_field(Nt,**kwargs):
                tgrid = construct_tgrid(Nt, **kwargs)
                Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
                return tgrid, Efield
            self.make_field = make_field
            
            def inputs_converter(*args,given_inps=['omega0','omegac','E0','phi0']):
                
                # omega0
                if ('omega0' in given_inps): omega0 = args[given_inps.index('omega0')]
                elif ('lambda0' in given_inps): omega0 = ConvertPhoton(args[given_inps.index('lambda0')], 'lambdaau', 'omegaau') 
                else: raise NotImplementedError("omega0 doesn't have input")
                
                # omegac
                if ('omegac' in given_inps): omegac = args[given_inps.index('omegac')]
                elif ('T_FWHM' in given_inps):
                    omegac = 2.0*np.arccos(2.0**(-0.25))/args[given_inps.index('T_FWHM')]
                elif ('Ncyc' in given_inps):
                    omegac = omega0/(2.0*args[given_inps.index('Ncyc')])
                else: raise NotImplementedError("omegac doesn't have input")
                
                # E0
                if ('E0' in given_inps): E0 = args[given_inps.index('E0')]
                else: raise NotImplementedError("E0 doesn't have input")
                
                # phi0
                if ('phi0' in given_inps): phi0 = args[given_inps.index('phi0')]
                else: raise NotImplementedError("phi0 doesn't have input")
                
                return omega0, omegac, E0, phi0
            self.inputs_converter = inputs_converter
            
            def construct_tgrid(*args, N_points_control='dt', duration_definition='T_FWHM'):
                
                if (duration_definition=='T_FWHM'):
                    omegac = 2.0*np.arccos(2.0**(-0.25))/args[1]
                elif (duration_definition=='omegac'):
                    omegac = args[1]
                elif (duration_definition=='Ncyc'):
                    omegac = args[1]/(2.0*args[2])
                else: raise NotImplementedError("srongly sepcified 'duration_definition'")
                
                if (N_points_control=='dt'):
                    tgrid = np.arange(0, (np.pi/omegac)+args[0], args[0])
                elif (N_points_control=='Nt'):
                    tgrid = np.linspace(0,np.pi/omegac,args[0])
                else: raise NotImplementedError("srongly sepcified 'N_points_control'")
                return tgrid
            self.construct_tgrid = construct_tgrid
            

        # elif (pulse_type == 'Gaussian'): # sin^2 - envelope
        #     def Gaussian_pulse(t,omega0,tFWHM,E0,phi0):
        #         return E0* np.exp(-(2.0*np.log(2.0)*t/tFWHM)**2)*  np.cos(omega0*t + phi0)
        #     self.pulse = Gaussian_pulse
        #     self.inputs_list_direct = ['omega0', 'T_FWHM', 'E0', 'phi0']
        #     self.inputs_list = ['lambda', 'T_FWHM', 'E0', 'phi0']
            
        #     def inputs_converter(lambdaSI,tFWHMSI,E0,phi0):
        #         omega0 = mn.ConvertPhoton(lambdaSI, 'lambdaSI', 'omegaau')
        #         return omega0, tFWHMSI/units.TIMEau, E0, phi0
        #     self.inputs_converter = inputs_converter
            
        #     def construct_tgrid(Nt, **kwargs):
        #         tmax = kwargs['t_expand'] * kwargs['tFWHMSI']/units.TIMEau
        #         tgrid = np.linspace(-0.5*tmax, 0.5*tmax,Nt)
        #         return tgrid
        #     self.construct_tgrid = construct_tgrid
            
        #     def make_field(Nt,**kwargs):
        #         tgrid = construct_tgrid(Nt, **kwargs)
        #         Efield = sin2pulse(tgrid,*inputs_converter(**kwargs))
        #         return tgrid, Efield
        #     self.make_field = make_field
                
        else: raise NotImplementedError('only sin2 implementented now')
        
        
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
  dset_id.attrs['units']=np.bytes_(unit)
  return dset_id


def addrealdataset_setprec(h_path, dset_name, dset_data, unit, precision):
  if ( precision == 'f'):
    dset_id = h_path.create_dataset(dset_name, data=dset_data, dtype='f')
    dset_id.attrs['units'] = np.bytes_(unit)
  elif (precision == 'd'):
    dset_id = h_path.create_dataset(dset_name, data=dset_data, dtype='d')
    dset_id.attrs['units'] = np.bytes_(unit)


def readscalardataset(file,path,type): # type is (S)tring or (N)umber
  if (type == 'N'): return file[path][()]
  elif (type == 'S'):
    try:
      return file[path][()].decode()
    except:
      try:
        return file[path][0].decode()
      except:
        raise TypeError('problem in decoding string')
    else:
      raise TypeError('problem in decoding string')
      
  else: raise TypeError('accepts only (N)umbers or (S)trings')


def h5_seek_for_scalar(file,dtype,*args):
    for path in args:
        try:
            return readscalardataset(file,path,dtype)
        except:
            pass
    raise ReferenceError('Dateset not found in args.')


## Other functions

def romberg_test(x_length,fx,eps,n0):
  N = len(fx)
  if ( not IsPowerOf2(N-1) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( not IsPowerOf2(n0) ): sys.exit("romberg: initial stepsize isn't 2**k")
  elif ( n0 > (N-1) ): sys.exit("romberg: initial number of points is larger than provided grid")
  dx = x_length/(N-1)
  step = (N-1)//n0 # adjust to n0 points, divisibility already checked
  k1 = 0
  I = [] # empty list of lists to store the integrals
  while (step >= 1):
    I.append([])
    indices = [k2 for k2 in range(0,N,step)]
    for k2 in range(k1+1):
      if (k2 == 0): value = integrate.trapz(fx[indices],dx=step*dx) # this is inefficient, we already keep the previous results, just refine them
      else: value = (4.0**k2 * I[k1][k2-1] - I[k1-1][k2-1]) / (4.0**k2-1.0)
      I[k1].append(value)

    if (k1>0):# convergence test
      Res = abs(I[k1][k1]-I[k1-1][k1-1])/abs(I[k1][k1])
      if (Res <= eps): return k1, I[k1][k1], Res, I

    step = step // 2
    k1 = k1+1

  return -1, I[-1][-1], Res, I # didn't converged in requested precision, returns the last value
  #  return [-1, I[-1][-1]] # didn't converged in requested precision, returns the last value


def rombergeff_test(x_length,fx,eps,n0):
  N = len(fx)
  if ( not IsPowerOf2(N-1) ): sys.exit("romberg: input isn't 2**k+1")
  elif ( not IsPowerOf2(n0) ): sys.exit("romberg: initial stepsize isn't 2**k")
  elif ( n0 > (N-1) ): sys.exit("romberg: initial number of points is larger than provided grid")
  dx = x_length/(N-1)
  step = (N-1)//n0 # adjust to n0 points, divisibility already checked
  k1 = 0
  I = [] # empty list of lists to store the integrals
  while (step >= 1):
    I.append([])
    indices = [k2 for k2 in range(0,N,step)]
    for k2 in range(k1+1):
      if (k2 == 0): value = integrate.trapz(fx[indices],dx=step*dx) # this is inefficient, we already keep the previous results, just refine them
      else: value = (4.0**k2 * I[k1][k2-1] - I[k1-1][k2-1]) / (4.0**k2-1.0)
      I[k1].append(value)

    if (k1>0):# convergence test
      Res = abs(I[k1][k1]-I[k1-1][k1-1])/abs(I[k1][k1])
      if (Res <= eps): return k1, I[k1][k1], Res, I

    step = step // 2
    k1 = k1+1

  return -1, I[-1][-1], Res, I # didn't converged in requested precision, returns the last value
  #  return [-1, I[-1][-1]] # didn't converged in requested precision, returns the last value


def gabor_transf(arr, t, t_min, t_max, N_steps, a, omegamax = -1.0):
  '''
  Gabor transform
  ===============


  gabor_transf(arr, t, t_min, t_max, N_steps = 400, a = 8)


  Computes Gabor transform for arbitrary number array 'arr' of 
  lenght N. 

      Parameters:
          arr (np.array or list of floats): input data for Gabor transform, length N
          t (np.array or list of floats): time domain for data, length N
          t_min (float): minimum time for Gabor transform domain
          t_max (float): maximum time for Gabor transform domain

      Optional parameters:
          N_steps (int): number of discretization points
          a (float): Gabor window parameter
  
      Returns:
          gabor_transf (np.array(N_steps, N)): Gabor transform of arr 

  Note: time domain 't' must correspond to the array 'arr'.

  Example:
      import numpy as np
      import random

      t_min = 0
      t_max = 1
      N_steps = 100

      x = np.linspace(0,1,100)
      y = [np.cos(2*np.pi*t) + np.sin(np.pi*t) + 0.1*random.randrange(-1,1) for t in x]

      z = gabor_transform(y, x, t_min, t_max, N_steps = N_steps)

  '''
  if (len(arr) != len(t)):
    raise ValueError('Arrays must have same dimension.')

  if (t_max < t_min):
    raise ValueError('Maximum time must be larger than minimum time')

  ### Time domain for Gabor transform
  t_0 = np.linspace(t_min, t_max, N_steps)

  ### Number of sample points for fft
  N = len(arr)

  ### Init np.array
  gabor_transf = np.zeros((N_steps, (N // 2) + 1))
  
  omega = np.linspace(0, (np.pi * (N - 1) / (t[-1] - t[0])), (N // 2) + 1)

  ### Compute Gabor transform (normalised)
  for i in range(0, N_steps):
    fft_loc = np.fft.fft(np.exp(-np.power((t-t_0[i])/a, 2))*arr[:])
    gabor_transf[i, :] = 2.0 / N * np.abs(fft_loc[0:((N // 2) + 1)])
    
  if (omegamax < 0.0):    
      return t_0, omega, gabor_transf
  else:
      k_omegamax = FindInterval(omega, omegamax)
      return t_0, omega[:k_omegamax], gabor_transf[:,:k_omegamax]


## SIGNAL PROCESSING
def identity(M):
    return np.ones(M)

def filter_box(M,xgrid,x_box,apply_blackman = False):
    if apply_blackman:
        k1 = FindInterval(xgrid, x_box[0])
        k2 = FindInterval(xgrid, x_box[1])
        # np.blackman(k2-k1)
        return np.concatenate((np.zeros(k1), np.blackman(k2-k1), np.zeros(M-k2)))
    else:
        return (x_box[0] <= xgrid)*(xgrid <= x_box[1])


def apply_filter(fx,filter_funct,*args,**kwargs):
    if (filter_funct == identity): return fx
    else:   return filter_funct(len(fx),*args,**kwargs) * fx


def clamp_array(array,filter_type,filter_threshold, replace_value = 0.0):
    if (filter_type is None): return array
    array_flat = array.flatten()
    if (filter_type == 'low_abs_cut'):
        dum = filter_threshold*np.max(np.abs(array_flat))
        if (replace_value == 'threshold'): replace_value = dum
        array_flat[np.abs(array_flat) < dum] = replace_value
        return array_flat.reshape(array.shape)
    elif ('low_abs2_cut'):
        dum = filter_threshold*np.max(np.abs(array_flat))**2
        if (replace_value == 'threshold'): replace_value = np.sqrt(dum)
        array_flat[np.abs(array_flat)**2 <= dum] = replace_value #np.sqrt(dum) # 1e-10
        return array_flat.reshape(array.shape)
    else:
        warnings.warn('filter wrongly specified: returned original')
        return array
    
## BEAM MEASURE
def measure_beam_RMS(x,fx):
    return np.sqrt(
        integrate.trapz( (x**2) * fx, x=x) / integrate.trapz(fx, x=x)
        )


def measure_beam_max_ratio_zeromax(x,fx,alpha):
    N = len(fx); fxmax = fx[0]; fxFWHM = alpha*fxmax
    
    for k1 in range(N):
        if (fxFWHM > fx[k1]): break
    
    if (k1 == (N-1)): return np.Inf
    
    x1 = x[k1-1]; x2 = x[k1]
    fx1 = fx[k1-1]; fx2 = fx[k1]
    dx = x2 - x1; dfx = fx2 - fx1

    if (abs(dfx) < np.spacing(fxmax)): return x1 # compare with an appropriate machine epsilon for flat cases
    else:   return (dx*fxFWHM + fx2*x1 - fx1*x2) / dfx
    
    
def measure_beam_FWHM_zeromax(x,fx):
    return measure_beam_max_ratio_zeromax(x,fx,0.5)
    

def measure_beam_E_alpha_zeromax(x,fx,alpha): #  taken from my Matlab implementation
    E_tot = integrate.trapz(fx, x=x)
    E_prim_norm = (1.0/E_tot) * integrate.cumtrapz(fx, x=x, initial=0)
    
    
    # k1 = FindInterval(E_prim_norm, alpha) # It'd be faster, need to deal with the boundaries.
    
    N = len(fx)
    for k1 in range(N):
        if (E_prim_norm[k1] > alpha): break
    
    if (k1 == (N-1)): return np.Inf
    elif (k1 == 0): return x[0] 

    x1 = x[k1-1]; x2 = x[k1]
    fx1 = E_prim_norm[k1-1]; fx2 = E_prim_norm[k1]
    dx = x2 - x1; dfx = fx2 - fx1

    if (abs(dfx) < np.spacing(1.0)): return x1 # compare with an appropriate machine epsilon for flat cases
    else:   return (dx*alpha + fx2*x1 - fx1*x2) / dfx
    
    
    # cumtrapz(x, initial=0)
    
    
    
# interpolation etc.
def interpolate_2D(x,y,fxy,Nx,Ny):
    f_interp = interpolate.interp2d(x,y,fxy)
    xnew = np.linspace(x[0], x[-1], Nx)
    ynew = np.linspace(y[0], y[-1], Ny)
    fxy_interp = f_interp(xnew, ynew)
    return xnew, ynew, fxy_interp




## inputs abstractions
def multiparameters_lines2dict(lines):
    processing_fixed_inputs = False;  processing_varying_inputs = False
    varying_params = []; fixed_params = []; params_dict = {}
    params_dict['_units']={}
    for line in lines:
        sep_line = line.split()  # separate the line
        
        if ((len(sep_line) == 0) or (sep_line[0] == '#') or (sep_line[0] == '##')):
            pass  # print('empty or commented line')
        else:
          if (sep_line[0] == '$fixed'):
              processing_fixed_inputs = True; processing_varying_inputs = False
          elif (sep_line[0] == '$varying'):
              processing_fixed_inputs = False; processing_varying_inputs = True
          elif processing_fixed_inputs: 
              fixed_params.append(sep_line[0])
              if (sep_line[2] == 'R'): params_dict[sep_line[0]] = float(sep_line[1])
              elif (sep_line[2] == 'I'): params_dict[sep_line[0]] = int(sep_line[1])
              elif (sep_line[2] == 'S'): params_dict[sep_line[0]] = sep_line[1]
              else: raise TypeError('line:' + line + '\n Specify datatype by R or I.')
              params_dict['_units'][sep_line[0]] = sep_line[3]
          elif processing_varying_inputs: 
              varying_params.append(sep_line[0]) 
              try:
                  params_dict[sep_line[0]] = [float(sep_line[2]), float(sep_line[3]),
                                          int(sep_line[4])]
              except:
                  params_dict[sep_line[0]] = [float(sep_line[2]), float(sep_line[3]),
                                          int(float(sep_line[4]))]
              params_dict['_units'][sep_line[0]] = sep_line[1]
          else: raise NotImplementedError('The input file must follow the $fixed - $varying structure now')
          
    return varying_params, fixed_params, params_dict


class parameters_selector:
    def __init__(self, 
                 varying_params, fixed_params, params_constructor,
                 dtypes={}, assumed_output_order=None):
        self.varying_params = varying_params
        self.fixed_params = fixed_params
        param_grids = {}
        self.N_combinations = 1; self.N_fixed = len(fixed_params); self.N_varying = len(varying_params)
        self.varying_params_lengths = []
        
        # param = varying_list[k1]
        # N_points = parameters[param][2] + 2
        # if (parameters[param][2] == -1): param_grids.append(
        #                                             np.asarray([parameters[param][1]])
        #                                             )
        # else: param_grids.append( np.linspace( *parameters[param][:2], N_points ) )
        # N *= N_points; param_dims.append(N_points)  
        
        for param in varying_params:   
            dtype = None if not(param in dtypes.keys()) else dtypes[param]
            if (params_constructor[param][2] == -1):                
                param_grids[param] = np.asarray([params_constructor[param][1]],dtype=dtype)
            else:
                N_points = params_constructor[param][2] + 2
                param_grids[param] = np.linspace( *params_constructor[param][:2], N_points, dtype=dtype)
                self.N_combinations *= N_points
            self.varying_params_lengths.append( len(param_grids[param]) )
        for param in fixed_params:
            param_grids[param] = params_constructor[param]
        self.param_grids = param_grids
                
                    
                

        self.assumed_output_order = (list(varying_params+fixed_params)
                                     if (assumed_output_order is None) 
                                     else list(assumed_output_order))
        
        if ('_units' in params_constructor.keys()): self.units = params_constructor['_units']
        
        
    def ret(self,N,variables=None):
        MultInd = np.unravel_index(N, self.varying_params_lengths)
        output_ordered = {}; ouput_required = []
        for k1 in range(self.N_varying):
            output_ordered[self.varying_params[k1]] = self.param_grids[self.varying_params[k1]][MultInd[k1]]
        if (variables == None): variables=self.assumed_output_order
        all_possible_ouputs = output_ordered|{key:self.param_grids[key] for key in self.fixed_params}
        for var in variables:
            ouput_required.append(all_possible_ouputs[var])
        return ouput_required
    
    def store_to_h5(self,h_path):
        for k1 in range(self.N_varying):
            try:
                adddataset(h_path, 'param_'+str(k1) , self.param_grids[self.varying_params[k1]] , '['+self.units[self.varying_params[k1]] +']' )
            except:
                adddataset(h_path, 'param_'+str(k1) , self.param_grids[self.varying_params[k1]] , '[?]' )
        h_path.create_dataset('varying_params',data=np.bytes_(self.varying_params))
        
        for k1 in range(self.N_fixed):
            try:
                adddataset(h_path, self.fixed_params[k1] , self.param_grids[self.fixed_params[k1]] , '['+self.units[self.fixed_params[k1]] +']' )
            except:
                adddataset(h_path, self.fixed_params[k1] , self.param_grids[self.fixed_params[k1]] , '[?]' )
        h_path.create_dataset('fixed_params',data=np.bytes_(self.fixed_params))


# array processing
def symmetrize_y(y,array_xy):
    
    """
    Designed to make symmetric plots of radially symmetric data. It flips 
    y-axis and concenate a duplicated part of `array_xy` along this axis.
    It removes the first row (assuming 0 included).
    
    Parameters
    ----------
    y : array-like
        y-axis, (expected to be the radial axis)
    array_xy : the array to be symmetrized
        DESCRIPTION.

    Returns
    -------
    y_sym   : the symmerrized y-axis
    arr_sym : the symmetrized array 

    """
    
    y_sym   = np.concatenate((np.flip(-y[1:]), y))
    arr_sym = np.concatenate((np.flip(array_xy[:,1:],axis=1),
                                      array_xy[:,: ]           
                             ), axis = 1)
    return y_sym, arr_sym
    
    
    
# def si_prefix(number, base_unit):
#     prefixes = {
#         -24: "y",  # yocto
#         -21: "z",  # zepto
#         -18: "a",  # atto
#         -15: "f",  # femto
#         -12: "p",  # pico
#         -9: "n",  # nano
#         -6: "\mu",  # micro
#         -3: "m",  # milli
#         0: "",    # no prefix
#         3: "k",   # kilo
#         6: "M",   # mega
#         9: "G",   # giga
#         12: "T",  # tera
#         15: "P",  # peta
#         18: "E",  # exa
#         21: "Z",  # zetta
#         24: "Y"   # yotta
#     }


#     sign = -1 if number < 0 else 1
#     number *= sign
#     exp = 0

#     while number >= 1000 and exp < 24:
#         number /= 1000.0
#         exp += 3

#     while number < 1 and exp > -24:
#         number *= 1000.0
#         exp -= 3

#     prefix = prefixes.get(exp, "?")
#     return rf"{sign * number:.2f}{prefix}"