from scipy import special
from scipy import integrate
import numpy as np
import struct
import array
import os
import time
# import ray
# import matlab.engine
# import string
import multiprocessing as mp
import math
# import joblib
# from mpi4py import MPI
# import oct2py
import shutil
import h5py
import sys
import units
import mynumerics as mn

def ComputeFieldsInRFromIntensityList(z_medium, rgrid, Hgrid, Nomega, LaserParams, Igrid, FSourceterm):
    print('Computing fields in the grid');
    Nz_medium = len(z_medium); Nr = len(rgrid);
    FField_r=np.empty([Nz_medium,Nomega,Nr], dtype=np.cdouble)
    for k3 in range(Nz_medium):
      for k1 in range(Nr): # We use linear interpolation using the intensity-grid at the instant
        I_r, phase_r = mn.GaussianBeam(rgrid[k1],z_medium[k3],0,LaserParams['I0']/units.INTENSITYau,LaserParams['w0'],1,LaserParams['lambda']) # a.u.
        phase_XUV = phase_r*Hgrid

        # find a proper interval in the Igrid, we use linear interp written by hand now
        k2 = mn.FindInterval(Igrid,I_r)
        weight1 = (Igrid[k2+1]-I_r)/(Igrid[k2+1]-Igrid[k2]); weight2 = (I_r-Igrid[k2])/(Igrid[k2+1]-Igrid[k2]);
        FField_r[k3, :, k1] = units.EFIELDau*np.exp(1j*phase_XUV)*(weight1*FSourceterm[k2,:]+weight2*FSourceterm[k2,:]); # free-form works? # not sure about the sign, if really (-i*(omega*t-kz-phiIR)), than should be "exp(+1j*...)"
    return FField_r


# define function to integrate, there are some global variables! ## THE OUTPUT IS IN THE MIX OF ATOMIC UNITS (field) and SI UNITS (radial coordinate + dr in the integral)
def FieldOnScreen(z_medium, omegagrid, omega_step, rgrid, FField_r, rgrid_anal, zgrid_anal, k_start, k_num, integrator):
# this function computes the Hankel transform of a given source term in omega-domain stored in FField_r
# all the grids are specified in the inputs, except the analysis in omega_anal, this is specified by 'k_start' and 'k_num', it is used in the multiprocessing scheme
    Nz_anal = np.asarray(zgrid_anal.shape); Nz_anal = Nz_anal[1];
    Nr_anal = len(rgrid_anal); Nr = len(rgrid); Nz_medium=len(z_medium);
    FHHGOnScreen = np.empty([Nz_medium, Nz_anal, k_num, Nr_anal], dtype=np.cdouble)
    for k7 in range(Nz_medium): # loop over different medium positions
        print('process starting at omega', k_start, ' started computation of zgrid', z_medium[k7])
        k4 = 0  # # of loops in omega
        for k1 in range(k_num):  # Nomega
            k5 = k_start + k1 * omega_step  # accesing the grid
            tic = time.process_time()
            for k6 in range(Nz_anal):
                for k2 in range(Nr_anal):  # Nr
                    k_omega = omegagrid[k5] / (units.TIMEau * units.c_light);  # omega divided by time: a.u. -> SI
                    integrand = np.empty([Nr], dtype=np.cdouble)
                    for k3 in range(Nr): integrand[k3] = np.exp(-1j*k_omega * (rgrid[k3] ** 2) / (2.0 * (zgrid_anal[k7,k6] - z_medium[k7]))) * rgrid[k3] * FField_r[k7, k5, k3] * special.jn(0, k_omega * rgrid[k3] * rgrid_anal[k2] / (zgrid_anal[k7,k6] - z_medium[k7]));  # rescale r to atomic units for spectrum in atomic units! (only scaling)
                    # for k3 in range(Nr): integrand[k3] = rgrid[k3] * FField_r[k7, k5, k3] * special.jn(0, k_omega * rgrid[k3] * rgrid_anal[k2] / (zgrid_anal[k7,k6] - z_medium[k7]));  # rescale r to atomic units for spectrum in atomic units! (only scaling)
                    if (integrator['method'] == 'Romberg'):
                        nint, value, err = mn.romberg(rgrid[-1]-rgrid[0],integrand,integrator['tol'],integrator['n0'])
                        # FHHGOnScreen[k7, k6, k4, k2] = (1.0 / (zgrid_anal[k7, k6] - z_medium[k7])) * value;
                        FHHGOnScreen[k7, k6, k4, k2] = value;
                    # elif (integrator['method'] == 'Trapezoidal'): FHHGOnScreen[k7,k6, k4, k2] = (1.0 / (zgrid_anal[k7,k6] - z_medium[k7])) * integrate.trapz(integrand, rgrid);
                    elif (integrator['method'] == 'Trapezoidal'): FHHGOnScreen[k7, k6, k4, k2] = integrate.trapz(integrand, rgrid);
                    # elif (integrator['method'] == 'Simpson'): FHHGOnScreen[k7,k6, k4, k2] = (1.0 / (zgrid_anal[k7,k6] - z_medium[k7])) * integrate.simps(integrand, rgrid);
                    elif (integrator['method'] == 'Simpson'): FHHGOnScreen[k7, k6, k4, k2] = integrate.simps(integrand, rgrid);
                    else: sys.exit('Wrong integrator')
                # k2 loop end
            # k6 loop end
            toc = time.process_time()
            #    print('cycle',k1,'duration',toc-tic)
            k4 = k4 + 1
        # k1 loop end
    # k7 loop end
    return (k_start, k_num, FHHGOnScreen)



# Optimal workload is obtained by the same amount of load for each woker if possible, eventually one extra task for the last worker. Otherwise (NOT OPTIMAL!!!), every worker takes more load and some workers may be eventually not used. An optimal routine would be either balance the load among the workers or employ some sophisticated  parallel scheduler.
def ObtainWorkload(Nomega_points,W):
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

  return W, N_PointsGrid, N_PointsForProcess
# optimal workload is now given by the number of processes


def CoalesceResults(results,Nz_medium,Nz_anal,Nomega_anal_start,Nomega_points,Nr_anal,W):
    FHHGOnScreen = np.empty([Nz_medium,Nz_anal, Nomega_points, Nr_anal, 2], dtype=np.double)
    for k1 in range(W):  # loop over unsorted results
        for k5 in range(Nz_medium):
            for k2 in range(results[k1][1]):  # # of omegas computed by this worker
                for k3 in range(Nr_anal):  # results in the radial grid
                    for k4 in range(Nz_anal):  # results in the z grid
                        FHHGOnScreen[k5, k4, results[k1][0] - Nomega_anal_start + k2, k3, 0] = results[k1][2][k5][k4][k2][k3].real  # we adjust the field index properly to the field it just sorts matrices the following way [A[1], A[2], ...], the indices are retrieved by the append mapping
                        FHHGOnScreen[k5, k4, results[k1][0] - Nomega_anal_start + k2, k3, 1] = results[k1][2][k5][k4][k2][k3].imag
        #      FHHGOnScreen[ results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k2][k3]/(r_Bohr**2) # eventually fully in atomic units for the integral, but there is still a prefactor of the integral!!!
    return FHHGOnScreen









# define function to integrate, there are some global variables! ## THE OUTPUT IS IN THE MIX OF ATOMIC UNITS (field) and SI UNITS (radial coordinate + dr in the integral)
def FieldOnScreen_singleplane(z_medium, omegagrid, omega_step, rgrid, FField_r, rgrid_anal, zgrid_anal, k_start, k_num, integrator):
# this function computes the Hankel transform of a given source term in omega-domain stored in FField_r
# all the grids are specified in the inputs, except the analysis in omega_anal, this is specified by 'k_start' and 'k_num', it is used in the multiprocessing scheme
    Nz_anal = len(zgrid_anal); Nr_anal = len(rgrid_anal); Nr = len(rgrid);
    FHHGOnScreen = np.empty([Nz_anal, k_num, Nr_anal], dtype=np.cdouble)
    k4 = 0  # # of loops in omega
    for k1 in range(k_num):  # Nomega
        k5 = k_start + k1 * omega_step  # accesing the grid
        tic = time.process_time()
        for k6 in range(Nz_anal):
            for k2 in range(Nr_anal):  # Nr
                k_omega = omegagrid[k5] / (units.TIMEau * units.c_light);  # omega divided by time: a.u. -> SI # use this as a prefactor to obtain curvature
                integrand = np.empty([Nr], dtype=np.cdouble)
                for k3 in range(Nr): integrand[k3] = np.exp(-(rgrid[k3] ** 2) / (2.0 * (zgrid_anal[k6] - z_medium))) * rgrid[k3] * FField_r[k5, k3] * special.jn(0, k_omega * rgrid[k3] * rgrid_anal[k2] / (zgrid_anal[k6] - z_medium));  # rescale r to atomic units for spectrum in atomic units! (only scaling)
                if (integrator == 'Trapezoidal'): FHHGOnScreen[k6, k4, k2] = (1.0 / (zgrid_anal[k6] - z_medium)) * integrate.trapz(integrand, rgrid);
                elif (integrator == 'Simpson'): FHHGOnScreen[k6, k4, k2] = (1.0 / (zgrid_anal[k6] - z_medium)) * integrate.simps(integrand, rgrid);
                else: sys.exit('Wrong integrator')
            # k2 loop end
        # k6 loop end
        toc = time.process_time()
        #    print('cycle',k1,'duration',toc-tic)
        k4 = k4 + 1
    return (k_start, k_num, FHHGOnScreen)



def CoalesceResults_serial(results,Nz_anal,Nomega_anal_start,Nomega_points,Nr_anal,W):
    FHHGOnScreen = np.empty([Nz_anal, Nomega_points, Nr_anal, 2], dtype=np.double)
    for k1 in range(W):  # loop over unsorted results
        for k2 in range(results[k1][1]):  # # of omegas computed by this worker
            for k3 in range(Nr_anal):  # results in the radial grid
                for k4 in range(Nz_anal):  # results in the z grid
                    FHHGOnScreen[k4, results[k1][0] - Nomega_anal_start + k2, k3, 0] = results[k1][2][k4][k2][k3].real  # we adjust the field index properly to the field it just sorts matrices the following way [A[1], A[2], ...], the indices are retrieved by the append mapping
                    FHHGOnScreen[k4, results[k1][0] - Nomega_anal_start + k2, k3, 1] = results[k1][2][k4][k2][k3].imag
    #      FHHGOnScreen[ results[k1][0]-Nomega_anal_start+k2 , k3 ] = results[k1][2][k2][k3]/(r_Bohr**2) # eventually fully in atomic units for the integral, but there is still a prefactor of the integral!!!
    return FHHGOnScreen






###########################################################
#  there is the part for phenomenological dipoles:        #
###########################################################
#

# this should be loaded from somewhere or computed, or whatever... NOT directly in the code!
omegawidth = 4.0/np.sqrt(4000.0**2); # roughly corresponds to 100 fs
PhenomParams = np.array([
[1, 2, 35, 39], # harmonics
[0, 0, 1775., 3600.], # alphas
# [1, 29, 35, 39], # harmonics
# [0, 500., 1775., 3600.], # alphas
[omegawidth, omegawidth, omegawidth, omegawidth]
])
NumHarm = 2; # number of harmonics

## define dipole function
def dipoleTimeDomainApp(tgrid,z,r,I0,PhenomParams,tcoeff,rcoeff,omega0,zR):
# !!!! THERE ARE SEVERAL POINTS NEEDED TO MENTION HERE
#  - We don't follow the notation of the complex field exp(-i*(omega*t-...)). The reason is that we can then use fft instead of ifft... Be careful with conversions.
#  - One should then check the sign of alpha as well...

#  tcoeff = 4.0*np.log(2.0)*TIMEau**2 / ( TFWHMSI**2 )
#  rcoeff = 2.0/(w0r**2)
  # we are in atomic units, we go to electric field instead of intensity
  E0 = np.sqrt(I0) / np.sqrt(1.0+(z/zR)**2) # reduced in z
  lambd = mn.ConvertPhoton(omega0,'omegaau','lambdaSI')
  res = []
  for k1 in range(len(tgrid)):
    res1 = 0.0*1j;
    # intens = E0*np.exp(-tcoeff*(tgrid[k1])**2 - rcoeff*r**2) # qeff here?
    for k2 in range(NumHarm):
        # if (k2 == 2): intens = E0 * np.exp(-tcoeff * (tgrid[k1]-900.0) ** 2 - rcoeff * r ** 2)  # qeff here?
        # else: intens = E0 * np.exp(-tcoeff * (tgrid[k1]) ** 2 - rcoeff * r ** 2)  # qeff here?
        intens = E0 * np.exp(-tcoeff * (tgrid[k1]) ** 2 - rcoeff * r ** 2)  # qeff here?
        alpha = PhenomParams[1, k2]
        order = PhenomParams[0, k2]
        k_omega_wave = 2.0*np.pi*order/lambd
        IR_induced_phase = mn.GaussianBeamCurvaturePhase(r,z,k_omega_wave,zR) # order included
        res1 = res1 + intens*np.exp(1j*(tgrid[k1]*omega0*order-IR_induced_phase-alpha*intens))
    res.append(res1); ## various points in time
  return np.asarray(res)



def ComputeFieldsPhenomenologicalDipoles(I0SI,omega0,TFWHM,w0,tgrid,omegagrid,rgrid,z_medium):
  print('Computing phenomenological dipoles: FFTs')
  Nomega = len(omegagrid); Nr = len(rgrid); Nz_medium = len(z_medium);
  FField_r=np.empty([Nz_medium,Nomega,Nr], dtype=np.cdouble)
#   FField_r=np.empty([Nomega,Nr], dtype=np.cdouble)

  # coeffs applied in intensity
  # tcoeff = 4.0*np.log(2.0)*units.TIMEau**2 / ( TFWHM**2 )
  # rcoeff = 2.0/(w0**2)

  # but we apply them in the field
  zR = mn.GaussianBeamRayleighRange(w0, mn.ConvertPhoton(omega0, 'omegaau', 'lambdaSI'))
  tcoeff = 2.0*np.log(2.0)*units.TIMEau**2 / ( TFWHM**2 )



  k1 = mn.FindInterval(tgrid,0.0)
  dt = tgrid[k1+1]-tgrid[k1]

  for k1 in range(Nz_medium):
    w = w0 * np.sqrt(1.0 + (z_medium[k1] / zR) ** 2)  # reduced in z
    rcoeff = 1.0 / (w ** 2)
    for k2 in range(Nr): # multiprocessing?
      print('kr',k2)
      # the expression of the dipole is tricky, see appendix C in Jan Vabek's diploma thesis (using "two-times" Fourier in omega)
      dum = dipoleTimeDomainApp(tgrid,z_medium[k1],rgrid[k2],I0SI/units.INTENSITYau,PhenomParams,tcoeff,rcoeff,omega0,zR) # compute "complex" dipole in t-domain
      dum = (dt/np.sqrt(2.0*np.pi))*np.fft.fft(dum) # fft with normalisation
      # FField_r[k1,:,k2] = dum # !! CANNOT BE DONE THIS WAY, WE HAVE EXTRA ASSUMPTION THAT OUR SIGNAL IS REAL fft: C -> C
      for k3 in range(Nomega): FField_r[k1,k3,k2] = dum[k3]
  print('FFT computed');
  return  FField_r


# def ComputeFieldsPhenomenologicalDipoles_mp(I0SI,omega0,TFWHM,w0,tgrid,omegagrid,rgrid,z_medium):
#   print('Computing phenomenological dipoles: FFTs')
#   Nomega = len(omegagrid); Nr = len(rgrid); Nz_medium = len(z_medium);
#   FField_r=np.empty([Nz_medium,Nomega,Nr], dtype=np.cdouble)
# #   FField_r=np.empty([Nomega,Nr], dtype=np.cdouble)
#
#   tcoeff = 4.0*np.log(2.0)*units.TIMEau**2 / ( TFWHM**2 )
#   rcoeff = 2.0/(w0**2)
#   k1 = mn.FindInterval(tgrid,0.0)
#   dt = tgrid[k1+1]-tgrid[k1]
#
#   for k1 in range(Nz_medium):
#     # # define output queue
#     # output = mp.Queue()
#     for k2 in range(Nr): # multiprocessing?
#       print('kr',k2)
#       # the expression of the dipole is tricky, see appendix C in Jan Vabek's diploma thesis (using "two-times" Fourier in omega)
#       dum = dipoleTimeDomainApp(tgrid,rgrid[k2],I0SI/units.INTENSITYau,PhenomParams,tcoeff,rcoeff,omega0) # compute "complex" dipole in t-domain
#       # add the extra phase
#       dum = (dt/np.sqrt(2.0*np.pi))*np.fft.fft(dum) # fft with normalisation
#       # FField_r[k1,:,k2] = dum # !! CANNOT BE DONE THIS WAY, WE HAVE EXTRA ASSUMPTION THAT OUR SIGNAL IS REAL fft: C -> C
#       for k3 in range(Nomega): FField_r[k1,k3,k2] = dum[k3]
#   print('FFT computed');
#   return  FField_r