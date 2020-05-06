import numpy as np
import os
import time
import multiprocessing as mp
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
# - we use precomputed dipoles now

#### THE MAIN PROGRAM #####
class LaserParamsClass:
  def __init__(self):
    pass

class NumericalParamsClass:
  def __init__(self):
    pass

LaserParams = LaserParamsClass() ## define macroscopic gaussian beam # try also fancy reading directly here
NumericalParams = NumericalParamsClass()

###################### THE PARAMETERS OF SIMULATION
ParamFile = 'results.h5'
ParamFile = h5py.File(ParamFile,'r')

# z_medium = np.array([-0.025, -0.02, -0.015, -0.01, -0.005, 0.0, 0.01])  #np.asarray([-0.025, -0.02, -0.015, -0.01, -0.005, 0.0, 0.01])  # np.array([-0.003, 0.0, 0.003]);

if ( 'array' == mn.readscalardataset(ParamFile,'inputs/'+'jetpositions_type','S') ):
  NumericalParams.z_medium = ParamFile['inputs/'+'jetpositions'][()]
elif ( 'grid' == mn.readscalardataset(ParamFile,'inputs/'+'jetpositions_type','S') ):
  NumericalParams.z_medium_max = mn.readscalardataset(ParamFile,'inputs/'+'jetpositions_max','N'); NumericalParams.z_medium_min = mn.readscalardataset(ParamFile, 'inputs/' + 'jetpositions_min', 'N'); NumericalParams.z_medium_Npoints = mn.readscalardataset(ParamFile, 'inputs/' + 'jetpositions_Npoints', 'N')
  NumericalParams.z_medium = np.linspace(NumericalParams.z_medium_min, NumericalParams.z_medium_max, NumericalParams.z_medium_Npoints);
else: sys.exit('wrongly specified jetpositions')
z_medium = NumericalParams.z_medium

NumericalParams.storing_source_terms = mn.readscalardataset(ParamFile,'inputs/'+'storing_source_terms','S')


NumericalParams.diffraction_integral = mn.readscalardataset(ParamFile,'inputs/'+'diffraction_integral','S')

LaserParams.I0 = mn.readscalardataset(ParamFile,'inputs/'+'I0','N')
LaserParams.w0 = mn.readscalardataset(ParamFile,'inputs/'+'w0','N')
LaserParams.lambd = mn.readscalardataset(ParamFile,'inputs/'+'lambda','N') #800.0e-9, # must correspond

LaserParams.omega0 = mn.ConvertPhoton(LaserParams.lambd ,'lambdaSI','omegaau') # find frequency in atomic units
LaserParams.zR = np.pi*(LaserParams.w0**2)/LaserParams.lambd
omega0 = LaserParams.omega0; zR = LaserParams.zR

## Optical system (part of laser)
LaserParams.r_pinhole = mn.readscalardataset(ParamFile,'inputs/'+'r_pinhole','N')
LaserParams.z_pinhole = mn.readscalardataset(ParamFile,'inputs/'+'z_pinhole','N')

# anlyses params # at the moment optimised for t he intensity list, change later


# rmax = 3.0*;
rmax = mn.readscalardataset(ParamFile,'inputs/'+'rmax_fact','N')*LaserParams.w0
# Nr = 300;
Nr = mn.readscalardataset(ParamFile,'inputs/'+'Nr_int','N') # 8193; # 2049; #1000; # 2049

rmax_anal = mn.readscalardataset(ParamFile,'inputs/'+'rmax_anal','N') # 3*1e-3 # [SI] on screen # 0.0001
Nr_anal = mn.readscalardataset(ParamFile,'inputs/'+'Nr_anal','N') # 25 #250 #750

zmin_anal = mn.readscalardataset(ParamFile,'inputs/'+'zmin_anal','N') # 0.05 # !!!!!! in the reference of the jet, the grid is then reshaped correctly
zmax_anal = mn.readscalardataset(ParamFile,'inputs/'+'zmax_anal','N') # 0.4
only_one_plane_bool = 1 == mn.readscalardataset(ParamFile,'inputs/'+'only_one_plane','N') # True # only zmax used # Nz_anal is overrun
fix_zmax_bool = 1 == mn.readscalardataset(ParamFile,'inputs/'+'fix_zmax','N') # True
Nz_anal = mn.readscalardataset(ParamFile,'inputs/'+'Nz_anal','N') # 100 #200

Hmin_anal = mn.readscalardataset(ParamFile,'inputs/'+'Hmin_anal','N') # 28.5 # 0.0 #28.5 np.nan
if ( Hmin_anal < 0 ): Hmin_anal = np.nan
Hmax_anal = mn.readscalardataset(ParamFile,'inputs/'+'Hmax_anal','N') # 29.5 # 2.5 #29.5 71.5
if ( Hmax_anal < 0 ): Hmax_anal = np.nan
omega_step = mn.readscalardataset(ParamFile,'inputs/'+'omega_step','N') # 1

# used only for phenomenological dipoles
tcoeff = 6.0; # extension of tgrid in the units of TFWHM
Nt = 1000;


## other parameters
integrator={}
integrator['method'] = mn.readscalardataset(ParamFile,'inputs/'+'I_method','S') #'Romberg'; # 'Trapezoidal', Simpson, Romberg
integrator['tol'] = mn.readscalardataset(ParamFile,'inputs/'+'I_tol','N') # 1e-2;
integrator['n0'] = mn.readscalardataset(ParamFile,'inputs/'+'I_n0','N') # 2;
dipole_model = mn.readscalardataset(ParamFile,'inputs/'+'dipole_model','S') # 'IntensityList' # 'IntensityList', Phenomenological

W = mn.readscalardataset(ParamFile,'inputs/'+'num_of_processes','N') # 4; # W = mp.cpu_count() # this is the number of workers


if ( 'local' == mn.readscalardataset(ParamFile,'inputs/'+'computer','S') ):
  outpath = os.path.join("/mnt", "c", "data", "ThinTargets_collab", "loc_tests")
  IntensityListFile = os.path.join("/mnt", "c", "data", "ThinTargets_collab", mn.readscalardataset(ParamFile, 'inputs/' + 'IntensityListFileName','S'))  # used only for the list
elif ( 'occigen' == mn.readscalardataset(ParamFile,'inputs/'+'computer','S') ):
  outpath = os.getcwd()
  IntensityListFile = os.path.join("/scratch", "cnt0025", "cli7594", "vabekjan", "ThinTargets_collab", "1DTDSE", "data", mn.readscalardataset(ParamFile, 'inputs/' + 'IntensityListFileName','S'))  # used only for the list
else: sys.exit('wrongly specified computer')

## output file specification
OutputFileName = mn.readscalardataset(ParamFile,'inputs/'+'OutputFileName','S') # "romb_iters_test.h5" # "results_phenom8.h5"
OutputFileAccessMethod = 'r+'

if ( 'single' == mn.readscalardataset(ParamFile,'inputs/'+'precision','S') ): precision = 'f'
elif ( 'double' == mn.readscalardataset(ParamFile,'inputs/'+'precision','S') ): precision = 'd'
else: sys.exit('precision wrongly specified')


## Cases of more advanced schemes
if NumericalParams.diffraction_integral == "CircularAperture_2D" :
  NumericalParams.Nr2 = mn.readscalardataset(ParamFile, 'inputs/' + 'Nr_int2', 'N')  # 8193; # 2049; #1000; # 2049
  NumericalParams.rgrid2 = np.linspace(0.0, LaserParams.r_pinhole, NumericalParams.Nr2)


ParamFile.close()
shutil.copy('results.h5', os.path.join(outpath,OutputFileName))


print(omega0,'omega0 in a.u.')


#quit()
########################################################## THE BODY OF THE PROGRAM


if (dipole_model == 'IntensityList'):
  # loading procedure
  file1 = h5py.File(IntensityListFile, 'r')
  NumericalParams.Igrid = file1['/Igrid'][:]
  NumericalParams.omegagrid = file1['/omegagrid'][:]
  NumericalParams.FSourceterm = file1['/FDipoleAccelerations'][:]
  NumericalParams.FSourceterm = np.squeeze(NumericalParams.FSourceterm[:,:,0] + 1j*NumericalParams.FSourceterm[:,:,1]) # convert to complex numbers

elif (dipole_model == 'Phenomenological'):
  tgrid = np.linspace(-tcoeff * 0.5 * LaserParams.TFWHM / units.TIMEau, tcoeff * 0.5 * LaserParams.TFWHM / units.TIMEau, Nt)
  Nomega = len(tgrid)//2 + 1
  NumericalParams.omegagrid = np.linspace(0,2.0*np.pi*Nomega/(tcoeff*LaserParams.TFWHM/units.TIMEau),Nomega)





## create grids
Nomega = len(NumericalParams.omegagrid)
Hgrid = np.empty([Nomega], dtype=np.double)
Hgrid[:] = NumericalParams.omegagrid[:]/LaserParams.omega0
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
  omegamin_anal = NumericalParams.omegagrid[0]
  Nomega_anal_start = 0
else:
  omegamin_anal = omega0*Hmin_anal
  Nomega_anal_start = mn.FindInterval(NumericalParams.omegagrid,omegamin_anal)

if (np.isnan(Hmax_anal)):
  omegamax_anal = NumericalParams.omegagrid[-1]
  Nomega_anal = len(NumericalParams.omegagrid)-1
else:
  omegamax_anal = omega0 * Hmax_anal
  Nomega_anal = mn.FindInterval(NumericalParams.omegagrid, omegamax_anal)

Nomega_points = mn.NumOfPointsInRange(Nomega_anal_start,Nomega_anal,omega_step); # Nomega_points is the number of simulations we want to perform



## compute fields in our rgrid
if ( (dipole_model == 'IntensityList') and (NumericalParams.storing_source_terms == 'table') ):
  NumericalParams.FField_r = Hfn.ComputeFieldsInRFromIntensityList(z_medium, rgrid, Hgrid, Nomega, LaserParams, NumericalParams.Igrid, NumericalParams.FSourceterm)
  del NumericalParams.FSourceterm; del NumericalParams.Igrid;  # this table may be big
elif (dipole_model == 'Phenomenological'):
  omegagrid = [];
  NumericalParams.FField_r = Hfn.ComputeFieldsPhenomenologicalDipoles(LaserParams.I0,omega0,LaserParams.TFWHM,LaserParams.w0,tgrid,omegagrid,rgrid,z_medium)




## print some analyses outputs
print('om_max', Nomega_anal)
print('om_min', Nomega_anal_start)
# print('tmax',tgrid[len(tgrid)-1])
print('omax',NumericalParams.omegagrid[-1])

  
dr = rgrid[1]-rgrid[0]






omegagrid_anal=[]

print('Nomega_points = ', Nomega_points);

  

#########################################################################################
#### MAIN INTEGRATION ####
tic1 = time.process_time()
ttic1 = time.time()

NumericalParams.integrator = integrator
NumericalParams.rgrid = rgrid
NumericalParams.rmax= rmax

NumericalParams.z_medium = z_medium
NumericalParams.omega_step = omega_step
NumericalParams.rgrid_anal = rgrid_anal
NumericalParams.zgrid_anal = zgrid_anal




# define output queue
output = mp.Queue()

# passing by reference is unPythonic, we define the extra function though
def FieldOnScreen_handle(k_start, k_num, NumericalParams, LaserParams):

  if NumericalParams.diffraction_integral == "CircularAperture_2D" :
    res = Hfn.FieldOnScreenApertured2D1(k_start, k_num, NumericalParams, LaserParams)
  elif NumericalParams.diffraction_integral == "CircularAperture_analytic" :
    res = Hfn.FieldOnScreenApertured1(k_start, k_num, NumericalParams, LaserParams)
  else: res = Hfn.FieldOnScreen(k_start, k_num, NumericalParams, LaserParams) # default

  output.put(res)


# split the workload evenly
W, N_PointsGrid, N_PointsForProcess = Hfn.ObtainWorkload(Nomega_points,W)

print('----')
print('process grids for workers: starting point + loads')
print(N_PointsGrid)
print(N_PointsForProcess)
print('----')

# now we need to retrieve the number of loops for each worker in general cases...


### we use multiprocessing by assigning each part of the load as one process
# define processes
processes = [mp.Process(target=FieldOnScreen_handle, args=(Nomega_anal_start+N_PointsGrid[k1], N_PointsForProcess[k1], NumericalParams, LaserParams)) for k1 in range(W)]

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
for k1 in range(Nomega_anal_start,Nomega_anal,omega_step): omegagrid_anal.append(NumericalParams.omegagrid[k1])
omegagrid_anal=np.asarray(omegagrid_anal);


# coalesce the results
FHHGOnScreen = Hfn.CoalesceResults(results,Nz_medium,Nz_anal,Nomega_anal_start,Nomega_points,Nr_anal,W)



#################################################################################
#### PROCESS OUTPUTS ####

    
#file1=open("Spectrum.dat","w")
# np.savetxt(os.path.join(outpath,"Spectrumreal.dat"),FHHGOnScreen[0,0,:,:,0].real,fmt="%e")
# np.savetxt(os.path.join(outpath,"Spectrumimag.dat"),FHHGOnScreen[0,0,:,:,1].imag,fmt="%e")
# np.savetxt(os.path.join(outpath,"omegagrid_anal.dat"),omegagrid_anal,fmt="%e")
# np.savetxt(os.path.join(outpath,"rgrid_anal.dat"),rgrid_anal,fmt="%e")
# np.savetxt(os.path.join(outpath,"zgrid_anal.dat"),zgrid_anal[0,:],fmt="%e")

#HDF5 results
f = h5py.File(os.path.join(outpath,OutputFileName),OutputFileAccessMethod)
grp = f.create_group('XUV')

h5spectrum_r = grp.create_dataset('Spectrum', data=FHHGOnScreen, dtype=precision)
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
mn.adddataset(grp,'w0',LaserParams.w0,'[SI]')
mn.adddataset(grp,'I0',LaserParams.I0,'[SI]')

## numerical parameters
grp = f.create_group('numerics')
mn.adddataset(grp,'integral_rmax',rgrid[-1],'[SI]')
mn.adddataset(grp,'integral_rmin',rgrid[0],'[SI]')
mn.adddataset(grp,'integral_dr',dr,'[SI]')
mn.adddataset(grp,'integral_points',Nr,'[-]')


f.close()

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

