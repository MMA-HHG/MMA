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

### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factor to atomic units
TIME = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);




params={
'Nsim': 4096, # simulations in r #4096
'w0z' : 96.0e-6,
'r_extend' : 4.0,
'E0' : 0.075, 
'z' : 0.05,
'lambda' : 810.0e-9,
'phase0' : 0.0, # initial CEP
'N1cycles' : 20.0
}
dt=1





## load z-grid
file1=open("zgrid.dat","r")
lines = file1.readlines()
k1=0; zgrid=[]
for line in lines: zgrid.append(float(line)); k1=k1+1
Nz = k1; zgrid=np.asarray(zgrid);  
file1.close()  

print(zgrid[0])

## create rgrid
Nr = params['Nsim']
rgrid=[]
for k2 in range(Nr): rgrid.append(k2*params['r_extend']*params['w0z']/Nr)
rgrid=np.asarray(rgrid);
  


## remove the directories eventually

simpath = 'sims11'

if os.path.exists(simpath) and os.path.isdir(simpath):
  shutil.rmtree(simpath)
  print('deleted sims')

os.mkdir(simpath)


print(params['lambda'])

## calculator of params
def FieldParams(r,z,p):
  kwave = 2.0*np.pi/p['lambda']
  zR = (np.pi*(p['w0z'])**2)/p['lambda']
  wz = p['w0z']*np.sqrt(1.0+(z/zR)**2)
  if (z == 0.0):
    invRz = 0.0
  else:
    invRz = 1.0/(z + (zR**2)/z)
  Erz = p['E0']*(p['w0z']/wz)*np.exp(-(r/wz)**2)
  phase = p['phase0'] - 0.5*(r**2)*kwave*invRz + np.arctan(z/zR) # be careful with notation of the wave E = e^(-i*(omega*t-k*z)), it means that for a fixed z, phi should be nagative to correspond with E~sin(omega*t+phi0)
  return Erz, phase
  
test = os.path.join('a','b')
print(test)


## loop for making directories
for k1 in range(Nz):
  fold1 = str(k1+1); fold1 = fold1.zfill(6);
  fold1 = os.path.join(simpath,'z_'+fold1)
  os.mkdir(fold1)
  for k2 in range(Nr):
    # folder to be created
    fold2 = str(k2+1); fold2 = fold2.zfill(6);
    fold = os.path.join(fold1,'r_'+fold2)
    os.mkdir(fold)
    file1=open( os.path.join(fold,'param.txt') ,"a")    
    print(fold)
    Erz, phase = FieldParams(rgrid[k2],zgrid[k1],params)
    content = "// definition of calculation conditions for a numerical field\n"
    content = content + "2 : Efield.fieldtype // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical\n"
    content = content + "\n"
    content = content + "// definition of parameters\n"
    content = content + "-1. : Eguess // Energy of the initial state\n"
    content = content + "64000 : num_r // Number of points of the initial spatial grid 16000\n"
    content = content + "0 : num_exp // Number of points of the spatial grid for the expansion\n"
    content = content + "0.2 : dx // resolution for the grid\n"
    content = content + "0 : InterpByDTorNT // refine resolution only for numerical fields (0 - by dt, 1 - by number of points)\n"
    content = content + str(dt)+" : dt // resolution in time\n"
    content = content + "1 : Ntinterp // number of intermediate points\n"
    content = content + "200. : textend // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700\n"
    content = content + "0 : analy.writewft // writewavefunction (1-writting every tprint)\n"
    content = content + "10. : analy.tprint // time spacing for writing the wavefunction\n"
    content = content + "2. : x_int // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields\n"
    content = content + "0 : PrintGaborAndSpectrum // print Gabor and partial spectra (1-yes)\n"
    content = content + "8 : a_Gabor // the parameter of the gabor window [a.u.]\n"
    content = content + "15. : omegaMaxGabor // maximal frequency (source term) in Gabor [a.u.]\n"
    content = content + "10 : dtGabor // spacing in Gabor (source term) [a.u.]\n"
    content = content + "2000. : tmin1window); // analyse 1st part of the dipole\n"
    content = content + "5000. : tmax1window); // analyse 1st part of the dipole\n"
    content = content + "5250. : tmin2window // analyse 2nd part of the dipole\n"
    content = content + "10000. : tmax2window // analyse 2nd part of the dipole\n"
    content = content + "2 : PrintOutputMethod // (0 - only text, 1 - only binaries, 2 - both)\n"
    content = content + "\n"
    content = content + "// TARGET definition:\n"
    content = content + "1.189 : trg.a\n"
    content = content + "\n"
    content = content + "// GAUGES\n"
    content = content + "0 : gauge // 0-length, otherwise velocity, velocity available only for analytic field (A needed)\n"
    content = content + "0 : transformgauge // 1 - transform also to another gauge during the calculation, (A needed)\n"
    content = content + "\n"
    content = content + "// PARAMETERS OF SIN2\n"
    content = content + str(Erz) + " : amplitude\n"
    omega_au = (2.0*np.pi*hbar*inverse_alpha_fine**2)/(params['lambda']*elmass*c_light)# find frequency in atomic units
    content = content + str(omega_au) + " : frequency\n"
    content = content + "0.0 : initial time\n"
    content = content + str(params['N1cycles']) + " : # of cycles\n"
    content = content + str(phase) + " : CEP [in radians, reference is a cosine pulse in A]\n"
    file1.write(content)
    

np.savetxt(os.path.join(simpath,'rgrid.dat'), rgrid, fmt="%e")
shutil.copyfile('zgrid.dat',os.path.join(simpath,'zgrid.dat'))


#// definition of calculation conditions for a numerical field
#0 : Efield.fieldtype // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical



#// definition of parameters
#-1. : Eguess // Energy of the initial state
#64000 : num_r // Number of points of the initial spatial grid 16000
#0 : num_exp // Number of points of the spatial grid for the expansion
#0.2 : dx // resolution for the grid
#0 : InterpByDTorNT // refine resolution only for numerical fields (0 - by dt, 1 - by number of points)
#2 : dt // resolution in time 0.0625
#1 : Ntinterp // number of intermediate points
#200. : textend // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
#0 : analy.writewft // writewavefunction (1-writting every tprint)
#10. : analy.tprint // time spacing for writing the wavefunction
#2. : x_int // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
#1 : PrintGaborAndSpectrum // print Gabor and partial spectra (1-yes)
#8 : a_Gabor // the parameter of the gabor window [a.u.]
#15. : omegaMaxGabor // maximal frequency (source term) in Gabor [a.u.]
#10 : dtGabor // spacing in Gabor (source term) [a.u.]
#2000. : tmin1window); // analyse 1st part of the dipole
#5000. : tmax1window); // analyse 1st part of the dipole
#5250. : tmin2window // analyse 2nd part of the dipole
#10000. : tmax2window // analyse 2nd part of the dipole
#2 : PrintOutputMethod // (0 - only text, 1 - only binaries, 2 - both)




#// definition of calculation conditions for a numerical field
#2 : Efield.fieldtype // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical

#// definition of parameters
#-1. : Eguess // Energy of the initial state
#64000 : num_r // Number of points of the initial spatial grid 16000
#0 : num_exp // Number of points of the spatial grid for the expansion
#0.2 : dx // resolution for the grid
#0 : InterpByDTorNT // refine resolution only for numerical fields (0 - by dt, 1 - by number of points)
#0.00125 : dt // resolution in time
#1 : Ntinterp // number of intermediate points
#200. : textend // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
#0 : analy.writewft // writewavefunction (1-writting every tprint)
#10. : analy.tprint // time spacing for writing the wavefunction
#2. : x_int // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields

#// TARGET definition:
#1.189 : trg.a

#// GAUGES
#0 : gauge // 0-length, otherwise velocity, velocity available only for analytic field (A needed)
#0 : transformgauge // 1 - transform also to another gauge during the calculation, (A needed)

#// PARAMETERS OF SIN2
#0.3 : amplitude
#0.057 : frequency
#0.0 : initial time
#90.0 : # of cycles
#0.0 : CEP [in radians, reference is a cosine pulse in A]



##zpadding=5;
#echo $Nzsim

#rm -r sims
#mkdir sims

## length of file to do correct padding

#cd sims
#k1=0;
#while read LINE; do
#let k1=k1+1; printf -v kdum "%0${zpadding}d" $k1;
#	mkdir z_$kdum
#	cd z_$kdum
#	z=${LINE//E/d};
#	echo 'value ' $z
#	for k2 in `seq -w 0 $Nsim`; do
#		mkdir r_$k2


################################################## INSIDE OF THE EVERY FOLDER MADE
#dum=$(./../../calculator.out << INPUTS
#$E0
#$lambda
#$z
#$w0z
#$r_extend
#$Nsim
#$k2
#$phase0
#INPUTS
#)
#	read -a c <<< $dum;
#	Erz=${c[0]}; phase=${c[1]}; ## input in c code not FORTRAN
##	echo ${c[*]}
##	echo $Erz
##	echo $phase


#{ echo 	"// definition of calculation conditions for a numerical field
#2 : Efield.fieldtype // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical

#// definition of parameters
#-1. : Eguess // Energy of the initial state
#64000 : num_r // Number of points of the initial spatial grid 16000
#0 : num_exp // Number of points of the spatial grid for the expansion
#0.2 : dx // resolution for the grid
#0.025 : dt // resolution in time
#200. : textend // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
#0 : analy.writewft // writewavefunction (1-writting every tprint)
#10. : analy.tprint // time spacing for writing the wavefunction
#2. : x_int // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields

#// TARGET definition:
#0.695 : trg.a

#// GAUGES
#0 : gauge // 0-length, otherwise velocity, velocity available only for analytic field (A needed)
#0 : transformgauge // 1 - transform also to another gauge during the calculation, (A needed)

#// PARAMETERS OF SIN2
#$Erz : amplitude
#0.057 : frequency
#0.0 : initial time
#10.0 : # of cycles
#$phase : CEP [in radians, reference is a cosine pulse in A]"
#} > r_$k2/param.txt



#cp ../../programs/submit_1DTDSE.slurm r_$k2/submit_1DTDSE.slurm
#cp ../../programs/TDSE/TDSE1D.out r_$k2/TDSE1D.out

## now submit the job

#cd r_$k2
#	mkdir results
#	sbatch -J TDSE_z_{$kdum}_r_${k2} submit_1DTDSE.slurm
#cd ..




################################################### END
#	done
#	cd ..
#done < ../zgrid.dat


#cd ..


##j=1
##while read LINE; do 
##	read -a run <<< $LINE;
##	echo "${run[*]}" >> "params_check.txt";
##	I0=${run[0]}; press=${run[1]}; w0=${run[2]}; focus=${run[3]}; tau=${run[4]}; n2atm=${run[5]}; lambda=${run[6]}; THEORY=${run[7]}; medium_length=${run[8]}; PRINTFIELDS=${run[9]};
##	export I0 press w0 focus tau n2atm lambda j THEORY medium_length PRINTFIELDS; ./script2.sh; let j++; 
##done < params.txt


###run=(3 4); let j++; n2atm=${run[0]}; n2atm=n2atm/2.0 | bc; export n2atm j; ./script2.sh;

##echo "test4"


##for fold in sim*; do echo $fold; echo $fold >> "ListOfSimulations.txt"; done;

### clearing folder

##mv ListOfSimulations.txt logfiles/
##mv params_check.txt logfiles/
##mv *.dat logfiles/
##rm *.out
##rm script2.sh

###if (( $(echo "$focus < 0" | bc -l) )); then
###  echo "no lense"
###fi



#echo "done" 
