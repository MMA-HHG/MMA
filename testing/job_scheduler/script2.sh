#!/bin/bash


#pass parameters to the calculator #available inps: n2atm peak_intens j


#work with directories
mkdir sim$j


#echo $n2atm
dum=$(./calculator.out << INPUTS
$I0
$press
$w0
$focus
$tau
$n2atm
$lambda
INPUTS
)

read -a c <<< $dum;

echo ${c[*]}

PinPcr=${c[0]};
w0z=${c[1]};
focal_length=${c[2]};

unset c;



# recalculate to input units
dum=$(./convert.out << INPUTS
$lambda
$w0z
$focal_length
$tau
INPUTS
)

echo $dum

read -a c <<< $dum;

#k1=0;
#while read LINE; do c[k1]=$LINE; let k1++; done < converted.dat

lambda=${c[0]//E/d};
w0z=${c[1]//E/d};
focal_length=${c[2]//E/d};
tau=${c[3]//E/d};

echo "lambda = $lambda";

if (( $(echo "$focus < 0" | bc -l) )); then
  focal_length="0.D0"
fi



cd sim$j

      # Generate the PARAMETRIC FILE
{ echo 	"	CODE PARAMETERS
Number of processors                                   : 32
Run time in hours                                      : 23.5d0
Length of window for t, normalized to pulse duration   : 8.d0
Number of points in t                                  : 2048
Length of window for r, normalized to beamwaist        : 4.d0
Number of points in r                                  : 2048
Number of absorber points in time                      : 16
Phase threshold for decreasing delta_z                 : 2.d-3
Physical distance of propagation (m)                   : ${medium_length//E/d}
Physical output distance for matlab files (m)          : 20.d-4
Output distance in z-steps for fluence and power       : 100
Radius for diagnostics (mm)                            : 0.15d0
Physical first stepwidth (mm)                          : 1.d-2
Operators T,T-1                                        : 2
[ 1- nothing / 2- Krausz / 3- uppe / 4- uppe without harmonics ]

	LASER PARAMETERS
Laser wavelength (cm)                                  : $lambda
Beamwaist (cm)                                         : $w0z
Degree of SuperGaussian                                : 1
Pulse duration in 1/e (fs)                             : $tau
Degree of SuperGaussian in time                        : 1
Ratio Pin/Pcr                                          : ${PinPcr//E/d}
Input                                                  : 1
[ 1- Gaussian temporal shape / 2- temporal shape from file / 3- CONTINUATION ]
filename for method 2                                  : nrl.dat
filename for method 3                                  : 000_000000
amplituderatio for method 3                            : 1.d0
spatial noise on the input shape                       : 0.D0
temporal noise on the input shape                      : 0.D0
noise on the input shape                               : 0.D0
focal length in the medium cm (0 for no lense)         : $focal_length
Initial chirp phase                                    : 0.d0      

        MEDIUM PARAMETERS
pressure in bar                                        : $press
SPECIFY ALL OTHER PARAMETERS AT 1 BAR, THE CODE WILL EVALUATE THE VALUES AT THE GIVEN PRESSURE

Type of dispersion law                                 : 3
[ 1- Taylor expansion / 2- Silica Agrawal / 3- Argon Dalgarno / 4- Miro data / 5- air / 6- xenon / 7- neon / 8- KDPo]
filename for method 4                                  : waterchi.tab
For method 1 only:
Linear refractive index                                : 1.45d0
inverse GV coefficient - n0/c (fs/cm)                  : 0.d0
GVD coefficient (fs2/cm)                               : -279.d0
Third order dispersion coefficient (fs3/cm)            : 1510.d0
Fourth order dispersion coefficient (fs4/cm)           : -4930.d0
Fifth order dispersion coefficient ( fs5/cm)           : 23245.d0

Nonlinear refractive index, Kerr coefficient (cm2/W)   : ${n2atm//E/d}
Type of Delayed Kerr response                          : 1
[ 1-No delayed Kerr / 2-Response in xdk*exp(t/tdk) / 3-Response in xdk*exp(t/tdk)*sin(wr t) ]
ratio of delayed Kerr, xdk                             : 0.5d0
Time of delayed Kerr, tdk (fs)                         : 77.d0
Frequency in delayed Kerr, wr (fs-1)                   : 1.6d-2

Chi5 coefficient (cm4/W2)                              : 0.D0

Effective density of neutral molecules (cm-3)          : 2.7d19
Ionization poential of neutral molecules (eV)          : 15.8d0
Initial electron density (cm-3)                        : 0.d0
Type of ionization method                              : $THEORY
[1-Euler method / 2- Semi-analytical / 3- PPT / 4-ADK-molecular / 5- Keldysh for crystal only / 6- KDP Guillaume / 7- air Pavel]
MPI cross section for method 1-2 (s-1cm2K/WK)          : 1.9d-120
angular momentum for method 3,7                        : 1
Effective residue charge for method 3-4,7              : 1.d0
Reduced mass of hole-electron for method 5             : 0.5d0
Number of photons to ionize from SLG1 (method 6)       : 3
cross section to ionize from SLG1 (m 6, s-1cm2Kp/WKp)  : 8.6d-27
density of defects SLG1 (cm-3)                         : 2.D17
cross section to ionize from SLG2 (m 6, s-1cm2/W)      : 2.d0
sigmacvref for I_ref (method 6, cm2/s)                 : 4.35d-7
reference intensity I_ref (method 6, W/cm^2)           : 42.32d12
reference exponent exp_ref (method 6)                  : -3.3d0
Number of photons to populate SLG2 (method 6)          : 2
cross section to populate SLG2 (m 6, s-1cm2Kpp/WKpp)   : 1.3d-12
saturation density for SLG2 (method 6,  cm-3)          : 2.D17
Effective density of neutral N2 molecules (m 7, cm-3)  : 2.2d19
Ionization poential of neutral N2 molecules (m 7, eV)  : 15.6d0
angular momentum N2 for method 7                       : 0
Effective residue charge N2 for method 7               : 0.9d0
Initial free electron temperature (m 7, eV)            : 0.025d0

electron colision time (method 1-6, fs)                : 190.D0
Linear recombination coefficient (fs-1)                : 0.D0
Linear recombination coefficient (6: SLG1 elec.) (fs-1): 3.3D-3
Linear recombination coefficient (6: holes) (fs-1)     : 1.D-3
quadratic recombination (gasses) fs-1cm3               : 0.D0

number of photons involved in the n-absorption         : 2  
the n-photon absoption cross section [s-1cm2N/Wn]      : 0.d0 
density of absorbing molecules                         : 0.d0" 
} > input.inp

cd ..
cp -r starting/* sim$j/
cp sim$j/input.inp input$j.inp

cp convert.out sim$j/convert.out


#for i in {1..5}
#do
#   let j=i
#   echo "loop $j"
#   # cp file2 inputs/input$j.dat
#   # cp -r direct direct$j
#   cd intensity_study
#   mkdir sim$j	
#   cd ..
#   cp -r starting/* intensity_study/sim$j/
#   # rm -r direct$j
#done


#echo $n2atm
#echo $j


