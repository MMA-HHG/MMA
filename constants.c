#include<time.h> 
#include<stdio.h>
#include <mpi.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "util.h"


void Init_constants(void){
    double pi = acos(-1.0);
    Ip_HeV = 27.21138602;
    hbar = 1.054571800e-34;
    alpha_fine = 1/137.035999139;
    c_light = 299792458;
    elcharge = 1.602176565e-19;
    elmass = 9.10938356e-31;
    mu0 = 4.0*pi*1e-7;
    eps0 = 1.0/(mu0*c_light*c_light);
    r_Bohr = 4.0*pi*eps0*hbar*hbar/(elmass*elcharge*elcharge);
}