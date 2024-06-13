/**
 * @file constants.c
 * @brief Contains physical constants.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#define _USE_MATH_DEFINES
#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif
#include<math.h>
#include<stdlib.h>
#include <string.h>
#include <stdio.h>

#include "constants.h"

const double Pi = M_PI;
double Ip_HeV, hbar, alpha_fine, c_light, elcharge, elmass, mu0, eps0, r_Bohr, 
    TIMEau, EFIELDau, k_Boltz, absolute_zero, torr2SI;

/**
 * @brief Inits physical constants
 */
void Init_constants(){
    Ip_HeV = 27.21138602;
    hbar = 1.054571800e-34;
    alpha_fine = 1/137.035999139;
    c_light = 299792458.;
    elcharge = 1.602176565e-19;
    elmass = 9.10938356e-31;
    mu0 = 4.0*Pi*1e-7;
    eps0 = 1.0/(mu0*c_light*c_light);
    r_Bohr = 4.0*Pi*eps0*hbar*hbar/(elmass*elcharge*elcharge);
    TIMEau = (elmass*r_Bohr*r_Bohr)/hbar;
    EFIELDau = hbar*hbar/(elmass*r_Bohr*r_Bohr*r_Bohr*elcharge);
    k_Boltz = 1.38064852e-23;
    absolute_zero = -273.15;
    torr2SI = 101325./760.;
}

/**
 * @brief Tabulated values for 1D-soft-Coulomb
 * 
 * see table 7.1 of https://theses.hal.science/tel-04192431v1/document
 * 
 */
void Soft_Coulomb_parameters(char * gas_type, double * a_value){

    printf("Setting gas preset\n");
    if (strcmp(gas_type, "He") == 0) 
    {
       *a_value = 0.6950;
    } 
    else if (strcmp(gas_type, "Ne") == 0)
    {
        *a_value = 0.8161;
    }
    else if (strcmp(gas_type, "Ar") == 0)
    {
        *a_value = 1.1893;
    }
    else if (strcmp(gas_type, "Kr") == 0)
    {
        *a_value = 1.3676;
    }
    else if (strcmp(gas_type, "Xe") == 0)
    {
        *a_value = 1.6171;
    }
    else
    {
        printf("Soft-coulomb parameter for this input not supported\n");
        exit(EXIT_FAILURE);
    }

}

/**
 * @brief Tabulated ground_states for 1D-soft-Coulomb
 * 
 * See table 7.1 of https://theses.hal.science/tel-04192431v1/document
 *  (It still needs to be refined for TDSE calculations.)
 * 
 */
void Soft_Coulomb_ground_states(char * gas_type, double * E_GS){

    printf("Setting gas preset\n");
    if (strcmp(gas_type, "He") == 0) 
    {
       *E_GS = 0.9036;
    } 
    else if (strcmp(gas_type, "Ne") == 0)
    {
        *E_GS = 0.7924;
    }
    else if (strcmp(gas_type, "Ar") == 0)
    {
        *E_GS = 0.5792;
    }
    else if (strcmp(gas_type, "Kr") == 0)
    {
        *E_GS = 0.5145;
    }
    else if (strcmp(gas_type, "Xe") == 0)
    {
        *E_GS = 0.4458;
    }
    else
    {
        printf("Soft-coulomb parameter for this input not supported\n");
        exit(EXIT_FAILURE);
    }

}
