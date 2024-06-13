/**
 * @file constants.h
 * @brief Header containing physical constants declaration and initialization.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef CONSTANTS_H
#define CONSTANTS_H

extern const double Pi;
extern double Ip_HeV, hbar, alpha_fine, c_light, elcharge, elmass, mu0, eps0, r_Bohr, 
    TIMEau, EFIELDau, k_Boltz, absolute_zero, torr2SI;
void Init_constants();
void Soft_Coulomb_parameters(char *, double *);
void Soft_Coulomb_ground_states(char *, double *);


#endif