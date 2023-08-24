/**
 * @file structures.h
 * @brief Header containing structures and their methods used throughout the code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef STRUCTURES
#define STRUCTURES

#include <stdlib.h>

// various analytic envelopes (not reintroduced yet)
typedef struct sin2_definition {
	double oc, o, A0, nc, phi, phi0, ti, E0;
} sin2_definition;

typedef struct Esin2_definition {
	double oc, o, A0, nc, phi, phi0, ti, E0;
} Esin2_definition;

typedef struct trap_definition {
	double omega,E0,phi,ton,toff;
	int nc;
} trap_definition;

typedef struct flattop1_definition {
	double ton,toff,T,o,phi,A,ti;
} flattop1_definition;

typedef struct flattop1chirp_definition {
	double ton,toff,T,o,phi,A,ti,b,c;
} flattop1chirp_definition;

/**
 * @brief Full specification of the input field
 * 
 */
typedef struct Efield_var {
	int fieldtype,fieldtype2;
	double *tgrid;
	double *Field;
	double dt;
	int Nt;
	struct trap_definition trap;
	struct sin2_definition *sin2;
	struct Esin2_definition *Esin2;
	struct flattop1_definition *flt1;
	struct flattop1chirp_definition *flt1ch;
	int Nflt1,Nsin2,NEsin2,Nflt1ch;
	double omega,E0,phi,ton,toff;
	int nc;
} Efield_var;

/**
 * @brief Microscopic target specification
 * 
 * @param a (double) Parameter of the soft core potential.
 */
typedef struct trg_def{
	double a;
} trg_def;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary structures for drivers of I/O
typedef struct analy_def{
	double tprint;
	int writewft;
} analy_def;

typedef struct output_print_def{
	int Efield, FEfield, sourceterm, Fsourceterm, FEfieldM2, FsourceTermM2, PopTot, tgrid, omegagrid, PopInt, expval_x;
} output_print_def;

////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O structures & functions operating on them
typedef struct inputs_def{
	struct trg_def trg;
	struct Efield_var Efield;
	double Eguess, Einit;
	double tmin;
	int Nt;
	int num_t;
	double dt;
	int num_r;
	int num_exp;
	double dx;
	double *psi0;
	double *psi;
	double *x;
	double ton;
	double toff;
	double *timet;
	double *dipole;
	int gauge;
	int transformgauge;
	double x_int;
	struct analy_def analy;	
	int InterpByDTorNT, Ntinterp, PrintGaborAndSpectrum, PrintOutputMethod;	
	double textend, dtGabor, tmin1window, tmin2window, tmax1window, tmax2window, a_Gabor, omegaMaxGabor;
	struct output_print_def Print;
	double CV;
	char precision[2];
} inputs_def;

typedef struct outputs_def{ // only * can be modified by direct inputs
	double *tgrid;
	double *tgrid_fftw;
	double *Efield;
	double *sourceterm;
	double *omegagrid;
	double **FEfield, *FEfield_data;
	double **Fsourceterm, *Fsourceterm_data;
	double *FEfieldM2;
	double *FsourcetermM2;
	double *PopTot;
	double *sourcetermfiltered;
	double *PopInt;
	double *expval;
	int Nt; 
	int Nomega;
/*	double *tmax;*/
} outputs_def;

// functions operating on them
output_print_def Initialise_Printing_struct(void);
void outputs_destructor(outputs_def *);
void inputs_destructor(inputs_def *);
output_print_def Initialise_Printing_struct(void);
output_print_def Set_all_prints(void);

#endif