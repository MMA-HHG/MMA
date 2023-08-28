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
 * @brief Full specification of the input field.
 */
typedef struct Efield_var {
	// Type of field – numerical/analytical.
	int fieldtype,fieldtype2;
	// Time grid
	double *tgrid;
	// Field array
	double *Field;
	// Timestep
	double dt;
	// Number of timesteps
	int Nt;
	// Trapezoidal envelope
	struct trap_definition trap;
	// Sine-squared envelope of vector potential
	struct sin2_definition *sin2;
	// Sine-squared envelope of electric field
	struct Esin2_definition *Esin2;
	// Flat-top envelope
	struct flattop1_definition *flt1;
	// Flat-top envelope with chirp
	struct flattop1chirp_definition *flt1ch;
	// Number of points in one cycle
	int Nflt1,Nsin2,NEsin2,Nflt1ch;
	// Field parameters
	double omega,E0,phi,ton,toff;
	// Number of points in the envelope
	int nc;
} Efield_var;

/**
 * @brief Microscopic target specification
 * 
 */
typedef struct trg_def{
	// Parameter of the soft core potential
	double a;
} trg_def;

/**
 * @brief Auxiliary structures for drivers of I/O
 * 
 */
typedef struct analy_def{
	// time spacing for writing the wavefunction
	double tprint;
	// writewavefunction (1-writting every tprint)
	int writewft;
} analy_def;

/**
 * @brief Auxiliary structures for drivers of I/O
 * 
 */
typedef struct output_print_def{
	int Efield, FEfield, sourceterm, Fsourceterm, FEfieldM2, FsourceTermM2, PopTot, tgrid, omegagrid, PopInt, expval_x;
} output_print_def;


/**
 * @brief Input structure
 * 
 * @details Full specification of the input parameters for the TDSE computation.
 * 
 */
typedef struct inputs_def {
	// Target info
	struct trg_def trg;
	// Field info
	struct Efield_var Efield;
	// Energy of the ground state
	double Eguess, Einit;
	// Start time
	double tmin;
	// Timesteps for TDSE
	int Nt;
	// Timestep grid size CUPRAD
	int num_t;
	// Timestep
	double dt;
	// Radial grid size for CUPRAD
	int num_r;
	// Number of points of the spatial grid for the expansion
	int num_exp;
	// Spatial step
	double dx;
	// Initial wavefunction
	double *psi0;
	// Wavefunction in time t
	double *psi;
	// Spatial grid for the TDSE
	double *x;
	// Start of the envelope
	double ton;
	// End of the envelope
	double toff;
	// Time grid TDSE
	double *timet;
	// TDSE dipole
	double *dipole;
	// Choice of gauge
	int gauge;
	// Gauge transform
	int transformgauge;
	// Integration limit for the ionization computation, note: 2 works fine with the lenth gauge and strong fields
	double x_int;
	// Specifies output dataset
	struct analy_def analy;	
	// Switches (1 = yes)
	int InterpByDTorNT, Ntinterp, PrintGaborAndSpectrum;
	// Switch (0 - only text, 1 - only binaries, 2 - both)
	int PrintOutputMethod;	
	// extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD
	double textend;
	// Gabor variables
	double dtGabor, tmin1window, tmin2window, tmax1window, tmax2window, a_Gabor, omegaMaxGabor;
	// I/O printing structure
	struct output_print_def Print;
	// Precision of the ground state energy
	double CV;
	char precision[2];
} inputs_def;

typedef struct outputs_def{ 
	// Temporal grid
	double *tgrid;
	// Temporal grid for FFTW
	double *tgrid_fftw;
	// Electric field
	double *Efield;
	// Source term (current)
	double *sourceterm;
	// Angular frequency grid
	double *omegagrid;
	// Field spectrum
	double **FEfield, *FEfield_data;
	// Source term spectrum
	double **Fsourceterm, *Fsourceterm_data;
	double *FEfieldM2;
	double *FsourcetermM2;
	// Population of the ground state
	double *PopTot;
	// Filtered source term
	double *sourcetermfiltered;
	// Ionization probability
	double *PopInt;
	// Expectation value of x
	double *expval;
	// Number of temporal gridpoints
	int Nt; 
	// Number of frequency gridpoints
	int Nomega;
} outputs_def;

// functions operating on them
output_print_def Initialise_Printing_struct(void);
void outputs_destructor(outputs_def *);
void inputs_destructor(inputs_def *);
output_print_def Initialise_Printing_struct(void);
output_print_def Set_all_prints(void);

#endif