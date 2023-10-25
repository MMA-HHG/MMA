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
	double oc;
	double o;
	double A0;
	double nc;
	double phi;
	double phi0;
	double ti;
	double E0;
} sin2_definition;

typedef struct Esin2_definition {
	double oc;
 	double o;
	double A0;
	double nc;
	double phi;
	double phi0;
	double ti;
	double E0;
} Esin2_definition;

typedef struct trap_definition {
	double omega;
	double E0;
	double phi;
	double ton;
	double toff;
	int nc;
} trap_definition;

typedef struct flattop1_definition {
	double ton;
	double toff;
	double T;
	double o;
	double phi;
	double A;
	double ti;
} flattop1_definition;

typedef struct flattop1chirp_definition {
	double ton;
	double toff;
	double T;
	double o;
	double phi;
	double A;
	double ti;
	double b;
	double c;
} flattop1chirp_definition;

/**
 * @brief Full specification of the input field.
 */
typedef struct Efield_var {
	// Type of field – numerical/analytical.
	int fieldtype;
	int fieldtype2;
	// Time grid
	double *tgrid;
	// Field array
	double *Field;
	// Timestep
	double dt;
	// Field parameters
	double omega;
	double E0;
	double phi;
	double ton;
	double toff;
	// Number of timesteps
	int Nt;
	// Number of points in one cycle
	int nc;
	int Nftlt1;
	int N;
	int NEsin2;
	int Nftlt1ch;
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
	// Number of points in the envelope
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
	int Efield;
	int FEfield;
	int sourceterm;
	int Fsourceterm;
	int FEfieldM2;
	int FsourceTermM2;
	int PopTot;
	int tgrid;
	int omegagrid;
	int PopInt;
	int expval_x;
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
	double Eguess;
	double Einit;
	// Start time
	double tmin;
	// Timesteps for TDSE
	int Nt;
	// Length of one cycle for 800 nm
	int num_t;
	// Timestep
	double dt;
	// Spatial grid size for TDSE
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
	int InterpByDTorNT;
	int Ntinterp; 
	// Switch (0 - only text, 1 - only binaries, 2 - both)
	int PrintOutputMethod;	
	// extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD
	double textend;
	// I/O printing structure
	struct output_print_def Print;
	// Precision of the ground state energy
	double CV;
	char precision[2];
} inputs_def;

typedef struct outputs_def{ 
	// Temporal grid
	double *tgrid;
	// Electric field
	double *Efield;
	// Source term (current)
	double *sourceterm;
	// Angular frequency grid
	double *omegagrid;
	// Field spectrum
	double **FEfield;
	double *FEfield_data;
	// Source term spectrum
	double **Fsourceterm;
	double *Fsourceterm_data;
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
	// Wavefunction in time t
	double **psi;
} outputs_def;

// functions operating on them
output_print_def Initialise_Printing_struct(void);
void outputs_destructor(outputs_def *);
void inputs_destructor(inputs_def *);
void free_mtrx(double ** buf, int N_rows);
output_print_def Initialise_Printing_struct(void);
output_print_def Set_all_prints(void);

#endif