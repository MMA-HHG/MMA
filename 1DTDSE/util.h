// This file contains the variables used by various part of the code.
// It should not contain any external constructions (fftw3,hdf5,mpi,...), so it can be used in any subtask using none of these.


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The field and atomic target:

// various analytic envelopes (not reintroduced yet)
struct sin2_definition {
	double oc, o, A0, nc, phi, phi0, ti, E0;
};
struct Esin2_definition {
	double oc, o, A0, nc, phi, phi0, ti, E0;
};
struct trap_definition {
	double omega,E0,phi,ton,toff;
	int nc;
};
struct flattop1_definition {
	double ton,toff,T,o,phi,A,ti;
};
struct flattop1chirp_definition {
	double ton,toff,T,o,phi,A,ti,b,c;
};

// full specification of the input field & functions operating on them
struct Efield_var {
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
};

// microscopic target:
struct trg_def{
	double a;
};

// functions
double potential(double, struct trg_def);
double gradpot(double, struct trg_def);
double Afieldflattop1(double, double, double , double , double , double , double , double );
double Afieldflattop1ch(double, double, double , double , double , double , double , double , double , double);
double smootherstep(double , double, double);
double clamp(double, double, double);

double Afieldsin2(double, double, double, double, double, double, double);
double dAfieldsin2(double, double, double, double, double, double, double);

double Primsin2cos(double , double , double , double , double );
double AfieldEsin2(double , double , double , double , double , double , double );

double EField(double,double,double,double,int,double,double);
double AField(struct Efield_var,double);
double dAField(struct Efield_var,double);

void define_analytical(struct Efield_var *, FILE *); // obsolete


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Auxiliary structures for drivers of I/O
struct analy_def{
	double tprint;
	int writewft;
};
struct output_print_def{
	int Efield, FEfield, sourceterm, Fsourceterm, FEfieldM2, FsourceTermM2, PopTot, tgrid, omegagrid;
};



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// I/O structures & functions operating on them
struct inputs_def{
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
};
struct outputs_def{ // only * can be modified by direct inputs
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
};

// functions operating on them
struct output_print_def Initialise_Printing_struct(void);
void outputs_destructor(struct outputs_def *);
void inputs_destructor(struct inputs_def *);
struct output_print_def Initialise_Printing_struct(void);
struct output_print_def Set_all_prints(void);
void Initialise_grid_and_ground_state(struct inputs_def *);
void Initialise_grid_and_D2(double, int, double **, double **, double **);
struct outputs_def call1DTDSE(struct inputs_def); // the caller



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The core fuctions of the code
double* propagation(struct trg_def, struct Efield_var,double,int,int,double,int,int,double,double *,double *,double *
				 ,FILE *,FILE *, double,double,double*,double*,int,int,double, struct analy_def, struct outputs_def);
void compute_population(struct trg_def, struct Efield_var,int,double *,int, double *, double, double *,double, double, double, double, double, struct outputs_def);



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FFTW3 drivers, that do not tneed fftw3 header
void printFFTW3(FILE *, FILE *, double *, int, double );
void print2FFTW3(FILE *, FILE *, double *, double *, int, double, double);
void print2FFTW3binary(FILE *, FILE *, FILE *, FILE *,FILE *, FILE *,FILE *, FILE *, FILE *, double *, double *, int, double, double);
void printGaborFFTW3(FILE *, FILE *, FILE *, FILE *, double *, int, double, double, double, double);
void printlimitedFFTW3(FILE *, double *, int, double, double, double);
void printGaborFFTW3binary(FILE *, FILE *, FILE *, FILE *, double *, int, double, double, double, double);

void calc2FFTW3(int, double, double, double *, double *, double **, double **, double **, double **, double **, double **, int *);
void calcFFTW3(int, double, double, double *, double **, double **, double **, double **, int *);



////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Algorithmical tools
void coarsen_grid_real(double *, int, double **, int *, int, int);


double* FourInterp(int , double * , int );
double interpolate( int , double , double* , double* );
double interpolate2n( int , double , double* , double* , int );
void findinterval(int , double , double* , int* , int* );
double findnextinterpolatedzero(int, double, double* , double* );

void nxtval_init(int, int *);
void nxtval_strided(int, int *);






////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Here be dragons

void tqli(double *,double *, int n,double *);
double Calculate_Shift(double, double, double); 
void Transform_Matrix(double*, double, double, int, int);
int QL_Tridiagonal_Symmetric_Matrix( double *, double *,double *, int, int);

void Inv_Tridiagonal_Matrix( double *, double *, double *, double *, double *, int);
void Inv_Tridiagonal_Matrix_complex( double *, double *, double *, double *, double *, int );
void Inv_Tridiagonal_Matrix_complex_Numerov( double *, double *, double *, double *, double*, int );


double norme(double *,int);
void normalise(double *,int);
double Einitialise(struct trg_def, double *,double *,double *,double *,double *,double,double ,int);
double E_calculation(double *,double *,double *,double *,int);
double E_calculation_numerov(struct trg_def, double *,double ,double *,int );





void window_analysis(struct trg_def,double,double,double,int,int,double,double*,double*,double*,double*,double*);
void dipole_analysis(double,double,double*,double*,int,int);
double* extend_grid(double *,int,int,int);


void printresults(struct trg_def, struct Efield_var,FILE *,int,double *,int, double *, double, double *,double, double, double, double, double, struct outputs_def);




void vander(double *, double *, double *, int);




double ** create_2Darray_accessor_real(int *, double *);



