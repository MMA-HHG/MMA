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

struct trg_def{
	double a;
};

struct analy_def{
	double tprint;
	int writewft;
};

struct outputs_def{
	double *tgrid;
	double *Efield;
	double *sourceterm;
	double *omegagrid;
	double *FEfield;
	double *Fsourceterm;
	double *PopTot;
	double *sourcetermfiltered;
/*	double *tmax;*/
};

void tqli(double *,double *, int n,double *);
double Calculate_Shift(double, double, double); 
void Transform_Matrix(double*, double, double, int, int);
int QL_Tridiagonal_Symmetric_Matrix( double *, double *,double *, int, int);
void Inv_Tridiagonal_Matrix( double *, double *, double *, double *, double *, int);
void Inv_Tridiagonal_Matrix_complex( double *, double *, double *, double *, double *, int );
void Inv_Tridiagonal_Matrix_complex_Numerov( double *, double *, double *, double *, double*, int );

double potential(double, struct trg_def);
double gradpot(double, struct trg_def);
double norme(double *,int);
void normalise(double *,int);
void Initialise(int);
double Einitialise(struct trg_def, double *,double *,double *,double *,double *,double,double ,int);
double E_calculation(double *,double *,double *,double *,int);
double E_calculation_numerov(struct trg_def, double *,double ,double *,int );
double EField(double,double,double,double,int,double,double);
double AField(struct Efield_var,double);
double dAField(struct Efield_var,double);
void gaussian(void);
void volkov_state(void);
void volkov_state_vg(void);


double* propagation(struct trg_def, struct Efield_var,double,int,int,double,int,int,double,double *,double *,double *
				 ,FILE *,FILE *, double,double,double*,double*,int,int,double, struct analy_def, struct outputs_def);
void window_analysis(struct trg_def,double,double,double,int,int,double,double*,double*,double*,double*,double*);
void dipole_analysis(double,double,double*,double*,int,int);
double* extend_grid(double *,int,int,int);
double* rmv_gs(double*,double*,double*,double);
double interpolate( int , double , double* , double* );
double interpolate2n( int , double , double* , double* , int );
void findinterval(int , double , double* , int* , int* );
double findnextinterpolatedzero(int, double, double* , double* );

void printresults(struct trg_def, struct Efield_var,FILE *,int,double *,int, double *, double, double *,double, double, double, double, double, struct outputs_def);
double Afieldflattop1(double, double, double , double , double , double , double , double );
double Afieldflattop1ch(double, double, double , double , double , double , double , double , double , double);
double smootherstep(double , double, double);
double clamp(double, double, double);

double Afieldsin2(double, double, double, double, double, double, double);
double dAfieldsin2(double, double, double, double, double, double, double);

double Primsin2cos(double , double , double , double , double );
double AfieldEsin2(double , double , double , double , double , double , double );


void vander(double *, double *, double *, int);

double* FourInterp(int , double * , int );
void printFFTW3(FILE *, FILE *, double *, int, double );
void print2FFTW3(FILE *, FILE *, double *, double *, int, double, double);
void print2FFTW3binary(FILE *, FILE *, FILE *, FILE *,FILE *, FILE *,FILE *, FILE *, FILE *, double *, double *, int, double, double);
void printGaborFFTW3(FILE *, FILE *, FILE *, FILE *, double *, int, double, double, double, double);
void printlimitedFFTW3(FILE *, double *, int, double, double, double);

void define_analytical(struct Efield_var *, FILE *);


// HDF5
int MPE_Counter_create(MPI_Comm, int, MPI_Win *);
int MPE_Counter_nxtval(MPI_Win, int, int *, int);
int MPE_Mutex_acquire(MPI_Win, int, int);
int MPE_Mutex_release(MPI_Win, int, int);


void addone(int *); // to test pointers
double readreal(hid_t, char *, herr_t *, double *);






















