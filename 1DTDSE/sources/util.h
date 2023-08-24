#include "structures.h"


// This file contains the variables used by various part of the code.
// It should not contain any external constructions (fftw3,hdf5,mpi,...), so it can be used in any subtask using none of these.


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// The field and atomic target:



// functions
double potential(double, trg_def);
double gradpot(double, trg_def);
double norme(double *,int);
void normalise(double *,int);
double Afieldflattop1(double, double, double , double , double , double , double , double );
double Afieldflattop1ch(double, double, double , double , double , double , double , double , double , double);
double smootherstep(double , double, double);
double clamp(double, double, double);

double Afieldsin2(double, double, double, double, double, double, double);
double dAfieldsin2(double, double, double, double, double, double, double);

double Primsin2cos(double , double , double , double , double );
double AfieldEsin2(double , double , double , double , double , double , double );

double EField(double,double,double,double,int,double,double);
double AField( Efield_var,double);
double dAField( Efield_var,double);

void define_analytical( Efield_var *, FILE *); // obsolete


void Initialise_grid_and_D2(double, int, double **, double **, double **);
void Initialise_grid_and_ground_state(inputs_def *);


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// FFTW3 drivers, that do not tneed fftw3 header


////////////////////////////////////////////////////////////////////////////////////////////////////////////
// Algorithmical tools
void coarsen_grid_real(double *, int, double **, int *, int, int);

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



double Einitialise( trg_def, double *,double *,double *,double *,double *,double,double ,int);
double E_calculation(double *,double *,double *,double *,int);
double E_calculation_numerov( trg_def, double *,double ,double *,int );





void window_analysis( trg_def,double,double,double,int,int,double,double*,double*,double*,double*,double*);
void dipole_analysis(double,double,double*,double*,int,int);
double* extend_grid(double *,int,int,int);


void printresults( trg_def,  Efield_var,FILE *,int,double *,int, double *, double, double *,double, double, double, double, double,  outputs_def);




void vander(double *, double *, double *, int);




double ** create_2Darray_accessor_real(int *, double *);



