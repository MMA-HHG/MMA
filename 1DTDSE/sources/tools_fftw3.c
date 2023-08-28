/**
 * @file tools_fftw3.c
 * @brief Contains functions employing spectral methods using FFTW3 library.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <math.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <fftw3.h>
#include "constants.h"
#include "tools_fftw3.h"


/**
 * @brief Takes signal and interpolates it using Fourier transform. 
 * 
 * @details The function takes the FFT of the input signal, extends the signal
 * to higher frequencies by adding zeros in the end of the signal and then 
 * creates the interpolated signal by IFFT and normalization.
 * 
 * @param k (int) Number of points to add between the two consecutive points in the original signal. 
 * @param signal (double*) Signal for the interpolation.
 * @param N (int) Number of points in the original signal.
 * @return double* Interpolated signal.
 */
double* FourInterp(int k, double *signal, int N)
{
	int Nc, N2, Nc2;
	fftw_complex *out, *in2;
	double *in, *out2;
	fftw_plan p;
	int k1;


	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal[k1];} // !!! REDUNDANT

	
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	fftw_destroy_plan(p); // deallocate plan memory
		
	// INVERSE + PADDING
	N2 = k*N;

	Nc2 = floor(((double)N2) / 2.); Nc2++;
	
	in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc2);
	out2 =  calloc(2*Nc2,sizeof(double));
	
	
	// zeroes
	for(k1 = 0; k1 <= (Nc2-1); k1++){in2[k1][0]=0.; in2[k1][1]=0.;}// zeroes initialised
	for(k1 = 0; k1 <= (Nc-1); k1++){in2[k1][0]=out[k1][0]/N; in2[k1][1]=out[k1][1]/N;} // previous data + normalisation	

	p = fftw_plan_dft_c2r_1d(N2, in2, out2, FFTW_ESTIMATE); // plan iFFTW
	fftw_execute(p); // run iFFTW

	fftw_destroy_plan(p); // deallocate plan memory

	
	free(in); fftw_free(in2); fftw_free(out);
	return out2;

}


/**
 * @brief Prints FFT of a signal.
 * 
 * @details Takes real signal spaced by given "dt" and it computes and prints its FFT.
 * 
 * @param sig 
 * @param fsig 
 * @param signal 
 * @param N 
 * @param dx 
 */
void printFFTW3(FILE *sig, FILE *fsig, double *signal, int N, double dx) //
{
	int Nc;
	fftw_complex *out;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signal
	for(k1 = 0; k1 <= (N-1); k1++){
					fprintf(sig,"%e\t%e\n", ((double)k1)*dx , signal[k1]);
					}


	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW


	dxi = 2.*Pi/(  ((double)N) * dx);

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out[k1][0], -coeff1*out[k1][1] , coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));
					}

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}

void print2FFTW3(FILE *sig, FILE *fsig, double *signal1, double *signal2, int N, double dx, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signal
	for(k1 = 0; k1 <= (N-1); k1++){
					fprintf(sig,"%e\t%e\t%e\n", ((double)k1)*dx , signal1[k1], signal2[k1]);
					}

	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal1[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out1, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT

	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW
	
	dxi = 2.*Pi/xmax; 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out1[k1][0], -coeff1*out1[k1][1] , coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]),coeff1*out2[k1][0], -coeff1*out2[k1][1] , coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]));
					}
	
	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


// PRINT Gabor of a signal
void printGaborFFTW3(FILE *Gsig, FILE *xgrid, FILE *xigrid, FILE *Gsigbin, double *signal, int N, double dx, double dxG, double a, double xiMaxPrint) // takes real signal speced by given "dx" and it computes and prints its Gabor transform, The parameters of the Gabor transform are new "dxG" (will be adjusted to a close one matching the points) and gabor parameter "a"
{
	int Nc, Ncprint;
	fftw_complex *out;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1,k2,kstep2;
	
	double *dumptr;

	a = 1.0/a;

	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);


	dxi = 2.*Pi/(  ((double)N) * dx); // dimension-preserving coefficients
	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);

	Ncprint = floor( xiMaxPrint/dxi ); Ncprint++; if (Ncprint > (Nc-1) ){Ncprint = Nc-1;} // maximum frequency

	
	if (dxG > dx){kstep2=floor(dxG/dx);}else{kstep2=1;}


	for(k2 = 0; k2 <= (N-1); k2 = k2 + kstep2) // gabor loop
	{

		for(k1 = 0; k1 <= (N-1); k1++){in[k1]= exp( -pow( a*dx*( ((double)k1)-((double)k2) ),2.)  ) * signal[k1];} // 
		
		p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
		fftw_execute(p); // run FFTW

		dumptr = calloc(1,sizeof(double));

		for(k1 = 0; k1 <= (Ncprint-1); k1++){
						dumptr[0] = sqrt( coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));
						fprintf(Gsig,"%e\t", dumptr[0]) ;

						fwrite( dumptr ,sizeof(double),1,Gsigbin);

						}
		dumptr[0] = sqrt( coeff2*(out[Ncprint][0]*out[Ncprint][0]+out[Ncprint][1]*out[Ncprint][1]));
		fprintf(Gsig,"%e\n", dumptr[0] ) ;
		fwrite( dumptr,sizeof(double),1,Gsigbin);

		

		fprintf(xgrid,"%e\n", ((double)k2)*dx );
	}
	for(k1 = 0; k1 <= Ncprint; k1++){fprintf(xigrid,"%e\n", ((double)k1)*dxi) ;}

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}

void print2FFTW3binary(FILE *xgrid, FILE *sig1, FILE *sig2, FILE *xigrid, FILE *fsig1, FILE *fsig2, FILE *fsig1M2, FILE *fsig2M2, FILE *GridDimensions, double *signal1, double *signal2, int N, double dx, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2;
	double *in, *dum;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signals
	fwrite(signal1,sizeof(double),N,sig1);
	fwrite(signal2,sizeof(double),N,sig2);
	dum = calloc(1,sizeof(double));
	for(k1 = 0; k1 <= (N-1); k1++){dum[0]=((double)k1)*dx; fwrite(dum,sizeof(double),1,xgrid);}

	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal1[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out1, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT

	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW
	
	
	dxi = 2.*Pi/xmax; 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);

	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=((double)k1)*dxi; fwrite(dum,sizeof(double),1,xigrid);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff1*out1[k1][0]; fwrite(dum,sizeof(double),1,fsig1);dum[0]=-coeff1*out1[k1][1]; fwrite(dum,sizeof(double),1,fsig1);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff1*out2[k1][0]; fwrite(dum,sizeof(double),1,fsig2);dum[0]=-coeff1*out2[k1][1]; fwrite(dum,sizeof(double),1,fsig2);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]); fwrite(dum,sizeof(double),1,fsig1M2);}	
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]); fwrite(dum,sizeof(double),1,fsig2M2);}

	fprintf(GridDimensions,"%i\n", N);
	fprintf(GridDimensions,"%i\n", Nc);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


// PRINT limited range FFTW3 of a signal (rectangular window), the result IS NOT PHASE SHIFTED !!!
void printlimitedFFTW3(FILE *fsig, double *signal, int N, double dx, double xmin, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc, NF;
	fftw_complex *out;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1,kmin,kmax;



	kmin = floor( xmin/dx ); // maximum frequency
	kmax = floor( xmax/dx ); kmax++; if (kmax > N-1){kmax = N;} // maximum frequency
	NF = kmax - kmin;


	Nc = floor(((double)NF) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (NF-1); k1++){in[k1]=signal[k1+kmin];} // !!! NEEDED
	
	p = fftw_plan_dft_r2c_1d(NF, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	dxi = 2.*Pi/(  ((double)NF) * dx); 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out[k1][0], -coeff1*out[k1][1] , coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));
					}
	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


void printGaborFFTW3binary(FILE *Gsize, FILE *xgrid, FILE *xigrid, FILE *Gsigbin, double *signal, int N, double dx, double dxG, double a, double xiMaxPrint) // takes real signal speced by given "dx" and it computes and prints its Gabor transform, The parameters of the Gabor transform are new "dxG" (will be adjusted to a close one matching the points) and gabor parameter "a"
{
	int Nc, Ncprint;
	fftw_complex *out;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1,k2,kstep2,NxG;
	
	double *dumptr1, *dumptr2;

	a = 1.0/a;

	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);


	dxi = 2.*Pi/(  ((double)N) * dx); // dimension-preserving coefficients
	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);

	Ncprint = floor( xiMaxPrint/dxi ); Ncprint++; if (Ncprint > (Nc-1) ){Ncprint = Nc-1;} // maximum frequency
	fprintf(Gsize,"%i\n", Ncprint+1);

	
	if (dxG > dx){kstep2=floor(dxG/dx);}else{kstep2=1;}
	

	NxG = 0; // counting the length of the loop
	for(k2 = 0; k2 <= (N-1); k2 = k2 + kstep2) // gabor loop
	{
		NxG++;
		for(k1 = 0; k1 <= (N-1); k1++){in[k1]= exp( -pow( a*dx*( ((double)k1)-((double)k2) ),2.)  ) * signal[k1];} // 
 
		dumptr1 = calloc( (Ncprint+1) ,sizeof(double));
		dumptr2 = calloc(1,sizeof(double));
		p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
		fftw_execute(p); // run FFTW

		for(k1 = 0; k1 <= Ncprint; k1++){dumptr1[k1] = sqrt( coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));}
		fwrite( dumptr1 ,sizeof(double), (Ncprint+1) ,Gsigbin);
		free(dumptr1);

		dumptr2[0] = ((double)k2)*dx;
		fwrite(dumptr2 ,sizeof(double),1,xgrid);
		free(dumptr2);
	}
	fprintf(Gsize,"%i\n", NxG);

	dumptr1 = calloc( (Ncprint+1) ,sizeof(double));
	for(k1 = 0; k1 <= Ncprint; k1++){dumptr1[k1] = ((double)k1)*dxi;}
	fwrite( dumptr1 ,sizeof(double), (Ncprint+1) ,xigrid);
	free(dumptr1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft
	
	return;

}

void calc2FFTW3(int N, double dx, double xmax, double *signal1, double *signal2, double **xgrid, double **xigrid, double **fsig1, double **fsig2, double **fsig1M2, double **fsig2M2, int *Nxi) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// DO THE TRANSFORMS
	
	Nc = floor(((double)N) / 2.); Nc++; // # of points
	in = calloc(2*Nc,sizeof(double));
	
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); 
	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal1[k1];} // !!! REDUNDANT	
	p = fftw_plan_dft_r2c_1d(N, in, out1, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p);
	fftw_destroy_plan(p); // deallocate plan memory

	out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT
	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p);
	fftw_destroy_plan(p); // deallocate plan memory


	//	RESCALE TO OUTPUTS
	dxi = 2.*Pi/xmax;
	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);

	*xgrid = (double*) calloc(N,sizeof(double));
	*xigrid = (double*) calloc(Nc,sizeof(double));
	*fsig1 = (double*) calloc(2*Nc,sizeof(double));
	*fsig2 = (double*) calloc(2*Nc,sizeof(double));
	*fsig1M2 = (double*) calloc(Nc,sizeof(double));
	*fsig2M2 = (double*) calloc(Nc,sizeof(double));

	// write results
	for(k1 = 0; k1 <= (N-1); k1++){(*xgrid)[k1]=((double)k1)*dx;}
	for(k1 = 0; k1 <= (Nc-1); k1++){
		(*xigrid)[k1] = ((double)k1)*dxi;
		(*fsig1)[2*k1] = coeff1*out1[k1][0]; (*fsig1)[2*k1+1] = - coeff1*out1[k1][1]; // !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft
		(*fsig2)[2*k1] = coeff1*out2[k1][0]; (*fsig2)[2*k1+1] = - coeff1*out2[k1][1];
		(*fsig1M2)[k1] = coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]);
		(*fsig2M2)[k1] = coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]);
	}

	fftw_free(out1); fftw_free(out2); free(in);
	*Nxi = Nc;
	return;

}



void calcFFTW3(int N, double dx, double xmax, double *signal, double **xgrid, double **xigrid, double **fsig, double **fsigM2, int *Nxi) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// DO THE TRANSFORMS
	
	Nc = floor(((double)N) / 2.); Nc++; // # of points
	in = calloc(2*Nc,sizeof(double));
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); 
	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal[k1];} // !!! REDUNDANT	
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p);

	//	RESCALE TO OUTPUTS
	dxi = 2.*Pi/xmax;
	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);


	// double *ptr1, *ptr2;
	*xgrid = (double*) calloc(N,sizeof(double));
	*xigrid = (double*) calloc(Nc,sizeof(double));
	*fsig = (double*) calloc(2*Nc,sizeof(double));
	*fsigM2 = (double*) calloc(Nc,sizeof(double));

	// write results
	for(k1 = 0; k1 <= (N-1); k1++){(*xgrid)[k1]=((double)k1)*dx;}
	for(k1 = 0; k1 <= (Nc-1); k1++){
		(*xigrid)[k1] = ((double)k1)*dxi;
		(*fsig)[2*k1] = coeff1*out[k1][0]; (*fsig)[2*k1+1] = - coeff1*out[k1][1]; // !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft
		(*fsigM2)[k1] = coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]);
	}

	fftw_destroy_plan(p); // deallocate plan memory
	fftw_free(out); free(in); 
	*Nxi = Nc;
	return;

}
