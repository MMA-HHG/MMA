/**
 * @file tools_fftw3.c
 * @brief Contains functions employing spectral methods using FFTW3 library.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "constants.h"
#include "tools_fftw3.h"

/**
 * @brief Takes signal and interpolates it using Fourier transform. 
 * 
 * @details The function takes the FFT of the input signal, extends the signal
 * to higher frequencies by adding zeros in the end of the signal and then 
 * creates the interpolated signal by IFFT and normalization.
 * 
 * @param k Number of points to add between the two consecutive points in the original signal. 
 * @param signal Signal for the interpolation.
 * @param N Number of points in the original signal.
 * @return Interpolated signal.
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
 * @brief Computes FFT of a signal and stores the FFTs.
 * 
 * @warning Our convention of FT is complex conjugate wrt FFTW3 DFT.
 * 
 * @param N Sample points.
 * @param dx Spacing of the samples.
 * @param xmax Maximum size of the sample grid.
 * @param signal Signal for the transform.
 * @param xigrid Frequencies.
 * @param fsig Positive part of the Fourier spectrum.
 * @param fsigM2 Absolute value of the spectrum.
 * @param Nxi Number of spectral points.
 */
void calcFFTW3(int N, double dx, double xmax, double *signal, double **xigrid,
			   double **fsig, double **fsigM2, int *Nxi) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	// Number of real spectral points
	int Nc;
	// Output spectrum array for FFTW
	fftw_complex *out;
	// Complex input array for storing the signal
	double *in;
	// Frequency stepsize
	double dxi;
	// FFT prefactors
	double coeff1, coeff2;
	// FFTW plan
	fftw_plan p;
	// Iterable
	int k1;

	// Half of the total number of points in the signal
	Nc = floor(((double)N) / 2.); 
	Nc++; 

	// Apply FFT prefactors according to the Parseval's theorem
	dxi = 2.*Pi/xmax;
	coeff1 = dx/ sqrt(2.*Pi); 
	coeff2 = dx*dx/(2.*Pi);
	
	// Allocate arrays
	*xigrid = (double*) calloc(Nc,sizeof(double));
	*fsig = (double*) calloc(2*Nc,sizeof(double));
	*fsigM2 = (double*) calloc(Nc,sizeof(double));
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	in = calloc(2*Nc,sizeof(double));
	// Store input signal into reduced array
	for(k1 = 0; k1 <= (N-1); k1++) {
		in[k1] = signal[k1];
	} 

	// Execute FFT
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE);
	fftw_execute(p);

	// write results
	for(k1 = 0; k1 <= (Nc-1); k1++) {
		(*xigrid)[k1] = ((double)k1)*dxi;
		(*fsig)[2*k1] = coeff1*out[k1][0]; 
		(*fsig)[2*k1+1] = - coeff1*out[k1][1]; 
		(*fsigM2)[k1] = coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]);
	}

	// Store number of frequency samples
	*Nxi = Nc;

	// Deallocate
	fftw_destroy_plan(p); 
	fftw_free(out); 
	free(in); 

	return;
}

/**
 * @brief Computes Gabor transform of signal.
 * 
 * @param signal Signal for computing Gabor.
 * @param dt Timestep of signal.
 * @param N Number of points in the signal.
 * @param N_freq Number of frequency points to include in the signal (relates to maximum omega).
 * @param N_G Number of temporal points for Gabor evaluation.
 * @param t_min Minimum time for Gabor evaluation.
 * @param t_max Maximum time for Gabor evaluation.
 * @param a Gabor parameter.
 * @return double** 
 */
double ** GaborTransform(double *signal, double dt, int N, int N_freq, int N_G, double t_min, double t_max, double a) 
{
	// Number of points in the FFT
	int Nc;
	// Output array for FFTW
	fftw_complex *out;
	// Input array for FFTW
	double *in;
	// Normalization coefficient
	double norm;
	// FFTW plan 
	fftw_plan p;
	// dt for Gabor
	double dt_G;
	// Gabor transform array
	double **gabor_transform;
	// Iterables
	int i, j;

	// Set variables
	a = 1.0/a;
	dt_G = (t_max - t_min)/N_G;
	Nc = floor(((double)N)/2.); 
	Nc++;

	// Allocate arrays
	gabor_transform = malloc(sizeof(double *) * N_G);
	for (j = 0; j < N_G; j++) {
		gabor_transform[j] = calloc(N_freq, sizeof(double));
	}
	in = calloc(2*Nc, sizeof(double));	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	// DFT normalization coefficient dt/sqrt(2 \pi)
	norm = dt/sqrt(2.*Pi);

	for (i = 0; i < N_G; i++) // gabor loop
	{

		for(j = 0; j < N; j++) {
			in[j]= exp(-pow(a * (((double)j)*dt - (t_min + ((double)i)*dt_G)), 2.)) * signal[j];
		} 
		
		// Plan FFTW - real to complex in 1D
		p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); 
		// Run FFTW
		fftw_execute(p); 

		for(j = 0; j < N_freq; j++) {
			gabor_transform[i][j] = norm * sqrt((out[j][0]*out[j][0] + out[j][1]*out[j][1]));
		}
	}

	return gabor_transform;

}
