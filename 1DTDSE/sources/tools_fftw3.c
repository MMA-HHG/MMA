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
#include "numerical_constants.h"
//#include "util.h"
#include "tools_fftw3.h"




// FOURIER INTERPOLATION
double* FourInterp(int k, double *signal, int N)
{
	FILE *file1, *file2, *newxgrid, *newygrid, *newgrid;
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

	/*
	// WRITE RESULT	
	file1 = fopen("ftransform.dat" , "w");	
	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(file1,"%e\t%e\n",out[k1][0],out[k1][1]);
					} // !!!!!!!!!!! WHY POINTERS THING DOESN'T WORK?!
	fclose(file1);
	*/

		
	// INVERSE + PADDING
	N2 = k*N;

	Nc2 = floor(((double)N2) / 2.); Nc2++;
	
	in2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc2);
	out2 =  calloc(2*Nc2,sizeof(double));
	
	
	// zeroes
	for(k1 = 0; k1 <= (Nc2-1); k1++){in2[k1][0]=0.; in2[k1][1]=0.;}// zeroes initialised


	//for(k1 = 0; k1 <= (1077-1); k1++){in2[k1][0]=out[k1][0]/N; in2[k1][1]=out[k1][1]/N;} // previous data + normalisation	!!! TEST FOR NOISE 0-PADDING
	for(k1 = 0; k1 <= (Nc-1); k1++){in2[k1][0]=out[k1][0]/N; in2[k1][1]=out[k1][1]/N;} // previous data + normalisation	

	p = fftw_plan_dft_c2r_1d(N2, in2, out2, FFTW_ESTIMATE); // plan iFFTW
	fftw_execute(p); // run iFFTW

	fftw_destroy_plan(p); // deallocate plan memory

	/*
	// WRITE THE RESULTS
	Nxnew = N2;
	// dx = (xgrid[Nx-1]-xgrid[0])/((double)Nxnew);

	// newxgrid = fopen( "newxgrid.dat", "w" );
	
	// newgrid = fopen( "newgrid.dat", "w" );

	// x = xgrid[7990]+dx;
	// printf("x,  %lf \n",x);

	
	printf("ftest1 \n");
	newygrid = fopen( "newygrid.dat", "w" );;
	for ( k1 = 0 ; k1 <= (N2-1); k1++)
	{
		fprintf(newygrid,"%e\n",out2[k1]);
	}

	fclose(newygrid);
	printf("writting done \n");
	*/
	
	free(in); fftw_free(in2); fftw_free(out);
	return out2;

}


// PRINT FFTW3 of a signal
void printFFTW3(FILE *sig, FILE *fsig, double *signal, int N, double dx) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out, *in2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signal
	for(k1 = 0; k1 <= (N-1); k1++){
					fprintf(sig,"%e\t%e\n", ((double)k1)*dx , signal[k1]);
					}
	//fclose(file1);


	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	
	// print fourier transform
	// file1 = fopen("ftransform.dat" , "w");

	dxi = 2.*Pi/(  ((double)N) * dx);

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out[k1][0], -coeff1*out[k1][1] , coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));
					}
	//fclose(file1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


void print2FFTW3(FILE *sig, FILE *fsig, double *signal1, double *signal2, int N, double dx, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2, *in2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signal
	for(k1 = 0; k1 <= (N-1); k1++){
					fprintf(sig,"%e\t%e\t%e\n", ((double)k1)*dx , signal1[k1], signal2[k1]);
					}
	//fclose(file1);


	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal1[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out1, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT

	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW
	
	// print fourier transform
	// file1 = fopen("ftransform.dat" , "w");

/*	dxi = 2.*Pi/(  ((double)N) * dx); */
	dxi = 2.*Pi/xmax; 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out1[k1][0], -coeff1*out1[k1][1] , coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]),coeff1*out2[k1][0], -coeff1*out2[k1][1] , coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]));
					}
	//fclose(file1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}






// PRINT Gabor of a signal
void printGaborFFTW3(FILE *Gsig, FILE *xgrid, FILE *xigrid, FILE *Gsigbin, double *signal, int N, double dx, double dxG, double a, double xiMaxPrint) // takes real signal speced by given "dx" and it computes and prints its Gabor transform, The parameters of the Gabor transform are new "dxG" (will be adjusted to a close one matching the points) and gabor parameter "a"
{
	int Nc, Ncprint;
	fftw_complex *out, *in2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1,k2,kstep2;
	
	double dum;
	double *dumptr;


/*	// print signal*/
/*	for(k1 = 0; k1 <= (N-1); k1++){*/
/*					fprintf(sig,"%e\t%e\n", ((double)k1)*dx , signal[k1]);*/
/*					}*/
/*	//fclose(file1);*/

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

		
		// print fourier transform
		// file1 = fopen("ftransform.dat" , "w");

	
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

		//fclose(file1);
	}
	for(k1 = 0; k1 <= Ncprint; k1++){fprintf(xigrid,"%e\n", ((double)k1)*dxi) ;}

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


void print2FFTW3binary(FILE *xgrid, FILE *sig1, FILE *sig2, FILE *xigrid, FILE *fsig1, FILE *fsig2, FILE *fsig1M2, FILE *fsig2M2, FILE *GridDimensions, double *signal1, double *signal2, int N, double dx, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2, *in2;
	double *in, *dum;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1;


	// print signals
	fwrite(signal1,sizeof(double),N,sig1);
	fwrite(signal2,sizeof(double),N,sig2);
	dum = calloc(1,sizeof(double));
	for(k1 = 0; k1 <= (N-1); k1++){dum[0]=((double)k1)*dx; fwrite(dum,sizeof(double),1,xgrid);}

/*	for(k1 = 0; k1 <= (N-1); k1++){*/
/*					fprintf(sig,"%e\t%e\t%e\n", ((double)k1)*dx , signal1[k1], signal2[k1]);*/
/*					}*/
/*	//fclose(file1);*/


	Nc = floor(((double)N) / 2.); Nc++;

	in = calloc(2*Nc,sizeof(double));	
	
	out1 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc); out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal1[k1];} // !!! REDUNDANT
	
	p = fftw_plan_dft_r2c_1d(N, in, out1, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW

	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT

	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p); // run FFTW
	
	// print fourier transform
	// file1 = fopen("ftransform.dat" , "w");

/*	dxi = 2.*Pi/(  ((double)N) * dx); */
	dxi = 2.*Pi/xmax; 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);

	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=((double)k1)*dxi; fwrite(dum,sizeof(double),1,xigrid);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff1*out1[k1][0]; fwrite(dum,sizeof(double),1,fsig1);dum[0]=-coeff1*out1[k1][1]; fwrite(dum,sizeof(double),1,fsig1);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff1*out2[k1][0]; fwrite(dum,sizeof(double),1,fsig2);dum[0]=-coeff1*out2[k1][1]; fwrite(dum,sizeof(double),1,fsig2);}
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]); fwrite(dum,sizeof(double),1,fsig1M2);}	
	for(k1 = 0; k1 <= (Nc-1); k1++){dum[0]=coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]); fwrite(dum,sizeof(double),1,fsig2M2);}

	fprintf(GridDimensions,"%i\n", N);
	fprintf(GridDimensions,"%i\n", Nc);

/*	for(k1 = 0; k1 <= (Nc-1); k1++){*/
/*					fprintf(fsig,"%e\t%e\t%e\t%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out1[k1][0], -coeff1*out1[k1][1] , coeff2*(out1[k1][0]*out1[k1][0]+out1[k1][1]*out1[k1][1]),coeff1*out2[k1][0], -coeff1*out2[k1][1] , coeff2*(out2[k1][0]*out2[k1][0]+out2[k1][1]*out2[k1][1]));*/
/*					}*/
	//fclose(file1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}




// PRINT limited range FFTW3 of a signal (rectangular window), the result IS NOT PHASE SHIFTED !!!
void printlimitedFFTW3(FILE *fsig, double *signal, int N, double dx, double xmin, double xmax) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc, NF;
	fftw_complex *out, *in2;
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

	
	// print fourier transform
	// file1 = fopen("ftransform.dat" , "w");

	dxi = 2.*Pi/(  ((double)NF) * dx); 

	coeff1 = dx/ sqrt(2.*Pi); coeff2 = dx*dx/(2.*Pi);	

	for(k1 = 0; k1 <= (Nc-1); k1++){
					fprintf(fsig,"%e\t%e\t%e\t%e\n",((double)k1)*dxi,coeff1*out[k1][0], -coeff1*out[k1][1] , coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));
					}
	//fclose(file1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	return;

}


void printGaborFFTW3binary(FILE *Gsize, FILE *xgrid, FILE *xigrid, FILE *Gsigbin, double *signal, int N, double dx, double dxG, double a, double xiMaxPrint) // takes real signal speced by given "dx" and it computes and prints its Gabor transform, The parameters of the Gabor transform are new "dxG" (will be adjusted to a close one matching the points) and gabor parameter "a"
{
	int Nc, Ncprint;
	fftw_complex *out, *in2;
	double *in;
	double dxi,coeff1,coeff2;
	fftw_plan p;
	int k1,k2,kstep2,NxG;
	
	double dum;
	double *dumptr1, *dumptr2;


/*	// print signal*/
/*	for(k1 = 0; k1 <= (N-1); k1++){*/
/*					fprintf(sig,"%e\t%e\n", ((double)k1)*dx , signal[k1]);*/
/*					}*/
/*	//fclose(file1);*/

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

		
		// print fourier transform
		// file1 = fopen("ftransform.dat" , "w");

	
/*		dumptr1 = calloc(1,sizeof(double));*/

		for(k1 = 0; k1 <= Ncprint; k1++){dumptr1[k1] = sqrt( coeff2*(out[k1][0]*out[k1][0]+out[k1][1]*out[k1][1]));}
		fwrite( dumptr1 ,sizeof(double), (Ncprint+1) ,Gsigbin);
		free(dumptr1);



/*		dumptr[0] = sqrt( coeff2*(out[Ncprint][0]*out[Ncprint][0]+out[Ncprint][1]*out[Ncprint][1]));*/
/*		fprintf(Gsig,"%e\n", dumptr[0] ) ;*/
/*		fwrite( dumptr,sizeof(double),1,Gsigbin);*/

		
		dumptr2[0] = ((double)k2)*dx;
		fwrite(dumptr2 ,sizeof(double),1,xgrid);
		free(dumptr2);
		 

		//fclose(file1);
	}
	fprintf(Gsize,"%i\n", NxG);

	dumptr1 = calloc( (Ncprint+1) ,sizeof(double));
	for(k1 = 0; k1 <= Ncprint; k1++){dumptr1[k1] = ((double)k1)*dxi;}
	fwrite( dumptr1 ,sizeof(double), (Ncprint+1) ,xigrid);
	free(dumptr1);

	// !!!!! OUR CONVENTION OF ft IS COMLEX CONJUGATE WRT dft

	
	return;

}


	/* case 3: // sin^2

		
		// printf("test2\n");
		// components of the vec.pot. for both pulses in all cases
		switch ( F.sin2.overlap ){

		case 1:
			if (t < F.sin2.tmax1)
			{
				//printf("case 1 \n");
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos(F.sin2.o1*t + F.sin2.phi1);
				A2=0.;
			} else if ((t >= F.sin2.tmax1) && (t < F.sin2.tmin2)){
				//printf("case 2 \n");
				A1=0.;
				A2=0.;
			} else if ((t >= F.sin2.tmin2) && (t <= F.sin2.tmax2)){
				//printf("case 3 \n");
				A1 = 0.;
				A2 = F.sin2.A2*( pow(sin(F.sin2.oc2*t+F.sin2.phi0),2.) ) * cos(F.sin2.o2*t + F.sin2.phi2);
			}
			// printf("test2 %lf \n", A1);
		break;

		case 2:
			if (t < F.sin2.tmin2)
			{
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos(F.sin2.o1*t + F.sin2.phi1);
				A2 = 0.;
			} else if ((t >= F.sin2.tmin2) && (t < F.sin2.tmax1)){
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos( F.sin2.o1*t + F.sin2.phi1);
				A2 = F.sin2.A2*( pow(sin(F.sin2.oc2*t + F.sin2.phi0),2.) ) * cos(F.sin2.o2*t + F.sin2.phi2);
			} else if ((t >= F.sin2.tmax1) && (t <= F.sin2.tmax2)){
				A1 = 0.;
				A2 = F.sin2.A2*( pow(sin(F.sin2.oc2*t + F.sin2.phi0),2.) ) * cos(F.sin2.o2*t + F.sin2.phi2);	
			}
		break;

		case 3:
			if (t < F.sin2.tmin2)
			{
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos(F.sin2.o1*t + F.sin2.phi1);
				A2 = 0.;
			} else if ((t >= F.sin2.tmin2) && (t < F.sin2.tmax2)){
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos(F.sin2.o1*t + F.sin2.phi1);
				A2 = F.sin2.A2*( pow(sin(F.sin2.oc2*t + F.sin2.phi0),2.) ) * cos(F.sin2.o2*t + F.sin2.phi2);

			} else if ((t >= F.sin2.tmax2) && (t <= F.sin2.tmax1)){
				A1 = F.sin2.A1*( pow(sin(F.sin2.oc1*t),2.) ) * cos(F.sin2.o1*t + F.sin2.phi1);
				A2 = 0.;
			}

		break;
		}


		// full vector potential
		// dum = A1 + A2;
		// printf("test2 %lf \n", dum);
		return A1+A2; 

	break;










// interpolate by 2n-th order interpolation, uses 0's outside the grid, assuming uniform grid for extrapolation
double interpolate2n( int n, double x, double *x_grid, double* y_grid, int n_int) //!inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
{
	int k1,k2,k3,k4;
	double dx,y;
	double *a, *xdum, *ydum;

	
	
	k1=0;
	k2=0;
	xdum = calloc(2*n_int,sizeof(double)); ydum = calloc(2*n_int,sizeof(double)); a = calloc(2*n_int,sizeof(double));
	findinterval(n, x, x_grid, &k1, &k2);
	// printf("\ninside interpolate \n");
	// printf("interval,  %i \n",k1);
	// printf("interval,  %i \n",k2);

	if ( (k1 >= n_interp) && (k2 <= n-n_interp) )
	{
		for (k3 = 0 , k3 <= 2*n_int , k3++)
		{
			xdum[k3] = x[k1-n_int+k3];
			ydum[k3] = y[k1-n_int+k3];
		} 	
	} else {
		dx = x[1] - x[0];
		for (k3 = 0 , k3 <= 2*n_int , k3++)
		{
			xdum[k3] = ( (k1-n_int+k3) >= 0 ) ? x[k1-n_int+k3] : x[k1]-((double)(k1-n_int+k3))*dx;
			ydum[k3] = ( (k1-n_int+k3) >= 0 ) ? y[k1-n_int+k3] : 0. ;
		} 		
	}
	
	vander(xdum,a,ydum,2*n_interp); // compute vandermond

	// evaluate polynomial
	y = a[0];
	for(k1 = 1, k1 <= 2*n_int , k1++)
		
	

	
	if( k1 == -1 )
	{
		y=y_grid[0];
	} else if( k2 == n+1 ){
		y=y_grid[n];
	} else {
		y=y_grid[k1]+(x-x_grid[k1])*(y_grid[k2]-y_grid[k1])/(x_grid[k2]-x_grid[k1]);	
	}
	
	//y = 0.;
	//printf("interpolated value,  %lf \n",y);
	//printf("\n");	
	return y;

}

void vander(double x[], double w[], double q[], int n)
{
	int i,j,k;
	double b,s,t,xx;
	double *c;

	c = calloc(n,sizeof(double)); //dvector(1,n); 
	if (n == 1) w[1]=q[1];
	else {
		for (i=1;i<=n;i++) c[i]=0.0;
		c[n] = -x[1];
		for (i=2;i<=n;i++) {
			xx = -x[i];
			for (j=(n+1-i);j<=(n-1);j++) c[j] += xx*c[j+1];
			c[n] += xx;
		}
		for (i=1;i<=n;i++) {
			xx=x[i];
			t=b=1.0;
			s=q[n];
			for (k=n;k>=2;k--) {
				b=c[k]+xx*b;
				s += q[k-1]*b;
				t=xx*t+b;
			}
			w[i]=s/t;
		}
	}
	free(c);//free_dvector(c,1,n);
}
*/







void calc2FFTW3(int N, double dx, double xmax, double *signal1, double *signal2, double **xgrid, double **xigrid, double **fsig1, double **fsig2, double **fsig1M2, double **fsig2M2, int *Nxi) //takes real signal speced by given "dt" and it computes and prints its FFTW3
{
	int Nc;
	fftw_complex *out1, *out2;
	double *in, *dum;
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


	// int size2D = sizeof(double *) * Nc + sizeof(double) * 2 * Nc; // Nc-rows - the size required for the array in the memory
	// double *ptr1, *ptr2;
	*xgrid = (double*) calloc(N,sizeof(double));
	*xigrid = (double*) calloc(Nc,sizeof(double));
	*fsig1 = (double*) calloc(2*Nc,sizeof(double));
	*fsig2 = (double*) calloc(2*Nc,sizeof(double));
	*fsig1M2 = (double*) calloc(Nc,sizeof(double));
	*fsig2M2 = (double*) calloc(Nc,sizeof(double));

	// ptr1 = (double *) ((*fsig1) + Nc); ptr2 = (double *) ((*fsig2) + Nc);
	// for(k1=0; k1 < Nc; k1++){(*fsig1)[k1] = ptr1 + 2 * k1; (*fsig2)[k1] = ptr2 + 2 * k1;}

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
	double *in, *dum;
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


	// int size2D = sizeof(double *) * Nc + sizeof(double) * 2 * Nc; // Nc-rows - the size required for the array in the memory
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
