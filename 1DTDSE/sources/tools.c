/**
 * @file tools.c
 * @brief Contains functions for the core TDSE code.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include <math.h>
#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>
#include <hdf5.h>
#include "constants.h"
#include "util.h"
#include "tools_hdf5.h"
#include "structures.h"
#include "tridiag.h"
#include "tools_algorithmic.h"
#include "tools.h"

clock_t start, finish;
clock_t start2, finish2;

extern double* timet,dipole;

void Initialise_grid_and_ground_state(inputs_def *in)
{
	/* 
	Comment to the choice of the CV criterion:
	This number has to be small enough to assure a good convregence of the wavefunction
	if it is not the case, then the saclar product of the the ground state and the excited states 
	is not quite 0 and those excited appears in the energy analysis of the gorund states, so the propagation !!
	CV = 1E-25 has been choosen to have a scalar product of 10^-31 with the third excited state for num_r = 5000 and dx=0.1
	*/

	int k1;
	int size = 2*((*in).num_r+1);
	double *off_diagonal, *diagonal;
	(*in).psi0 = calloc(size,sizeof(double));
	for(k1=0;k1<=(*in).num_r;k1++){(*in).psi0[2*k1] = 1.0; (*in).psi0[2*k1+1] = 0.;}
	Initialise_grid_and_D2((*in).dx, (*in).num_r, &((*in).x), &diagonal, &off_diagonal); // !!!! dx has to be small enough, it doesn't converge otherwise
	(*in).Einit = Einitialise((*in).trg, (*in).psi0, off_diagonal, diagonal, off_diagonal, (*in).x, (*in).Eguess, (*in).CV, (*in).num_r); // originally, some possibility to have also excited state
	free(diagonal); free(off_diagonal);
}

void Initialise_grid_and_D2(double dx, int num_r, double **x, double **diagonal, double **off_diagonal) // Initialise ground-state
{
    int k1;
    double xmax = 0.5*num_r*dx;
	// printf("numr : %i \n",num_r);
	*x = calloc((num_r+1),sizeof(double));
	*off_diagonal = calloc(2*(num_r+1),sizeof(double));
	*diagonal = calloc(2*(num_r+1),sizeof(double));	


	//Initialisation Matrix corresponding to D2
	for(k1=0;k1<=num_r;k1++)
	{
		(*x)[k1] = (double)k1*dx-xmax;
		(*off_diagonal)[2*k1] = -0.5/(dx*dx); (*off_diagonal)[2*k1 + 1] = 0.;
		(*diagonal)[2*k1] = 1./(dx*dx); (*diagonal)[2*k1 + 1] = 0.;
	}	
}

//// POTENTIAL SPECIFICATION

double potential(double x,  trg_def trg)
{
	// Soft core potential
	return -1.0/sqrt(trg.a*trg.a+x*x);

	//return 0.025*x*x;
	
	//double a0 = 0.695, a1 = 0.5, a2 = 1., a3 = 1.;  
	//return -2.0/sqrt(a0*a0+x*x) + a1*exp(-pow( (x-a2)/a3 ,2));

	//double q = 3.;
	//return -pow((pow(trg.a,q)+pow(abs(x),q)), 1./q );
	
}


double gradpot(double x,  trg_def trg)
{
  
  return x*pow(trg.a*trg.a+x*x,-1.5);


  //return x*pow(c*c+(x-R*0.5)*(x-R*0.5),-1.5)+x*pow(c*c+(x+R*0.5)*(x+R*0.5),-1.5);

  
  //return 0.05*x;

  //double a0 = 0.695, a1 = 0.5, a2 = 1., a3 = 1.;  
  //return 2.0*x*pow(a0*a0+x*x,-1.5) - 2.*a0*(x-a2)*exp(-pow( (x-a2)/a3 ,2))/(a3*a3) ;
  //double q = 3.;
  //if ( x >= 0.){ return pow(abs(x),q-1) * pow((pow(trg.a,q)+pow(abs(x),q)), (1.+q)/q );}else{return -pow(abs(x),q-1) * pow((pow(trg.a,q)+pow(abs(x),q)), (1.+q)/q );}

}

double norme(double *x,int Num_r)
{
	int i;
	double sum = 0.;
	for(i=0;i<=Num_r;i++){sum = sum + x[2*i]*x[2*i] + x[2*i+1]*x[2*i+1];}
	return sum;
}
void normalise(double *x,int Num_r)
{
	int i;
	double sum = 0.;
	for(i=0;i<=Num_r;i++){sum = sum + x[2*i]*x[2*i] + x[2*i+1]*x[2*i+1];}
	sum = sqrt(sum);
	for(i=0;i<=Num_r;i++){x[2*i]/=sum;x[2*i+1]/=sum;}
}

double * extend_grid(double *pold,int size,int oldsize,int shift)
{
	int i;
	double *pnew;

	pnew = calloc(2*(size+oldsize+1),sizeof(double));

	if ((shift >= 0) && (size >= shift))
	{

	 for(i=oldsize+1;i<=oldsize+size;i++) {pnew[2*i] = 0.; pnew[2*i+1] = 0.;}

	 for(i=0;i<=oldsize;i++) 
	 {
		pnew[2*shift+2*(oldsize-i)] = pold[2*(oldsize-i)]; 
		pnew[2*shift+2*(oldsize-i)+1] = pold[2*(oldsize-i)+1]; 
	 }

	for(i=0;i<shift;i++) {pnew[2*i] = 0.; pnew[2*i+1] = 0.;}

	}
	else
	{printf("\nExpansion of the grid incorrect : shift <0 or size < shift \n\n");}
	
	free(pold);

	return pnew;
}


void window_analysis( trg_def trg, double dE,double Estep,double E_start,int num_E,int num_r,double dx,double *psi,double *dinf,double *d,double *dsup,double *x)
{

	
	double *dnew,*dnew2,*dinfnew,*dsupnew,*res,*res2,*psi2;
	double *dinfnew2,*dsupnew2;
	double prob;
	int i,j;
	FILE *fel; 


	dnew = calloc(2*(num_r+1),sizeof(double));
	dnew2 = calloc(2*(num_r+1),sizeof(double)); 
	dinfnew = calloc(2*(num_r+1),sizeof(double)); 
	dsupnew = calloc(2*(num_r+1),sizeof(double));
	dinfnew2 = calloc(2*(num_r+1),sizeof(double));
	dsupnew2 = calloc(2*(num_r+1),sizeof(double));
	res = calloc(2*(num_r+1),sizeof(double));
	res2 = calloc(2*(num_r+1),sizeof(double));
	psi2 = calloc(2*(num_r+1),sizeof(double));



	fel = fopen("electron_spectrum.dat","w");
	if(fel == NULL){ printf("Cannot open electron_spectrum.dat"); exit(1);} 
	printf("Working on the bin 00000");



	
/*
	for(j = 0 ; j<= num_r ; j++) 
	{	
			dinfnew[2*j] = -dinf[2*j]; dinfnew[2*j+1] = -dinf[2*j+1];
			dsupnew[2*j] = -dsup[2*j]; dsupnew[2*j+1] = -dsup[2*j+1];
	}
*/





	for(i=0; i<= num_E ; i++)
	{
	 		
	  for(j = 0 ; j<= num_r ; j++) 
	  {	
			//dnew[2*j] = E_start+Estep*i-d[2*j]-potential(x[j],trg)-dE/sqrt(8); dnew[2*j+1] = -d[2*j+1]-dE/sqrt(8);
			//dnew2[2*j] = E_start+Estep*i-d[2*j]-potential(x[j],trg)+dE/sqrt(8); dnew2[2*j+1] = -d[2*j+1]+dE/sqrt(8);

		  dnew[2*j] = 10*(E_start+Estep*i-dE/sqrt(8)-potential(x[j],trg))/12. - d[2*j]; dnew[2*j+1] = -d[2*j+1]-10*dE/(12*sqrt(8));
		  dnew2[2*j] = 10*(E_start+Estep*i+dE/sqrt(8)-potential(x[j],trg))/12. - d[2*j]; dnew2[2*j+1] = -d[2*j+1]+10*dE/(12*sqrt(8));
		    
		  dinfnew[2*j] = (E_start+Estep*i-dE/sqrt(8)-potential(x[j],trg))/12.-dinf[2*j]; dinfnew[2*j+1] = -dinf[2*j+1]-dE/(12*sqrt(8));
		  dsupnew[2*j] = (E_start+Estep*i-dE/sqrt(8)-potential(x[j+1],trg))/12.-dsup[2*j]; dsupnew[2*j+1] = -dsup[2*j+1]-dE/(12*sqrt(8));	 

		  dinfnew2[2*j] = (E_start+Estep*i+dE/sqrt(8)-potential(x[j],trg))/12.-dinf[2*j]; dinfnew2[2*j+1] = -dinf[2*j+1]+dE/(12*sqrt(8));
		  dsupnew2[2*j] = (E_start+Estep*i+dE/sqrt(8)-potential(x[j+1],trg))/12.-dsup[2*j]; dsupnew2[2*j+1] = -dsup[2*j+1]+dE/(12*sqrt(8));

	  }



/*
	  Inv_Tridiagonal_Matrix_complex(dinfnew,dnew,dsupnew,psi,res,num_r+1);
		
	  for(j=0; j<= num_r ; j++) {psi2[2*j] = res[2*j];psi2[2*j+1] = res[2*j+1];}

	  Inv_Tridiagonal_Matrix_complex(dinfnew,dnew2,dsupnew,psi2,res2,num_r+1);
*/




	  Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew,dnew,dsupnew,psi,res,num_r);
		
	  for(j=0; j<= num_r ; j++) {psi2[2*j] = res[2*j];psi2[2*j+1] = res[2*j+1];}

	  Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew2,dnew2,dsupnew2,psi2,res2,num_r);



	  prob = norme(res2,num_r);
      	  prob = prob*dx*pow(dE,4.);

	  fprintf(fel,"%e\t%e\n",E_start+Estep*i,prob);
	  
	  printf("\b\b\b\b\b%5d", i); fflush(stdout); 

	}

	printf("\n");



	free(dnew); free(dnew2); free(dinfnew);
	free(res);free(dsupnew);
	free(dinfnew2);free(dsupnew2);
	free(res2);free(psi2); 
	fclose(fel);

}


void dipole_analysis(double num_w,double dw,double *timet,double *dipole, int nc, int num_t)
{

	int i,j;
	double FFT_re,FFT_im,omega;
	FILE *fhhg;
	
	fhhg = fopen("HHG_spectrum.dat", "w" );
	printf("Working on the bin 00000");

	
	omega = 0;
	for (i = 0 ; i < num_w ; i++)
	{
		
		FFT_re = 0; FFT_im = 0;
		for( j = 0 ; j < nc*num_t ; j++)
		{
			FFT_re += dipole[2*j]*cos(omega*timet[j])-dipole[2*j+1]*sin(omega*timet[j]);
			FFT_im += dipole[2*j]*sin(omega*timet[j])+dipole[2*j+1]*cos(omega*timet[j]);
		}
	
	   fprintf(fhhg,"%e\t%e\t%e\n",omega,FFT_re,FFT_im);
	   printf("\b\b\b\b\b%5d", i); fflush(stdout);
		
	   omega += dw;
	}

    printf("\n"); fclose(fhhg);

}



void projection_analysis(double Estep,double E_start,int num_E,int num_r,double dx,double *psi,double *x)
{


	double prob,prob_re,prob_im,k,delta;
	FILE *fel;
	int i,j,jmin,jmax;


	fel = fopen("electron_spectrum.dat", "w" );
	printf("Working on the bin 00000");

	num_E = 1000;
	Estep = 0.001;	jmin = 7000; jmax = 8000; delta = (double) jmax- (double) jmin;
		
	for(i=0; i<=  num_E ; i++)
	{
	 		
	  k = sqrt(2.0*Estep*i); 
	  prob_re = 0.0; prob_im = 0.0;
          for(j = jmin ; j<= jmax ; j++) 
	  {	
		prob_re +=  psi[2*j]*cos(k*x[j])-psi[2*j+1]*sin(k*x[j]);  
		prob_im +=  psi[2*j]*sin(k*x[j])+psi[2*j+1]*cos(k*x[j]);

	  }

	  prob_re *= dx/delta; prob_im *= dx/delta; 	


	  prob = (prob_re*prob_re+prob_im*prob_im)/(2.0*Pi);


	  fprintf(fel,"%e\t%e\n",Estep*i,prob);
	  printf("\b\b\b\b\b%5d", i); fflush(stdout); 

	}

	printf("\n");

	fclose(fel);


}

void projection_analysis_EV( trg_def trg, double dE,double Estep,double E_start,int num_E,int num_r,double dx,double *psi,double *dinf,double *d,double *dsup,double *x) // inti procedure incompatible
{


	double prob,CV,Eguess,E_previous,E,ps_re,ps_im;
	double *psi_EV;
	FILE *fel;
	int i,j;

	psi_EV = calloc(2*(num_r+1),sizeof(double));


	fel = fopen("electron_spectrum.dat", "w" );
	printf("Working on the bin 00000");

	num_E = 2000;
	Estep = 0.001;	


	CV = 1E-10; // CV criteria  
	E_previous = 0;
		
	for(i=0; i<= num_E ; i++)
	{

		for(j=0;j<=num_r;j++) {psi_EV[2*j] = 1; psi_EV[2*j+1] = 0.;}
		normalise(psi_EV,num_r); // Initialise psi_EV for Einitialise

	 	Eguess = i*Estep; 
		
		E = Einitialise(trg,psi_EV,dinf,d,dsup,x,Eguess,CV,num_r);
		printf("%e\t%e\n",Eguess,E);


		ps_re = 0; ps_im = 0;
		for(j=0;j<=num_r;j++)
		{
			ps_re += (psi[2*j]*psi_EV[2*j]-psi[2*j+1]*psi_EV[2*j+1])*x[j];
			ps_im += (psi[2*j]*psi_EV[2*j+1]+psi[2*j+1]*psi_EV[2*j])*x[j];
			//fprintf(fel,"%e\t%e\t%e\n",x[j],psi_EV[2*j],psi_EV[2*j+1]);
		}

	 prob = ps_re*ps_re + ps_im*ps_im; 
	 // if ( Eguess != E_previous)
	  //{	
	    fprintf(fel,"%e\t%e\n",E,prob); //E_previous = E;
	  //}
	  //printf("\b\b\b\b\b%5d", i); fflush(stdout); 

	}

	printf("\n");

	fclose(fel); free(psi_EV);


}


