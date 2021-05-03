#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<fftw3.h>
#include "hdf5.h"
#include "mpi.h"

#include "numerical_constants.h"
#include "util.h"

clock_t start, finish;
clock_t start2, finish2;

extern double* timet,dipole;


void printresults(struct trg_def trg, struct Efield_var Efield, FILE *file1, int k, double *psi, int num_r, double *psi0, double tt, double *x, double dx, double Field, double Apot, double x_int, double dip_pre, struct outputs_def outputs)
{
	double dip,pop_re,pop_im,pop_tot,current,position,ion_prob2, dum;
	int j,k1,k2,k3,k4;

		// printf("test2,\t%i\n",k);
		// the population in the ground state (gauge dependent)
		 pop_re=0.; pop_im=0.;
		for(j = 0 ; j <= num_r ; j++) {pop_re = pop_re + psi[2*j]*psi0[2*j] + psi[2*j+1]*psi0[2*j+1]; pop_im = pop_im + psi[2*j]*psi0[2*j+1] - psi[2*j+1]*psi0[2*j];}
		pop_tot = pop_re*pop_re + pop_im*pop_im;


		// the gauge independent probability of the electron being between -x_int and x_int
		k1 = 0; k2 = 0; k3 = 0; k4 = 0;
		findinterval(num_r, -x_int, x, &k1, &k2);
		findinterval(num_r, x_int, x, &k3, &k4);
		ion_prob2 = 0;
		for(j=k1;j<=k4;j++){ion_prob2 = ion_prob2 + psi[2*j]*psi[2*j] + psi[2*j+1]*psi[2*j+1];}

		
		/*
		// calculation of <-grad(V)> (gauge independent dipole) !!!!!!!! SIGN?
		dip=0.; // dip_im=0.;
		
		for(j = 0 ; j <= num_r ; j++) {dip = dip + (psi[2*j]*psi[2*j] + psi[2*j+1]*psi[2*j+1])*gradpot(x[j],trg); 
		 // dip_im = dip_im + (psi[2*j]*psi[2*j+1] - psi[2*j+1]*psi[2*j])*gradpot(x[j]);
		 }
		*/
		dip = dip_pre;
		

		// calculation of <x> (gauge independent)
		position=0.;
		
		for(j = 0 ; j <= num_r ; j++)
		{
			position = position + (psi[2*j]*psi[2*j] + psi[2*j+1]*psi[2*j+1])*x[j]; 
		}
	       	

		current = 0.; // (gauge dependent, current+Apot is gauge-independent) 
		for(j = 1 ; j<= num_r-1 ; j++) 
		{
		current = current + psi[2*j]*(psi[2*(j+1)+1]-psi[2*(j-1)+1]) - psi[2*j+1]*(psi[2*(j+1)]-psi[2*(j-1)]);                 
		}
		current = current*0.5/dx;



		dum = 0.; // test
		for(j = 1 ; j<= num_r-1 ; j++) 
		{
		dum = dum + psi[2*j]*(psi[2*(j+1)+1]-psi[2*(j-1)+1]) + psi[2*j+1]*(psi[2*(j+1)]-psi[2*(j-1)]);                 
		}
		dum = dum*0.5/dx;

	

		//printf("test4,\t%i\n",k);
		//dum = tt/41.34144728;

		
		
		switch (Efield.fieldtype){
		case 0: case 1:
			fprintf(file1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",tt,Field,pop_tot,dip,-dip+Field,current,position,ion_prob2,(1.-pop_tot)*Field);
		break;
		case 2: case 3:
			fprintf(file1,"%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\t%e\n",tt,Field,Apot,pop_tot,dip,-dip+Field,current,position,ion_prob2,(1.-pop_tot)*Field,dum);
		break;
		}

		// save to outputs
		outputs.PopTot[k+1]=pop_tot;		
}