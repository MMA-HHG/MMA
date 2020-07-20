#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#include "numerical_constants.h"
#include "util.h"

clock_t start, finish;
clock_t start2, finish2;

extern double* timet,dipole;

void Initialise_grid_and_ground_state(struct inputs_def *in)
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



struct output_print_def Initialise_Printing_struct(void) // Initialise ground-state
{
	struct output_print_def res;

	res.Efield = 0;
	res.FEfield = 0;
	res.sourceterm = 0;
	res.Fsourceterm = 0;
	res.FEfieldM2 = 0;
	res.FsourceTermM2 = 0;
	res.PopTot = 0;
	res.tgrid = 0;
	res.omegagrid = 0;

	return res;
}

struct output_print_def Set_all_prints(void) // Initialise ground-state
{
	struct output_print_def res;

	res.Efield = 1;
	res.FEfield = 1;
	res.sourceterm = 1;
	res.Fsourceterm = 1;
	res.FEfieldM2 = 1;
	res.FsourceTermM2 = 1;
	res.PopTot = 1;
	res.tgrid = 1;
	res.omegagrid = 1;

	return res;
}

// void outputs_constructor(struct outputs_def *outputs, int Nt) 
// {
// 	(*outputs).tgrid = calloc((Nt+1),sizeof(double));
// 	(*outputs).Efield = calloc((Nt+1),sizeof(double));
// 	(*outputs).sourceterm = calloc((Nt+1),sizeof(double));
// 	(*outputs).PopTot = calloc((Nt+1),sizeof(double));
// 	(*outputs).Nt = Nt+1;
// }

void outputs_destructor(struct outputs_def *outputs) // frees memory allocated for outputs
{
	free((*outputs).tgrid);
	free((*outputs).omegagrid);
	free((*outputs).tgrid_fftw);

	free((*outputs).Efield);
	free((*outputs).sourceterm);
	free((*outputs).PopTot);

	free((*outputs).FEfield_data);
	free((*outputs).Fsourceterm_data);

	free((*outputs).FEfieldM2);
	free((*outputs).FsourcetermM2);

	(*outputs).Nt = 0;
	(*outputs).Nomega = 0;
}

void inputs_destructor(struct inputs_def *in) // frees memory allocated for outputs
{
	free((*in).psi0);
	free((*in).x);
}

void Initialise_grid_and_D2(double dx, int num_r, double **x, double **diagonal, double **off_diagonal) // Initialise ground-state
{
    int k1;
    double xmax = 0.5*num_r*dx;
	printf("numr : %i \n",num_r);
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






void compute_population(struct trg_def trg, struct Efield_var Efield, int k, double *psi, int num_r, double *psi0, double tt, double *x, double dx, double Field, double Apot, double x_int, double dip_pre, struct outputs_def outputs)
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

	
		// save to outputs
		outputs.PopTot[k+1]=pop_tot;		
}



//// POTENTIAL SPECIFICATION

double potential(double x, struct trg_def trg)
{
	return -1.0/sqrt(trg.a*trg.a+x*x);

	//return 0.025*x*x;
	
	//double a0 = 0.695, a1 = 0.5, a2 = 1., a3 = 1.;  
	//return -2.0/sqrt(a0*a0+x*x) + a1*exp(-pow( (x-a2)/a3 ,2));

	//double q = 3.;
	//return -pow((pow(trg.a,q)+pow(abs(x),q)), 1./q );
	
}


double gradpot(double x, struct trg_def trg)
{
  
  return x*pow(trg.a*trg.a+x*x,-1.5);


  //return x*pow(c*c+(x-R*0.5)*(x-R*0.5),-1.5)+x*pow(c*c+(x+R*0.5)*(x+R*0.5),-1.5);

  
  //return 0.05*x;

  //double a0 = 0.695, a1 = 0.5, a2 = 1., a3 = 1.;  
  //return 2.0*x*pow(a0*a0+x*x,-1.5) - 2.*a0*(x-a2)*exp(-pow( (x-a2)/a3 ,2))/(a3*a3) ;
  //double q = 3.;
  //if ( x >= 0.){ return pow(abs(x),q-1) * pow((pow(trg.a,q)+pow(abs(x),q)), (1.+q)/q );}else{return -pow(abs(x),q-1) * pow((pow(trg.a,q)+pow(abs(x),q)), (1.+q)/q );}

}