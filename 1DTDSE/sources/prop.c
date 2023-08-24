#include "tools_MPI-RMA.h"
#include "numerical_constants.h"
#include "prop.h"
#include "tools_fftw3.h"
#include "tridiag.h"
#include "structures.h"
#include "tools.h"
#include "tools_algorithmic.h"

clock_t start, finish;
clock_t start2, finish2;
double MPI_clock_start, MPI_clock_finish;

extern double* timet,dipole;

// extern  Efield_var;


double* propagation( trg_def trg,  Efield_var Efield, double tmin, int Nt, int num_t,double dt,int num_r
,int num_exp,double dx,double *psi0,double *psi,double *x
,FILE *timef,FILE *timef2,double ton,double toff, double *timet, double *dipole, int gauge, int transformgauge, double x_int,  analy_def analy,  outputs_def outputs)
{	
	
	double *res1,*dnew1,*dinfnew1,*dsupnew1,*psi_inter1;
	double *res2,*dnew2,*dinfnew2,*dsupnew2,*psi_inter2;
	double Field,tt,coef,Apot,atten;
	int j,save,comp,index,k,ip,k1,k2,k3,k4;	
	
	double cpot;
	double dip;
	double ion_prob2, ttrig; 
	

	double t_zero1,t_zero2;

	
	
	double *psi2;	

	psi_inter1 = calloc(2*(num_r+1),sizeof(double));
	res1 = calloc(2*(num_r+1),sizeof(double));
	dnew1 = calloc(2*(num_r+1),sizeof(double));
	dinfnew1 = calloc(2*(num_r+1),sizeof(double));
	dsupnew1 = calloc(2*(num_r+1),sizeof(double));
	psi_inter2 = calloc(2*(num_r+1),sizeof(double));
	res2 = calloc(2*(num_r+1),sizeof(double));
	dnew2 = calloc(2*(num_r+1),sizeof(double));
	dinfnew2 = calloc(2*(num_r+1),sizeof(double));
	dsupnew2 = calloc(2*(num_r+1),sizeof(double));

	psi2 = calloc(2*(num_r+1),sizeof(double));


	tt = tmin; ttrig=tmin;

	comp = 100; index = 0; save = 1;

		// the gauge independent probability of the electron being between -x_int and x_int for psi0
		k1 = 0; k2 = 0; k3 = 0; k4 = 0;
		findinterval(num_r, -x_int, x, &k1, &k2);
		findinterval(num_r, x_int, x, &k3, &k4);
		ion_prob2 = 0;
		for(j=k1;j<=k4;j++){ion_prob2 = ion_prob2 + psi0[2*j]*psi0[2*j] + psi0[2*j+1]*psi0[2*j+1];}



		Field = 0;
		// variables for apodization
		t_zero1 = tmin;
		t_zero2 = findnextinterpolatedzero(Efield.Nt-1, t_zero1 + Efield.dt, Efield.tgrid, Efield.Field);		
	ip = 0;

	cpot = 1.;

	
	outputs.tgrid[0] = tt, outputs.sourceterm[0] = 0.; outputs.Efield[0]=Field; outputs.PopTot[0]=1.0; // wrong, first values shall be computed
	outputs.PopInt[0]=ion_prob2; outputs.expval[0]=0.0;


	for(j = 0 ; j<= num_r ; j++) {psi[2*j] = psi0[2*j]; psi[2*j+1] = psi0[2*j+1];}

	
	start2 = clock();
	MPI_clock_start = MPI_Wtime(); 
	
	int do_zeroing = 0;
	for(k = 0 ; k < Nt ; k++)
	{

		// printf("tcycle %i \n",k); fflush(NULL);	
		if( k%num_t == 0 )
		{
		start = clock();	
		// printf("Cycle number : %i ; size of the box : %i ; progress %i/%i \n",(k/num_t)+1,num_r,k,Nt); fflush(NULL);
		}


		tt = tt + dt;
		
		coef = 0.5*dt/(dx*dx);

		
		
		comp = 100; index = 0; save = 1;

		if(do_zeroing == 0){
			if(Efield.Field[k]*Efield.Field[k+1] <= 0.0){do_zeroing = 1;}
			Field = 0.;
		}else{
			Field = Efield.Field[k]; 
		}
		
		// if(tt <= t_zero2)
		// {
		// 	Field = 0.;
		// }else{
		// 	Field = Efield.Field[k];
		// }
		


		for(j = 0 ; j<= num_r ; j++) 
		{	
			dinfnew1[2*j] = 1/12.; dinfnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[j],trg));
			dnew1[2*j] = 10/12.; dnew1[2*j+1] = 0.5*dt*( 1./(dx*dx) )+0.5*dt*10/12.*(cpot*potential(x[j],trg));
			//dsupnew1[2*j] = 1/12.; dsupnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[j+1],trg));			
		}
		for(j = 0 ; j<num_r ; j++) { dsupnew1[2*j] = 1/12.; dsupnew1[2*j+1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[j+1],trg));}
		dsupnew1[2*num_r ] = 1/12.; dsupnew1[2*num_r +1] = 0.5*dt*( -0.5/(dx*dx) )+0.5*dt*1/12.*(cpot*potential(x[num_r]+dx,trg));	
		
		// first part of the evolution (H0+V)

		psi_inter1[0] = (10/12.)*psi[0]+coef*psi[1]+1/12.*psi[2]-0.5*coef*psi[3];
		psi_inter1[0] = psi_inter1[0]+0.5*dt*((10/12.)*psi[1]*(cpot*potential(x[0],trg))
						+(1/12.)*psi[3]*(cpot*potential(x[1],trg)));

		psi_inter1[1] = (10/12.)*psi[1]-coef*psi[0]+1/12.*psi[3]+0.5*coef*psi[2];	
		psi_inter1[1] = psi_inter1[1]-0.5*dt*((10/12.)*psi[0]*(cpot*potential(x[0],trg))
						+(1/12.)*psi[2]*(cpot*potential(x[1],trg)));

		for(j = 1 ; j< num_r ; j++)
		{

			psi_inter1[2*j] = (10/12.)*psi[2*j]+coef*psi[2*j+1]+1/12.*psi[2*(j+1)]
							  +1/12.*psi[2*(j-1)]-0.5*coef*(psi[2*(j-1)+1]+psi[2*(j+1)+1]);
			psi_inter1[2*j] = psi_inter1[2*j]+0.5*dt*((10/12.)*psi[2*j+1]*(cpot*potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)+1]*(cpot*potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)+1]*(cpot*potential(x[j+1],trg)));

			psi_inter1[2*j+1] = (10/12.)*psi[2*j+1]-coef*psi[2*j]+1/12.*psi[2*(j+1)+1]
							  +1/12.*psi[2*(j-1)+1]+0.5*coef*(psi[2*(j-1)]+psi[2*(j+1)]);
			psi_inter1[2*j+1] = psi_inter1[2*j+1]-0.5*dt*((10/12.)*psi[2*j]*(cpot*potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)]*(cpot*potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)]*(cpot*potential(x[j+1],trg)));

		}

		psi_inter1[2*num_r] = (10/12.)*psi[2*num_r]+coef*psi[2*num_r+1]+1/12.*psi[2*(num_r-1)]-0.5*coef*psi[2*(num_r-1)+1];
		psi_inter1[2*num_r] = psi_inter1[2*num_r]+0.5*dt*((10/12.)*psi[2*num_r+1]*(cpot*potential(x[num_r],trg))
							  +(1/12.)*psi[2*(num_r-1)+1]*(cpot*potential(x[num_r-1],trg)));

		psi_inter1[2*num_r+1] = (10/12.)*psi[2*num_r+1]-coef*psi[2*num_r]+1/12.*psi[2*(num_r-1)+1]+0.5*coef*psi[2*(num_r-1)];
		psi_inter1[2*num_r+1] = psi_inter1[2*num_r+1]-0.5*dt*((10/12.)*psi[2*num_r]*(cpot*potential(x[num_r],trg))
							    +(1/12.)*psi[2*(num_r-1)]*(cpot*potential(x[num_r-1],trg)));


		Inv_Tridiagonal_Matrix_complex(dinfnew1,dnew1,dsupnew1,psi_inter1,res1,num_r+1);


		// second part of the evolution (Hint)

		if( gauge == 0 )
		{
			for(j = 0 ; j<= num_r ; j++) 
			{
				psi[2*j] = cos(Field*dt*x[j])*res1[2*j]-sin(Field*dt*x[j])*res1[2*j+1]; 
				psi[2*j+1] = cos(Field*dt*x[j])*res1[2*j+1]+sin(Field*dt*x[j])*res1[2*j];
			}
		}
		else // velocity gauge (A has to be available)
		{

			for(j = 0 ; j<= num_r ; j++) 
			{			
				dinfnew2[2*j] = 1/6.+0.5*dt*Apot*0.5/dx; dinfnew2[2*j+1] = 0;
				dnew2[2*j] = 4/6.; dnew2[2*j+1] = 0;
				dsupnew2[2*j] = 1/6.-0.5*dt*Apot*0.5/dx; dsupnew2[2*j+1] = 0;
			}

			psi_inter2[0] = 4/6.*res1[0]+(1/6.+0.5*dt*Apot*0.5/dx)*res1[2];
			psi_inter2[1] = 4/6.*res1[1]+(1/6.+0.5*dt*Apot*0.5/dx)*res1[3];

			for(j = 1 ; j< num_r ; j++)
			{

				psi_inter2[2*j] = 4/6.*res1[2*j] + (1/6. + 0.5*dt*Apot*0.5/dx)*res1[2*(j+1)];
				psi_inter2[2*j] = psi_inter2[2*j] + (1/6. - 0.5*dt*Apot*0.5/dx)*res1[2*(j-1)];
				psi_inter2[2*j+1] = 4/6.*res1[2*j+1] + (1/6. + 0.5*dt*Apot*0.5/dx)*res1[2*(j+1)+1];
				psi_inter2[2*j+1] = psi_inter2[2*j+1] + (1/6. - 0.5*dt*Apot*0.5/dx)*res1[2*(j-1)+1];

			}

			psi_inter2[2*num_r] = 4/6.*res1[2*num_r]+(1/6.-0.5*dt*Apot*0.5/dx)*res1[2*(num_r-1)];
			psi_inter2[2*num_r+1] = 4/6.*res1[2*num_r+1]+(1/6.-0.5*dt*Apot*0.5/dx)*res1[2*(num_r-1)+1];	

			Inv_Tridiagonal_Matrix_complex(dinfnew2,dnew2,dsupnew2,psi_inter2,res2,num_r+1);
			for(j = 0 ; j<= num_r ; j++) 
			{
				atten = 1.;		
		        	psi[2*j] = atten*res2[2*j]; 
				psi[2*j+1] = atten*res2[2*j+1];	
			}

		}
			

///////// ABSORBING WALL (REINTRODUCE!)


		// // absorbing wall
		// for(j = 0 ; j<= num_r ; j++) 
		// {
		// 	if ( j >= (num_r-300) ) 
		// 	{	phiabs = x[j]-x[num_r-300];
		// 		phiabs *= 1./(x[num_r]-x[num_r-300]);
		// 		atten = 1-phiabs;
		// 		phiabs *= Pi*0.5;
		// 		atten *= cos(phiabs);
		// 	}
		// 	else {atten = 1;}

		// 	if ( j <= 300 ) 
		// 	{	phiabs = x[j]-x[300];
		// 		phiabs *= 1./(x[0]-x[300]);
		// 		atten = 1-phiabs;
		// 		phiabs *= Pi*0.5;
		// 		atten *= cos(phiabs);
		// 	}
		// 	else {if (j < (num_r-300)) atten = 1;}
		// }


		// PRINTING
		dip=0.; for(k1 = 0 ; k1 <= num_r ; k1++) {dip = dip + (psi[2*k1]*psi[2*k1] + psi[2*k1+1]*psi[2*k1+1])*gradpot(x[k1],trg);};
		outputs.tgrid[k+1] = tt, outputs.sourceterm[k+1] = -dip+Field; outputs.Efield[k+1]=Field;
		//outputs.tgrid[k+1] = tt, outputs.sourceterm[k+1] = Field; outputs.Efield[k+1]=Field;


		// printresults(trg,Efield, timef,k,psi,num_r,psi0,tt,x,dx,Field,Apot,x_int,dip,outputs); population was computed there
		compute_population(trg,Efield,k,psi,num_r,psi0,tt,x,dx,Field,Apot,x_int,dip,outputs);

		if (ip == (int)floor(num_t/20.))
		{
			// printf("*/"); fflush(stdout);
			ip = 0;
		}
		ip++;


	if( ( k%num_t == num_t-1) && ( k != 0 ) )
	{
	finish = clock();
	}

	

	

	} // end of the main loop

	finish2 = clock();
	MPI_clock_finish = MPI_Wtime(); 
	// printf("\nDuration of calculation for the whole problem %f sec, MPI time %f sec\n\n",(double)(finish2 - start2) / CLOCKS_PER_SEC, MPI_clock_finish-MPI_clock_start);

	
	free(psi_inter1);free(res1);free(dnew1);free(dinfnew1);free(dsupnew1);
    	free(psi_inter2);free(res2);free(dnew2);free(dinfnew2);free(dsupnew2);
	free(psi2);

	return psi;
}

void compute_population( trg_def trg,  Efield_var Efield, int k, double *psi, int num_r, double *psi0, double tt, double *x, double dx, double Field, double Apot, double x_int, double dip_pre,  outputs_def outputs)
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

		outputs.expval[k+1]=position;
		outputs.PopInt[k+1]=ion_prob2;	
}

