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

// extern struct Efield_var;


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



double* propagation(struct trg_def trg, struct Efield_var Efield, double tmin, int Nt, int num_t,double dt,int num_r
,int num_exp,double dx,double *psi0,double *psi,double *x
,FILE *timef,FILE *timef2,double ton,double toff, double *timet, double *dipole, int gauge, int transformgauge, double x_int, struct analy_def analy, struct outputs_def outputs)
{	
	
	double *res1,*dnew1,*dinfnew1,*dsupnew1,*psi_inter1;
	double *res2,*dnew2,*dinfnew2,*dsupnew2,*psi_inter2;
	double Field,tt,coef,Apot,atten,phiabs;
	int i,j,save,comp,index,k,ip,shift,k1,k2,k3,k4,k5,k6;	
	FILE *wf,*zplot,*phase,*wft,*testfile;
	char name[20];
	double cpot,c0re,c0im;
	double dip,dip_re,dip_im,pop_re,pop_im;
        double current, ion_prob2, ttrig; 
	double av_x_re,av_px_re,av_x_im,av_px_im,derpsi_re,derpsi_im;

	double t_zero1,t_zero2;

	double dum, dum2;
	
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

	
	outputs.tgrid[0] = tt, outputs.sourceterm[0] = 0.; outputs.Efield[0]=Field; outputs.PopTot[0]=1.0;

	for(j = 0 ; j<= num_r ; j++) {psi[2*j] = psi0[2*j]; psi[2*j+1] = psi0[2*j+1];}

	
	start2 = clock();
	
	for(k = 0 ; k < Nt ; k++)
	{

		//printf("tcycle %i \n",k); fflush(NULL);	
		if( k%num_t == 0 )
		{
		start = clock();	
		printf("Cycle number : %i ; size of the box : %i ; progress %i/%i \n",(k/num_t)+1,num_r,k,Nt); fflush(NULL);
		}


		tt = tt + dt;
		
		coef = 0.5*dt/(dx*dx);

		
		
		comp = 100; index = 0; save = 1;

		Field = Efield.Field[k];


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
	printf("\nDuration of calculation for the whole problem %f sec\n\n",(double)(finish2 - start2) / CLOCKS_PER_SEC);

	
	free(psi_inter1);free(res1);free(dnew1);free(dinfnew1);free(dsupnew1);
    	free(psi_inter2);free(res2);free(dnew2);free(dinfnew2);free(dsupnew2);
	free(psi2);

	return psi;
}


void window_analysis(struct trg_def trg, double dE,double Estep,double E_start,int num_E,int num_r,double dx,double *psi,double *dinf,double *d,double *dsup,double *x)
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

void projection_analysis_EV(struct trg_def trg, double dE,double Estep,double E_start,int num_E,int num_r,double dx,double *psi,double *dinf,double *d,double *dsup,double *x) // inti procedure incompatible
{


	double prob,prob_re,prob_im,k,delta,CV,Eguess,E_previous,E,ps_re,ps_im;
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


