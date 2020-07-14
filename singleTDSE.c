#include<time.h>
#include<stdio.h>
#include<string.h> 
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"
#include "mpi.h"



#include"util.h"

struct Efield_var Efield;
struct trg_def trg;
struct analy_def analy;
// struct outputs_def outputs;
double *y,*x,*a,*c,sum,*diagonal,*off_diagonal,*eigenvector,*u,*r,*vector;
double *timet,*dipole;
double dx,xmax,Eguess,Einit,CV,phi,omega,E0,period,Pi,tfinal,alpha,v,mod1,mod2,dE,Estep,norm_gauss,x_int,textend;

int gauge,transformgauge,fieldinau,input0,Ntinterp,InterpByDTorNT;

double dt, tmax,tmin;
int Nt;

double dum,dum1,dum2;
double *psi0,*psi,Einit2,ps_re,ps_im,*psi_rmv_gs;
double E_start,ton,toff,dw;
int num_E,num_exp,num_w,N_t,dumint;


double a_Gabor, omegaMaxGabor, dtGabor, tmin1window, tmax1window, tmin2window, tmax2window, IonFilterThreshold;
int PrintGaborAndSpectrum, IonisationFilterForTheSourceTerm;

int PrintOutputMethod;


double *test_expand;

double *dinf,*dsup,*d,*u1,*res,*t;

int i,j,k,l,m,num_r,num_t,err,max_iteration_count,size,nc,k1,k2,k2,k3,k5;
clock_t start, finish;

char ch;

int size_exp,shift;

FILE *newygrid,*eingenvaluef,*eingenvectorf,*timef,*timef2,*gaussianwp,*volkovwp,*param,*pot,*file1,*file3,*file2,*file4,*file5,*file6,*file7,*file8,*file9;

char filename1[25], filename2[25];




struct outputs_def call1DTDSE(struct inputs_def inputs) // this is a wrapper that will by bypassed. It's here due to the design of the original code.
{
	// declarations
	struct outputs_def outputs;	
	Pi = acos(-1.);

	///////////////////////////////////////////////
	// local copies of variables given by inputs //
	///////////////////////////////////////////////

	// (bypass ?)

	Efield.fieldtype = 0; // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical

	Eguess = inputs.Eguess; // Energy of the initial state
	num_r = inputs.num_r; // Number of points of the initial spatial grid 16000
	num_exp = inputs.num_exp; // Number of points of the spatial grid for the expansion
	dx = inputs.dx; // resolution for the grid
	InterpByDTorNT = inputs.InterpByDTorNT; // Number of points of the spatial grid for the expansion
	dt = inputs.dt; // resolution in time
	Ntinterp = inputs.Ntinterp; // Number of points of the spatial grid for the expansion
	textend = inputs.textend; // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
	analy.writewft = inputs.analy.writewft; // writewavefunction (1-writting every tprint)
	analy.tprint = inputs.analy.tprint; // time spacing for writing the wavefunction	
	x_int = inputs.x_int; // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
	PrintGaborAndSpectrum = inputs.PrintGaborAndSpectrum; // print Gabor and partial spectra (1-yes)
	a_Gabor = inputs.a_Gabor; // the parameter of the gabor window [a.u.]
	omegaMaxGabor = inputs.omegaMaxGabor; // maximal frequency in Gabor [a.u.]
	dtGabor = inputs.dtGabor; // spacing in Gabor
	tmin1window = inputs.tmin1window; // analyse 1st part of the dipole
	tmax1window = inputs.tmax1window; // analyse 1st part of the dipole
	tmin2window = inputs.tmin2window; // analyse 2nd part of the dipole
	tmax2window = inputs.tmax2window; // analyse 2nd part of the dipole
	PrintOutputMethod = inputs.PrintOutputMethod; // (0 - only text, 1 - only binaries, 2 - both)
	// IonisationFilterForTheSourceTerm = inputs.IonisationFilterForTheSourceTerm; // filter source term by high-ionisation components (1-yes)
	// IonFilterThreshold = inputs.IonFilterThreshold; // threshold for the ionisation [-]
	trg.a = inputs.trg.a; // the limit of the integral for the ionisation //2 2 works fine with the
 
	Efield = inputs.Efield;


	gauge = 1;
	transformgauge = 0;
	input0 = 1;
	printf("tgrid,  %e, %e \n",Efield.tgrid[0],Efield.tgrid[1]);


	////////////////////////////////
	// PREPARATIONAL COMPUTATIONS //
	////////////////////////////////

	// find dt from the grid around 0	
	switch ( input0 ){case 0: dumint = 0; break; case 1: dumint = round(Efield.Nt/2.); /* field centered around 0 */ break;} // choosing the best resolution	
	Efield.dt = Efield.tgrid[dumint+1]-Efield.tgrid[dumint]; // 
	
	tmax = Efield.tgrid[Efield.Nt]-Efield.tgrid[0]; // total length of the grid
   
	// refine dt either by given number of points or by required dt
	if (InterpByDTorNT == 1){k1 = Ntinterp + 1;} else { k1 = floor(Efield.dt/dt); Ntinterp = k1; k1++; }
	dt = Efield.dt/((double)k1); // redefine dt properly


	Efield.Field = FourInterp(k1, Efield.Field, Efield.Nt); // make the interpolation !!!!!! tgrid does not correspond any more
	Nt = k1*Efield.Nt + 1; // not very nice to compute it, FourInterp should do it

	// PRINT new field and its transform == assign to outputs
	// file1 = fopen("inputs/InputField2.dat" , "w"); file2 = fopen("inputs/InputFField2.dat" , "w"); printFFTW3(file1, file2, Efield.Field, Nt, dt); fclose(file1); fclose(file2);

	
	num_t = floor((2*Pi)/(0.057*dt)); num_t++;  // the length of one cycle for 800 nm (i.e. omega=0.057) 
	printf("Efield.dt,  %e \n",Efield.dt);  fflush(NULL);
	printf("dt,  %e \n",dt);

	size = 2*(num_r+1);// for complex number


	// Allocation memory

	//printf("test1 %i \n",size);

	// ALLOCATE MEMORY, COPY INITIAL ARRAYS AND PREPARE THEM FOR THE PROPAGATOR

	x = malloc((num_r+1)*sizeof(double)); // we keep this construction and not use directly initial x due to the extensibility of the grid
	memcpy(x,inputs.x,(num_r+1)*sizeof(double));

	psi0 = malloc(size*sizeof(double));	
	memcpy(psi0,inputs.psi0,(num_r+1)*sizeof(double));

	psi = calloc(size,sizeof(double));
	t = calloc(Nt,sizeof(double));
	timet = calloc(Nt,sizeof(double));
	dipole = calloc(2*Nt,sizeof(double));



	// prepare outputs (there should be written a constructor depending on required values and a destructor on allocated memory should be called in the main code)
	// ineficient allocate every time... Should be done by reference as well.
	// outputs_constructor(outputs,Nt);
	outputs.tgrid = calloc((Nt+1),sizeof(double)); outputs.Efield = calloc((Nt+1),sizeof(double)); outputs.sourceterm = calloc((Nt+1),sizeof(double)); outputs.PopTot = calloc((Nt+1),sizeof(double));
	outputs.Nt = (Nt+1);
	 


	//printf("\n");	
	//printf("Propagation procedure ...\n");
	//printf("\n");	



	start = clock();

	psi = propagation(trg,Efield,tmin,Nt,num_t,dt,num_r,num_exp,dx,psi0,psi,x,timef,timef2,ton,toff,timet,dipole,gauge,transformgauge,x_int,analy,outputs);
	printf("TDSE done ...\n");


	
	// TEST filtering for high ionisation // MORE EFFICIENT WOULD BE FILTER WHILE ASSIGNING VALUE
	if(IonisationFilterForTheSourceTerm == 1){
		outputs.sourcetermfiltered = calloc((Nt+1),sizeof(double));
		for(k1 = 0; k1 <= Nt; k1++){ outputs.sourcetermfiltered[k1] = outputs.sourceterm[k1];}
		for(k1 = 0; k1 <= Nt; k1++){
			if( outputs.PopTot[k1] < IonFilterThreshold){
				for(k2 = k1; k2 <= Nt; k2++){outputs.sourcetermfiltered[k2]=0.0;}
				break;
			}
		}
	}

	// PRINT field and source terms in both domains // ADD switch to do one of them or both of them
	


	// case 1:
	// 	file1 = fopen("results/tgrid.bin","wb"); file2 = fopen("results/Efield.bin","wb"); file3 = fopen("results/SourceTerm.bin","wb");
	// 	file4 = fopen("results/omegagrid.bin","wb"); file5 = fopen("results/FEfield.bin","wb"); file6 = fopen("results/FSourceTerm.bin","wb");
	// 	file7 = fopen("results/Spectrum2Efield.bin","wb"); file8 = fopen("results/Spectrum2SourceTerm.bin","wb"); file9 = fopen("results/GridDimensionsForBinaries.dat","w");
	
	
	// 	print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourceterm, (Nt+1), dt, outputs.tgrid[Nt]);
	// 	fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);

	// 	file1 = fopen("results/GS_population.bin","wb");
	// 	fwrite(outputs.PopTot,sizeof(double),(Nt+1),file1);
	// 	fclose(file1);

    //             if(IonisationFilterForTheSourceTerm == 1){
	// 	file1 = fopen("results/tmp1.bin","wb"); file2 = fopen("results/tmp2.bin","wb"); file3 = fopen("results/SourceTermFiltered.bin","wb"); // We just use the function as it is and remove redundant files... not optimal
	// 	file4 = fopen("results/tmp3.bin","wb"); file5 = fopen("results/tmp4.bin","wb"); file6 = fopen("results/FSourceTermFiltered.bin","wb");
	// 	file7 = fopen("results/tmp5.bin","wb"); file8 = fopen("results/Spectrum2SourceTermFiltered.bin","wb"); file9 = fopen("results/tmp1.dat","w");
	// 	print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourcetermfiltered, (Nt+1), dt, outputs.tgrid[Nt]);
	// 	fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);
	// 	dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp2.bin"); dumint=remove("results/tmp3.bin"); dumint=remove("results/tmp4.bin");
	// 	dumint=remove("results/tmp5.bin"); dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp1.dat");
	// 	}

	// TO COMPARE
		file1 = fopen("results/TimeDomain.dat" , "w"); file2 = fopen("results/OmegaDomain.dat" , "w");
		print2FFTW3(file1, file2, outputs.Efield, outputs.sourceterm, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2);
		file1 = fopen("results/GS_population.dat" , "w");
		for(k1 = 0; k1 <= Nt; k1++){fprintf(file1,"%e\t%e\n", outputs.tgrid[k1] , outputs.PopTot[k1]);}
		fclose(file1);

                if(IonisationFilterForTheSourceTerm == 1){
		file1 = fopen("results/TimeDomainFiltered.dat" , "w"); file2 = fopen("results/OmegaDomainFiltered.dat" , "w");
		print2FFTW3(file1, file2, outputs.Efield, outputs.sourcetermfiltered, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2);
		}


	// SAVE THE RESULTS
	calc2FFTW3(outputs.Nt, dt, tmax, outputs.Efield, outputs.sourceterm, &outputs.tgrid_fftw, &outputs.omegagrid, &outputs.FEfield,
				&outputs.Fsourceterm, &outputs.FEfieldM2, &outputs.FsourcetermM2, &outputs.Nomega); //takes real signal speced by given "dt" and it computes and prints its FFTW3

	printf("fftw_computed in single\n");  fflush(NULL);


	// print Gabor and partial spectra
	if (PrintGaborAndSpectrum == 1){
	file1 = fopen("results/OmegaDipolewindow1.dat" , "w"); 
	printlimitedFFTW3(file1, outputs.sourceterm, (Nt+1), dt, tmin1window, tmax1window);
	fclose(file1);

	file1 = fopen("results/OmegaDipolewindow2.dat" , "w"); 
	printlimitedFFTW3(file1, outputs.sourceterm, (Nt+1), dt, tmin2window, tmax2window);
	fclose(file1);

	file1 = fopen("results/GaborDimensions.dat" , "w"); file2 = fopen("results/GaborDipole_tgrid.bin" , "wb"); file3 = fopen("results/GaborDipole_omegagrid.bin" , "wb"); file4 = fopen("results/GaborDipole.bin" , "wb");
	printGaborFFTW3binary(file1, file2, file3, file4, outputs.sourceterm, (Nt+1), dt, dtGabor, a_Gabor, omegaMaxGabor);
	fclose(file1); fclose(file2); fclose(file3); fclose(file4);
	}


	//  Remove the Ground state from Psi
	// psi_rmv_gs = calloc(2*(num_r+1),sizeof(double));
	// psi_rmv_gs = rmv_gs(psi0,psi,x,num_r);




	finish = clock();


	
	//fclose(eingenvectorf);

	
	//printf("\n");
	// printf("Duration of calculation %f sec\n",(double)(finish - start) / CLOCKS_PER_SEC);
	//printf("\n");
 
  //printf("%e \n",outputs.Efield[0]);
 // printf("address %p \n",outputs.Efield);
	
	//printf("Calculation terminated ; good analysis\n");

	
return outputs;

/*
free(psi); free(psi_rmv_gs); free(psi0); free(x); free(off_diagonal);
free(diagonal); free(vector);
free(t); free(timet); free(dipole); 

fclose(eingenvaluef); fclose(eingenvectorf); fclose(pot); 
*/

}



void Initialise(int num_r)
{
    double xmax = 0.5*num_r*dx;
	x = calloc((num_r+1),sizeof(double));
	off_diagonal = calloc(2*(num_r+1),sizeof(double));
	diagonal = calloc(2*(num_r+1),sizeof(double));	

	//Initialisation Matrix corresponding to D2
	for(i=0;i<=num_r;i++)
	{
		x[i] = i*dx-xmax; 
		off_diagonal[2*i] = -0.5/(dx*dx); off_diagonal[2*i + 1] = 0.;
		diagonal[2*i] = 1./(dx*dx); diagonal[2*i + 1] = 0.;
	}

	
}




void gaussian(void)
{

 double phi1,phi2,phi3,mod1,mod2,psigaussian_re,psigaussian_im;

 

 gaussianwp = fopen("results/gaussianwp.dat", "w" );

for(i=0;i<=num_r;i++)
{
  phi1 = v*x[i] - 0.5*v*v*tfinal;
  
  phi2 = 0.125*tfinal*pow(pow(alpha,4)+tfinal*tfinal/4.,-1);
  phi2=phi2*(x[i]-v*tfinal)*(x[i]-v*tfinal);

  phi3 = -0.5*atan(0.5*tfinal/(alpha*alpha));

  mod1 = -0.25*alpha*alpha*pow(pow(alpha,4)+tfinal*tfinal/4.,-1);
  mod1 = mod1*pow(x[i]-v*tfinal,2);
  mod2 = -0.25*log(pow(alpha,4)+tfinal*tfinal/4);

  psigaussian_re = norm_gauss*sqrt(Pi)*exp(mod1+mod2)*cos(phi1+phi2+phi3);
  psigaussian_im = norm_gauss*sqrt(Pi)*exp(mod1+mod2)*sin(phi1+phi2+phi3);

  fprintf(gaussianwp,"%f\t%f\t%f\n",x[i],psigaussian_re,psigaussian_im);
}


 fclose(gaussianwp);

}



void volkov_state(void)
{
 
double phi1,phi2,phi3,phi4,mod1,mod2,psivolkov_re,psivolkov_im;
double At,intA,intA2;
double beta_re,beta_im,num; 


volkovwp = fopen("volkovwp.dat", "w" );


At = E0/omega*(cos(omega*tfinal)-1);
intA = E0/omega*(sin(omega*tfinal)/omega-tfinal);
intA2 = (E0/omega)*(E0/omega)*(1.5*tfinal+0.25*sin(2*omega*tfinal)/omega-2*sin(omega*tfinal)/omega);


for(i=0;i<=num_r;i++)
{

  phi1 = -x[i]*At;
  
  phi2 = -0.5*intA2;

  phi3 = -0.5*atan(0.5*tfinal/(alpha*alpha));

  beta_re = 2*alpha*alpha*v;
  beta_im = x[i]+intA;

  num = pow(alpha,4)+tfinal*tfinal/4.;	

  phi4 = 0.25*(2*alpha*alpha*beta_re*beta_im-0.5*tfinal*(beta_re*beta_re-beta_im*beta_im))/num;

  mod1 = 0.25*(alpha*alpha*(beta_re*beta_re-beta_im*beta_im)+tfinal*beta_im*beta_re)/num;

  mod2 = -0.25*log(num);

  psivolkov_re = norm_gauss*sqrt(Pi)*exp(mod1+mod2)*cos(phi1+phi2+phi3+phi4)*exp(-v*v*alpha*alpha);
  psivolkov_im = norm_gauss*sqrt(Pi)*exp(mod1+mod2)*sin(phi1+phi2+phi3+phi4)*exp(-v*v*alpha*alpha);

  fprintf(volkovwp,"%f\t%f\t%f\n",x[i],psivolkov_re,psivolkov_im);
}


 fclose(volkovwp);

}

void volkov_state_vg(void)
{
 
double phi1,phi2,phi3,mod1,mod2,mod3,psivolkov_re,psivolkov_im,xp;
double intA;
double num; 


volkovwp = fopen("volkovwp_vg.dat", "w" );

printf("TIME FOR VOLKOV : %f \n",tfinal);

intA = E0/omega*(-cos(omega*tfinal)/omega+1./omega);

norm_gauss = pow(2*alpha/Pi,0.25);

for(i=0;i<=num_r;i++)
{
  xp = x[i] + intA;
  num = pow(alpha,4)+tfinal*tfinal/4.;	

  mod1 = (4*pow(alpha,4)*v*v - pow(xp,2.))*alpha*alpha;
  mod1 = 0.25*mod1/num;

  mod2 = 2*xp*alpha*alpha*tfinal*v;
  mod2 = 0.25*mod2/num;

  phi1 = 4*pow(alpha,4)*xp*v;
  phi1 = 0.25*phi1/num;

  phi2 = -0.5*tfinal*(4*pow(alpha,4)*v*v-xp*xp);
  phi2 = 0.25*phi2/num;

  mod3 = -0.25*log(num);	
  phi3 = -0.5*atan(0.5*tfinal/(alpha*alpha));

  psivolkov_re = norm_gauss*sqrt(Pi)*exp(mod1+mod2+mod3)*cos(phi1+phi2+phi3)*exp(-v*v*alpha*alpha);
  psivolkov_im = norm_gauss*sqrt(Pi)*exp(mod1+mod2+mod3)*sin(phi1+phi2+phi3)*exp(-v*v*alpha*alpha);


  psi0[2*i] = psivolkov_re;
  psi0[2*i + 1] = psivolkov_im;


  fprintf(volkovwp,"%f\t%f\t%f\n",x[i],psivolkov_re,psivolkov_im);
}

 fclose(volkovwp);

}


double* rmv_gs(double *psi0,double *psi, double *x, double num_r)
{

  double *psi_new;
  double c_gs_re,c_gs_im;
  int j;

  psi_new = calloc(2*(num_r+1),sizeof(double));  

  c_gs_re = 0.0; c_gs_im = 0.0;
  for(j = 0 ; j<= num_r ; j++) 
  {
     c_gs_re += psi[2*j]*psi0[2*j];
     c_gs_im += psi[2*j+1]*psi0[2*j];
  }

  for(j = 0 ; j<= num_r ; j++) 
  {
     psi_new[2*j] =  psi[2*j] - c_gs_re*psi0[2*j];
     psi_new[2*j+1] =  psi[2*j+1] - c_gs_im*psi0[2*j];
  }

  return psi_new;

  free(psi_new);

}
