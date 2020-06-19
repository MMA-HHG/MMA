#include<time.h> 
#include<stdio.h>
#include<stdlib.h>
#include<malloc.h>
#include<math.h>
#include "hdf5.h"



#include"util.h"
#include"util2.h"

struct Efield_var Efield;
struct trg_def trg;
struct analy_def analy;
struct outputs_def outputs;
double *y,*x,*a,*c,sum,*diagonal,*off_diagonal,*eigenvector,*u,*r,*vector;
double *timet,*dipole;
double dx,xmax,Eguess,Einit,CV,phi,omega,E0,period,Pi,tfinal,alpha,v,mod1,mod2,dE,Estep,norm_gauss,x_int,textend;

int gauge,transformgauge,fieldinau,input0,Ntinterp,InterpByDTorNT;

double dt, tmax,tmin;
int Nt;

double dum,dum1,dum2;
double *psi0,*psi,*psi2,Einit2,ps_re,ps_im,*psiexc,*psi_rmv_gs;
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

//#pragma warning( disable : 4996 ) // warning for fopen in visual 2005



int main(void)
{	
	Pi = acos(-1.);

	// TESTING VARIABLES
	int kz = 0;
	int kr = 0;

	// Open the param.txt file for intialisation of the parameter
	param = fopen("param.txt" , "r");
	if(param == NULL) {printf("DATA could not be found in param.txt file\n");}

	dumint=fscanf(param, "%*[^\n]\n", NULL);
	dumint=fscanf(param,"%i %*[^\n]\n",&Efield.fieldtype); // 0-numerical, loaded in femtosecons, 1-numerical, loaded in atomic units in whole grid, 2-analytical


	dumint=fscanf(param, "%*[^\n]\n", NULL);

	dumint=fscanf(param,"%lf %*[^\n]\n",&Eguess); // Energy of the initial state
	dumint=fscanf(param,"%i %*[^\n]\n",&num_r); // Number of points of the initial spatial grid 16000
	dumint=fscanf(param,"%i %*[^\n]\n",&num_exp); // Number of points of the spatial grid for the expansion
	dumint=fscanf(param,"%lf %*[^\n]\n",&dx); // resolution for the grid
	dumint=fscanf(param,"%i %*[^\n]\n",&InterpByDTorNT); // Number of points of the spatial grid for the expansion
	dumint=fscanf(param,"%lf %*[^\n]\n",&dt); // resolution in time
	dumint=fscanf(param,"%i %*[^\n]\n",&Ntinterp); // Number of points of the spatial grid for the expansion
	dumint=fscanf(param,"%lf %*[^\n]\n",&textend); // extension of the calculation after the last fields ends !!! NOW ONLY FOR ANALYTICAL FIELD //700
	dumint=fscanf(param,"%i %*[^\n]\n",&analy.writewft); // writewavefunction (1-writting every tprint)
	dumint=fscanf(param,"%lf %*[^\n]\n",&analy.tprint); // time spacing for writing the wavefunction	
	dumint=fscanf(param,"%lf %*[^\n]\n",&x_int); // the limit of the integral for the ionisation //2 2 works fine with the lenth gauge and strong fields
	dumint=fscanf(param,"%i %*[^\n]\n",&PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)
	dumint=fscanf(param,"%lf %*[^\n]\n",&a_Gabor); // the parameter of the gabor window [a.u.]
	dumint=fscanf(param,"%lf %*[^\n]\n",&omegaMaxGabor); // maximal frequency in Gabor [a.u.]
	dumint=fscanf(param,"%lf %*[^\n]\n",&dtGabor); // spacing in Gabor
	dumint=fscanf(param,"%lf %*[^\n]\n",&tmin1window); // analyse 1st part of the dipole
	dumint=fscanf(param,"%lf %*[^\n]\n",&tmax1window); // analyse 1st part of the dipole
	dumint=fscanf(param,"%lf %*[^\n]\n",&tmin2window); // analyse 2nd part of the dipole
	dumint=fscanf(param,"%lf %*[^\n]\n",&tmax2window); // analyse 2nd part of the dipole
	dumint=fscanf(param,"%i %*[^\n]\n",&PrintOutputMethod); // (0 - only text, 1 - only binaries, 2 - both)
	dumint=fscanf(param,"%i %*[^\n]\n",&IonisationFilterForTheSourceTerm); // filter source term by high-ionisation components (1-yes)
	dumint=fscanf(param,"%lf %*[^\n]\n",&IonFilterThreshold); // threshold for the ionisation [-]

	dumint=fscanf(param, "%*[^\n]\n", NULL);

	dumint=fscanf(param,"%lf %*[^\n]\n",&trg.a); // the limit of the integral for the ionisation //2 2 works fine with the

	dumint=fscanf(param, "%*[^\n]\n", NULL);


	// paramters either for numerical field or analytic
	// we extend it for now as new functionality and then resort it
	switch (Efield.fieldtype){
	case 3:
		// HDF5, filename given directly when loading
	break;
	case 2:
		// GAUGES
		dumint=fscanf(param,"%i %*[^\n]\n",&gauge); // 0-length, otherwise velocity, velocity available only for analytic field (A needed)
		dumint=fscanf(param,"%i %*[^\n]\n",&transformgauge); // 1 - transform also to another gauge during the calculation, (A needed)
/*		printf("gauge,  %i \n",gauge);*/
	break;
	case 0:	case 1:
		// FILENAMES

		dumint=fscanf(param,"%s %*[^\n]\n",filename1); // filename1
		dumint=fscanf(param,"%s %*[^\n]\n",filename2); // filename2
		dumint=fscanf(param, "%*[^\n]\n", NULL);
		dumint=fscanf(param,"%i %*[^\n]\n",&fieldinau); // 0-input inatomic units, 1 - in femto and GV/m
		dumint=fscanf(param,"%i %*[^\n]\n",&input0); // 0 field starts with 0, 1 in the middle of the file

		printf("filename1  %s \n",filename1);
	break;

	}

	

	num_t = floor((2.*Pi)/(0.057*dt)); // the length of one cycle for 800 nm (i.e. omega=0.057) 
	

	// define the properties of the temporal grid
	switch (Efield.fieldtype){
	case 2:
		printf("Analytical fields are used\n");
/*		dumint=fscanf(param, "%*[^\n]\n", NULL); // move in file*/
/*		dumint=fscanf(param,"%lf %*[^\n]\n",&dum); printf("E0,  %f \n",dum);*/

		define_analytical(&Efield, param);
	break;
	case 0:
		printf("Numerical field 1 (in femtoseconds) is used\n");
	break;
	case 1:
		printf("Numerical field 2 (in atomic units) is used\n");
	break;

	}

	fclose(param);


	// !!!!! loading procedure
	switch (Efield.fieldtype){

	case 3:
		// OPEN HDF5 file	
        file_id = H5Fopen ("results.h5", H5F_ACC_RDWR, H5P_DEFAULT); // open file

		// find dimensions	
		dataset_id = H5Dopen2 (file_id, "IRprop/tgrid ", H5P_DEFAULT); // open dataset	     
        dspace_id=h5dget_space(dset_id); // Get the dataspace ID     
        h5sget_simple_extent_dims(dspace_id, dims, maxdims)  //Getting dims from dataspace

		// allocate fields

		// set hyperslab
		dataset_id = H5Dopen2 (file_id, "IRprop/Fields_rzt", H5P_DEFAULT); // open dataset

		// load fields

		// close field
		
		Efield.Nt = 0;	
		while((ch=fgetc(file1))!=EOF)
		{if(ch == '\n'){Efield.Nt++;}}
		rewind(file1);
		printf("Efield.Nt,  %i \n",Efield.Nt);		
			
		Efield.tgrid = calloc(Efield.Nt,sizeof(double)); 
		Efield.Field = calloc(Efield.Nt,sizeof(double)); 

		// LOAD FILES		
		for(k1 = 0 ; k1 <= Efield.Nt-1 ; k1++) //Efield.Nt-1
		{
			dumint = fscanf(file1,"%lf",&Efield.tgrid[k1]); // file in femtoseconts and GV/m;
			dumint = fscanf(file2,"%lf",&Efield.Field[k1]);
			if (!(fieldinau == 0))
			{
				Efield.tgrid[k1] = Efield.tgrid[k1]*41.34144728; // timegrid in a.u.
				Efield.Field[k1] = Efield.Field[k1]*0.001944689151; // corresponding field
			}
		}	
		fclose(file1);
		fclose(file2);		

		// k1 = 0; k2 = 0;	findinterval(Efield.Nt, 0., Efield.tgrid, &k1, &k2);// find zero of the grid, the best resolution is around
		switch ( input0 ){case 0: dumint = 0; break; case 1: dumint = round(Efield.Nt/2.); /* field centered around 0 */ break;} // original definition
	
		Efield.dt = Efield.tgrid[dumint+1]-Efield.tgrid[dumint]; // Efield.dt = Efield.tgrid[1+round(Efield.Nt/2.)]-Efield.tgrid[round(Efield.Nt/2.)];
		tmax = Efield.tgrid[Efield.Nt]-Efield.tgrid[0];

		// PRINT field and its transform
		file1 = fopen("inputs/InputField.dat" , "w"); file2 = fopen("inputs/InputFField.dat" , "w"); printFFTW3(file1, file2, Efield.Field, Efield.Nt, Efield.dt); fclose(file1); fclose(file2);
		
		// find the padding we need (now we assume coarser grid in the input, need an "if" otherwise)
		if (InterpByDTorNT == 1)
		{
			k1 = Ntinterp + 1;
		} else {
			k1 = floor(Efield.dt/dt); Ntinterp = k1; k1++;
		}

		dt = Efield.dt/((double)k1); // redefine dt properly

		Efield.Field = FourInterp(k1, Efield.Field, Efield.Nt); // make the interpolation !!!!!! tgrid does not correspond any more

		Nt = k1*Efield.Nt + 1;

		// PRINT new field and its transform
		file1 = fopen("inputs/InputField2.dat" , "w"); file2 = fopen("inputs/InputFField2.dat" , "w"); printFFTW3(file1, file2, Efield.Field, Nt, dt); fclose(file1); fclose(file2);

		
		num_t = floor((2*Pi)/(0.057*dt)); num_t++; 
		printf("Efield.dt,  %lf \n",Efield.dt);
		printf("points per interval  %i \n",num_t);
		printf("number of interpolated points per interval  %i \n",Ntinterp);
		
	break;

	case 0:
		// LOAD THE FILES		
		
		file1 = fopen(filename1 , "r"); if(file1 == NULL) {printf("timegrid file %s doesn't exist\n",filename1);}
		file2 = fopen(filename2 , "r"); if(file2 == NULL) {printf("Field file %s doesn't exist\n", filename2);}	
		
		
		Efield.Nt = 0;	
		while((ch=fgetc(file1))!=EOF)
		{if(ch == '\n'){Efield.Nt++;}}
		rewind(file1);
		printf("Efield.Nt,  %i \n",Efield.Nt);		
			
		Efield.tgrid = calloc(Efield.Nt,sizeof(double)); 
		Efield.Field = calloc(Efield.Nt,sizeof(double)); 

		// LOAD FILES		
		for(k1 = 0 ; k1 <= Efield.Nt-1 ; k1++) //Efield.Nt-1
		{
			dumint = fscanf(file1,"%lf",&Efield.tgrid[k1]); // file in femtoseconts and GV/m;
			dumint = fscanf(file2,"%lf",&Efield.Field[k1]);
			if (!(fieldinau == 0))
			{
				Efield.tgrid[k1] = Efield.tgrid[k1]*41.34144728; // timegrid in a.u.
				Efield.Field[k1] = Efield.Field[k1]*0.001944689151; // corresponding field
			}
		}	
		fclose(file1);
		fclose(file2);		

		// k1 = 0; k2 = 0;	findinterval(Efield.Nt, 0., Efield.tgrid, &k1, &k2);// find zero of the grid, the best resolution is around
		switch ( input0 ){case 0: dumint = 0; break; case 1: dumint = round(Efield.Nt/2.); /* field centered around 0 */ break;} // original definition
	
		Efield.dt = Efield.tgrid[dumint+1]-Efield.tgrid[dumint]; // Efield.dt = Efield.tgrid[1+round(Efield.Nt/2.)]-Efield.tgrid[round(Efield.Nt/2.)];
		tmax = Efield.tgrid[Efield.Nt]-Efield.tgrid[0];

		// PRINT field and its transform
		file1 = fopen("inputs/InputField.dat" , "w"); file2 = fopen("inputs/InputFField.dat" , "w"); printFFTW3(file1, file2, Efield.Field, Efield.Nt, Efield.dt); fclose(file1); fclose(file2);
		
		// find the padding we need (now we assume coarser grid in the input, need an "if" otherwise)
		if (InterpByDTorNT == 1)
		{
			k1 = Ntinterp + 1;
		} else {
			k1 = floor(Efield.dt/dt); Ntinterp = k1; k1++;
		}

		dt = Efield.dt/((double)k1); // redefine dt properly

		Efield.Field = FourInterp(k1, Efield.Field, Efield.Nt); // make the interpolation !!!!!! tgrid does not correspond any more

		Nt = k1*Efield.Nt + 1;

		// PRINT new field and its transform
		file1 = fopen("inputs/InputField2.dat" , "w"); file2 = fopen("inputs/InputFField2.dat" , "w"); printFFTW3(file1, file2, Efield.Field, Nt, dt); fclose(file1); fclose(file2);

		
		num_t = floor((2*Pi)/(0.057*dt)); num_t++; 
		printf("Efield.dt,  %lf \n",Efield.dt);
		printf("points per interval  %i \n",num_t);
		printf("number of interpolated points per interval  %i \n",Ntinterp);
		
	break;
	case 1:
		// scan field

		file1 = fopen("inputs/timegrid2.dat" , "r");
		file2 = fopen("inputs/Field.dat" , "r");

		Efield.Nt = 0;	
		while((ch=fgetc(file1))!=EOF)
		{if(ch == '\n'){Efield.Nt++;}}
		rewind(file1);

		Efield.tgrid = calloc(Efield.Nt,sizeof(double));
		Efield.Field = calloc(Efield.Nt,sizeof(double));
 
		for(k1 = 0 ; k1<= Efield.Nt-1 ; k1++) //Efield.Nt-1
			{
				dumint = fscanf(file1,"%lf",&Efield.tgrid[k1]); // timegrid in a.u.
				//Efield.tgrid[k1] = Efield.tgrid[k1]*41.34144728;
				dumint = fscanf(file2,"%lf",&Efield.Field[k1]); // corresponding field
				//Efield.Field[k1] = Efield.Field[k1]*0.001944689151;
				// printf("row %i\n",k1);
			}	
		fclose(file1);
		fclose(file2);		
		
		tmin = Efield.tgrid[0]; //Efield.tgrid[Efield.Nt-1]; //field is long, try to take only fewer points, given by hand
		// tmin = findnextinterpolatedzero(Efield.Nt-1, tmin, Efield.tgrid, Efield.Field); // start in nearest exact zero
		printf("tmin is,  %lf \n",tmin);
		tmax = Efield.tgrid[Efield.Nt-1];
		printf("tmax is,  %lf \n",tmax);
		Nt = floor((tmax-tmin)/dt);
		tmax = tmin+dt*(double)Nt;
		
		Efield.dt = Efield.tgrid[1]-Efield.tgrid[0]; // not used

		num_t = floor((2*Pi)/(Efield.omega*dt)); num_t++; 
		printf("Efield.dt,  %lf \n",Efield.dt);
		printf("points per interval  %i \n",num_t);	

	break;
	case 2:
		//---------------------------- THERE IS DEFINED ANALYTICAL FIELD

		printf("\ntest analytic\n");
		// finding the maximal time & fit the field specification
		//for(k1 = 0 ; k1 <= F.Nflt1 ; k1++)
		tmax = 0.;
		if (Efield.Nsin2 > 0){
			for(k1 = 0 ; k1 <= (Efield.Nsin2-1) ; k1++){
				
				Efield.sin2[k1].A0 = Efield.sin2[k1].E0/Efield.sin2[k1].o;
				Efield.sin2[k1].oc = Efield.sin2[k1].o/(2.*Efield.sin2[k1].nc);			
				Efield.sin2[k1].phi0 = Pi/2.0 - Efield.sin2[k1].oc * ( Efield.sin2[k1].ti + 0.5 * Pi/Efield.sin2[k1].oc ); //adjust of the envelope phase
				Efield.sin2[k1].phi = Efield.sin2[k1].phi - Efield.sin2[k1].o * Efield.sin2[k1].ti - Pi * Efield.sin2[k1].nc; // recalculate CEP


				tmax = fmax(tmax, Efield.sin2[k1].ti + Pi / Efield.sin2[k1].oc );
			}
		}

// ANOTHER FIELD TYPES
		if (Efield.Nflt1 > 0){for(k1 = 0 ; k1 <= (Efield.Nflt1-1) ; k1++){tmax = fmax(tmax, Efield.flt1[k1].ti + Efield.flt1[k1].ton + Efield.flt1[k1].toff + Efield.flt1[k1].T );}}
		if (Efield.Nflt1ch > 0){for(k1 = 0 ; k1 <= (Efield.Nflt1ch-1) ; k1++){tmax = fmax(tmax, Efield.flt1ch[k1].ti + Efield.flt1ch[k1].ton + Efield.flt1ch[k1].toff + Efield.flt1ch[k1].T );}}
		
		//printf("test2 \n");		
		if (Efield.NEsin2 > 0){
			
			for(k1 = 0 ; k1 <= (Efield.NEsin2-1) ; k1++){
				
				Efield.Esin2[k1].A0 = Efield.Esin2[k1].E0/Efield.Esin2[k1].o;
				Efield.Esin2[k1].oc = Efield.Esin2[k1].o/(2.*Efield.Esin2[k1].nc);			
				Efield.Esin2[k1].phi0 = Pi/2.0 - Efield.Esin2[k1].oc * ( Efield.Esin2[k1].ti + 0.5 * Pi/Efield.Esin2[k1].oc ); //adjust of the envelope phase
				Efield.Esin2[k1].phi = Efield.Esin2[k1].phi - Efield.Esin2[k1].o * Efield.Esin2[k1].ti - Pi * Efield.Esin2[k1].nc; // recalculate CEP


				tmax = fmax(tmax, Efield.Esin2[k1].ti + Pi / Efield.Esin2[k1].oc );
			}
		}

		tmax = tmax + textend;
		Nt = floor( tmax/dt ); Nt++;
		// Nt = Nt + 1000;
		/*Nt = Efield.trap.nc*(num_t+1);
		period = 2.*Pi/Efield.trap.omega;
		dt = period/num_t;
		Efield.trap.phi = Efield.trap.phi*Pi*0.5;
		*/
	break;
	}

	



	// TESTING WORKSPACE

/*	printf("Test started \n");*/
/*	double *testarray;*/
/*	dum = 5.0;*/
/*	testarray = calloc(1,sizeof(double));*/
/*	testarray[0] = dum;*/

/*	file1 = fopen("results/binaryfile1.bin" , "wb");*/
/*	fwrite(testarray,sizeof(double),1,file1);*/
/*	*/
/*	fclose(file1);*/
/*	*/
/*	printf("test2 \n");*/

/*	double *testarray2;*/
/*	testarray2 = calloc(1,sizeof(double));*/

/*	file2 = fopen("results/binaryfile1.bin" , "rb");*/
/*	fread(testarray2,sizeof(double),1,file2);*/
/*	fclose(file2);*/

/*	printf("elem1  %lf \n",testarray2[0]);*/



/*	printf("Test finished \n");*/
/*	exit(0);*/

/*	printf("Test started \n");*/
/*	double *testarray;*/
/*	testarray = calloc(3,sizeof(double));*/
/*	testarray[0] = 0.5;  testarray[2] = 1.; testarray[3] = 1.5;*/

/*	file1 = fopen("results/binaryfile1.bin" , "wb");*/
/*	fwrite(testarray,sizeof(double),3,file1);*/
/*	*/
/*	testarray[1] = 2.0;  testarray[2] = 2.5; testarray[3] = 3.0;*/
/*	fwrite(testarray,sizeof(double),3,file1);*/

/*	fclose(file1);*/
/*	*/
/*	printf("test2 \n");*/

/*	double *testarray2;*/
/*	testarray2 = calloc(6,sizeof(double));*/


/*	file2 = fopen("results/binaryfile1.bin" , "rb");*/
/*	fread(testarray2,sizeof(double),6,file2);*/
/*	fclose(file2);*/

/*	printf("elem1  %lf \n",testarray2[0]);*/

/*	printf("elem4  %lf \n",testarray2[4]);*/


/*	printf("Test finished \n");*/
/*	exit(0);*/

	
	
	

	
	//printf("Efield.tgrid[0],  %lf \n",Efield.tgrid[0]);
	//printf("Efield.tgrid[Efield.Nt-1],  %lf \n",Efield.tgrid[Efield.Nt-1]);

	// dum = (double)Nt;
	printf("Implicitly in atomic units \n");

	printf("\ntime properties \n");	
	printf("dt,  %lf \n",dt);
	//printf("points per cycle,  %i \n",num_t);
	printf("total points,  %i \n",Nt);
	// printf("Nt_old,  %i \n",nc*(num_t+1));

	/*
	printf("\nField properties \n");	
	printf("Amplitude,  %lf \n",Efield.trap.E0);
	printf("omega0,  %lf \n",Efield.trap.omega);
	*/

	printf("\nspace properties\n");	
	printf("dx,  %lf \n",dx);
	printf("intial nmax : %i \n",num_r);
	printf("nextend : %i \n",num_exp);	
	printf("initial xmax : %lf \n",num_r*dx/2.);

	printf("\n\n");	

	size = 2*(num_r+1);// for complex number


	// Allocation memory

	printf("test1 %i \n",size);

	x = calloc((num_r+1),sizeof(double));
	printf("test2 \n"); 
	off_diagonal = calloc(size,sizeof(double));
	printf("test3 \n"); 
	diagonal = calloc(size,sizeof(double));
	vector = calloc(size,sizeof(double));
	psi0 = calloc(size,sizeof(double));
	psi2 = calloc(size,sizeof(double));
	psi = calloc(size,sizeof(double));
	psiexc = calloc(size,sizeof(double));
	
	t = calloc(Nt,sizeof(double));

	printf("test4 \n");

	timet = calloc(Nt,sizeof(double));
	dipole = calloc(2*Nt,sizeof(double));

	eingenvaluef = fopen("results/eingenvalue.dat", "w" );
	eingenvectorf = fopen("results/eingenvector.dat", "w" );
	pot = fopen("results/latice.dat", "w" );

	if(eingenvaluef==NULL || eingenvectorf==NULL)
	{printf("pb d'ouverture de fichier");}
	

	printf("Initialisation \n");

	// Initialise vectors and Matrix 
	Initialise(num_r);
	for(i=0;i<=num_r;i++){psi0[2*i] = 1.0; psi0[2*i+1] = 0.; psiexc[2*i] = 1; psiexc[2*i+1] = 0.;}
	//normalise(psi0,num_r); // Initialise psi0 for Einitialise
	//normalise(psiexc,num_r);

	CV = 1E-20; // CV criteria

	/* This number has to be small enough to assure a good convregence of the wavefunction
	if it is not the case, then the saclar product of the the ground state and the excited states 
	is not quite 0 and those excited appears in the energy analysis of the gorund states, so the propagation !!
	CV = 1E-25 has been choosen to have a scalar product of 10^-31 with the third excited state for num_r = 5000 and dx=0.1
	*/

	printf("Calculation of the energy of the ground sate ; Eguess : %f\n",Eguess);

	Einit = Einitialise(trg,psi0,off_diagonal,diagonal,off_diagonal,x,Eguess,CV,num_r);
	for(i=0;i<=num_r;i++) {fprintf(eingenvectorf,"%f\t%e\t%e\n",x[i],psi0[2*i],psi0[2*i+1]); fprintf(pot,"%f\t%e\n",x[i],potential(x[i],trg));}

	 
	printf("Initial energy is : %1.12f\n",Einit);
	printf("first excited energy is : %1.12f\n",Einit2);

	printf("\n");	
	printf("Propagation procedure ...\n");
	printf("\n");	


	outputs.tgrid = calloc((Nt+1),sizeof(double)); outputs.Efield = calloc((Nt+1),sizeof(double)); outputs.sourceterm = calloc((Nt+1),sizeof(double)); outputs.PopTot = calloc((Nt+1),sizeof(double));

	start = clock();


	psi = propagation(trg,Efield,tmin,Nt,num_t,dt,num_r,num_exp,dx,psi0,psi,x,timef,timef2,ton,toff,timet,dipole,gauge,transformgauge,x_int,analy,outputs);

/*	printf("\ntmax test\n");	*/
/*	printf("tmax,  %lf \n",*outputs.tmax);*/

//	volkov_state_vg();

	
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
	

	switch (PrintOutputMethod){
	printf("\nPrintingResults\n");
	case 0:
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
	break;
	case 1:
		file1 = fopen("results/tgrid.bin","wb"); file2 = fopen("results/Efield.bin","wb"); file3 = fopen("results/SourceTerm.bin","wb");
		file4 = fopen("results/omegagrid.bin","wb"); file5 = fopen("results/FEfield.bin","wb"); file6 = fopen("results/FSourceTerm.bin","wb");
		file7 = fopen("results/Spectrum2Efield.bin","wb"); file8 = fopen("results/Spectrum2SourceTerm.bin","wb"); file9 = fopen("results/GridDimensionsForBinaries.dat","w");
		print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourceterm, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);

		file1 = fopen("results/GS_population.bin","wb");
		fwrite(outputs.PopTot,sizeof(double),(Nt+1),file1);
		fclose(file1);

                if(IonisationFilterForTheSourceTerm == 1){
		file1 = fopen("results/tmp1.bin","wb"); file2 = fopen("results/tmp2.bin","wb"); file3 = fopen("results/SourceTermFiltered.bin","wb"); // We just use the function as it is and remove redundant files... not optimal
		file4 = fopen("results/tmp3.bin","wb"); file5 = fopen("results/tmp4.bin","wb"); file6 = fopen("results/FSourceTermFiltered.bin","wb");
		file7 = fopen("results/tmp5.bin","wb"); file8 = fopen("results/Spectrum2SourceTermFiltered.bin","wb"); file9 = fopen("results/tmp1.dat","w");
		print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourcetermfiltered, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);
		dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp2.bin"); dumint=remove("results/tmp3.bin"); dumint=remove("results/tmp4.bin");
		dumint=remove("results/tmp5.bin"); dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp1.dat");
		}
	break;
	case 2:
		file1 = fopen("results/tgrid.bin","wb"); file2 = fopen("results/Efield.bin","wb"); file3 = fopen("results/SourceTerm.bin","wb");
		file4 = fopen("results/omegagrid.bin","wb"); file5 = fopen("results/FEfield.bin","wb"); file6 = fopen("results/FSourceTerm.bin","wb");
		file7 = fopen("results/Spectrum2Efield.bin","wb"); file8 = fopen("results/Spectrum2SourceTerm.bin","wb"); file9 = fopen("results/GridDimensionsForBinaries.dat","w");
		print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourceterm, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);

		file1 = fopen("results/GS_population.bin","wb");
		fwrite(outputs.PopTot,sizeof(double),(Nt+1),file1);
		fclose(file1);

                if(IonisationFilterForTheSourceTerm == 1){
		file1 = fopen("results/tmp1.bin","wb"); file2 = fopen("results/tmp2.bin","wb"); file3 = fopen("results/SourceTermFiltered.bin","wb"); // We just use the function as it is and remove redundant files... not optimal
		file4 = fopen("results/tmp3.bin","wb"); file5 = fopen("results/tmp4.bin","wb"); file6 = fopen("results/FSourceTermFiltered.bin","wb");
		file7 = fopen("results/tmp5.bin","wb"); file8 = fopen("results/Spectrum2SourceTermFiltered.bin","wb"); file9 = fopen("results/tmp1.dat","w");
		print2FFTW3binary(file1, file2, file3, file4, file5, file6, file7, file8, file9, outputs.Efield, outputs.sourcetermfiltered, (Nt+1), dt, outputs.tgrid[Nt]);
		fclose(file1); fclose(file2); fclose(file3); fclose(file4); fclose(file5); fclose(file6); fclose(file7); fclose(file8); fclose(file9);
		dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp2.bin"); dumint=remove("results/tmp3.bin"); dumint=remove("results/tmp4.bin");
		dumint=remove("results/tmp5.bin"); dumint=remove("results/tmp1.bin"); dumint=remove("results/tmp1.dat");
		}

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
	break;
	}




	printf("\n");
	printf("Calculation of the HHG spectrum \n");


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

/*	if (PrintInputGabor==1){*/
/*	file1 = fopen("results/GaborInput.dat" , "w"); file2 = fopen("results/GaborInput_tgrid.dat" , "w"); file3 = fopen("results/GaborInput_omegagrid.dat" , "w"); file4 = fopen("results/GaborDipole.bin" , "wb");*/
/*	printGaborFFTW3(file1, file2, file3, file4, outputs.sourceterm, (Nt+1), dt, dtGabor, a_Gabor, omegaMaxGabor);*/
/*	fclose(file1); fclose(file2); fclose(file3); fclose(file4);		*/
/*	}*/


/*	dumint=fscanf(param,"%i %*[^\n]\n",&PrintGaborAndSpectrum); // print Gabor and partial spectra (1-yes)*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&a_Gabor); // the parameter of the gabor window [a.u.]*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&omegaMaxGabor); // maximal frequency in Gabor*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&dtGabor); // spacing in Gabor*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&tmin1window); // analyse 1st part of the dipole*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&tmax1window); // analyse 1st part of the dipole*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&tmin2window); // analyse 2nd part of the dipole*/
/*	dumint=fscanf(param,"%lf %*[^\n]\n",&tmax2window); // analyse 2nd part of the dipole*/


/// howto write binary

/*	printf("Test started \n");*/
/*	double *testarray;*/
/*	testarray = calloc(3,sizeof(double));*/
/*	testarray[1] = 0.5;  testarray[2] = 1.; testarray[3] = 1.5;*/

/*	file1 = fopen("results/binaryfile1.bin" , "wb");*/
/*	fwrite(testarray,sizeof(double),3,file1);*/
/*	*/
/*	testarray[1] = 2.0;  testarray[2] = 2.5; testarray[3] = 3.0;*/
/*	fwrite(testarray,sizeof(double),3,file1);*/

/*	fclose(file1);*/
/*	*/
/*	printf("test2 \n");*/

/*	double *testarray2;*/
/*	testarray2 = calloc(6,sizeof(double));*/

/*	file2 = fopen("results/binaryfile1.bin" , "rb");*/
/*	fread(testarray2,sizeof(double),6,file2);*/
/*	fclose(file2);*/

/*	printf("elem1  %lf \n",testarray2[1]);*/

/*	printf("elem4  %lf \n",testarray2[4]);*/


/*	printf("Test finished \n");*/
/*	exit(0);*/




	// Remove the Ground state from Psi
	psi_rmv_gs = calloc(2*(num_r+1),sizeof(double));
	psi_rmv_gs = rmv_gs(psi0,psi,x,num_r);




	finish = clock();


	
	//fclose(eingenvectorf);

	
	printf("\n");
	printf("Duration of calculation %f sec\n",(double)(finish - start) / CLOCKS_PER_SEC);
	printf("\n");
	
	printf("Calculation terminated ; good analysis\n");

	
return 0;

/*
free(psi); free(psi_rmv_gs); free(psi0); free(x); free(off_diagonal);
free(diagonal); free(vector); free(psi2); free(psiexc); 
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







































	/*case 3:

		// pulses boundaries
		Efield.sin2.tmax1 = 2.*Pi*Efield.sin2.nc1/Efield.sin2.o1;
		Efield.sin2.tmin2 = Efield.sin2.delay;
		Efield.sin2.tmax2 = Efield.sin2.delay + 2.*Pi*Efield.sin2.nc2/Efield.sin2.o2;


		// envelope frequencies and phases
		Efield.sin2.oc1 = Efield.sin2.o1/(2.0*Efield.sin2.nc1); // 1st envelope frequency
		Efield.sin2.oc2 = Efield.sin2.o2/(2.0*Efield.sin2.nc2); // 2nd envelope frequency
		Efield.sin2.Tc2 = 2.0*Pi*Efield.sin2.nc2 / Efield.sin2.o2; //2nd pulse duration

		Efield.sin2.phi0 = Pi/2.0 - Efield.sin2.oc2*(Efield.sin2.delay + 0.5*Efield.sin2.Tc2); //adjust of the envelope phase

		// adjust of the field phases
		Efield.sin2.phi1 = Efield.sin2.phi1 - Efield.sin2.nc1*Pi; // 1st field phase adjusted for sin^2 pulse



		// Efield.sin2.phi2 = Efield.sin2.phi2 - Efield.sin2.o2*(Efield.sin2.delay + 0.5*Efield.sin2.Tc2);	// 2nd field phase adjusted for the delayed sin^2 pulse
		Efield.sin2.phi2 = Efield.sin2.phi1; // common phase for both pulses




		if (Efield.sin2.tmax1 >= Efield.sin2.tmax2){Efield.sin2.tmax = Efield.sin2.tmax1;}
		else{Efield.sin2.tmax = Efield.sin2.tmax2;}

		// THREE DIFFERENT CASES FOR THE OVERLAP
		if (Efield.sin2.tmax1<=Efield.sin2.tmin2){Efield.sin2.overlap = 1;}
		else if ( ( Efield.sin2.tmin2<Efield.sin2.tmax1 ) && ( Efield.sin2.tmax1<Efield.sin2.tmax2) ){Efield.sin2.overlap = 2;}
		else if ( (Efield.sin2.tmin2<=Efield.sin2.tmax1) && ( Efield.sin2.tmax2<=Efield.sin2.tmax1) ){Efield.sin2.overlap = 3;}	

		// printf("overlap %i\n", Efield.sin2.overlap);	
		printf("tmax %lf\n", Efield.sin2.tmax);	

		dum = (Efield.sin2.tmax*Efield.sin2.o1)/(2.*Pi); // virtual number of cycles according to the first field frequency
		Nt = floor( dum*(num_t+1) ); Nt++;

		period = 2.*Pi/Efield.sin2.o1;

		dt = period/num_t;	


	break;*/

