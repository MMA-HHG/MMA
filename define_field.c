#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>

#include"util.h"

#define Pi acos(-1.)
clock_t start, finish;
clock_t start2, finish2;
#pragma warning( disable : 4996 ) // warning for fopen in visual 2005

extern double* timet,dipole;

// extern struct Efield_var;


void define_analytical(struct Efield_var *Efield, FILE *param)
{

double E0, dum, dum1, dum2;
int dumint;

// defined in vector potential

		// FLAT TOPS

		// be careful if  omega = 2*omega_envelope, problem in analytic expression
		

		// SIN^2 IN VECTOR POTENTIAL PULSES (SATISFYING A(tf) = 0 by definition)
		printf("test1\n");
		dumint=fscanf(param, "%*[^\n]\n", NULL); // move in file
		printf("test2\n");
		Efield->Nsin2 = 1;
		if( Efield->Nsin2 > 0){
			Efield->sin2 = calloc(Efield->Nsin2,sizeof(struct sin2_definition));
			
			// 1st field
			dumint=fscanf(param,"%lf %*[^\n]\n",&Efield->sin2[0].E0);
			dumint=fscanf(param,"%lf %*[^\n]\n",&Efield->sin2[0].o);
			dumint=fscanf(param,"%lf %*[^\n]\n",&Efield->sin2[0].ti); 
			dumint=fscanf(param,"%lf %*[^\n]\n",&Efield->sin2[0].nc);  // # of cycles
			dumint=fscanf(param,"%lf %*[^\n]\n",&Efield->sin2[0].phi); // CEP [in radians, reference is a cosine pulse in A]
		}
}
