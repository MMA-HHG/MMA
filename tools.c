#include<math.h>
#include<malloc.h>
#include<stdlib.h>
#include<stdio.h>
#include<time.h>
#include<fftw3.h>
#include "hdf5.h"
#include "mpi.h"

#include"util.h"

clock_t start, finish;
clock_t start2, finish2;

extern double* timet,dipole;

// extern struct Efield_var;

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

void PrintOutputs(hid_t file_id, struct inputs_def *in, struct outputs_def *out) // Initialise ground-state
{
	if ( (*in).Print.Efield == 1 ){}
	res.Efield = 0;
	res.FEfield = 0;
	res.sourceterm = 0;
	res.Fsourceterm = 0;
	res.FEfieldM2 = 0;
	res.FsourceTermM2 = 0;
	res.PopTot = 0;
	res.tgrid = 0;
	res.omegagrid = 0;
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



double AField(struct Efield_var F, double t) // ANAlytic field is -dA/dt, the sign!
{
	
	double omegap = F.trap.omega/((double)F.trap.nc),ts,Tf,A0;

	double a,b;
	double A,A1,A2,dum;
	int k1;

switch (F.fieldtype){

	case 2:  // analytic 
		A = 0.;
		//for(k1 = 0 ; k1 <= F.Nflt1 ; k1++)
		if (F.Nflt1 > 0)
		{
			for(k1 = 0 ; k1 <= (F.Nflt1-1) ; k1++)
			{
				A = A + Afieldflattop1(t, F.flt1[k1].ti, F.flt1[k1].ton, F.flt1[k1].toff, F.flt1[k1].T, F.flt1[k1].o, F.flt1[k1].phi, F.flt1[k1].A);
			}
		}

		if (F.Nsin2 > 0)
		{
			for(k1 = 0 ; k1 <= (F.Nsin2-1) ; k1++)
			{
				// Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
				A = A + Afieldsin2(t , F.sin2[k1].ti , F.sin2[k1].A0 , F.sin2[k1].oc , F.sin2[k1].phi0 , F.sin2[k1].o , F.sin2[k1].phi);
			}
		}

		if (F.NEsin2 > 0)
		{
			for(k1 = 0 ; k1 <= (F.NEsin2-1) ; k1++)
			{
				// Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
				A = A + AfieldEsin2(t , F.Esin2[k1].ti , F.Esin2[k1].A0 , F.Esin2[k1].oc , F.Esin2[k1].phi0 , F.Esin2[k1].o , F.Esin2[k1].phi);
			}
		}

		if (F.Nflt1ch > 0)
		{
			for(k1 = 0 ; k1 <= (F.Nflt1ch-1) ; k1++)
			{
				A = A + Afieldflattop1ch(t, F.flt1ch[k1].ti, F.flt1ch[k1].ton, F.flt1ch[k1].toff, F.flt1ch[k1].T, F.flt1ch[k1].o, F.flt1ch[k1].phi, F.flt1ch[k1].A, F.flt1ch[k1].b, F.flt1ch[k1].c);
			}
		}

		return A;	
	break;


	}


}


double AfieldEsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
{
	if ( (t <= ti) )
	{
		return 0.;
	} else if ( t >= (ti+pi/oc) ){
		return Primsin2cos(oc, o, phi, phi0, (ti+pi/oc) ) - Primsin2cos(oc, o, phi, phi0, ti);
	} else {
		return Primsin2cos(oc, o, phi, phi0, t) - Primsin2cos(oc, o, phi, phi0, ti);
	}

	// oc, o, A0, nc1, nc2, phi, phi0, ti;

	
}


double Primsin2cos(double a, double b, double c, double d, double t) // based on the antiderivative used in the SFA, not consistent notation... 
{
/*real*8 FUNCTION aDA(a,b,c,d,t);
	real*8, intent(in) :: a,b,c,d,t;
*/

// antiderivative of sin^2 in 't' (normalized by A0) antiderivative of sin^2(a*t+d)*cos(b*t+c)
	double aDA;

	aDA = (2. * cos(b*t)*sin(c))/b;
	aDA = aDA + (2.*cos(c)*sin(b*t))/b; 
	aDA = aDA + sin(c - 2.*d - 2.*a*t + b*t)/(2.*a - b);
	aDA = aDA - sin(c + 2.*d + 2.*a*t + b*t)/(2.*a + b);
	aDA = aDA/4.;

	return aDA;

// end function aDA; 
}


double Afieldflattop1(double t, double ti, double ton, double toff, double T, double o, double phi, double A0)
{
	double envelope;

	//oc2 = pi/(2.*toff); phienvel = pi - oc2*(T+ton+toff); tend = T+ton+toff;
	//return A0*pow(sin(oc2*t+phienvel),2.)*sin(o*t+phi);
	
	
	// envelope
	envelope = smootherstep(ti,ti+ton,t)*smootherstep(0.,toff,ti+ton+T+toff-t);

	return A0*envelope*sin(o*t+phi);	
	
}


double Afieldflattop1ch(double t, double ti, double ton, double toff, double T, double o, double phi, double A0, double b , double c)
{
	double envelope;

	//oc2 = pi/(2.*toff); phienvel = pi - oc2*(T+ton+toff); tend = T+ton+toff;
	//return A0*pow(sin(oc2*t+phienvel),2.)*sin(o*t+phi);
	
	
	// envelope
	envelope = smootherstep(ti,ti+ton,t)*smootherstep(0.,toff,ti+ton+T+toff-t);

	return (A0 + c*t) * envelope*sin(o*t+b*t*t+phi);	
	
}


double Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
{
	if ( (t <= ti) || ( t >= (ti+pi/oc) ) )
	{
		return 0.;
	} else {
		return A0*( pow(sin(oc*t + phi0),2.) ) * cos(o*t + phi);
	}

	// oc, o, A0, nc1, nc2, phi, phi0, ti;

	
}

double dAField(struct Efield_var F, double t) // ANAlytic field is -dA/dt, the sign!
{
	
	double omegap = F.trap.omega/((double)F.trap.nc),ts,Tf,A0;

	double a,b;
	double A,A1,A2,dum;
	int k1;

switch (F.fieldtype){

	case 2:  // analytic 
		A = 0.;
		//for(k1 = 0 ; k1 <= F.Nflt1 ; k1++)
		if (F.Nflt1 > 0)
		{
			perror("derivative of vectpot not implemented");
/*			for(k1 = 0 ; k1 <= (F.Nflt1-1) ; k1++)*/
/*			{*/
/*				A = A + Afieldflattop1(t, F.flt1[k1].ti, F.flt1[k1].ton, F.flt1[k1].toff, F.flt1[k1].T, F.flt1[k1].o, F.flt1[k1].phi, F.flt1[k1].A);*/
/*			}*/
		}

		if (F.Nsin2 > 0)
		{
			for(k1 = 0 ; k1 <= (F.Nsin2-1) ; k1++)
			{
				// Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
				A = A + dAfieldsin2(t , F.sin2[k1].ti , F.sin2[k1].A0 , F.sin2[k1].oc , F.sin2[k1].phi0 , F.sin2[k1].o , F.sin2[k1].phi);
			}
		}

		if (F.NEsin2 > 0)
		{
			perror("derivative of vectpot not implemented");
/*			for(k1 = 0 ; k1 <= (F.NEsin2-1) ; k1++)*/
/*			{*/
/*				// Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)*/
/*				A = A + AfieldEsin2(t , F.Esin2[k1].ti , F.Esin2[k1].A0 , F.Esin2[k1].oc , F.Esin2[k1].phi0 , F.Esin2[k1].o , F.Esin2[k1].phi);*/
/*			}*/
		}

		if (F.Nflt1ch > 0)
		{
			perror("derivative of vectpot not implemented");
/*			for(k1 = 0 ; k1 <= (F.Nflt1ch-1) ; k1++)*/
/*			{*/
/*				A = A + Afieldflattop1ch(t, F.flt1ch[k1].ti, F.flt1ch[k1].ton, F.flt1ch[k1].toff, F.flt1ch[k1].T, F.flt1ch[k1].o, F.flt1ch[k1].phi, F.flt1ch[k1].A, F.flt1ch[k1].b, F.flt1ch[k1].c);*/
/*			}*/
		}

		return A;	
	break;


	}


}


double dAfieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
{
	if ( (t <= ti) || ( t >= (ti+pi/oc) ) )
	{
		return 0.;
	} else {
		return A0*( oc * sin(2.*(oc*t + phi0)) * cos(o*t + phi) - o * pow(sin(oc*t + phi0),2.) * sin(o*t + phi));
	}

	// oc, o, A0, nc1, nc2, phi, phi0, ti;

	
}


// pi/2.0d0 - wc2*(delay+0.5d0*Tc2);

// attempt to do exact envelope in the electric field
/*
double Afieldflattop3(double t, double ton, double toff, double T, double o, double phi, double A0)
{
	 double oc1, oc2, phienvel, A, dum1, dum2, tend;


	if( (t >= ton ) && ( t <= (ton+T) ) )
	{
		oc1 = pi/(2.*ton);
		dum2 = (4*o*oc1*cos(phi + o*ton)*sin(2*oc1*ton) + 2*(pow(o,2) - 4*pow(oc1,2) - pow(o,2)*cos(2*oc1*ton))*sin(phi + o*ton))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // A(ton)
		dum2 = o*dum2;
		return A0*(sin(o*t+phi)-sin(o*ton+phi)+dum2);
	} else if (t < ton) {
		oc1 = pi/(2.*ton);

		A = (2*(-pow(o,2) + 4*pow(oc1,2) + pow(o,2)*cos(2*oc1*t))*cos(phi + o*t) + 4*o*oc1*sin(2*oc1*t)*sin(phi + o*t))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // undefinite integral result
		dum1 = (2*pow(oc1,2)*cos(phi))/(pow(o,3) - 4*o*pow(oc1,2)); // A(0)
		// dum2 = (4*o*oc1*cos(phi + o*ton)*sin(2*oc1*ton) + 2*(pow(o,2) - 4*pow(oc1,2) - pow(o,2)*cos(2*oc1*ton))*sin(phi + o*ton))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // A(ton)
		return o*A0*(A-dum1); // shift to be 0 in 0
	} else {
		oc2 = pi/(2.*toff); phienvel = (pi/2) - oc2*(T+ton); tend = T+ton+toff;

		A = ((-2*cos(phi + o*t))/o + cos(phi - 2*phienvel + o*t - 2*oc2*t)/(o - 2*oc2) + cos(phi + 2*phienvel + o*t + 2*oc2*t)/(o + 2*oc2))/4.; // undefinite integral result
		dum1 = ((-2*cos(phi + o*tend))/o + cos(phi - 2*phienvel + o*tend - 2*oc2*tend)/(o - 2*oc2) + cos(phi + 2*phienvel + o*tend + 2*oc2*tend)/(o + 2*oc2))/4.; // A(tend)
		return A0*(A-dum1);
	}
}
*/


// from wiki
double smootherstep(double edge0, double edge1, double x) {
  // Scale, and clamp x to 0..1 range
  x = clamp((x - edge0) / (edge1 - edge0), 0.0, 1.0);
  // Evaluate polynomial
  return x * x * x * (x * (x * 6 - 15) + 10);
}

double clamp(double x, double lowerlimit, double upperlimit) {
  if (x < lowerlimit)
  {
    x = lowerlimit;
  } else if (x > upperlimit) {
    x = upperlimit;
  }
  return x;
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





// MANIPULATION WITH DATA

void findinterval(int n, double x, double *x_grid, int *k1, int *k2) //! returns interval where is placed x value, if it is out of the range, 0 is used
{
//intervals are ordered: <..>(..>(..>...(..>
	int i;
	

	// printf("x_grid[0],  %lf \n",x_grid[0]);
	
	if( x < x_grid[0] )
	{
			*k1 = -1;
			*k2= 0;
			return;
	}

	for(i=0;i< n;i++)
	{
		if ( x <= x_grid[i+1] )
		{
			*k1 = i;
			*k2= i+1;
			// printf("interval,  %i \n",*k1);
			// printf("interval,  %i \n",*k2);
			return;
		}
	}
	
 	*k1=n; *k2=n+1;
	// !write(*,*) "error in the interval subroutine"

}





// NUMERICS

double interpolate( int n, double x, double *x_grid, double* y_grid) //!inputs: # of points, x(n), y(x(n)), x, returns y(x) (linearinterpolation), extrapolation by the boundary values
{
	int k1,k2;
	double y;
	
	k1=0;
	k2=0;
	findinterval(n, x, x_grid, &k1, &k2);
	// printf("\ninside interpolate \n");
	// printf("interval,  %i \n",k1);
	// printf("interval,  %i \n",k2);
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


double findnextinterpolatedzero(int n, double x, double* x_grid, double* y_grid) // it next zero according to an input value
{
	int k1,k2,k3,k4;
	double x_root,x1,x2,y1,y2;
	
	k1=0;
	k2=0;
	k3=0;
	findinterval(n, x, x_grid, &k1, &k2);
	// printf("\ninside interpolate \n");
	// printf("interval,  %i \n",k1);
	// printf("interval,  %i \n",k2);
	if( ( k1 == -1 ) ||  ( k2 == n+1 ))
	{
		printf("Cannot find interpolated zero: out of range\n");
	} else {
		k4 = k1;
		while (k4 < n)
		{
		if ( ( y_grid[k4]*y_grid[k4+1] ) < 0  )
		{
			x1 = x_grid[k4];
			x2 = x_grid[k4+1];
			y1 = y_grid[k4];
			y2 = y_grid[k4+1];
			x_root = (y2*x1-y1*x2)/(y2-y1);
			return x_root;
		}
		k4 = k4+1;
		}
		printf("There is no zero \n");	
	}
	
	//y = 0.;
	//printf("interpolated value,  %lf \n",y);
	//printf("\n");	
}



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

	dxi = 2.*pi/(  ((double)N) * dx);

	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);	

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

/*	dxi = 2.*pi/(  ((double)N) * dx); */
	dxi = 2.*pi/xmax; 

	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);	

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


	dxi = 2.*pi/(  ((double)N) * dx); // dimension-preserving coefficients
	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);

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

/*	dxi = 2.*pi/(  ((double)N) * dx); */
	dxi = 2.*pi/xmax; 

	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);

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

	dxi = 2.*pi/(  ((double)NF) * dx); 

	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);	

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


	dxi = 2.*pi/(  ((double)N) * dx); // dimension-preserving coefficients
	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);

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

	out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * Nc);
	for(k1 = 0; k1 <= (N-1); k1++){in[k1]=signal2[k1];} // !!! REDUNDANT
	p = fftw_plan_dft_r2c_1d(N, in, out2, FFTW_ESTIMATE); //fftw_plan_dft_r2c_1d(int n, double *in, fftw_complex *out, unsigned flags); // plan FFTW
	fftw_execute(p);


	//	RESCALE TO OUTPUTS

	dxi = 2.*pi/xmax;
	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);


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
	dxi = 2.*pi/xmax;
	coeff1 = dx/ sqrt(2.*pi); coeff2 = dx*dx/(2.*pi);


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

	fftw_free(out); free(in);
	*Nxi = Nc;
	return;

}


double ** create_2Darray_accessor_real(int * dims, double *array_data) //takes a contiguous block of memory and reconstruct a 2D arroy from that !!! generalise it to nD using void* (see discussion)
{
	int k1;
	double **array_accessor;
	array_accessor = (double**) malloc(dims[0]*sizeof(double));
	for(k1 = 0; k1 < dims[0];k1++){array_accessor[k1] = &array_data[dims[1]*k1];}	
	return array_accessor;
}
// how-to generalise: there could be a problem to declare a correct number of '*', chain somehow voids (seems to be possible)? or hot-fix it by log if?


// void * create_nDarray_accessor(int ndims, int * dims, void *array_data) //takes a contiguous block of memory and reconstruct a 2D array from that !!! generalise it to nD using void* (see discussion)
// { // it's quite intersting exercise, but the gain is small anyway, write rather just index-mapping function
// 	int k1;
// 	size_t size = sizeof(&array_data[0]); 
// 	// we find respective sizes of pointers
// 	size_t * sizes;
// 	void * ptr;
// 	sizes[0] = sizeof(&array_data[0]); ptr = malloc(sizes[0]);
// 	for(k1=1; k1<ndims; k1++){sizes[k1]=sizeof(&ptr[k1-1]); free(ptr); ptr = malloc(sizes[k1]);}
// 	void * accesor;
// 	size_t accessor_size;
// 	accesor_size = dims[0]*sizes[0];
// 	for(k1=1; k1 < ndims; k1++){accesor_size+=dims[k1]*sizes[k1];} // product shoul be involved, not only sum
// 	accesor = malloc(accesor_size); // much bigger draw it by multiplications, aftermost pointers should point to every line (i.e. N1*N2*...*N(n-1))-pointers
// 	// fill pointers here // most of them points within the structure and only last ones outside to the array
	
// 	// double **array_accessor;
// 	// array_accessor = (double**) malloc(dims[0]*sizeof(double));
// 	// for(k1 = 0; k1 < dims[0];k1++){array_accessor[k1] = &array_data[dims[1]*k1];}	
// 	// return array_accessor;
// 	return accessor;
// }