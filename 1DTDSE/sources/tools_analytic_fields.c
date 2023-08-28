#include "constants.h"
#include "structures.h"
#include "tools_analytic_fields.h"

clock_t start, finish;
clock_t start2, finish2;

extern double* timet,dipole;

// extern struct Efield_var;

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
	} else if ( t >= (ti+Pi/oc) ){
		return Primsin2cos(oc, o, phi, phi0, (ti+Pi/oc) ) - Primsin2cos(oc, o, phi, phi0, ti);
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

	//oc2 = Pi/(2.*toff); phienvel = Pi - oc2*(T+ton+toff); tend = T+ton+toff;
	//return A0*pow(sin(oc2*t+phienvel),2.)*sin(o*t+phi);
	
	
	// envelope
	envelope = smootherstep(ti,ti+ton,t)*smootherstep(0.,toff,ti+ton+T+toff-t);

	return A0*envelope*sin(o*t+phi);	
	
}


double Afieldflattop1ch(double t, double ti, double ton, double toff, double T, double o, double phi, double A0, double b , double c)
{
	double envelope;

	//oc2 = Pi/(2.*toff); phienvel = Pi - oc2*(T+ton+toff); tend = T+ton+toff;
	//return A0*pow(sin(oc2*t+phienvel),2.)*sin(o*t+phi);
	
	
	// envelope
	envelope = smootherstep(ti,ti+ton,t)*smootherstep(0.,toff,ti+ton+T+toff-t);

	return (A0 + c*t) * envelope*sin(o*t+b*t*t+phi);	
	
}


double Afieldsin2(double t, double ti, double A0, double oc, double phi0, double o, double phi)
{
	if ( (t <= ti) || ( t >= (ti+Pi/oc) ) )
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
	if ( (t <= ti) || ( t >= (ti+Pi/oc) ) )
	{
		return 0.;
	} else {
		return A0*( oc * sin(2.*(oc*t + phi0)) * cos(o*t + phi) - o * pow(sin(oc*t + phi0),2.) * sin(o*t + phi));
	}

	// oc, o, A0, nc1, nc2, phi, phi0, ti;

	
}


// Pi/2.0d0 - wc2*(delay+0.5d0*Tc2);

// attempt to do exact envelope in the electric field
/*
double Afieldflattop3(double t, double ton, double toff, double T, double o, double phi, double A0)
{
	 double oc1, oc2, phienvel, A, dum1, dum2, tend;


	if( (t >= ton ) && ( t <= (ton+T) ) )
	{
		oc1 = Pi/(2.*ton);
		dum2 = (4*o*oc1*cos(phi + o*ton)*sin(2*oc1*ton) + 2*(pow(o,2) - 4*pow(oc1,2) - pow(o,2)*cos(2*oc1*ton))*sin(phi + o*ton))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // A(ton)
		dum2 = o*dum2;
		return A0*(sin(o*t+phi)-sin(o*ton+phi)+dum2);
	} else if (t < ton) {
		oc1 = Pi/(2.*ton);

		A = (2*(-pow(o,2) + 4*pow(oc1,2) + pow(o,2)*cos(2*oc1*t))*cos(phi + o*t) + 4*o*oc1*sin(2*oc1*t)*sin(phi + o*t))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // undefinite integral result
		dum1 = (2*pow(oc1,2)*cos(phi))/(pow(o,3) - 4*o*pow(oc1,2)); // A(0)
		// dum2 = (4*o*oc1*cos(phi + o*ton)*sin(2*oc1*ton) + 2*(pow(o,2) - 4*pow(oc1,2) - pow(o,2)*cos(2*oc1*ton))*sin(phi + o*ton))/(4.*(pow(o,3) - 4*o*pow(oc1,2))); // A(ton)
		return o*A0*(A-dum1); // shift to be 0 in 0
	} else {
		oc2 = Pi/(2.*toff); phienvel = (Pi/2) - oc2*(T+ton); tend = T+ton+toff;

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