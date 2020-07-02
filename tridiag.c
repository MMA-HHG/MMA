#include<math.h>
#include<float.h>                // required for DBL_EPSILON
#include<malloc.h>
#include<stdlib.h>
#include<stdio.h>

#include "hdf5.h"
#include "mpi.h"
#include"util.h"





//                        Internally Defined Routines 

////////////////////////////////////////////////////////////////////////////////
//  int QL_Tridiagonal_Symmetric_Matrix(double diagonal[],                    //
//          double off_diagonal[], double *U, int n, int max_iteration_count) //
//                                                                            //
//  Description:                                                              //
//     This function implements the QL algorithm with implicit shifts of the  //
//     of the origin for a symmetric tridiagonal matrix.  If the original     //
//     symmetric matrix was tridiagonalized using Householder transformations //
//     then the tridiagonal matrix should be organized so that the largest    //
//     elements are at the bottom right-hand of the matrix.  Also, the process//
//     of tridiagonalizing the matrix will introduce round-off errors, so that//
//     the estimates of the eigenvalues and eigenvectors may be improved by   //
//     calling the inverse_power_method.                                      //
//                                                                            //
//  Arguments:                                                                //
//     double diagonal[]                                                      //
//            On input, the diagonal of the tridiagonal symmetric matrix.     //
//            On output, the eigenvalues.                                     //
//     double off_diagonal[]                                                  //
//            The subdiagonal and superdiagonal of the tridiagonal symmetric  //
//            matrix.  The subdiagonal and superdiagonal are (off_diagonal[1],//
//            ,,, off_diagonal[n-1]).                                         //
//     double *U                                                              //
//            On input, U[n][n] is the identity matrix if the tridiagonal     //
//            matrix is the primary data or is the transformation matrix      //
//            of Householder transformations which is the output of           //
//            tridiagonalize_symmetric_matrix.  On output, U[n][n] is the     //
//            matrix of eigenvectors, the ith column being the eigenvector    //
//            corresponding to the eigenvalue d[i].  These are the eigen-     //
//            vectors of the full matrix, if on input U is not the identity.  //
//     int    n                                                               //
//            The number of rows or columns of the symmetric matrix A and of  //
//            the matrix U.                                                   //
//     int    max_iteration_count                                             //
//            The maximum number of iterations to try to annihilate an off-   //
//            diagonal elements.                                              //
//                                                                            //
//  Return Values:                                                            //
//     0  Success                                                             //
//    -1  Failure - Failed to converge within max_iteration_count iterations. //
//                                                                            //
//  Example:                                                                  //
//     #define N                                                              //
//     double U[N][N], diagonal[N], double off_diagonal[N];                   //
//     int err;                                                               //
//     int max_iteration_count = 30;                                          //
//                                                                            //
//     (your code to create the vectors diagonal and off_diagonal)            //
//     if (QL_Tridiagonal_Symmetric_Matrix(diagonal, off_diagonal, (double*)U,//
//         N, max_iteration_count) ) printf("Failed to converge\n");          //
//     else printf("Success\n");                                              //
////////////////////////////////////////////////////////////////////////////////
//  
                                                                          //
int QL_Tridiagonal_Symmetric_Matrix( double *diagonal, double *off_diagonal,
                                   double *U, int n, int max_iteration_count)
{
   int i, j, k;
   int iteration;
   double epsilon; 
   double *p_off = off_diagonal + 1;
   double s,c,g,h,q;
   double shift;
   double dum;

   p_off[n-1] = 0.0;
   for (i = 0; i < n; i++) {
      for (iteration = 0; iteration < max_iteration_count; iteration++) { 
         for (k = i; k < (n - 1); k++) {
            epsilon = DBL_EPSILON * ( fabs(diagonal[k]) + fabs(diagonal[k+1]) );
            if ( fabs(p_off[k]) <= epsilon ) break;
         }
         if ( k == i ) break;
         shift = Calculate_Shift(diagonal[i+1], diagonal[i], p_off[i]);
         q = diagonal[k] - shift;
         c = 1.0;
         s = 1.0;
         shift = 0.0;
         for (j = k - 1; j >= i; j--) {
            h = c * p_off[j];
            g = s * p_off[j];
            if ( fabs( g ) >= fabs( q ) ) {
               c = q / g;
               dum = sqrt( c * c + 1.0 );
               p_off[j+1] = g * dum;
               s = 1.0 / dum;
               c *= s;
            }
            else {
               s = g / q;
               dum = sqrt( s * s + 1.0 );
               p_off[j+1] = q * dum;
               c = 1.0 / dum;
               s *= c;
            }
            q = diagonal[j + 1] - shift;
            dum = s * (diagonal[j] - q) + 2.0 * c * h;
            shift = s * dum;
            diagonal[j+1] = q + shift;
            q = c * dum - h;
            Transform_Matrix(U, s, c, j, n);
         }
         diagonal[i] -= shift;
         p_off[i] = q;
         p_off[k] = 0.0;
      }
   }
   if ( iteration >= max_iteration_count ) return -1;    
   return 0;
}

double Calculate_Shift(double d2, double d1, double off) 
{
   double h;

   h = ( d2 - d1 ) / ( off + off );
   off = sqrt( h * h + 1.0 );
   d2 = ( h < 0.0 ) ? h - off : h + off;
   return d1 - off / d2;
}

void Transform_Matrix(double* U, double s, double c, int col, int n)
{
   double x;
   int i;
   int col1 = col + 1;

   for (i = 0; i < n; i++, U += n) {
      x = *(U + col1);
      *(U + col1) = s * *(U + col) + c * x;
      *(U + col) *= c;
      *(U + col) -= s * x;
   }
}



//////////////////////////////////////////////////////////////////////////////////////////////
//
//				Calculation of the vector u such as M*u = r where M is tridiagonal NOT necessarly symatric.
//               M is build up of diagonal elements b[1...n]
//					the upper element is c[1...n-1]
//					the lower element is a[1...n-1]
//
//
//					ATTENNTION : THE FIRST ELEMENT b[0] should be != 0	
//////////////////////////////////////////////////////////////////////////////////////////////
void Inv_Tridiagonal_Matrix( double *a, double *b, double *c, double *r, double *u, int n)
{
	int j;
	double bet,*gam;

	gam = (double *)calloc(n+1,sizeof(double));

	if(b[0] == 0.) printf("Error 1 in Inv_Tridiagonal_Matrix");
	u[0]=r[0]/(bet=b[0]);

	printf("init ! %f \n",u[0]);
	for(j=1;j<=n-1;j++)
	{
		gam[j]=c[j-1]/bet;
		bet=b[j]-a[j-1]*gam[j];
		printf("beta ! %f \n",bet);
		if(bet == 0.) printf("Error 2 in Inv_Tridiagonal_Matrix");
		u[j] = (r[j]-a[j-1]*u[j-1])/bet;

		printf("boucle ! %f \n",u[j]);
	}

	for(j=(n-2);j>=0;j--) u[j] -= gam[j+1]*u[j+1];

	free(gam);

}

// For deltails see Inv_Tridiagonal_Matrix same dut for complex
void Inv_Tridiagonal_Matrix_complex( double *a, double *b, double *c, double *r, double *u, int n)
{
	int j;
	double bet_re,bet_im,bet_norm,*gam;

	gam = (double *)calloc(2*(n+1),sizeof(double));

	if((b[0] == 0.) && (b[1] == 0.)) printf("Error 1 in Inv_Tridiagonal_Matrix_complex");
	u[0]=(r[0]*b[0]+r[1]*b[1])/(b[0]*b[0]+b[1]*b[1]) /*real*/;u[1]=(r[1]*b[0]-r[0]*b[1])/(b[0]*b[0]+b[1]*b[1]);/*imag*/
	bet_re = b[0];bet_im = b[1]; bet_norm = bet_re*bet_re + bet_im*bet_im;

	for(j=1;j<=n-1;j++)
	{
		gam[2*j]=(c[2*(j-1)]*bet_re + c[2*(j-1)+1]*bet_im)/bet_norm;
		gam[2*j+1]=(c[2*(j-1)+1]*bet_re - c[2*(j-1)]*bet_im)/bet_norm;

		bet_re=b[2*j] - (a[2*(j-1)]*gam[2*j]-a[2*(j-1)+1]*gam[2*j+1]);
		bet_im=b[2*j+1] - (a[2*(j-1)+1]*gam[2*j]+a[2*(j-1)]*gam[2*j+1]);
		bet_norm = bet_re*bet_re + bet_im*bet_im;
		
		if(bet_norm == 0.) printf("Error 2 in Inv_Tridiagonal_Matrix_complex");

		u[2*j] = (r[2*j]-a[2*(j-1)]*u[2*(j-1)]+a[2*(j-1)+1]*u[2*(j-1)+1])*bet_re;
		u[2*j] += (r[2*j+1]-a[2*(j-1)]*u[2*(j-1)+1]-a[2*(j-1)+1]*u[2*(j-1)])*bet_im;
		u[2*j] = u[2*j]/bet_norm;
		u[2*j+1] = -(r[2*j]-a[2*(j-1)]*u[2*(j-1)]+a[2*(j-1)+1]*u[2*(j-1)+1])*bet_im;
		u[2*j+1] += (r[2*j+1]-a[2*(j-1)]*u[2*(j-1)+1]-a[2*(j-1)+1]*u[2*(j-1)])*bet_re;
		u[2*j+1] = u[2*j+1]/bet_norm;
	}



	for(j=(n-2);j>=0;j--) 
	{
		u[2*j] -= gam[2*(j+1)]*u[2*(j+1)]-gam[2*(j+1)+1]*u[2*(j+1)+1];
		u[2*j+1] -= gam[2*(j+1)]*u[2*(j+1)+1]+gam[2*(j+1)+1]*u[2*(j+1)];
	}
	
	
	free(gam);

}


double Einitialise(struct trg_def trg, double *psi0,double *dinf,double *d,double *dsup,double *x,double Eguess,double CV,int num_r)
{
    	double *res,*dnew,*diag,*dinfnew,*dsupnew;
	double Energy,test,Eold;
	int i,size = 2*(num_r+1);
	double sum;

	res = (double *)calloc(size,sizeof(double));
	dnew = (double *)calloc(size,sizeof(double));
	dinfnew = (double *)calloc(size,sizeof(double));
	dsupnew = (double *)calloc(size,sizeof(double));
	diag = (double *)calloc(size,sizeof(double));


	printf("%e \n",psi0[0]);
printf("%e \n",psi0[0]);
printf("%e \n",dinf[0]);
printf("%e \n",d[0]);
printf("%e \n",dsup[0]);
printf("%e \n",x[0]);
printf("%e \n",CV);


	printf("%e \n",psi0[3]);
printf("%e \n",psi0[3]);
printf("%e \n",dinf[3]);
printf("%e \n",d[3]);
printf("%e \n",dsup[3]);
printf("%e \n",x[3]);
printf("%e \n",CV);

	
	  for(i=0;i<=num_r;i++)
	  {
		  dinfnew[2*i] = dinf[2*i] - Eguess/12. + potential(x[i],trg)/12.; dinfnew[2*i+1] = dinf[2*i+1];		  
		  dnew[2*i] = 10*potential(x[i],trg)/12.+ d[2*i] - 10*Eguess/12.; dnew[2*i+1] = d[2*i+1];
		  dsupnew[2*i] = dsup[2*i] - Eguess/12. + potential(x[i+1],trg)/12.; dsupnew[2*i+1] = dsup[2*i+1];
		  diag[2*i] = potential(x[i],trg)+ d[2*i] ; diag[2*i+1] = d[2*i+1];
	  }




	Eold = Eguess;
	printf(" Energy : %e\n",Eold);
	do
	{
	
	  Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew,dnew,dsupnew,psi0,res,num_r); // WARNING here num_r
	  //res = Inv_Tridiagonal_Matrix_complex(dinfnew,dnew,dsupnew,psi0,res,num_r+1);

	  normalise(res,num_r); for(i=0;i<=num_r;i++){psi0[2*i] = res[2*i]; psi0[2*i+1] = res[2*i+1];}

	  //Energy = E_calculation(res,dinf,diag,dsup,num_r);
	  Energy = E_calculation_numerov(trg, res,d[0],x,num_r);

	  test = sqrt((Energy-Eold)*(Energy-Eold));
	  Eold = Energy;

	  printf("cv : %e",test);
	  printf(" Energy : %e\n",Eold);
	
	}
	while(test > CV);
	
	normalise(psi0,num_r);

	//sum = 0;
	//for(i=0;i<num_r;i++){sum = sum + psi0[2*i]*psi0[2*i]*potential(x[i],trg);}
	//printf("Average value of the potential : %e\n",sum);
	
	
	return Energy;

	free(dnew); free(res); free(diag);free(dinfnew);free(dsupnew);
}



double E_calculation(double *psi0,double *dinf,double *d,double *dsup,int num_r)
{
    double sum_re = 0,sum_im = 0,sum_re_1,sum_re_2,sum_re_3,sum_im_1,sum_im_2,sum_im_3;
    int i;

	
	sum_re_2 = d[0]*psi0[0] - d[1]*psi0[1];
	sum_im_2 = d[1]*psi0[0] + d[0]*psi0[1];
	
	sum_re_2 = sum_re_2*psi0[0] + sum_im_2*psi0[1]; // WARNING Here it is the complex conjugate
	sum_im_2 = -sum_re_2*psi0[1] + sum_im_2*psi0[0]; 

	sum_re_3 = dsup[0]*psi0[2] - dsup[1]*psi0[3];
	sum_im_3 = dsup[1]*psi0[2] + dsup[0]*psi0[3];

	sum_re_3 = sum_re_3*psi0[0] + sum_im_3*psi0[1]; 
	sum_im_3 = -sum_re_3*psi0[1] + sum_im_3*psi0[0]; 

	sum_re = sum_re_2 + sum_re_3; sum_im = sum_im_2 + sum_im_3;
	
	for (i=1;i<=num_r-1;i++)
	{
      sum_re_1 = dinf[2*(i-1)]*psi0[2*(i-1)] - dinf[2*(i-1)+1]*psi0[2*(i-1)+1];
	  sum_im_1 = dinf[2*(i-1)+1]*psi0[2*(i-1)] + dinf[2*(i-1)]*psi0[2*(i-1)+1];

	  sum_re_1 = sum_re_1*psi0[2*i] + sum_im_1*psi0[2*i+1]; 
	  sum_im_1 = -sum_re_1*psi0[2*i+1] + sum_im_1*psi0[2*i]; 
	 
	  sum_re_2 = d[2*i]*psi0[2*i] - d[2*i+1]*psi0[2*i+1];
	  sum_im_2 = d[2*i+1]*psi0[2*i] + d[2*i]*psi0[2*i+1];

	  sum_re_2 = sum_re_2*psi0[2*i] + sum_im_2*psi0[2*i+1]; 
	  sum_im_2 = -sum_re_2*psi0[2*i+1] + sum_im_2*psi0[2*i]; 

	  sum_re_3 = dsup[2*(i+1)]*psi0[2*(i+1)] - dsup[2*(i+1)+1]*psi0[2*(i+1)+1];
	  sum_im_3 = dsup[2*(i+1)+1]*psi0[2*(i+1)] + dsup[2*(i+1)]*psi0[2*(i+1)+1];

	  sum_re_3 = sum_re_3*psi0[2*i] + sum_im_3*psi0[2*i+1]; 
	  sum_im_3 = -sum_re_3*psi0[2*i+1] + sum_im_3*psi0[2*i]; 
	
	  sum_re = sum_re + sum_re_1 + sum_re_2 + sum_re_3; sum_im = sum_im + sum_im_1 + sum_im_2 + sum_im_3;

	}

	  sum_re_2 = d[2*(num_r)]*psi0[2*(num_r)] - d[2*(num_r)+1]*psi0[2*(num_r)+1];
	  sum_im_2 = d[2*(num_r)+1]*psi0[2*(num_r)] + d[2*(num_r)]*psi0[2*(num_r)+1];

	  sum_re_2 = sum_re_2*psi0[2*(num_r)] + sum_im_2*psi0[2*(num_r)+1]; 
	  sum_im_2 = -sum_re_2*psi0[2*(num_r)+1] + sum_im_2*psi0[2*(num_r)]; 

	  sum_re_3 = dinf[2*(num_r-1)]*psi0[2*(num_r-1)] - dinf[2*(num_r-1)+1]*psi0[2*(num_r-1)+1];
	  sum_im_3 = dinf[2*(num_r-1)+1]*psi0[2*(num_r-1)] + dinf[2*(num_r-1)]*psi0[2*(num_r-1)+1];

	  sum_re_3 = sum_re_3*psi0[2*(num_r)] + sum_im_3*psi0[2*(num_r)+1]; 
	  sum_im_3 = -sum_re_3*psi0[2*(num_r)+1] + sum_im_3*psi0[2*(num_r)]; 
	
	  sum_re = sum_re + sum_re_2 + sum_re_3; sum_im = sum_im + sum_im_2 + sum_im_3;

	  return sum_re;  // for an hermitian hamiltonian, sum_im is 0


}


double E_calculation_numerov(struct trg_def trg, double *psi,double dx,double *x,int num_r)
{
    int j;
	double *psi_inter1,*psi_inter2;
	double coef,E,E1,E2,norme;
	
	psi_inter1 = calloc(2*(num_r+1),sizeof(double));
	psi_inter2 = calloc(2*(num_r+1),sizeof(double));




	coef = dx; // that is in fact 1/dx^2

		psi_inter1[0] = coef*psi[0]-0.5*coef*psi[2];
		psi_inter1[0] = psi_inter1[0]+((10/12.)*psi[0]*(potential(x[0],trg))
						+(1/12.)*psi[2]*(potential(x[1],trg)));
		psi_inter1[1] = coef*psi[1]-0.5*coef*psi[3];
		psi_inter1[1] = psi_inter1[1]+((10/12.)*psi[1]*(potential(x[0],trg))
						+(1/12.)*psi[3]*(potential(x[1],trg)));
		
		psi_inter2[0] = 10*psi[0]/12. + psi[2]/12.;
		psi_inter2[1] = 10*psi[1]/12. + psi[3]/12.;

		for(j = 1 ; j< num_r ; j++)
		{

			psi_inter1[2*j] = coef*psi[2*j]-0.5*coef*(psi[2*(j-1)]+psi[2*(j+1)]);
			psi_inter1[2*j] = psi_inter1[2*j]+((10/12.)*psi[2*j]*(potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)]*(potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)]*(potential(x[j+1],trg)));

			psi_inter1[2*j+1] = coef*psi[2*j+1]-0.5*coef*(psi[2*(j-1)+1]+psi[2*(j+1)+1]);
			psi_inter1[2*j+1] = psi_inter1[2*j+1]+((10/12.)*psi[2*j+1]*(potential(x[j],trg))
							  +(1/12.)*psi[2*(j-1)+1]*(potential(x[j-1],trg))
							  +(1/12.)*psi[2*(j+1)+1]*(potential(x[j+1],trg)));

			psi_inter2[2*j] = 10*psi[2*j]/12. + (psi[2*(j+1)]+psi[2*(j-1)])/12.;
			psi_inter2[2*j+1] = 10*psi[2*j+1]/12. + (psi[2*(j+1)+1]+psi[2*(j-1)+1])/12.;

		}

		psi_inter1[2*num_r] = coef*psi[2*num_r]-0.5*coef*psi[2*(num_r-1)];
		psi_inter1[2*num_r] = psi_inter1[2*num_r]+((10/12.)*psi[2*num_r]*(potential(x[num_r],trg))
							  +(1/12.)*psi[2*(num_r-1)]*(potential(x[num_r-1],trg)));
		psi_inter1[2*num_r+1] = coef*psi[2*num_r+1]-0.5*coef*psi[2*(num_r-1)+1];
		psi_inter1[2*num_r+1] = psi_inter1[2*num_r+1]+((10/12.)*psi[2*num_r+1]*(potential(x[num_r],trg))
							  +(1/12.)*psi[2*(num_r-1)+1]*(potential(x[num_r-1],trg)));
			
		psi_inter2[2*num_r] = 10*psi[2*num_r]/12. + psi[2*(num_r-1)]/12.;
		psi_inter2[2*num_r+1] = 10*psi[2*num_r+1]/12. + psi[2*(num_r-1)+1]/12.;


		E = 0;norme = 0;
		for(j = 0 ; j<= num_r ; j++)
		{
			E1 = psi_inter1[2*j]*psi_inter2[2*j]+psi_inter1[2*j+1]*psi_inter2[2*j+1];
			E2 = -psi_inter1[2*j]*psi_inter2[2*j+1]+psi_inter1[2*j+1]*psi_inter2[2*j];
			E+=E1+E2;
			
			norme += psi_inter2[2*j]*psi_inter2[2*j]+psi_inter2[2*j+1]*psi_inter2[2*j+1];
		}


		free(psi_inter1);free(psi_inter2);

		return E/norme;

}


void Inv_Tridiagonal_Matrix_complex_Numerov( double *dinf, double *d, double *dsup, double *psi, double *res,int num_r)
{

	// This function calculate (H)^-1*M2*Psi  where H is defined by dinf,d,dsup

        double *psi_inter;
	int i,size = 2*(num_r+1);


	psi_inter = calloc(size,sizeof(double));

	
	  psi_inter[0] = 10*psi[0]/12. + psi[2]/12.; psi_inter[1] = 10*psi[1]/12. + psi[3]/12.;
	  for(i=1;i<num_r;i++)
	  {
		  psi_inter[2*i] = 10*psi[2*i]/12. + (psi[2*(i+1)]+psi[2*(i-1)])/12.; 
		  psi_inter[2*i+1] = 10*psi[2*i+1]/12. + (psi[2*(i+1)+1]+psi[2*(i-1)+1])/12.;		  
	  }
	  psi_inter[2*num_r] = 10*psi[2*num_r]/12. + psi[2*(num_r-1)]/12.; 
	  psi_inter[2*num_r+1] = 10*psi[2*num_r+1]/12. + psi[2*(num_r-1)+1]/12.;
	 

	  Inv_Tridiagonal_Matrix_complex(dinf,d,dsup,psi_inter,res,num_r+1);

	  free(psi_inter);

}
