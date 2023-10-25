/**
 * @file tridiag.c
 * 
 * @brief Contains methods for computing systems with tridiagonal matrices.
 * 
 * @details The complex variables are stored in a 1D double array of size 2N, 
 * i.e., Re(psi[0]) = psi[0], Im(psi[0]) = psi[1] etc.
 */
#include "tridiag.h"
#include "structures.h"
#include "tools.h"

/**
 * @brief Computes inverse of complex tridiagonal matrix
 * 
 * @details The procedure computes complex system of tridiagonal matrix defined
 * via its sub/super- diagonal.
 * 
 * @param a Subdiagonal.
 * @param b Diagonal.
 * @param c Superdiagonal.
 * @param r RHS vector.
 * @param u Solution vector.
 * @param n Rank of the matrix.
 */
void Inv_Tridiagonal_Matrix_complex(double *a, double *b, double *c, double *r, double *u, int n)
{
	int j;
	double bet_re,bet_im,bet_norm,*gam;

	gam = (double *)calloc(2*(n+1),sizeof(double));

	if((b[0] == 0.) && (b[1] == 0.)) printf("Error 1 in Inv_Tridiagonal_Matrix_complex");
	/*real*/
	u[0]=(r[0]*b[0]+r[1]*b[1])/(b[0]*b[0]+b[1]*b[1]);
	/*imag*/
	u[1]=(r[1]*b[0]-r[0]*b[1])/(b[0]*b[0]+b[1]*b[1]);
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

/**
 * @brief Computes energy of the ground state using the resolvent method.
 * 
 * @param trg Target info.
 * @param psi0 Initial wavefunction.
 * @param dinf Subdiagonal of hamiltonian matrix.
 * @param d Diagonal of hamiltonian matrix.
 * @param dsup Superdiagonal of hamiltonian matrix.
 * @param x Spatial grid.
 * @param Eguess Initial guess of the ground state energy.
 * @param CV Error threshold for energy computation.
 * @param num_r Number of points in the wavefunction.
 * @return double 
 */
double Einitialise(trg_def trg, double *psi0, double *dinf, double *d, double *dsup,
	double *x, double Eguess, double CV, int num_r)
{
    double *res,*dnew,*diag,*dinfnew,*dsupnew;
	double Energy,test,Eold,dx;
	int i,size = 2*(num_r+1);

	res = (double *)calloc(size,sizeof(double));
	dnew = (double *)calloc(size,sizeof(double));
	dinfnew = (double *)calloc(size,sizeof(double));
	dsupnew = (double *)calloc(size,sizeof(double));
	diag = (double *)calloc(size,sizeof(double));


	// Initialization of tridiagonal matrix of the Hamiltonian with Numerov
	for(i = 0; i <= num_r; i++)
	{
		dinfnew[2*i] = dinf[2*i] - Eguess/12. + potential(x[i],trg)/12.; 
		dinfnew[2*i+1] = dinf[2*i+1];		  
		dnew[2*i] = 10*potential(x[i],trg)/12.+ d[2*i] - 10*Eguess/12.; 
		dnew[2*i+1] = d[2*i+1];
		diag[2*i] = potential(x[i],trg)+ d[2*i]; 
		diag[2*i+1] = d[2*i+1];
	}
	for(i = 0; i < num_r; i++) {
		dsupnew[2*i] = dsup[2*i] - Eguess/12. + potential(x[i+1],trg)/12.; 
		dsupnew[2*i+1] = dsup[2*i+1];
	}
	dx = x[num_r] - x[num_r-1];
	dsupnew[2*num_r] = dsup[2*num_r] - Eguess/12. + potential(x[num_r]+dx,trg)/12.; dsupnew[2*num_r+1] = dsup[2*num_r+1];


	Eold = Eguess;
	// Application of resolvent until the most precise guess of the ground state is found
	do {
		// Inverse Hamiltonian with Numerov
		Inv_Tridiagonal_Matrix_complex_Numerov(dinfnew, dnew, dsupnew, psi0, res, num_r); 
		// Normalise wavefunction
		normalise(res, num_r); 
		// Set new initial state
		for(i=0; i <= num_r; i++) {
			psi0[2*i] = res[2*i]; 
			psi0[2*i+1] = res[2*i+1];
		}
		// Compute average value of Hamiltonian - Energy
		Energy = E_calculation_numerov(trg, res, d[0], x, num_r);

		// Compute stopping criterium
		test = sqrt((Energy-Eold)*(Energy-Eold));
		Eold = Energy;

	}
	while(test > CV);
	
	// Normalise the ground state
	normalise(psi0, num_r);

	// Free memory
	free(dnew); 
	free(res); 
	free(diag);
	free(dinfnew);
	free(dsupnew);
	return Energy;	
}

/**
 * @brief Computes dot product of the tridiagonal matrix solution with the diffused wavefunction 
 * and returns the energy E.
 * 
 * @param trg Target information.
 * @param psi Solution of tridiagonal matrix.
 * @param dx Spatial step.
 * @param x Spatial grid.
 * @param num_r Number of gridpoints.
 * @return double Energy
 */
double E_calculation_numerov(trg_def trg, double *psi, double dx, double *x, int num_r)
{
    int j;
	double *psi_inter1, *psi_inter2;
	double coef, E, E1, E2, norme;
	
	psi_inter1 = calloc(2*(num_r+1),sizeof(double));
	psi_inter2 = calloc(2*(num_r+1),sizeof(double));

	coef = dx; // that is in fact 1/dx^2

	// Compute dot product
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


	E = 0;
	norme = 0;
	// Compute norm and energy
	for(j = 0 ; j<= num_r ; j++)
	{
		E1 = psi_inter1[2*j]*psi_inter2[2*j]+psi_inter1[2*j+1]*psi_inter2[2*j+1];
		E2 = -psi_inter1[2*j]*psi_inter2[2*j+1]+psi_inter1[2*j+1]*psi_inter2[2*j];
		E+=E1+E2;
		
		norme += psi_inter2[2*j]*psi_inter2[2*j]+psi_inter2[2*j+1]*psi_inter2[2*j+1];
	}

	free(psi_inter1);
	free(psi_inter2);

	return E/norme;
}

/**
 * @brief Computes (H)^-1*M2*Psi where H (Hamiltonian) is defined by subdiagonal, 
 * diagonal and superdiagonal.
 * 
 * @param dinf Subdiagonal.
 * @param d Diagonal.
 * @param dsup Superdiagonal.
 * @param psi Wavefunction.
 * @param res Result of the inverse.
 * @param num_r Matrix rank.
 */
void Inv_Tridiagonal_Matrix_complex_Numerov(double *dinf, double *d, double *dsup, 
	double *psi, double *res, int num_r)
{
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
	
	// Compute inverse matrix
	Inv_Tridiagonal_Matrix_complex(dinf,d,dsup,psi_inter,res,num_r+1);

	free(psi_inter);
}
