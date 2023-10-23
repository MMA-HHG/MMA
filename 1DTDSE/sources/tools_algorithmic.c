/**
 * @file tools_algorithmic.c
 * @brief Contains interpolation tools.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "tools_algorithmic.h"

/**
 * @brief Finds value in an interval within indices k1 and k2.
 * 
 * @param n Grid size.
 * @param x Searched value.
 * @param x_grid Array to be searched.
 * @param k1 Index before x.
 * @param k2 Index after x.
 */
void findinterval(int n, double x, double *x_grid, int *k1, int *k2) //! returns interval where is placed x value, if it is out of the range, 0 is used
{
	// intervals are ordered: <..>(..>(..>...(..>
	int i;
	
	if(x < x_grid[0])
	{
		*k1 = -1;
		*k2= 0;
		return;
	}

	for(i = 0; i < n; i++)
	{
		if (x <= x_grid[i+1])
		{
			*k1 = i;
			*k2= i+1;
			return;
		}
	}
	
 	*k1=n; 
	*k2=n+1;
}

/**
 * @brief Coarsens a real array.
 * 
 * @param in_array Input array to be coarsened.
 * @param out_array Coarsened output array.
 * @param length_out Output array length after coarsening.
 * @param k_step Number of steps to skip.
 * @param N_max Maximum number of points.
 */
void coarsen_grid_real(double *in_array, double **out_array, int *length_out, int k_step, int N_max)
{
	*length_out = N_max/k_step;
	*out_array = calloc(*length_out,sizeof(double));
	int k1;
	for(k1=0; k1 < *length_out; k1++)
	{
		(*out_array)[k1] = in_array[k1*k_step];
	}
	return;

}

/**
 * @brief Returns a linear interpolation of a point x, i.e. given x it returns 
 * y(x).
 * 
 * @details If x outside the boundary given by x_grid, it returns the boundary values
 * of y_grid. 
 * 
 * @param n Grid size.
 * @param x Interpolated point.
 * @param x_grid X-axis where the point x is located.
 * @param y_grid Y(X) array.
 * @return double y_interp(x)
 */
double interpolate(int n, double x, double *x_grid, double* y_grid) 
{
	int k1,k2;
	double y;
	
	k1=0;
	k2=0;
	findinterval(n, x, x_grid, &k1, &k2);
	if( k1 == -1 )
	{
		y=y_grid[0];
	} else if( k2 == n+1 ){
		y=y_grid[n];
	} else {
		y=y_grid[k1]+(x-x_grid[k1])*(y_grid[k2]-y_grid[k1])/(x_grid[k2]-x_grid[k1]);	
	}
	
	return y;
}

/**
 * @brief Finds next closest zero, i.e. x for which y(x) = 0, according to an input array.
 * 
 * @param n Grid size.
 * @param x Value from which to start searching for zero.
 * @param x_grid X-grid where the zeros are located.
 * @param y_grid Y(X) array.
 * @return double 
 */
double findnextinterpolatedzero(int n, double x, double * x_grid, double * y_grid) 
{
	int k1,k2,k4;
	double x_root,x1,x2,y1,y2;
	
	k1=0;
	k2=0;
	
	findinterval(n, x, x_grid, &k1, &k2);

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
	
	return 0.;
}
