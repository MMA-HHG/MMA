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

void nxtval_init(int init_offset, int *val)
{
	*val = init_offset;
}

void nxtval_strided(int stride, int *val)
{
	printf("nxtval1 %i, %i \n",*val, stride); fflush(NULL);
	*val = *val + stride;
	printf("nxtval2 %i, %i \n",*val, stride); fflush(NULL);
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


void coarsen_grid_real(double *in_array, int length_in, double **out_array, int *length_out, int k_step, int N_max)
{
	*length_out = N_max/k_step;
	printf("coarsing, length %i \n",*length_out); fflush(NULL);
	*out_array = calloc(*length_out,sizeof(double));
	int k1;
	for(k1=0; k1 < *length_out; k1++)
	{
		(*out_array)[k1] = in_array[k1*k_step];
	}
	return;

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
