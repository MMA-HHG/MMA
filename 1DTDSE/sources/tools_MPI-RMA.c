/**
 * @file tools_MPI-RMA.c
 * @brief Contains MPI routines for MPI tasks management.
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#include "tools_MPI-RMA.h"

void nxtval_init(int init_offset, int *val)
{
	*val = init_offset;
}

void nxtval_strided(int stride, int *val)
{
	*val = *val + stride;
}