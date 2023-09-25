/**
 * @file tools_MPI-RMA.h
 * @brief Header containing MPI routines. 
 * 
 * @copyright Copyright (c) 2023
 * 
 */
#ifndef TOOLS_MPI_RMA_H
#define TOOLS_MPI_RMA_H

#include <time.h> 
#include <stdio.h>
#include <mpi.h>
#include <stdlib.h>
#include <math.h>

void MPE_Counter_create(MPI_Comm, int, MPI_Win *);
int MPE_Counter_nxtval(MPI_Win, int, int *);
int MPE_Mutex_acquire(MPI_Win, int, int);
int MPE_Mutex_release(MPI_Win, int, int);
void nxtval_init(int, int *);
void nxtval_strided(int, int *);

#endif