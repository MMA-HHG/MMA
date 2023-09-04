### Compile 1D TDSE

CC=mpicc
CFLAGS=-Wall -std=c99 -g -O0 
LIBS=-L/usr/local/lib -lhdf5 -lhdf5_hl -lfftw3 -lm
INCLUDE=-I/usr/local/include
SOURCE=sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c \
		sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c sources/structures.c \
		sources/tools.c sources/prop.c sources/singleTDSE.c \
		sources/strided_scheduler_separated_outs_coarsen.c
SOURCE_STRIDE=sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c \
		sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c sources/structures.c \
		sources/tools.c sources/prop.c sources/singleTDSE.c \
		sources/single_stride.c
OBJ=constants.o tools_algorithmic.o tools_MPI-RMA.o \
		tools_hdf5.o tridiag.o tools_fftw3.o structures.o \
		tools.o prop.o singleTDSE.o \
		strided_scheduler_separated_outs_coarsen.o
OBJ_STRIDE=constants.o tools_algorithmic.o tools_MPI-RMA.o \
		tools_hdf5.o tridiag.o tools_fftw3.o structures.o \
		tools.o prop.o singleTDSE.o \
		single_stride.o
TDSE=TDSE.e
TDSE_STRIDE=TDSE_stride.e

all: $(TDSE) $(TDSE_STRIDE)

$(TDSE): $(SOURCE)
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SOURCE)
	$(CC) $(CFLAGS) $(LIBS) -o $(TDSE) $(OBJ)

$(TDSE_STRIDE): $(SOURCE_STRIDE)
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SOURCE_STRIDE)
	$(CC) $(CFLAGS) $(LIBS) -o $(TDSE_STRIDE) $(OBJ_STRIDE)

clean_obj: *.o
	rm -f *.o

clean:
	rm *.o $(TDSE)