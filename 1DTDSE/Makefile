### Compile 1D TDSE

CC=mpicc
CFLAGS=-Wall -std=c99 -O3
LIBS=-L/usr/local/lib -lhdf5 -lhdf5_hl -lfftw3 -lm
INCLUDE=-I/usr/local/include
SOURCE=sources/constants.c sources/tools_algorithmic.c sources/tools_MPI-RMA.c \
		sources/tools_hdf5.c sources/tridiag.c sources/tools_fftw3.c \
		sources/tools.c sources/prop.c sources/singleTDSE.c \
		sources/strided_scheduler_separated_outs_coarsen_test.c
OBJ=*.o
TDSE=TDSE.e

all: $(TDSE)

$(TDSE): $(SOURCE)
	$(CC) $(CFLAGS) $(INCLUDE) -c $(SOURCE)
	$(CC) $(CFLAGS) $(LIBS) -o $(TDSE) $(OBJ)

clean_obj: *.o
	rm -f *.o

clean:
	rm *.o $(TDSE)