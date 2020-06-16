#
#       Definitions
#
#CINES correction   MPIF90 = mpif90
#CINES commentaire : mpif90 est A utiliser avec openmpi
#CINES commentaire : mpiifort est A utiliser avec intelmpi
#  MPIF90 = mpiifort
  MPIF90 = mpif90
#  F90    = ifort
  F90    = gfortran
#  FFLAGS =  -check all -warn all -traceback -fpe0
#  FFLAGS =  -Wall -O2 -ffree-line-length-none -fbacktrace -fcheck=all -v -g3 # -g[0-3] is for debuging
  FFLAGS =  -Wall -O0 -ffree-line-length-none -fbacktrace -fcheck=all -Wuninitialized -v -g -fbounds-check
#  FFLAGS =  -O3 -ansi_alias 
#  LFLAGS = -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lz 
#  LFLAGS = -lfftw3_mpi -lfftw3 -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lz -lm
  LDLIBS = -lfftw3_mpi -lfftw3 -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lz -lm
#  FFTPATH    = $(MKL_LIBS_MT)
  FFTPATH    = /gpfs/apps/tools/modulefiles/software/FFTW/3.3.5-OpenMPI-1.10.2-GCC-6.1.0-psm-dbg/lib/libfftw3.so
  H5INCLUDE  = /gpfs/apps/tools/modulefiles/software/HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm/include
  H5LIBS = /gpfs/apps/tools/modulefiles/software/HDF5/1.8.13-OpenMPI-1.10.2-GCC-6.1.0-psm/lib
#  FFT     = -lfftw3
  PROPAGATION       = cuprad_occigen.e
  S_PROPAGATION       =  modules.f90 code_ionisation.f90 fft_rk.f90 finalize.f90 \
                        firststep.f90 longstep_rk.f90 output.f90 outputHDF5_occigen.f90 cuprad.f90
  MAKE_START	   = make_start_occigen.e
  S_MAKE_START       = write_start.f90 calc_start.f90 normalisation.f90 \
                       write_listing.f90 make_start.f90
  CLEAN       = clean_occigen.e
  S_CLEAN	= clean.f90
  MERGE       = merge_occigen.e
  S_MERGE	= merge.f90
  CALC_FT	= calc_ft_occigen.e
  S_CALC_FT	  = modules.f90 fft.f90 calc_ft.f90


#
#       Make definition
#
$(PROPAGATION) : $(S_PROPAGATION)
#	${MPIF90} -o $(PROPAGATION) $(S_PROPAGATION) ${FFTPATH} ${FFT} ${FFLAGS} ${LFLAGS}
	${MPIF90} -I${H5INCLUDE} -L${H5LIBS} -o $(PROPAGATION) $(S_PROPAGATION) ${FFTPATH} ${FFT} ${FFLAGS} ${LFLAGS} ${LDLIBS}
#	${MPIF90} -o $(PROPAGATION) $(S_PROPAGATION) ${FFTPATH} ${FFT} ${FFLAGS}
$(PROPAGATION_DEFOK) : $(S_PROPAGATION_DEFOK)
#	${MPIF90} -o $(PROPAGATION_DEFOK) $(S_PROPAGATION_DEFOK) ${FFTPATH} ${FFT} ${FFLAGS}
	${MPIF90} ${FFLAGS} -o $(PROPAGATION_DEFOK) $(S_PROPAGATION_DEFOK) ${FFTPATH} ${FFT} ${FFLAGS}
$(PROPAGATION_MICKAEL) : $(S_PROPAGATION_MICKAEL)
	${MPIF90} -o $(PROPAGATION_MICKAEL) $(S_PROPAGATION_MICKAEL) ${FFTPATH} ${FFT} ${FFLAGS}
$(MAKE_START) : $(S_MAKE_START)
#	${F90} -o $(MAKE_START) $(S_MAKE_START) ${FFTPATH} ${FFT} ${FMFLAGS}
	${F90} ${FFLAGS} -o $(MAKE_START) $(S_MAKE_START) ${FFTPATH} ${FFT} ${FMFLAGS}
$(MAKE_START_NU0) : $(S_MAKE_START_NU0)
	${F90} -o $(MAKE_START_NU0) $(S_MAKE_START_NU0) ${FFTPATH} ${FFT} ${FFLAGS}
$(CLEAN) : $(S_CLEAN)
	${F90} -o $(CLEAN) $(S_CLEAN) ${FFLAGS}
$(MERGE) : $(S_MERGE)
	${F90} -o $(MERGE) $(S_MERGE) ${FFLAGS}
$(CALC_FT) : $(S_CALC_FT)
	${MPIF90} -o $(CALC_FT) $(S_CALC_FT) ${FFTPATH} ${FFT} ${FFLAGS}
#
#      Utilities
#
clean :
	/bin/rm -f *.e
	/bin/rm -f *.o
	/bin/rm -f *.mod
#
binary :
	/bin/mv *.e ../binary
#
code : $(CALC_FT) $(MAKE_START) $(MAKE_START_4WM) $(PROPAGATION) $(PROPAGATION_DEFOK) $(PROPAGATION_MICKAEL) $(CLEAN) $(MERGE) binary clean
#
util : $(CALC_FT) $(MAKE_START) $(MAKE_START_NU0) $(CLEAN) $(MERGE) binary clean
