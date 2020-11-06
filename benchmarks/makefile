#
#       Definitions

  MPIF90 = mpiifort

  F90    = ifort

  FFLAGS =  -check all -warn all -traceback -fpe0

  LFLAGS = -lhdf5 -lhdf5_hl -lhdf5_fortran -lhdf5hl_fortran -lz 

  FFTPATH    = $(MKL_LIBS_MT)
  

  TEST_ARRAY	   = test_array.e
  S_TEST_ARRAY       = ../sources/array_helper.f90 test_arrays.f90

  TEST_ARRAY2	   = test_interp.e
  S_TEST_ARRAY2       = ../sources/array_helper.f90 test_interps.f90


#       Make definition
#
linked_list : linked_list_module.f90
	$(F90) -c linked_list_module.f90


$(TEST_ARRAY) : $(S_TEST_ARRAY)
	${F90} -o $(TEST_ARRAY) $(S_TEST_ARRAY) ${FMFLAGS} ${LFLAGS}

$(TEST_ARRAY2) : $(S_TEST_ARRAY2)
	${F90} -o $(TEST_ARRAY2) $(S_TEST_ARRAY2) ${FMFLAGS} ${LFLAGS}


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
clean_preprocessor :
	rm -f listing STOP *.DAT *.LOG
#
test_find : $(TEST_ARRAY) # binary # clean

test_interp : $(TEST_ARRAY2) # binary # clean
