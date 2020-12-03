PROGRAM test_modules

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  integer :: iarray4(4) = (/1, 2, 3, 4/)
  logical :: dumlog




  print *,"Test any mask on integers"

  
  print *, .true.

  dumlog = ANY(iarray4 == 1)

  print *,"arr ni 1", dumlog

  dumlog = ANY(1 == iarray4 )

  print *,"1 in arr", dumlog

  dumlog = ANY(iarray4 == 5)

  print *,"arr ni 5", dumlog

  dumlog = ANY(5 == iarray4 )

  print *,"5 in arr", dumlog  


END PROGRAM test_modules
