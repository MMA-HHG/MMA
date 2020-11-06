PROGRAM test_modules
  USE array_helper

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  real(8), parameter :: default_array(Np) = (/ (k1, k1=1,Np) /)
  real(8), parameter :: default_array3(3) = (/ (k1, k1=1,3) /), default_array4(4) = (/ (k1, k1=1,4) /)



  print *,"Series of tests of the array procedures"

  k1 = 0
  print *, default_array
  call findinterval(k1,56.5D0,default_array,Np,57)
  print *, k1

END PROGRAM test_modules
