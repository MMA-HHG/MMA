PROGRAM test_modules
  USE array_helper

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  real(8), parameter :: default_array(Np) = (/ (k1, k1=1,Np) /)
  real(8), parameter :: default_array3(3) = (/ (k1, k1=1,3) /), default_array4(4) = (/ (k1, k1=1,4) /)



  print *,"Series of tests of the array procedures"

  
  print *, default_array
  print *,"notip"
  call findinterval(k1,-1.0D0,default_array,Np)
  print *, k1
  call findinterval(k1,1.0D0,default_array,Np)
  print *, k1
  call findinterval(k1,1.5D0,default_array,Np)
  print *, k1
  call findinterval(k1,2.0D0,default_array,Np)
  print *, k1
  call findinterval(k1,2.5D0,default_array,Np)
  print *, k1
  call findinterval(k1,2.999D0,default_array,Np)
  print *, k1
  call findinterval(k1,3.5D0,default_array,Np)
  print *, k1
  call findinterval(k1,57.5D0,default_array,Np)
  print *, k1
  call findinterval(k1,99.1D0,default_array,Np)
  print *, k1
  call findinterval(k1,100.0D0,default_array,Np)
  print *, k1
  call findinterval(k1,101.0D0,default_array,Np)
  print *, k1

  print *,"tip"
  call findinterval(k1,-1.0D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,1.0D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,1.5D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,2.0D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,2.5D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,2.999D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,3.5D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,57.5D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,99.1D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,100.0D0,default_array,Np, k_tip = 50)
  print *, k1
  call findinterval(k1,101.0D0,default_array,Np, k_tip = 50)
  print *, k1

  print *,"2D"
  call findinterval(k1,k2,1.5D0,3.5D0,default_array,default_array,Np,Np)
  print *, k1, k2
  call findinterval(k1,k2,1.5D0,3.5D0,default_array,default_array,Np,Np,kx_tip=25)
  print *, k1, k2
  call findinterval(k1,k2,1.5D0,3.5D0,default_array,default_array,Np,Np,ky_tip=25)
  print *, k1, k2
  call findinterval(k1,k2,1.5D0,3.5D0,default_array,default_array,Np,Np,kx_tip=25,ky_tip=75)
  print *, k1, k2


END PROGRAM test_modules
