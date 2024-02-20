PROGRAM test_modules
  USE array_helper

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  real(8), parameter :: default_array(Np) = (/ (k1, k1=1,Np) /)
  real(8), parameter :: default_array0(Np) = (/ (k1, k1=0,Np-1) /)
  real(8), parameter :: default_array3(3) = (/ (k1, k1=1,3) /), default_array4(4) = (/ (k1, k1=1,4) /)

  real(8), parameter :: x_test(11) = (/ -1.0D0, 1.0D0, 1.5D0, 2.0D0, 2.5D0, 2.999D0, 3.5D0, 57.5D0, 99.1D0, 100.0D0, 101.0D0 /)

  real(8) :: array2d(Np,Np), x, y, fx, fxy

  do k1 = 0, (Np-1)
    array2d(:,k1+1) = k1 + default_array0
  enddo



  print *,"Series of tests of the array procedures"

  
  print *, default_array

  print *,"notip"

  do k1 = 1,11
    call interpolate_lin(x_test(k1),fx,default_array,default_array+1,Np)
    print *, x_test(k1), fx  
  enddo

  print *,"tip"
  do k1 = 1,11
    call interpolate_lin(x_test(k1),fx,default_array,default_array+2,Np, k_guess = 50)
    print *, x_test(k1), fx  
  enddo


  print *,"2D"
  call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np)
  print *, 1.5D0, 1.5D0, fxy 

  call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, kx_guess = 50)
  print *, 1.5D0, 1.5D0, fxy 

  call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, ky_guess = 50)
  print *, 1.5D0, 1.5D0, fxy 

  call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, kx_guess = 50, ky_guess = 50)
  print *, 1.5D0, 1.5D0, fxy 

  call interpolate_lin(-1.0D0,-1.0D0,fxy,default_array0,default_array0,array2d,Np,Np)
  print *, -1.0D0, -1.0D0, fxy 

  call interpolate_lin(-1.0D0, 2.5D0,fxy,default_array0,default_array0,array2d,Np,Np)
  print *, -1.0D0, 2.5D0, fxy 

  call interpolate_lin( 2.5D0,-1.0D0,fxy,default_array0,default_array0,array2d,Np,Np)
  print *, 2.5D0, -1.0D0, fxy 


END PROGRAM test_modules
