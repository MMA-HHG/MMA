PROGRAM test_modules
  USE array_helper
  USE HDF5
  USE HDF5_helper_serial

  IMPLICIT NONE
  
  integer, parameter :: Np = 100;
  integer :: k1,k2
  real(8), parameter :: default_array(Np) = (/ (k1, k1=1,Np) /)
  real(8), parameter :: default_array0(Np) = (/ (k1, k1=0,Np-1) /)
  real(8), parameter :: default_array3(3) = (/ (k1, k1=1,3) /), default_array4(4) = (/ (k1, k1=1,4) /)

  real(8), parameter :: x_test(11) = (/ -1.0D0, 1.0D0, 1.5D0, 2.0D0, 2.5D0, 2.999D0, 3.5D0, 57.5D0, 99.1D0, 100.0D0, 101.0D0 /)



  real(8), dimension(:), allocatable      :: rgrid
  real(8), dimension(:), allocatable      :: zgrid
  real(8), dimension(:,:), allocatable    :: density_profile_matrix

  real(8), dimension(:), allocatable      :: dumr_arr_1D

  real(8) :: array2d(Np,Np), x, y, fx, fxy

  character(*), parameter   ::  main_h5_fname = "results.h5" , density_mod_grpname =   "density_mod"

  real(8), parameter :: zmax_default = 1., rmax_default = 1.

  real(8)   :: r, z, density_dum
  integer   :: Nr, Nz, kr, kz
  integer   :: h5err

  logical zgrid_exists, rgrid_exists

  integer(hid_t)              :: file_id



  CALL h5open_f(h5err)
  CALL h5fopen_f (main_h5_fname, H5F_ACC_RDONLY_F, file_id, h5err)
  call h5lexists_f(file_id, density_mod_grpname//'/zgrid',zgrid_exists,h5err)
  call h5lexists_f(file_id, density_mod_grpname//'/rgrid',rgrid_exists,h5err)

  print *, 'zgrid exists', zgrid_exists, 'rgrid exists', rgrid_exists
  
  ! We create a matrix in all the cases. If only one grid is provided, max-z (-r) values are used

  if (zgrid_exists .and. rgrid_exists) then
      call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
      call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)

      allocate(rgrid(Nr), zgrid(Nz))

      call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
      call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)

    !   rgrid = rgrid/w0m ! convert units [m -> C.U.]  
    !   zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

      call read_dset(file_id, density_mod_grpname//'/table', density_profile_matrix, Nr, Nz)

  elseif (zgrid_exists) then
      call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
      print *, 'zsize is', Nz
      allocate(zgrid(Nz))
      call read_array_real_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)
      print *, 'zgrid is', zgrid


      

    !   print *, 'zgrid created, myrank', my_rank
    !   print *, '(density mod) normalisation in z:', four_z_Rayleigh

    !   zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

      rgrid = (/0.d0, rmax_default/); Nr = 2

      allocate(density_profile_matrix(2,Nz),dumr_arr_1D(Nz))
      call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nz)
      density_profile_matrix(1,:) = dumr_arr_1D
      density_profile_matrix(2,:) = dumr_arr_1D
      deallocate(dumr_arr_1D)

  elseif (rgrid_exists) then
      call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
      allocate(rgrid(Nr))
      call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
    !   rgrid = rgrid/w0m ! convert units [m -> C.U.]

      zgrid = (/0.d0, zmax_default/)

      allocate(density_profile_matrix(Nr,2),dumr_arr_1D(Nr))
      call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nr)
      density_profile_matrix(:,1) = dumr_arr_1D
      density_profile_matrix(:,2) = dumr_arr_1D
      deallocate(dumr_arr_1D)

  else
      error stop 'unable to read data'
  endif

  CALL h5fclose_f(file_id, h5err)
  CALL h5close_f(h5err)

  print *, '--------------------------------------------------------------'


  do while (.true.)
    print *, '--------------------------------------------------------------'
    print *, 'z-value:'
    read *, z
    print *, 'r-value:'
    read *, r
    call findinterval(kz,z,zgrid,Nz)
    call findinterval(kr,r,rgrid,Nr)

    call interpolate2D_decomposed_eq(kr,kz,r,z,density_dum,rgrid,zgrid,density_profile_matrix,Nr,Nz)
    print *, 'the indexes are', (/kr, kz/)
    print *, 'the input is', (/r,z/)
    print *, 'the value is', density_dum

  enddo


!   do k1 = 0, (Np-1)
!     array2d(:,k1+1) = k1 + default_array0
!   enddo



!   print *,"Series of tests of the array procedures"

  
!   print *, default_array

!   print *,"notip"

!   do k1 = 1,11
!     call interpolate_lin(x_test(k1),fx,default_array,default_array+1,Np)
!     print *, x_test(k1), fx  
!   enddo

!   print *,"tip"
!   do k1 = 1,11
!     call interpolate_lin(x_test(k1),fx,default_array,default_array+2,Np, k_guess = 50)
!     print *, x_test(k1), fx  
!   enddo


!   print *,"2D"
!   call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np)
!   print *, 1.5D0, 1.5D0, fxy 

!   call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, kx_guess = 50)
!   print *, 1.5D0, 1.5D0, fxy 

!   call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, ky_guess = 50)
!   print *, 1.5D0, 1.5D0, fxy 

!   call interpolate_lin(1.5D0,1.5D0,fxy,default_array0,default_array0,array2d,Np,Np, kx_guess = 50, ky_guess = 50)
!   print *, 1.5D0, 1.5D0, fxy 

!   call interpolate_lin(-1.0D0,-1.0D0,fxy,default_array0,default_array0,array2d,Np,Np)
!   print *, -1.0D0, -1.0D0, fxy 

!   call interpolate_lin(-1.0D0, 2.5D0,fxy,default_array0,default_array0,array2d,Np,Np)
!   print *, -1.0D0, 2.5D0, fxy 

!   call interpolate_lin( 2.5D0,-1.0D0,fxy,default_array0,default_array0,array2d,Np,Np)
!   print *, 2.5D0, -1.0D0, fxy 


END PROGRAM test_modules
