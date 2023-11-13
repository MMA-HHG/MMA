!> @author Jan Vábek
!! @brief module for pre-ionisation
!!
!! This module contains ...
!! 
module pre_ionised
! It was developed by Jan Vabek
! 
! method
! 1 - constant (ratio compared to neutrals)
! 2 - full-table
! 3 - r-table
! 4 - z-table
! It replaces the original rho0 in the code

! units
! 1 - relative
! 2 - cm-3

use HDF5
use HDF5_helper
use array_helper
use parameters
use fields
use normalization
use h5namelist

use density_module

use mpi_stuff ! only for testing ot print procnumber etc., remove after

implicit none
private
public  :: init_pre_ionisation, initial_electron_density !, initial_electron_density_guess


logical, public                         :: apply_pre_ionisation
integer                                 :: method_geometry, method_units
integer                                 :: Nr, Nz
integer                                 :: h5err
real(8), dimension(:), allocatable      :: rgrid
real(8), dimension(:), allocatable      :: zgrid
real(8), dimension(:,:), allocatable    :: initial_electrons_ratio_matrix
real(8)                                 :: initial_electrons_ratio_scalar
integer, parameter, dimension(3)        :: table_geometries = (/2,3,4/)
integer, parameter, dimension(2)        :: table_1D_geometries = (/3,4/)
logical                                 :: scalar_case = .false.



CONTAINS
! preparation
subroutine init_pre_ionisation(file_id)
    integer(hid_t)              :: file_id
    real(8)                     :: dumr

    real(8), dimension(:), allocatable      :: dumr_arr_1D

    logical                                 :: rgrid_exists, zgrid_exists, initial_electrons_ratio_exists

    print *, 'pre-inoisation accessed, proc', my_rank


    call h5lexists_f(file_id, density_mod_grpname//'/zgrid',zgrid_exists,h5err)
    call h5lexists_f(file_id, density_mod_grpname//'/rgrid',rgrid_exists,h5err)
    call h5lexists_f(file_id, density_mod_grpname//'/initial_electrons_ratio',initial_electrons_ratio_exists,h5err)
    
    ! We create a matrix in all the cases. If only one grid is provided, max-z (-r) values are used

    if (zgrid_exists.and. rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        rgrid = rgrid/w0m ! convert units [m -> C.U.]
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)
        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        call read_dset(file_id, density_mod_grpname//'/initial_electrons_ratio', initial_electrons_ratio_matrix, Nr, Nz)

        initial_electrons_ratio_matrix = initial_electrons_ratio_matrix/rhoat_inv ! convert units


    elseif (zgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)
        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        rgrid = (/0.d0, delta_r*dim_r/); Nr = 2

        allocate(initial_electrons_ratio_matrix(2,Nz))
        call read_dset(file_id, density_mod_grpname//'/initial_electrons_ratio', dumr_arr_1D, Nz)
        initial_electrons_ratio_matrix(1,:) = dumr_arr_1D/rhoat_inv
        initial_electrons_ratio_matrix(2,:) = dumr_arr_1D/rhoat_inv

    elseif (rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        rgrid = rgrid/w0m ! convert units [m -> C.U.]

        zgrid = (/0.d0, proplength/); Nz = 2

        allocate(initial_electrons_ratio_matrix(Nr,2))
        call read_dset(file_id, density_mod_grpname//'/initial_electrons_ratio', dumr_arr_1D, Nr)
        initial_electrons_ratio_matrix(:,1) = dumr_arr_1D/rhoat_inv
        initial_electrons_ratio_matrix(:,2) = dumr_arr_1D/rhoat_inv

    elseif (initial_electrons_ratio_exists) then        
        call read_dset(file_id, pre_ionised_grpname//'/initial_electrons_ratio',initial_electrons_ratio_scalar)
        initial_electrons_ratio_scalar = initial_electrons_ratio_scalar/rhoat_inv
        scalar_case = .true.
        
    else
        error stop 'pre-ionisation group present, but inputs therein are wrongly specified'
    endif

end subroutine init_pre_ionisation


function initial_electron_density(r,z,kr_actual,kr_first) ! already rescaled to C.U.
    real(8)                    :: initial_electron_density
    real(8)                    :: r,z 
    integer                    :: kr_actual, kr_first
    integer                    :: kr,kz    

    logical, save              :: first_run = .true.
    integer, save              :: my_kr_start, kr_guess, kz_guess = 1

    if (scalar_case) then
        initial_electron_density = density_mod(kr_actual)*initial_electrons_ratio_scalar
        return
    endif

    if (first_run) then
        if (kr_actual == kr_first) then
            call findinterval(my_kr_start,r,rgrid,Nr) ! obtain guess, assume it's first called at right place
        else
            error stop "initial_electron_density requires to be firstly called for kr_actual = kr_first"
        endif
        first_run = .false.
    endif

    if (kr_actual == kr_first) then
        kr = my_kr_start
    else
        call findinterval(kr,r,rgrid,Nr,k_guess = kr_guess)
    endif

    kr_guess = kr ! bookkeeping: remember kr from the last run

    call findinterval(kz,z,zgrid,Nz,k_guess = kz_guess)
    kz_guess = kz

    call interpolate2D_decomposed_eq(kr,kz,r,z,initial_electron_density,rgrid,zgrid,initial_electrons_ratio_matrix,Nr,Nz)
    initial_electron_density = density_mod(kr_actual)*initial_electron_density
    return

end function


end module pre_ionised
