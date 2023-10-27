!>
!!
!!
!!
!!



module density_module
use HDF5
use HDF5_helper
use array_helper
use normalization
use h5namelist
use parameters
use fields

use mpi_stuff ! only for testing ot print procnumber etc., remove after

implicit none
private

! public
public  :: init_density_mod,calc_density_mod

REAL(8), dimension(:), ALLOCATABLE, public    :: density_mod
logical, public                         :: apply_density_mod,is_density_changed

! module internals
integer                                 :: method_geometry
integer                                 :: Nr, Nz
integer                                 :: h5err
real(8), dimension(:), allocatable      :: rgrid
real(8), dimension(:), allocatable      :: zgrid
real(8), dimension(:,:), allocatable    :: table_2D
real(8), dimension(:,:), allocatable    :: density_profile_matrix


integer, parameter, dimension(3)        :: table_geometries = (/2,3,4/)
integer, parameter, dimension(2)        :: table_1D_geometries = (/3,4/)

CONTAINS

! preparation
subroutine init_density_mod(file_id)
    integer(hid_t)                          :: file_id
    real(8)                                 :: dumr
    real(8), dimension(:), allocatable      :: dumr_arr_1D
    logical                                 :: rgrid_exists, zgrid_exists

    print *, 'density-mod acessed, proc', my_rank

    call h5lexists_f(file_id, density_mod_grpname//'/zgrid',zgrid_exists,h5err)
    call h5lexists_f(file_id, density_mod_grpname//'/rgrid',rgrid_exists,h5err)
    
    ! We create a matrix in all the cases. If only one grid is provided, max-z (-r) values are used

    if (zgrid_exists.and. rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        rgrid = rgrid/w0m ! convert units [m -> C.U.]
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)
        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        call read_dset(file_id, density_mod_grpname//'/table', density_profile_matrix, Nr, Nz)

    elseif (zgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)
        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        rgrid = (/0.d0, delta_r*dim_r/) ! automatic allocation test

        allocate(density_profile_matrix(2,Nz))
        call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nz)
        density_profile_matrix(1,:) = dumr_arr_1D
        density_profile_matrix(2,:) = dumr_arr_1D

    elseif (rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        rgrid = rgrid/w0m ! convert units [m -> C.U.]

        zgrid = (/0.d0, proplength/) ! automatic allocation test

        allocate(density_profile_matrix(Nr,2))
        call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nr)
        density_profile_matrix(:,1) = dumr_arr_1D
        density_profile_matrix(:,2) = dumr_arr_1D

    else
        error stop 'density_mod group present but grids are wrongly specified'
    endif

    ! select case (method_geometry)
    !     case (1)
    !         if (method_units == 1) then
    !             call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
    !             call read_dset(file_id, 'pre_ionised/initial_electrons_ratio',rho0_loc)
    !             rho0_loc = rho0_loc/dumr ! C.U.
    !             print *, 'the initial atom density is',1.0D0/dumr, 'proc', my_rank
    !             print *, 'the initial density is',rho0_loc, 'proc', my_rank
    !             return
    !         elseif (method_units == 2) then
    !             call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
    !             call read_dset(file_id, 'pre_ionised/initial_electron_density',rho0_loc) ! in cm-3
    !             rho0_loc = dumr*rho0_loc ! C.U.
    !             return
    !         endif
    !     case (2)
    !         call ask_for_size_1D(file_id, 'pre_ionised/rgrid', Nr)
    !         call read_dset(file_id, 'pre_ionised/rgrid', rgrid, Nr)
    !         rgrid = rgrid/w0m ! m -> C.U.
    !         call ask_for_size_1D(file_id, 'pre_ionised/zgrid', Nz)
    !         call read_dset(file_id, 'pre_ionised/zgrid', zgrid, Nz)
    !         zgrid = zgrid/four_z_Rayleigh ! m -> C.U.
    !         call read_dset(file_id, 'pre_ionised/table', table_2D, Nr, Nz)
    !     case (3)
    !         call ask_for_size_1D(file_id, 'pre_ionised/rgrid', Nr)
    !         call read_dset(file_id, 'pre_ionised/rgrid', rgrid, Nr)
    !         rgrid = rgrid/w0m ! m -> C.U.
    !         call read_dset(file_id, 'pre_ionised/table', table_1D, Nr)
    !     case (4)
    !         call ask_for_size_1D(file_id, 'pre_ionised/zgrid', Nz)
    !         call read_dset(file_id, 'pre_ionised/zgrid', zgrid, Nz)
    !         zgrid = zgrid/four_z_Rayleigh ! m -> C.U.
    !         call read_dset(file_id, 'pre_ionised/table', table_1D, Nz)
    ! end select

    !if (method_units == 1)  then ! any( table_geometries == method_geometry) eventual condition for extended prescriptions
    !    call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
    !    dumr = 1.D0/dumr ! conversion factor (- -> C.U.)
    !elseif (method_units == 2) then
    !    call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
    !endif
end subroutine init_density_mod

subroutine calc_density_mod(z)

    real(8) :: z
    integer, save       :: kz_guess = 1
    integer             :: k1, kr, kz, kr_guess
    real(8)             :: r, density_dum
    logical             :: first_iteration


    is_density_changed = .false.
    first_iteration = .true.

    kr_guess = 1
    do k1 = 1, dim_r
        r=(k1-1)*delta_r
        call findinterval(kr,kz,r,z,rgrid,zgrid,Nr,Nz,kx_guess=kr_guess,ky_guess=kz_guess)
        if (first_iteration) then
            first_iteration = .false.
            kz_guess = kz
        endif

        call interpolate2D_decomposed_eq(kr,kz,r,z,density_dum,rgrid,zgrid,density_profile_matrix,Nr,Nz)
        if (density_dum /= density_mod(k1)) then
            density_mod(k1) = density_dum
            is_density_changed = .true.
        endif
        kr_guess = kr
    enddo


! NOTE: The procedure is written universally to always interpolate. Can be optimised, e.g. do not recompute at all in the case of purely r-modulation.
end subroutine calc_density_mod


end module density_module
