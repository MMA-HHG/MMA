!> @brief This module maintains the density modulation in the case there's
!! a non-trivial density profile. It contains subroutines to initialise and
!! compute the local densitive relative to the reference density \ref
!! parameters::rhoat_inv "rhoat_inv".
!!
!! @author Jan Vábek
!! @author Stefan Skupin
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

REAL(8), dimension(:), allocatable, public    :: density_mod            !< vector containg the radial modulation of the density (respective to the reference density \ref parameters::rhoat_inv "rhoat_inv")
logical, public                               :: apply_density_mod      !< TRUE: the density modulation is applied; FALSE: default value \ref parameters::rhoat_inv "rhoat_inv" is used. Depends on \ref h5namelist::density_mod_grpname "density_mod_grpname".
logical, public                               :: is_density_changed     !< maintans if the density is changed from one z-step to another


!> @cond INCLUDE_DENSITY_MOD_INTERNALS
! module internals
integer                                 :: h5err
real(8), dimension(:), allocatable      :: rgrid, zgrid                 ! local grids for the interpolation of the density modulation and their dimensions ↓↓
integer                                 :: Nr, Nz
real(8), dimension(:,:), allocatable    :: density_profile_matrix       ! matrix to store the density modulation on the input grids ↑↑
!> @endcond 


CONTAINS

!> @brief This subroutine initialises the interpolation procedure te retrieve the density modulation based
!! on the table stored in 'file_id'. It is called if \apply_density_mod in the initial phase of the code.
!!
!! @param[in]       file_id          The id of the main HDF5 file. 
subroutine init_density_mod(file_id)
    integer(hid_t)                          :: file_id
    real(8)                                 :: dumr
    real(8), dimension(:), allocatable      :: dumr_arr_1D
    logical                                 :: rgrid_exists, zgrid_exists

    print *, 'density-mod acessed, proc', my_rank

    call h5lexists_f(file_id, density_mod_grpname//'/zgrid',zgrid_exists,h5err)
    call h5lexists_f(file_id, density_mod_grpname//'/rgrid',rgrid_exists,h5err)
    
    ! We create a matrix in all the cases. If only one grid is provided, max-z (-r) values are used

    if (zgrid_exists .and. rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)

        allocate(rgrid(Nr), zgrid(Nz))

        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)

        rgrid = rgrid/w0m ! convert units [m -> C.U.]  
        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        call read_dset(file_id, density_mod_grpname//'/table', density_profile_matrix, Nr, Nz)

    elseif (zgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/zgrid', Nz)
        allocate(zgrid(Nz))
        call read_dset(file_id, density_mod_grpname//'/zgrid', zgrid, Nz)

        print *, 'zgrid created, myrank', my_rank
        print *, '(density mod) normalisation in z:', four_z_Rayleigh

        zgrid = zgrid/four_z_Rayleigh ! convert units [m -> C.U.]

        rgrid = (/0.d0, delta_r*dim_r/); Nr = 2

        allocate(density_profile_matrix(2,Nz),dumr_arr_1D(Nz))
        call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nz)
        density_profile_matrix(1,:) = dumr_arr_1D
        density_profile_matrix(2,:) = dumr_arr_1D
        deallocate(dumr_arr_1D)

    elseif (rgrid_exists) then
        call ask_for_size_1D(file_id, density_mod_grpname//'/rgrid', Nr)
        allocate(rgrid(Nr))
        call read_dset(file_id, density_mod_grpname//'/rgrid', rgrid, Nr)
        rgrid = rgrid/w0m ! convert units [m -> C.U.]

        zgrid = (/0.d0, proplength/); Nz = 2

        allocate(density_profile_matrix(Nr,2),dumr_arr_1D(Nr))
        call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nr)
        density_profile_matrix(:,1) = dumr_arr_1D
        density_profile_matrix(:,2) = dumr_arr_1D
        deallocate(dumr_arr_1D)

    else
        error stop 'density_mod group present but grids are wrongly specified'
    endif

end subroutine init_density_mod



!> @brief This subroutine computes the radial density modulation for a given 'z'. The result is
!! stored in \ref density_mod
!!
!! @param[in]       z          The local 'z'-coordinate.
subroutine calc_density_mod(z)

    real(8) :: z
    integer, save       :: kz_guess = 1
    integer             :: k1, kr, kz, kr_guess
    real(8)             :: r, density_dum
    logical, save       :: first_call = .false.

    ! testing
    real(8), allocatable :: density_mod_compare(:)
    integer(4)           :: mpi_snd, mpi_rcv


    is_density_changed = .false.

    if (first_call) then
        density_mod = 1.d0
        is_density_changed = .true.
        first_call = .false.
    endif

    call findinterval(kz,z,zgrid,k_guess=kz_guess)
    kz_guess = kz ! see the save attribute
    
    kr_guess = 1
    do k1 = 1, dim_r
        r = (k1-1)*delta_r
        call findinterval(kr,r,rgrid,k_guess=kr_guess)

        call interpolate2D_lin(r,z,density_dum,rgrid,zgrid,density_profile_matrix,Nr,Nz,kx_known=kr,ky_known=kz)
        if (density_dum /= density_mod(k1)) then
            density_mod(k1) = density_dum
            is_density_changed = .true.
        endif
        kr_guess = kr
    enddo

    !!! TEST MPI, compare accross workers
    allocate(density_mod_compare(dim_r))

    call MPI_SENDRECV( &
                        density_mod,         dim_r, MPI_DOUBLE_PRECISION, mod(my_rank            + 1, num_proc), 0, & ! send
                        density_mod_compare, dim_r, MPI_DOUBLE_PRECISION, mod(my_rank + num_proc - 1, num_proc), 0, & ! receive
                        MPI_COMM_WORLD, MPI_STATUS_IGNORE, ierr)

    if ( any( density_mod /= density_mod_compare ) ) then
        print *, "density modulation does not exacly match for rank", my_rank, "and z[m] =", z*four_z_Rayleigh
        error stop "ERROR IN THE DENSITY MODULATION"
    else
        if (my_rank == 0) then
            print *, "density modulation test passed for z[m] =", z*four_z_Rayleigh
            print *, "density modulation (first, middle, last)", (/ density_mod(1), density_mod(dim_r/2), density_mod(dim_r) /)
        endif
    endif



    deallocate(density_mod_compare)

! NOTE: The procedure is written universally to always interpolate. Can be optimised, e.g. do not recompute at all in the case of purely r-modulation.
end subroutine calc_density_mod


end module density_module