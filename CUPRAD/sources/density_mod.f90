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

        rgrid = (/0.d0, delta_r*dim_r/)

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

        zgrid = (/0.d0, proplength/)

        allocate(density_profile_matrix(Nr,2),dumr_arr_1D(Nr))
        call read_dset(file_id, density_mod_grpname//'/table', dumr_arr_1D, Nr)
        density_profile_matrix(:,1) = dumr_arr_1D
        density_profile_matrix(:,2) = dumr_arr_1D
        deallocate(dumr_arr_1D)

    else
        error stop 'density_mod group present but grids are wrongly specified'
    endif

end subroutine init_density_mod

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

    call findinterval(kr,kz,r,z,rgrid,zgrid,Nr,Nz,kx_guess=kr_guess,ky_guess=kz_guess)
    kz_guess = kz ! see the save attribute
    
    kr_guess = 1
    do k1 = 1, dim_r
        r=(k1-1)*delta_r
        call findinterval(kr,r,rgrid,Nr,k_guess=kr_guess)

        call interpolate2D_decomposed_eq(kr,kz,r,z,density_dum,rgrid,zgrid,density_profile_matrix,Nr,Nz)
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
        print *, "density modulation does not exacly match for rank", my_rank, "and z[m] =", z/four_z_Rayleigh
        error stop "ERROR IN THE DENSITY MODULATION"
    else
        if (my_rank == 0) then
            print *, "density modulation test passed for z[m] =", z/four_z_Rayleigh
            print *, "density modulation (first, middle, last)", (/ density_mod(1), density_mod(dim_r/2), density_mod(dim_r) /)
        endif
    endif



    deallocate(density_mod_compare)

! NOTE: The procedure is written universally to always interpolate. Can be optimised, e.g. do not recompute at all in the case of purely r-modulation.
end subroutine calc_density_mod


end module density_module

! call mpi_isend(density_mod, ierr)
! call mpi_recv(density_mod_compare, ierr)

! call mpi_wait()
! call mpi_wait()
