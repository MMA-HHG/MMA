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

module density_module
use HDF5
use HDF5_helper
use array_helper
use parameters
use normalization

use mpi_stuff ! only for testing ot print procnumber etc., remove after

implicit none
private
public  :: init_density_mod,calc_density_mod


logical, public                         :: apply_density_mod,is_density_changed
integer                                 :: method_geometry, method_units
integer                                 :: Nr, Nz
real(8), dimension(:), allocatable      :: rgrid
real(8), dimension(:), allocatable      :: zgrid
real(8), dimension(:,:), allocatable    :: table_2D
real(8), dimension(:), allocatable      :: table_1D
REAL(8), dimension(:), ALLOCATABLE, public    :: density_mod
real(8)                                 :: rho0_loc ! there is rho0 in the global scope (parameters)
integer, parameter, dimension(3)        :: table_geometries = (/2,3,4/)
integer, parameter, dimension(2)        :: table_1D_geometries = (/3,4/)

CONTAINS
! preparation
subroutine init_density_mod(file_id)
    integer(hid_t)              :: file_id
    real(8)                     :: dumr

    print *, 'pre-inoisation accessed, proc', my_rank

    call read_dset(file_id, 'pre_ionised/method_geometry',method_geometry)
    call read_dset(file_id, 'pre_ionised/method_units',method_units)
    select case (method_geometry)
        case (1)
            if (method_units == 1) then
                call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
                call read_dset(file_id, 'pre_ionised/initial_electrons_ratio',rho0_loc)
                rho0_loc = rho0_loc/dumr ! C.U.
                print *, 'the initial atom density is',1.0D0/dumr, 'proc', my_rank
                print *, 'the initial density is',rho0_loc, 'proc', my_rank
                return
            elseif (method_units == 2) then
                call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
                call read_dset(file_id, 'pre_ionised/initial_electron_density',rho0_loc) ! in cm-3
                rho0_loc = dumr*rho0_loc ! C.U.
                return
            endif
        case (2)
            call ask_for_size_1D(file_id, 'pre_ionised/rgrid', Nr)
            call read_dset(file_id, 'pre_ionised/rgrid', rgrid, Nr)
            rgrid = rgrid/w0m ! m -> C.U.
            call ask_for_size_1D(file_id, 'pre_ionised/zgrid', Nz)
            call read_dset(file_id, 'pre_ionised/zgrid', zgrid, Nz)
            zgrid = zgrid/four_z_Rayleigh ! m -> C.U.
            call read_dset(file_id, 'pre_ionised/table', table_2D, Nr, Nz)
        case (3)
            call ask_for_size_1D(file_id, 'pre_ionised/rgrid', Nr)
            call read_dset(file_id, 'pre_ionised/rgrid', rgrid, Nr)
            rgrid = rgrid/w0m ! m -> C.U.
            call read_dset(file_id, 'pre_ionised/table', table_1D, Nr)
        case (4)
            call ask_for_size_1D(file_id, 'pre_ionised/zgrid', Nz)
            call read_dset(file_id, 'pre_ionised/zgrid', zgrid, Nz)
            zgrid = zgrid/four_z_Rayleigh ! m -> C.U.
            call read_dset(file_id, 'pre_ionised/table', table_1D, Nz)
    end select

    !if (method_units == 1)  then ! any( table_geometries == method_geometry) eventual condition for extended prescriptions
    !    call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
    !    dumr = 1.D0/dumr ! conversion factor (- -> C.U.)
    !elseif (method_units == 2) then
    !    call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
    !endif
end subroutine init_density_mod

subroutine calc_density_mod(z)

REAL(8) :: z

end subroutine calc_density_mod


end module density_module
