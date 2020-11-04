! It was developed by Jan Vabek
! 
! method
! 1 - constant (ratio compared to neutrals)
! 2 - full-table
! 3 - r-table
! 4 - z-table
! It replaces the original rho0 in the code

module pre_ionised
use HDF5, HDF5_helper
use parameters, normalization

implicit none
private
public  :: init_pre_ionisation, initial_electron_density


logical, public                         :: apply_pre_ionisation
integer                                 :: method_geometry, method_units
integer                                 :: Nr, Nz
real(8), dimension(:), allocatable      :: rgrid
real(8), dimension(:), allocatable      :: zgrid
real(8), dimension(:,:), allocatable    :: table_2D
real(8), dimension(:), allocatable      :: table_1D
real(8)                                 :: rho0
integer, parameter, dimension(3)        :: table_geometries = (\2,3,4\)
integer, parameter, dimension(2)        :: table_1D_geometries = (\3,4\)

CONTAINS
! preparation
subroutine init_pre_ionisation(file_id)
    integer(hid_t)              :: file_id
    real(8)                     :: dumr

    call read_dset(file_id, 'pre_ionised/method_geometry',method_geometry)
    select case (method_geometry)
        case (1)
            if (method_units == 1) then
                call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
                call read_dset(file_id, 'pre_ionised/initial_electrons_ratio',rho0)
                rho0 = rho0/dumr ! C.U.
                return
            elseif (method_units == 2) then
                call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
                call read_dset(file_id, 'pre_ionised/initial_electron_density',rho0) ! in cm-3
                rho0 = dumr*rho0 ! C.U.
                return
            endif
        case (2)
            call ask_for_size_1D(group_id, 'pre_ionised/rgrid', Nr)
            call read_dset(group_id, 'pre_ionised/rgrid', rgrid, Nr)
            rgrid = rgrid/w0m ! m -> C.U.
            call ask_for_size_1D(group_id, 'pre_ionised/zgrid', Nz)
            call read_dset(group_id, 'pre_ionised/zgrid', zgrid, Nz)
            zgrid = zgrid/four_z_Rayleigh ! m -> C.U.
            call read_dset(group_id, 'pre_ionised/table', table_2D, Nr, Nz)
    end select

    if (method_units == 1)  then ! any( table_geometries == method_geometry) eventual condition for extended prescriptions
        call read_dset(file_id, 'pre-processed/rhoat_inv',dumr) ! inverse of the neutrals density, C.U.
        dumr = 1.D0/dumr ! conversion factor (- -> C.U.)
    elseif (method_units == 2) then
        call read_dset(file_id, 'pre-processed/density_normalisation_factor',dumr) ! conversion factor (cm-3 -> C.U.)
    endif
end subroutine init_pre_ionisation

!
function initial_electron_density(r,z) ! already rescaled to C.U.
    real(8)                    :: initial_electron_density
    real(8)                    :: r,z 
    ! each worker should know its first interpolation guess

    select case (method_geometry)
        case (1)
            initial_electron_density = rho0;
            return
    end select
end function


end module pre_ionised