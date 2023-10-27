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

module pre_ionised
use HDF5
use HDF5_helper
use array_helper
use parameters
use normalization

use mpi_stuff ! only for testing ot print procnumber etc., remove after

implicit none
private
public  :: init_pre_ionisation, initial_electron_density, initial_electron_density_guess


logical, public                         :: apply_pre_ionisation
integer                                 :: method_geometry, method_units
integer                                 :: Nr, Nz
real(8), dimension(:), allocatable      :: rgrid
real(8), dimension(:), allocatable      :: zgrid
real(8), dimension(:,:), allocatable    :: table_2D
real(8), dimension(:), allocatable      :: table_1D
real(8)                                 :: rho0_loc ! there is rho0 in the global scope (parameters)
integer, parameter, dimension(3)        :: table_geometries = (/2,3,4/)
integer, parameter, dimension(2)        :: table_1D_geometries = (/3,4/)

CONTAINS
! preparation
subroutine init_pre_ionisation(file_id)
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
end subroutine init_pre_ionisation


!
function initial_electron_density(r,z,reset_r_guess) ! already rescaled to C.U.
    real(8)                    :: initial_electron_density
    real(8)                    :: r,z 
    logical, optional          :: reset_r_guess
    integer                    :: kr,kz    

    logical, save              :: first_run = .true.
    integer, save              :: my_first_kr_guess
    integer, save              :: kr_guess, kz_guess

    if (method_geometry == 1) then
        initial_electron_density = rho0_loc;
        return
    endif

    if (first_run) then
        if ( (method_geometry == 2) .or. (method_geometry == 3) ) then
            call findinterval(my_first_kr_guess,r,rgrid,Nr) ! obtain guess, assume it's first called at right place
            kr_guess = my_first_kr_guess
        endif
        if ( (method_geometry == 2) .or. (method_geometry == 4) ) kz_guess = 1
        first_run = .false.
    endif

    if (present(reset_r_guess)) then
        if ( ( (method_geometry == 2) .or. (method_geometry == 3) ) .and. (reset_r_guess) ) kr_guess = my_first_kr_guess
    endif

    select case (method_geometry)
        case (2)
            call findinterval(kr,kz,r,z,rgrid,zgrid,Nr,Nz,kx_guess=kr_guess,ky_guess=kz_guess)
            kr_guess = kr; kz_guess = kz
            call interpolate2D_decomposed_eq(kr,kz,r,z,initial_electron_density,rgrid,zgrid,table_2D,Nr,Nz)
            return
        case (3)  
            call findinterval(kr,r,rgrid,Nr,k_guess=kr_guess)
            kr_guess = kr;
            call interpolate1D_decomposed_eq(kr,r,initial_electron_density,rgrid,table_1D,Nr)
            return
        case (4)   
            call findinterval(kz,z,zgrid,Nz,k_guess=kz_guess)
            kz_guess = kz;
            call interpolate1D_decomposed_eq(kz,z,initial_electron_density,zgrid,table_1D,Nz)
            return   
    end select
end function

function initial_electron_density_guess(r,z,k_actual,k_first) 
    real(8)                    :: initial_electron_density_guess
    real(8)                    :: r,z 
    integer                    :: k_actual,k_first

    if (k_actual == k_first) then
        initial_electron_density_guess = initial_electron_density(r,z,reset_r_guess=.TRUE.)
    else
        initial_electron_density_tip = initial_electron_density(r,z)
    endif

end function


end module pre_ionised
