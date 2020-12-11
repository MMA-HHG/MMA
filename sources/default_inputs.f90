! It was developed by Jan Vabek

module default_inputs
use write_listing

IMPLICIT NONE
real(8) :: Intensity_entry, Intensity_focus, waist_focus, Curvature_radius_entry, focus_position
character(15)   ::  gas_preset

integer                 :: k1
integer, parameter      :: N_tests = 9
character(*), parameter :: available_tests(N_tests) = (/"test", "test2", "GfP", "GfI", "GfFWHME", "GfFWHMI", "GfH5w", "PI", "PIPPT"/) ! "GfH5w_pre_ionised_PPT"
! integer, parameter      :: test_numbers(N_tests) =  (k1, k1=1,N_tests)

CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!    DEFAULT VALUES FOR NOBLE GASES

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


subroutine preset_parameters_gas
integer         :: switch_ionisation, switch_atom

    ! switch_ionisation
    ! 1 - PPT, 2 - external

    select case (gas_preset)
    case ('Ar_PPT')
        switch_ionisation = 1
        switch_atom = 1
    case ('Ar_ext')
        switch_ionisation = 2
        switch_atom = 1
    case ('Ne_PPT')
        switch_ionisation = 1
        switch_atom = 2
    case ('Ne_ext')
        switch_ionisation = 2
        switch_atom = 2
    case ('Xe_PPT')
        switch_ionisation = 1
        switch_atom = 3
    case ('Xe_ext')
        switch_ionisation = 2
        switch_atom = 3
    case ('Kr_PPT')
        switch_ionisation = 1
        switch_atom = 4
    case ('Kr_ext')
        switch_ionisation = 2
        switch_atom = 4
    case default
        print *, 'wrong preset gas model entry'
    end select


    select case(switch_ionisation)
    case(1)
        switch_rho = 3
    case(2)
        switch_rho = 8
    end select


    select case(switch_atom)
    case (1) ! Argon
        Ui_eV_phys =            15.76D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     3
        n2_phys =               1.d-19 ! Kerr as in n_2*I (cm2/W)


    case (2) ! Neon
        Ui_eV_phys =            21.56D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     7
        n2_phys =               1.d-19 ! Kerr as in n_2*I (cm2/W)

    case (3) ! Xenon
        Ui_eV_phys =            12.13D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     6
        n2_phys =               1.d-19 ! Kerr as in n_2*I (cm2/W)

    case (4) ! Krypton
        Ui_eV_phys =            14.00D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     9
        n2_phys =               1.d-19 ! Kerr as in n_2*I (cm2/W)
    end select



    ! shared default values
    rhont_cm3_phys = 2.7d19! effective density of neutral molecules, it is the density of an ideal gas for 1 bar 0 °C in cm-3 (https://en.wikipedia.org/wiki/Number_density#Units)


    ! generally target-specific values, but kept as these estimates for our results:

    ions_Kerr_ratio = 1.D0/3.D0
    n4_phys = 0.d0 ! Kerr
    switch_dKerr = 1
    
    sigmak_phys = 1.9d-120  ! for only some ionisation models

    tauc_fs_phys = 190.d0 ! THIS IS APPLIED

    alpha_fs_phys = 0.d0
    alphaquad_fscm3_phys = 0.d0 !1.d-3

    NN = 2
    sigman_phys = 0.d0
    rhoabs_cm3_phys = 0.d0

end subroutine preset_parameters_gas


subroutine preset_laser
    super_N = 1
    super_t = 1
    switch_start = 1

    chirp_factor = 0.d0
end subroutine preset_laser


subroutine preset_numerics

    ! field numerical properties (noise, etc.)
    inputfilename_t = 'nrl.dat'
    ! inputfilename_c = '000_000000'
    restartamp = 1.d0
    noise_s = 0.d0
    noise_t = 0.d0
    noise = 0.d0

end subroutine preset_numerics



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!     FUNCTIONS TO RECALCULATE VARIOUS INPUT FORMS 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

real(8) function e_inv2FWHM(e_inv)
    real(8) :: e_inv
    e_inv2FWHM = 2.d0*sqrt(log(2.d0))*e_inv
end function e_inv2FWHM

real(8) function FWHM2e_inv(FWHM)
    real(8) :: FWHM
    FWHM2e_inv = FWHM/2.d0*sqrt(log(2.d0))
end function FWHM2e_inv

real(8) function Convert_pulse_duration(t_in, type_in, type_out, type2_in, type2_out)
    ! 1/e is the reference
    real(8)                 :: t_in, dum
    character(*)            :: type_in, type_out
    character(*), optional  :: type2_in, type2_out

    select case(type_in)
    case('1/e')
        dum = t_in
    case('FWHM')
        dum = t_in/(2.d0*sqrt(log(2.d0)))
    case('rms')
        dum = sqrt(2.d0)*t_in
    case default
        print *, 'wrong input of pulse duration, nothing done'
        Convert_pulse_duration = t_in
        return
    end select

    if (present(type2_in)) then
        print *, type2_in
        select case(type2_in)
        case('Efield')
            dum = dum/sqrt(2.d0)
        case('Intensity')
            ! nothing
        case default
            print *, 'wrong input2 of pulse duration, nothing done'
            Convert_pulse_duration = t_in
            return
        end select
    endif

    select case(type_out)
    case('1/e')
        Convert_pulse_duration = dum
    case('FWHM')
        Convert_pulse_duration = 2.d0*sqrt(log(2.d0))*dum
    case('rms')
        Convert_pulse_duration = dum/sqrt(2.d0)
    case default
        print *, 'wrong output of pulse duration, nothing done'
        Convert_pulse_duration = t_in
    end select

    if (present(type2_out)) then
        print *, type2_out, Convert_pulse_duration 
        select case(type2_out)
        case('Efield')
            Convert_pulse_duration = sqrt(2.d0)*Convert_pulse_duration
        case('Intensity')
            ! nothing
        case default
            print *, 'wrong output2 of pulse duration, nothing done'
            Convert_pulse_duration = t_in
            return
        end select
        print *, type2_out, Convert_pulse_duration 
    endif

end function Convert_pulse_duration

real(8) function ratio_Pin_Pcr_entry2I_entry(Pin_Pcr,wz,n2p,lambda)
    real(8) :: Pin_Pcr,wz,n2p,lambda
    ratio_Pin_Pcr_entry2I_entry = (Pin_Pcr*(lambda/(PI*wz))**2) / n2p
end function ratio_Pin_Pcr_entry2I_entry

real(8) function I_entry2ratio_Pin_Pcr_entry(I_entry,wz,n2p,lambda)
    real(8) :: I_entry,wz,n2p,lambda
    I_entry2ratio_Pin_Pcr_entry = n2p*I_entry*(PI*wz/lambda)**2
end function I_entry2ratio_Pin_Pcr_entry

real(8) function Energy2ratio_Pin_Pcr_entry(Energy,n2p,lambda,t0)
    real(8) :: Energy,n2p,lambda,t0
    Energy2ratio_Pin_Pcr_entry = 2.d0*sqrt(2.d0*PI)*n2p*Energy/(t0*lambda**2)
end function Energy2ratio_Pin_Pcr_entry

real(8) function ratio_Pin_Pcr_entry2Energy(Pin_Pcr,n2p,lambda,t0)
    real(8) :: Pin_Pcr,n2p,lambda,t0
    ratio_Pin_Pcr_entry2Energy = (t0*lambda**2)*Pin_Pcr/(2.d0*sqrt(2.d0*PI)*n2p)
end function ratio_Pin_Pcr_entry2Energy

subroutine Gaussian_focus2Gaussian_entry(I0,w0,z,Iz,wz,Rz,lambda)
    real(8) :: I0,w0,z,Iz,wz,Rz,lambda,zR
    zR = (PI*w0**2)/lambda
    wz = w0 * SQRT(1.D0+(z/zR)**2)
    Iz = I0 * (w0/wz)**2
    Rz = z + (zR**2) / z
end subroutine Gaussian_focus2Gaussian_entry

subroutine Gaussian_entry2Gaussian_focus(Iz,wz,Rz,I0,w0,focus,lambda)
    ! for small z's (large R(z)), this could be critical
    real(8)     :: Iz,wz,Rz,I0,w0,focus,lambda
    w0 = wz**2 / (1.D0+((PI*wz**2)/(lambda*Rz))**2)
    !zR = (PI * lambda * Rz**2 * wz**2) / (lambda**2 * Rz**2 + PI**2 * wz**4)
    focus = -Rz/(1.D0+(lambda*Rz/(PI*wz**2))**2)
    I0 = Iz * (wz/w0)**2
end subroutine Gaussian_entry2Gaussian_focus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!           PRESET VALUES FOR TESTING MODE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function get_test_number(testname)
    character(*)    :: testname

    do k1 = 1, N_tests
        if (testname == available_tests(k1)) then
            get_test_number = k1
            return
        endif
    enddo
    stop "wrong test name"

end function get_test_number

subroutine preset_dispersion_tests(test_number)
    integer :: test_number
    ! Dispersion law
    dispfilename = 'waterchi.tab'

    ! The following values are applied only for Taylorised dispersion law
    !n0 = 1.45d0
    !delta_k_p_fs_per_cm_phys = 0.d0
    !k_pp_fs2_per_cm_phys = -279.d0
    !k_ppp_fs3_per_cm_phys = 1510.d0
    !k_pppp_fs4_per_cm_phys = -4930.d0
    !k_ppppp_fs5_per_cm_phys = 23245.d0

end subroutine preset_dispersion_tests

subroutine preset_delayed_Kerr_tests(test_number)
    integer :: test_number

    ! Delayed Kerr
    
    !xdk = 0.5d0
    !tdk_fs_phys = 77.d0
    !raman_phys = 1.6d-2

end subroutine preset_delayed_Kerr_tests

subroutine preset_numerics_tests(test_number)
    integer :: test_number

    num_proc = 32
    time_limit = 0.48d0

    ! time
    lt = 8.d0
    dim_t = 2048 ! asymmetric
    absorb = 16

    ! space
    lr = 4.d0
    dim_r = 1024
    

    ! propagation & adaptive steps
    delta_z_mm_phys = 1.d-2
    decrease = 2.d-3
    switch_T = 2 ! operator

    ! writing
    rfil_mm_phys = 0.1d0
    rhodist = 100
    outlength_m_phys = 0.001d0

    select case(test_number)
    case(1:6)
        outlength_Efield_m_phys = outlength_m_phys
    case(7:N_tests)
        outlength_Efield_m_phys = 0.0005d0        
    end select
    call save_or_replace(file_id, 'inputs/numerics_physical_output_distance_for_Efield_only', outlength_Efield_m_phys, error, units_in = '[m]')
    
end subroutine preset_numerics_tests

subroutine preset_physics(test_number)
    integer :: test_number

!---------------------------------------------------------------------------------------------------------------------!    
    lambda0_cm_phys = 8.d-5

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1,9)
        gas_preset = 'Ar_PPT'
    case(2:8)
        gas_preset = 'Ar_ext'
    end select
    call save_or_replace(file_id, 'inputs/gas_preset', gas_preset, error, units_in = '[-]')

!---------------------------------------------------------------------------------------------------------------------!
    proplength_m_phys = 0.005d0

!---------------------------------------------------------------------------------------------------------------------!
    w0_cm_phys = 0.1d0
    call save_or_replace(file_id, 'inputs/laser_beamwaist', w0_cm_phys, error, units_in = '[cm]')

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:3,5:N_tests)
        numcrit = 2.0d0
        call save_or_replace(file_id, 'inputs/laser_ratio_pin_pcr', numcrit, error, units_in = '[-]')
    case(4)
        Intensity_entry = 1.d18
        call save_or_replace(file_id, 'inputs/laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    end select
    

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:4,7:N_tests)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
    case(5)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Efield', tp_fs_phys, error, units_in = '[fs]')
    case(6)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Intensity', tp_fs_phys, error, units_in = '[fs]')
    end select   
    
    !call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')


!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1, 2)
        f_cm_phys = 50.d0 ! THIS IS SOMETHING TO COMPUTE
    case(3:N_tests)
        f_cm_phys = 0.d0
    end select

!---------------------------------------------------------------------------------------------------------------------!
    pressure = 1.d0


!---------------------------------------------------------------------------------------------------------------------!
! pre-ionized
    if ( any(test_number == (/8, 9/)) ) then
        call h5gcreate_f(file_id, 'pre_ionised', group_id2, error)
        call save_or_replace(group_id2, 'method_geometry', 1, error, units_in = '[-]')
        call save_or_replace(group_id2, 'method_units', 1, error, units_in = '[-]')
        call save_or_replace(group_id2, 'initial_electrons_ratio', 0.04d0, error, units_in = '[-]')
        call h5gclose_f(group_id2, error)   
    endif    

end subroutine preset_physics

subroutine testing_values(test_number) ! set values for testing
    use HDF5
    use HDF5_helper
    integer :: test_number

    call preset_numerics_tests(test_number)
    call preset_physics(test_number)
end subroutine testing_values

end module default_inputs
