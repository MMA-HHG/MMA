! It was developed by Jan Vabek

module default_inputs
use write_listing

IMPLICIT NONE
real(8) :: Intensity_entry, Intensity_focus, waist_focus, Curvature_radius_entry, invCurvature_radius_entry, focus_position
character(15)   ::  gas_preset

integer                 :: k1
integer, parameter      :: N_tests = 119
!character, parameter :: available_tests(N_tests) = (/ "test", "test2", "GfP", "GfI", "GfFWHME", "GfFWHMI", "GfH5w", &
!                                                        "PI", "PIPPT", "pressure", "ELI1", "ELI1ppt", "ELI2", "ELI3", &
!                                                        "ELI4", "ELI_PI_PPT_Kr", "Ar_vacuum1", "Ar_vacuum1_long", &
!                                                        "Ar_vacuum2_foc_entry", "Ar_vacuum2_foc_half", "Ar_vacuum2_foc_end", &
!                                                        "Ar_vacuum2_f_half", "TDSE1", "TDSE1_long", "TDSE1_long2", &
!                                                        "TDSE1_10mm", "TDSE1_10mm_one_node", "TDSE1_10mm_sparse", &
!                                                        "TDSE1_1mm", "TDSE1_03mm", "run_test" /) ! "GfH5w_pre_ionised_PPT"
character(len=255), parameter :: available_tests(N_tests) = [character(len=255) :: "test", "test2", "GfP", "GfI", "GfFWHME",& 
                                                        "GfFWHMI", "GfH5w", "PI", "PIPPT", "pressure", "ELI1", "ELI1ppt", &
                                                        "ELI2", "ELI3", "ELI4", "ELI_PI_PPT_Kr", "Ar_vacuum1", "Ar_vacuum1_long", &
                                                        "Ar_vacuum2_foc_entry", "Ar_vacuum2_foc_half", "Ar_vacuum2_foc_end", &
                                                        "Ar_vacuum2_f_half", "TDSE1", "TDSE1_long", "TDSE1_long2", &
                                                        "TDSE1_10mm", "TDSE1_10mm_one_node", "TDSE1_10mm_sparse", &
                                                        "TDSE1_1mm", "TDSE1_03mm", "run_test", "test_modulations_base", &
                                                        "test_modulations_base_half_dens", "test_modulations_base_half_dens_table", &
                                                        "test_modulations_base2", "test_modulations_base_half_dens2", &
                                                        "test_modulations_base_half_dens_table2", "test_modulations_base_100dens", &
                                                        "test_modulations_base_100dens_table", "test_modulations_base_100dens2", &
                                                        "test_modulations_base_100dens_table2", "test_modulations_base_100dens_table2_short", &
                                                        "test_modulations_base_100dens2_short", &
                                                        ("undefined", k1 = 44, 100), &
                                                        "test1_modT1", "test1_mod_incT1", "test1_mod_decT1", &
                                                        "test1_modT2", "test1_mod_incT2", "test1_mod_decT2", &
                                                        "test1_modT3", "test1_mod_incT3", "test1_mod_decT3", &
                                                        "undefined", &
                                                        "test2_modT1", "test2_mod_incT1", "test2_mod_decT1", &
                                                        "test2_modT2", "test2_mod_incT2", "test2_mod_decT2", &
                                                        "test2_modT3", "test2_mod_incT3", "test2_mod_decT3"] ! "GfH5w_pre_ionised_PPT"
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

    ! ionization potentials https://en.wikipedia.org/wiki/Ionization_energies_of_the_elements_(data_page)

    !----------------------------------------------------------------------------------------------------------------------
    ! Kerr effect:
    ! There are various sources of the Kerr effect constants, some of them are overviewed in 

    ! [*]: C. Bree, A. Demircan and G. Steinmeyer, "Method for Computing the Nonlinear Refractive Index via Keldysh Theory," 
    ! in IEEE Journal of Quantum Electronics, vol. 46, no. 4, pp. 433-437, April 2010, doi: 10.1109/JQE.2009.2031599.

    ! The following table contains optimised values for 800-nm field in 1e-19 cm2/W
    
    !     Keldysh | ref [8]   | ESHG  |
    ! He    0.034 |   0.041   | 0.037 |
    ! Ne    0.085 |   0.074   | 0.094 |
    ! Ar    0.796 |   1.04    | 1.09  |
    ! Kr    1.89  |   2.94    | 2.47  |
    ! Xe    5.48  |   9.35    | 6.39  |

    ! ref [8] is H.J.Lehmeier, W.Leupacher and A.Penzkofer, https://doi.org/10.1016/0030-4018(85)90069-0
    ! ESHG is compiled from various sorces referred in [*]
    ! [*] also provide the formula (12) for scaling for different wave-lengths

    ! We use ref [8] values as pre-set values
    !----------------------------------------------------------------------------------------------------------------------

    

    case (1) ! Argon
        Ui_eV_phys =            15.75962D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     3
        n2_phys =               1.04d-19 ! Kerr as in n_2*I (cm2/W)


    case (2) ! Neon
        Ui_eV_phys =            21.5646D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     7
        n2_phys =               0.074d-19 ! Kerr as in n_2*I (cm2/W)

    case (3) ! Xenon
        Ui_eV_phys =            12.1298D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     6
        n2_phys =               9.35d-19 ! Kerr as in n_2*I (cm2/W)

    case (4) ! Krypton
        Ui_eV_phys =            13.99961D0 ! ionisation potential (eV)
        angular_momentum  =     1 ! (-)
        residue_charge =        1.0D0 ! (-)

        switch_dispersion =     9
        n2_phys =               2.94d-19 ! Kerr as in n_2*I (cm2/W)
    end select

    ! Helium 24.58738D0


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

    NN = 0
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

    print *, 'conversion to pin/pcr'
    print *, 'n2p', n2p
    print *, 'I_entry', I_entry
    print *, 'wz', wz
    print *, 'lambda', lambda

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

subroutine Gaussian_focus2Gaussian_entry(I0,w0,z,Iz,wz,invRz,lambda)
    real(8) :: I0,w0,z,Iz,wz,invRz,lambda,zR
    zR = (PI*w0**2)/lambda
    wz = w0 * SQRT(1.D0+(z/zR)**2)
    Iz = I0 * (w0/wz)**2
    invRz = z/(z**2 + zR**2)
end subroutine Gaussian_focus2Gaussian_entry

subroutine Gaussian_entry2Gaussian_focus(Iz,wz,invRz,I0,w0,focus,lambda)
    ! for small z's (large R(z)), this could be critical
    real(8)     :: Iz,wz,invRz,I0,w0,focus,lambda
    w0 = sqrt( wz**2 / (1.D0+((PI*invRz*wz**2)/lambda)**2) )
    !zR = (PI * lambda * Rz**2 * wz**2) / (lambda**2 * Rz**2 + PI**2 * wz**4)
    ! focus = -Rz/(1.D0+(lambda*Rz/(PI*wz**2))**2)
    focus = -invRz/(invRz**2 + (lambda/(PI*wz**2))**2)
    I0 = Iz * (wz/w0)**2
end subroutine Gaussian_entry2Gaussian_focus




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!           PRESET VALUES FOR TESTING MODE 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

integer function get_test_number(testname)
    character(*)    :: testname

    do k1 = 1, N_tests
        if (testname.EQ.available_tests(k1)) then
            if (testname.eq.'undefined') stop "udenfined test"
            print *, "TEST:", testname
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


    ! computer time & procs
    select case(test_number)
    case(1:22)
        num_proc = 32
        time_limit = 0.48d0
    case(23,24)
        num_proc = 32
        time_limit = 1.98d0 
    case(25,26,28:30, 38:43)  
        num_proc = 32
        time_limit = 23.95d0 
    case(27)  
        num_proc = 16
        time_limit = 23.95d0 
    case(31:37)
        num_proc = 4
        time_limit = 23.95d0

    case(101:109, 111:119)  
        num_proc = 32
        time_limit = 19.9d0 
    end select

    ! time
    select case(test_number)
    case(1:14,17,18:30,32:43)
        lt = 8.d0
        dim_t = 2048 ! asymmetric
        absorb = 16
    case(16)
        lt = 12.d0
        dim_t = 2048 ! asymmetric
        absorb = 16     
    case(31)
        lt = 12.d0
        !dim_t = 2048 ! asymmetric
        dim_t = 16 ! asymmetric
        absorb = 16   

    case(101:109, 111:119)
        lt = 8.d0
        dim_t = 2048 ! asymmetric
        absorb = 16  
    end select


    ! space
    lr = 4.d0
    dim_r = 1024
    !dim_r = 8
    

    ! propagation & adaptive steps
    decrease = 2.d-3

    select case(test_number)
    case(1:43)
        switch_T = 2 ! operator

    case(101:103, 111:113)
        switch_T = 1
    case(104:106, 114:116)
        switch_T = 2
    case(107:109, 117:119)
        switch_T = 3
    end select

    select case(test_number)
    case(1:22)
        delta_z_mm_phys = 1.d-2  
    case(23:43)
        delta_z_mm_phys = 1.d-3    

    case(101:109,111:119)
        delta_z_mm_phys = 1.d-3   
    end select

    ! writing
    rfil_mm_phys = 0.1d0
    rhodist = 100

    select case(test_number)
    case(1:6)
        outlength_m_phys = 0.001d0
        outlength_Efield_m_phys = outlength_m_phys
    case(7:14)
        outlength_m_phys = 0.001d0
        outlength_Efield_m_phys = 0.0005d0 
    case(15,16)
        outlength_m_phys = 0.000125d0
        outlength_Efield_m_phys = 0.075d0    
    case(17)
        outlength_m_phys = 0.0005d0
        outlength_Efield_m_phys = 0.075d0
    case(18:22)
        outlength_m_phys = 0.001d0
        outlength_Efield_m_phys = 0.075d0     
    case(23)
        outlength_m_phys = 0.000001d0  
        outlength_Efield_m_phys = 0.075d0    
    case(24:27,29,30,31)
        outlength_m_phys = 0.000002d0  
        outlength_Efield_m_phys = 0.075d0 
    case(28)
        outlength_m_phys = 0.0005d0  
        outlength_Efield_m_phys = 0.075d0    
    !case(31)
    !    outlength_m_phys = 0.000005d0  
    !    outlength_Efield_m_phys = 0.00075d0    
    case(32:41)
        outlength_m_phys = 1.d-4  
        outlength_Efield_m_phys = 0.075d0    
    case(42,43)
        outlength_m_phys = 1.d-5  
        outlength_Efield_m_phys = 0.075d0  
        
    case(101:109)
        outlength_m_phys = 1.d-4  
        outlength_Efield_m_phys = 0.075d0
    case(111:119)
        outlength_m_phys = 0.5d-4  
        outlength_Efield_m_phys = 0.075d0

    end select
    call save_or_replace(file_id, 'inputs/numerics_physical_output_distance_for_Efield_only', &
                         outlength_Efield_m_phys, error, units_in = '[m]')
    
end subroutine preset_numerics_tests

subroutine preset_physics(test_number)
    integer :: test_number

!---------------------------------------------------------------------------------------------------------------------!    
    select case(test_number)
    case(1:10,17,18:22,32:43)
        lambda0_cm_phys = 8.d-5
    case(11:16,23:30)
        lambda0_cm_phys = 7.92d-5
    case(31)
        lambda0_cm_phys = 7.92d-5

    case(101:109,111:119)
        lambda0_cm_phys = 8.d-5
    end select

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1,9,12,15,17,18:22,32:43)
        gas_preset = 'Ar_PPT'
    case(2:8,10,11,13,14)
        gas_preset = 'Ar_ext'
    case(16,23:30)
        gas_preset = 'Kr_PPT'
    case(31)
        gas_preset = 'Ar_PPT'

    case(101:109,111:119)
        gas_preset = 'Ar_PPT'
    end select
    call save_or_replace(file_id, 'inputs/gas_preset', gas_preset, error, units_in = '[-]')

!---------------------------------------------------------------------------------------------------------------------!
 
     select case(test_number)
    case(1:11,13,14)
        proplength_m_phys = 0.005d0
    case(12)
        proplength_m_phys = 0.0025d0
    case(16,17)
        proplength_m_phys = 0.015d0
    case(18)
        proplength_m_phys = 0.06d0
    case(19:22)
        proplength_m_phys = 0.015d0
    case(23)
        proplength_m_phys = 0.0002d0
    case(24,25)
        proplength_m_phys = 0.002d0
    case(26,27,28,32:41)
        proplength_m_phys = 0.01d0
    case(42,43)
        proplength_m_phys = 2.5d-4
    case(29)
        proplength_m_phys = 0.001d0
    case(30)
        proplength_m_phys = 0.0003d0
    case(31)
        proplength_m_phys = 0.0003d0

    case(101:109)
        proplength_m_phys = 0.01d0
    case(111:119)
        proplength_m_phys = 0.005d0
    end select   

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:10)
        w0_m_phys = 0.001d0      ! m
    case(11:16)
        w0_m_phys = 0.00011d0    ! m
    case(17,18,22)
        w0_m_phys = 0.0001d0     ! m
    case(19:21)
        waist_focus = 100.d-6   ! m
    case(23:30,32:43)
        waist_focus = 110.d-6   ! m
    case(31)
        w0_m_phys = 0.001d0      ! m
        !waist_focus = 110.d-6   ! m

    case(101:109,111:119)
        waist_focus = 110.d-6   ! m
    end select

    select case(test_number)
    case(1:18,22)    
        call save_or_replace(file_id, 'inputs/laser_beamwaist_entry', w0_m_phys, error, units_in = '[m]')
    case(31)    
        call save_or_replace(file_id, 'inputs/laser_beamwaist_entry', w0_m_phys, error, units_in = '[m]')
    case(19:21,23:30,32:43)
        call save_or_replace(file_id, 'inputs/laser_focus_beamwaist_Gaussian', waist_focus, error, units_in = '[m]')

    case(101:109,111:119)
        call save_or_replace(file_id, 'inputs/laser_focus_beamwaist_Gaussian', waist_focus, error, units_in = '[m]')
    end select

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:3,5:10)
        numcrit = 2.0d0
        call save_or_replace(file_id, 'inputs/laser_ratio_pin_pcr', numcrit, error, units_in = '[-]')
    case(4,15,16,17,18,22)
        Intensity_entry = 1.d18
        call save_or_replace(file_id, 'inputs/laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    case(19:21)
        Intensity_focus = 1.d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    case(31)
        Intensity_focus = 1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    case(11,12)
        Intensity_entry = 1.129755554227896d19
        call save_or_replace(file_id, 'inputs/laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    case(13,14)
        Intensity_entry = 1.5d18
        call save_or_replace(file_id, 'inputs/laser_intensity_entry', Intensity_entry, error, units_in = '[SI]')
    case(23:30,32:34)
        Intensity_focus = 1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    case(35:39)
        Intensity_focus = 2.d0*1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    case(40:43)
        Intensity_focus = 5.d0*1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')

    case(101:109)
        Intensity_focus = 2.d0*1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    case(111:119)
        Intensity_focus = 5.d0*1.8d18
        call save_or_replace(file_id, 'inputs/laser_focus_intensity_Gaussian', Intensity_focus, error, units_in = '[SI]')
    end select
    

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:4,7:10,17,18:30,32:43)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
    case(5)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Efield', tp_fs_phys, error, units_in = '[fs]')
    case(6)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Intensity', tp_fs_phys, error, units_in = '[fs]')
    case(31)
        tp_fs_phys = 10.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Intensity', tp_fs_phys, error, units_in = '[fs]')
    case(11:16)
        tp_fs_phys = 35.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_FWHM_Intensity', tp_fs_phys, error, units_in = '[fs]')

    case(101:109,111:119)
        tp_fs_phys = 50.d0
        call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')
    end select   
    
    !call save_or_replace(file_id, 'inputs/laser_pulse_duration_in_1_e_Efield', tp_fs_phys, error, units_in = '[fs]')


!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1, 2)
        f_cm_phys = 50.d0 ! THIS IS SOMETHING TO COMPUTE
    case(22)
        f_cm_phys = 0.75d0
    case(3:18)
        f_cm_phys = 0.d0
    case(31)
        focus_position = 0.0d0
    case(19,23:30,32:43)
        focus_position = 0.0d0
    case(20)
        focus_position = 0.0075d0
    case(21)
        focus_position = 0.015d0

    case(101:109,111:119)
        focus_position = 0.0d0
    end select

!---------------------------------------------------------------------------------------------------------------------!
    select case(test_number)
    case(1:9)
        pressure = 1.d0
    case(10)
        pressure = 0.5d0
    case(11:14)
        pressure = 0.035d0
    case(15,16)
        pressure = 0.015d0
    case(17, 18, 19:22)
        pressure = 0.001d0
    case(23:30)
        pressure = 0.025d0
    case(31)
        pressure = 0.025d0
    case(32,34,35,37)
        pressure = 0.05d0
    case(33,36)
        pressure = 0.5d0*0.05d0
    case(38,40,43)
        pressure = 0.05d0
    case(39,41:42)
        pressure = 1.0d2 * 0.05d0

    case(101,104,107,111,114,117)
        pressure = 0.05d0
    case(102,105,108,112,115,118)
        pressure = 1.0d2 * 0.05d0
    case(103,106,109,113,116,119)
        pressure = 1.0d-2 * 0.05d0
    end select
    


!---------------------------------------------------------------------------------------------------------------------!
! pre-ionized
    if ( any(test_number == (/8, 9, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23, 24, 25, 26, 27, 28, 29, 30/)) ) then
        call h5gcreate_f(file_id, 'pre_ionised', group_id2, error)
        call save_or_replace(group_id2, 'method_geometry', 1, error, units_in = '[-]')
        call save_or_replace(group_id2, 'method_units', 1, error, units_in = '[-]')
        if ( any(test_number == (/8, 9, 14, 16/)) ) then
            call save_or_replace(group_id2, 'initial_electrons_ratio', 0.04d0, error, units_in = '[-]')
        else
            call save_or_replace(group_id2, 'initial_electrons_ratio', 0.0d0, error, units_in = '[-]')
        endif
        call h5gclose_f(group_id2, error)   
    endif    


    !---------------------------------------------------------------------------------------------------------------------!
    ! density modulation
    if ( any(test_number == (/34,37,39,41,42, &
                              102,103,105,106,108,109, &
                              112,113,115,116,118,119/)) ) then
        call h5gcreate_f(file_id, 'density_mod', group_id2, error)

        call create_dset(group_id2, 'zgrid', (/ 0.d0 , proplength_m_phys /), 2)
        call h5_add_units_1D(group_id2, 'zgrid', '[m]')

        if ( any(test_number == (/34,37/)) ) call create_dset(group_id2, 'table', (/ 0.5d0 , 0.5d0 /), 2)

        if ( any(test_number == (/39,41,42/)) ) call create_dset(group_id2, 'table', (/ 1.0d-2 , 1.0d-2 /), 2)

        if ( any(test_number == (/102,105,108,112,115,118/)) ) call create_dset(group_id2, 'table', (/ 1.0d-2 , 1.0d-2 /), 2)
        if ( any(test_number == (/103,106,109,113,116,119/)) ) call create_dset(group_id2, 'table', (/ 1.0d2  , 1.0d2 /), 2)

        call h5_add_units_1D(group_id2, 'table', '[-]')

        call h5gclose_f(group_id2, error)   
    endif    

end subroutine preset_physics

subroutine testing_values(test_number) ! set values for testing
    use HDF5
    use HDF5_helper_serial
    integer :: test_number

    call preset_numerics_tests(test_number)
    call preset_physics(test_number)
end subroutine testing_values

end module default_inputs
