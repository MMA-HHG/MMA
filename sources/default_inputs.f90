! It was developed by Jan Vabek

! module tabulated_parameters
! implicit none


! end type atom_info

! ! Argon
! type Argon_table
!     real(8), parameter      :: ionization_potential_eV = 15.8D0
!     integer
! end

! end module tabulated_parameters

module default_inputs
use write_listing

implicit none
! private
! public  :: a

! INTERFACE a
!     procedure a1
! END INTERFACE



CONTAINS



subroutine preset_parameters
integer         :: switch_ionisation, switch_atom


    switch_ionisation = 1; ! PPT
    !switch_ionisation = 2; ! External

    switch_atom = 1;


    ! shared default values
    rhont_cm3_phys = 2.7d19! effective density of neutral molecules, it is the density of an ideal gas for 1 bar 0 °C in cm-3 (https://en.wikipedia.org/wiki/Number_density#Units)

    select case(switch_ionisation)
    case(1)
        switch_rho = 3
    case(2)
        switch_rho = 8
    ! case default
    !     print *, 'unsupported'
    !     stop
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

    ! Preset values for models that should be partially optional (see 'argon.inp' for details)

    ! field numerical properties (noise, etc.)
    inputfilename_t = 'nrl.dat'
    inputfilename_c = '000_000000'
    restartamp = 1.d0
    noise_s = 0.d0
    noise_t = 0.d0
    noise = 0.d0
    f_cm_phys = 50.d0 ! THIS IS SOMETHING TO COMPUTE
    chirp_factor = 0.d0



    ! Dispersion law
    dispfilename = 'waterchi.tab'

    ! The following values are applied only for Taylorised dispersion law
    n0 = 1.45d0
    delta_k_p_fs_per_cm_phys = 0.d0
    k_pp_fs2_per_cm_phys = -279.d0
    k_ppp_fs3_per_cm_phys = 1510.d0
    k_pppp_fs4_per_cm_phys = -4930.d0
    k_ppppp_fs5_per_cm_phys = 23245.d0

    ! Dealyed Kerr + chi5
    switch_dKerr = 1
    xdk = 0.5d0
    tdk_fs_phys = 77.d0
    raman_phys = 1.6d-2
    n4_phys = 0.d0

    ! variables for only some ionisation models
    sigmak_phys = 1.9d-120

    reduced_mass = 0.5d0
    KKp = 3
    sigmakp_phys = 8.6d-27
    rhoslg1_phys = 2.d17
    sigma_phys = 2.d0
    sigmacv_ref_phys = 4.35d-7
    I_ref_phys = 42.32d12
    exp_ref = -3.3d0
    KKpp = 2
    sigmakpp_phys = 1.3d-12
    rhosat_phys = 2.d17
    rhont_N2_cm3_phys = 2.2d19
    Ui_N2_eV_phys = 15.6d0
    angular_momentum_N2 = 0
    residue_charge_N2 = 0.9d0
    T_init_eV_phys = 0.025d0

    tauc_fs_phys = 190.d0 ! THIS IS APPLIED
    alpha_fs_phys = 0.d0
    alpha1_fs_phys = 3.3d-3
    alphah_fs_phys = 1.d-3
    alphaquad_fscm3_phys = 1.d-3

    NN = 2
    sigman_phys = 0.d0
    rhoabs_cm3_phys = 0.d0

end subroutine preset_parameters

end module default_inputs
