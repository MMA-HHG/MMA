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

implicit none
private
public  :: a

INTERFACE a
    procedure a1
END INTERFACE



CONTAINS



subroutine preset_parameters
integer         :: switch_ionisation, switch_atom


    switch_ionisation = 1; ! PPT
    switch_ionisation = 2; ! External

    switch_atom = 1;

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


end subroutine preset_parameters

end module default_inputs
