!> @brief This module contains phyisical cnostants and the function \ref
!! ConvertPhoton that links various electromagnetic-field descriptions
!! (wavelength, frequency, photon's energy, ...).
!!
!! @author Jan Vábek
module constants
    real(8), parameter  :: PI = 4.d0*atan(1.d0)

    ! physical constants in SI units
    real(8), parameter  :: echarge = 1.60217662d-19                         !< electron's charge [SI]
    real(8), parameter  :: emass = 9.10938356d-31                           !< electron's mass [SI]
    real(8), parameter  :: c_light = 299792458.d0                           !< the speed of light in vacuum [SI]
    real(8), parameter  :: mu_0 = 4.d0*PI*1.d-7                             !< vacuum permeability [SI]
    real(8), parameter  :: eps0 = 1.d0/(mu_0*c_light**2)                    !< vacuum permittivity [SI]
    real(8), parameter  :: hbar = 1.054571800d-34                           !< reduced Planck's constant [SI]
    real(8), parameter  :: alpha_fine_inv = 137.035999139d0                 !< inverse fine-structure constant (\f$1/\alpha_{\text{fine}}\f$) 
    real(8), parameter  :: alpha_fine = 1.d0/alpha_fine_inv                 !< fine-structure constant
    real(8), parameter  :: r_Bohr = hbar*alpha_fine_inv/(c_light*emass)     !< Bohr's radius [SI]
    real(8), parameter  :: Ip_HeV = 27.21138602d0                           !< The Hydrogen's ionisation potential [eV]

    real(8), parameter  :: TIMEau = (emass*r_Bohr**2) / hbar                !< one atomic unit of time (time[SI]=time[a.u.]*TIMEau)

    CONTAINS

    !! An alternative for a fast use: call once with x=1. to get conversion numerical factor.
    !> @brief converts various electromagnetic-field descriptions (wavelength, frequency, photon's energy, ... ).
    !! Available inputs/outputs are: \n
    !! * `omegaau`  - frequency/energy in atomic units
    !! * `omegaSI`  - (angular) frequency in SI
    !! * `lambdaSI` - wavelength in SI
    !! * `lambdaau` - wavelength in atomic units
    !! * `T0SI`     - period in SI
    !! * `T0au`     - period in atomic units
    !! * `eV`       - energy in eV
    !! * `Joule`    - period in atomic units
    !!
    !! Note: If needed to use for an array/in a loop, consider calling once with x=1 to get conversion numerical factor.
    !!
    !! @param[in]       x           input
    !! @param[in]       inp         specifier of the input
    !! @param[in]       outp        specifier of the output
    !! @return                      value specified by **outp**
    function ConvertPhoton(x,inp,outp)
    ! It is not optimalised (case constructs etc.). 
        real(8)         :: ConvertPhoton

        real(8)         :: x, omega
        character(*)    :: inp, outp

        select case (inp)
            case ('omegaau')
                omega = x
            case ('lambdaSI')
                omega = 2.d0 * PI* hbar / (x * emass * c_light * alpha_fine**2)
            case ('lambdaau')
                omega = 2.d0 * PI/(alpha_fine*x)
            case ('omegaSI')
                omega = x * TIMEau
            case ('eV')
                omega = x * echarge/(emass*alpha_fine**2*c_light**2)
            case ('T0SI')
                omega = TIMEau*2.d0*PI/x
            case ('T0au')
                omega = 2.d0*PI/x
            case ('Joule')
                omega = x / (emass*alpha_fine**2 * c_light**2)
            case default
                print *, 'Wrong input in ConvertPhoton'
        end select

        select case (outp)
            case ('omegaau')
                ConvertPhoton = omega
            case ('lambdaSI')
                ConvertPhoton = 2.d0*PI*hbar/(omega*emass*c_light*alpha_fine**2)
            case ('lambdaau')
                ConvertPhoton = 2.d0*PI/(alpha_fine*omega)
            case ('omegaSI')
                ConvertPhoton = omega/TIMEau
            case ('eV')
                ConvertPhoton = omega/(echarge/(emass*alpha_fine**2 * c_light**2))
            case ('T0SI')
                ConvertPhoton = TIMEau*2.d0*PI/omega
            case ('T0au')
                ConvertPhoton = 2.d0*PI/omega
            case ('Joule')
                ConvertPhoton = omega*(emass*alpha_fine**2 * c_light**2)
            case default
                print *, 'Wrong output in ConvertPhoton'
        end select

    end function
end module constants
