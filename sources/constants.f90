module constants
    real(8), parameter  :: PI = 4.d0*atan(1.d0)

    ! physical constants in SI units
    real(8), parameter  :: echarge = 1.60217662d-19
    real(8), parameter  :: emass = 9.10938356d-31
    real(8), parameter  :: c_light = 299792458.d0
    real(8), parameter  :: mu_0 = 4.d0*PI*1.d-7
    real(8), parameter  :: eps0 = 1.d0/(mu_0*c_light**2)
    real(8), parameter  :: hbar = 1.054571800d-34
    real(8), parameter  :: alpha_fine_inv = 137.035999139d0
    real(8), parameter  :: alpha_fine = 1.d0/alpha_fine_inv
    real(8), parameter  :: r_Bohr = hbar*alpha_fine_inv/(c_light*emass)
    real(8), parameter  :: Ip_HeV = 27.21138602d0

    real(8), parameter  :: TIMEau = (emass*r_Bohr**2) / hbar

    CONTAINS

    function ConvertPhoton(x,inp,outp)
    ! It is not optimalised (case constructs etc.). An alternative for a fast use: call once with x=1. to get conversion numerical factor.
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
                omega = x * units.TIMEau
            case ('eV')
                omega = x * np.elcharge/(emass*alpha_fine**2*c_light**2)
            case ('T0SI')
                omega = units.TIMEau*2.d0*PI/x
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
                ConvertPhoton = omega/units.TIMEau
            case ('eV')
                ConvertPhoton = omega/(units.elcharge/(emass*alpha_fine**2 * c_light**2))
            case ('T0SI')
                ConvertPhoton = units.TIMEau*2.d0*PI/omega
            case ('T0au')
                ConvertPhoton = 2.d0*PI/omega
            case ('Joule')
                ConvertPhoton = omega*(emass*alpha_fine**2 * c_light**2)
            case default
                print *, 'Wrong output in ConvertPhoton'
        end select

    end function
end module constants
