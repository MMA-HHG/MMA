module constants
    real(8), parameter  :: PI = 4.d0*acos(1.d0)

    ! physical constants in SI units
    real(8), parameter  :: echarge = 1.60217662d-19
    real(8), parameter  :: emass = 9.10938356d-31
    real(8), parameter  :: c_light = 299792458d0
    real(8), parameter  :: mu_0 = 4.d0*PI*1.d-7
    real(8), parameter  :: eps0 = 1.d0/(mu_0*c_light**2)
    real(8), parameter  :: hbar = 1.054571800d-34
    real(8), parameter  :: alpha_fine_inv = 137.035999139d0
    real(8), parameter  :: alpha_fine = 1.d0/alpha_fine_inv
    real(8), parameter  :: Ip_HeV = 27.21138602d0
end module constants