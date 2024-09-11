import numpy as np
### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458.0; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass); alpha_fine = 1.0/inverse_alpha_fine;
Ip_HeV = 27.21138602;
mu0 = 4.0*np.pi*1e-7; eps0 = 1.0/(mu0*c_light**2);
r_electron_classical = r_Bohr*(alpha_fine**2)

Boltzmann_constant = 1.380649e-23
Avogadro_constant = 6.02214076e23
universal_gas_constant = Boltzmann_constant*Avogadro_constant

N_ref_default = 1e5/(Boltzmann_constant*(273.15+20.)) # reference gas number density (p = 1 bar & T = 20 °C)

# conversion factor to atomic units
TIMEau = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);
INTENSITYau = (inverse_alpha_fine/(8.0*np.pi))*(hbar**3)/((elmass**2)*(r_Bohr**6));
ENERGYau = hbar**2/(elmass*r_Bohr**2);
TIMEau = (elmass*r_Bohr**2)/hbar;
EFIELDau = hbar*hbar/(elmass*r_Bohr*r_Bohr*r_Bohr*elcharge)
LENGTHau = r_Bohr