import numpy as np
### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458.0; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass); alpha_fine = 1.0/inverse_alpha_fine;

# conversion factor to atomic units
TIMEau = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);
INTENSITYau = (inverse_alpha_fine/(8.0*np.pi))*(hbar**3)/((elmass**2)*(r_Bohr**6));
ENERGYau = hbar^2/(elmass*r_Bohr**2);
TIMEau = (elmass*r_Bohr**2)/hbar;
