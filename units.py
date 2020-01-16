import numpy as np
### physical constants
hbar=1.0545718e-34; inverse_alpha_fine=137.035999139; c_light=299792458; elcharge=1.602176565e-19; elmass=9.10938356e-31;
r_Bohr = hbar*inverse_alpha_fine/(c_light*elmass);

# conversion factor to atomic units
TIMEau = (inverse_alpha_fine**2)*hbar/(elmass*c_light**2);
INTENSITYau = (inverse_alpha_fine/(8.0*np.pi))*(hbar**3)/((elmass**2)*(r_Bohr**6));
