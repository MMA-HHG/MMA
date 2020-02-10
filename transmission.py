import numpy as np
import math
from scipy import special

wir = 100 #in micron waist

I = 2; # in 10^14 W/cm²


z_jet = -58; #in mm
z_fente = 370; #in mm


R_fente = 70; #in micron


q = 29;
alphaq = 3; # in 10^-14 cm²/W
qeff = 4.7;




# /// Def of parameters
wir = wir*1e-6
z_fente = z_fente*1e-3;
R_fente = R_fente*1e-6;
z_jet = z_jet*1e-3;

zrir = wir**2*math.pi/(800e-9)

Rir = z_jet + zrir**2/z_jet





phaseq = alphaq*I;
k0 = 2.*math.pi/(800e-9);


Ratom = q*k0*wir**2*(1.+(z_jet/zrir)**2)**2/(4*phaseq); 

## Calculates the XUV parameters

Ruvx = (1/Rir + 1/Ratom)**(-1.);
kuvx = q*k0;



wuvx = wir*(1.+(z_jet/zrir)**2)**0.5/(np.sqrt(qeff));


#####

z_sol = wuvx**(4)*Ruvx/(wuvx**4 + 4*Ruvx**2/kuvx**2)

zruvx = np.sqrt(Ruvx*z_sol-z_sol**2)

wuvx_fente = np.sqrt(2*((-(z_jet-z_sol) + z_fente)**2 + zruvx**2)/(kuvx*zruvx))

w0uvx = np.sqrt(2*zruvx/(q*k0))

T = 1.-math.exp(-2.0*(R_fente/(wuvx_fente))**2.) 

print(Rir*100,'cm')
print(zrir*100,'cm')
print((z_jet-z_sol)*1000,'mm')
print(w0uvx*1e6,'micron')
print(wuvx_fente*1e6,'micron')
print(T*100,'Transmission in %')








