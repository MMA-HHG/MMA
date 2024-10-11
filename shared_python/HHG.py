# import numpy as np
import units

# Some HHG characteristics

def ComputeCutoff(Intensity, omega, Ip):
  '''
  It computes the cutoff from the formula E_cutoff = (3.17Up+Ip) Energy [a.u.] and harmonic order
  inputs are intensity [a.u.], omega [a.u.] and Ip [a.u.]
  It returns [Energy, order]
  '''

  Up = Intensity/(4.0*omega**2)
  Energy = 3.17*Up + Ip
  return Energy, Energy/omega

def ComputeInvCutoff(order,omega,Ip):
  '''
  order [-]
  omega [a.u.]
  Ip [a.u.]
  returns intensity in a.u.
  
  The related formula is 'E_cutoff = (3.17Up+Ip)'
  '''
  
  return (4.0*omega**2) * (omega*order - Ip)/3.17

# Ionisation potentials in atomic units
Ip_list={'H' : 0.5,
         'He': 0.9036,
         'Ne': 0.7924,
         'Ar': 0.5792,
         'Kr': 0.5145,
         'Xe': 0.4458}

def ComputeCutoff_gas(Intensity, omega, gas = 'H'):
  '''
  It computes the cutoff from the formula E_cutoff = (3.17Up+Ip) Energy [a.u.] and harmonic order
  inputs are intensity [a.u.], omega [a.u.] and gas (available 'H', 'He', 'Ne', 'Ar', 'Kr', 'Xe')
  It returns [Energy, order]
  '''

  Up = Intensity/(4.0*omega**2)
  Energy = 3.17*Up + Ip_list[gas]
  return Energy, Energy/omega


def ComputeInvCutoff_gas(order,omega,gas = 'H'):
  '''
  order [-]
  omega [a.u.]
  gas (available 'H', 'He', 'Ne', 'Ar', 'Kr', 'Xe')
  returns intensity in a.u.
  
  The related formula is 'E_cutoff = (3.17Up+Ip)'
  '''
  
  return (4.0*omega**2) * (omega*order - Ip_list[gas])/3.17


def Critical_ATI_intensity_rough(Ip):
    '''
    It returns a very rough estimate of the above-threshold-ionisation intensity [a.u.]
    from the ionisation potential Ip [a.u.]. It assumes the suppression of the Coulomb
    potential [-1/x] without any further corrections.
    '''
    return (Ip**2/(4.))**2


def eta_opt(omega_driver_SI, XUV_polarisability, IR_polarisability):
    """
    Compute the optimal ionisation degree for phase matching.
    (It compensates the dispersion of both IR and XUV.)

    Parameters
    ----------
    xxx

    Returns
    -------
    optimal ionisation degree [-]
    """
        
    delta_polarisability = IR_polarisability - XUV_polarisability

    return (omega_driver_SI**2) *units.eps0*units.elmass*delta_polarisability/(units.elcharge**2)
    