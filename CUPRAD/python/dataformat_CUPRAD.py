"""Utilities for loading and analysing CUPRAD output archives.

The module provides a lightweight data container for CUPRAD simulation results,
including propagated electric fields, plasma profiles, fluence, spectra,
ionisation tables, and selected simulation parameters useful for plotting and
post-processing. It also provides some methods that can be applied to perform
usual analyses of the data.
"""
import numpy as np
import units
import mynumerics as mn
import MMA_administration as MMA

# def load_data

class empty_class:
    pass


class get_data:
    """Container and post-processing interface for CUPRAD simulation data.

    The class loads CUPRAD field outputs, grids, medium parameters, laser
    parameters, and selected diagnostics from an HDF5 archive. The main electric
    field is stored as ``E_zrt``, using the ordering ``z, r, t``.

    Parameters
    ----------
    InputArchive : h5py-like file object
        Open archive containing CUPRAD inputs, outputs, logs, and diagnostics.
    r_resolution : list, optional
        Radial loading specification. ``[True]`` loads the full radial grid.
        Otherwise, the expected form is ``[False, dr, rmax]``, where ``dr`` is
        the requested radial spacing and ``rmax`` is the maximal radius to read.

    Attributes
    ----------
    E_zrt : ndarray
        Electric field array ordered as ``z, r, t``.
    tgrid, rgrid, zgrid : ndarray
        Temporal, radial, and propagation grids.
    omega0 : float
        Central angular frequency of the driving laser in SI units.
    k0_wave : float
        Central vacuum wave number of the driving laser.
    VG_IR : float
        Infrared group velocity inferred from the CUPRAD log.

    Available methods
    -----------------
    vacuum_shift
        Shift the electric field to coordinates moving at the vacuum speed of
        light.
    complexify_envel
        Construct the complex envelope of the electric field by removing the
        carrier oscillation.
    get_Fluence
        Load or compute the radial fluence diagnostic.
    get_plasma
        Load the plasma-density output.
    compute_spectrum
        Fourier-transform the temporal electric field and optionally compute a
        radially integrated spectrum.
    get_ionisation_model
        Load the tabulated ionisation model.
    co_moving_t_grid
        Return a time grid shifted to the co-moving frame. This helper is
        attached during initialisation.

    Notes
    -----
    The class assumes the archive follows the MMA/CUPRAD HDF5 layout defined by
    ``MMA.paths``.
    """
    def __init__(self,InputArchive,r_resolution=[True]):
        """Initialise the CUPRAD data object from an archive.

        The constructor loads grids, electric fields, laser parameters, medium
        parameters, optional density modulation data, and selected long-step CUPRAD
        diagnostics. If requested, the radial grid and field data are down-sampled while
        loading.

        Parameters
        ----------
        InputArchive : h5py-like file object
            Open archive containing the CUPRAD simulation data.
        r_resolution : list, optional
            Radial resolution selector. ``[True]`` keeps the full radial resolution;
            otherwise use ``[False, dr, rmax]`` to load approximately every ``dr`` up to
            radius ``rmax``.
        """
        full_resolution = (r_resolution[0] is True)
        self.omega0 = mn.ConvertPhoton(1e-2*mn.readscalardataset(InputArchive,
                      MMA.paths['CUPRAD_inputs'] +'/laser_wavelength','N'),'lambdaSI','omegaSI')
        self.k0_wave = 2.0*np.pi/mn.ConvertPhoton(self.omega0,'omegaSI','lambdaSI')
        self.tgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/tgrid'][:]; Nt = len(self.tgrid)
        rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]; Nr = len(rgrid)            
        self.zgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/zgrid'][:]; Nz = len(self.zgrid)            
        if full_resolution:
            kr_step = 1; Nr_max = Nr
        else:
            dr = r_resolution[1]; rmax = r_resolution[2]
            dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
            rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid) 
        
        # CUPRAD ouputs (z,t,r) (c-like, original Fortran is reversed)
        
        self.E_zrt = np.transpose(InputArchive[MMA.paths['CUPRAD_outputs'] +'/output_field'][:Nz,:,0:Nr_max:kr_step],
                                  axes=(0,2,1))# Arrays may be over-allocated by CUPRAD
        
        # hot-fix case of underallocated array, happens rarely
        if (self.E_zrt.shape[0] < Nz):
            Nz = self.E_zrt.shape[0]
            self.zgrid = self.zgrid[:Nz]
        

        self.inverse_GV = InputArchive[MMA.paths['CUPRAD_logs'] +'/inverse_group_velocity_SI'][()]
        self.VG_IR = 1.0/self.inverse_GV               
        self.rho0_init = 1e6 * mn.readscalardataset(InputArchive,
                         MMA.paths['CUPRAD_inputs'] +'/calculated/medium_effective_density_of_neutral_molecules','N')
        self.Ip_eV = InputArchive[MMA.paths['CUPRAD_inputs'] +'/ionization_ionization_potential_of_neutral_molecules'][()]
        self.pressure_mbar = 1e3*InputArchive[MMA.paths['CUPRAD_inputs'] +'/medium_pressure_in_bar'][()]; self.pressure_string = "{:.1f}".format(self.pressure_mbar)+' mbar'
        try:
            self.preionisation_ratio = InputArchive[MMA.paths['global_inputs'] +'/pre_ionised/initial_electrons_ratio'][()]
        except:
            self.preionisation_ratio = 0
        self.preionisation_string = "{:.1f}".format(100*self.preionisation_ratio) + ' %'
        
        self.effective_neutral_particle_density = 1e6*mn.readscalardataset(InputArchive, MMA.paths['CUPRAD_inputs'] +'/calculated/medium_effective_density_of_neutral_molecules','N')
        
        if 'density_mod' in InputArchive[MMA.paths['global_inputs']].keys():
            self.density_mod_profile_relative = InputArchive[MMA.paths['global_inputs']+'/density_mod/table'][:]
            self.density_mod_profile_mbar = self.pressure_mbar*self.density_mod_profile_relative
            try: self.density_mod_zgrid = InputArchive[MMA.paths['global_inputs']+'/density_mod/zgrid'][:]
            except: pass
        
            try: self.density_mod_rgrid = InputArchive[MMA.paths['global_inputs']+'/density_mod/rgrid'][:]
            except: pass
        
        self.w0_entry = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_beamwaist_entry',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_beamwaist_entry')
        
        self.pulse_duration_entry = mn.h5_seek_for_scalar(InputArchive,'N',
                                        MMA.paths['CUPRAD_inputs'] +'/laser_pulse_duration_in_1_e_Efield',
                                        MMA.paths['CUPRAD_inputs'] +'/calculated/laser_pulse_duration_in_1_e_Efield')

        try:
            self.gas_type = mn.h5_seek_for_scalar(InputArchive,'S',
                                        MMA.paths['global_inputs'] + '/gas_preset',
                                        MMA.paths['CUPRAD_inputs'] + '/gas_preset')
        except:
            self.gas_type = 'not specified'

        self.rgrid = rgrid
        self.Nr = Nr; self.Nt = Nt; self.Nz = Nz
        
        # Further analyses that may be stored in various directions
        self.Gaussian_focus = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_focus_position_Gaussian',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_focus_position_Gaussian')
        self.Gaussian_focus_string = "{:.1f}".format(1e3*self.Gaussian_focus) + ' mm' # 'z='+"{:.1f}".format(1e3*res.zgrid[0])
        
        self.Intensity_entry = mn.h5_seek_for_scalar(InputArchive,'N',
                              MMA.paths['CUPRAD_inputs'] +'/laser_intensity_entry',
                              MMA.paths['CUPRAD_inputs'] +'/calculated/laser_intensity_entry')
        self.Intensity_entry_string = "{:.1f}".format(1e-18*self.Intensity_entry) + ' 1e18 W/m2'
        
        # Idealised beam Rayleigh range (defined only if Gaussian waist is available)
        
        self.energy = InputArchive[MMA.paths['CUPRAD'] +'/longstep/energy'][:]
        self.energy_zgrid = InputArchive[MMA.paths['CUPRAD'] +'/longstep/z_buff'][:]
            
        try:
            Gaussian_w0 = mn.readscalardataset(InputArchive, MMA.paths['CUPRAD_inputs'] +'/laser_focus_beamwaist_Gaussian','N')
            self.Gaussian_zR = np.pi*(Gaussian_w0**2)/(1e-2*mn.readscalardataset(InputArchive,MMA.paths['CUPRAD_inputs'] +'/laser_wavelength','N'))
        except:
            self.Gaussian_zR = np.NaN

        try:
            self.Intensity_Gaussian_focus = mn.readscalardataset(InputArchive, MMA.paths['CUPRAD_inputs'] +'/laser_focus_intensity_Gaussian','N')
            self.Intensity_Gaussian_focus_string = "{:.1f}".format(1e-18*self.Intensity_Gaussian_focus) + ' 1e18 W/m2'
        except:
            self.Intensity_Gaussian_focus = np.NaN   
            self.Intensity_Gaussian_focus_string = "xxx"
            
            
        def co_moving_t_grid(z_def):
            """Return the time grid shifted to a frame co-moving with the pulse.

            Parameters
            ----------
            z_def : int or float
                Propagation position. If an integer is given, it is interpreted as an index
                into ``self.zgrid``. Otherwise, it is interpreted directly as a propagation
                distance.

            Returns
            -------
            ndarray
                Time grid corrected by the difference between propagation at the infrared
                group velocity and propagation at the vacuum speed of light.
            """
            if isinstance(z_def, int): z = self.zgrid[z_def]
            else: z = z_def
                
            delta_t_lab = z/self.VG_IR
            delta_t_vac  = z/units.c_light
            return self.tgrid + delta_t_lab - delta_t_vac
        
        self.co_moving_t_grid = co_moving_t_grid
        
    def vacuum_shift(self,output='replace'):
        """Shift the electric field to coordinates moving at the vacuum speed of light.

        For each ``z`` and ``r`` position, the temporal field is Fourier-transformed,
        phase-shifted according to the difference between the group-velocity frame and
        the vacuum-light frame, and transformed back to the time domain.

        Parameters
        ----------
        output : {'replace', 'return', 'add'}, optional
            Output mode. ``'replace'`` overwrites ``self.E_zrt`` with the shifted field,
            ``'return'`` returns the shifted field, and ``'add'`` stores it as
            ``self.E_zrt_vac``.

        Returns
        -------
        ndarray or None
            Shifted electric field if ``output='return'``; otherwise ``None``.

        Raises
        ------
        ValueError
            If ``output`` is not one of the supported modes.
        """
        E_vac = np.zeros(self.E_zrt.shape)   
        for k1 in range(self.Nz):
            delta_z = self.zgrid[k1] # local shift
            delta_t_lab = self.inverse_GV*delta_z # shift to the laboratory frame
            delta_t_vac = delta_t_lab - delta_z/units.c_light # shift to the coordinates moving by c.
            for k2 in range(self.Nr):
                ogrid_nn, FE_s, NF = mn.fft_t_nonorm(self.tgrid, self.E_zrt[k1,k2,:]) # transform to omega space        
                FE_s = np.exp(1j*ogrid_nn*delta_t_vac) * FE_s # phase factor        
                tnew, E_s = mn.ifft_t_nonorm(ogrid_nn,FE_s,NF)
                E_vac[k1,k2,:] = E_s.real
        
        if (output == 'replace'):      self.E_zrt = E_vac 
        elif (output == 'return'):     return E_vac 
        elif (output == 'add'):        self.E_zrt_vac = E_vac 
        else: raise ValueError('wrongly specified output for the vacuum shift.')
        

    def complexify_envel(self,output='return'):
        """Compute the complex envelope of the propagated electric field.

        The method constructs the analytic signal of each temporal field trace and
        removes the fast carrier oscillation at ``self.omega0``.

        Parameters
        ----------
        output : {'return', 'add'}, optional
            Output mode. ``'return'`` returns the complex envelope, while ``'add'``
            stores it as ``self.E_zrt_cmplx_envel``.

        Returns
        -------
        ndarray or None
            Complex envelope array with the same shape as ``self.E_zrt`` if
            ``output='return'``; otherwise ``None``.

        Raises
        ------
        ValueError
            If ``output`` is not one of the supported modes.
        """
        E_zrt_cmplx_envel = np.zeros(self.E_zrt.shape,dtype=complex)
        rem_fast_oscillations = np.exp(-1j*self.omega0*self.tgrid)
            
        for k1 in range(self.Nz):
            for k2 in range(self.Nr):
                E_zrt_cmplx_envel[k1,k2,:] = rem_fast_oscillations*mn.complexify_fft(self.E_zrt[k1,k2,:])
        
        if (output == 'return'):     return E_zrt_cmplx_envel
        elif (output == 'add'):      self.E_zrt_cmplx_envel = E_zrt_cmplx_envel
        else: raise ValueError('wrongly specified output for the vacuum shift.') 
        
    def get_Fluence(self, InputArchive, fluence_source='file'):
        """Load or compute the radial fluence diagnostic.

        The result is stored in ``self.Fluence`` as a lightweight container with fields
        ``value``, ``zgrid``, ``rgrid``, and ``units``.

        Parameters
        ----------
        InputArchive : h5py-like file object
            Open archive containing CUPRAD data and optional stored fluence diagnostics.
        fluence_source : {'file', 'computed'}, optional
            Source of the fluence. ``'file'`` reads the CUPRAD long-step fluence output.
            ``'computed'`` integrates ``|E(z, r, t)|**2`` over time using the loaded
            electric field.

        Notes
        -----
        The computed fluence is evaluated as ``c * eps0 * integral |E|**2 dt`` and is
        stored in SI units of ``J/m2``.
        """
        self.Fluence = empty_class()
        if (fluence_source == 'file'):
            self.Fluence.value = InputArchive[MMA.paths['CUPRAD'] +'/longstep/fluence'][:,:]
            self.Fluence.zgrid = InputArchive[MMA.paths['CUPRAD'] +'/longstep/zgrid_analyses2'][:]
            self.Fluence.rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]
            self.Fluence.units = 'C.U.'                
        
        elif (fluence_source == 'computed'):                
            self.Fluence.zgrid = self.zgrid
            self.Fluence.rgrid = self.rgrid
            self.Fluence.value = np.zeros((self.Nr, self.Nz))
            for k1 in range(self.Nz):
                for k2 in range(self.Nr):
                    # self.Fluence.value[k2, k1] = sum(abs(self.E_trz[:, k2, k1])**2)
                    self.Fluence.value[k2, k1] = units.c_light*units.eps0 * np.trapz(abs(self.E_zrt[k1, k2, :])**2,self.tgrid)
            self.Fluence.units = 'J/m2'

    def get_plasma(self, InputArchive, r_resolution=[True]): # analogy to the fields
        """Load the plasma-density output from the CUPRAD archive.

        The plasma data are stored in ``self.plasma`` as a lightweight container with
        the array ``value_zrt`` and the corresponding grids and dimensions.

        Parameters
        ----------
        InputArchive : h5py-like file object
            Open archive containing CUPRAD plasma output.
        r_resolution : list, optional
            Radial loading specification. ``[True]`` loads the full radial grid.
            Otherwise, use ``[False, dr, rmax]`` to down-sample the radial direction
            while loading.

        Attributes Created
        ------------------
        self.plasma.value_zrt : ndarray
            Plasma output array ordered as ``z, r, t``.
        self.plasma.tgrid, self.plasma.rgrid, self.plasma.zgrid : ndarray
            Temporal, radial, and propagation grids for the plasma output.
        """
        full_resolution = r_resolution[0]
        self.plasma = empty_class()
        
        self.plasma.tgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/tgrid'][:]; Nt = len(self.tgrid)
        rgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/rgrid'][:]; Nr = len(rgrid)            
        self.plasma.zgrid = InputArchive[MMA.paths['CUPRAD_outputs'] +'/zgrid'][:]; Nz = len(self.zgrid)            
        if full_resolution:
            kr_step = 1; Nr_max = Nr
        else:
            dr = r_resolution[1]; rmax = r_resolution[2]
            dr_file = rgrid[1]-rgrid[0]; kr_step = max(1,int(np.floor(dr/dr_file))); Nr_max = mn.FindInterval(rgrid, rmax)
            rgrid = rgrid[0:Nr_max:kr_step]; Nr = len(rgrid) 
            
        self.plasma.value_zrt = np.transpose(InputArchive[MMA.paths['CUPRAD_outputs'] +'/output_plasma'][:Nz,:,0:Nr_max:kr_step],
                                             axes=(0,2,1)) # [:,0:Nr_max:kr_step,:Nz] # Arrays may be over-allocated by CUPRAD
        
        
        # self.plasma.value_trz = self.plasma.value_trz.transpose(1,2,0) # (1,2,0) # hot-fix to reshape
        
        self.plasma.rgrid = rgrid
        self.plasma.Nr = Nr; self.plasma.Nt = Nt; self.plasma.Nz = Nz

    def compute_spectrum(self,output='add',compute_dE_domega = False):
        """Compute temporal spectra of the propagated electric field.

    The method Fourier-transforms each temporal trace of ``self.E_zrt`` and stores
    or returns the resulting spectral field ordered as ``z, r, omega``. Optionally,
    it also computes a radially integrated spectral quantity for each propagation
    position.

    Parameters
    ----------
    output : {'add', 'return'}, optional
        Output mode. ``'add'`` stores the spectrum as ``self.FE_zrt`` and, if
        requested, ``self.dE_domega``. ``'return'`` returns the computed arrays.
    compute_dE_domega : bool, optional
        If ``True``, also compute the radial integral of ``abs(FE_zrt)**2`` for each
        ``z`` and frequency.

    Returns
    -------
    ndarray or tuple of ndarray or None
        If ``output='return'`` and ``compute_dE_domega`` is ``False``, returns
        ``FE_zrt``. If ``compute_dE_domega`` is ``True``, returns
        ``(FE_zrt, dE_domega)``. Otherwise returns ``None``.

    Raises
    ------
    ValueError
        If ``output`` is not one of the supported modes.
    """
        self.ogrid = mn.fft_t(self.tgrid, self.E_zrt[0,0,:])[0]
        
        No = len(self.ogrid); Nr = len(self.rgrid); Nz = len(self.zgrid)
                
        FE_zrt = np.zeros((Nz,Nr,No),dtype=complex)
        
        for k1 in range(Nz):
            for k2 in range(Nr):
                FE_zrt[k1,k2,:] = mn.fft_t(self.tgrid, self.E_zrt[k1,k2,:])[1]
  
        
        if compute_dE_domega:
            dE_domega = np.empty((Nz,No))
            
            for k1 in range(Nz):
                for k2 in range(No):
                    dE_domega[k1,k2] = np.trapz(np.abs(FE_zrt[k1,:,k2])**2,self.rgrid)
            
            if (output == 'return'):     return FE_zrt, dE_domega
            elif (output == 'add'):      self.FE_zrt = FE_zrt; self.dE_domega = dE_domega
            else: raise ValueError('wrongly specified output for the vacuum shift.') 
        else:
            if (output == 'return'):     return FE_zrt
            elif (output == 'add'):      self.FE_zrt = FE_zrt
            else: raise ValueError('wrongly specified output for the vacuum shift.') 

    def get_ionisation_model(self, InputArchive):
        """Load the tabulated ionisation model from the CUPRAD archive.

        The ionisation data are stored in ``self.ionisation_model`` as a lightweight
        container with the electric-field grid and the corresponding ionisation rates.

        Parameters
        ----------
        InputArchive : h5py-like file object
            Open archive containing the CUPRAD ionisation-model group.
        """
        self.ionisation_model = empty_class()
        self.ionisation_model.Egrid =            InputArchive[MMA.paths['CUPRAD_ionisation_model'] +'/Egrid'][:]
        self.ionisation_model.ionisation_rates = InputArchive[MMA.paths['CUPRAD_ionisation_model'] +'/ionisation_rates'][:]
        
        
def add_print_parameter(parameter,data):
    """Return a formatted string for a selected simulation parameter.

    Parameters
    ----------
    parameter : str
        Parameter selector. Supported values are ``'pressure'``,
        ``'preionisation'``, ``'focus_in_medium'``, ``'intensity_entry'``, and
        ``'intensity_Gaussian_focus'``.
    data : get_data-like object
        Data object containing the formatted parameter strings.

    Returns
    -------
    str
        Formatted parameter string, or an empty string for an unsupported selector.
    """
    if (parameter=='pressure'): return data.pressure_string
    elif (parameter=='preionisation'): return data.preionisation_string
    elif (parameter=='focus_in_medium'): return data.Gaussian_focus_string
    elif (parameter=='intensity_entry'): return data.Intensity_entry_string
    elif (parameter=='intensity_Gaussian_focus'): return data.Intensity_Gaussian_focus_string
    else: return ''


def create_param_string(params,data):
    res = ''
    for param in params:
        curr = add_print_parameter(param,data)
        if not(len(curr)==0):
            res = res + ', ' + curr     
    return res   


def measure_beam(grid, beam, measure, *args, measured_axis = 0):
    N0, N1 = beam.shape
    if (measured_axis == 0):
        radius = np.zeros(N1)
        for k1 in range(N1):
            radius[k1] = measure(grid,beam[:,k1],*args)
    if (measured_axis == 1):
        radius = np.zeros(N0)
        for k1 in range(N0):
            radius[k1] = measure(grid,beam[k1,:],*args)
    
    return radius


