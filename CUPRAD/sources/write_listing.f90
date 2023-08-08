MODULE write_listing
  USE normalisation

  INTEGER(4) p,switch_dispersion,switch_start
  CHARACTER(LEN=15)      :: indexfile
  
CONTAINS

  SUBROUTINE write_listingfile
    IMPLICIT NONE

    OPEN(unit=100,FILE='listing',form='formatted',STATUS='NEW')

    WRITE(100,'(t5,a,t20,a/)') 'NUMERICAL PARAMETERS'
    WRITE(100,'(a,t50,i3)') 'number of processors' ,num_proc
    WRITE(100,'(a,t50,i5)') 'number of points in t',dim_t
    WRITE(100,'(a,t50,i5)') 'number of points in r',dim_r
    WRITE(100,'(a,t50,es12.4)') 'length of window for t, normalized to t_p',lt
    WRITE(100,'(a,t50,es12.4)') 'length of window for r, normalized to w_0',lr
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'step width in time',lt*tp_fs_phys/dim_t,'fs'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'x-step width',lr*w0_cm_phys*10000/dim_r,'mum'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'physical distance of propagation',proplength_m_phys,'m'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'physical output distance for matlab files',outlength_m_phys,'m'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'radius for diagnostics',rfil_mm_phys,'mm'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'physical first stepwidth (code might reduce it)',delta_z_mm_phys,'mm'
    WRITE(100,'(a,t50,i3)') 'number of absorber points in time',absorb
    WRITE(100,'(a,t50,es12.4)') 'phase threshold for increasing delta_z',increase
    WRITE(100,'(a,t50,es12.4)') 'phase threshold for decreasing delta_z',decrease
    WRITE(100,'(a,t50,i4)') 'output distance in z-steps for plasma-channel',rhodist
    WRITE(100,'(a,t50,es12.4/)') 'approx runtime in hours',time_limit

    WRITE(100,'(t5,a,t20,a/)') 'PHYSICAL PARAMETERS'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'wavelenght lambda0',lambda0_cm_phys,'cm'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'beamwaist',w0_cm_phys,'cm'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'focal length in the medium (0 = no lense)',f_cm_phys,'cm'
    WRITE(100,'(a,t50,es12.4)') 'lense parameter exp(-L((x/w0)^2+(y/w0)^2))',lense_factor
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'pulse duration',tp_fs_phys,'fs'
    SELECT CASE (switch_start)
    CASE(1) 
       WRITE(100,'(a,t50,es12.4)') 'chirp parameter exp(-C(t/t_p)^2)',chirp_factor
       WRITE(100,'(a,t50,es12.4)') 'ratio Pin/Pcr',numcrit
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise spatial (amplitude)',int( noise_s * 100.d0),'%'
       WRITE(100,'(a,t50,i3)') 'SuperGaussian-in-space' ,super_N
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise temporal (amplitude)',int( noise_t * 100.d0),'%'
       WRITE(100,'(a,t50,i3)') 'SuperGaussian-in-time' ,super_t
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise (amplitude)',int( noise * 100.d0),'%'
    CASE(2) 
       WRITE(100,'(a,t50,es12.4)') 'chirp parameter exp(-C(t/t_p)^2)',chirp_factor
       WRITE(100,'(a,t50,es12.4)') 'ratio Pin/Pcr',numcrit
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise spatial (amplitude)',int( noise_s * 100.d0),'%'
       WRITE(100,'(a,t50,i3)') 'SuperGaussian-in-space' ,super_N
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise temporal (amplitude)',int( noise_t * 100.d0),'%'
       WRITE(100, '(a,t50)') 'input temporal profile from file '//inputfilename_t
       WRITE(100,'(a,t50,i4, 1x,a)') 'field noise (amplitude)',int( noise * 100.d0),'%'
    CASE(3)
       WRITE(100, '(a,t50)') 'CONTINUATION!!!!!!!!!!!!'
       ! WRITE(100, '(a,t50)') 'initial condition from files '//inputfilename_c
       WRITE(100,'(a,t50,es12.4)') 'amplitude factor for medium change',restartamp
    CASE(4)
       WRITE(100, '(a,t50)') 'CONTINUATION WITH FLAT PHASE!!!!!!!!!!!!'
       ! WRITE(100, '(a,t50)') 'initial condition from files '//inputfilename_c
       WRITE(100,'(a,t50,es12.4)') 'amplitude factor for medium change',restartamp
    END SELECT
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'pressure',pressure,'bar'
    WRITE(100,'(a,t50,es12.4)') 'refractive index at lambda0',n0
    SELECT CASE (switch_dispersion)
    CASE(1) 
       WRITE(100, '(a,t50)') 'Taylor expansion for the dispersion law'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'First order dispersion term - n0/c', delta_k_p_fs_per_cm_phys, 'fs/cm'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'First order dispersion term', k_p_fs_per_cm_phys, 'fs/cm'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'Second order dispersion term', k_pp_fs2_per_cm_phys, 'fs2/cm'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'Third order dispersion term',  k_ppp_fs3_per_cm_phys , 'fs3/cm'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'Fourth order dispersion term', k_pppp_fs4_per_cm_phys, 'fs4/cm'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'Fifth order dispersion term', k_ppppp_fs5_per_cm_phys, 'fs5/cm'
    CASE(2)
       WRITE(100, '(a,t50)') 'Sellmeier formula for silica, Agrawal' 
    CASE(3)
       WRITE(100, '(a,t50)') 'dispersion for argon, Dalgarno and Kingston' 
    CASE(4)
       WRITE(100, '(a,t50)') 'chi taken from file '//dispfilename
    CASE(5)
       WRITE(100, '(a,t50)') 'dispersion for air, Peck and Reeder'
    CASE(6)
       WRITE(100, '(a,t50)') 'dispersion for xenon, Dalgarno and Kingston' 
    CASE(7)
       WRITE(100, '(a,t50)') 'dispersion for neon, Dalgarno and Kingston' 
    CASE(8)
       WRITE(100, '(a,t50)') 'dispersion for KDPo, Handbook of Optics' 
    END SELECT
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'Kerr coefficient',n2_phys,'cm2/W'
    SELECT CASE (switch_dKerr)
    CASE(1)
       WRITE(100, '(a,t50)') 'Only instantaneous Kerr response'
    CASE(2)
       WRITE(100, '(a,t50)') 'Delayed Kerr response in exp[(T-t)/tk]'
       WRITE(100,'(a,t50,es12.4)') 'ratio delayed kerr /instantaneous kerr',xdk
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'time delayed kerr coefficient tk',tdk_fs_phys,'fs'
    CASE(3)
       WRITE(100, '(a,t50)') 'Delayed Kerr response in exp[(T-t)/tk]*sin[wr(T-t)]'
       WRITE(100,'(a,t50,es12.4)') 'ratio delayed kerr / instantaneous kerr',xdk
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'time delayed kerr coefficient tk',tdk_fs_phys,'fs'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'frequency in raman response wr',raman_phys,'fs-1'
    END SELECT
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'chi5 coefficient',n4_phys,'cm4/W2'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'critical plasma density',rhoc_cm3_phys,'cm-3'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'effective density of neutral molecules',rhont_cm3_phys,'cm-3'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'physical initial electron density',rho0_phys,'cm-3'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'gap potential for ionization of oxygen molecules',Ui_eV_phys,'eV'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'gap potential for ionization of oxygen molecules',Ui_au_phys,'au'
    SELECT CASE (switch_rho)
    CASE(8)
       WRITE(100,'(a,t50,i3)') 'CPR',switch_rho
       WRITE(100,'(a,t50,es12.4)') 'complex rotation used'
    CASE(1)
       WRITE(100,'(a,t50,i3)') 'Plasma equation solved with euler method',switch_rho
       WRITE(100,'(a,t50,i3)') 'number of photon needed to extract electrons',KK
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'coefficient for MPI',sigmak_phys,'s-1cm2K/W2'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'coefficient of multiphoton absorption',betak_phys,'cm2K-3/WK-1'
    CASE(2)
       WRITE(100,'(a,t50,i3)') 'Plasma equation solved semi-analyticaly',switch_rho
       WRITE(100,'(a,t50,i3)') 'number of photon needed to extract electrons',KK
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'coefficient for MPI',sigmak_phys,'s-1cm2K/W2'
       WRITE(100,'(a,t50,es12.4, 1x,a)') 'coefficient of multiphoton absorption',betak_phys,'cm2K-3/WK-1'
    CASE(3)
       WRITE(100,'(a,t50,i3)') 'PPT',switch_rho
       WRITE(100,'(a,t50,es12.4)') 'residue charge', residue_charge
       WRITE(100,'(a,t50,i3)') 'angular momentum',angular_momentum
    END SELECT
    SELECT CASE (switch_T)
    CASE(1)
       WRITE(100,'(a,t50)') 'no self-steepening and space-time focusing'
    CASE(2)
       WRITE(100,'(a,t50)') 'Operators T, T^-1 included'
    CASE(3)
       WRITE(100,'(a,t50)') 'UPPE terms included'
    CASE(4)
       WRITE(100,'(a,t50)') 'UPPE terms without harmonics included'
    END SELECT
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'electron collision time',tauc_fs_phys,'fs'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'cross section for inverse bremsstrahlung',sigma_cm2_phys,'cm2'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'rec. coefficient (linear) for MPI e',alpha_fs_phys,'fs-1'


    WRITE(100,'(a,t50,es12.4, 1x,a)') 'recombination coefficent (quadratic)',alphaquad_fscm3_phys,'fs-1cm3'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'photon energy',photon_energy_au_phys,'au'
    WRITE(100,'(a,t50,i3)') 'Number of photons for n-photon absoprtion',NN
    WRITE(100,'(a,t50,es12.4, 1x,a)')'The n-photon absorption cross-section', sigman_phys, 's-1cm2N/WN'
    WRITE(100,'(a,t50,es12.4, 1x,a/)')'Density of absorbent molecules', rhoabs_cm3_phys, 'cm-3'

    WRITE(100,'(t5,a,t20,a/)') 'CARACTERISTIC VALUES'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'Rayleigh length',z_rayleigh_cm_phys,'cm'
    WRITE(100,'(a,t50,es12.4, 1x,a)') 'critical power',Pcr_phys,'W'
    WRITE(100,'(a,t50,es12.4,1x,a/)') 'central wave number in vacuum',k0_phys,'cm-1'

    WRITE(100,'(t5,a,t20,a/)') 'ADIMMENSIONNED COEFFICIENTS'
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned distance of propagation',proplength
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned output distance matlab files',outlength
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned first stepwidth',delta_z
    WRITE(100,'(a,t50,es12.4)') 'adimensionned central wavenumber in the medium',rek0
    WRITE(100,'(a,t50,es12.4)') 'adimensionned inverse group velocity',rekp
    IF (switch_T.EQ.3) THEN
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned phase mismatch THG',deltak3omega(1)
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned part phase mismatch THG lin ',deltak3omega(2)
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned part phase mismatch THG exp ',deltak3omega(3)
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned phase mismatch FHG',deltak5omega(1)
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned part phase mismatch FHG lin ',deltak5omega(2)
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned part phase mismatch FHG exp ',deltak5omega(3)
    ENDIF
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for Kerr effect',c3
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for chi5 effect',c5
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for losses due to AI',gamma1
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for plasma defocusing',gamma2
    IF ((switch_rho.EQ.1).OR.(switch_rho.EQ.2)) THEN
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for MPA effect',muk
       WRITE(100,'(a,t50,es12.4)') 'adimmensionned MPI coefficient',beta_inv_2KK
    ENDIF
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned initial electron density',rho0
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned avalanche coefficient',nu
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned recombination coefficient',alpha

    WRITE(100,'(a,t50,es12.4)') 'adimensionned recombination coefficient (quadratic)',alphaquad
    WRITE(100,'(a,t50,es12.4)') 'adimmensionned inverse density of neutral molecules',rhoat_inv
    IF ((switch_dkerr.EQ.2).OR.(switch_dkerr.EQ.3)) WRITE(100,'(a,t50,es12.4)') 'adimmensionned coefficient for delay time',tdk
    IF (switch_dkerr.EQ.3) WRITE(100,'(a,t50,es12.4)') 'adimensionned raman frequency',raman
    WRITE(100,'(a,t50,es12.4)') 'adimensioned radius for diagnostics',rfil
    WRITE(100,'(a,t50,es12.4)') 'adimensionned central frequency',omega
    WRITE(100,'(a,t50,es12.4)') 'adimensionned central frequency of the box',omega_uppe
    WRITE(100,'(a,t50,es12.4)') 'adimensionned nonlinear losses coefficient', eta1
    WRITE(100,'(a,t50,es12.4/)') 'adimensionned coeff. for excited moelcules', eta2

    WRITE(100, '(a,t50)') 'index distribution taken from file '//indexfile
    CLOSE(100)

    RETURN
  END SUBROUTINE write_listingfile

END MODULE write_listing
