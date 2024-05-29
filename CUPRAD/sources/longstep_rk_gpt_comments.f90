MODULE long_step
    USE constants
    USE parameters
    USE density_module
    REAL(8) rhompi, rho1, rho2, rhoth, rhotr, rhofh, rhoslg2, rhoav
  CONTAINS
  
    ! \brief Subroutine for refractive index lookup and interpolation
    ! \param phase_index Output interpolated phase index
    ! \param r Radial position
    SUBROUTINE index_interpolation(phase_index, r)
      USE fields
      IMPLICIT NONE
  
      INTEGER(4) i_x, i_z, help
      REAL(8) phase_index, r
  
      IF ((z.GE.zz(1)).AND.(z.LE.zz(i_z_max)).AND.(r.LE.xx(i_x_max))) THEN
         IF (z.LT.zz(i_z_old)) THEN
            i_z_old = 2
         ENDIF
         IF (r.LT.xx(i_x_old)) THEN
            i_x_old = 2
         ENDIF
         DO i_z = i_z_old, i_z_max
            IF (z.LE.zz(i_z)) THEN
               help = i_z
               EXIT
            ENDIF
         ENDDO
         i_z_old = help   
         DO i_x = i_x_old, i_x_max
            IF (r.LE.xx(i_x)) THEN
               help = i_x
               EXIT
            ENDIF
         ENDDO
         i_x_old = help
         IF ((z.LT.zz(i_z_old-1)).OR.(z.GT.zz(i_z_old))) THEN
            print*, 'Error in z index interpolation', r, z
            READ(5)
         ENDIF
         IF ((r.LT.xx(i_x_old-1)).OR.(r.GT.xx(i_x_old))) THEN
            print*, 'Error in r index interpolation', r, z
            READ(5)
         ENDIF              
         phase_index = (zz(i_z_old) - z) * (xx(i_x_old) - r) * Indice_norm(i_x_old-1, i_z_old-1)
         phase_index = phase_index + (z - zz(i_z_old-1)) * (xx(i_x_old) - r) * Indice_norm(i_x_old-1, i_z_old)
         phase_index = phase_index + (zz(i_z_old) - z) * (r - xx(i_x_old-1)) * Indice_norm(i_x_old, i_z_old-1)
         phase_index = phase_index + (z - zz(i_z_old-1)) * (r - xx(i_x_old-1)) * Indice_norm(i_x_old, i_z_old)
         phase_index = phase_index / ((zz(i_z_old) - zz(i_z_old-1)) * (xx(i_x_old) - xx(i_x_old-1)))
      ELSE
         phase_index = 0.D0
      ENDIF
  
      RETURN
    END SUBROUTINE index_interpolation
  
    
    ! \brief Subroutine to apply absorption using cosh in time on the boundary
    SUBROUTINE absorbation
      USE fields
      USE mpi_stuff
      IMPLICIT NONE
  
      INTEGER(4) j, l
  
      DO l = dim_r_start(num_proc), dim_r_end(num_proc)
         DO j = 1, dim_t
            e(j, l) = e(j, l) * bound_t(j)
         ENDDO
      ENDDO
  
      RETURN
    END SUBROUTINE absorbation
  
  
    ! \brief Subroutine to solve a tridiagonal system of equations
    ! \param dim Dimension of the system
    ! \param DLn Lower diagonal element at boundary
    ! \param D Diagonal elements
    ! \param DL Lower diagonal elements
    ! \param DU Upper diagonal elements
    ! \param B Right-hand side vector
    SUBROUTINE tridag(dim, DLn, D, DL, DU, B)
      IMPLICIT NONE
  
      INTEGER(4) dim
      COMPLEX(8) DLn
      COMPLEX(8) D(dim), DL(2:dim-1), DU(dim-1), B(dim)
      INTEGER(4) j
      COMPLEX(8) bet, gam(2:dim)
  
      bet = D(1)
      B(1) = B(1) / bet
      DO j = 2, dim-1
         gam(j) = DU(j-1) / bet
         bet = D(j) - DL(j) * gam(j)
         B(j) = (B(j) - DL(j) * B(j-1)) / bet
      ENDDO
      gam(dim) = DU(dim-1) / bet
      bet = D(dim) - DLn * gam(dim)
      B(dim) = (B(dim) - DLn * B(dim-1)) / bet
      DO j = dim-1, 1, -1
         B(j) = B(j) - gam(j+1) * B(j+1)
      ENDDO
  
      RETURN
    END SUBROUTINE tridag
  
    ! \brief Subroutine to perform the Crank-Nicolson method for a single time step
    ! \param B Right-hand side vector
    ! \param k Spatial index
    SUBROUTINE cn(B, k)
      USE fields
      USE parameters
      IMPLICIT NONE
  
      INTEGER(4) k
      COMPLEX(8) B(dim_r)
      INTEGER(4) j
      COMPLEX(8) DLn, help_1, help_2
  
      help_1 = B(1)
      B(1) = (CMPLX(1.D0, 0.D0, 8) + CMPLX(0.D0, -2.D0, 8) * delta_rel(1, k)) * B(1) + CMPLX(0.D0, 2.D0, 8) * delta_rel(1, k) * B(2)
      DO j = 2, dim_r-1
         help_2 = help_1
         help_1 = B(j)
         B(j) = -DL(j-1, k) * help_2 + (CMPLX(1.D0, 0.D0, 8) + CMPLX(0.D0, -1.D0, 8) * delta_rel(j, k)) * B(j) - DU(j, k) * B(j+1)
      ENDDO
      B(dim_r) = CMPLX(0.D0, 0.D0, 8) * B(dim_r)
  
      IF (help_2.NE.(0.D0, 0.D0)) THEN
         DLn = -help_1 / help_2
         IF (AIMAG(DLn).GT.0.D0) DLn = CMPLX(-1.D0, 0.D0, 8)
      ELSE
         DLn = (-1.D0, 0.D0)
      ENDIF
  
      CALL tridag(dim_r, DLn, D(:, k), DL(:, k), DU(:, k), B)
  
      RETURN
    END SUBROUTINE cn
  
    ! \brief Subroutine to perform the linear propagation in z-direction
    SUBROUTINE mult_propagator
      USE fields
      USE mpi_stuff
      IMPLICIT NONE
  
      INTEGER(4) k, l
  
      DO k = dim_t_start(num_proc), dim_t_end(num_proc)
         CALL cn(efft(1:dim_r, k), k)
         DO l = 1, dim_r
           efft(l, k) = efft(l, k) * p_t(l, k) ! Apply the dispersion
         ENDDO
      ENDDO
  
      RETURN
    END SUBROUTINE mult_propagator
  
    ! \brief Subroutine to calculate N-photon absorption (to be removed)
    ! \param rhoabs Absorption factor
    ! \param mediumabs Medium absorption
    ! \param eti Initial electric field
    ! \param etip1 Final electric field
    SUBROUTINE calc_absorption(rhoabs, mediumabs, eti, etip1)
      REAL(8) :: rhoabs, mediumabs, eti, etip1
      REAL(8) :: intF
  
      IF (eta1.NE.0.D0) THEN
         intF = -eta2 * 0.5D0 * (eti**NN + etip1**NN) * delta_t  
         rhoabs = 1.D0 - (1.D0 - rhoabs) * EXP(intF)
         IF (rhoabs.LT.0.D0) rhoabs = 0.D0
         mediumabs = eta1 * etip1**(NN-1) * (1.D0 - rhoabs)
      ENDIF
  
      RETURN
    END SUBROUTINE calc_absorption
  
    ! \brief Subroutine to calculate the density of ionized particles
    ! \param rho Density of ionized particles
    ! \param mpa Multi-photon absorption coefficient
    ! \param eti Initial electric field
    ! \param etip1 Final electric field
    ! \param l Radial index
    SUBROUTINE calc_rho(rho, mpa, eti, etip1, l)
      USE PPT
      USE External_ionisation_table
      IMPLICIT NONE
  
      INTEGER(4) :: l
      REAL(8) :: rho, mpa
      REAL(8) :: eti, etip1
      REAL(8) intF, var1, rhosave
  
      rhosave = rho
      SELECT CASE (switch_rho)
      CASE(8)
         CALL interpolate_ext(var1, mpa, etip1)
         intF = (nu * 0.5D0 * (etip1 + eti) * density_mod(l) - rhoat_inv * var1 - alpha) * delta_t
         rho = rho * EXP(intF) + var1 * density_mod(l) * delta_t
         mpa = mpa * (density_mod(l) - rhosave * rhoat_inv)
      CASE(1)
         rho = rho + (nu * rho * eti * density_mod(l) + beta_inv_2KK * eti**KK * (density_mod(l) - rho * rhoat_inv) - alpha * rho) * delta_t
         mpa = muk * etip1**(KK-1) * (density_mod(l) - rhosave * rhoat_inv)
      CASE(2)
         var1 = 0.5D0 * beta_inv_2KK * (etip1**KK + eti**KK)
         intF = (nu * 0.5D0 * (etip1 + eti) * density_mod(l) - rhoat_inv * var1 - alpha) * delta_t
         rho = rho * EXP(intF) + var1 * density_mod(l) * delta_t
         mpa = muk * etip1**(KK-1) * (density_mod(l) - rhosave * rhoat_inv)
      CASE(3)
         CALL interpolate_ppt(var1, mpa, etip1)
         intF = (nu * 0.5D0 * (etip1 + eti) * density_mod(l) - rhoat_inv * var1 - alpha) * delta_t
         rho = rho * EXP(intF) + var1 * density_mod(l) * delta_t
         mpa = mpa * (density_mod(l) - rhosave * rhoat_inv)
      END SELECT
      rho = rho - alphaquad * rhosave**2 * delta_t
      IF (rho.LT.0.D0) rho = 0.D0
  
      IF (density_mod(l) .GT. 0.D0) THEN
        IF ((rho * rhoat_inv / density_mod(l)) .GT. 1.D0) rho = density_mod(l) / rhoat_inv
      ENDIF
  
      IF (rhompi.LT.0.D0) rhompi = 0.D0
      IF (rho1.LT.0.D0) rho1 = 0.D0
      IF (rho2.LT.0.D0) rho2 = 0.D0
      IF (rhoth.LT.0.D0) rhoth = 0.D0
      IF (rhotr.LT.0.D0) rhotr = 0.D0
      IF (rhofh.LT.0.D0) rhofh = 0.D0
      IF (rhoslg2.LT.0.D0) rhoslg2 = 0.D0
      IF (rhoav.LT.0.D0) rhoav = 0.D0
  
      RETURN
    END SUBROUTINE calc_rho
  
    ! \brief Subroutine to calculate the Kerr effect
    ! \param delkerr Kerr effect term
    ! \param delkerrp Previous Kerr effect term
    ! \param eti Initial electric field
    ! \param etip1 Final electric field
    SUBROUTINE calc_delkerr(delkerr, delkerrp, eti, etip1)
      IMPLICIT NONE
  
      REAL(8) :: delkerr, delkerrp
      REAL(8) :: eti
      REAL(8) :: etip1
  
      SELECT CASE (switch_dKerr)
      CASE(1)
         CONTINUE
      CASE(2)
         delkerr = expt1 * delkerr + expt2 * eti + expt3 * etip1
      CASE(3)
         delkerr = expt1 * delkerr + expt2 * delkerrp + expt3 * eti + expt4 * etip1
         delkerrp = expt1p * delkerr + expt2p * delkerrp + expt3p * eti + expt4p * etip1
      END SELECT
  
      RETURN
    END SUBROUTINE calc_delkerr
  
    ! \brief The main computational subroutine
    SUBROUTINE mult_phase
      USE fields
      USE fft
      USE mpi_stuff
      USE HDF5
      USE HDF5_helper
      USE longstep_vars
      USE linked_list
      USE ll_data
      USE ppt
      USE normalization
      USE pre_ionised
      USE h5namelist
      IMPLICIT NONE
  
      INTEGER(4) j, l, k, k1
      REAL(8) phase, maxphase_part, peakmax_part, energy_part, energy_fil_part, rhomax_part, delkerr, delkerrp, rhotemp, mpa, r, phase_p, phase_j, losses_j, phase_index
      REAL(8) mediumabs, rhoabs_max_part, rhoabstemp
  
      ! For storing in HDF5
      INTEGER(HID_T) :: file_id       !< File identifier 
      INTEGER(HID_T) :: group_id      !< Group identifier 
      INTEGER :: error                !< HDF error
      LOGICAL :: group_status         !< Boolean indicating if the group exists
      
      ! Names of datasets and groups
      CHARACTER(*), PARAMETER :: groupname = longstep_grpname     
      CHARACTER(*), PARAMETER :: rhoabs_max_dset_name = longstep_grpname // "/rhoexcmax_max"
      CHARACTER(*), PARAMETER :: peakmax_dset_name = longstep_grpname // "/peakmax"
      CHARACTER(*), PARAMETER :: energy_dset_name = longstep_grpname // "/energy"
      CHARACTER(*), PARAMETER :: energy_fil_dset_name = longstep_grpname // "/energy_fil"
      CHARACTER(*), PARAMETER :: rhomax_dset_name = longstep_grpname // "/rhomax"
      CHARACTER(*), PARAMETER :: powmax_dset_name = longstep_grpname // "/powmax"
      CHARACTER(*), PARAMETER :: z_buff_dset_name = longstep_grpname // "/z_buff"
      CHARACTER(*), PARAMETER :: every_rhodist_z_dset_name = longstep_grpname // "/zgrid_analyses2"
      CHARACTER(*), PARAMETER :: onax_t_dset_name = longstep_grpname // "/Efied_onaxis"
      
      ! Arrays for temporary storage and dataset writing
      INTEGER(HSIZE_T), DIMENSION(1) :: new_dims, memspace_dims, offset, hyperslab_size
      REAL, DIMENSION(1,2) :: powmax_data
      REAL, ALLOCATABLE :: onax_t_data(:,:)
  
      ! Variables for linked list buffering
      REAL, ALLOCATABLE, TARGET :: fluence_data(:)
      REAL, ALLOCATABLE, TARGET :: plasma_channel_data(:)
      REAL, ALLOCATABLE, TARGET :: losses_plasma_data(:)
      REAL, ALLOCATABLE, TARGET :: losses_ionization_data(:)
     
      ! Initialization phase
      rho = 0.D0
      rhoabs = 0.D0 
      rhoabs_max_part = 0.D0
      losses_plasma = 0.D0
      losses_ionization = 0.D0
      peakmax_part = 0.D0
      rhomax_part = 0.D0
      energy_part = 0.D0
      energy_fil_part = 0.D0
      maxphase_part = 0.D0
  
      DO l = dim_r_start(num_proc), dim_r_end(num_proc)
         r = REAL(l-1) * delta_r
         e_2 = ABS(e(1:dim_t, l))
         peakmax_part = MAX(peakmax_part, MAXVAL(e_2))
         e_2 = e_2 ** 2
         fluence(l) = SUM(e_2)
         energy_part = energy_part + fluence(l) * REAL(l-1, 8)
         IF (rfil.GT.r) energy_fil_part = energy_fil_part + fluence(l) * REAL(l-1, 8)
         delkerr = 0.D0
         delkerrp = 0.D0
        IF (apply_pre_ionisation) THEN
           rhotemp = initial_electron_density(r, z, l, dim_r_start(num_proc))
        ELSE
           rhotemp = 0.D0
        ENDIF
  
         rhompi = 0.D0
         rho1 = 0.D0
         rho2 = 0.D0
         rhoth = 0.D0
         rhotr = 0.D0
         rhofh = 0.D0
         rhoslg2 = 0.D0
         rhoav = 0.D0
  
         rhoabstemp = 0.D0
         mediumabs = 0.D0
  
         ! Physical effects
         CALL index_interpolation(phase_index, r)
         phase_index = phase_index * delta_zh
         CALL calc_rho(rhotemp, mpa, 0.D0, e_2(1), l)
         CALL calc_absorption(rhoabstemp, mediumabs, 0.D0, e_2(1))
  
         DO j = 1, dim_t
            phase_p = (c3i * e_2(j) + c3d * delkerr - c5 * e_2(j)**2) * ((density_mod(l) - rhotemp * rhoat_inv) + ions_Kerr_ratio * rhotemp * rhoat_inv) * delta_zh
            phase_j = -gamma2 * rhotemp * delta_zh
            losses_j = -gamma1 * density_mod(l) * rhotemp * delta_zh
            rho(l) = MAX(rho(l), rhotemp)
            rhoabs(l) = MAX(rhoabs(l), rhoabstemp)
            losses_plasma(l) = losses_plasma(l) + 2.D0 * e_2(j) * gamma1 * density_mod(l) * rhotemp
            losses_ionization(l) = losses_ionization(l) + 2.D0 * e_2(j) * mpa
            phase = phase_p + phase_j + phase_index
            maxphase_part = MAX(maxphase_part, ABS(phase))
  
            SELECT CASE (switch_T)
            CASE(1)
               etemp(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_zh * (mpa + mediumabs), phase, 8))
            CASE(2)
               ptemp(j, l) = e(j, l) * CMPLX(0.D0, phase_p)
               jtemp(j, l) = e(j, l) * CMPLX(0.D0, phase_j)
               etemp(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_zh * (mpa + mediumabs), phase_index, 8))
            CASE(3)
               ptemp(j, l) = e(j, l) * CMPLX(0.D0, phase_p) + (hfac(j, 4) * CONJG(hfac(j, 0) * e(j, l)) * CMPLX(0.D0, phase_p) &
                    + (hfac(j, 1) * (hfac(j, 0) * e(j, l))**3 + hfac(j, 2) * (hfac(j, 0) * e(j, l))**3 * e_2(j)) &
                    + (hfac(j, 3) * (hfac(j, 0) * e(j, l))**5)) * ((density_mod(l) - rhotemp * rhoat_inv) + ions_Kerr_ratio * rhotemp * rhoat_inv)
               jtemp(j, l) = e(j, l) * CMPLX(0.D0, phase_j)
               etemp(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_zh * (mpa + mediumabs), phase_index, 8))
            CASE(4)
               ptemp(j, l) = e(j, l) * CMPLX(0.D0, phase_p)
               jtemp(j, l) = e(j, l) * CMPLX(0.D0, phase_j)
               etemp(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_zh * (mpa + mediumabs), phase_index, 8))
            END SELECT
            IF (j.NE.dim_t) THEN
               CALL calc_rho(rhotemp, mpa, e_2(j), e_2(j+1), l)
               CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
               CALL calc_delkerr(delkerr, delkerrp, e_2(j), e_2(j+1))
            ENDIF
         ENDDO
      ENDDO
  
      ! Application of the propagator (physics in omega-domain)
      SELECT CASE (switch_T)
      CASE(1)
         CONTINUE
      CASE(2)
         CALL dfftw_execute(plan_forward_erk)
         CALL dfftw_execute(plan_p)
         CALL dfftw_execute(plan_j)
         DO k = dim_r_start(num_proc), dim_r_end(num_proc)
            DO l = 1, dim_t
               etemp(l, k) = etemp(l, k) + op_t(l, k) * ptemp(l, k) + op_t_inv(l, k) * jtemp(l, k)
            ENDDO
         ENDDO
         CALL dfftw_execute(plan_backward_erk)
         etemp = diminv * etemp
      CASE(3)
         CALL dfftw_execute(plan_forward_erk)
         CALL dfftw_execute(plan_p)
         CALL dfftw_execute(plan_j)
         DO k = dim_r_start(num_proc), dim_r_end(num_proc)
            DO l = dim_th+1, dim_t
               etemp(l, k) = etemp(l, k) + op_t(l, k) * ptemp(l, k) + op_t_inv(l, k) * jtemp(l, k)
            ENDDO
         ENDDO
         CALL dfftw_execute(plan_backward_erk)
         etemp = diminv * etemp
      CASE(4)
         CALL dfftw_execute(plan_forward_erk)
         CALL dfftw_execute(plan_p)
         CALL dfftw_execute(plan_j)
         DO k = dim_r_start(num_proc), dim_r_end(num_proc)
            DO l = 1, dim_t
               etemp(l, k) = etemp(l, k) + op_t(l, k) * ptemp(l, k) + op_t_inv(l, k) * jtemp(l, k)
            ENDDO
         ENDDO
         CALL dfftw_execute(plan_backward_erk)
         etemp = diminv * etemp
      END SELECT
  
      ! Physics in time domain
      DO l = dim_r_start(num_proc), dim_r_end(num_proc)
         r = REAL(l-1) * delta_r
         e_2 = ABS(etemp(1:dim_t, l))
         e_2 = e_2 ** 2
         delkerr = 0.D0
         delkerrp = 0.D0
        IF (apply_pre_ionisation) THEN
           rhotemp = initial_electron_density(r, z, l, dim_r_start(num_proc))
        ELSE
           rhotemp = 0.D0
        ENDIF
         rhompi = 0.D0
         rho1 = 0.D0
         rho2 = 0.D0
         rhoth = 0.D0
         rhotr = 0.D0
         rhofh = 0.D0
         rhoslg2 = 0.D0
         rhoav = 0.D0
  
         rhoabstemp = 0.D0
         mediumabs = 0.D0
         CALL index_interpolation(phase_index, r)
         phase_index = phase_index * delta_z
         CALL calc_rho(rhotemp, mpa, 0.D0, e_2(1), l)
         CALL calc_absorption(rhoabstemp, mediumabs, 0.D0, e_2(1))
         DO j = 1, dim_t
            phase_p = (c3i * e_2(j) + c3d * delkerr - c5 * e_2(j)**2) * ((density_mod(l) - rhotemp * rhoat_inv) + ions_Kerr_ratio * rhotemp * rhoat_inv) * delta_z
            phase_j = -gamma2 * rhotemp * delta_z
            losses_j = -gamma1 * density_mod(l) * rhotemp * delta_z
            phase = phase_p + phase_j + phase_index
            SELECT CASE (switch_T)
            CASE(1)
               e(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_z * (mpa + mediumabs), phase, 8))
            CASE(2)
               ptemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_p)
               jtemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_j)
               e(j, l) = e(j, l) * EXP(CMPLX(losses_j - delta_z * (mpa + mediumabs), phase_index, 8))
            CASE(3)
               ptemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_p) + (hfac(j, 4) * CONJG(hfac(j, 0) * etemp(j, l)) * CMPLX(0.D0, phase_p) &
                    + 2.D0 * (hfac(j, 1) * (hfac(j, 0) * etemp(j, l))**3 + hfac(j, 2) * (hfac(j, 0) * etemp(j, l))**3 * e_2(j)) &
                    + 2.D0 * (hfac(j, 3) * (hfac(j, 0) * etemp(j, l))**5)) * ((density_mod(l) - rhotemp * rhoat_inv) + ions_Kerr_ratio * rhotemp * rhoat_inv)
               jtemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_j)
               e(j, l) = etemp(j, l) * EXP(CMPLX(losses_j - delta_z * (mpa + mediumabs), phase_index, 8))
            CASE(4)
               ptemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_p)
               jtemp(j, l) = etemp(j, l) * CMPLX(0.D0, phase_j)
               e(j, l) = etemp(j, l) * EXP(CMPLX(losses_j - delta_z * (mpa + mediumabs), phase_index, 8))
            END SELECT
            IF (j.NE.dim_t) THEN
               CALL calc_rho(rhotemp, mpa, e_2(j), e_2(j+1), l)
               CALL calc_absorption(rhoabstemp, mediumabs, e_2(j), e_2(j+1))
               CALL calc_delkerr(delkerr, delkerrp, e_2(j), e_2(j+1))
            ENDIF
         ENDDO
      ENDDO
  
      ! Save outputs
      CALL MPI_REDUCE(maxphase_part, maxphase, 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_BCAST(maxphase, 1, MPI_DOUBLE_PRECISION, 0, MPI_COMM_WORLD, ierr)
      count = count + 1
      rhomax_part = MAXVAL(rho(dim_r_start(num_proc):dim_r_end(num_proc)))
      rhoabs_max_part = MAXVAL(rhoabs(dim_r_start(num_proc):dim_r_end(num_proc)))
      CALL MPI_REDUCE(peakmax_part, peakmax(count), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(rhomax_part, rhomax(count), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(rhoabs_max_part, rhoabs_max(count), 1, MPI_DOUBLE_PRECISION, MPI_MAX, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(energy_part, energy(count), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
      CALL MPI_REDUCE(energy_fil_part, energy_fil(count), 1, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
  
      IF (my_rank.EQ.0) z_buff(count) = four_z_Rayleigh * z
      IF (count.GE.rhodist .OR. z.LE.delta_z) THEN
        ALLOCATE(onax_t_data(1,1:dim_t))
        ALLOCATE(fluence_data(1:dim_r/num_proc))
        ALLOCATE(plasma_channel_data(1:dim_r/num_proc))
        ALLOCATE(losses_plasma_data(1:dim_r/num_proc))
        ALLOCATE(losses_ionization_data(1:dim_r/num_proc))
      ENDIF
  
      IF (count.GE.rhodist) THEN
         e_2 = 0.D0
         DO l = dim_r_start(num_proc), dim_r_end(num_proc)
            e_2 = e_2 + ABS(e(1:dim_t, l))**2 * REAL(l-1, 8)
         ENDDO
         CALL MPI_REDUCE(e_2(1:dim_t), e_2KK(1:dim_t), dim_t, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         IF (my_rank.EQ.0) THEN
            CALL h5open_f(error)
            CALL h5fopen_f(main_h5_fname, H5F_ACC_RDWR_F, file_id, error)
            CALL h5lexists_f(file_id, longstep_grpname, group_status, error)
            IF (.NOT. group_status) THEN
              CALL h5gcreate_f(file_id, longstep_grpname, group_id, error)
              CALL h5gclose_f(group_id, error)
            ENDIF
            IF (longstep_write_count.EQ.0) THEN
              CALL create_1D_dset_unlimited(file_id, z_buff_dset_name, REAL(z_buff(1:rhodist), 4), rhodist)
              CALL create_1D_dset_unlimited(file_id, rhoabs_max_dset_name, REAL(rhoabs_max(1:rhodist), 4), rhodist)
              CALL create_1D_dset_unlimited(file_id, TRIM(rhoabs_max_dset_name) // "_normalised", REAL(rhoabs_max(1:rhodist), 4), rhodist)
              CALL create_1D_dset_unlimited(file_id, peakmax_dset_name, REAL(peakmax(1:rhodist), 4), rhodist)
              peakmax(1:rhodist) = critical_power / (1000.D0 * 4.D0 * PI * beam_waist**2) * peakmax(1:rhodist) * peakmax(1:rhodist)
              CALL create_1D_dset_unlimited(file_id, TRIM(peakmax_dset_name) // "_normalised", REAL(peakmax(1:rhodist), 4), rhodist)
              CALL create_1D_dset_unlimited(file_id, rhomax_dset_name, REAL(rhomax(1:rhodist), 4), rhodist)
              energy(1:rhodist) = 1.D-6 * critical_power * 1.D9 / (4.D0 * PI) * pulse_duration * 1.D-15 * 1.D6 * energy(1:rhodist)
              CALL create_1D_dset_unlimited(file_id, energy_dset_name, REAL(2.D0 * PI * energy(1:rhodist) * delta_t * delta_r**2, 4), rhodist)
              CALL h5_add_units_1D(file_id, energy_dset_name, '[J]')
              energy_fil(1:rhodist) = 1.D-6 * critical_power * 1.D9 / (4.D0 * PI) * pulse_duration * 1.D-15 * 1.D6 * energy_fil(1:rhodist)
              CALL create_1D_dset_unlimited(file_id, energy_fil_dset_name, REAL(2.D0 * PI * energy_fil(1:rhodist) * delta_t * delta_r**2, 4), rhodist)
              CALL h5_add_units_1D(file_id, energy_fil_dset_name, '[J]')
              DO k1 = 1, dim_t
                onax_t_data(1, k1) = REAL(efield_factor * REAL(efield_osc(k1) * e(k1, 1)), 4)
              ENDDO
              IF (dset_write_count.EQ.0) THEN
                CALL create_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)/), 1)
                CALL h5_add_units_1D(file_id, powmax_dset_name, '[C.U.]')
                CALL create_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, dim_t)
                CALL h5_add_units_1D(file_id, onax_t_dset_name, '[V/m]')
                CALL create_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh * z, 4)/), 1)
                CALL h5_add_units_1D(file_id, every_rhodist_z_dset_name, '[m]')
              ELSE
                CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)/), &
                new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
                CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
                  new_dims = (/INT(dset_write_count + 1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                  memspace_dims = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                  offset = (/INT(dset_write_count, HSIZE_T), INT(0, HSIZE_T)/), & 
                  hyperslab_size = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/))
                CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh * z, 4)/), &
                  new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                  offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
              ENDIF
              dset_write_count = dset_write_count + 1
              longstep_write_count = longstep_write_count + 1
              original_rhodist = rhodist
            ELSE
              IF (rhodist.NE.original_rhodist) THEN
                new_dims = (/INT(longstep_write_count * original_rhodist + rhodist, HSIZE_T)/)
              ELSE
                new_dims = (/INT((longstep_write_count + 1) * rhodist, HSIZE_T)/)
              ENDIF
              memspace_dims = (/INT(rhodist, HSIZE_T)/)
              offset = (/INT(longstep_write_count * original_rhodist, HSIZE_T)/)
              hyperslab_size = (/INT(rhodist, HSIZE_T)/)
              IF (count.NE.0) THEN
                CALL extend_1D_dset_unlimited(file_id, z_buff_dset_name, & 
                  REAL(z_buff(1:count), 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                CALL extend_1D_dset_unlimited(file_id, rhoabs_max_dset_name, & 
                  REAL(rhoabs_max(1:count), 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                CALL extend_1D_dset_unlimited(file_id, peakmax_dset_name, & 
                  REAL(peakmax(1:count), 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                peakmax(1:rhodist) = critical_power / (1000.D0 * 4.D0 * PI * beam_waist**2) * peakmax(1:rhodist) * peakmax(1:rhodist)
                CALL extend_1D_dset_unlimited(file_id, TRIM(peakmax_dset_name) // "_normalised", & 
                  REAL(peakmax(1:count), 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                CALL extend_1D_dset_unlimited(file_id, rhomax_dset_name, & 
                  REAL(rhomax(1:count), 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                energy(1:rhodist) = 1.D-6 * critical_power * 1.D9 / (4.D0 * PI) * pulse_duration * 1.D-15 * 1.D6 * energy(1:rhodist)
                CALL extend_1D_dset_unlimited(file_id, energy_dset_name, & 
                  REAL(2.D0 * PI * energy(1:count) * delta_t * delta_r**2, 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
                energy_fil(1:rhodist) = 1.D-6 * critical_power * 1.D9 / (4.D0 * PI) * pulse_duration * 1.D-15 * 1.D6 * energy_fil(1:rhodist)
                CALL extend_1D_dset_unlimited(file_id, energy_fil_dset_name, & 
                  REAL(2.D0 * PI * energy_fil(1:count) * delta_t * delta_r**2, 4), new_dims, &
                  memspace_dims, offset, hyperslab_size)
              ENDIF
              DO k1 = 1, dim_t
                onax_t_data(1, k1) = REAL(efield_factor * REAL(efield_osc(k1) * e(k1, 1)), 4)
              ENDDO
              CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)/), &
                new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
              CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
                new_dims = (/INT(dset_write_count + 1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                memspace_dims = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                offset = (/INT(dset_write_count, HSIZE_T), INT(0, HSIZE_T)/), & 
                hyperslab_size = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/))
              CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh * z, 4)/), &
                new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
              dset_write_count = dset_write_count + 1
              longstep_write_count = longstep_write_count + 1
            ENDIF
            CALL h5fclose_f(file_id, error)
            CALL h5close_f(error)
         ENDIF
         fluence_data(:) = REAL(fluence * delta_t, 4)
         plasma_channel_data(:) = REAL(rho, 4)
         losses_plasma_data(:) = REAL(losses_plasma * delta_t, 4)
         losses_ionization_data(:) = REAL(losses_ionization * delta_t, 4)
         ptr_f => fluence_data
         ptr_p => plasma_channel_data
         ptr_lp => losses_plasma_data
         ptr_li => losses_ionization_data
         IF (length_of_linked_list .EQ. 0) THEN
           CALL list_init(fluence_ll, DATA = TRANSFER(ptr_f, list_data))
           CALL list_init(plasma_channel_ll, DATA = TRANSFER(ptr_p, list_data))
           CALL list_init(losses_plasma_ll, DATA = TRANSFER(ptr_lp, list_data))
           CALL list_init(losses_ionization_ll, DATA = TRANSFER(ptr_li, list_data))
         ELSE
           CALL list_append(fluence_ll, DATA = TRANSFER(ptr_f, list_data))
           CALL list_append(plasma_channel_ll, DATA = TRANSFER(ptr_p, list_data))
           CALL list_append(losses_plasma_ll, DATA = TRANSFER(ptr_lp, list_data))
           CALL list_append(losses_ionization_ll, DATA = TRANSFER(ptr_li, list_data))
         ENDIF
         NULLIFY(ptr_f, ptr_p, ptr_lp, ptr_li)
         length_of_linked_list = length_of_linked_list + 1
         count = 0
      ELSE IF (z.LE.delta_z) THEN
         fluence_data(:) = REAL(fluence * delta_t, 4)
         plasma_channel_data(:) = REAL(rho, 4)
         losses_plasma_data(:) = REAL(losses_plasma * delta_t, 4)
         losses_ionization_data(:) = REAL(losses_ionization * delta_t, 4)
         ptr_f => fluence_data
         ptr_p => plasma_channel_data
         ptr_lp => losses_plasma_data
         ptr_li => losses_ionization_data
         IF (length_of_linked_list .EQ. 0) THEN
           CALL list_init(fluence_ll, DATA = TRANSFER(ptr_f, list_data))
           CALL list_init(plasma_channel_ll, DATA = TRANSFER(ptr_p, list_data))
           CALL list_init(losses_plasma_ll, DATA = TRANSFER(ptr_lp, list_data))
           CALL list_init(losses_ionization_ll, DATA = TRANSFER(ptr_li, list_data))
         ELSE
           CALL list_append(fluence_ll, DATA = TRANSFER(ptr_f, list_data))
           CALL list_append(plasma_channel_ll, DATA = TRANSFER(ptr_p, list_data))
           CALL list_append(losses_plasma_ll, DATA = TRANSFER(ptr_lp, list_data))
           CALL list_append(losses_ionization_ll, DATA = TRANSFER(ptr_li, list_data))
         ENDIF
         NULLIFY(ptr_f, ptr_p, ptr_lp, ptr_li)
         length_of_linked_list = length_of_linked_list + 1
         e_2 = 0.D0
         DO l = dim_r_start(num_proc), dim_r_end(num_proc)
            e_2 = e_2 + ABS(e(1:dim_t, l))**2 * REAL(l-1, 8)
         ENDDO
         CALL MPI_REDUCE(e_2(1:dim_t), e_2KK(1:dim_t), dim_t, MPI_DOUBLE_PRECISION, MPI_SUM, 0, MPI_COMM_WORLD, ierr)
         IF (my_rank.EQ.0) THEN
            CALL h5open_f(error)
            CALL h5fopen_f(main_h5_fname, H5F_ACC_RDWR_F, file_id, error)
            CALL h5lexists_f(file_id, longstep_grpname, group_status, error)
            IF (.NOT. group_status) THEN
              CALL h5gcreate_f(file_id, longstep_grpname, group_id, error)
              CALL h5gclose_f(group_id, error)
            ENDIF
            powmax_data(1, 1) = REAL(z, 4)
            powmax_data(1, 2) = REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)
            DO k1 = 1, dim_t
              onax_t_data(1, k1) = REAL(efield_factor * REAL(efield_osc(k1) * e(k1, 1)), 4)
            ENDDO
            IF (dset_write_count.EQ.0) THEN
              CALL create_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)/), 1)
              CALL create_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, dim_t)
              CALL create_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh * z, 4)/), 1)
              CALL h5_add_units_1D(file_id, powmax_dset_name, '[C.U.]2')
              CALL h5_add_units_1D(file_id, onax_t_dset_name, '[V/m]?2')
              CALL h5_add_units_1D(file_id, every_rhodist_z_dset_name, '[m]2')
            ELSE
              CALL extend_1D_dset_unlimited(file_id, powmax_dset_name, (/REAL(2.D0 * PI * MAXVAL(e_2KK) * delta_r**2, 4)/), &
                new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
              CALL extend_2D_dset_unlimited(file_id, onax_t_dset_name, onax_t_data, & 
                new_dims = (/INT(dset_write_count + 1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                memspace_dims = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/), & 
                offset = (/INT(dset_write_count, HSIZE_T), INT(0, HSIZE_T)/), & 
                hyperslab_size = (/INT(1, HSIZE_T), INT(dim_t, HSIZE_T)/))
              CALL extend_1D_dset_unlimited(file_id, every_rhodist_z_dset_name, (/REAL(four_z_Rayleigh * z, 4)/), &
                new_dims = (/INT(dset_write_count + 1, HSIZE_T)/), memspace_dims = (/INT(1, HSIZE_T)/), &
                offset = (/INT(dset_write_count, HSIZE_T)/), hyperslab_size = (/INT(1, HSIZE_T)/))
            ENDIF
            dset_write_count = dset_write_count + 1
            CALL h5fclose_f(file_id, error)
            CALL h5close_f(error)
         ENDIF
         DEALLOCATE(fluence_data, plasma_channel_data, losses_plasma_data, losses_ionization_data)
      ENDIF
      RETURN
    END SUBROUTINE mult_phase
  
    ! \brief Subroutine to drive the single-step calculation
    SUBROUTINE propagation
      USE fields
      USE fft
      USE mpi_stuff
  
      IMPLICIT NONE
      
      CALL fft_forward_inplace(.FALSE.)
      CALL mult_propagator
      CALL fft_backward2_inplace
      z = z + delta_zh
      CALL mult_phase
      CALL fft_forward_inplace(.TRUE.)
      CALL mult_propagator
      CALL fft_backward2_inplace
      CALL absorbation
      z = z + delta_zh
  
      RETURN
    END SUBROUTINE propagation
  
  END MODULE long_step
  