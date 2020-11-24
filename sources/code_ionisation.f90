! There are various models of ionisation implemented, these modules are accessed
! from the main code for both computing and then evaluating the ionisation rates.

! It is contributed by more authors using ionisation of their needs.
! Some of them are: Rachel Nuter (original PPT-procedure)
!                   Jan Vabek (loading from an external table)

!=========================================================================

!                 MODULE FOR THE MEDIA PARAMETERS

!==========================================================================

MODULE MEDIUM_PARAMETER

  IMPLICIT NONE

  DOUBLE PRECISION, SAVE  :: ionisation_potential
  DOUBLE PRECISION, SAVE  :: residue_charge
  DOUBLE PRECISION, SAVE  :: n0_indice
  DOUBLE PRECISION, SAVE  :: critical_density
  DOUBLE PRECISION, SAVE  :: atomic_density
  DOUBLE PRECISION, SAVE  :: reduced_mass
  INTEGER, SAVE  :: angular_momentum

END MODULE MEDIUM_PARAMETER

! MODULE MEDIUM_PARAMETER_N2

!   IMPLICIT NONE

!   DOUBLE PRECISION, SAVE  :: ionisation_potential
!   DOUBLE PRECISION, SAVE  :: residue_charge
!   DOUBLE PRECISION, SAVE  :: n0_indice
!   DOUBLE PRECISION, SAVE  :: critical_density
!   DOUBLE PRECISION, SAVE  :: atomic_density
!   DOUBLE PRECISION, SAVE  :: reduced_mass
!   INTEGER, SAVE  :: angular_momentum

! END MODULE MEDIUM_PARAMETER_N2

!==========================================================================

!                      MODULE FOR THE LASER FIELD

!==========================================================================
MODULE LASER_PARAMETER

  ! This module is used to read and save the Laser Parameters

  IMPLICIT NONE

  DOUBLE PRECISION, SAVE     :: photon_energy
  DOUBLE PRECISION, SAVE     :: pulse_duration
  DOUBLE PRECISION, SAVE     :: critical_power
  DOUBLE PRECISION, SAVE     :: beam_waist

END MODULE LASER_PARAMETER

!==========================================================================

!                      MODULE FOR ATOMS - PPT METHOD

!==========================================================================

MODULE PPT_PARAMETER

  ! All data are evaluated in Atomic Unit
  ! Modified 19/04/04

  USE MEDIUM_PARAMETER

  IMPLICIT NONE

  DOUBLE PRECISION, SAVE :: n_number                         
  DOUBLE PRECISION, SAVE :: coulomb_field
  DOUBLE PRECISION, SAVE :: C_factor
  DOUBLE PRECISION, SAVE :: C_ADK_factor
  DOUBLE PRECISION, SAVE :: f_factor


CONTAINS

  SUBROUTINE ATTRIBUTE_DATA_PPT                 ! This routine calculates the parameters

    IMPLICIT NONE

    DOUBLE PRECISION                :: n_star, l_star

    INTERFACE
       FUNCTION gamm(xx)
         DOUBLE PRECISION    :: xx
         DOUBLE PRECISION                :: gamm
       END FUNCTION gamm
    END INTERFACE


    n_number             = residue_charge / sqrt(2.d0 * ionisation_potential)
    coulomb_field        = (2.d0 * ionisation_potential)**(1.5d0)

    n_star = n_number
    l_star = n_star - 1.d0
    C_factor = 2.d0**(2.d0*n_star) / (n_star * gamm(n_star+l_star+1.d0)* &
         gamm(n_star - l_star) )
    SELECT CASE(angular_momentum)
    CASE(0) 
      f_factor = 1.d0
    CASE(1)
      f_factor = 3.d0
    CASE(2)
      f_factor = 5.D0
    END SELECT
     
         

  END SUBROUTINE ATTRIBUTE_DATA_PPT

!---------------------------------------------------------------------------
  SUBROUTINE ATTRIBUTE_DATA_ADK                 ! This routine calculates the parameters

    IMPLICIT NONE

    DOUBLE PRECISION :: C_l_2, C_l_4

    C_l_2   = 0.683d0
    C_l_4   = 0.033d0

    C_ADK_factor = 0.5 * ( C_l_2 * sqrt(30.d0) + C_l_4 * sqrt(180.d0) )**2


  END SUBROUTINE ATTRIBUTE_DATA_ADK


END MODULE PPT_PARAMETER

! MODULE PPT_PARAMETER_N2

!   ! All data are evaluated in Atomic Unit
!   ! Modified on the 19/04/04

!   USE MEDIUM_PARAMETER_N2

!   IMPLICIT NONE

!   DOUBLE PRECISION, SAVE :: n_number                         
!   DOUBLE PRECISION, SAVE :: coulomb_field
!   DOUBLE PRECISION, SAVE :: C_factor
!   DOUBLE PRECISION, SAVE :: C_ADK_factor
!   DOUBLE PRECISION, SAVE :: f_factor


! CONTAINS

!   SUBROUTINE ATTRIBUTE_DATA_PPT                 ! This routine calculates the parameters

!     IMPLICIT NONE

!     DOUBLE PRECISION                :: n_star, l_star

!     INTERFACE
!        FUNCTION gamm(xx)
!          DOUBLE PRECISION   :: xx
!          DOUBLE PRECISION                :: gamm
!        END FUNCTION gamm
!     END INTERFACE


!     n_number             = residue_charge / sqrt(2.d0 * ionisation_potential)
!     coulomb_field        = (2.d0 * ionisation_potential)**(1.5d0)

!     n_star = n_number
!     l_star = n_star - 1.d0
!     C_factor = 2.d0**(2.d0*n_star) / (n_star * gamm(n_star+l_star+1.d0)* &
!          gamm(n_star - l_star) )
!     SELECT CASE(angular_momentum)
!     CASE(0) 
!       f_factor = 1.d0
!     CASE(1)
!       f_factor = 3.d0
!     CASE(2)
!       f_factor = 5.D0
!     END SELECT
     
         

!   END SUBROUTINE ATTRIBUTE_DATA_PPT

! END MODULE PPT_PARAMETER_N2

!===================================================================================================

MODULE PPT

  ! This module is used for PPT. It does the following procedure
  ! 1- Read the input parameters
  ! 2- Find the intensity range to build the table
  ! 3- Create the table
  ! $- A routine to serach the intensity

  USE MEDIUM_PARAMETER
  USE LASER_PARAMETER
  USE PPT_PARAMETER

  IMPLICIT NONE
  DOUBLE PRECISION, SAVE   :: INTENSITY_STEP,intensity_step_inv
  INTEGER, SAVE            :: DIMENSION_PPT
  INTEGER, SAVE            :: JLOTI, JLOTIP1                          ! Integer for Intensity Index
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPT_TABLE 



CONTAINS

  !----------------------------------------
  !1- Read the input parameters

  SUBROUTINE INITIALISE_PPT(THEORY)

    IMPLICIT NONE
    CHARACTER(LEN=3), INTENT(IN)  :: THEORY
    
    IF ( THEORY == "PPT") CALL ATTRIBUTE_DATA_PPT
    IF ( THEORY == "ADK") THEN
       CALL ATTRIBUTE_DATA_PPT
       CALL ATTRIBUTE_DATA_ADK
    ENDIF   
    JLOTI=1
    JLOTIP1=1


  END SUBROUTINE INITIALISE_PPT

  !---------------------------------------
  ! 2- Find the intensity range to build the table

  SUBROUTINE FIND_INTENSITY_AREA(THEORY)

    ! This function solves the equation n2*I = rho / 2*rho_c
    ! The intensities are done in W/cm2

    IMPLICIT NONE
    CHARACTER(LEN=3), INTENT(IN)  :: THEORY
    DOUBLE PRECISION, PARAMETER   :: PI = 3.14159265d0
    DOUBLE PRECISION, PARAMETER   :: epsilon = 0.1d0
    DOUBLE PRECISION              :: INTENSITY_MAX
    DOUBLE PRECISION              :: factor
    DOUBLE PRECISION              :: intensity
    DOUBLE PRECISION              :: ionisation_rate
    DOUBLE PRECISION              :: delta

    INTERFACE
     FUNCTION IONISATION_RATE_PPT(intensity)
         DOUBLE PRECISION, INTENT(IN) :: intensity
         DOUBLE PRECISION             :: IONISATION_RATE_PPT
     END FUNCTION IONISATION_RATE_PPT
    !  FUNCTION IONISATION_RATE_ADK(intensity)
    !    DOUBLE PRECISION, INTENT(IN) :: intensity
    !    DOUBLE PRECISION             :: IONISATION_RATE_ADK
    !  END FUNCTION IONISATION_RATE_ADK
    !  FUNCTION IONISATION_RATE_CRISTAL(intensity)
    !    DOUBLE PRECISION, INTENT(IN) :: intensity
    !    DOUBLE PRECISION             :: IONISATION_RATE_CRISTAL
    !  END FUNCTION
    END INTERFACE


    IF ( (THEORY == "PPT").OR.(THEORY == "ADK") ) factor = ( PI * 1.0d25 /( 45.5635d0**2 * 2.41889)) * critical_power * photon_energy**2 *&
         ( atomic_density / critical_density ) *  MIN(pulse_duration,1000.d0)

    IF (THEORY == "IRC") factor = ( PI * 1.0d52 /( 45.5635d0**2 * 2.41889 * 5.29177**3)) * critical_power * photon_energy**2 *&
         ( MIN(pulse_duration,1000.d0) / critical_density ) 


    intensity = 1.d10
    delta = 1.d0

    DO WHILE( (delta > epsilon).AND.(intensity <= 1.d16) )
       intensity = intensity + 10.d0**(Int(log10(intensity))) /50.0
     IF ( THEORY == "PPT") THEN
        ionisation_rate = IONISATION_RATE_PPT(intensity)
    !  ELSE IF( THEORY == "ADK" ) THEN
    !     ionisation_rate = IONISATION_RATE_ADK(intensity)
    !  ELSE IF( THEORY == "IRC") THEN
    !     ionisation_rate = IONISATION_RATE_CRISTAL(intensity)
     ELSE
        print *, 'The THEORY FLAG is wrong, choose ADK, PPT or IRC'
        print *, 'The code is stopped'
        STOP
     ENDIF
       delta = ABS( (intensity - ionisation_rate * factor) / (ionisation_rate * factor)  )
    ENDDO

    INTENSITY_STEP = intensity / 5000.d0
    INTENSITY_MAX  = intensity * 100.d0
    DIMENSION_PPT  = Int(intensity_max / intensity_step) +1

  END SUBROUTINE FIND_INTENSITY_AREA

  !----------------------------------------
  ! 3- Create the table

  SUBROUTINE FILL_TABLE(THEORY)
    USE mpi_stuff
    USE HDF5
    USE HDF5_helper
    IMPLICIT NONE
    CHARACTER(len = 3),INTENT(IN):: THEORY
    DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
    DOUBLE PRECISION, PARAMETER :: field_intensity_au = 3.50944758d16 !unit of the field intensity in atomic units
    DOUBLE PRECISION            :: intensity_factor
    DOUBLE PRECISION            :: rate_factor
    DOUBLE PRECISION            :: MPA_factor
    DOUBLE PRECISION            :: intensity
    DOUBLE PRECISION            :: efield
    DOUBLE PRECISION            :: ionisation_rate
    DOUBLE PRECISION            :: o_atom_dens, o_crit_dens, o_beam_waist, o_photenergy, o_pulse_duration, o_n0, o_ionpot
    INTEGER                     :: i,error
    INTEGER(HID_T)              :: file_id, group_id
    CHARACTER(LEN=25)           :: filename = "PPT_tables.h5", groupname = "PPT"
    REAL(8), ALLOCATABLE        :: rates_table(:,:), reference_table(:,:)
    LOGICAL                     :: file_exists, just_read = .FALSE.
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: Egrid, ionisation_rates
    CHARACTER(LEN=25)           :: outfilename = "results.h5", outgroupname="ionisation_model"

    INTERFACE
     FUNCTION IONISATION_RATE_PPT(intensity)
         DOUBLE PRECISION, INTENT(IN) :: intensity
         DOUBLE PRECISION             :: IONISATION_RATE_PPT
     END FUNCTION IONISATION_RATE_PPT
    !  FUNCTION IONISATION_RATE_ADK(intensity)
    !      DOUBLE PRECISION, INTENT(IN) :: intensity
    !      DOUBLE PRECISION             :: IONISATION_RATE_ADK
    !  END FUNCTION IONISATION_RATE_ADK
    !  FUNCTION IONISATION_RATE_CRISTAL(intensity)
    !    DOUBLE PRECISION, INTENT(IN) :: intensity
    !    DOUBLE PRECISION             :: IONISATION_RATE_CRISTAL
    !  END FUNCTION
    END INTERFACE

    IF (THEORY == "PPT") THEN
    PRINT*, 'PPT-fill-table is acessed, proc', my_rank
    ENDIF


    ALLOCATE(PPT_TABLE(DIMENSION_PPT, 3))
    ! Normalised factors
    intensity_factor = 4.d0 * PI * beam_waist**2 * 1.d-9 / critical_power
 
    IF ( (THEORY == "PPT") .OR. (THEORY == "ADK") ) THEN
    rate_factor      =(4.d0 * PI**2 * 1.d16 / (45.5635**2 * 2.41889)) * &
         (atomic_density / critical_density) *&
         (photon_energy**2) * (beam_waist**2) * pulse_duration
    MPA_factor       = n0_indice * critical_density * ionisation_potential * 45.5635 * 4.359d-10 / &
         (photon_energy * pulse_duration * 2.d0 * PI)
    ELSE IF ( THEORY == "IRC") THEN
    rate_factor      =(4.d0 * PI**2 * 1.d43 / (45.5635**2 * 2.41889 * 5.29177**3)) * &
         (pulse_duration / critical_density) * photon_energy**2 * beam_waist**2
    MPA_factor       = (2.d0 * PI * 4.359d33 / (45.5635 * 2.41889 * 5.29177**3)) * &
         n0_indice * beam_waist**2 * photon_energy * ionisation_potential 
    ENDIF


    ! Fill the Table
    PPT_TABLE(1, 1) = 0.d0                   
    PPT_TABLE(1, 2) = 0.d0
    PPT_TABLE(1, 3) = 0.d0

    IF (THEORY == "PPT") THEN

        ALLOCATE(Egrid(DIMENSION_PPT),ionisation_rates(DIMENSION_PPT))
        Egrid(1) = 0.0d0; ionisation_rates(1) = 0.0d0 ! 0 imposed for 0-field

        DO i = 2, dimension_PPT
         intensity = (i-1) * intensity_step
         ionisation_rate = ionisation_rate_PPT(intensity)
         Egrid(i) = sqrt(intensity/field_intensity_au)
         ionisation_rates(i) = ionisation_rate
         PPT_TABLE(i, 1) = intensity * intensity_factor                    ! normalised intensity (C.U.)
         PPT_TABLE(i, 2) = ionisation_rate * rate_factor                   ! ionisation rate (C.U.)
         PPT_TABLE(i, 3) = MPA_factor * ( ionisation_rate  * rate_factor/ intensity )    ! Normalised MPA
        ENDDO

        
        ! The ionisation table is provided in the outputs
        IF (my_rank.EQ.0) THEN 
          CALL h5open_f(error)
          CALL h5fopen_f(outfilename, H5F_ACC_RDWR_F, file_id, error)
          CALL h5gcreate_f(file_id, outgroupname, group_id, error)
          CALL create_dset(group_id, 'Egrid', Egrid, DIMENSION_PPT)
          CALL create_dset(group_id, 'ionisation_rates', ionisation_rates, DIMENSION_PPT)
          CALL h5_add_units_1D(group_id, 'Egrid', '[a.u.]')
          CALL h5_add_units_1D(group_id, 'ionisation_rates', '[a.u.]')
          CALL h5gclose_f(group_id, error)
          CALL h5fclose_f(file_id, error)
        ENDIF

	DEALLOCATE(Egrid,ionisation_rates)


    ! these models don't print reference, should we keep them or take as obsolete non-canonical cases?
    ! ELSE IF (THEORY == "ADK") THEN
    ! DO i = 2, dimension_PPT
    !    intensity = (i-1) * intensity_step
    !    ionisation_rate = ionisation_rate_ADK(intensity)
    !    PPT_TABLE(i, 1) = intensity * intensity_factor                    ! Normalised Intensity
    !    PPT_TABLE(i, 2) = ionisation_rate * rate_factor                   ! Gamma
    !    PPT_TABLE(i, 3) = MPA_factor * ( ionisation_rate  * rate_factor/ intensity )    ! Normalised MPA
    ! ENDDO
    ! ELSE IF (THEORY == "IRC") THEN
    ! DO i = 2, dimension_PPT
    !    intensity = (i-1) * intensity_step
    !    ionisation_rate = ionisation_rate_cristal(intensity)
    !    PPT_TABLE(i, 1) = intensity * intensity_factor                    ! Normalised Intensity
    !    PPT_TABLE(i, 2) = ionisation_rate * rate_factor                   ! Normalised Gamma
    !    PPT_TABLE(i, 3) = MPA_factor * ionisation_rate  / intensity       ! Normalised MPA
    ! ENDDO    
    ENDIF

    IF(my_rank.EQ.0) THEN
       OPEN (UNIT = 2, FILE = 'ionisation_table.dat', STATUS = 'unknown')
       DO i = 1, dimension_PPT
          WRITE(2, '(3(2x, e12.5))') PPT_TABLE(i, 1),PPT_TABLE(i, 2),PPT_TABLE(i, 3)
       ENDDO
       CLOSE(2); close(4); close(7);
    ENDIF



    PPT_TABLE(dimension_PPT, 1)=PPT_TABLE(dimension_PPT, 1)*1.d99

    OPEN(UNIT=2,FORM='unformatted',FILE='table_binary.dat')
    WRITE(2) PPT_TABLE
    CLOSE(2)
    intensity_step_inv = 1.d0 /(intensity_step * intensity_factor)

  END SUBROUTINE FILL_TABLE

  !----------------------------------------
  ! 4- Find the intensitty and Interpolate

  SUBROUTINE INTERPOLATE_PPT(var1,mpa,intensity)

    IMPLICIT NONE

    REAL(8) var1,mpa,intensity
    INTEGER(4) jmin,jmax

    REAL(8) intensity_step

    INTERFACE
       SUBROUTINE HUNT(xx, n, x, jlo)
         INTEGER                       :: n
         DOUBLE PRECISION, DIMENSION(n) :: XX
         DOUBLE PRECISION             :: x
         INTEGER                      :: jlo
       END SUBROUTINE HUNT

    END INTERFACE


    CALL hunt(ppt_table,dimension_ppt,intensity,jlotip1)
    IF (jlotip1.GE.DIMENSION_PPT-1) THEN 
!       print *, "Ionisation table exceeded"
!       print *, intensity
!       READ(5,*)
       jlotip1=DIMENSION_PPT-1
    ENDIF
    IF (jlotip1.LT.1) jlotip1=1
    intensity_step = ppt_table(jlotip1+1,1) - ppt_table(jlotip1,1)
    mpa=        ((intensity-ppt_table(jlotip1,1))*ppt_table(jlotip1+1,3)+(ppt_table(jlotip1+1,1)-intensity)*ppt_table(jlotip1,3))/intensity_step
!    mpa=        ((intensity-ppt_table(jlotip1,1))*ppt_table(jlotip1+1,3)+(ppt_table(jlotip1+1,1)-intensity)*ppt_table(jlotip1,3))*intensity_step_inv
    jmin=MIN(jloti,jlotip1)
    jmax=MAX(jloti,jlotip1)
    var1=SUM(ppt_table(jmin:jmax+1,2),DIM=1)/REAL(jmax-jmin+2,8)
    jloti=jlotip1

    RETURN
  END SUBROUTINE INTERPOLATE_PPT


END MODULE PPT


!=================================================================================================== Load from external table

module libraries
IMPLICIT NONE;
CONTAINS

  integer FUNCTION filelength(filename); 
    character(*), intent(in) :: filename;

    integer :: IO, N_load;

    open(UNIT=1,FILE=filename,FORM="FORMATTED",action='read');
    N_load = 0; IO = 0;
    do while(IO >= 0)
      read(1,*,IOSTAT = IO); N_load = N_load+1; 
    enddo
    close(1);
    N_load = N_load - 1; !# of lines in the loaded file
    filelength = N_load;
  end function filelength; 
end module libraries

MODULE External_ionisation_table

  ! This module loads an external ionisation table from pre-calculated HDF5 archive (or from an ASCII file-testing case). It also provides the rate by interpolation lookup.
   
  ! routines
  ! RESCALE_TABLE_EXT: it loads the table and transform it into C.U.
  ! INTERPOLATE_EXT: table lookup

  USE MEDIUM_PARAMETER
  USE LASER_PARAMETER
  USE PPT_PARAMETER

  IMPLICIT NONE
  DOUBLE PRECISION, SAVE   :: INTENSITY_STEP,intensity_step_inv
  INTEGER, SAVE            :: DIMENSION_EXT
  INTEGER, SAVE            :: JLOTI, JLOTIP1                          ! Integer for Intensity Index
  DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: EXT_TABLE 

CONTAINS
  SUBROUTINE RESCALE_TABLE_EXT
  ! It takes the ionsation table from the pre-computed HDF5 archive (default option).
  ! There is also possibility to load it from ASCII-formated 'rates_atomic.dat' file, this is mainly for testing. 
  ! The ionization is stored as follows:
  ! || electric field ||  ionisation rate  ||, both in atomic units, corresponding datasets are:
  ! ||     Egrid      ||  ionisation_rates ||, stored in the 'ionisation_model' group. ASCII file contains thesw two columns
  ! It is required to include 0 (the first row is omitted).
    
    USE mpi_stuff
    USE libraries
    USE HDF5
    USE HDF5_helper
    IMPLICIT NONE
    DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
    DOUBLE PRECISION, PARAMETER :: field_intensity_au = 3.50944758d16 !field intensity in atomic units
    DOUBLE PRECISION            :: intensity_factor
    DOUBLE PRECISION            :: rate_factor
    DOUBLE PRECISION            :: MPA_factor
    DOUBLE PRECISION            :: intensity
    DOUBLE PRECISION            :: ionisation_rate
    INTEGER                     :: i
    DOUBLE PRECISION, DIMENSION(:), ALLOCATABLE :: dumvect, Egrid, ionisation_rates
    LOGICAL                     :: file_exists
    CHARACTER(LEN=25)           :: filename = "calculated_tables.h5", groupname = "ionisation_model"
    CHARACTER(LEN=25)           :: outfilename = "results.h5", outgroupname="ionisation_model"
    INTEGER                     :: error
    INTEGER(HID_T)              :: file_id, group_id
    REAL(8), ALLOCATABLE        :: rates_atomic(:,:)
    INTEGER(HSIZE_T), DIMENSION(2) :: dims_ext
    !INTEGER(HSIZE_T), DIMENSION(1) :: h5dim1
    


    ! INTERFACE
    !  FUNCTION IONISATION_RATE_ADK(intensity)
    !      DOUBLE PRECISION, INTENT(IN) :: intensity
    !      DOUBLE PRECISION             :: IONISATION_RATE_ADK
    !  END FUNCTION IONISATION_RATE_ADK
    !  FUNCTION IONISATION_RATE_CRISTAL(intensity)
    !    DOUBLE PRECISION, INTENT(IN) :: intensity
    !    DOUBLE PRECISION             :: IONISATION_RATE_CRISTAL
    !  END FUNCTION
    ! END INTERFACE

    ! Factors to rescale in computational units
    intensity_factor = 4.d0 * PI * beam_waist**2 * 1.d-9 / critical_power    
    rate_factor      =(4.d0 * PI**2 * 1.d16 / (45.5635**2 * 2.41889)) * &
         (atomic_density / critical_density) *&
         (photon_energy**2) * (beam_waist**2) * pulse_duration
    MPA_factor       = n0_indice * critical_density * ionisation_potential * 45.5635 * 4.359d-10 / &
         (photon_energy * pulse_duration * 2.d0 * PI)

    ! default option if HDF5-archive exists
    INQUIRE(FILE="calculated_tables.h5", EXIST=file_exists) ! there are only these tables within now, if extended, use h5lexists
    IF (file_exists) THEN
        CALL h5open_f(error)
        CALL h5fopen_f(filename, H5F_ACC_RDONLY_F, file_id, error) ! all workers do
        CALL h5gopen_f(file_id, groupname, group_id, error)
        CALL ask_for_size_1D(group_id, 'Egrid', DIMENSION_EXT)

        ALLOCATE(EXT_TABLE(DIMENSION_EXT, 3))
        ALLOCATE(Egrid(DIMENSION_EXT),ionisation_rates(DIMENSION_EXT))
        
        CALL read_dset(group_id, 'Egrid', Egrid, DIMENSION_EXT)
        CALL read_dset(group_id, 'ionisation_rates', ionisation_rates, DIMENSION_EXT)

        CALL h5gclose_f(group_id, error)
        CALL h5fclose_f(file_id, error)
        CALL h5close_f(error)

        open(UNIT=4,FILE='reference_table.dat',FORM="FORMATTED",action='write'); ! still kept, remove before final release
        DO i = 2, DIMENSION_EXT
          intensity = Egrid(i)
          intensity = intensity * intensity * field_intensity_au !rescale intensity intensity from atomic units to W/cm^2 (not SI)
          ionisation_rate = ionisation_rates(i)
          
          WRITE(4, '(3(2x, e12.5))') sqrt(intensity/field_intensity_au), intensity, ionisation_rate
          EXT_TABLE(i, 1) = intensity * intensity_factor   ! normalised intensity (C.U.)
          EXT_TABLE(i, 2) = ionisation_rate * rate_factor  ! ionisation rate (C.U.)
          IF (ionisation_rate.EQ.0.D0) THEN ! ionisation for small fields can be zero in stored precision
             EXT_TABLE(i, 3) = 0.D0
          ELSE
             EXT_TABLE(i, 3) = MPA_factor * ( ionisation_rate  * rate_factor/ intensity )  ! Normalised MPA
          ENDIF
       ENDDO

       EXT_TABLE(1, 1) = 0.d0; EXT_TABLE(1, 2) = 0.d0; EXT_TABLE(1, 3) = 0.d0 ! first line is imposed as 0-ionisation for 0-field
       
       close(4)

    ELSE ! As supposed as a testing case, there are some non-hdf5 outputs
      IF(my_rank.EQ.0) DIMENSION_EXT = filelength('rates_atomic.dat')
      CALL MPI_BCAST(DIMENSION_EXT,1,MPI_INTEGER,0,MPI_COMM_WORLD,ierr)
      ALLOCATE(EXT_TABLE(DIMENSION_EXT, 3),dumvect(DIMENSION_EXT))
      ALLOCATE(Egrid(DIMENSION_EXT),ionisation_rates(DIMENSION_EXT))
      ! Rescale the Table (done only by one worker)
      IF(my_rank.EQ.0) THEN
        open(UNIT=3,FILE='rates_atomic.dat',FORM="FORMATTED",action='read');    
        open(UNIT=4,FILE='reference_table.dat',FORM="FORMATTED",action='write');

        DO i = 2, DIMENSION_EXT
          read(unit=3,fmt=*) intensity, ionisation_rate !the organisation:  ||field amplitude imag. part of eigenerg. || all in atomic units
          Egrid(i) = intensity; ionisation_rates(i) = ionisation_rate
          intensity = intensity * intensity * field_intensity_au !rescale intensity intensity from atomic units to W/cm^2 (not SI)

          WRITE(4, '(3(2x, e12.5))') sqrt(intensity/field_intensity_au), intensity, ionisation_rate
         !intensity = (i-1) * intensity_step
         !ionisation_rate = ionisation_rate_PPT(intensity)
          EXT_TABLE(i, 1) = intensity * intensity_factor     ! normalised intensity (C.U.)
          EXT_TABLE(i, 2) = ionisation_rate * rate_factor    ! ionisation rate (C.U.)
          IF (ionisation_rate.EQ.0.D0) THEN
            EXT_TABLE(i, 3) = 0.D0
          ELSE
            EXT_TABLE(i, 3) = MPA_factor * ( ionisation_rate  * rate_factor/ intensity )  ! Normalised MPA
          ENDIF
        ENDDO
        EXT_TABLE(1, 1) = 0.d0   ! first row in input stores the result of EXT for zero field, not meaningful + division                
        EXT_TABLE(1, 2) = 0.d0
        EXT_TABLE(1, 3) = 0.d0
        Egrid(1) = 0; ionisation_rates(1) = 0;
 
        close(3); close(4);
      ENDIF   
   
      ! Distribute the table among all workers
      IF(my_rank.EQ.0) dumvect = EXT_TABLE(:, 1)
      CALL MPI_BCAST(dumvect,DIMENSION_EXT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      IF(my_rank.NE.0)  EXT_TABLE(:, 1) = dumvect

      IF(my_rank.EQ.0) dumvect = EXT_TABLE(:, 2)
      CALL MPI_BCAST(dumvect,DIMENSION_EXT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      IF(my_rank.NE.0)  EXT_TABLE(:, 2) = dumvect

      IF(my_rank.EQ.0) dumvect = EXT_TABLE(:, 3)
      CALL MPI_BCAST(dumvect,DIMENSION_EXT,MPI_DOUBLE_PRECISION,0,MPI_COMM_WORLD,ierr)
      IF(my_rank.NE.0)  EXT_TABLE(:, 3) = dumvect
      deallocate(dumvect)
    ENDIF

    ! The ionisation table is provided in the outputs
    IF (my_rank.EQ.0) THEN 
      CALL h5open_f(error)
      CALL h5fopen_f(outfilename, H5F_ACC_RDWR_F, file_id, error)
      CALL h5gcreate_f(file_id, outgroupname, group_id, error)
      CALL create_dset(group_id, 'Egrid', Egrid, DIMENSION_EXT)
      CALL create_dset(group_id, 'ionisation_rates', ionisation_rates, DIMENSION_EXT)
      CALL h5_add_units_1D(group_id, 'Egrid', '[a.u.]')
      CALL h5_add_units_1D(group_id, 'ionisation_rates', '[a.u.]')
      CALL h5gclose_f(group_id, error)
      CALL h5fclose_f(file_id, error)
    ENDIF

    EXT_TABLE(dimension_EXT, 1)=EXT_TABLE(dimension_EXT, 1)*1.d99 ! for diverging fields
!    intensity_step_inv = 1.d0 /(intensity_step * intensity_factor)

    JLOTI=1
    JLOTIP1=1
    DEALLOCATE(Egrid,ionisation_rates); 

  END SUBROUTINE RESCALE_TABLE_EXT

  !--------------------------------------------

  SUBROUTINE INTERPOLATE_EXT(var1,mpa,intensity)

    IMPLICIT NONE

    REAL(8) var1,mpa,intensity
    INTEGER(4) jmin,jmax

    REAL(8) intensity_step

    INTERFACE
       SUBROUTINE HUNT(xx, n, x, jlo)
         INTEGER                       :: n
         DOUBLE PRECISION, DIMENSION(n) :: XX
         DOUBLE PRECISION             :: x
         INTEGER                      :: jlo
       END SUBROUTINE HUNT

    END INTERFACE


    CALL hunt(ext_table,dimension_ext,intensity,jlotip1)
    IF (jlotip1.GE.DIMENSION_EXT-1) THEN 
!       print *, "Ionisation table exceeded"
!       print *, intensity
!       READ(5,*)
       jlotip1=DIMENSION_EXT-1
    ENDIF
    IF (jlotip1.LT.1) jlotip1=1
    intensity_step = ext_table(jlotip1+1,1) - ext_table(jlotip1,1)
    mpa=        ((intensity-ext_table(jlotip1,1))*ext_table(jlotip1+1,3)+(ext_table(jlotip1+1,1)-intensity)*ext_table(jlotip1,3))/intensity_step
!    mpa=        ((intensity-ppt_table(jlotip1,1))*ppt_table(jlotip1+1,3)+(ppt_table(jlotip1+1,1)-intensity)*ppt_table(jlotip1,3))*intensity_step_inv
    jmin=MIN(jloti,jlotip1)
    jmax=MAX(jloti,jlotip1)
    var1=SUM(ext_table(jmin:jmax+1,2),DIM=1)/REAL(jmax-jmin+2,8)
    jloti=jlotip1

    RETURN
  END SUBROUTINE INTERPOLATE_EXT


END MODULE External_ionisation_table





! The following models are not modified during the last development, will they be released?

!========================================== PPT N2
! MODULE PPT_N2

!   ! This module is used for PPT. It does the following procedure
!   ! 1- Read the input parameters
!   ! 2- Find the intensity range to build the table
!   ! 3- Create the table
!   ! $- A routine to serach the intensity

!   USE MEDIUM_PARAMETER_N2
!   USE LASER_PARAMETER
!   USE PPT_PARAMETER_N2

!   IMPLICIT NONE
!   DOUBLE PRECISION, SAVE   :: INTENSITY_STEP,intensity_step_inv
!   INTEGER, SAVE            :: DIMENSION_PPT
!   INTEGER, SAVE            :: JLOTI, JLOTIP1                          ! Integer for Intensity Index
!   DOUBLE PRECISION, DIMENSION(:,:), ALLOCATABLE :: PPT_TABLE 



! CONTAINS

!   !----------------------------------------
!   !1- Read the input parameters

!   SUBROUTINE INITIALISE_PPT(THEORY)

!     IMPLICIT NONE
!     CHARACTER(LEN=3), INTENT(IN)  :: THEORY
    
!     IF ( THEORY == "PPT") THEN
!        CALL ATTRIBUTE_DATA_PPT
!     ELSE
!        print *, 'The THEORY FLAG for N2 is wrong, choose PPT'
!        print *, 'The code is stopped'
!        STOP      
!     ENDIF   
!     JLOTI=1
!     JLOTIP1=1


!   END SUBROUTINE INITIALISE_PPT

!   !---------------------------------------
!   ! 2- Find the intensity range to build the table

!   SUBROUTINE FIND_INTENSITY_AREA(THEORY)

!     ! This function solves the equation n2*I = rho / 2*rho_c
!     ! The intensities are done in W/cm2

!     IMPLICIT NONE
!     CHARACTER(LEN=3), INTENT(IN)  :: THEORY
!     DOUBLE PRECISION, PARAMETER   :: PI = 3.14159265d0
!     DOUBLE PRECISION, PARAMETER   :: epsilon = 0.1d0
!     DOUBLE PRECISION              :: INTENSITY_MAX
!     DOUBLE PRECISION              :: factor
!     DOUBLE PRECISION              :: intensity
!     DOUBLE PRECISION              :: ionisation_rate
!     DOUBLE PRECISION              :: delta

!     INTERFACE
!      FUNCTION IONISATION_RATE_PPT_N2(intensity)
!          DOUBLE PRECISION, INTENT(IN) :: intensity
!          DOUBLE PRECISION             :: IONISATION_RATE_PPT_N2
!      END FUNCTION IONISATION_RATE_PPT_N2
!     END INTERFACE


!     IF ( (THEORY == "PPT").OR.(THEORY == "ADK") ) factor = ( PI * 1.0d25 /( 45.5635d0**2 * 2.41889)) * critical_power * photon_energy**2 *&
!          ( atomic_density / critical_density ) *  MIN(pulse_duration,1000.d0)


!     intensity = 1.d10
!     delta = 1.d0

!     DO WHILE( (delta > epsilon).AND.(intensity <= 1.d16) )
!        intensity = intensity + 10.d0**(Int(log10(intensity))) /50.0
!      IF ( THEORY == "PPT") THEN
!         ionisation_rate = IONISATION_RATE_PPT_N2(intensity)
!      ELSE
!         print *, 'The THEORY FLAG is wrong, choose ADK, PPT or IRC'
!         print *, 'The code is stopped'
!         STOP
!      ENDIF
!        delta = ABS( (intensity - ionisation_rate * factor) / (ionisation_rate * factor)  )
!     ENDDO

!     INTENSITY_STEP = intensity / 5000.d0
!     INTENSITY_MAX  = intensity * 100.d0
!     DIMENSION_PPT  = Int(intensity_max / intensity_step) +1

!   END SUBROUTINE FIND_INTENSITY_AREA

!   !----------------------------------------
!   ! 3- Create the table

!   SUBROUTINE  FILL_TABLE(THEORY)
!     USE mpi_stuff
!     IMPLICIT NONE
!     CHARACTER(len = 3),INTENT(IN):: THEORY
!     DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
!     DOUBLE PRECISION            :: intensity_factor
!     DOUBLE PRECISION            :: rate_factor
!     DOUBLE PRECISION            :: MPA_factor
!     DOUBLE PRECISION            :: intensity
!     DOUBLE PRECISION            :: ionisation_rate
!     INTEGER                     :: i

!     INTERFACE
!      FUNCTION IONISATION_RATE_PPT_N2(intensity)
!          DOUBLE PRECISION, INTENT(IN) :: intensity
!          DOUBLE PRECISION             :: IONISATION_RATE_PPT_N2
!      END FUNCTION IONISATION_RATE_PPT_N2
!     END INTERFACE


!     ALLOCATE(PPT_TABLE(DIMENSION_PPT, 3))

!     ! Normalised factors
!     intensity_factor = 4.d0 * PI * beam_waist**2 * 1.d-9 / critical_power 
!     IF ( (THEORY == "PPT") .OR. (THEORY == "ADK") ) THEN
!     rate_factor      =(4.d0 * PI**2 * 1.d16 / (45.5635**2 * 2.41889)) * &
!          (atomic_density / critical_density) *&
!          (photon_energy**2) * (beam_waist**2) * pulse_duration
!     MPA_factor       = n0_indice * critical_density * ionisation_potential * 45.5635 * 4.359d-10 / &
!          (photon_energy * pulse_duration * 2.d0 * PI)
!     ENDIF


!     ! Fill the Table
!     PPT_TABLE(1, 1) = 0.d0                   
!     PPT_TABLE(1, 2) = 0.d0
!     PPT_TABLE(1, 3) = 0.d0


!     IF (THEORY == "PPT") THEN
!     DO i = 2, dimension_PPT
!        intensity = (i-1) * intensity_step
!        ionisation_rate = ionisation_rate_PPT_N2(intensity)
!        PPT_TABLE(i, 1) = intensity * intensity_factor                    ! Normalised Intensity
!        PPT_TABLE(i, 2) = ionisation_rate * rate_factor                   ! Gamma
!        PPT_TABLE(i, 3) = MPA_factor * ( ionisation_rate  * rate_factor/ intensity )    ! Normalised MPA
!     ENDDO
!     ENDIF

!     IF(my_rank.EQ.0) THEN
!        OPEN (UNIT = 2, FILE = 'ionisation_table_N2.dat', STATUS = 'unknown')
!        DO i = 1, dimension_PPT
!           WRITE(2, '(3(2x, e12.5))') PPT_TABLE(i, 1),PPT_TABLE(i, 2),PPT_TABLE(i, 3)
!        ENDDO
!        CLOSE(2)
!     ENDIF



!     PPT_TABLE(dimension_PPT, 1)=PPT_TABLE(dimension_PPT, 1)*1.d99
!     intensity_step_inv = 1.d0 /(intensity_step * intensity_factor)

!   END SUBROUTINE FILL_TABLE

!   !----------------------------------------
!   ! 4- Find the intensitty and Interpolate

!   SUBROUTINE INTERPOLATE_PPT(var1,mpa,intensity)

!     IMPLICIT NONE

!     REAL(8) var1,mpa,intensity
!     INTEGER(4) jmin,jmax

!     INTERFACE
!        SUBROUTINE HUNT(xx, n, x, jlo)
!          INTEGER                        :: n
!          DOUBLE PRECISION, DIMENSION(n) :: XX
!          DOUBLE PRECISION               :: x
!          INTEGER                      :: jlo
!        END SUBROUTINE HUNT

!     END INTERFACE


!     CALL hunt(ppt_table,dimension_ppt,intensity,jlotip1)
!     IF (jlotip1.GE.DIMENSION_PPT-1) THEN 
! !       print *, "Ionisation table N2 exceeded"
! !       print *, intensity
! !       READ(5,*)
!        jlotip1=DIMENSION_PPT-1
!     ENDIF
!     IF (jlotip1.LT.1) jlotip1=1
!     mpa=        ((intensity-ppt_table(jlotip1,1))*ppt_table(jlotip1+1,3)+(ppt_table(jlotip1+1,1)-intensity)*ppt_table(jlotip1,3))*intensity_step_inv
!     jmin=MIN(jloti,jlotip1)
!     jmax=MAX(jloti,jlotip1)
!     var1=SUM(ppt_table(jmin:jmax+1,2),DIM=1)/REAL(jmax-jmin+2,8)
!     jloti=jlotip1

!     RETURN
!   END SUBROUTINE INTERPOLATE_PPT


! END MODULE PPT_N2


!===================================================================================
!===================================================================================

!                            The following code is needed for the PPT MODULE


!===================================================================================
!===================================================================================

DOUBLE PRECISION FUNCTION IONISATION_RATE_PPT(intensity)

  ! Function which evaluates cthe ionisation rate in atomic unit in the case of gas
  ! Modified 19/04/04

  USE MEDIUM_PARAMETER
  USE LASER_PARAMETER
  USE PPT_PARAMETER

  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: electric_field_au = 5.14224d11
  DOUBLE PRECISION, PARAMETER :: epsilon_0 = 8.85419d-12
  DOUBLE PRECISION, PARAMETER :: light_velocity = 3.0d8
  DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
  DOUBLE PRECISION, INTENT(IN):: intensity 
  DOUBLE PRECISION            :: electric_field
  DOUBLE PRECISION            :: gamma
  DOUBLE PRECISION            :: G_factor


  INTERFACE

     FUNCTION A_FACTOR(gamma)
       double precision, intent(in)   :: gamma
       double precision               :: A_FACTOR
     END FUNCTION A_FACTOR

  END INTERFACE


  ! Initialise the parameters
  ELECTRIC_FIELD = sqrt(2.d0 * intensity * 1.0d4 /(light_velocity*epsilon_0*n0_indice)) / &
       electric_field_au
  GAMMA = photon_energy * sqrt(2.d0 * ionisation_potential) / electric_field 

  G_factor = (3.d0 /(2.d0 * gamma) ) * ( (1.d0 + 1.d0/(2.d0 * gamma**2) ) * &
       log( gamma+sqrt(1.d0+gamma**2) ) - sqrt(1.d0 + gamma**2) / (2.d0 * gamma))

  IONISATION_RATE_PPT  = sqrt(6.d0 / PI) * C_factor * f_factor * ionisation_potential * &
       (2.d0*coulomb_field/(electric_field*sqrt(1.d0+gamma**2)))**(2.d0*n_number-1.5d0) * &
       A_factor(gamma) * exp(-2.d0 * coulomb_field*G_factor / (3.d0 * electric_field) )



END FUNCTION IONISATION_RATE_PPT

! DOUBLE PRECISION FUNCTION IONISATION_RATE_PPT_N2(intensity)

!   ! Function which evaluates cthe ionisation rate in atomic unit in the case of gas
!   ! Modified on the 19/04/04

!   USE MEDIUM_PARAMETER_N2
!   USE LASER_PARAMETER
!   USE PPT_PARAMETER_N2

!   IMPLICIT NONE
!   DOUBLE PRECISION, PARAMETER :: electric_field_au = 5.14224d11
!   DOUBLE PRECISION, PARAMETER :: epsilon_0 = 8.85419d-12
!   DOUBLE PRECISION, PARAMETER :: light_velocity = 3.0d8
!   DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
!   DOUBLE PRECISION, INTENT(IN):: intensity 
!   DOUBLE PRECISION            :: electric_field
!   DOUBLE PRECISION            :: gamma
!   DOUBLE PRECISION            :: G_factor


!   INTERFACE

!      FUNCTION A_FACTOR_N2(gamma)
!        double precision, intent(in)   :: gamma
!        double precision               :: A_FACTOR_N2
!      END FUNCTION A_FACTOR_N2

!   END INTERFACE


!   ! Initialise the parameters
!   ELECTRIC_FIELD = sqrt(2.d0 * intensity * 1.0d4 /(light_velocity*epsilon_0*n0_indice)) / &
!        electric_field_au
!   GAMMA = photon_energy * sqrt(2.d0 * ionisation_potential) / electric_field 

!   G_factor = (3.d0 /(2.d0 * gamma) ) * ( (1.d0 + 1.d0/(2.d0 * gamma**2) ) * &
!        log( gamma+sqrt(1.d0+gamma**2) ) - sqrt(1.d0 + gamma**2) / (2.d0 * gamma))

!   IONISATION_RATE_PPT_N2  = sqrt(6.d0 / PI) * C_factor * f_factor * ionisation_potential * &
!        (2.d0*coulomb_field/(electric_field*sqrt(1.d0+gamma**2)))**(2.d0*n_number-1.5d0) * &
!        A_factor_N2(gamma) * exp(-2.d0 * coulomb_field*G_factor / (3.d0 * electric_field) )



! END FUNCTION IONISATION_RATE_PPT_N2


!==========================================================================

! DOUBLE PRECISION FUNCTION IONISATION_RATE_ADK(intensity)

!   ! Function which evaluates cthe ionisation rate in atomic unit in the case of gas
!   ! Modified on the 19/04/04

!   USE MEDIUM_PARAMETER
!   USE LASER_PARAMETER
!   USE PPT_PARAMETER

!   IMPLICIT NONE
!   DOUBLE PRECISION, PARAMETER :: electric_field_au = 5.14224d11
!   DOUBLE PRECISION, PARAMETER :: epsilon_0 = 8.85419d-12
!   DOUBLE PRECISION, PARAMETER :: light_velocity = 3.0d8
!   DOUBLE PRECISION, PARAMETER :: PI = 3.14159265d0
!   DOUBLE PRECISION, INTENT(IN):: intensity 
!   DOUBLE PRECISION            :: electric_field
!   DOUBLE PRECISION            :: gamma
!   DOUBLE PRECISION            :: G_factor


!   INTERFACE

!      FUNCTION A_FACTOR(gamma)
!        double precision, intent(in)   :: gamma
!        double precision               :: A_FACTOR
!      END FUNCTION A_FACTOR

!   END INTERFACE

!   ! Initialise the parameters
!   ELECTRIC_FIELD = sqrt(2.d0 * intensity * 1.0d4 /(light_velocity*epsilon_0*n0_indice)) / &
!        electric_field_au
!   GAMMA = photon_energy * sqrt(2.d0 * ionisation_potential) / electric_field 

!   G_factor = (3.d0 /(2.d0 * gamma) ) * ( (1.d0 + 1.d0/(2.d0 * gamma**2) ) * &
!        log( gamma+sqrt(1.d0+gamma**2) ) - sqrt(1.d0 + gamma**2) / (2.d0 * gamma))

!   IONISATION_RATE_ADK  = sqrt(6.d0 / PI) * C_ADK_factor * ionisation_potential * &
!        (2.d0*coulomb_field/(electric_field*sqrt(1.d0+gamma**2)))**(2.d0*n_number-2.5d0) * &
!        A_factor(gamma) * exp(-2.d0 * coulomb_field*G_factor / (3.d0 * electric_field) )

! END FUNCTION IONISATION_RATE_ADK


!==========================================================================

DOUBLE PRECISION FUNCTION A_FACTOR(gamma)

  ! Fonction qui calcule le facteur A(gamma)

  USE MEDIUM_PARAMETER
  USE LASER_PARAMETER

  IMPLICIT NONE

  INTEGER, PARAMETER            :: Nmax = 5000
  DOUBLE PRECISION, PARAMETER   :: PI = 3.14159265d0
  DOUBLE PRECISION, PARAMETER   :: tiny = 1.d-13
  DOUBLE PRECISION, intent(in)  :: gamma 
  DOUBLE PRECISION              :: last_value, delta, alpha, beta, nu 
  INTEGER                       :: n

  INTERFACE

     FUNCTION  DAWSON(x)
       DOUBLE PRECISION, intent(in) :: x
       DOUBLE PRECISION             :: DAWSON
     END FUNCTION DAWSON

  END INTERFACE

  last_value = 0.d0
  delta      = 1.d0
  A_FACTOR   = 0.d0
  alpha      = 2.d0 * ( log(gamma + sqrt(1.d0 + gamma**2)) - gamma / sqrt(1.d0 + gamma**2) )
  beta       = 2.d0 * gamma / sqrt(1.d0 + gamma**2)
  nu         =  (ionisation_potential / photon_energy) * (1.d0 + 1.d0 / (2.d0 * gamma**2) )
  n          = Int(nu+1.d0)


  DO WHILE( (delta>tiny).AND.(n<Nmax) )
     A_FACTOR = A_FACTOR + DAWSON(sqrt(beta*(n-nu))) * exp(-1.d0*alpha*(n-nu))
     delta = A_FACTOR - last_value
     last_value = A_FACTOR
     n = n+1
  END DO

  IF ( n>Nmax) THEN 
     print *, " The iteration in A_factor was stopped because Nmax was reached"
     print *, delta
  ENDIF

  A_FACTOR = A_FACTOR * (4.D0 / SQRT(3.D0*PI)) * gamma**2 / (1.d0 + gamma**2)

END FUNCTION A_FACTOR

! DOUBLE PRECISION FUNCTION A_FACTOR_N2(gamma)

!   ! Fonction qui calcule le facteur A(gamma)

!   USE MEDIUM_PARAMETER_N2
!   USE LASER_PARAMETER

!   IMPLICIT NONE

!   INTEGER, PARAMETER            :: Nmax = 5000
!   DOUBLE PRECISION, PARAMETER   :: PI = 3.14159265d0
!   DOUBLE PRECISION, PARAMETER   :: tiny = 1.d-13
!   DOUBLE PRECISION, intent(in)  :: gamma 
!   DOUBLE PRECISION              :: last_value, delta, alpha, beta, nu 
!   INTEGER                       :: n

!   INTERFACE

!      FUNCTION  DAWSON(x)
!        DOUBLE PRECISION, intent(in) :: x
!        DOUBLE PRECISION             :: DAWSON
!      END FUNCTION DAWSON

!   END INTERFACE

!   last_value = 0.d0
!   delta      = 1.d0
!   A_FACTOR_N2   = 0.d0
!   alpha      = 2.d0 * ( log(gamma + sqrt(1.d0 + gamma**2)) - gamma / sqrt(1.d0 + gamma**2) )
!   beta       = 2.d0 * gamma / sqrt(1.d0 + gamma**2)
!   nu         =  (ionisation_potential / photon_energy) * (1.d0 + 1.d0 / (2.d0 * gamma**2) )
!   n          = Int(nu+1.d0)


!   DO WHILE( (delta>tiny).AND.(n<Nmax) )
!      A_FACTOR_N2 = A_FACTOR_N2 + DAWSON(sqrt(beta*(n-nu))) * exp(-1.d0*alpha*(n-nu))
!      delta = A_FACTOR_N2 - last_value
!      last_value = A_FACTOR_N2
!      n = n+1
!   END DO

!   IF ( n>Nmax) THEN 
!      print *, " The iteration in A_factor was stopped because Nmax was reached"
!      print *, delta
!   ENDIF

!   A_FACTOR_N2 = A_FACTOR_N2 * (4.D0 / SQRT(3.D0*PI)) * gamma**2 / (1.d0 + gamma**2)

! END FUNCTION A_FACTOR_N2

!==========================================================================

! DOUBLE PRECISION FUNCTION IONISATION_RATE_CRISTAL(intensity)

! ! This function is based on the Keldysh theory to evaluate 
! ! the ionization rate for the crystal.
! ! The ionization rate is expressed in atomic unit

! USE MEDIUM_PARAMETER
! USE LASER_PARAMETER

! IMPLICIT NONE
! DOUBLE PRECISION, INTENT(IN)  :: intensity

! DOUBLE PRECISION, PARAMETER   :: ELECTRIC_FIELD_AU = 5.14224d11
! DOUBLE PRECISION, PARAMETER   :: LENGTH_AU = 5.29177e-9
! DOUBLE PRECISION, PARAMETER   :: TIME_AU = 2.41889e-17
! DOUBLE PRECISION, PARAMETER   :: ENERGY_AU = 4.35981e-18
! DOUBLE PRECISION, PARAMETER   :: INTENSITY_AU = 3.509338d16
! DOUBLE PRECISION, PARAMETER   :: epsilon_0 = 8.854188d-12
! DOUBLE PRECISION, PARAMETER   :: light_velocity = 2.997925d8
! DOUBLE PRECISION              :: PI 

! DOUBLE PRECISION              :: electric_field, gamma, gamma1, delta_omega


! INTERFACE
!    FUNCTION DELK(x)
!      DOUBLE PRECISION, intent(in) :: x
!      DOUBLE PRECISION             :: delk
!    END FUNCTION DELK

!    FUNCTION DELE(x)
!      DOUBLE PRECISION, intent(in) :: x
!      DOUBLE PRECISION             :: dele 
!    END FUNCTION DELE

!    FUNCTION Q_FACTOR(x,y)
!      DOUBLE PRECISION, intent(in) :: x,y
!      DOUBLE PRECISION             :: Q_FACTOR 
!    END FUNCTION Q_FACTOR

!    FUNCTION rate_tunnel(x,y)
!      DOUBLE PRECISION, intent(in) :: x,y
!      DOUBLE PRECISION             :: rate_tunnel
!    END FUNCTION

!    SUBROUTINE SIGMA_MPI_AU(sigma)
!      DOUBLE PRECISION, intent(out) :: sigma
!    END SUBROUTINE SIGMA_MPI_AU
! END INTERFACE

! PI = 4.d0 * Atan(1.d0)



! ELECTRIC_FIELD = sqrt(2.d0 * intensity * 1.d4 / (light_velocity*epsilon_0*n0_indice)) / ELECTRIC_FIELD_AU

! GAMMA = photon_energy * sqrt(reduced_mass * ionisation_potential) / ELECTRIC_FIELD
! GAMMA1 = GAMMA**2 / (1.d0 + GAMMA**2)

! delta_omega = 2.d0 * ionisation_potential* DELE(1.d0-gamma1) / (PI*photon_energy*sqrt(gamma1))

! IONISATION_RATE_CRISTAL = (2.d0 * photon_energy / (9.d0 * PI)) * ( reduced_mass * photon_energy / sqrt(gamma1) )**1.5 * &
!                          Q_FACTOR(gamma1, delta_omega) * &
!                          exp(-1.d0 * PI * Int(delta_omega+1.d0) *(DELK(GAMMA1)-DELE(GAMMA1)) / DELE(1.d0-GAMMA1))


! END FUNCTION IONISATION_RATE_CRISTAL

!====================================================================================

! DOUBLE PRECISION FUNCTION Q_FACTOR(x,y)

! ! Function which evaluates the Q function

!   IMPLICIT NONE
!   INTEGER, PARAMETER          :: Nmax = 500
!   INTEGER                     :: n
!   DOUBLE PRECISION, intent(in):: x,y
!   DOUBLE PRECISION            :: PI
!   DOUBLE PRECISION, PARAMETER :: TINY = 1.d-30
!   DOUBLE PRECISION            :: ARG_DAWSON, ARG_EXP, ALPHA
!   DOUBLE PRECISION            :: last_value, delta

! INTERFACE
!    FUNCTION DELK(x)
!      DOUBLE PRECISION, intent(in) :: x
!      DOUBLE PRECISION             :: delk
!    END FUNCTION DELK

!    FUNCTION DELE(x)
!      DOUBLE PRECISION, intent(in) :: x
!      DOUBLE PRECISION             :: dele 
!    END FUNCTION DELE

!    FUNCTION DAWSON(x)
!      DOUBLE PRECISION, intent(in) :: x
!      DOUBLE PRECISION             :: DAWSON
!    END FUNCTION

! END INTERFACE

! PI = 4.d0 * Atan(1.d0)

! ARG_DAWSON = PI**2 / (4.d0 * DELK(1.d0-x) * DELE(1.d0-x))
! ARG_EXP = -1.d0 * PI * ( DELK(x)-DELE(x) ) / DELE(1.d0-x)
! ALPHA = 2.d0 * ( Int(y+1.d0)-y ) 

! delta = 1.d0
! n = 0
! Q_FACTOR = 0.d0
! LAST_VALUE = Q_FACTOR




! DO WHILE( (delta>TINY).AND.(n<Nmax) )
!    Q_FACTOR = Q_FACTOR + EXP(ARG_EXP * n) * DAWSON( sqrt(ARG_DAWSON * ( ALPHA+n ) ) )
!    delta = Q_FACTOR - LAST_VALUE
!    LAST_VALUE = Q_FACTOR
!    n = n+1
! ENDDO

! IF (n>Nmax) THEN
!    PRINT *, ' The iteration in Q_factor function was stopped before convergence'
!    PRINT *, ' The precision of result is', delta 
! ENDIF

! Q_FACTOR = Q_FACTOR / sqrt(2.d0 * DELK(1.d0-x) / PI)

! END FUNCTION Q_FACTOR

!=====================================================================================================


! Routine to search the intensity in the code from numnerical recipes
SUBROUTINE hunt(xx,n,x,jlo)
  INTEGER(4) jlo,n
  REAL(8) x,xx(n)
  INTEGER(4) inc,jhi,jm
  LOGICAL ascnd
  ascnd=xx(n).gt.xx(1)
  if(jlo.le.0.or.jlo.gt.n)then
     jlo=0
     jhi=n+1
     goto 3
  endif
  inc=1
  if(x.ge.xx(jlo).eqv.ascnd)then
1    jhi=jlo+inc
     if(jhi.gt.n)then
        jhi=n+1
     else if(x.ge.xx(jhi).eqv.ascnd)then
        jlo=jhi
        inc=inc+inc
        goto 1
     endif
  else
     jhi=jlo
2    jlo=jhi-inc
     if(jlo.lt.1)then
        jlo=0
     else if(x.lt.xx(jlo).eqv.ascnd)then
        jhi=jlo
        inc=inc+inc
        goto 2
     endif
  endif
3 if(jhi-jlo.eq.1)return
  jm=(jhi+jlo)/2
  if(x.gt.xx(jm).eqv.ascnd)then
     jlo=jm
  else
     jhi=jm
  endif
  goto 3
END SUBROUTINE hunt

      FUNCTION dawson(x)
      INTEGER NMAX
      DOUBLE PRECISION dawson,H,A1,A2,A3
      DOUBLE PRECISION, intent(in) :: x
      PARAMETER (NMAX=6,H=0.4d0,A1=2.d0/3.d0,A2=0.4d0,A3=2.d0/7.d0)
      INTEGER i,init,n0
      DOUBLE PRECISION d1,d2,e1,e2,sum,x2,xp,xx,c(NMAX)
      SAVE init,c
      DATA init/0/
      if(init.eq.0)then
        init=1
        do 11 i=1,NMAX
          c(i)=exp(-((2.d0* dble(i)-1.d0)*H)**2)
11      continue
      endif
      if(abs(x).lt.0.2d0)then
        x2=x**2
        dawson=x*(1.d0-A1*x2*(1.d0-A2*x2*(1.d0-A3*x2)))
      else
        xx=abs(x)
        n0=2*nint(0.5d0*xx/H)
        xp=xx- dble(n0)*H
        e1=exp(2.d0*xp*H)
        e2=e1**2
        d1= dble(n0+1)
        d2=d1-2.d0
        sum=0.d0
        do 12 i=1,NMAX
          sum=sum+c(i)*(e1/d1+1.d0/(d2*e1))
          d1=d1+2.d0
          d2=d2-2.d0
          e1=e2*e1
12      continue
        dawson=0.5641895835d0*sign(exp(-xp**2),x)*sum
      endif
      return
      END

     FUNCTION gamm(xx)
      DOUBLE PRECISION gamm,xx
      INTEGER j
      DOUBLE PRECISION ser,stp,tmp,x,y,cof(6)
      SAVE cof,stp
      DATA cof,stp/76.18009172947146d0,-86.50532032941677d0,24.01409824083091d0,-1.231739572450155d0,.1208650973866179d-2,-.5395239384953d-5,2.5066282746310005d0/
      x=xx
      y=x
      tmp=x+5.5d0
      tmp=(x+0.5d0)*log(tmp)-tmp
      ser=1.000000000190015d0
      do 11 j=1,6
        y=y+1.d0
        ser=ser+cof(j)/y
11    continue
      gamm=tmp+log(stp*ser/x)
      gamm = exp(gamm)
      return
      END


!=============================================================================================================
! ELLIPTIC INTEGRAL COMPLETE used for the IONISATION_RATE_CRISTAL function
!=============================================================================================================

!  DOUBLE PRECISION FUNCTION DELK(xx)
! ! Compute the complete elliptic integral of the first kind. x must be lower than 1
!   IMPLICIT NONE
!   DOUBLE PRECISION, intent(in) :: xx
  
!   INTERFACE
!      FUNCTION rf(x,y,z)
!        DOUBLE PRECISION, intent(in) :: x, y, z
!        DOUBLE PRECISION             :: rf
!      END FUNCTION rf
!   END INTERFACE

!   DELK = rf(0.d0, 1.d0 - xx, 1.d0)

!   END FUNCTION DELK

!  DOUBLE PRECISION FUNCTION DELE(xx)
! ! Compute the complete elliptic integral of the second kind. x must be lower than 1
!   IMPLICIT NONE
!   DOUBLE PRECISION, intent(in) :: xx
  
!   INTERFACE
!      FUNCTION rf(x,y,z)
!        DOUBLE PRECISION, intent(in) :: x, y, z
!        DOUBLE PRECISION             :: rf
!      END FUNCTION rf
!      FUNCTION rd(x,y,z)
!        DOUBLE PRECISION, intent(in) :: x, y, z
!        DOUBLE PRECISION             :: rd
!      END FUNCTION rd
     
!   END INTERFACE

!   DELE = rf(0.d0, 1.d0 - xx, 1.d0) - xx * rd(0.d0, 1.d0 - xx, 1.d0) / 3.d0

!   END FUNCTION DELE  

!  DOUBLE PRECISION FUNCTION rf(x,y,z)
! ! Compute Carlson's elliptic integral of the first kind, Rf(x,y,z). x,y and z must be nonnegative, and
! ! at most one can be zero. TINY must be at least 5 times the machine underflow limit, BIG at most
! ! one fifth the machine overflow limit.
!    IMPLICIT NONE
!    DOUBLE PRECISION, intent(in) :: x, y, z
!    DOUBLE PRECISION, PARAMETER  :: ERRTOL=0.0025,TINY=1.5e-38, BIG=3.E37
!    DOUBLE PRECISION, PARAMETER  :: THIRD=1./3.,C1=1./24.,C2=.1,C3=3./44.,C4=1./14. 
!    DOUBLE PRECISION             :: alamb,ave,delx,dely,delz,e2,e3,sqrtx,sqrty,sqrtz,xt,yt,zt

!    if(min(x,y,z).lt.0..or.min(x+y,x+z,y+z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rf'
!       xt=x
!       yt=y
!       zt=z
! 1     continue
!         sqrtx=sqrt(xt)
!         sqrty=sqrt(yt)
!         sqrtz=sqrt(zt)
!         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
!         xt=.25*(xt+alamb)
!         yt=.25*(yt+alamb)
!         zt=.25*(zt+alamb)
!         ave=THIRD*(xt+yt+zt)
!         delx=(ave-xt)/ave
!         dely=(ave-yt)/ave
!         delz=(ave-zt)/ave
!       if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
!       e2=delx*dely-delz**2
!       e3=delx*dely*delz
!       rf=(1.+(C1*e2-C2-C3*e3)*e2+C4*e3)/sqrt(ave)
!       return
!  END FUNCTION rf

!  DOUBLE PRECISION FUNCTION rd(x,y,z)
! ! Compute Carlson's elliptic integral of the second kind, Rd(x,y,z), x and y must be nonnegative, and 
! ! at most one can be zero. z must be positive. TINY must be at least twice the negative 2/3 power of the
! ! machine overflow limit. BIG must be at most 0.1 * ERRTOL times the negative 2/3 power of the machine
! ! underflow limit.
!    IMPLICIT NONE
!    DOUBLE PRECISION, intent(in)  :: x,y,z
!    DOUBLE PRECISION, PARAMETER   :: ERRTOL=.0015,TINY=1.e-25,BIG=4.5E21
!    DOUBLE PRECISION, PARAMETER   :: C1=3./14.,C2=1./6.,C3=9./22.,C4=3./26.,C5=.25*C3,C6=1.5*C4
!    DOUBLE PRECISION              :: alamb,ave,delx,dely,delz,ea,eb,ec,ed,ee,fac,sqrtx,sqrty,sqrtz,sum,xt,yt,zt
!    if(min(x,y).lt.0..or.min(x+y,z).lt.TINY.or.max(x,y,z).gt.BIG)pause 'invalid arguments in rd'
!       xt=x
!       yt=y
!       zt=z
!       sum=0.
!       fac=1.
! 1     continue
!         sqrtx=sqrt(xt)
!         sqrty=sqrt(yt)
!         sqrtz=sqrt(zt)
!         alamb=sqrtx*(sqrty+sqrtz)+sqrty*sqrtz
!         sum=sum+fac/(sqrtz*(zt+alamb))
!         fac=.25*fac
!         xt=.25*(xt+alamb)
!         yt=.25*(yt+alamb)
!         zt=.25*(zt+alamb)
!         ave=.2*(xt+yt+3.*zt)
!         delx=(ave-xt)/ave
!         dely=(ave-yt)/ave
!         delz=(ave-zt)/ave
!       if(max(abs(delx),abs(dely),abs(delz)).gt.ERRTOL)goto 1
!       ea=delx*dely
!       eb=delz*delz
!       ec=ea-eb
!       ed=ea-6.*eb
!       ee=ed+ec+ec
!       rd=3.*sum+fac*(1.+ed*(-C1+C5*ed-C6*delz*ee)+delz*(C2*ee+delz*(-C3*ec+delz*C4*ea)))/(ave*sqrt(ave))
!       return
!  END FUNCTION rd

!  SUBROUTINE INITIALIZE_PPT_N2(ionisation_potential_N2,residue_charge_N2,n0_indice_N2,critical_density_N2,atomic_density_N2,reduced_mass_N2,angular_momentum_N2)
!     USE PPT_N2
!     IMPLICIT NONE

!     DOUBLE PRECISION ionisation_potential_N2
!     DOUBLE PRECISION residue_charge_N2
!     DOUBLE PRECISION n0_indice_N2
!     DOUBLE PRECISION critical_density_N2
!     DOUBLE PRECISION atomic_density_N2
!     DOUBLE PRECISION reduced_mass_N2
!     INTEGER angular_momentum_N2

!     ionisation_potential = ionisation_potential_N2
!     residue_charge = residue_charge_N2
!     n0_indice = n0_indice_N2 
!     critical_density = critical_density_N2 
!     atomic_density = atomic_density_N2
!     reduced_mass = reduced_mass_N2
!     angular_momentum = angular_momentum_N2

!     CALL INITIALISE_PPT('PPT')
!     CALL FIND_INTENSITY_AREA('PPT')
!     CALL FILL_TABLE('PPT')

!  END SUBROUTINE INITIALIZE_PPT_N2

!  SUBROUTINE INTERPOLATE_PPT_N2(var1,mpa,intensity)
!     USE PPT_N2
!     IMPLICIT NONE

!     REAL(8) var1,mpa,intensity

!     CALL interpolate_ppt(var1,mpa,intensity)

!     RETURN
!  END SUBROUTINE INTERPOLATE_PPT_N2
