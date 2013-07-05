!====================================================================
!>
!! CRTM AMSU-B module for STMAS assimilation
!!
!! \history \n
!! Creation: Yuanfu Xie 2013
!
!====================================================================
MODULE crtm_kmatrix

  USE CRTM_MODULE
  USE prmtrs_stmas_cloud, only : satellite_obs

  IMPLICIT NONE

  CHARACTER(*), PARAMETER :: program_name = 'lvd_da'

  ! Parameters which could be ported to sat_da.nl later.
  INTEGER, PARAMETER :: num_lvd_sensor = 1
  INTEGER, PARAMETER :: num_lvd_absorber = 2
  INTEGER, PARAMETER :: num_lvd_cloud = 0
  INTEGER, PARAMETER :: num_lvd_aerosol = 0
  REAL(FP),    PARAMETER :: lvd_zenith_angle=-34.04_fp, lvd_scan_angle=-30.20_fp

  INTEGER :: num_lvd_channels

  type(CRTM_ChannelInfo_Type), allocatable:: lvd_ChannelInfo(:)
  type(CRTM_Geometry_Type),    allocatable:: lvd_GeometryInfo(:)
  type(CRTM_Atmosphere_Type),  allocatable:: lvd_Atm(:)     ! Profiles
  type(CRTM_Atmosphere_Type),  allocatable:: lvd_Atm_K(:,:) ! Channel x profiles
  type(CRTM_Surface_Type),     allocatable:: lvd_Sfc(:)     ! Profiles
  type(CRTM_Surface_Type),     allocatable:: lvd_Sfc_K(:,:) ! Channel x profile
  ! RTSolution and RTSolution_K: Channels x profiles
  type(CRTM_RTSolution_Type),  allocatable:: lvd_RTSolution(:,:), &
                                             lvd_RTSolution_K(:,:)
CONTAINS

!====================================================================
!> 
!! CRTM configuration for GOES.
!!
!! \history \n
!! Creation: Yuanfu Xie 2013, S Albers (modified for GOES)
!
!====================================================================

SUBROUTINE conf4lvd(lvd,numgrid,sensor_id,coe_path,moist_unit,&
                      badflag)

  IMPLICIT NONE

  CHARACTER(*), INTENT(IN) :: sensor_id,coe_path,moist_unit
  INTEGER, INTENT(IN) :: numgrid(4)
  REAL,    INTENT(IN) :: badflag
  TYPE(SATELLITE_OBS), INTENT(IN) :: lvd
  REAL(fp) :: angle(4,lvd%numpro)

  ! Local variables:
  CHARACTER(256) :: message
  INTEGER :: istatus,np

  ! Allocate memory for AMSU-B assimilation:
  ALLOCATE(lvd_ChannelInfo(num_lvd_sensor), &
           lvd_GeometryInfo(lvd%numpro), &
           STAT=istatus)

  istatus = CRTM_Init( (/sensor_id/), lvd_ChannelInfo, &
                        EmisCoeff_File='Wu-Smith.EmisCoeff.bin', &
                        File_Path=TRIM(coe_path)//'Coefficient_Data/', &
                        Quiet = .TRUE.)  
  IF ( istatus /= SUCCESS ) THEN
    message = 'Error initializing CRTM'
    CALL Display_Message( PROGRAM_NAME, message, FAILURE )
    STOP
  ELSE
    message = 'CRTM initialized'
    call Display_Message( program_name, message, SUCCESS )
  ENDIF
  num_lvd_channels = SUM(lvd_ChannelInfo%n_channels)
  IF (num_lvd_channels .NE. lvd%numchs) THEN
    PRINT*,'Numbers of AMSUB channels do not match!'
    STOP
  ENDIF

  ! Allocate memory for CRTM arrays:
  ALLOCATE(lvd_RTSolution(lvd%numchs,lvd%numpro), &
           lvd_RTsolution_K(lvd%numchs,lvd%numpro), &
           lvd_Atm(lvd%numpro), &
           lvd_Sfc(lvd%numpro), &
           lvd_Atm_K(lvd%numchs,lvd%numpro), &
           lvd_Sfc_K(lvd%numchs,lvd%numpro), &
           STAT=istatus)

  ! Create CRTM:
  CALL CRTM_Atmosphere_Create(lvd_Atm,numgrid(3),num_lvd_absorber, &
                              num_lvd_cloud,num_lvd_aerosol)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(lvd_Atm)) ) THEN
     message = 'Error allocating CRTM Atmosphere structures'
     CALL Display_Message( PROGRAM_NAME, message, FAILURE )
     STOP
  ENDIF
  CALL CRTM_Atmosphere_Create(lvd_Atm_K,numgrid(3),num_lvd_absorber, &
                              num_lvd_cloud,num_lvd_aerosol)
  IF ( ANY(.NOT. CRTM_Atmosphere_Associated(lvd_Atm)) ) THEN
     message = 'Error allocating CRTM Atmosphere_K structures'
     CALL Display_Message( PROGRAM_NAME, message, FAILURE )
     STOP
  ENDIF

  ! Data type conversion: real to real(fp) of CRTM:
  angle = lvd%angles

  ! Geometry setting:
  CALL CRTM_Geometry_SetValue( lvd_GeometryInfo, &
                               iFOV = lvd%numpro, &
                               Sensor_Zenith_Angle  = angle(1,1:lvd%numpro), & 
                               Sensor_Azimuth_Angle = angle(2,1:lvd%numpro), &
                               Source_Zenith_Angle  = angle(3,1:lvd%numpro), &
                               Source_Azimuth_Angle = angle(4,1:lvd%numpro) ) 

  ! Pass the non-control profile values to CRTM:
  DO np=1,lvd%numpro

    IF (lvd%lndcvr(np) .GT. 0) THEN
      lvd_Sfc(np)%Land_Coverage  = 1.0_fp
      lvd_Sfc(np)%Water_Coverage  = 0.0_fp
    ELSE
      lvd_Sfc(np)%Land_Coverage  = 0.0_fp
      lvd_Sfc(np)%Water_Coverage  = 1.0_fp
    ENDIF
    lvd_Sfc(np)%Land_type      = 1 ! COMPACTED SOIL, ACCORDING CRTM 2.1.1
    lvd_Sfc(np)%Water_type     = 1 ! SEA_WATER
    lvd_Atm(np)%Absorber_Id    = (/ H2O_ID , O3_ID /)
    lvd_Atm(np)%Absorber(:,2)  = 2.53e-3  ! Fake value as MHS does not use O3

    ! Currently, use mass mixing ratio only CRTM user guide 2.1.1 section 4.6:
    IF (moist_unit .EQ. 'MASS_MIXING_RATIO_UNITS') THEN
      lvd_Atm(np)%Absorber_Units = (/ MASS_MIXING_RATIO_UNITS,  &
                                        VOLUME_MIXING_RATIO_UNITS /)      
    ELSE
      PRINT*,'Not implemented, check the namelist, sat_da.nl'
      STOP
    ENDIF
    
  ENDDO

END SUBROUTINE conf4lvd

!====================================================================
!>
!! CRTM forward and K-matrix calculation for AMSU-B.
!!
!! \history \n
!! Creation: Yuanfu Xie 2013
!
!====================================================================

SUBROUTINE lvd_costgrad(numgrid,u_sfc,v_sfc,pres,sh,temp,lvd,tf, &
                          badflag,cost,grad,lvd_channels_used, &
                          iu,iv,it,iq)

  IMPLICIT NONE

  TYPE(satellite_obs), INTENT(IN) :: lvd

  INTEGER, INTENT(IN) :: numgrid(4),lvd_channels_used(*)
  INTEGER, INTENT(IN) :: tf,iu,iv,it,iq  ! Time frame and indices for control variables
  REAL,    INTENT(IN) :: pres(numgrid(3))
  REAL,    INTENT(IN) ::   sh(lvd%numpro,numgrid(3))
  REAL,    INTENT(IN) :: temp(lvd%numpro,numgrid(3))
  REAL,    INTENT(IN) :: u_sfc(lvd%numpro)
  REAL,    INTENT(IN) :: v_sfc(lvd%numpro)
  REAL,    INTENT(IN) :: badflag
  REAL,   INTENT(OUT) :: cost,grad(numgrid(1),numgrid(2),numgrid(3),numgrid(4),*)

  ! Local variables:
  REAL, PARAMETER :: pi=3.1415926
  CHARACTER(256) :: message
  DOUBLE PRECISION :: f
  INTEGER :: np,nl,nc,istatus
  REAL    :: bt,spd,weight

  PRINT*,'AMSUB cost function and gradient calculation...'

  DO np=1,lvd%numpro

    ! 2D fields:
    lvd_Sfc(np)%Land_Temperature = temp(np,1)
    lvd_Sfc(np)%Wind_Speed = SQRT(u_sfc(np)**2+v_sfc(np)**2)
    lvd_Sfc(np)%Wind_Direction = (2.0*pi+ASIN(v_sfc(np)/lvd_Sfc(np)%Wind_Speed))*180.0/pi

    ! 3D fields:
    ! Note: future improvement would be using terrain-following coordinate:
    ! Currently, we do not treat terrain!!!
    lvd_Atm(np)%Level_Pressure(0) = 1.5*pres(numgrid(3))-0.5*pres(numgrid(3)-1)
    DO nl=1,numgrid(3)-1
      lvd_Atm(np)%Level_Pressure(nl) = 0.5*(pres(numgrid(3)-nl)+pres(numgrid(3)-nl+1))
      lvd_Atm(np)%Pressure(nl)       = pres(numgrid(3)-nl+1)
    ENDDO
    lvd_Atm(np)%Level_Pressure(numgrid(3)) = 1.5*pres(1)-0.5*pres(2)
    lvd_Atm(np)%Pressure(numgrid(3)) = pres(1)

    ! Assign meteorological profiles to CRTM structured data:

    ! Convert specific humidity to mass mixing ratio:
    ! From the book of "Fundamentals of Atmospheric Modeling by Mark Z. Jacobson,
    ! mass mixing ratio (mmr) = rho_v/rho_d (equation 2.26) and 
    ! specific humidity (sh)  = rho_v/(rho_d+rho_v) (equation 2.27). Thus
    ! mmr = sh/(1-sh).
    DO nl=1,numgrid(3)
      ! CRTM mass mixing ratio is in g/kg: Table 4.7 CRTM user guide 2.1.1:
      lvd_Atm(np)%Absorber(nl,1) = 1.0e3*sh(np,numgrid(3)-nl+1)/(1.0e3-sh(np,numgrid(3)-nl+1))
      IF (lvd_Atm(np)%Absorber(nl,1) .lt. 0.0) then
        PRINT*,'Negative Vapor? ',lvd_Atm(np)%Absorber(nl,1),nl,np
        lvd_Atm(np)%Absorber(nl,1) = 0.0
      ENDIF

      lvd_Atm(np)%Temperature(nl) = temp(np,numgrid(3)-nl+1)
    ENDDO
    
  ENDDO

  ! Clean up arrays for K-matrix:
  CALL CRTM_Atmosphere_Zero( lvd_Atm_k )
  CALL CRTM_Surface_Zero( lvd_Sfc_K )
  lvd_RTSolution_K%brightness_Temperature = ONE
  lvd_RTSolution_K%Radiance = ZERO

  print*,'Calling k_matrix...'

  ! K-Matrix:
  istatus = CRTM_K_Matrix( lvd_Atm         , &
                           lvd_Sfc         , &
                           lvd_RTSolution_K, &
                           lvd_GeometryInfo, &
                           lvd_ChannelInfo , &
                           lvd_Atm_K       , &
                           lvd_Sfc_K       , &
                           lvd_RTSolution  )
  IF ( istatus /= SUCCESS ) THEN
    Message = 'Error in CRTM K_Matrix Model'
    CALL Display_Message( program_name, message, FAILURE )
    STOP
  ENDIF

  ! Cost and grad:
  f = 0.0D0
  ! Break function and gradient for calculating weight:
  weight = 0.0
  DO np=1,lvd%numpro

    DO nc=1,lvd%numchs
      IF (lvd_channels_used(nc) .NE. 1) CYCLE

      bt = lvd_RTSolution(nc,np)%Brightness_Temperature

      ! Cost:
      f = f+(bt-lvd%satval(nc,np))**2
      weight = weight+0.5
    ENDDO

  ENDDO
  ! Pass the double precision cost function value to the output variable:
  ! weight = SQRT(weight)

  cost = 1.0*f/weight

  DO np=1,lvd%numpro

    DO nc=1,lvd%numchs
      IF (lvd_channels_used(nc) .NE. 1) CYCLE

      bt = lvd_RTSolution(nc,np)%Brightness_Temperature

      ! Grad:
      ! Surface:
      spd = SQRT(u_sfc(np)**2+v_sfc(np)**2)
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,it) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,it)+(bt-lvd%satval(nc,np))/weight* &
        lvd_Sfc_K(nc,np)%Land_Temperature
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu)+(bt-lvd%satval(nc,np))/weight* &
        (lvd_Sfc_K(nc,np)%Wind_Speed*     u_sfc(np)/spd - &
         lvd_Sfc_K(nc,np)%Wind_Direction*SIGN(1.0,u_sfc(np))*v_sfc(np)/(spd*spd))*180.0/pi
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iu) = &
      grad(lvd%grdidx(1,np),lvd%grdidx(2,np),1,tf,iv)+(bt-lvd%satval(nc,np))/weight* &
        (lvd_Sfc_K(nc,np)%Wind_Speed*     v_sfc(np)/spd + &
         lvd_Sfc_K(nc,np)%Wind_Direction*ABS(u_sfc(np))/(spd*spd))*180.0/pi

      ! 3D profiles:
      DO nl=1,numgrid(3)
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,it) = &
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,it) + &
          (bt-lvd%satval(nc,np))/weight*lvd_Atm_K(nc,np)%Temperature(nl)
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,iq) = &
        grad(lvd%grdidx(1,np),lvd%grdidx(2,np),numgrid(3)-nl+1,tf,iq) + &
          (bt-lvd%satval(nc,np))/weight*lvd_Atm_K(nc,np)%Absorber(nl,1)* &
          1.0e6/(1.0e3-sh(np,numgrid(3)-nl+1))**2
      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE lvd_costgrad

!====================================================================
!>
!! Destroy the CRTM for AMSU-B.
!!
!! \history \n
!! Creation: Yuanfu Xie 2013
!
!====================================================================

SUBROUTINE dstrylvd

  ! Local variables:
  INTEGER :: istatus

  istatus = CRTM_Destroy(lvd_ChannelInfo)

END SUBROUTINE dstrylvd

END MODULE crtm_kmatrix

