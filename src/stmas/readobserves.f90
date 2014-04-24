MODULE READOBSERVES
!*************************************************
! READ IN RADAR OBSERVATIONS
! HISTORY: JANUARY 2008, CODED by ZHONGJIE HE.
!*************************************************
  USE PRMTRS_STMAS
  USE GENERALTOOLS, ONLY : VRTCLPSTN, VRTCLPSTN8, GETOBDATE, INTERPLTN, INTERPLTN_XIE, ZPCONVERT, OBSTOGRID
  USE READ_BACKGRD, ONLY : BK0, X00, Y00, P00, Z00, T00, HEIGHTU, HEIGHTL

  PRIVATE  HANDLEOBS, HT_TO_PRS, LURAO, LUNDX, LUNEW
  PUBLIC   RDRADROBS, RDBUFROBS, RDBUFROBS_XIE, ALLOCATOB, DEALOCTOB, RDOBSTEST
  PUBLIC   NOBSMAX, OBP, OBS, OBE, NST, OBA      !  , NALLOBS

  PUBLIC   X_RADAR, Y_RADAR                      ! JUST FOR TEST

  INTEGER  , PARAMETER :: LURAO=11,LUNDX=21,LUNEW=50,NOBSMAX=10000000
  REAL     , ALLOCATABLE :: OBP(:,:,:), OBS(:,:), OBE(:,:), OBA(:,:)
  INTEGER  , ALLOCATABLE :: NST(:)

!**************************************************
!COMMENT:
!   THIS MODULE IS MAINLY USED BY INPUT_BG_AND_OBS.F90.
!   SUBROUTINES:
!      ALLOCATOB: MEMORY ALLOCATE FOR OBP, OBS, OBE, OBA AND STT.
!      DEALOCTOB: MEMORY DEALLOCATE FOR OBP, OBS, OBE, AND STT.
!      RDRADROBS: READ IN DATA OF RADIAL WIND FROM LAPS.
!      RDBUFROBS: READ IN CONVENTIONAL DATA FROM LAPS.
!      HANDLEOBS: DETERMINE THE POSITION OF THE OBSERVATION IN THE BACKGROUND FIELD.
!      HT_TO_PRS: TRANSLATE THE HEIGHT OF THE OBSERVATION TO PRESURE ACCORDING TO THE BACKGROUND, USED FOR THE CASE OF PRESURE COORDINATE.
!      RDOBSTEST: JUST A TEST SUBROUTINE FOR READING SOME IDEAL DATAS FOR A TEST CASE.
!
!   ARRAYS:
!      OBP: POSITION OF EACH OBSERVATION IN THE BACKGROUND.
!      OBS: VALUE OF OBSERVATION.
!      OBE: ERROR OF OBSERVATION.
!      NST: NUMBER OF EACH VARIABLE OBSERVATION DATAS.
!      OBA: AZIMUTH AND TILT ANGLES OF OBSERVATION.
!
!   VARIABLE:
!       X_RADAR: THE LONGITUDE OF RADAR, JUST USED BY output_analysis.f90 TO DRAW PICTURES IN THE TEST CASE
!       Y_RADAR: THE LATITUDE OF RADAR, JUST USED BY output_analysis.f90 TO DRAW PICTURES IN THE TEST CASE
!**************************************************
CONTAINS

SUBROUTINE ALLOCATOB
!*************************************************
! MEMORY ALLOCATE FOR OBP, OBS, OBE, OBA AND STT
! HISTORY: SEPTEMBER 2007, CODED BY WEI LI
! HISTORY: JANUARY 2008, DEPARTED FORM THE MODULE 'INPUT_BG_OBS' by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: S,ER
! --------------------

  ALLOCATE(OBP(NUMDIMS,NOBSMAX,NUMSTAT+3),STAT=ER)
  IF(ER.NE.0)STOP 'OBP ALLOCATE WRONG'
  ALLOCATE(OBS(NOBSMAX,NUMSTAT+3),STAT=ER)
  IF(ER.NE.0)STOP 'OBS ALLOCATE WRONG'
  ALLOCATE(OBE(NOBSMAX,NUMSTAT+3),STAT=ER)
  IF(ER.NE.0)STOP 'OBE ALLOCATE WRONG'
  ALLOCATE(NST(NUMSTAT+3),STAT=ER)
  IF(ER.NE.0)STOP 'NST ALLOCATE WRONG'
  DO S=1,NUMSTAT+3
    NST(S)=0
  ENDDO
  ALLOCATE(OBA(NOBSMAX,3),STAT=ER)
  IF(ER.NE.0)STOP 'OBA ALLOCATE WRONG'

  ! ALLOCATE GRID SPACE FOR LAPS RADIAL WIND:
  ! ALLOCATE(IDXRADWN(4,FCSTGRD(1)*FCSTGRD(2)*FCSTGRD(3)*FCSTGRD(4)), &
  !            LAPSRADW(FCSTGRD(1)*FCSTGRD(2)*FCSTGRD(3)*FCSTGRD(4)), STAT=ER)
  ! IF (ER .NE. 0) THEN
  !   print*,'Size too large: ',FCSTGRD(1)*FCSTGRD(2)*FCSTGRD(3)*FCSTGRD(4)
  !   print*,'Size too large: ',FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4)
  !   PRINT*,'ALLOCATOB: cannot allocate memory for LAPS radial wind'
  !   STOP
  ! ENDIF

  RETURN
END SUBROUTINE ALLOCATOB


SUBROUTINE DEALOCTOB
!*************************************************
! MEMORY DEALLOCATE FOR OBP, OBS, OBE, AND STT
! HISTORY: SEPTEMBER 2007, CODED BY WEI LI
! HISTORY: JANUARY 2008, DEPARTED FORM THE MODULE 'INPUT_BG_OBS' by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------

  DEALLOCATE(OBP)
  DEALLOCATE(OBS)
  DEALLOCATE(OBE)
  DEALLOCATE(NST)
  DEALLOCATE(OBA)

  ! DEALLOCATE LAPS RADIAL WIND:
  ! DEALLOCATE(IDXRADWN, LAPSRADW)

  RETURN
END SUBROUTINE DEALOCTOB

SUBROUTINE RDLAPSRDR
!*************************************************
!  THIS ROUTINE READS GRIDDED RADAR DATA FROM LAPS
!  REMAPPED DATASET USING LAPS:
!               GET_MULTIRADAR_VEL AND
!               QC_RADAR_OBS
!       ASSUMING THE BACKGROUND HAS BEEN READ IN.
!
!  HISTORY: MARCH 2008, CODED BY YUANFU XIE.
!*************************************************

  USE MEM_NAMELIST              ! LAPS WIND PARAMETER MODULE

  IMPLICIT NONE

  ! LOCAL VARIABLES:
  CHARACTER*31 :: RADEXT(MAX_RADARS)    ! POSSIBLE RADAR NAME EXTENSIONS
  CHARACTER*4  :: RADNAM(MAX_RADARS)    ! RADAR STATION NAMES
  CHARACTER*8  :: SID			! RADAR STATION NAME IN 8 BYTE
  INTEGER      :: NRADAR                ! NUMBER OF RADAR AVAILABLE
  INTEGER      :: STTRAD,STTNQY         ! RADAR AND ITS NYQUIST STATUS
  INTEGER      :: NGRDRD(MAX_RADARS)    ! NUMBER OF GRIDPOINTS WITH MEASURABLE VEL
  INTEGER      :: RADTIM(MAX_RADARS)    ! RADAR OBSERVATION TIME
  INTEGER      :: RADIDS(MAX_RADARS)    ! RADAR IDX
  INTEGER      :: IOFFSET(MAX_RADARS),JOFFSET(MAX_RADARS) ! OFFSET FOR THE NEW get_multiradar_vel
  INTEGER      :: I,J,K,L,M,IX0,IY0,IX1,IY1	! GRID INDICES, NUMBER OF RADAR, TIME FRAME
  LOGICAL      :: CLUTTR                ! .TRUE. -- REMOVE 3D RADAR CLUTTER
  REAL         :: RADVEL(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),MAX(3,MAX_RADARS))  ! 3 To allow cloud liquid read in
                  ! RADAR 4D VELOCITY GRID
  REAL         :: RADNQY(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),MAX(3,MAX_RADARS))  ! 3 To allow cloud ice read in
                  ! RADAR 4D NYQUIST VELOCITY
  REAL         :: UVZERO(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),3)    ! INCREASED TO 3 SHARING WITH CLOUD READ
                  ! ZERO UV GRIDS USED FOR CALLING LAPS QC_RADAR_OBS
  REAL         :: VOLNQY(MAX_RADARS)            ! VOLUME NYQUIST VELOCITY
  REAL         :: RADLAT(MAX_RADARS),RADLON(MAX_RADARS),RADHGT(MAX_RADARS)
  REAL         :: UVGRID(2)

  INTEGER      :: TOLTIM,INC        ! RADAR DATA TOLERATE WINDOW, TIME INCREMENT
  REAL         :: XRADAR,YRADAR,ZRADAR  ! RADAR STATION LOCATION IN ANALYSIS GRID
  REAL         :: HEIGHT_TO_ZCOORD2     ! LAPS ROUTINE FOR CONVERTING HEIGHT TO GRID
  REAL         :: XSPACE,YSPACE ! SPACINGS

  ! VARIABLES FOR USING HANDLE OBS ROUTINE:
  INTEGER      :: IP
  REAL         :: OP(4),PRSLVL(FCSTGRD(3)),AZ,EA,OB,OE
  REAL         :: HEIGHT_OF_LEVEL,OH,SR		! OH OBS HEIGHT, SR RADAR SLANT RANGE

  !====== USED BY THE MODIFICATION OF ZHONGJIE HE
  INTEGER      :: I0,J0        ! INDEX OF THE GRID POINT AT THE WEST-SOUTH CORNER TO RADAR STATION
  REAL         :: XRDR,YRDR    ! LONGITUDE AND LATITUDE OF RADAR STATION USED TO CALCULATE AZ AND EA
  REAL         :: RE           ! RADIUS OF EARTH
  !====== END MODIFICATION OF ZHONGJIE HE
!jhui
  INTEGER :: TT,TT1,INC1,nn
  REAL    :: REFSCL
  REAL    :: lat(FCSTGRD(1),FCSTGRD(2))
  REAL    :: lon(FCSTGRD(1),FCSTGRD(2))
  REAL    :: topo(FCSTGRD(1),FCSTGRD(2))
  REAL    :: rlaps_land_frac(FCSTGRD(1),FCSTGRD(2))
  REAL    :: grid_spacing_cen_m
  INTEGER :: istatus, i4_tol,i4_ret,iqc_2dref
  CHARACTER :: units*10,comment*125,radar_name*4,iext*31,c_filespec*255
  REAL :: heights_3d(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL :: radar_ref_3d(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4))
  REAL :: closest_radar(FCSTGRD(1),FCSTGRD(2))
  REAL ::  rlat_radar(5), rlon_radar(5),rheight_radar(5)
  INTEGER :: n_ref_grids,n_2dref,n_3dref
!  INTEGER :: istat_radar_2dref_a,istat_radar_3dref_a
!!Variable's defination are not same with fountion read_multiradar_3dref ,modified by shuyuan 20100525
  INTEGER :: istat_radar_2dref_a(FCSTGRD(1),FCSTGRD(2))
  INTEGER :: istat_radar_3dref_a(FCSTGRD(1),FCSTGRD(2))
!! liu 20100525
  character*40 c_vars_req
  character*180 c_values_req
  INTEGER :: i4time_radar
  REAL     :: tempref,make_ssh,make_td,tw,make_rh

! add the following definition since build on 5/23/2011 failed. HJ 5/23/2011
  integer nx_r, ny_r 

  INCLUDE 'main_sub.inc'
  INCLUDE 'laps_cloud.inc'
  REAL :: CLD_PRS(KCLOUD)
  REAL :: CLOUD3D(FCSTGRD(1),FCSTGRD(2),KCLOUD)

  IF(IFPCDNT.NE.1) THEN              !   BY ZHONGJIEHE
    PRINT*, 'ERROR! LAPS IS PRESSURE COORDINATE! PLEASE SET IFPCDNT TO 1!'
    STOP
  ENDIF                              !   END OF MODIFIED BY ZHONGJIE HE

  ! CLUTTR = .TRUE.               ! TRUE. -- REMOVE 3D RADAR CLUTTER for old get_multiradar_vel
  TOLTIM = 600		! DEFAULT 10 MINUTE TOLERATE TIME WINDOW
  RE=6365.0E3                   ! ADDED BY ZHONGJIE HE
  IF (FCSTGRD(4) .GT. 1) TOLTIM = T00(FCSTGRD(4))-T00(1)  
  ! GET PRESSURE LEVELS:
  CALL GET_PRES_1D(LAPSI4T,FCSTGRD(3),PRSLVL,STTRAD)

  ! READ RADAR DATA CLOSE TO EACH TIME FRAME:
  INC = TOLTIM/(FCSTGRD(4)-1)

  DO M=1,FCSTGRD(4)
    ! GET UNFOLDED RADAR THROUGH LAPS:
    CALL GET_MULTIRADAR_VEL(ITIME2(1)+(M-1)*INC,INC/2,RADTIM,MAX_RADARS, &
                            NRADAR,RADEXT,RMISSING, &
                            FCSTGRD(1),FCSTGRD(2),FCSTGRD(3), &
                            LATITUDE,LONGITUD, &
                            FCSTGRD(1),FCSTGRD(2),0, &
                            RADVEL,RADNQY,RADIDS,VOLNQY, &
                            IOFFSET,JOFFSET,CLUTTR,NGRDRD, &
                            RADLAT,RADLON,RADHGT,RADNAM,STTRAD,STTNQY)

    ! SET UV ZERO GRIDS:
    UVZERO = 0.0

    ! NYQUIST UNFOLDING:
    DO L=1,NRADAR
      ! QC and unfolding radar nyquist:
      nx_r = FCSTGRD(1)
      ny_r = FCSTGRD(2)
      ioffset = 0
      joffset = 0
      CALL QC_RADAR_OBS(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3), &
                         RMISSING,nx_r,ny_r,ioffset,joffset,RADVEL(1,1,1,L),RADNQY(1,1,1,L),NGRDRD(L), &
                         LATITUDE,LONGITUD,RADLAT(L),RADLON(L),RADHGT(L), &
                         UVZERO(1,1,1,1),UVZERO(1,1,1,2),BK0(1,1,1,1,1),BK0(1,1,1,1,2), &
                         VOLNQY(L),L_CORRECT_UNFOLDING,L_GRID_NORTH,STTRAD)
      PRINT*,'STATUS OF QC_RADAR_OBS: ',STTRAD,' FOR TIME FRAME: ',M

      ! ASSIGN AZIMUTH AND ELEVATION ANGLES: BK0 STORES UVZT
      CALL LATLON_TO_RLAPSGRID(RADLAT(L),RADLON(L),LATITUDE,LONGITUD,&
                                FCSTGRD(1),FCSTGRD(2),XRADAR,YRADAR,STTRAD)

      ! Yuanfu found ZRADAR is not being used: turn off for now 07/2009
      ! For radar stations outside domain, this caused problem as 
      ! HEIGHT_TO_ZCOORD2 uses height(xradar,yradar...):
      !ZRADAR = HEIGHT_TO_ZCOORD2(RADHGT(L),BK0(1,1,1,1,3), &
      !                           FCSTGRD(1),FCSTGRD(2),FCSTGRD(3), &
      !                           NINT(XRADAR),NINT(YRADAR),STTRAD)

      ! LAPS GRID SPACING:
      CALL GET_GRID_SPACING_ACTUAL_XY(RADLAT(L),RADLON(L),XSPACE,YSPACE,STTRAD)

      DO K=1,FCSTGRD(3)
        ! NOTE LAPS REMAP ROUTINE CURRENTLY USES STANDARD ATMOSPHERE
        ! FOR OBS HEIGHT. WHEN IT CHANGES, THE HEIGHT FED TO 
        ! LATLON_TO_RADAR HAS TO CHANGED:
        OH = HEIGHT_OF_LEVEL(k)
        DO J=1,FCSTGRD(2)
          DO I=1,FCSTGRD(1)

            IF (RADVEL(I,J,K,L) .NE. RMISSING) THEN

!               PRINT*,'RDLAPSRDR: --RADIAL WIND: ', &
!		 RADVEL(I,J,K,L),RADNQY(I,J,K,L),I,J,K,L,XRADAR,YRADAR,NGRDRD(L),VOLNQY(L)

	       ! COMPUTE AZIMUTH AND ELEVATION ANGLES USING LAPS ROUTINE
	       ! LATLON_TO_RADAR.
               CALL LATLON_TO_RADAR(LATITUDE(I,J),LONGITUD(I,J),OH,AZ,SR,EA, &
                                     RADLAT(L),RADLON(L),RADHGT(L))

               OP(2) = LATITUDE(I,J)      ! OP(1): LONGITUDE; OP(2): LATITUDE
               OP(1) = LONGITUD(I,J)
               OP(3) = PRSLVL(K)          ! IN PASCAL

               ! FOR GIVEN TOLERATED TIME WINDOW, ASSUME AVAILABLE RADAR IS
               ! AT THE TIME FRAME:
               OP(4) = RADTIM(L)-ITIME2(1)	! ACTUAL RADAR TIME
               ! OP(4) = T00(M)

               IP = 1
               OB=RADVEL(I,J,K,L)		! POSITIVE WIND IS TOWARD THE STATION BY YUANFU
               OE=0.5
!jhui
!               OE=0.3
               SID(1:4) = RADNAM(L)

               CALL HANDLEOBS_SIGMA(OP,OB,OE,NUMSTAT+1,NALLOBS,IP,AZ,EA,SID)

           ENDIF
          ENDDO	! I
        ENDDO		! J
      ENDDO		! K
    ENDDO		! L -- RADAR LEVELS

  ENDDO		! M -- TIME FRAMES
  PRINT*,'NUMBER OF RADAR RADIAL WIND TAKEN AS OBS: ',NST(NUMSTAT+1)


  call get_laps_domain_95(FCSTGRD(1),FCSTGRD(2),lat,lon,topo &
                ,rlaps_land_frac,grid_spacing_cen_m,istatus)
  if(istatus .ne. 1)then
       write(6,*)' Error getting LAPS domain'
       return
  endif
  write(6,*)' Actual grid spacing in domain center = ',grid_spacing_cen_m

!=======reflectivity==for time cycle ,read multitime data file *.vrz======
!=========added by shuyuan liu 20100830================== 
  
  i4_tol = 900
  i4_ret = 0
  ref_base = -10
  nn =0
  REFSCL =0.0  
  iext="vrz"
  BK0(:,:,:,:,NUMSTAT+1) = 0.0
  radar_ref_3d = 0.0
  ! Get radar at current and previous time:
  DO L=1,2  !for L   time          
     call get_filespec(iext(1:3),2,c_filespec,istatus)
     call get_file_time(c_filespec,LAPSI4T+(L-2)*ICYCLE,i4time_radar)

     call read_multiradar_3dref(LAPSI4T+(L-2)*ICYCLE,i4_tol,i4_ret,&!I
                   .true.,ref_base,&                              ! I
                   FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),iext, &   ! I 
                   lat,lon,topo,.false.,.false., heights_3d, &  ! I
                   radar_ref_3d(1,1,1,L), &                                      ! O
                   rlat_radar,rlon_radar,rheight_radar,radar_name, &  ! O
                   iqc_2dref,closest_radar, &                          ! O
                   n_ref_grids,n_2dref,n_3dref,istat_radar_2dref_a, &  ! O
                   istat_radar_3dref_a)                              ! O
     DO K = 1, FCSTGRD(3)
      DO J = 1, FCSTGRD(2)
       DO I = 1, FCSTGRD(1) 
       IF(radar_ref_3d(I,J,K,L) .GT. 5. .AND. radar_ref_3d(I,J,K,L) .LT.100) THEN
            ! modified shuyuan 20100719
                        
            tempref=(radar_ref_3d(I,J,K,L)-43.1)/17.5
            tempref=(10.**tempref)       !g/m3
             ! shuyuan 20100719
            REFSCL = REFSCL + tempref**2             
                 
            ! CHECK IF RAIN AND SNOW IS ANALYZED:
            IF (NUMSTAT .LE. 5) GOTO 555 
            nn = nn + 1
            OP(2) = LATITUDE(I,J)      ! OP(1): LONGITUDE; OP(2): LATITUDE
            OP(1) = LONGITUD(I,J)
            OP(3) = PRSLVL(K)          ! IN PASCAL
!           OP(4) = RADTIM(L)-ITIME2(1)    ! ACTUAL RADAR TIME
            OP(4) = T00(L)
            OB= tempref         
            OE=0.01  ! shuyuan   test 0.1 0.01 1 
            SID(1:3) = "vrz"
            CALL HANDLEOBS_SIGMA(OP,OB,OE,NUMSTAT+3,NALLOBS,IP,AZ,EA,SID) 
            ! SKIP REFLECTIVITY:
555         CONTINUE 

            ! Modified by Yuanfu Xie Nov. 2011 for adding radar reflectivity generated SH obs:
            !NST(HUMIDITY) = NST(HUMIDITY)+1
            !OBP(1,NST(HUMIDITY),HUMIDITY) = I-1
            !OBP(2,NST(HUMIDITY),HUMIDITY) = J-1
            !OBP(3,NST(HUMIDITY),HUMIDITY) = K-1
            !OBP(4,NST(HUMIDITY),HUMIDITY) = L-1
            !OBS(NST(HUMIDITY),HUMIDITY) = MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,0.75,-132.0)
            !OBE(NST(HUMIDITY),HUMIDITY) = 1.0
            !NALLOBS = NALLOBS+1
        ENDIF
       ENDDO
      ENDDO
     ENDDO
  ENDDO  ! for L
  ! INTERPOLATION TO THE THREE TIME FRAMES OF STMAS ANALYSIS:
  radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),3) = &
    1.5*radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)- &
    0.5*radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1) = &
    0.5*radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)+ &
    0.5*radar_ref_3d(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)

  ! Calculate the SH bounds from radar reflectivity:
  DO L=1,FCSTGRD(4)
    DO K=1,FCSTGRD(3)
      DO J=1,FCSTGRD(2)
        DO I=1,FCSTGRD(1)
          IF(radar_ref_3d(I,J,K,L) .GT. 5. .AND. radar_ref_3d(I,J,K,L) .LT.100) THEN

          ! ASSUME 75% satured RH where reflectivity present as SH lower bounds:
          IF (BK0(I,J,K,L,TEMPRTUR) .GT. 273.15) THEN
            IF (radar_ref_3d(I,J,K,L) .GT. 45.0) THEN
              BK0(I,J,K,L,NUMSTAT+1)= &
                MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,1.0,0.0)
            ELSE
              BK0(I,J,K,L,NUMSTAT+1)= &
                MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,0.1+0.9*radar_ref_3d(I,J,K,L)/45.0,0.0)
            ENDIF
          ELSE
            IF (radar_ref_3d(I,J,K,L) .GT. 30.0) THEN
              BK0(I,J,K,L,NUMSTAT+1)= &
                MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,1.0,0.0)
            ELSE
              BK0(I,J,K,L,NUMSTAT+1)= &
                MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,0.2+0.8*radar_ref_3d(I,J,K,L)/30.0,0.0)
            ENDIF
          ENDIF

          ENDIF
        ENDDO 
      ENDDO 
    ENDDO 
  ENDDO 
  PRINT*,'Range of reflectivity derived bound: ', &
   maxval(BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2,NUMSTAT+1)), &
   minval(BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2,NUMSTAT+1))

  ! Interpolate RADREF to GRDBKGD0 as SH lower bound:
  IF (MAXGRID(3) .NE. FCSTGRD(3) .OR. MAXGRID(4) .NE. FCSTGRD(4)) THEN
    PRINT*,'The analysis grid does not match multigrid in Z or T, please check!'
    PRINT*,'The SH lower bound calculation assumes they are the same'
    STOP
  ENDIF

  ! READ IN LAPS CLOUD LIQUID AND ICE:
  iext = "lwc"
  DO L=1,2
    CALL GET_LAPS_3DGRID(LAPSI4T+(L-2)*ICYCLE,i4_tol,i4_ret, &
               FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),iext,iext, &
               units,comment,RADVEL(1,1,1,L),STTRAD)
  ENDDO
  ! INTERPOLATION TO THE THREE TIME FRAMES OF STMAS ANALYSIS:
  RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),3) = &
    1.5*RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)- &
    0.5*RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1) = &
    0.5*RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)+ &
    0.5*RADVEL(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  PRINT*,'LAPS cloud liquid: ',minval(RADVEL(:,:,:,2))*1000.0,maxval(RADVEL(:,:,:,2))*1000.0
  DO L=1,2
    CALL GET_LAPS_3DGRID(LAPSI4T+(L-2)*ICYCLE,i4_tol,i4_ret, &
               FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),iext,"ice", &
               units,comment,RADNQY(1,1,1,L),STTRAD)
  ENDDO
  ! INTERPOLATION TO THE THREE TIME FRAMES OF STMAS ANALYSIS:
  RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),3) = &
    1.5*RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)- &
    0.5*RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1) = &
    0.5*RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)+ &
    0.5*RADNQY(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  PRINT*,'LAPS cloud ice: ',minval(RADNQY(:,:,:,2))*1000.0,maxval(RADNQY(:,:,:,2))*1000.0

  ! READ IN LAPS CLOUD ANALYSIS FOR BOUNDS OF SH: USE OF ARRAY UZERO
  UVZERO = 0.0
  iext = "lc3"
  DO L=1,2
    CALL GET_CLOUDS_3DGRID(LAPSI4T+(L-2)*ICYCLE,J,FCSTGRD(1),FCSTGRD(2), &
           KCLOUD,iext,CLOUD3D,CLD_HTS,CLD_PRS,I)
    IF (J .NE. LAPSI4T+(L-2)*ICYCLE) THEN
      PRINT*,'lc3 time does not match the analysis, skip'
      cycle
    ELSE
      PRINT*,'readobservs: found file -- ',J
    ENDIF

    DO K=1,FCSTGRD(3)
      DO M=1,KCLOUD-1
        IF (PRSLVL(K) .LE. CLD_PRS(M) .AND. (PRSLVL(K) .GE. CLD_PRS(M+1))) THEN
          XSPACE = (LOG(CLD_PRS(M))-LOG(PRSLVL(K)))/(LOG(CLD_PRS(M))-LOG(CLD_PRS(M+1)))
          YSPACE = 1.0-XSPACE
          DO J=1,FCSTGRD(2)
          DO I=1,FCSTGRD(1)
            UVZERO(I,J,K,L) = YSPACE*CLOUD3D(I,J,M)+XSPACE*CLOUD3D(I,J,M+1)
          ENDDO
          ENDDO
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  ! INTERPOLATION TO THE THREE TIME FRAMES OF STMAS ANALYSIS:
  UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),3) = &
    1.5*UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)- &
    0.5*UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)
  UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1) = &
    0.5*UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2)+ &
    0.5*UVZERO(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),1)

  ! CONVERT CLOUD TO SH BOUNDS:
  DO L=1,3
    DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      ! XSPACE = 0.0
      ! IF (UVZERO(I,J,K,L) .GT. 0.1) THEN
      !   XSPACE = MAKE_SSH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15, &
      !                     0.8*(UVZERO(I,J,K,L)**0.2),0.0)
      ! ENDIF
      ! IF (XSPACE .GT. BK0(I,J,K,L,NUMSTAT+1)) BK0(I,J,K,L,NUMSTAT+1)=XSPACE

      ! Adjust SH bounds based on cloud fraction:
      IF ( (UVZERO(I,J,K,L) .LT. 0.1) .AND. (radar_ref_3d(I,J,K,L) .GT. 5.0) ) THEN
        ! XSPACE = MAKE_TD(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15, & ! DEW
        !                  BK0(I,J,K,L,NUMSTAT+1),0.0)
        ! YSPACE = TW(BK0(I,J,K,L,TEMPRTUR)-273.15,XSPACE,PRSLVL(K)/100.0) ! WET BULB T
        ! XSPACE = MAKE_RH(PRSLVL(K)/100.0,BK0(I,J,K,L,TEMPRTUR)-273.15,& ! RH
        !                  BK0(I,J,K,L,NUMSTAT+1),0.0)
        BK0(I,J,K,L,NUMSTAT+1) = 0.0 !MAKE_SSH(PRSLVL(K)/100.0,YSPACE,XSPACE,0.0)
      ENDIF
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  PRINT*,'Range of reflectivity and cloud derived bound: ', &
   maxval(BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2,NUMSTAT+1)), &
   minval(BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),2,NUMSTAT+1))

  ! ASSIGN REFLECTIVITY DERIVED BOUNDS TO GRDBKGD0:
  DO I=1,MAXGRID(1)
    IX0 = FLOAT(I-1)/FLOAT(MAXGRID(1)-1)*(FCSTGRD(1)-1)+1
    IX1 = MIN(IX0+1,FCSTGRD(1))
    DO J=1,MAXGRID(2)
      IY0 = FLOAT(J-1)/FLOAT(MAXGRID(2)-1)*(FCSTGRD(2)-1)+1
      IY1 = MIN(IY0+1,FCSTGRD(2))
      DO K=1,MAXGRID(3)
        DO L=1,MAXGRID(4)
          ! Simple shift instead of interpolation:
          GRDBKGD0(I,J,K,L,NUMSTAT+1) = 0.25*(BK0(IX0,IY0,K,L,NUMSTAT+1)+ &
            BK0(IX1,IY0,K,L,NUMSTAT+1)+BK0(IX0,IY1,K,L,NUMSTAT+1)+BK0(IX1,IY1,K,L,NUMSTAT+1))
          GRDBKGD0(I,J,K,L,NUMSTAT+2) = 1000.0 !GRDBKGD0(I,J,K,L,NUMSTAT+1) ! TEST UPPER BOUND
        ENDDO
      ENDDO
    ENDDO
  ENDDO
  PRINT*,'Max dBZ over finest grid: ',maxval(GRDBKGD0(:,:,:,:,NUMSTAT+1)), &
                                      minval(GRDBKGD0(:,:,:,:,NUMSTAT+1))
END SUBROUTINE RDLAPSRDR



SUBROUTINE RDRADROBS
!*************************************************
! READ IN RADAR OBSERVATIONS
! HISTORY: JANUARY 2008, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER ,PARAMETER :: NH=12,NR=7,ML=1500
  CHARACTER(LEN=8)       :: SS
  INTEGER         :: L,NL,IR,DT,NC,NW,O,OS,IS,IP
  REAL            :: HD(NH)
  REAL            :: RA(NR,ML)
  REAL            :: X,Y,P,T,UV,ZZ,TT,UE,MS,EA,AZ,ZE
  REAL            :: OP(NUMDIMS),OB,OE

  REAL            :: X0,Y0,T0,SH,AA,ID,TD
  INTEGER         :: YR,MN,DY,HR
  INTEGER         :: II

  INTEGER     , EXTERNAL :: IREADSB,IREADMG,I4DY

!  CHARACTER*150          :: OD
!  CHARACTER*9            :: A9
!  INTEGER                :: LD,I4,N4,ST

  CHARACTER*80            :: HDRSTR, DATSTR

  INTEGER          :: FG

  DATA HDRSTR / 'SSTN CLAT CLON SELV ANAZ ANEL YEAR MNTH DAYS HOUR MINU MGPT' /
  DATA DATSTR / 'STDM SUPLAT SUPLON HEIT RWND RWAZ RSTD' /


!  CALL GET_SYSTIME(I4,A9,ST)
!  CALL GET_FILESPEC('bufr',2,OD,ST)
!  CALL GET_FILE_TIME(OD,I4,N4)
!  CALL S_LEN(OD,LD)
!  CALL MAKE_FNAM_LP(N4,A9,ST)
!  OD(LD+4:LD+9) = OD(LD-4:LD)
!  OD(LD-5:LD+3) = A9
!  LD = LD+9


  PRINT*, 'READING RADAR BUFR FILE...........'
  OPEN(UNIT=LURAO,FILE='./radarII.bufr',FORM='UNFORMATTED'   &
       ,STATUS='OLD',ACTION='READ')
  CALL OPENBF(LURAO,'IN',LURAO)

  FG=0
  NC=0
  NW=0
!  CALL GET_CONFIG(IS)
!  IF(IS.NE.1)STOP 'LAPS PARAMETERS ARE WRONG!!!'
  O=NALLOBS
  DO WHILE(IREADMG(LURAO,SS,DT).EQ.0)

!    WRITE(6,*) 'READ_RADAR: BUFR FILE DATE IS ',DT

    DO WHILE(IREADSB(LURAO).EQ.0)

      IF(FG .EQ.1) EXIT             ! JUST FOR TEST, BY ZHONGJIE HE

!     READ HEADER. EXTRACT STATION INFORMATION
      CALL UFBINT(LURAO,HD,NH,1,IR,HDRSTR)

!      iyr = hd(7)
!      imo = hd(8)
!      idy = hd(9)
!      ihr = hd(10)
!      imn = hd(11)
!      isc = izero
      T0=HD(10)*60+HD(11) ! BASE TIME   (HOURS)
      X0=HD(3)            ! STATION LON (DEGREES)
      Y0=HD(2)            ! STATION LAT (DEGREES)
      SH=HD(4)            ! STATION ELEVATION
      AA=HD(5)            ! AZIMUTH OF RADIA (DEGREES)
      EA=HD(6)            ! ELEVATION ANGLE    (DEGREES)
      IF(X0.LT.0.0)   X0=X0+360.0
      IF(X0.GE.360.0) X0=X0-360.0

!     Go THROUGH THE DATA LEVELS
      CALL UFBINT(LURAO,RA,NR,ML,NL,DATSTR)

! ===============just for test by zhongjie he
      PRINT*, 'LON=',X0,'LAT=',Y0
      IF(X0-150.GE.118 .AND. X0-150.LE.124 .AND. Y0-20.GE.21. .AND. Y0-20.LE.25 .AND. NL.GE.100) THEN
        FG=1
      ELSE
        CYCLE
      ENDIF

      X_RADAR=X0-150
      Y_RADAR=Y0-150
!===============

      IF(NL>ML) THEN
        WRITE(6,*)'READ_RADAR:  ***ERROR*** INCREASE READ_RADAR BUFR SIZE SINCE', & 
             'NUMBER OF LEVS=', NL,' > MAXLEVS=',ML
        STOP
      ENDIF

      DO L=1,NL
        NW=NW+1

! FOR T LOCATION
        T=(T0+RA(1,L))/60.         ! UNIT IS HOURS
! FOR X LOCATION
        X=RA(3,L)

        X=X-150                    ! JUST FOR A TEST
!        PRINT*, 'X=X-150 IS JUST FOR TEST ----------------------'

        IF(X.LT.0.0)   X=X+360.0
        IF(X.GE.360.0) X=X-360.0        
! FOR Y LOCATION
        Y=RA(2,L)

        Y=Y-20                     ! JUST FOR A TEST
!        PRINT*, 'Y=Y-20 IS JUST FOR TEST ---------------------'

! FOR ZZ OBSERVATION
        ZZ=RA(4,L)
! FOR AZIMUTH OF WIND
        AZ=RA(6,L)

! TRANSFORM
        OP(1)=X
        OP(2)=Y
        IP=1
        CALL HT_TO_PRS(OP(1),OP(2),ZZ,P,IS)
        IF(IS.NE.1)CYCLE
        OP(3)=P
! FOR RADIAL WIND OBSERVATION
        UV=RA(5,L)
        UE=RA(7,L)
! OUTPUT
        ZE=2.0
        OB=UV
        OE=UE
        OS=NUMSTAT+1
        CALL HANDLEOBS(OP,OB,OE,OS,O,IP,AZ,EA)

        NC=NC+1

      ENDDO
    ENDDO
  ENDDO
  NALLOBS=O
  PRINT*,'THE NUMBER OF RADAR DATA IS',NW,'AND',NC,'AVAILABLE'
  CALL CLOSBF(LURAO)
  CLOSE(LURAO)
!  CLOSE(LUNDX)
!  CLOSE(LUNEW)
  PRINT*,'NALLOBS=',NALLOBS
  RETURN
END SUBROUTINE RDRADROBS

SUBROUTINE RDBUFROBS_XIE
!*************************************************
! READ CONVENSIONAL OBS FROM A BUFR FILE, REPLACE
! THE OLD RDBUFROBS DEVELOPED BY WEI LI FOR MORE
! EFFICIENCY.
!
! HISTORY: JAN. 2009 BY YUANFU XIE
!          APR. 2009 BY YUANFU XIE FOR USING LAPS
!          ROUTINE FOR HEIGHT COMPUTATION.
!*************************************************

  IMPLICIT NONE
! --------------------
  INTEGER ,PARAMETER  :: NH=6,NR=16,LN=255
  CHARACTER(LEN=8)    :: SS,SID,TYP
  INTEGER             :: L,NL,IR,DT,NC,NW,O,OS,IS,IP,I,J,K,KK
  INTEGER             :: IDX(2,NUMDIMS)	! INTERPOLATION INDICES
  REAL                :: COE(2,NUMDIMS)	! INTERPOLATION COEFFICENTS
  REAL                :: BD,P1
  REAL(KIND=8)        :: RA(NR,LN),HD(NH)
  REAL                :: X,Y,Z,P,T,DEL,UU,VV,TT,ZZ,ZE,OP(NUMDIMS),DI,SP
  INTEGER, EXTERNAL :: IREADSB,IREADMG,I4DY
  REAL                :: HEIGHT_TO_PRESSURE,RLEVEL_OF_FIELD	! LAPS FUNCTIONS

  REAL                :: AZ,EA            ! ADDED BY ZHONGJIE HE

  ! VARIABLES FOR LAPS OBS: BY YUANFU
  CHARACTER*150       :: OD
  CHARACTER*9         :: A9,WFO_FNAME13_TO_FNAME9
  CHARACTER*13	      :: YYYYMMDD_HHMM
  INTEGER             :: LD,I4,N4,ST
!jhui
  INTEGER :: data_acar
  INTEGER :: tw_u, tw_t, tw_sh, tw_p
  EQUIVALENCE(SID,HD(5))
!jhui

  ! Include a statement function for converting sh to 'rh' = sh/s2r(p) by Yuanfu Xie:
  include 'sh2rh.inc'

  data_acar =0
  tw_u=0
  tw_t=0
  tw_sh=0
  tw_p=0


  CALL GET_FILESPEC('bufr',2,OD,ST)	! GET PATH TO BUFR DIRECTORY
  CALL GET_FILE_TIME(OD,LAPSI4T,N4)	! GET NEAREST I4TIME INTO N4
  CALL S_LEN(OD,LD)
  CALL MAKE_FNAM_LP(N4,A9,ST)
  OD(LD+4:LD+9) = OD(LD-4:LD)
  OD(LD-5:LD+3) = A9
  LD = LD+9

  ! OPEN THE BUFR FILE:
  OPEN(UNIT=LURAO,FILE=OD(1:LD),FORM='UNFORMATTED',STATUS='OLD',ACTION='READ')

  CALL OPENBF(LURAO,'IN',LURAO)	! OPEN BUFR IN CHANNEL

  NC=0		! NUMBER OF VALID OBS
  NW=0		! NUMBER OF OBS IN BUFR

  O=0
  DO WHILE(IREADMG(LURAO,SS,DT).EQ.0)		! OPEN A BUFR MESSAGE
    DO WHILE(IREADSB(LURAO).EQ.0)		! OPEN A SUBSET

      CALL UFBINT(LURAO,HD,NH,1,IR,'XOB YOB ELV DHR SID TYP')	! READ HEADER
!jhui
      IF ( SID =="ACAR    ") THEN
      data_acar = data_acar +1
      END IF
      ! READ OBS:
      CALL UFBINT(LURAO,RA,NR,LN,NL, &
           'XDR YDR PRSS POB POE HRDR UOB VOB WOE ZOB ZOE TOB TOE QOB QOE PMO')

      DO L=1,NL		! FOR ALL LEVELS
        NW=NW+1

        ! FOR X LONGITUDE
        OP(1)=HD(1)
        IF(RA(1,L) .LT. BUFRMISS) OP(1)=RA(1,L)	! BALLOON DRIFT X
        IF(OP(1) .GE. BUFRMISS) CYCLE			! INVALID OBS LONGITUDE
        IF(OP(1) .LT. 0.0) OP(1) = OP(1)+360.0	! EASTWARD LONGITUDE
        ! FOR Y LATITUDE
        OP(2) = HD(2)
        IF(RA(2,L) .LT. BUFRMISS) OP(2) = RA(2,L)	! BALLOON DRIFT Y
        IF(OP(2) .GE. BUFRMISS) CYCLE			! INVALID OBS LATITUDE

        ! CHECK WHETHER THE OBS IN THE HORIZONTAL DOMAIN:
        IF(IF_TEST.NE.1) THEN
          CALL LATLON_TO_RLAPSGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
        ELSE
          CALL OBSTOGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
        ENDIF
        IF(X .LT. 1.0 .OR. Y .LT. 1.0 .OR. & 
	   X .GT. FCSTGRD(1) .OR. Y .GT. FCSTGRD(2) .OR. &
           IS .NE. 1) CYCLE				! OUTSIDE HORIZONTAL DOMAIN

        ! FOR T LOCATION
        ! FIT INTO LAPS TIME FRAMES BY YUANFU:
        I4 = MOD(DT,100)*100
        WRITE(YYYYMMDD_HHMM,1) 20000000+DT/100,'_',I4
1	FORMAT(I8,A1,I4)
        A9 = WFO_FNAME13_TO_FNAME9(YYYYMMDD_HHMM)
        CALL CV_ASC_I4TIME(A9,I4,ST)
        ! IF (I4 .NE. LAPSI4T) THEN
        !   WRITE(*,*) 'Analysis time does not match: ',I4,ITIME2(2)
        !   STOP
        ! ENDIF
        T = I4-ITIME2(1)+3600*HD(4)
        IF (T .LT. T00(1) .OR. T .GT. T00(FCSTGRD(4))) CYCLE	! OUT OF TIME WINDOW

        ! FOR P LOCATION
        P = BUFRMISS
        ! USE OBSERVED PRESSURE IN PASCAL AS BUFR SAVES IN MB:
        IF(RA(4,L).LT. BUFRMISS) P = RA(4,L)*100.0D0	
        ! FOR ZZ OBSERVATION
        ZZ=BUFRMISS					! INITIAL
        IF(RA(10,L).LT. BUFRMISS) ZZ = RA(10,L) 	! USE OBSERVED HEIGHT
        ZE=RA(11,L)					! OBS ERROR IN HEIGHT
        IF(P.GE. BUFRMISS .AND. ZZ .GE. BUFRMISS) CYCLE ! INVALID VERTICAL OBS

        ! INTERPOLATION: BACKGROUND TO OBS:
        IDX(1,1) = INT(X)
        IDX(1,2) = INT(Y)
        IDX(2,1:2) = MIN0(IDX(1,1:2)+1,FCSTGRD(1:2))
        COE(1,1) = IDX(2,1)-X
        COE(1,2) = IDX(2,2)-Y
        COE(2,1:2) = 1.0-COE(1,1:2)	! ASSUME HORIZONTAL GRID DISTANCE IS 1
        IF (FCSTGRD(4) .EQ. 1) THEN
          PRINT*,'RDBUFROBS -- ERROR: TEMPORAL GRID HAS ONE GRIDPOINT!'
          STOP
        ENDIF
        DEL = (ITIME2(2)-ITIME2(1))/(FCSTGRD(4)-1)	! DELTA TIME
        IDX(1,4) = INT(T/DEL)+1
        IDX(2,4) = MIN0(IDX(1,4)+1,FCSTGRD(4))
        COE(2,4) = T/DEL+1-IDX(1,4)
        COE(1,4) = 1.0-COE(2,4)

        ! FOR MISSING PRESSURE, CONVERT FROM BACKGROUND HEIGHT:
        IP=0		! DIRECT PRESSURE OBS; 1 BECOMES HYDROSTATIC OBS

        ! FILL MSLP IF AVAILABLE:
        IF (ZZ .LE. 10.0 .AND. RA(16,L) .NE. BUFRMISS) P = RA(16,L)*100.0D0

        IF(P .GE. BUFRMISS) THEN	! USE HYDROSTATIC RELATION	
          IP=1		! HYDROSTATIC OBS
          P = HEIGHT_TO_PRESSURE(ZZ,BK0(1,1,1,IDX(1,4),PRESSURE),P00,&
                    FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),IDX(1,1),IDX(1,2))
        ENDIF

        ! FIND THE VERTICAL LEVEL:
        KK = RLEVEL_OF_FIELD(P,P00,1,1,FCSTGRD(3),1,1,IS)
        IF (IS .NE. 1) CYCLE	! OUT OF VERTICAL DOMAIN

        ! VERTICAL INTERPOLATION:
        IF (KK .EQ. FCSTGRD(3)) THEN
          IDX(1,3) = KK-1
          COE(1,3) = 0.0
        ELSE
          IDX(1,3) = KK
          COE(1,3) = (P-P00(KK+1))/(P00(KK)-P00(KK+1))
        ENDIF
        IDX(2,3) = IDX(1,3)+1
        COE(2,3) = 1.0-COE(1,3)
        ! VERTICAL GRID POSITIN:
        Z = IDX(1,3)+COE(2,3)

        ! FOR WIND OBSERVATION
        IF (RA(9,L) .LT. BUFRMISS .AND. &
          RA(7,L) .LT. BUFRMISS .AND. RA(8,L) .LT. BUFRMISS) THEN

          ! SAVE OBS INTO DATA ARRAYS:
          O = O+2	! TOTAL SINGLE OBS COUNT
          NST(U_CMPNNT) = NST(U_CMPNNT)+1
          NST(V_CMPNNT) = NST(V_CMPNNT)+1	! WIND OBS COUNTS
          OBP(1,NST(U_CMPNNT),U_CMPNNT) = X-1
          OBP(2,NST(U_CMPNNT),U_CMPNNT) = Y-1
          OBP(3,NST(U_CMPNNT),U_CMPNNT) = Z-1
          OBP(4,NST(U_CMPNNT),U_CMPNNT) = T/(T00(2)-T00(1)) ! WIND OBS LOCATIONS
          OBP(1,NST(V_CMPNNT),V_CMPNNT) = X-1
          OBP(2,NST(V_CMPNNT),V_CMPNNT) = Y-1
          OBP(3,NST(V_CMPNNT),V_CMPNNT) = Z-1
          OBP(4,NST(V_CMPNNT),V_CMPNNT) = T/(T00(2)-T00(1))
!jhui
          IF ( OBP(1,NST(U_CMPNNT),U_CMPNNT) .GT. 38 .AND. &
               OBP(1,NST(U_CMPNNT),U_CMPNNT) .LT. 112 .AND. &
               OBP(2,NST(U_CMPNNT),U_CMPNNT) .GT. 33 .AND. &
               OBP(2,NST(U_CMPNNT),U_CMPNNT) .LT. 107 ) THEN
          tw_u = tw_u + 1
          ENDIF

          ! SAVE INCREMENT:
          UU = 0.0
          VV = 0.0
          DO KK=1,2
            DO K=1,2
              DO J=1,2
                DO I=1,2
                  UU = UU+COE(I,1)*COE(J,2)*COE(K,3)*COE(KK,4)* &
                       BK0(IDX(I,1),IDX(J,2),IDX(K,3),IDX(KK,4),U_CMPNNT)
                  VV = VV+COE(I,1)*COE(J,2)*COE(K,3)*COE(KK,4)* &
                       BK0(IDX(I,1),IDX(J,2),IDX(K,3),IDX(KK,4),V_CMPNNT)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          OBS(NST(U_CMPNNT),U_CMPNNT) =RA(7,L)-UU
          OBS(NST(V_CMPNNT),V_CMPNNT) =RA(8,L)-VV
          OBE(NST(U_CMPNNT),U_CMPNNT) =0.5	!RA(9,L)	! WIND ERROR
          OBE(NST(V_CMPNNT),V_CMPNNT) =0.5	!RA(9,L)
!jhui	
!          OBE(NST(U_CMPNNT),U_CMPNNT) =0.3	!RA(9,L)	! WIND ERROR
!          OBE(NST(V_CMPNNT),V_CMPNNT) =0.3	!RA(9,L)	
          ! SAVE OBS INTO PIG FILE OF LAPS AS WIND DIRECTION AND SPEED:
          ! UU = RA(7,L)
          ! VV = RA(8,L)
          ! CALL UV_TO_DISP(UU,VV,DI,SP)
          ! WRITE(PIGOBS_CHANNEL,*) X-1,Y-1,Z-1,DI,SP,SID

          ! Add a threshold check: Yuanfu June 2010
          IF ((ABS(OBS(NST(U_CMPNNT),U_CMPNNT)) .GT. 20.0) .OR. &
              (ABS(OBS(NST(V_CMPNNT),V_CMPNNT)) .GT. 20.0) ) THEN
            ! Over the threshold, remove this data:
            O = O-2
            NST(U_CMPNNT) = NST(U_CMPNNT)-1
            NST(V_CMPNNT) = NST(V_CMPNNT)-1
          ENDIF
        ENDIF

        ! FOR TEMPERATURE OBSERVATION
        IF (RA(12,L) .LT. BUFRMISS .AND. RA(13,L) .LT. BUFRMISS ) THEN
          ! SAVE OBS INTO DATA ARRAYS:
          O = O+1	! TOTAL SINGLE OBS COUNT
          NST(TEMPRTUR) = NST(TEMPRTUR)+1	! TEMPERATURE OBS COUNTS
          OBP(1,NST(TEMPRTUR),TEMPRTUR) = X-1
          OBP(2,NST(TEMPRTUR),TEMPRTUR) = Y-1
          OBP(3,NST(TEMPRTUR),TEMPRTUR) = Z-1
          OBP(4,NST(TEMPRTUR),TEMPRTUR) = T/(T00(2)-T00(1)) ! TEMPERATUER OBS LOCATION
          OBS(NST(TEMPRTUR),TEMPRTUR) =RA(12,L)+273.15D0	! IN KELVIN
          OBE(NST(TEMPRTUR),TEMPRTUR) =0.5	!RA(13,L)			! OBS ERROR
!jhui
!          OBE(NST(TEMPRTUR),TEMPRTUR) =0.3	!RA(13,L)			! OBS ERROR
!jhui
          IF ( OBP(1,NST(TEMPRTUR),TEMPRTUR) .GT. 38 .AND. &
               OBP(1,NST(TEMPRTUR),TEMPRTUR) .LT. 112 .AND. &
               OBP(2,NST(TEMPRTUR),TEMPRTUR) .GT. 33 .AND. &
               OBP(2,NST(TEMPRTUR),TEMPRTUR) .LT. 107 ) THEN
          tw_t = tw_t + 1
          ENDIF

          ! SAVE OBS INNOVATION:
          TT = 0.0
          DO KK=1,2
            DO K=1,2
              DO J=1,2
                DO I=1,2
                  TT = TT+COE(I,1)*COE(J,2)*COE(K,3)*COE(KK,4)* &
                       BK0(IDX(I,1),IDX(J,2),IDX(K,3),IDX(KK,4),TEMPRTUR)
                ENDDO
              ENDDO
            ENDDO
          ENDDO

          OBS(NST(TEMPRTUR),TEMPRTUR) = OBS(NST(TEMPRTUR),TEMPRTUR)-TT

          ! Add a threshold check: Yuanfu June 2010
          IF (ABS(OBS(NST(TEMPRTUR),TEMPRTUR)) .GT. 10.0) THEN
            ! Over the threshold, remove this data:
            O = O-1
            NST(TEMPRTUR) = NST(TEMPRTUR)-1
          ENDIF

        ENDIF

        ! FOR SPECIFIC HUMIDITY OBSERVATION
        IF (RA(14,L) .LT. BUFRMISS .AND. RA(15,L) .LT. BUFRMISS ) THEN
          ! SAVE OBS INTO DATA ARRAYS:
          O = O+1	! TOTAL SINGLE OBS COUNT
          NST(HUMIDITY) = NST(HUMIDITY)+1	! HUMIDITY OBS COUNTS
          OBP(1,NST(HUMIDITY),HUMIDITY) = X-1
          OBP(2,NST(HUMIDITY),HUMIDITY) = Y-1
          OBP(3,NST(HUMIDITY),HUMIDITY) = Z-1
          OBP(4,NST(HUMIDITY),HUMIDITY) = T/(T00(2)-T00(1)) ! HUMIDITY OBS LOCATION

!jhui
          IF ( OBP(1,NST(HUMIDITY),HUMIDITY) .GT. 38 .AND. &
               OBP(1,NST(HUMIDITY),HUMIDITY) .LT. 112 .AND. &
               OBP(2,NST(HUMIDITY),HUMIDITY) .GT. 33 .AND. &
               OBP(2,NST(HUMIDITY),HUMIDITY) .LT. 107 ) THEN
          tw_sh = tw_sh + 1
          ENDIF
          ! SAVE OBS INNOVATION:
          TT = 0.0
          DO KK=1,2
            DO K=1,2
              DO J=1,2
                DO I=1,2
                  TT = TT+COE(I,1)*COE(J,2)*COE(K,3)*COE(KK,4)* &
                       BK0(IDX(I,1),IDX(J,2),IDX(K,3),IDX(KK,4),HUMIDITY)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          ! Convert SH obs to 'RH' = SH/s2r(p) by Yuanfu Xie: Dec. 2013
          OBS(NST(HUMIDITY),HUMIDITY) =RA(14,L)/1000.0/s2r(p/100.0)-TT	! HUMIDITY OBS: BUFR USES MG/KG
          OBE(NST(HUMIDITY),HUMIDITY) = 1.0	!RA(15,L)	! OBS ERROR

          ! Add a threshold check: Yuanfu June 2010
          IF (ABS(OBS(NST(HUMIDITY),HUMIDITY)) .GT. 10.0) THEN
            ! Over the threshold, remove this data:
            O = O-1
            NST(HUMIDITY) = NST(HUMIDITY)-1
          ENDIF

        ENDIF
        ! FOR HEIGHT OBSERVATION: IF PRESSURE IS NOT DERIVED FROM HEIGHT:
        IF (IP .EQ. 0 .AND. ZZ .LT. BUFRMISS .AND. ZE .LT. BUFRMISS ) THEN

          ! REMOVE ALL HIGH LEVEL HEIGHT OBS:
          !IF (L .GT. 3) CYCLE

          ! SAVE OBS INTO DATA ARRAYS:
          O = O+1	! TOTAL SINGLE OBS COUNT
          NST(PRESSURE) = NST(PRESSURE)+1	! HEIGHT OBS COUNTS
          OBP(1,NST(PRESSURE),PRESSURE) = X-1
          OBP(2,NST(PRESSURE),PRESSURE) = Y-1
          OBP(3,NST(PRESSURE),PRESSURE) = Z-1
          OBP(4,NST(PRESSURE),PRESSURE) = T/(T00(2)-T00(1)) ! HEIGHT OBS LOCATION

!jhui
          IF ( OBP(1,NST(PRESSURE),PRESSURE) .GT. 38 .AND. &
               OBP(1,NST(PRESSURE),PRESSURE) .LT. 112 .AND. &
               OBP(2,NST(PRESSURE),PRESSURE) .GT. 33 .AND. &
               OBP(2,NST(PRESSURE),PRESSURE) .LT. 107 ) THEN
          tw_p = tw_p + 1
          ENDIF

          ! SAVE OBS INNOVATION:
          TT = 0.0
          DO KK=1,2
            DO K=1,2
              DO J=1,2
                DO I=1,2
                  TT = TT+COE(I,1)*COE(J,2)*COE(K,3)*COE(KK,4)* &
                       BK0(IDX(I,1),IDX(J,2),IDX(K,3),IDX(KK,4),PRESSURE)
                ENDDO
              ENDDO
            ENDDO
          ENDDO
          OBS(NST(PRESSURE),PRESSURE) =ZZ-TT		! HEIGHT OBS
          OBE(NST(PRESSURE),PRESSURE) = 2.0	!ZE	! OBS ERROR

          ! Add a threshold check: Yuanfu June 2010
          IF (ABS(ZZ-TT) .GT. 50) THEN ! Use 50m as the threshold value now.
            ! Over the threshold, remove this data:
            O = O-1
            NST(PRESSURE) = NST(PRESSURE)-1
          ENDIF

        ENDIF

        NC=NC+1			! NUMBER OF VALID OBS COUNTS

      ENDDO
    ENDDO
  ENDDO
  NALLOBS=O
  PRINT*,'THE NUMBER OF LOCATION IS',NW,'AND',NC,'AVAILABLE'
  PRINT*,'RDBUFROBS: NUMBER OF ALL OBS = ',NALLOBS
  CALL CLOSBF(LURAO)
  CLOSE(LURAO)
!  CLOSE(LUNDX)
!  CLOSE(LUNEW)

   
  RETURN
END SUBROUTINE RDBUFROBS_XIE

SUBROUTINE HANDLEOBS(OP,OB,OE,OS,O,IP,AZ,EA)
!*************************************************
! CALCULATE THE DIFFERENCE BETWEEN OBSERVATION AND BACKGROUND
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,M,N,NP(MAXDIMS),NN(MAXDIMS),IS,O,OS,IP
  REAL  :: X,Y,P
  REAL  :: DI,SP	! YUANFU FOR OUTPUT PIG FILE
  REAL  :: AC(NUMDIMS,NGPTOBS),OC(NUMDIMS),CO(NGPTOBS),HT
  REAL  :: OP(NUMDIMS),OB,OE

  INTEGER  :: UU,VV
  REAL     :: AZ,TU,TV,EA               ! ADDED BY ZHONGJIE HE
  REAL     :: D2R                       ! CONVERTER FROM DEG TO RADIAN

  D2R = 3.14159/180.0

! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
  CALL LATLON_TO_RLAPSGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
!  CALL OBSTOGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)

  IF(X.LT.1.0.OR.Y.LT.1.0.OR.X.GT.FCSTGRD(1).OR.Y.GT.FCSTGRD(2).OR.IS.NE.1)RETURN
  CALL VRTCLPSTN(FCSTGRD(3),LIMIT_3,P00,OP(3),P,IS)

  IF(IS.NE.1)RETURN
  DO N=1,MAXDIMS
    NP(N)=1
  ENDDO
  NP(1)=INT(X)
  NP(2)=INT(Y)
  NP(3)=INT(P)
  DO N=1,MAXDIMS
    IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
  ENDDO
  OC(1)=X
  OC(2)=Y
  OC(3)=P
  M=0
!==============================================================
!  DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
!  DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
!  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
!  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
!    NN(1)=I
!    NN(2)=J
!    NN(3)=K
!    NN(4)=T
!    M=M+1
!    DO N=1,NUMDIMS
!      AC(N,M)=NN(N)*1.0D0
!    ENDDO
!  ENDDO
!  ENDDO
!  ENDDO
!  ENDDO
!======================= modified by zhongjie he ==============
  DO T=NP(4),NP(4)+1
  DO K=NP(3),NP(3)+1
  DO J=NP(2),NP(2)+1
  DO I=NP(1),NP(1)+1
    NN(1)=MIN0(I,FCSTGRD(1))
    NN(2)=MIN0(J,FCSTGRD(2))
    NN(3)=MIN0(K,FCSTGRD(3))
    NN(4)=MIN0(T,FCSTGRD(4))
    M=M+1
    DO N=1,NUMDIMS
      AC(N,M)=NN(N)*1.0D0
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  ENDDO
!==============================================================
  CALL INTERPLTN(NUMDIMS,NGPTOBS,CO,AC,OC)
  !CALL INTERPLTN_XIE(NUMDIMS,NGPTOBS,CO,AC,OC,3,NUMGRID(3),PPM)
  HT=0.0
  M=0
  TU=0.0
  TV=0.0
  IF(OS.LE.NUMSTAT) THEN
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      HT=HT+CO(M)*BK0(I,J,K,T,OS)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ELSEIF(OS.EQ.NUMSTAT+1) THEN        ! FOR RADAR DATA , BY ZHONGJIE HE
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      TU=TU+CO(M)*BK0(I,J,K,T,UU)
      TV=TV+CO(M)*BK0(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=(TU*SIN(D2R*AZ)+TV*COS(D2R*AZ))*COS(D2R*EA)
  ELSEIF(OS.EQ.NUMSTAT+2) THEN        ! FOR SFMR DATA , BY ZHONGJIE HE
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      TU=TU+CO(M)*BK0(I,J,K,T,UU)
      TV=TV+CO(M)*BK0(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT(TU*TU+TV*TV)
  ENDIF

  IF(OS.EQ.4.AND.ABS(OB-HT).GE.8.0)RETURN
  IF(OS.EQ.3.AND.ABS(OB-HT).GE.50.0)RETURN
  IF(IP.EQ.1.AND.OS.EQ.3)RETURN
  OB=OB-HT

  CALL VRTCLPSTN8(MAXGRID(3),LIMIT_3,PP0,OP(3),P,IS)

  IF(IS.NE.1)RETURN

  OC(3)=P
  O=O+1
  NST(OS)=NST(OS)+1
  IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
  DO N=1,NUMDIMS
    OBP(N,NST(OS),OS)=OC(N)-1.0D0
  ENDDO
  OBS(NST(OS),OS)=OB
  OBE(NST(OS),OS)=OE
  IF(OS.EQ.NUMSTAT+1) THEN
    OBA(NST(OS),1)=AZ
    OBA(NST(OS),2)=EA
  ENDIF


  RETURN
END SUBROUTINE HANDLEOBS

SUBROUTINE HT_TO_PRS(X0,Y0,Z,P,IS)
!*************************************************
! CONVERT HEIGHT TO PRESSURE TO USE AS MANY DATA AS POSSIBLE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: IS,I,J,K,T,N,M,NP(2),ZZ
  INTEGER  :: LM(2)
  REAL  :: X,Y,X0,Y0
  REAL  :: Z,P,AC(2,4),OC(2),CO(4),HT(FCSTGRD(3))
! --------------------
  ZZ=PRESSURE                                                        !  BY ZHONGJIE HE
  CALL LATLON_TO_RLAPSGRID(Y0,X0,Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
!   CALL OBSTOGRID(Y0,X0,Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)        !  BY ZHONGJIE HE

  IF(X.LT.1.0.OR.Y.LT.1.0.OR.X.GT.FCSTGRD(1).OR.Y.GT.FCSTGRD(2).OR.IS.NE.1)THEN
    IS=0
    RETURN
  ENDIF
  OC(1)=X
  OC(2)=Y
  NP(1)=INT(X)
  NP(2)=INT(Y)
  DO N=1,2
    IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
  ENDDO
  M=0
!==============================================
!  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
!  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
!    M=M+1
!    AC(1,M)=I*1.0D0
!    AC(2,M)=J*1.0D0
!  ENDDO
!  ENDDO
!============== modified by zhongjie he =======
  DO N=1,2
    LM(N)=NP(N)+1
    IF(NUMDIMS.LT.N) LM(N)=NP(N)
  ENDDO
  DO J=NP(2),LM(2)
  DO I=NP(1),LM(1)
    M=M+1
    AC(1,M)=MIN0(I,FCSTGRD(1))*1.0
    AC(2,M)=MIN0(J,FCSTGRD(2))*1.0
  ENDDO
  ENDDO
!==============================================
  CALL INTERPLTN(2,4,CO,AC,OC)
  !CALL INTERPLTN_XIE(2,4,CO,AC,OC,3,NUMGRID(3),PPM)
  T=1
  DO K=1,FCSTGRD(3)
    HT(K)=0.0
    M=0
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      HT(K)=HT(K)+CO(M)*BK0(I,J,K,T,ZZ)
    ENDDO
    ENDDO
  ENDDO
  CALL ZPCONVERT(FCSTGRD(3),HT,P00,Z,P,IS)
  RETURN
END SUBROUTINE HT_TO_PRS


SUBROUTINE RDOBSTEST
!*************************************************
! READ IN TEST DATA
! HISTORY: FEBRUARY 2008, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  CHARACTER*8     :: SID
  INTEGER         :: L,O,OS,IS,IP
  REAL            :: X,Y,P,UV,ZZ,UE,EA,AZ,ZE,T
  REAL            :: OP(NUMDIMS),OB,OE
  INTEGER         :: NC,NW,LMAX
  REAL            :: MF                 ! NUMBER FLAG OF MISSED DATA

  MF=-9999.0
  NC=0
  NW=0
  SID='TEST OBS'
  O=NALLOBS
  OPEN(11,FILE='OBS_RADAR_CNVTN.DAT',ACTION='READ')
!  OPEN(11,FILE='OBS_RADAR.DAT',ACTION='READ')
    X_RADAR=350.
    Y_RADAR=150
    READ(11,*) LMAX
    DO L=1,LMAX
      NW=NW+1
      READ(11,'(8F10.3,I3,F12.3)') X,Y,ZZ,T,AZ,EA,UV,UE,OS,P
!      READ(11,*) X,Y,ZZ,T,AZ,EA,UV,UE,OS,P
!      READ(11,'(7F10.3,I3,F12.3)') X,Y,ZZ,AZ,EA,UV,UE,OS,P      !  FOR PRESURE COORDINATE AND P IS CONTENTED IN THE FILE OF OBS
      IF(ABS(X-MF).LE.1E-5) CYCLE
      IF(ABS(Y-MF).LE.1E-5) CYCLE
      IF(ABS(AZ-MF).LE.1E-5) CYCLE
      IF(ABS(EA-MF).LE.1E-5) CYCLE
      IF(ABS(UV-MF).LE.1E-5) CYCLE

!================= for test ===============
!      if(mod(l,11).ne.1) cycle
!      if(os.eq.numstat+1 ) cycle
!================= for test ===============

! TRANSFORM
      OP(1)=X
      OP(2)=Y
      IP=0
      IF(OS.EQ.NUMSTAT+1) IP=1

      IF(IFPCDNT.EQ.1) THEN
        IF(ABS(P-MF).LE.1E-5) THEN
          CALL HT_TO_PRS(OP(1),OP(2),ZZ,P,IS)
          IP=1
          IF(IS.NE.1)CYCLE
        ENDIF
        OP(3)=P
        IF(NUMDIMS.GE.4) OP(4)=T
        ZE=2.0
        OB=UV
        OE=UE
        CALL HANDLEOBS_SIGMA(OP,OB,OE,OS,O,IP,AZ,EA,SID)
        NC=NC+1

        OB=ZZ
        OE=ZE
        OS=PRESSURE
        CALL HANDLEOBS_SIGMA(OP,OB,OE,OS,O,IP,AZ,EA,SID)
        NC=NC+1
      ELSE
        OP(3)=ZZ
        IF(NUMDIMS.GE.4) OP(4)=T
        ZE=2.0
        OB=UV
        OE=UE
        CALL HANDLEOBS_SIGMA(OP,OB,OE,OS,O,IP,AZ,EA,SID)
        NC=NC+1
      ENDIF

    ENDDO
    NALLOBS=O
  CLOSE(11)

  DO OS=1,NUMSTAT+2
    IF(OS.EQ.U_CMPNNT) PRINT*, 'THE NUMBER OF OBSERVED U DATA IS:',NST(OS)
    IF(OS.EQ.V_CMPNNT) PRINT*, 'THE NUMBER OF OBSERVED V DATA IS:',NST(OS)
    IF(OS.EQ.W_CMPNNT) PRINT*, 'THE NUMBER OF OBSERVED W DATA IS:',NST(OS)
    IF(OS.EQ.PRESSURE) PRINT*, 'THE NUMBER OF OBSERVED PRESSURE DATA IS:',NST(OS)
    IF(OS.EQ.TEMPRTUR) PRINT*, 'THE NUMBER OF OBSERVED TEMPRTUR DATA IS:',NST(OS)
    IF(OS.EQ.NUMSTAT+1 .AND. NST(NUMSTAT+1).GE.1) PRINT*, 'THE NUMBER OF RADAR DATA IS:',NST(OS)
    IF(OS.EQ.NUMSTAT+2 .AND. NST(NUMSTAT+2).GE.1) PRINT*, 'THE NUMBER OF SFMR DATA IS:',NST(OS)
  ENDDO

  RETURN
END SUBROUTINE RDOBSTEST

SUBROUTINE HANDLEOBS_SIGMA(OP,OB,OE,OS,O,IP,AZ,EA,SID)
!*************************************************
! CALCULATE THE DIFFERENCE BETWEEN OBSERVATION AND BACKGROUND
! HISTORY: MARCH 2007, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  CHARACTER*8, INTENT(IN) :: SID	! STATION NAME
  INTEGER  :: I,J,K,T,M,N,NP(MAXDIMS),NN(MAXDIMS),IS,O,OS,IP
  INTEGER  :: LM(MAXDIMS)
  REAL  :: X,Y,P,TM
  REAL  :: DI,SP		! YUANFU FOR OUTPUT PIG FILE
  REAL  :: AC(NUMDIMS,NGPTOBS),OC(NUMDIMS),CO(NGPTOBS),HT
  REAL  :: OP(NUMDIMS),OB,OE

  INTEGER  :: UU,VV,WW
  REAL     :: AZ,TU,TV,EA,HS,HE,TW
  REAL     :: D2R
! --------------------

  D2R = 3.14159/180.0

  UU=U_CMPNNT
  VV=V_CMPNNT
  WW=W_CMPNNT

  IF(IP.EQ.1.AND.OS.EQ.3)RETURN

  IF(IF_TEST.NE.1) THEN
    CALL LATLON_TO_RLAPSGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
  ELSE
    CALL OBSTOGRID(OP(2),OP(1),Y00,X00,FCSTGRD(1),FCSTGRD(2),X,Y,IS)
  ENDIF

  IF(X.LT.1.0.OR.Y.LT.1.0.OR.X.GT.FCSTGRD(1).OR.Y.GT.FCSTGRD(2).OR.IS.NE.1)RETURN

  IF(IFPCDNT.EQ.1) THEN                      ! FOR PRESSURE COORDINATE
    CALL VRTCLPSTN(FCSTGRD(3),LIMIT_3,P00,OP(3),P,IS)
  ELSEIF(IFPCDNT.EQ.2) THEN                  ! FOR HEIGHT COORDINATE
    CALL VRTCLPSTN(FCSTGRD(3),LIMIT_3,Z00,OP(3),P,IS)
  ELSEIF(IFPCDNT.EQ.0) THEN                  ! FOR SIGMA COORDINATE
    DO N=1,2
      NP(N)=1
    ENDDO
    NP(1)=INT(X)
    NP(2)=INT(Y)
    NP(3)=1
    NP(4)=1
    DO N=1,MAXDIMS
      IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
    ENDDO
    OC(1)=X
    OC(2)=Y
    OC(3)=1
    OC(4)=1
    M=0
!============================================================
!    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
!    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
!    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
!    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
!      NN(1)=I
!      NN(2)=J
!      NN(3)=K
!      NN(4)=T
!      M=M+1
!      DO N=1,NUMDIMS
!        AC(N,M)=NN(N)*1.0D0
!      ENDDO
!    ENDDO
!    ENDDO
!    ENDDO
!    ENDDO
!=====================MODIFIED BY ZHONGJIE HE ==============
    DO N=1,MAXDIMS
      LM(N)=NP(N)+1
      IF(NUMDIMS.LT.N) LM(N)=NP(N)
    ENDDO
    DO T=NP(4),LM(4)
    DO K=NP(3),LM(3)
    DO J=NP(2),LM(2)
    DO I=NP(1),LM(1)
      NN(1)=MIN0(I,FCSTGRD(1))
      NN(2)=MIN0(J,FCSTGRD(2))
      NN(3)=MIN0(K,FCSTGRD(3))
      NN(4)=MIN0(T,FCSTGRD(4))
      M=M+1
      DO N=1,NUMDIMS
        AC(N,M)=NN(N)*1.0D0
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
!============================================================
    ! CALL INTERPLTN_XIE(NUMDIMS,NGPTOBS,CO,AC,OC,3,NUMGRID(3),PPM)
    CALL INTERPLTN(NUMDIMS,NGPTOBS,CO,AC,OC)
    M=0
    HS=0.0
    HE=0.0
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      HS=HS+CO(M)*HEIGHTL(I,J)
      HE=HE+CO(M)*HEIGHTU(I,J)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    OP(3)=(OP(3)-HS)/(HE-HS)
    CALL VRTCLPSTN(FCSTGRD(3),LIMIT_3,Z00,OP(3),P,IS)
  ENDIF
!jhui
!  FOR RADAR INGEST DATA that outside 0-3600
  IF ( OS .GT. NUMSTAT) THEN
  IF ( OP(4).GT. 3600)  OP(4) = 3600 
  IF ( OP(4).LT. 0)  OP(4) = 0 
  IF(NUMDIMS.GE.4) CALL VRTCLPSTN(FCSTGRD(4),LIMIT_4,T00,OP(4),TM,IS)
!  ELSE
!  IF(NUMDIMS.GE.4) CALL VRTCLPSTN(FCSTGRD(4),LIMIT_4,T00,OP(4),TM,IS)
  ENDIF

  IF(IS.NE.1)RETURN
  DO N=1,MAXDIMS
    NP(N)=1
  ENDDO
  NP(1)=INT(X)
  NP(2)=INT(Y)
  NP(3)=INT(P)
  IF(NUMDIMS.GE.4) NP(4)=INT(TM)
  DO N=1,MAXDIMS
    IF(NP(N).EQ.FCSTGRD(N).AND.FCSTGRD(N).NE.1)NP(N)=FCSTGRD(N)-1
  ENDDO
  OC(1)=X
  OC(2)=Y
  OC(3)=P
  IF(NUMDIMS.GE.4) OC(4)=TM
  M=0
!===================================================
!  DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
!  DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
!  DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
!  DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
!    NN(1)=I
!    NN(2)=J
!    NN(3)=K
!    NN(4)=T
!    M=M+1
!    DO N=1,NUMDIMS
!      AC(N,M)=NN(N)*1.0D0
!    ENDDO
!  ENDDO
!  ENDDO
!  ENDDO
!  ENDDO
!=================== modified by zhongjie he========
  DO N=1,MAXDIMS
    LM(N)=NP(N)+1
    IF(NUMDIMS.LT.N) LM(N)=NP(N)
  ENDDO
  DO T=NP(4),LM(4)
  DO K=NP(3),LM(3)
  DO J=NP(2),LM(2)
  DO I=NP(1),LM(1)
    NN(1)=MIN0(I,FCSTGRD(1))
    NN(2)=MIN0(J,FCSTGRD(2))
    NN(3)=MIN0(K,FCSTGRD(3))
    NN(4)=MIN0(T,FCSTGRD(4))
    M=M+1
    DO N=1,NUMDIMS
      AC(N,M)=NN(N)*1.0D0
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  ENDDO
!===================================================
  ! CALL INTERPLTN_XIE(NUMDIMS,NGPTOBS,CO,AC,OC,3,NUMGRID(3),PPM)
  CALL INTERPLTN(NUMDIMS,NGPTOBS,CO,AC,OC)
  HT=0.0
  M=0
  TU=0.0
  TV=0.0
  TW=0.0
  IF(OS.LE.NUMSTAT) THEN
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      HT=HT+CO(M)*BK0(I,J,K,T,OS)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ELSEIF(OS.EQ.NUMSTAT+1) THEN        ! FOR RADAR DATA , BY ZHONGJIE HE
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      TU=TU+CO(M)*BK0(I,J,K,T,UU)
      TV=TV+CO(M)*BK0(I,J,K,T,VV)
      IF(WW.NE.0) TW=TW+CO(M)*BK0(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=(TU*SIN(D2R*AZ)+TV*COS(D2R*AZ))*COS(D2R*EA)+TW*SIN(D2R*EA)
  ELSEIF(OS.EQ.NUMSTAT+2) THEN        ! FOR SFMR DATA , BY ZHONGJIE HE
    DO T=NP(4),MIN0(NP(4)+1,FCSTGRD(4))
    DO K=NP(3),MIN0(NP(3)+1,FCSTGRD(3))
    DO J=NP(2),MIN0(NP(2)+1,FCSTGRD(2))
    DO I=NP(1),MIN0(NP(1)+1,FCSTGRD(1))
      M=M+1
      TU=TU+CO(M)*BK0(I,J,K,T,UU)
      TV=TV+CO(M)*BK0(I,J,K,T,VV)
      IF(WW.NE.0) TW=TW+CO(M)*BK0(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT(TU*TU+TV*TV+TW*TW)
  ELSEIF(OS.EQ.NUMSTAT+3) THEN        
    HT =0.0
  ENDIF

  IF(OS.EQ.TEMPRTUR.AND.ABS(OB-HT).GE.8.0)RETURN
  IF(OS.EQ.PRESSURE.AND.ABS(OB-HT).GE.50.0)RETURN
  OB=OB-HT

  IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2) THEN              ! FOR SIGMA AND HEIGHT COORDINATE
    CALL VRTCLPSTN8(MAXGRID(3),LIMIT_3,ZZB,OP(3),P,IS)
  ELSE                               ! FOR PRESURE COORDINATE
    CALL VRTCLPSTN8(MAXGRID(3),LIMIT_3,PP0,OP(3),P,IS)
  ENDIF

  IF(IS.NE.1)RETURN

  OC(3)=P
  O=O+1
  NST(OS)=NST(OS)+1
  IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
  DO N=1,NUMDIMS
    OBP(N,NST(OS),OS)=OC(N)-1.0D0
  ENDDO
  OBS(NST(OS),OS)=OB
  OBE(NST(OS),OS)=OE
  IF(OS.EQ.NUMSTAT+1) THEN
    OBA(NST(OS),1)=AZ
    OBA(NST(OS),2)=EA
  ENDIF
  
  RETURN
END SUBROUTINE HANDLEOBS_SIGMA

END MODULE READOBSERVES
