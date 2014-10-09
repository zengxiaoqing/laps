MODULE INPUT_BG_OBS

  USE PRMTRS_STMAS
  USE GENERALTOOLS
  USE READOBSERVES, ONLY : RDRADROBS, RDBUFROBS, RDBUFROBS_XIE, ALLOCATOB, &
             DEALOCTOB, NOBSMAX, OBP, OBS, OBE, NST, OBA, RDOBSTEST, RDLAPSRDR
  USE READ_BACKGRD, ONLY : RDLAPSBKG, RDBCKGRND, ALLOCTBKG, DEALCTBKG, BK0, C00, D00, X00, Y00, P00,  &
                           DX0, DY0, Z00, DZ0, DT0, HEIGHTU, HEIGHTL, RDBKGTEST

  ! IPW data ingest:
  USE lapsdata

  IMPLICIT NONE
  INTEGER  ,ALLOCATABLE :: GRIDMASK(:,:,:,:,:)
  REAL     ,ALLOCATABLE :: OBRADIUS(:,:)
  
  ! For reading in IPW data namelist:
  TYPE(lapsData_t) :: ipws
  CLASS(lapsData_t), POINTER :: idat

  PUBLIC           BKGRNDOBS
  PRIVATE          READNMLST, BKGMEMALC, OBSMEMALC, MEMORYALC, MEMORYRLS, ADDBKGRND, GETBKGRND_NEW

!**************************************************
!COMMENT:
!   THIS MODULE IS USED BY MAIN.F90 TO GET INFORMATIONS OF BACKGROUND AND OBSERVATIONS.
!   SUBROUTINES:
!      BKGRNDOBS: MAIN ROUTINE OF THIS MODULE, GET IN INFORMATIONS OF BACKGROUND AND OBSERVATIONS.
!      READNMLST: READ NAME LIST FILE. 
!      BKGMEMALC: ALLOCATE MEMORY FOR BACKGROUND ARRAYS.
!      OBSMEMALC: ALLOCATE OBSERVATION MEMORY AND SAVE THE ABSERVATION INFORMATIONS.
!      MEMORYALC: GET INFORMATIONS FROM THE MODEL OR SOME INITIAL FILES, ALLOCATE MEMORY FOR THE RELEVANT ARRAYS. 
!      MEMORYRLS: MEMORY RELEASE. 
!      GETBKGRND: INTERPOLATE THE BACKGROUND FIELD TO ANALYSIS GRID POINTS. 
!      ADDBKGRND: SET BACKGROUND AS OBSERVATIONS AT THE GRID POINTS WHERE OBSERVATIONS CAN NOT AFFECT. 
!      GETBKGRND_NEW: MODIFIED FROM GETBKGRND.
!
!   ARRAYS:
!      GRIDMASK: MASK TO DEFINE WHETHERE THE GRID POINT IS AFFECTED BY OBSERVATIONS.  
!      OBRADIUS: THIS ARRAY DEFINES THE DISTANCE IN WHICH THE FIELD CAN BE AFFECT BY THE OBSERVATIONS.
!**************************************************

CONTAINS

!doc==================================================================
!
!>
!! This is a function of input_bg_obs module and included in
!! input_bg_obs.f90 for reading in STMAS 3D analysis namelist.
!! All of the namelist variables are defined in prmtrs_stmas.f90.
!!
!! \author Yuanfu Xie
!!
!! \b History \n
!! Created: Dec. 2013
!!
!endoc================================================================

!doc==================================================================
!
!>
!! This is a function of input_bg_obs module and included in
!! input_bg_obs.f90 for reading in STMAS 3D analysis namelist.
!! All of the namelist variables are defined in prmtrs_stmas.f90.
!!
!! \author Yuanfu Xie
!!
!! \b History \n
!! Created: Dec. 2013
!!
!endoc================================================================

SUBROUTINE read_namelist

  IMPLICIT NONE

  NAMELIST /number_states/numstat

  NAMELIST /stmas3d/ifbkgnd,ifbound,ifpcdnt, &
                    penal0x,penal0y,penal0z,penal0t, &
                    pnlt0pu,pnlt0pv,fnstgrd, &
                    numdims,numgrid,maxgrid, &
                    u_cmpnnt,v_cmpnnt,w_cmpnnt,pressure,temprtur,humidity, &
                    raincont,snowcont,grapcont,cloudice,cloudwat, &
                    cosstep,midgrid,finstep,pnlt0hy,taul_hy, &
                    endhylv,endgslv

  CHARACTER(LEN=256) ::filename
  INTEGER :: n,nm,ns,ierr,limgrid(2)
  REAL :: ratio

  ! Get namelist file from LAPS static:
  CALL get_directory('static',filename,n)
  filename(n:n+11) = '/stmas3d.nl'

  ! Open file for read:
  OPEN(11,file=filename(1:n+11))
  READ(11,NML=number_states,IOSTAT=ierr)

  ! Allocate memory:
  ALLOCATE(sl0(numstat),penal0x(numstat),penal0y(numstat),penal0z(numstat), &
           penal0t(numstat),penal_x(numstat),penal_y(numstat),penal_z(numstat), &
           penal_t(numstat),obradius(maxdims,numstat), STAT=ierr)

  READ(11,NML=stmas3d,IOSTAT=ierr)
  CLOSE(11)

  ! Default scaling:
  sl0 = 1.0

  ! Radar data influence radius: Option later for readin from namelist
  DO ns=1,numstat
    obradius(1:2,ns) = 200000.0
    obradius(3,ns) = 50000.0
    obradius(4,ns) = 0.0
  ENDDO

  ! Initial grid positions:
  oripstn = 0

  ! Coordinate indices:
  xsl = 1
  ysl = 2
  psl = 3
  csl = 4
  dsl = 5

  ! Unit conversion to meters:
  xytrans = 1.0
  z_trans = 1.0

  ! For testing: if_test = 1; otherwise,
  if_test = 0

  ! Multigrid V cycle option: 1 full cycle (coarse to fine repeated); 
  !                           0 half cycle (one coarse to fine)
  ifrepet = 0
  itrepet = 0

  ! Threshold values QC observation in the vertical and time:
  limit_3 = 10
  limit_4 = 1

  ! Initializing STMAS grids:
  ngptobs = 2**numdims        ! Number of observation grid indices
  nallobs = 0                ! Total number of all observations

  ! Multigrid setup:
  IF (maxgrid(1) .EQ. 0) THEN

    ! Using a default multigrid setup based on the fcstgrd:
    fnstgrd = 4       ! Default levels of multigrid

    ! Get LAPS fcstgrd:
    CALL get_grid_dim_xy(fcstgrd(1),fcstgrd(2),ierr)
    CALL get_LAPS_dimensions(fcstgrd(3),ierr)
    fcstgrd(4) = 3    ! Hardcode for now

    ! For X and Y directions:
    ! Limit of 401 maxgrid in x and y directions:
    ratio = FLOAT(fcstgrd(1)-1)/FLOAT(fcstgrd(2)-1)
    limgrid(1) = MIN(fcstgrd(1),401)
    limgrid(2) = MIN(fcstgrd(2),401)
    IF (ratio .GE. 1.0) THEN
      limgrid(2) = INT((limgrid(1)-1)/ratio)+1
    ELSE
      limgrid(1) = INT((limgrid(2)-1)/ratio)+1
    ENDIF
    DO n=1,2
       numgrid(n) = INT((limgrid(n)-1)/2**(fnstgrd-1))
       maxgrid(n) = 2**(fnstgrd-1)*numgrid(n)+1
       numgrid(n) = numgrid(n)+1
    ENDDO

    ! For Z direction:
    nm = 0
    numgrid(3) = fcstgrd(3)-1
    DO n=1,fnstgrd
      IF (MOD(numgrid(3),2) .EQ. 0) THEN
        numgrid(3) = numgrid(3)/2
        nm = nm+1
      ELSEIF (nm .EQ. 0) THEN
        PRINT*,'Currently, the number of analysis vertical levels must be odd!'
        STOP
      ELSE
        EXIT      ! Use current numgrid(3) to start multigrid
      ENDIF
    ENDDO
    maxgrid(3) = numgrid(3)*2**(nm-1)+1
    numgrid(3) = numgrid(3)+1

    ! For T direction:
    numgrid(4) = 2
    maxgrid(4) = 3

  ELSE

    ! Using maxgrid to setup multigrid:
    DO n=1,numdims
      IF(maxgrid(n) .GT. 1 .AND. numgrid(n) .GT. 1) THEN
        IF(MOD(maxgrid(n)-1,numgrid(n)-1) .EQ. 0) THEN
          nm=(maxgrid(N)-1)/(numgrid(N)-1)
          ns=1
          DO WHILE(nm .GE. 2)
            IF(MOD(nm,2) .EQ. 0) THEN
              nm=nm/2
              ns=ns+1
            ELSE
              PRINT*, 'MAXGRID SHOULD BE (NUMGRID-1)*2**N+1'
              STOP
            ENDIF
          ENDDO
        ELSE
          PRINT*, 'MAXGRID SHOULD BE (NUMGRID-1)*2**N+1'
          STOP
        ENDIF
        IF(ns .GT. fnstgrd) THEN
          maxgrid(n)=(numgrid(n)-1)*2**(fnstgrd-1)+1
        ENDIF
      ENDIF
    ENDDO

  ENDIF ! End multigrid setup

  ! Maxgrid in time is the same as final analysis:
  fcstgrd(4) = maxgrid(4)

  inigrid = numgrid          ! Save initial start multigrid numbers

  ! Initial vertical temporarl gridspacing
  grdspac(3:4) = 0.0
  IF (maxgrid(3) .GT. 1) grdspac(3) = (maxgrid(3)-1)/FLOAT(numgrid(3)-1)
  
  ! IPW data namelist:
  CALL ipws%getNamelist

  PRINT*,''
  PRINT*,'STMAS namelist has been read with'
  WRITE(*,1) numgrid, maxgrid
1 FORMAT(' numgrid: ',4i4,' maxgrid: ',4i4)
  PRINT*,''

END SUBROUTINE read_namelist

! INCLUDE 'laps_configs.f90'
SUBROUTINE LAPS_CONFIG

!********************************************************************
!  READ IN LAPS CONFIGURATION PARAMETERS.
!
!  HISTORY: JAN. 2009 BY YUANFU XIE.
!
!           MODIFIED DEC. 2013 BY YUANFU XIE FOR P_SFC_F USED FOR TPW
!           CALCULATION
!********************************************************************

  USE PRMTRS_STMAS

  IMPLICIT NONE

  INTEGER :: N,ISTATE

  ! SYSTEM TIME:
  CALL GET_SYSTIME(LAPSI4T,LAPSAST,ISTATE)

  ! LAPS CYCLE TIME:
  CALL GET_LAPS_CYCLE_TIME(ICYCLE,ISTATE)

  ! I4TIME OF JAN 1, 1970
  CALL CV_ASC_I4TIME('700010000',I4T_GPS)

  ! MISSING VALUE IN REAL:
  CALL GET_R_MISSING_DATA(RMISSING,ISTATE)
  ! BAD SURFACE DATA FLAG:
  CALL GET_SFC_BADFLAG(BADSFCDT,ISTATE)

  ! SPATIAL X Y DIMENSIONS:
  CALL GET_GRID_DIM_XY(FCSTGRD(1),FCSTGRD(2),ISTATE)
  CALL GET_LAPS_DIMENSIONS(FCSTGRD(3),ISTATE)
  PRINT*,(FCSTGRD(N),N=1,NUMDIMS)

  ! LAPS ARRAY:
  GRDSPAC(1)=((FCSTGRD(1)-1)*1.0)/((NUMGRID(1)-1)*1.0)
  GRDSPAC(2)=((FCSTGRD(2)-1)*1.0)/((NUMGRID(2)-1)*1.0)
  IF(NUMGRID(4).GE.2) THEN
    GRDSPAC(4)=(FCSTGRD(4)-1.0)/(NUMGRID(4)-1.0)
  ELSE
    GRDSPAC(4)=0.0
  ENDIF
  ALLOCATE(GRIDMASK(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4),NUMSTAT),STAT=ISTATE)
  IF(ISTATE.NE.0)STOP 'GRIDMASK ALLOCATE WRONG'

  ALLOCATE(Z_FCSTGD(FCSTGRD(3)),STAT=ISTATE)
  IF(ISTATE.NE.0)STOP 'Z_FCSTGD ALLOCATE WRONG'
  ALLOCATE(Z_MAXGID(MAXGRID(3)),STAT=ISTATE)
  IF(ISTATE.NE.0)STOP 'Z_MAXGID ALLOCATE WRONG'

  CALL ALLOCTBKG
  CALL ALLOCATOB

  ! LAPS LAT/LON/TOPOGRAPHY:
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'LAT',LATITUDE,ISTATE)
  IF (ISTATE .NE. 1) THEN
    WRITE(6,*) 'RDLAPSRDR: error get LAPS LAT'
    STOP
  ENDIF
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'LON',LONGITUD,ISTATE)
  IF (ISTATE .NE. 1) THEN
    WRITE(6,*) 'RDLAPSRDR: error get LAPS LON'
    STOP
  ENDIF
  CALL READ_STATIC_GRID(FCSTGRD(1),FCSTGRD(2),'AVG',TOPOGRPH,ISTATE)
  IF (ISTATE .NE. 1) THEN
    WRITE(6,*) 'RDLAPSRDR: error get LAPS AVG'
    STOP
  ENDIF

  ! ALLOCATE MEMORY FOR SURFACE PRESSURE:
  ALLOCATE (P_SFC_F(FCSTGRD(1),FCSTGRD(2)),stat=istate)

END SUBROUTINE LAPS_CONFIG


SUBROUTINE BKGRNDOBS
!*************************************************
! MAIN ROUTINE OF THIS PREPARATION CODE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE

  INTEGER      :: ISTATE,i
  
  INTEGER      :: nfiles
  CLASS(*), POINTER :: t
  TYPE(list), POINTER :: c

  PRINT*,'READNMLST'
  ! CALL READNMLST
  CALL read_namelist

  PRINT*,'BKGMEMALC'
  CALL BKGMEMALC
  PRINT*,'RDBCKGRND'
  IF(IF_TEST.NE.1) THEN
    ! LAPS CONFIGURATION:
    CALL LAPS_CONFIG
    !  CALL RDBCKGRND
    CALL RDLAPSBKG
  ELSE
    OPEN(11,FILE='fort.11',STATUS='OLD',ACTION='READ')
      READ(11,*)FCSTGRD(1:NUMDIMS)
    CLOSE(11)
    PRINT*,'MEMORYALC'
    CALL MEMORYALC
    CALL RDBKGTEST
  ENDIF
  PRINT*,'GETBKGRND_NEW'
  CALL GETBKGRND_NEW
  PRINT*,'RDBUFROBS'
  IF(IF_TEST.NE.1) THEN
    ! CALL OPEN_LAPSPRD_FILE(TMGOBS_CHANNEL,LAPSI4T,'tmg',ISTATE)
    ! CALL OPEN_LAPSPRD_FILE(PIGOBS_CHANNEL,LAPSI4T,'pig',ISTATE)
    CALL RDBUFROBS_XIE
    ! CLOSE(TMGOBS_CHANNEL)
    ! CLOSE(PIGOBS_CHANNEL)
  !  CALL RDRADROBS
    ! CALL RDLAPSRDR
    call read_laps_radar  ! Switch to a new routine using less memory by Yuanfu
    ! CALL GPSWDELAY
    
    ! IPW data: by Yuanfu Xie 2014: 
    ! numins: number of instruments from lapsdata module:
    IF (numins .GT. 0) THEN
      ALLOCATE(ipwdata,STAT=istate)
      CALL ipwdata%initial

      DO i=1,numins
print*,'Ins: ',names(i),trim(paths(i)),' ',trim(patterns(i)),' ',trim(xpatterns(i))
        ! GPS is located at index 1: only gps data ingested
        SELECT CASE (TRIM(names(i)))
        CASE ('gps')
          ALLOCATE(gpsTPW_t :: idat)
        CASE ('snd')
          ALLOCATE(sndIPW_t :: idat)
        CASE DEFAULT
          PRINT*,'Undefined IPW data'
          STOP
        END SELECT
           
        idat%a9time = ipws%a9time
        idat%timeWindow = ipws%timeWindow
        idat%numFiles = ipws%numFiles
        idat%filenames = ipws%filenames
        CALL idat%getFilePath(i)
    
        CALL idat%readData(i,ipwdata,p_sfc_f,latitude,longitud,fcstgrd)
        
        DEALLOCATE(idat)
      ENDDO
      
    ENDIF
    
  ELSE 
    CALL RDOBSTEST
  ENDIF
  IF(IFBKGND.EQ.1)THEN
    PRINT*,'ADDBKGRND'
    CALL ADDBKGRND
  ENDIF
  PRINT*,'OBSMEMALC'
  CALL OBSMEMALC
  PRINT*,'MEMORYRLS'
  CALL MEMORYRLS
  RETURN
END SUBROUTINE BKGRNDOBS

SUBROUTINE READNMLST
!*************************************************
! READ NAME LIST
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: S,N,ER,NU,LN,NM,NS
  CHARACTER(LEN=200) :: DR
  CHARACTER(LEN=256) :: ST

  ! GET NAMELIST FROM STATIC:
  CALL GET_DIRECTORY('static',ST,N)
! replacing the file name from stmas3d.nl to stmas3d.txt, and also increase the number by 1. N+10 to N+11. HJ 6/28/2011 
!  DR = ST(1:N)//'stmas3d.nl'
  DR = ST(1:N)//'stmas3d.txt'
! --------------------
  NU=2
  !OPEN(NU,FILE='namelist.txt',ACTION='READ',STATUS='OLD')
! replacing the file name from stmas3d.nl to stmas3d.txt, and also increase the number by 1. N+10 to N+11. HJ 6/28/2011 
!  OPEN(NU,FILE=DR(1:N+10),ACTION='READ',STATUS='OLD')
  OPEN(NU,FILE=DR(1:N+11),ACTION='READ',STATUS='OLD')


  READ(NU,*)IFBKGND
  READ(NU,*)IFBOUND
  READ(NU,*)IFPCDNT

  IF(IFPCDNT.EQ.1) THEN
    PRINT*,'IFPCDNT=',IFPCDNT,'FOR PRESURE COORDINATE'
  ELSEIF(IFPCDNT.EQ.2) THEN
    PRINT*,'IFPCDNT=',IFPCDNT,'FOR HEIGHT COORDINATE'
  ELSEIF(IFPCDNT.EQ.0) THEN
    PRINT*,'IFPCDNT=',IFPCDNT,'FOR SIGMA COORDINATE'
  ENDIF

  READ(NU,*)NUMSTAT
  ALLOCATE(SL0(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'SL0 ALLOCATE WRONG'
  DO S=1,NUMSTAT
    READ(NU,*)SL0(S)
  ENDDO
  ALLOCATE(PENAL0X(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0X ALLOCATE WRONG'
  ALLOCATE(PENAL0Y(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0Y ALLOCATE WRONG'
  ALLOCATE(PENAL0Z(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0Z ALLOCATE WRONG'
  ALLOCATE(PENAL0T(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL0T ALLOCATE WRONG'
  ALLOCATE(PENAL_X(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_X ALLOCATE WRONG'
  ALLOCATE(PENAL_Y(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_Y ALLOCATE WRONG'
  ALLOCATE(PENAL_Z(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_Z ALLOCATE WRONG'
  ALLOCATE(PENAL_T(NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'PENAL_T ALLOCATE WRONG'
  DO S=1,NUMSTAT
    READ(NU,*)PENAL0X(S)
    READ(NU,*)PENAL0Y(S)
    READ(NU,*)PENAL0Z(S)
    READ(NU,*)PENAL0T(S)
    PRINT*,'SMOOTHING VAR: ',S,' WITH: ',PENAL0X(S),PENAL0Y(S),PENAL0Z(S),PENAL0T(S)
  ENDDO
  READ(NU,*)PNLT0PU
  READ(NU,*)PNLT0PV
  READ(NU,*)NUMDIMS
  READ(NU,*)NUMGRID(1)
  READ(NU,*)NUMGRID(2)
  READ(NU,*)NUMGRID(3)
  READ(NU,*)NUMGRID(4)
  READ(NU,*)FCSTGRD(4)		! FOR LAPS INGEST: READ INTO LAPS TIME FRAME FOR TEMPORAL ANALYSIS YUANFU
  READ(NU,*)GRDSPAC(4)
  READ(NU,*)ORIPSTN(1)
  READ(NU,*)ORIPSTN(2)
  READ(NU,*)ORIPSTN(3)
  READ(NU,*)ORIPSTN(4)
  READ(NU,*)MAXGRID(1)
  READ(NU,*)MAXGRID(2)
  READ(NU,*)MAXGRID(3)
  READ(NU,*)MAXGRID(4)
  READ(NU,*)FNSTGRD
  ALLOCATE(OBRADIUS(MAXDIMS,NUMSTAT),STAT=ER)
  IF(ER.NE.0)STOP 'OBRADIUS ALLOCATE WRONG'
  DO S=1,NUMSTAT
    DO N=1,MAXDIMS
      READ(NU,*)OBRADIUS(N,S)
    ENDDO
  ENDDO
  READ(NU,*)U_CMPNNT
  READ(NU,*)V_CMPNNT
  READ(NU,*)W_CMPNNT
  READ(NU,*)PRESSURE         ! 'PRESSURE' IS FOR Z COORDINATE, FOR PRESURE COORDINATE IT IS HEIGHT
  READ(NU,*)TEMPRTUR
  READ(NU,*)HUMIDITY
  IF (NUMSTAT .GT. 5) &
    READ(NU,*) ROUR_CMPNNT  ! added by shuyuan20100721 for density(rou)*rain water mixing ratio
  IF (NUMSTAT .GT. 6) &
    READ(NU,*) ROUS_CMPNNT  ! added by shuyuan20100721 for density(rou)*snow water mixing ratio
  READ(NU,*)XSL
  READ(NU,*)YSL
  READ(NU,*)PSL
  READ(NU,*)CSL
  READ(NU,*)DSL
  READ(NU,*)COSSTEP
  READ(NU,*)MIDGRID
  READ(NU,*)FINSTEP

  READ(NU,*)XYTRANS             ! COEFFICIENT USED TO TRANSLATE THE X AND Y COORDINATE TO METERS
  READ(NU,*)Z_TRANS             ! COEFFICIENT USED TO TRANSLATE THE Z COORDINATE TO METERS
  READ(NU,*)IF_TEST             ! WHETHER RUN FOR TEST CASE, 1 IS FOR THE TEST CASE
  READ(NU,*)IFREPET             ! WHETHERE RUN THE MUTIGRID FRAME IN MONOTONOUS OR REPEATEDLY, 0 FOR MONOTONOUS, 1 FOR REPEATEDLY
  READ(NU,*)ITREPET             ! THE TIMES TO REPEAT
  IF(IFREPET.EQ.0) ITREPET=0

  READ(NU,*)PNLT0HY             ! INITIAL HYDROSTATIC CONDITION PENALTY COEFFICENT
  READ(NU,*)TAUL_HY             ! REDUCTION COEFFICENT OF HYDROSTATIC CONDITION PENALTY TERM
  READ(NU,*)ENDHYLV             ! THE GRID LEVEL AFTER WHICH THE HYDROSTATIC CONDITION PENALTY TERM IS OMMITTED
  READ(NU,*)ENDGSLV             ! THE GRID LEVEL AFTER WHICH THE GEOSTROPHIC BALANCE PENALTY TERM IS OMMITTED

print*,'HHH: ',humidity,midgrid,endgslv,endhylv
  READ(NU,*)LIMIT_3             ! THE LIMITATION OF HEIGHT OR PRESSURE TO DECIDE IN WHICH RANGE THE OBSERVATION IS AVIABLE. 
  READ(NU,*)LIMIT_4             ! THE LIMITATION OF TIME TO DECIDE IN WHICH RANGE THE OBSERVATION IS AVIABLE.

  CLOSE(NU)
  NGPTOBS=2**NUMDIMS
  NALLOBS=0                     ! INITIALIZING THE NUMBER OF OBSERVATION, ADDIED BY ZHONGJIE HE 
  GRDSPAC(3)=0.
  IF(MAXGRID(3).NE.1) GRDSPAC(3)=(MAXGRID(3)-1.0)/(NUMGRID(3)-1.0)      ! BY ZHONGJIE HE

  DO N=1,NUMDIMS
    IF(MAXGRID(N).GT.1 .AND. NUMGRID(N).GT.1) THEN
      IF(MOD(MAXGRID(N)-1,NUMGRID(N)-1).EQ.0) THEN
        NM=(MAXGRID(N)-1)/(NUMGRID(N)-1)
        NS=1
        DO WHILE(NM.GE.2)
          IF(MOD(NM,2).EQ.0) THEN
            NM=NM/2
            NS=NS+1
          ELSE
            PRINT*, 'MAXGRID SHOULD BE (NUMGRID-1)*2**N+1'
            STOP
          ENDIF
        ENDDO
      ELSE
        PRINT*, 'MAXGRID SHOULD BE (NUMGRID-1)*2**N+1'
        STOP
      ENDIF
      IF(NS.GT.FNSTGRD) THEN
        MAXGRID(N)=(NUMGRID(N)-1)*2**(FNSTGRD-1)+1
      ENDIF
    ENDIF
  ENDDO

  DO N=1,NUMDIMS
    INIGRID(N)=NUMGRID(N)
  ENDDO

  RETURN
END SUBROUTINE READNMLST

SUBROUTINE BKGMEMALC
!*************************************************
! ALLOCATE MEMORY FOR BACKGROUND ARRAY
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: ER
! --------------------
  ! TWO ADDITIONAL SPACES OF GRDBKGD0 FOR SAVING LOWER AND UPPER BOUNDS OF SH:
  ALLOCATE(GRDBKGD0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4),NUMSTAT+2),STAT=ER)
  IF(ER.NE.0)STOP 'GRDBKGD0 ALLOCATE WRONG'
  ALLOCATE(XX0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'XX0 ALLOCATE WRONG'
  ALLOCATE(YY0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'YY0 ALLOCATE WRONG'
  ALLOCATE(CR0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'CR0 ALLOCATE WRONG'
  ALLOCATE(DG0(MAXGRID(1),MAXGRID(2)),STAT=ER)
  IF(ER.NE.0)STOP 'DG0 ALLOCATE WRONG'
  ALLOCATE(DN0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4)),STAT=ER)
  IF(ER.NE.0)STOP 'DN0 ALLOCATE WRONG'
  IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN 
    ALLOCATE(ZZ0(MAXGRID(1),MAXGRID(2),MAXGRID(3),MAXGRID(4)),STAT=ER)
    IF(ER.NE.0)STOP 'ZZ0 ALLOCATE WRONG'
    ALLOCATE(ZZB(MAXGRID(3)),STAT=ER)
    IF(ER.NE.0)STOP 'ZZB ALLOCATE WRONG'
  ELSEIF(IFPCDNT.EQ.1)THEN
    ALLOCATE(PP0(MAXGRID(3)),PPM(MAXGRID(3)),STAT=ER)
    IF(ER.NE.0)STOP 'PP0 AND PPM ALLOCATE WRONG'
  ENDIF
  RETURN
END SUBROUTINE BKGMEMALC

SUBROUTINE OBSMEMALC
!*************************************************
! ALLOCATE OBSERVATION MEMORY AND READ IN DATA AND SCALE
! HISTORY: AUGUST 2007, CODED by WEI LI.
!        : MARCH 2008, MODIFIED BY ZHONGJIE HE
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: N,O,ER,LN,S,NO
  CHARACTER(LEN=200) :: DR
! --------------------
  IF(NALLOBS.EQ.0)RETURN
  ALLOCATE(OBSPOSTN(NUMDIMS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSPOSTN ALLOCATE WRONG'
  ALLOCATE(OBSCOEFF(NGPTOBS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSCOEFF ALLOCATE WRONG'
  ALLOCATE(OBSIDXPC(MAXDIMS,NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSIDXPC ALLOCATE WRONG'
  ALLOCATE(OBSVALUE(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSVALUE ALLOCATE WRONG'
  ALLOCATE(OBSERROR(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSERROR ALLOCATE WRONG'
!jhui
  ALLOCATE(NOBSTAT(NUMSTAT+3),STAT=ER)
  IF(ER.NE.0)STOP 'NOBSTAT ALLOCATE WRONG'
  ALLOCATE(OBSEINF1(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSEINF1 ALLOCATE WRONG'
  ALLOCATE(OBSEINF2(NALLOBS),STAT=ER)
  IF(ER.NE.0)STOP 'OBSEINF2 ALLOCATE WRONG'
!  ALLOCATE(OBSEINF3(NALLOBS),STAT=ER)
!  IF(ER.NE.0)STOP 'OBSEINF3 ALLOCATE WRONG'

  O=0
  DO S=1,NUMSTAT+3
    DO NO=1,NST(S)
      O=O+1
      DO N=1,NUMDIMS
        OBSPOSTN(N,O)=OBP(N,NO,S)
      ENDDO
      OBSVALUE(O)=OBS(NO,S)
      OBSERROR(O)=OBE(NO,S)
    ENDDO
    NOBSTAT(S)=NST(S)
  ENDDO
  S=NUMSTAT+1
  DO NO=1,NST(S)
    OBSEINF1(NO)=OBA(NO,1)
    OBSEINF2(NO)=OBA(NO,2)
  ENDDO
!jhui
!  S=NUMSTAT+3
!  DO NO=1,NST(S)
!   WRITE(306,*) OBS(NO,S)
!  ENDDO

  RETURN
END SUBROUTINE OBSMEMALC

SUBROUTINE MEMORYALC
!*************************************************
! MEMORY ALLOCATE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: N,ER
  INTEGER :: ST		! STATUS
! --------------------

  GRDSPAC(1)=((FCSTGRD(1)-1)*1.0)/((NUMGRID(1)-1)*1.0)
  GRDSPAC(2)=((FCSTGRD(2)-1)*1.0)/((NUMGRID(2)-1)*1.0)
  IF(NUMGRID(4).GE.2) THEN
    GRDSPAC(4)=(FCSTGRD(4)-1.0)/(NUMGRID(4)-1.0)
  ELSE
    GRDSPAC(4)=0.0
  ENDIF
  ALLOCATE(GRIDMASK(FCSTGRD(1),FCSTGRD(2),MAXGRID(3),FCSTGRD(4),NUMSTAT),STAT=ER) !!!!!ATTENTION
  IF(ER.NE.0)STOP 'GRIDMASK ALLOCATE WRONG'

  ALLOCATE(Z_FCSTGD(FCSTGRD(3)),STAT=ER)
  IF(ER.NE.0)STOP 'Z_FCSTGD ALLOCATE WRONG'
  ALLOCATE(Z_MAXGID(MAXGRID(3)),STAT=ER)
  IF(ER.NE.0)STOP 'Z_MAXGID ALLOCATE WRONG'

  CALL ALLOCTBKG
  CALL ALLOCATOB
  RETURN
END SUBROUTINE MEMORYALC

SUBROUTINE MEMORYRLS
!*************************************************
! MEMORY RELEASE
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
  DEALLOCATE(GRIDMASK)
  DEALLOCATE(OBRADIUS)
  CALL DEALCTBKG
  CALL DEALOCTOB
  RETURN
END SUBROUTINE MEMORYRLS

SUBROUTINE GETBKGRND_NEW
!*************************************************
! GET DATA OF BACKGROUND FOR ANALYSIS, ADAPT TO PRESSURE HEIGHT AND SIGMA COORDINATE.
! HISTORY: SEPTEMBER 2007, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,LN,NM
  REAL     :: XB(MAXGRID(1)),YB(MAXGRID(2)),ZB(MAXGRID(3)),TB(MAXGRID(4))
  REAL     :: XF(FCSTGRD(1)),YF(FCSTGRD(2)),ZF(FCSTGRD(3)),TF(FCSTGRD(4))
  INTEGER  :: FG(MAXDIMS),MG(MAXDIMS)
  REAL     :: Z1(1),T1(1),Z2(1),T2(1),DX,DY,DT
  CHARACTER(LEN=200) :: DR

  REAL     :: SH          ! COEFFICIENT TO TRANSLATE THE ZZ0 TO LENGTH WITH UNIT OF METERS.
  LOGICAL  :: LEXIST      ! USED TO DECIDE WHETHER THE LEVEL DATA FILE IS EXIST
  INTEGER :: CLOCK_COUNT0,CLOCK_RATE,CLOCK_MAX,CLOCK_COUNT1

  integer :: init_timer,ishow_timer
! --------------------

  SH=Z_TRANS

  DO I=1,FCSTGRD(1)
    XF(I)=(I-1)*1.0D0
  ENDDO
  DO J=1,FCSTGRD(2)
    YF(J)=(J-1)*1.0D0
  ENDDO
  DX=((FCSTGRD(1)-1)*1.0D0)/((MAXGRID(1)-1)*1.0D0)
  DO I=1,MAXGRID(1)
    XB(I)=XF(1)+(I-1)*DX
  ENDDO
  DY=((FCSTGRD(2)-1)*1.0D0)/((MAXGRID(2)-1)*1.0D0)
  DO J=1,MAXGRID(2)
    YB(J)=YF(1)+(J-1)*DY
  ENDDO

  IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
    DO K=1,FCSTGRD(3)
      ZF(K)=Z00(K)
    ENDDO
  ELSE                                         ! FOR PRESSURE COORDINATE
    DO K=1,FCSTGRD(3)
      ZF(K)=P00(K)
    ENDDO
  ENDIF

!  OPEN(2,FILE='P_LEVEL.DAT',STATUS='OLD',ACTION='READ')
!    DO K=1,MAXGRID(3)
!      READ(2,*)ZB(K)
!    ENDDO
!  CLOSE(2)
  IF(MAXGRID(3).LT.FCSTGRD(3)) THEN
    PRINT*, 'READ THE VERTICAL LEVELS FROM DATA FILE: LEVEL.DAT'
    INQUIRE(FILE='LEVEL.DAT',EXIST=LEXIST)
    IF(LEXIST) THEN
      OPEN(2,FILE='LEVEL.DAT',STATUS='OLD',ACTION='READ')
        DO K=1,MAXGRID(3)
          READ(2,*)ZB(K)
        ENDDO
      CLOSE(2)
    ELSE
      PRINT*, 'THE FILE LEVEL.DAT DOES NOT EXIST'
      STOP
    ENDIF
  ELSEIF(MAXGRID(3).GT.2*FCSTGRD(3)-1) THEN
    PRINT*, 'ERROR! MAXGRID(3) SHOULD BE SMALLER THAN 2*FCSTGRD(3)-1!'
    STOP
  ELSEIF(MAXGRID(3).EQ.FCSTGRD(3)) THEN
    DO K=1,MAXGRID(3)
      ZB(K)=ZF(K)
    ENDDO
  ELSE
    NM=MAXGRID(3)-FCSTGRD(3)
    DO K=1,NM
      ZB(2*K-1)=ZF(K)
      ZB(2*K)=0.5*(ZF(K)+ZF(K+1))
    ENDDO
    DO K=NM+1,FCSTGRD(3)
      ZB(K+NM)=ZF(K)
    ENDDO
!======
!    DO K=1,FCSTGRD(3)-1
!      ZB(K)=ZF(K)
!    ENDDO
!    DO K=FCSTGRD(3),MAXGRID(3)
!      ZB(K)=ZF(FCSTGRD(3)-1)+(ZF(FCSTGRD(3))-ZF(FCSTGRD(3)-1))/FLOAT(MAXGRID(3)-FCSTGRD(3)+1)*(K-FCSTGRD(3)+1)
!    ENDDO
!======
  ENDIF

!  DO K=1,MAXGRID(3)
!    ZB(K)=ZF(1)+(ZF(FCSTGRD(3))-ZF(1))/FLOAT(MAXGRID(3)-1)*(K-1)
!  ENDDO

!===============
  DO K=1,MAXGRID(3)
    Z_MAXGID(K)=ZB(K)
  ENDDO
!===============

  DO K=1,MAXGRID(3)
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2) THEN       ! FOR SIGMA AND HEIGHT COORDINATE
      ZZB(K)=ZB(K)
    ELSE                                          ! FOR PRESSURE COORDINATE
      PP0(K)=ZB(K)
    ENDIF
  ENDDO
  ! MULTIGRID PRESSURE COORDINATE:
  I = (MAXGRID(3)-1)/(NUMGRID(3)-1)
  DO K=1,NUMGRID(3)
    PPM(K) = PP0((K-1)*I+1)
  ENDDO

  IF(IFPCDNT.EQ.0) THEN                ! FOR SIGMA COORDINATE
    DO T=1,MAXGRID(4)
    DO J=1,MAXGRID(2)
    DO I=1,MAXGRID(1)
    DO K=1,MAXGRID(3)
      ZZ0(I,J,K,T)=(ZB(K)*(HEIGHTU(I,J)-HEIGHTL(I,J))+HEIGHTL(I,J))*SH   ! TRANSLATE THE UNIT TO METERS
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ELSEIF(IFPCDNT.EQ.2) THEN            ! FOR HEIGHT COORDINATE
    DO T=1,MAXGRID(4)
    DO J=1,MAXGRID(2)
    DO I=1,MAXGRID(1)
    DO K=1,MAXGRID(3)
      ZZ0(I,J,K,T)=ZB(K)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDIF

!=======================================

  DO T=1,FCSTGRD(4)
    TF(T)=(T-1)*1.0
  ENDDO
  IF(MAXGRID(4).GE.2) DT=((FCSTGRD(4)-1)*1.0)/((MAXGRID(4)-1)*1.0)
  DO T=1,MAXGRID(4)
    TB(T)=TF(1)+(T-1)*DT
  ENDDO

!  CALL SYSTEM_CLOCK(CLOCK_COUNT0,CLOCK_RATE,CLOCK_MAX)
!  CALL FCST2BKGD(NUMDIMS,NGPTOBS,NUMSTAT,FCSTGRD, &
!                 XF,YF,ZF,TF,MAXGRID,XB,YB,ZB,TB,BK0,GRDBKGD0)

  PRINT*,'Starting interpolating background onto maxgrid: ', init_timer()
  IF (UNIFORM .EQ. 0) THEN
    CALL BKGTOFINE(NUMSTAT,FCSTGRD,XF,YF,ZF,TF,MAXGRID,XB,YB,ZB,TB,BK0,GRDBKGD0)
  ELSE ! USE YUANFU'S UNIFORM INTERPOLATION ROUTINE:
    DO S=1,NUMSTAT
    DO T=1,FCSTGRD(4)
      CALL uniform_interpolation3(fcstgrd,maxgrid,BK0(:,:,:,T,S),GRDBKGD0(:,:,:,T,S))
    ENDDO
    ENDDO
  ENDIF
  PRINT*,'Ending : interpolating background onto maxgrid: ', ishow_timer()

!  CALL SYSTEM_CLOCK(CLOCK_COUNT1,CLOCK_RATE,CLOCK_MAX)
!  write ( *, '(a)' ) '  SYSTEM_CLOCK with INTEGER arguments reports: in subroutine of <FCST2BKGD> '
!  write ( *, '(a,i12)' ) '    The current clock count is    ', clock_count1-clock_count0
!  write ( *, '(a,i12)' ) '    The clock count per second is ', clock_rate
!  write ( *, '(a,i12)' ) '    The maximum clock count is    ', clock_max

  ! YUANFU REMOVE S00 ARRAY AS IT WAS SET TO 1 AND IT WAS PASSED TO DN0
  ! CALL BKGTOFINE(1,FCSTGRD,XF,YF,ZF,TF,MAXGRID,XB,YB,ZB,TB,S00,DN0)
  DN0 = 1.0 ! DENSITY ON THE MULTIGRID. IT IS SUPPOSED TO USE REAL DENSITY BUT CONSTANT FOR NOW.

  FG(1)=FCSTGRD(1)
  FG(2)=FCSTGRD(2)
  FG(3)=1
  FG(4)=1
  Z1(1)=0.0
  T1(1)=0.0
  MG(1)=MAXGRID(1)
  MG(2)=MAXGRID(2)
  MG(3)=1
  MG(4)=1
  Z2(1)=0.0
  T2(1)=0.0

!  CALL FCST2BKGD(2,4,1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,C00,CR0)
  IF (UNIFORM .EQ. 0) THEN
    CALL BKGTOFINE(1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,C00,CR0)
  ELSE
    CALL uniform_interpolation2(FG,MG,C00,CR0)
  ENDIF

!  CALL FCST2BKGD(2,4,1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,D00,DG0)
  ! YUANFU SKIPS THIS CALL FOR DG0 SINCE D00 IS SET TO ZERO, A CONSTANT IN READ_BACKGRD:
  ! CALL BKGTOFINE(1,FG,XF,YF,Z1,T1,MG,XB,YB,Z2,T2,D00,DG0)
  DG0 = 0.0

  DX=((FCSTGRD(1)-1)*DX0)/((MAXGRID(1)-1)*1.0D0)
  DY=((FCSTGRD(2)-1)*DY0)/((MAXGRID(2)-1)*1.0D0)
  IF(MAXGRID(4).GE.2) DT=((FCSTGRD(4)-1)*DT0)/((MAXGRID(4)-1)*1.0)
  DO I=1,MAXGRID(1)
  DO J=1,MAXGRID(2)
    XX0(I,J)=0.0D0+(I-1)*DX
    YY0(I,J)=0.0D0+(J-1)*DY
  ENDDO
  ENDDO

  RETURN
END SUBROUTINE GETBKGRND_NEW

SUBROUTINE ADDBKGRND
!*************************************************
! SET BACKGROUND AS OBSERVATIONS AT THE GRID POINTS WHERE OBSERVATIONS CAN NOT AFFECT
! HISTORY: OCTOBER 2007, CODED by WEI LI.
!          MARCH 2008, MODIFIED BY ZHONGJIE HE
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,N,RX,RY,RT,NP(MAXDIMS)
  REAL     :: RH,RZ,OE(NUMSTAT),OC(NUMDIMS),P,PP(MAXGRID(3))
  INTEGER  :: SS,OS                                             ! ADDED BY ZHONGJIE HE
! --------------------
  IF(NALLOBS.EQ.0)RETURN
  DO S=1,NUMSTAT
    DO T=1,FCSTGRD(4)
    DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
    DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      GRIDMASK(I,J,K,T,S)=1
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  DO S=1,NUMSTAT
    OE(S)=0.0D0
  ENDDO
  DO S=1,NUMSTAT
    DO O=1,NST(S)
      DO N=1,NUMDIMS
        OC(N)=OBP(N,O,S)+1.0
      ENDDO
      IF(IFPCDNT.EQ.1)THEN       ! FOR PRESSURE COORDINATE
        DO K=1,MAXGRID(3)
          PP(K)=PP0(K)
        ENDDO
        IF(MAXGRID(3).GE.2)THEN
          K=INT(OC(3))
          P=(OC(3)-K)*(PP(K+1)-PP(K))+PP(K)
        ENDIF
      ELSE                      ! FOR SIGMA AND HEIGHT COORDINATE 
        DO K=1,MAXGRID(3)
          PP(K)=Z00(K)
        ENDDO
        IF(MAXGRID(3).GE.2)THEN
          K=INT(OC(3))
          P=(OC(3)-K)*(PP(K+1)-PP(K))+PP(K)
        ENDIF
      ENDIF
      RX=0
      RY=0
      RT=0
      IF(FCSTGRD(1).GE.2)RX=INT(OBRADIUS(1,S)/DX0)+1
      IF(FCSTGRD(2).GE.2)RY=INT(OBRADIUS(2,S)/DY0)+1
      IF(FCSTGRD(4).GE.2)RT=INT(OBRADIUS(4,S)/DT0)+1
      DO T=MAX0(INT(OC(4))-RT,1),MIN0(INT(OC(4))+RT+1,FCSTGRD(4))
      DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
      DO J=MAX0(INT(OC(2))-RY,1),MIN0(INT(OC(2))+RY+1,FCSTGRD(2))
      DO I=MAX0(INT(OC(1))-RX,1),MIN0(INT(OC(1))+RX+1,FCSTGRD(1))
        RH=SQRT(((OC(1)-I)*DX0)**2+((OC(2)-J)*DY0)**2)
        IF(IFPCDNT.EQ.1 .OR. IFPCDNT.EQ.2)THEN       ! FOR PRESSURE AND HEIGHT COORDINATE
          RZ=ABS((PP(K)-P))
        ELSEIF(IFPCDNT.EQ.0)THEN                     ! FOR SIGMA COORDINATE
          RZ=ABS((PP(K)-P))*(HEIGHTU(I,J)-HEIGHTL(I,J))
        ENDIF
        IF(RH.LE.OBRADIUS(1,S).AND.RZ.LE.OBRADIUS(3,S))GRIDMASK(I,J,K,T,S)=0
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      OE(S)=MAX(OE(S),OBE(O,S)*1.0)
    ENDDO
  ENDDO

  O=NALLOBS
  DO S=1,NUMSTAT
    DO T=1,FCSTGRD(4)
    DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
    DO J=1,FCSTGRD(2),3
    DO I=1,FCSTGRD(1),3
      IF(GRIDMASK(I,J,K,T,S).EQ.1.AND.NST(S).GE.1)THEN
        NP(1)=I-1
        NP(2)=J-1
        NP(3)=K-1
        NP(4)=T-1
        O=O+1
        NST(S)=NST(S)+1
        IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
        DO N=1,NUMDIMS
          OBP(N,NST(S),S)=NP(N)
        ENDDO
        OBS(NST(S),S)=0.0D0
        OBE(NST(S),S)=OE(S)
      ENDIF
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  NALLOBS=O

  IF(.TRUE.) THEN                      ! THE FOLLOWING IS ADDED BY ZHONGJIE HE
    DO OS=NUMSTAT+1,NUMSTAT+2
    DO O=1,NST(OS)
      DO N=1,NUMDIMS
        OC(N)=OBP(N,O,OS)+1.0
      ENDDO
      IF(IFPCDNT.EQ.1)THEN        ! FOR PRESSURE COORDINATE
        DO K=1,MAXGRID(3)
          PP(K)=PP0(K)
        ENDDO
        IF(MAXGRID(3).GE.2)THEN
          K=INT(OC(3))
          P=(OC(3)-K)*(PP(K+1)-PP(K))+PP(K)
        ENDIF
      ELSE                       ! FOR SIGMA AND HEIGHT COORDINATE
        DO K=1,MAXGRID(3)
          PP(K)=Z00(K)
        ENDDO
        IF(MAXGRID(3).GE.2)THEN
          K=INT(OC(3))
          P=(OC(3)-K)*(PP(K+1)-PP(K))+PP(K)
        ENDIF
      ENDIF
      RX=0
      RY=0
      RT=0
      IF(FCSTGRD(1).GE.2)RX=INT(MAX(OBRADIUS(1,U_CMPNNT),OBRADIUS(1,V_CMPNNT))/DX0)+1
      IF(FCSTGRD(2).GE.2)RY=INT(MAX(OBRADIUS(2,U_CMPNNT),OBRADIUS(2,V_CMPNNT))/DY0)+1
      IF(FCSTGRD(4).GE.2)RT=INT(MAX(OBRADIUS(4,U_CMPNNT),OBRADIUS(4,V_CMPNNT))/DT0)+1
      DO T=MAX0(INT(OC(4))-RT,1),MIN0(INT(OC(4))+RT+1,FCSTGRD(4))
      DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
      DO J=MAX0(INT(OC(2))-RY,1),MIN0(INT(OC(2))+RY+1,FCSTGRD(2))
      DO I=MAX0(INT(OC(1))-RX,1),MIN0(INT(OC(1))+RX+1,FCSTGRD(1))
        RH=SQRT(((OC(1)-I)*DX0)**2+((OC(2)-J)*DY0)**2)
        IF(IFPCDNT.EQ.1 .OR. IFPCDNT.EQ.2)THEN       ! FOR PRESSURE AND HEIGHT COORDINATE
          RZ=ABS((PP(K)-P))
        ELSEIF(IFPCDNT.EQ.0)THEN                     ! FOR SIGMA COORDINATE
          RZ=ABS((PP(K)-P))*(HEIGHTU(I,J)-HEIGHTL(I,J))
        ENDIF
        DO SS=1,2
          IF(SS.EQ.1) S=U_CMPNNT
          IF(SS.EQ.2) S=V_CMPNNT
          IF(RH.LE.OBRADIUS(1,S).AND.RZ.LE.OBRADIUS(3,S))GRIDMASK(I,J,K,T,S)=0
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      OE(U_CMPNNT)=MAX(OE(U_CMPNNT),OBE(O,OS)*1.0)
      OE(V_CMPNNT)=MAX(OE(V_CMPNNT),OBE(O,OS)*1.0)
      ENDDO
    ENDDO

    O=NALLOBS
    IF(NST(NUMSTAT+1)+NST(NUMSTAT+2).GE.1) THEN
      DO T=1,FCSTGRD(4)
      DO K=1,MAXGRID(3)      !!!!!!!ATTENTION, DUE TO SPECIAL CHARACTER OF VERTICAL COORDINATE
      DO J=1,FCSTGRD(2),3
      DO I=1,FCSTGRD(1),3
        DO SS=1,2
        IF(SS.EQ.1) S=U_CMPNNT
        IF(SS.EQ.2) S=V_CMPNNT
        IF(GRIDMASK(I,J,K,T,S).EQ.1)THEN
          NP(1)=I-1
          NP(2)=J-1
          NP(3)=K-1
          NP(4)=T-1
          O=O+1
          NST(S)=NST(S)+1
          IF(O.GT.NOBSMAX)STOP 'NUMBER OF OBSERVATIONS EXCEEDED'
          DO N=1,NUMDIMS
            OBP(N,NST(S),SS)=NP(N)
          ENDDO
          OBS(NST(S),S)=0.0D0
          OBE(NST(S),S)=OE(S)
        ENDIF
        ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF
    NALLOBS=O
  ENDIF
  RETURN
  
  PRINT*, 'AFTER ADDING BACKGROUND, NALLOBS=',NALLOBS

END SUBROUTINE ADDBKGRND


END MODULE INPUT_BG_OBS
