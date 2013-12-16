MODULE OUTPUT_ANALS
!*************************************************
! OUTPUT THE ANALYSIS RESULT
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!          FEBRUARY 2008,  ZHONGJIE HE.
!
!          DECEMBER 2013, YUANFU XIE:
!          a) CHANGE ANA DIMENSION FOR SAVING MEMORY
!          b) MODIFY THE TPW CALCULATION
!*************************************************

  USE PRMTRS_STMAS
  USE GENERALTOOLS
  
  USE READ_BACKGRD   , ONLY : ORI_LON, ORI_LAT, END_LON, END_LAT,BK0
  USE DRAWCOUNTOUR   , ONLY : DRCONTOUR, DRCONTOUR_2D
  USE READOBSERVES   , ONLY : X_RADAR, Y_RADAR

  PUBLIC    OUTPTLAPS, OUTPUTANA, TMPOUTPUT
  PRIVATE   BKGMEMRLS, DRCONTOUR_0, RADIALWND
!***************************************************
!!COMMENT:
!   THIS MODULE IS USED BY THE MAIN PROGRAM, TO OUTPUT THE ANALIZED FIELDS.
!   SUBROUTINES:
!      OUTPTLAPS : OUTPUT THE ANALIZED FIELDS TO LAPS.
!      OUTPUTANA : OUTPUT THE ANALIZED FILEDS TO SOME FILES AND DRAW SOME PICTURES TO CHECK.
!      TMPOUTPUT : OUTPUT THE ANALIZED FILEDS TO SOME FILES AND DRAW SOME PICTURES TO CHECK.
!      BKGMEMRLS : RELEASE THE MEMORYS.
!      DRCONTOUR_0: WRITE SOME FIELDS INTO THE FILES WITH THE FORMAT OF SURFER( A SOFTWARE TO DRAW PICTURES).
!      RADIALWND : TRANSLATE THE U, V FIELDS TO RADIAL WIND FIELD, USED TO DRAW PICTURE.

CONTAINS

SUBROUTINE OUTPTLAPS
!*************************************************
! GET BACKGROUND FOR ANALYSIS
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,LN
  CHARACTER(LEN=200) :: DR
  REAL     :: XB(MAXGRID(1)),YB(MAXGRID(2)),ZB(MAXGRID(3)),TB(MAXGRID(4))
  REAL     :: XF(FCSTGRD(1)),YF(FCSTGRD(2)),ZF(FCSTGRD(3)),TF(FCSTGRD(4))
  REAL, ALLOCATABLE :: ANA(:,:,:,:)
  INTEGER  :: FG(MAXDIMS),MG(MAXDIMS)
  REAL     :: Z1(1),T1(1),Z2(1),T2(1),DX,DY

  ! PARAMETERS FOR CONVERTING Q, P AND T TO RH (COPIED FROM lib/degrib/rrpr.F90:
  real, parameter :: svp1=611.2
  real, parameter :: svp2=17.67
  real, parameter :: svp3=29.65
  real, parameter :: svpt0=273.15
  real, parameter :: eps = 0.622
  real, parameter :: r_dry = 287.0
  real            :: tmp,sph,density

  ! VARIABLES NEEDED FOR LAPS: YUANFU
  CHARACTER*9   :: A9
  CHARACTER*3   :: WN(3)=(/'U3','V3','OM'/)	! WIND NAMES
  CHARACTER*3   :: SN(2)=(/'SU','SV'/)		! WIND NAMES
  CHARACTER*4   :: WU(3)=(/'M/S ','M/S ','PA/S'/) ! WIND UNITS
  ! added by shuyuan 20100722
  CHARACTER*3   :: QW(2)=(/'RAI','SNO'/)       ! rain water content   snow water content
  CHARACTER*3   :: RE(1)=(/'REF'/)		! reflectivity  dbz
  CHARACTER*125 :: QWC(2)=(/'ROUR','ROUS'/) ! QW COMMENTS
  CHARACTER*125 :: RC(1)=(/'reflectivity'/) ! reflectivity COMMENTS 
  character*10   :: units_3D(2)=(/'kg/m**3','kg/m**3'/)
  !----------------------------------------------------------
  CHARACTER*125 :: WC(3)=(/'3DWIND','3DWIND','3DWIND'/) ! WIND COMMENTS
  CHARACTER*125 :: SC(2)=(/'SFCWIND','SFCWIND'/) ! SFC WIND COMMENTS
  INTEGER       :: I4,ST,IFRAME,ILV(FCSTGRD(3))  ! PRESSURE IN MB NEEDED USING HUMID WRITE ROUTINES 
  ! REAL          :: LA(FCSTGRD(1),FCSTGRD(2)),LO(FCSTGRD(1),FCSTGRD(2))
  ! REAL          :: TP(FCSTGRD(1),FCSTGRD(2)),DS,LV(FCSTGRD(3))
  REAL          :: DS,LV(FCSTGRD(3))

  ! YUANFU: MAKE ALL LARGE ARRAYS INTO ALLOCATABLE ONES FROM AUTOMATIC:
  REAL, ALLOCATABLE :: HT(:,:,:),RH(:,:,:),T3(:,:,:),W3(:,:,:,:),SF(:,:,:),SP(:,:),TPW(:,:)

  REAL          :: HEIGHT_TO_ZCOORD3,SSH2,MAKE_RH,MAKE_SSH,RM,A,DLNP
!ADDED BY SHUYUAN 20100722 FOR REFLECTIVITY
  REAL          :: REF_OUT(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  integer       :: istatus  ,N_3D_FIELDS
!Yuanfu test of cloud ice and liquid for temperature:
character :: ext*31,unit*10,comment*30 
integer   :: i4_tol,i4_ret,iqc,nref,n2d,n3d
integer   :: istatus2d(FCSTGRD(1),FCSTGRD(2)),istatus3d(FCSTGRD(1),FCSTGRD(2))
integer   :: CLOUD_BASE,CLOUD_TOP,K400
REAL :: closest_radar(FCSTGRD(1),FCSTGRD(2)),AT
REAL ::  rlat,rlon,rhgt
REAL :: FRACTION
i4_tol=900
i4_ret=0
! --------------------

  ! ALLOCATE MEMORY:
  ALLOCATE(HT(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)),RH(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)), &
           T3(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3)),W3(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),2), &
           SF(FCSTGRD(1),FCSTGRD(2),2),SP(FCSTGRD(1),FCSTGRD(2)),TPW(FCSTGRD(1),FCSTGRD(2)), &
          STAT=ST)
  
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
  
!===========
  DO K=1,FCSTGRD(3)
    ZF(K)=Z_FCSTGD(K)
  ENDDO
  DO K=1,MAXGRID(3)
    ZB(K)=Z_MAXGID(K)
  ENDDO
!=================
!  OPEN(2,FILE=DR(1:LN)//'P_LEVEL.DAT',STATUS='OLD',ACTION='READ')
!  DO K=1,MAXGRID(3)
!    READ(2,*)ZB(K)
!  ENDDO
!  CLOSE(2)
!==================
  DO T=1,FCSTGRD(4)
    TF(T)=(ITIME2(2)-ITIME2(1))*(T-1)/(FCSTGRD(4)-1)
  ENDDO
  DO T=1,MAXGRID(4)
    TB(T)=(ITIME2(2)-ITIME2(1))*(T-1)/(MAXGRID(4)-1)
  ENDDO

  ! HEIGHT FROM HYDROSTATIC:
  ! WRITE(15,*) MAXGRID,GRDBKGD0(1:MAXGRID(1),1:MAXGRID(2),1:MAXGRID(3),2,3)
  ! DO K=2,MAXGRID(3)
  !   DO J=1,MAXGRID(2)
   !    DO I=1,MAXGRID(1)
   !      GRDBKGD0(I,J,K,1:MAXGRID(4),3) = GRDBKGD0(I,J,K-1,1:MAXGRID(4),3)-287.0/9.8*0.5* &
! 	  (GRDBKGD0(I,J,K,1:MAXGRID(4),4)+GRDBKGD0(I,J,K-1,1:MAXGRID(4),4))* &
   !        (LOG(PPM(K))-LOG(PPM(K-1)))
   !    ENDDO
   !  ENDDO
  ! ENDDO
  ! WRITE(15,*) MAXGRID,GRDBKGD0(1:MAXGRID(1),1:MAXGRID(2),1:MAXGRID(3),2,3)
  print*,'Increment 1: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),1:3,1)))
  print*,'Increment 2: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),1:3,2)))
  print*,'Increment 3: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),1:3,3)))
  print*,'Increment 4: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),1:3,4)))
  print*,'Increment 5: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),1:maxgrid(3),1:3,5)))
  IF (NUMSTAT .GT. 5) print*,'Increment 6: ',maxval(ABS(grdbkgd0(1:maxgrid(1),1:maxgrid(2),maxgrid(3),1:3,6)))

  ! Yuanfu: Change Ana to a 4 dimensional array to save space:
  ! ALLOCATE MEMORY:
  ALLOCATE(ANA(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4)),STAT=ST)

  ! LOOP THROUGH ALL ANALYSIS VARIABLES:
  DO S=1,NUMSTAT
    print * ,'------------------------------------------------------'
    ANA = 0.0
    CALL BKGTOFINE(1,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,GRDBKGD0(:,:,:,:,S),ANA)

    ! MAKE SURE INTERPOLATED SH ANALYSIS POSITIVE:
    IF (S .EQ. HUMIDITY) THEN
      DO T=1,FCSTGRD(4)
      DO K=1,FCSTGRD(3)
      DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
        ANA(I,J,K,T) = &
          MAX(-BK0(I,J,K,T,HUMIDITY),ANA(I,J,K,T))
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDIF

    ! ADD INCREMENT ANA TO BK0:
    DO T=1,FCSTGRD(4)
    DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      BK0(I,J,K,T,S) = BK0(I,J,K,T,S)+ANA(I,J,K,T)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    print*,'bko_max=',maxval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1:2,S)),S
    print*,'        ',maxval(ANA(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1:2))
    print*,'bko_min=',minval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1:2,S)),S
    print*,'        ',minval(ANA(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),1:2))
  ENDDO
print*,'Specific humidity low bound: ',minval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2,5))
print*,'Specific humidity upp bound: ',maxval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2,5))

  CALL GET_GRID_SPACING_ACTUAL(LATITUDE((FCSTGRD(1)-1)/2+1,(FCSTGRD(2)-1)/2+1), &
                                LONGITUD((FCSTGRD(1)-1)/2+1,(FCSTGRD(2)-1)/2+1),DS,ST)
  ! PRESSURE LEVELS:
  CALL GET_PRES_1D(LAPSI4T,FCSTGRD(3),LV,ST)
  ILV = LV/100.0   ! INTEGER PRESSURES NEEDED IN HUMID LH3 OUTPUT
  CALL GET_R_MISSING_DATA(RM,ST)

  ! OUTPUT: WIND BY YUANFU --

  ! INTERPOLATE SURFACE WIND:
  ! HT = BK0(:,:,:,1,3)
  ! W3 = BK0(:,:,:,1,1:2)
  ! DO J=1,FCSTGRD(2)
  ! DO I=1,FCSTGRD(1)
  !   SP(I,J) = HEIGHT_TO_ZCOORD3(TP(I,J),HT,LV,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),I,J,ST)
  !   K = INT(SP(I,J))
  !   A = SP(I,J)-K
  !   SF(I,J,1:2) = (1.0-A)*W3(I,J,K,1:2)+A*W3(I,J,K+1,1:2)
  ! ENDDO
  ! ENDDO

  ! WRITE SURFACE WIND (THIS IS ALSO REQIRED TO PLOT 3D WIND ON-THE-FLY):
  ! CALL PUT_LAPS_MULTI_2D(I4,'lwm',SN,WU,SC,SF,FCSTGRD(1),FCSTGRD(2),2,ST)
  ! CALL WIND_POST_PROCESS(I4,'lw3',WN,WU,WC,W3(1,1,1,1),W3(1,1,1,2), &
  !                        FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),3,SF(1,1,1),SF(1,1,2),TP, &
  !                        LA,LO,DS,SP,RM,.TRUE.,ST)

  ! WRITE TEMPERATURE: W3 ARRAY SERVES AS WORKING ARRAY:
  ! T3(:,:,:) = BK0(:,:,:,1,4)
  ! CALL WRITE_TEMP_ANAL(I4,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),T3,HT,W3,A9,ST)

  ! OUTPUT TIME FRAME:
  IFRAME = 2	! TEMPORARY TESTING BY YUANFU XIE

  ! USE LAPS WRITE BALANCED FIELD:
  ! VERTICAL VELOCITY:
  T3 = 0.0 ! ALSO TEST FOR X Y BOUNDARIES: EXTRAPOLATE LATER YUANFU
  ! CENTER FINITE DIFFERNCE REQUIRE 2*DELTA X (or Y):
  DS = 2.0*DS
  DO K=2,FCSTGRD(3)
    DO J=2,FCSTGRD(2)-1
      DO I=2,FCSTGRD(1)-1
        ! V3(K)=V3(K-1)-DZ*(UX+VY):
        T3(I,J,K) = T3(I,J,K-1) - 0.5* &
          (LV(K)-LV(K-1))*( &
          ! (BK0(I,J,K,FCSTGRD(4),3)-BK0(I,J,K-1,FCSTGRD(4),3))*( &
          (BK0(I+1,J,K,IFRAME,1)-BK0(I-1,J,K,IFRAME,1))/DS+ &
          (BK0(I,J+1,K,IFRAME,2)-BK0(I,J-1,K,IFRAME,2))/DS )
        ! BK0(I,J,K,FCSTGRD(4),2) =  1.0e5*(BK0(I+1,J,K,FCSTGRD(4),3)-BK0(I-1,J,K,FCSTGRD(4),3))/DS
        ! BK0(I,J,K,FCSTGRD(4),1) = -1.0e5*(BK0(I,J+1,K,FCSTGRD(4),3)-BK0(I,J-1,K,FCSTGRD(4),3))/DS

	! Test vorticity:
	! T3(I,J,K) = (BK0(I+1,J,K,FCSTGRD(4),2)-BK0(I-1,J,K,FCSTGRD(4),2))/DS- &
        !             (BK0(I,J+1,K,FCSTGRD(4),1)-BK0(I,J-1,K,FCSTGRD(4),1))/DS
        ! T3(I,J,K) = T3(I,J,K)*1.0e4
      ENDDO
    ENDDO
  ENDDO

  ! According a discussion with Dan, specific humidity is adjusted by Q_r (rain content) using ssh2 routine of LAPS:
GOTO 111
  DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
        ! Assume saturation: TD = T:
        IF (BK0(I,J,K,IFRAME,6) .GT. 0.0) &
          BK0(I,J,K,IFRAME,5) = SSH2(LV(k)/100.0,BK0(I,J,K,IFRAME,4)-273.15,BK0(I,J,K,IFRAME,4)-273.15,-132.0)
      ENDDO
    ENDDO
  ENDDO
111 continue ! skip SH adjustment according to q_r

  ! Adjust SH by bounds:
goto 222
  do k=1,fcstgrd(3)
  do j=1,fcstgrd(2)
  do i=1,fcstgrd(1)
    bk0(i,j,k,iframe,5) = max(bk0(i,j,k,iframe,5),bk0(i,j,k,iframe,6))
  enddo
  enddo
  enddo
222 continue

  ! CONVERTED FROM P, Q T: NOTE: Q is in g/kg
  DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
        ! sph = BK0(I,J,K,IFRAME,5)*0.001
        ! tmp = BK0(I,J,K,IFRAME,4)
        ! RH(I,J,K) = 1.E2 * (LV(k)*sph/(sph*(1.-eps) + eps))/(svp1*exp(svp2*(tmp-svpt0)/(tmp-svp3)))
        IF (BK0(I,J,K,IFRAME,5) .GE. 0.0) THEN
          RH(I,J,K) = MAKE_RH(LV(k)/100.0,BK0(I,J,K,IFRAME,4)-273.15,BK0(I,J,K,IFRAME,5),-132.0)

          ! Ensure no RH greater than 1.0:
          IF (RH(I,J,K) .GT. 1.0) THEN
            RH(I,J,K) = 1.0
            BK0(I,J,K,IFRAME,5) = &
              MAKE_SSH(LV(k)/100.0,BK0(I,J,K,IFRAME,4)-273.15,1.0,-132.0)
          ENDIF

          ! Percentage:
          RH(I,J,K) = RH(I,J,K)*100.0
        ELSE
          RH(I,J,K) = 0.0
        ENDIF
      ENDDO
    ENDDO
  ENDDO
  print*,'Max/Min RH: ',maxval(RH),minval(RH)

  ! HEIGHT FROM HYDROSTATIC:
  ! BK0(1:FCSTGRD(1),1:FCSTGRD(2),1,IFRAME,3) = 0.0
  goto 1
  DO K=4,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
        BK0(I,J,K,IFRAME,3) = BK0(I,J,K-1,IFRAME,3)-287.0/9.8*0.5* &
	  (BK0(I,J,K,IFRAME,4)+BK0(I,J,K-1,IFRAME,4))* &
          (LOG(LV(K))-LOG(LV(K-1)))
      ENDDO
    ENDDO
  ENDDO
1 continue

  ! Temperature adjustment according to rain and snow: Testing by YUANFU
  ! SKIP FOR NOW AS CWB DOMAIN FOR MORAKOT, THIS CAUSES VERY LARGE TEMP INCREMENT
  GOTO 11
  DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
        ! Dry density:
        density = LV(K)/r_dry/BK0(I,J,K,IFRAME,4)
        ! Adjusted temperature:
        BK0(I,J,K,IFRAME,4) = LV(K)/r_dry/ &
          (density+BK0(I,J,K,IFRAME,6)*0.001+BK0(I,J,K,IFRAME,7)*0.001)
      ENDDO
    ENDDO
  ENDDO
11 CONTINUE
  ! Temperature adjustment according to cloud ice and liquid: Test by Yuanfu
  ext = "vrz"
  BK0(:,:,:,1,NUMSTAT+1) = 0.0      ! Reflectivity
  call read_multiradar_3dref(lapsi4t, &
      900,0,      & ! 900 tolerate
      .true.,-10.0, & ! apply_map: true; ref missing data value: 10
      fcstgrd(1),fcstgrd(2),fcstgrd(3),ext,latitude,longitud,topogrph, &
      .false.,.false., & ! l_low_fill: false; l_high_fill: false
      bk0(1,1,1,iframe,3),bk0(1,1,1,1,numstat+1),rlat,rlon,rhgt,unit,iqc,closest_radar, &
      nref,n2d,n3d,istatus2d,istatus3d)
  ext = "lwc"
  BK0(:,:,:,IFRAME,NUMSTAT+1) = 0.0 ! Cloud liquid
  BK0(:,:,:,IFRAME,NUMSTAT+2) = 0.0 ! Cloud ice
  CALL GET_LAPS_3DGRID(LAPSI4T,i4_tol,i4_ret, &
               FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),ext,"lwc", &
               unit,comment,BK0(1,1,1,IFRAME,NUMSTAT+1),ST)
  CALL GET_LAPS_3DGRID(LAPSI4T,i4_tol,i4_ret, &
               FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),ext,"ice", &
               unit,comment,BK0(1,1,1,IFRAME,NUMSTAT+2),ST)
  print*,'Max cloud liquid: ',maxval(BK0(:,:,:,IFRAME,NUMSTAT+1))
  print*,'Max cloud ice   : ',maxval(BK0(:,:,:,IFRAME,NUMSTAT+2))
 goto 333
  ! Adjust temperature according to cloud:
  DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      AT=0.0
      CLOUD_BASE = FCSTGRD(3)+1
      CLOUD_TOP = 0
      DO K=1,FCSTGRD(3)
        IF (LV(K) .EQ. 40000) THEN 
          IF (BK0(I,J,K,1,NUMSTAT+1) .GT. 5 .AND. &
              BK0(I,J,K,1,NUMSTAT+1) .LT.100.0 ) &
            AT=BK0(I,J,K,1,NUMSTAT+1)/10.0
          K400 = K
        ENDIF
        IF (BK0(I,J,K,IFRAME,NUMSTAT+1) .GT. 0.0 .OR. &
            BK0(I,J,K,IFRAME,NUMSTAT+2) .GT. 0.0) THEN
          CLOUD_BASE = MIN0(K,CLOUD_BASE)
          CLOUD_TOP  = MAX0(K,CLOUD_TOP)
        ENDIF
      ENDDO

      ! Only adjust when cloud top is above 400mb
      IF (CLOUD_TOP .LT. K400 .OR. AT .LE. 0.0) cycle

      IF (CLOUD_BASE .LE. K400) THEN
        DO K=CLOUD_BASE,K400
          BK0(I,J,K,IFRAME,4) = BK0(I,J,K,IFRAME,4)+ &
            AT*(K-CLOUD_BASE)/FLOAT(MAX(K400-CLOUD_BASE,1))
        ENDDO
      ENDIF
      DO K=K400+1,CLOUD_TOP
        BK0(I,J,K,IFRAME,4) = BK0(I,J,K,IFRAME,4)+ &
          AT*(CLOUD_TOP-K)/FLOAT(MAX(CLOUD_TOP-K400,1))
      ENDDO
    ENDDO
  ENDDO
333 continue
 
  ! TOTAL PRECIPITABLE WATER:
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    TPW(I,J) = 0.0
    DO K=1,FCSTGRD(3)-1
      IF (LV(K)/100.0 .LE. p_sfc_f(i,j)) THEN
        ! above topography: summed up
        TPW(I,J) = TPW(I,J) + 0.5* &
                 (BK0(I,J,K,IFRAME,5)+BK0(I,J,K+1,IFRAME,5))* &
                 (LV(K)-LV(K+1))/100.0 ! PRESSURE IN MB

      ELSEIF (LV(K+1)/100.0 .LE. p_sfc_f(i,j)) THEN
        TPW(I,J) = TPW(I,J) + BK0(I,J,K+1,IFRAME,5)* &
                   (p_sfc_f(i,j)-LV(k+1)/100.0)
      ENDIF
    ENDDO
    ! FROM G/KG TO CM:
    TPW(I,J) = TPW(I,J)/100.0/9.8 ! FOLLOWING DAN'S INT_IPW.F ROUTINE
  ENDDO
  ENDDO
  PRINT*,'TPW max/min: ', maxval(TPW),minval(TPW)

  ! KG/KG FOR BALANCE NETCDF:
  BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),IFRAME,5) = 0.001* &
    BK0(1:FCSTGRD(1),1:FCSTGRD(2),1:FCSTGRD(3),IFRAME,5)
  CALL WRITE_BAL_LAPS(LAPSI4T,BK0(1,1,1,IFRAME,3),BK0(1,1,1,IFRAME,1), &
                          BK0(1,1,1,IFRAME,2),BK0(1,1,1,IFRAME,4), &
                          T3,RH,BK0(1,1,1,IFRAME,5),FCSTGRD(1),FCSTGRD(2), &
                          FCSTGRD(3),LV,ST)

  ! WRITE TEMPERATURE AND WIND INTO LAPSPRD:
  CALL WRITE_TEMP_ANAL(LAPSI4T,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3), &
       BK0(1,1,1,IFRAME,4),BK0(1,1,1,IFRAME,3),DR,ST)

  ! WRITE WIND BY COMPONENTS:
  CALL PUT_LAPS_MULTI_3D_JACKET(LAPSI4T,'lw3',WN,WU,WC, &
                     BK0(1,1,1,IFRAME,1),BK0(1,1,1,IFRAME,2),T3, &
                     FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),3,ST)

  ! WRITE SP TO FILE FOLLOWING DAN'S LQ3_DRIVE1A.f:
  WRITE(*,*) 'WRITE-FILE...'
  CALL WRITEFILE(LAPSI4T,'STMAS SH analysis',ILV,BK0(1,1,1,IFRAME,5), &
                 FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),ST)

  ! WRITE TOTAL PRECIPITABLE WATER TO lh4 FILE:
  WRITE(*,*) 'WRITE-LH4...'
  CALL WRITE_LH4(LAPSI4T,TPW,1.0,FCSTGRD(1),FCSTGRD(2),ST)

  ! WRITE RH TO LH3:
  WRITE(*,*) 'WRITE-LH3...'
  CALL LH3_COMPRESS(BK0(1,1,1,IFRAME,5),BK0(1,1,1,IFRAME,4),LAPSI4T,ILV, &
                    -132.0,FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),0,ST)

  ! END OF OUTPUT TO LAPS BY YUANFU

! --------ADDED BY SHUYUAN 201007---------------------------------
! CALCULATE AND OUTPUT REFLECTIVITY
  ! CHECK IF RAIN AND SNOW IS ANALYZED:
  IF (NUMSTAT .LE. 5) GOTO 555 ! BY YUANFU
   T=2
   DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
      DO I=1,FCSTGRD(1)
      REF_OUT(I,J,K)= 0.
      if(BK0(I,J,K,T,6).GT. 0.0 ) then 
       !20100907
      ! REF_OUT(I,J,K)=REF_OUT(I,J,K)+(BK0(I,J,K,T,6))**1.75*17300.
        REF_OUT(I,J,K)=REF_OUT(I,J,K)+43.1+17.5*ALOG10(BK0(I,J,K,T,6))           
      else
       ! LAPS uses -10 as base value for reflectivity:
       REF_OUT(I,J,K)=-10.
      endif 
 
      if( REF_OUT(I,J,K) .LT. 0.) then
        ! LAPS uses -10 as base value for reflectivity:
        REF_OUT(I,J,K)=-10.
      endif    
      ENDDO
    ENDDO
   ENDDO 
   print*,'ref_max=',maxval(REF_OUT(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3)))
   print*,'ref_min=',minval(REF_OUT(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3)))
   call put_laps_3d(LAPSI4T,'lps',RE,'dBZ',RC,REF_OUT(1,1,1),FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
             !RC,BK0(1,1,1,T,10),FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
             
  !RAIN CONTENT(AIR DENSITY* RAIN WATER MIXING RATIO(kg/m3),
  !SNOW CONTENT(AIR DENSITY * SNOW WATER MIXING RATIO)
   N_3D_FIELDS=2   ! variable number
   istatus=0 
  !convert rc sc to rour rous 
   DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
     DO I=1,FCSTGRD(1)
     if(BK0(I,J,K,T,6) .NE. 0)then
     !BK0(I,J,K,T,6)=(BK0(I,J,K,T,6)/17300.)**(4./7.)/1000
      BK0(I,J,K,T,6)=BK0(I,J,K,T,6)/1000.  !20100907  kg/m3
     endif
     if(BK0(I,J,K,T,7) .NE. 0)then
     !BK0(I,J,K,T,7)=(BK0(I,J,K,T,7)/38000.)**(5./11.)/1000
      BK0(I,J,K,T,7)=BK0(I,J,K,T,7)/1000.   !20100907
     endif
     ENDDO
    ENDDO
   ENDDO  
   print*,'bko6_max=',maxval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2,6))
   print*,'bko7_max=',maxval(BK0(1:fcstgrd(1),1:fcstgrd(2),1:fcstgrd(3),2,7))
   ! OUT PUT SNOW CONTENT(SNO) AND RAI   
   call put_laps_3d_multi_R(LAPSI4T,'lwc',QW,units_3D,QWC ,  &
              BK0(1,1,1,T,6),BK0(1,1,1,T,7),     &                                                             
              FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),            &
              FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),N_3D_FIELDS,istatus)         
!--END OUTPUT REF RAI SNO  BY SHUYUAN---------------------

  ! SKIP OUTPUT RAIN AND SNOW IF NOT ANALYZED:
555 CONTINUE

  CALL BKGMEMRLS
  DEALLOCATE(Z_FCSTGD)
  ! FINALLY RELEASE THE BACKGROUND ARRAYS: BY YUANFU
  DEALLOCATE(BK0)
  ! FINALLY RELEASE MEMORY FOR LAT/LON/TOPO:
  DEALLOCATE(LATITUDE, LONGITUD, TOPOGRPH)
  RETURN
END SUBROUTINE OUTPTLAPS

SUBROUTINE OUTPUTANA
!*************************************************
! GET BACKGROUND FOR ANALYSIS
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!          FEBRUARY 2008, BY ZHONGJIE HE
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,L,T,S,LN
  CHARACTER(LEN=200) :: DR
  REAL     :: XB(MAXGRID(1)),YB(MAXGRID(2)),ZB(MAXGRID(3)),TB(MAXGRID(4))
  REAL     :: XF(FCSTGRD(1)),YF(FCSTGRD(2)),ZF(FCSTGRD(3)),TF(FCSTGRD(4))
  REAL     :: BK1(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4),NUMSTAT)
  INTEGER  :: FG(MAXDIMS),MG(MAXDIMS)
  REAL     :: Z1(1),T1(1),Z2(1),T2(1),DX,DY,DT

! ADDED BY ZHONGJIE HE==========
  REAL  :: UU(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL  :: VV(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL  :: WW(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL  :: UV(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3))
  REAL  :: W0(FCSTGRD(1),FCSTGRD(2),FCSTGRD(3),FCSTGRD(4))

  REAL  :: D1,D2,DD
!===============================

! --------------------
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
  CALL GET_DIRECTORY('static',DR,LN)
!================
  DO K=1,FCSTGRD(3)
    ZF(K)=Z_FCSTGD(K)
  ENDDO
  DO K=1,MAXGRID(3)
    ZB(K)=Z_MAXGID(K)
  ENDDO
!=============== JUST FOR TEST DATA 
!  OPEN(2,FILE=DR(1:LN)//'P_LEVEL.DAT',STATUS='OLD',ACTION='READ')
!  DO K=1,MAXGRID(3)
!    READ(2,*)ZB(K)
!  ENDDO
!  CLOSE(2)
  
!  DO K=1,MAXGRID(3)
!    ZB(K)=ZF(1)+(K-1)*(ZF(FCSTGRD(3))-ZF(1))/(MAXGRID(3)-1.0)
!  ENDDO
! ===============

  DO T=1,FCSTGRD(4)
    TF(T)=(T-1)*1.0
  ENDDO
  IF(MAXGRID(4).GE.2) DT=((FCSTGRD(4)-1)*1.0)/(MAXGRID(4)-1.0)
  DO T=1,MAXGRID(4)
    TB(T)=TF(1)+(T-1)*DT
  ENDDO

! ADDED BY ZHONGJIE HE==========

!  CALL FCST2BKGD(NUMDIMS,NGPTOBS,NUMSTAT,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,DIFFTOUT,BK1)
  CALL BKGTOFINE(NUMSTAT,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,GRDBKGD0,BK1)
!  CALL FCST2BKGD(NUMDIMS,NGPTOBS,1,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,WWWOUT,W0)
  CALL BKGTOFINE(1,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,WWWOUT,W0)

  DEALLOCATE(DIFFTOUT)
  DEALLOCATE(WWWOUT)

  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    UU(I,J,K)=BK1(I,J,K,1,U_CMPNNT)
    VV(I,J,K)=BK1(I,J,K,1,V_CMPNNT)
!    WW(I,J,K)=W0(I,J,K,1)
!    D1=ORI_LON+(END_LON-ORI_LON)/FLOAT(FCSTGRD(1)-1)*(I-1)
!    D1=D1-X_RADAR
!    D2=ORI_LAT+(END_LAT-ORI_LAT)/FLOAT(FCSTGRD(2)-1)*(J-1)
!    D2=D2-Y_RADAR
!    DD=SQRT(D1*D1+D2*D2)
!    UV(I,J,K)=UU(I,J,K)*D1/DD+VV(I,J,K)*D2/DD
  ENDDO
  ENDDO
  ENDDO

  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,UU,'UD.GRD')
  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,VV,'VD.GRD')
!  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,WW,'WD.GRD')
!  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,UV,'RD.GRD')

!===============================

!  CALL FCST2BKGD(NUMDIMS,NGPTOBS,NUMSTAT,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,GRDBKGD0,BK1)
  CALL BKGTOFINE(NUMSTAT,MAXGRID,XB,YB,ZB,TB,FCSTGRD,XF,YF,ZF,TF,GRDBKGD0,BK1)

! ADDED BY ZHONGJIE HE==========
  DO K=1,FCSTGRD(3)
  DO J=1,FCSTGRD(2)
  DO I=1,FCSTGRD(1)
    UU(I,J,K)=BK1(I,J,K,1,U_CMPNNT)
    VV(I,J,K)=BK1(I,J,K,1,V_CMPNNT)
    WW(I,J,K)=W0(I,J,K,1)
!    D1=ORI_LON+(END_LON-ORI_LON)/FLOAT(FCSTGRD(1)-1)*(I-1)
!    D1=D1-X_RADAR
!    D2=ORI_LAT+(END_LAT-ORI_LAT)/FLOAT(FCSTGRD(2)-1)*(J-1)
!    D2=D2-Y_RADAR
!    DD=SQRT(D1*D1+D2*D2)
!    UV(I,J,K)=UU(I,J,K)*D1/DD+VV(I,J,K)*D2/DD
  ENDDO
  ENDDO
  ENDDO

  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,UU,'UA.GRD')
  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,VV,'VA.GRD')
  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,WW,'WA.GRD')
!  CALL DRCONTOUR(ORI_LON,END_LON,ORI_LAT,END_LAT,NUMDIMS,FCSTGRD,UV,'RA.GRD')
!===============================

  CALL BKGMEMRLS
  DEALLOCATE(Z_FCSTGD)
  DEALLOCATE(Z_MAXGID)


!============= to check the analysis field======
  OPEN(1,FILE='STMAS_ANAL.DAT',ACTION='WRITE')
    WRITE(1,*) FCSTGRD(1:4)
    DO L=1,FCSTGRD(4)
    DO K=1,FCSTGRD(3)
    DO J=1,FCSTGRD(2)
    DO I=1,FCSTGRD(1)
      WRITE(1,*) BK1(I,J,K,L,U_CMPNNT),BK1(I,J,K,L,V_CMPNNT),BK1(I,J,K,L,TEMPRTUR),W0(I,J,K,L)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  CLOSE(1)
!===========================================


  RETURN
END SUBROUTINE OUTPUTANA

SUBROUTINE TMPOUTPUT
!*************************************************
! OUTPUT ANALYSIS
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S
  INTEGER  :: IC,JC
! --------------------
  OPEN(2,FILE='ANALYSIS.DAT')
  DO S=1,NUMSTAT
    DO T=1,1      ! NUMGRID(4)
    DO K=1,1      ! NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      WRITE(2,*) GRDBKGD0(I,J,K,T,S)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  CLOSE(2)
!  CALL DRCONTOUR_0(PRESSURE,'Z.GRD')
  CALL DRCONTOUR_0(U_CMPNNT,'U.GRD')
  CALL DRCONTOUR_0(V_CMPNNT,'V.GRD')
!  CALL DRCONTOUR_0(4,'T.GRD')
  IC=(NUMGRID(1)+1)/2
  JC=(NUMGRID(2)+1)/2
  CALL RADIALWND(IC,JC)
  CALL BKGMEMRLS
  RETURN
END SUBROUTINE TMPOUTPUT

SUBROUTINE BKGMEMRLS
!*************************************************
! RELEASE MEMORY FOR BACKGROUND ARRAY
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
  DEALLOCATE(GRDBKGD0)
  DEALLOCATE(XX0)
  DEALLOCATE(YY0)
  DEALLOCATE(CR0)
  DEALLOCATE(DG0)
  DEALLOCATE(DN0)
  IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN
    DEALLOCATE(ZZ0)
    DEALLOCATE(ZZB)
  ELSEIF(IFPCDNT.EQ.1)THEN
    DEALLOCATE(PP0,PPM)
  ENDIF
  RETURN
END SUBROUTINE BKGMEMRLS

SUBROUTINE DRCONTOUR_0(S,FN)
!*************************************************
! DRAW THE ANALYSIS (AFFILIATE)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER   :: I,J,K,T,S,IU
  REAL      :: MX,MN
  CHARACTER(LEN=5) :: FN
! --------------------
  K=2 !(NUMGRID(3)+1)/2 !NUMGRID(3)
  T=(NUMGRID(4)+1)/2 !NUMGRID(4) !1
  IU=2
  OPEN(IU,FILE=FN,STATUS='UNKNOWN')
  MX=-1000.0
  MN=1000.0
  DO I=1,NUMGRID(1)
  DO J=1,NUMGRID(2)
    IF(GRDBKGD0(I,J,K,T,S).LT.9000.0)THEN
      IF(GRDBKGD0(I,J,K,T,S).GT.MX)MX=GRDBKGD0(I,J,K,T,S)
      IF(GRDBKGD0(I,J,K,T,S).LT.MN)MN=GRDBKGD0(I,J,K,T,S)
    ENDIF
  ENDDO
  ENDDO
  WRITE(IU,'(A4)')'DSAA'
  WRITE(IU,'(2I4)')NUMGRID(1),NUMGRID(2)
  WRITE(IU,*)ORIPSTN(1),ORIPSTN(1)+(NUMGRID(1)-1)*GRDSPAC(1)
  WRITE(IU,*)ORIPSTN(2),ORIPSTN(2)+(NUMGRID(2)-1)*GRDSPAC(2)
  WRITE(IU,*)MN,MX
  DO I=1,NUMGRID(1)
  DO J=1,NUMGRID(2)
    IF(GRDBKGD0(I,J,K,T,S).GT.9000.0)GRDBKGD0(I,J,K,T,S)=2.E38
  ENDDO
  ENDDO
  DO J=1,NUMGRID(2)
    WRITE(IU,*)(GRDBKGD0(I,J,K,T,S),I=1,NUMGRID(1))
  ENDDO
  CLOSE(IU)
  RETURN
END SUBROUTINE DRCONTOUR_0

SUBROUTINE RADIALWND(IC,JC)
!*************************************************
! DRAW RADIAL WIND (AFFILIATE)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: IU,IC,JC,I,J,K,T
  REAL     :: DX,DY,DD,MN,MX
  REAL     :: UV(NUMGRID(1),NUMGRID(2),NUMGRID(3),NUMGRID(4))
! --------------------
  DO T=1,NUMGRID(4)
  DO K=1,NUMGRID(3)
  DO J=1,NUMGRID(2)
  DO I=1,NUMGRID(1)
    IF(I.EQ.IC.AND.J.EQ.JC)THEN
      UV(I,J,K,T)=2.E38
      CYCLE
    ENDIF
    DX=1.0D0*(I-IC)
    DY=1.0D0*(J-JC)
    DD=SQRT(DX*DX+DY*DY)
    UV(I,J,K,T)=GRDBKGD0(I,J,K,T,U_CMPNNT)*DX/DD+GRDBKGD0(I,J,K,T,V_CMPNNT)*DY/DD
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  K=(NUMGRID(3)+1)/2 !NUMGRID(3)
  T=(NUMGRID(4)+1)/2 !NUMGRID(4) !1
  IU=2
  OPEN(IU,FILE='RADIAL.GRD',STATUS='UNKNOWN')
  MX=-1000.0
  MN=1000.0
  DO I=1,NUMGRID(1)
  DO J=1,NUMGRID(2)
    IF(UV(I,J,K,T).LT.9000.0)THEN
      IF(UV(I,J,K,T).GT.MX)MX=UV(I,J,K,T)
      IF(UV(I,J,K,T).LT.MN)MN=UV(I,J,K,T)
    ENDIF
  ENDDO
  ENDDO
  WRITE(IU,'(A4)')'DSAA'
  WRITE(IU,'(2I4)')NUMGRID(1),NUMGRID(2)
  WRITE(IU,*)ORIPSTN(1),ORIPSTN(1)+(NUMGRID(1)-1)*GRDSPAC(1)
  WRITE(IU,*)ORIPSTN(2),ORIPSTN(2)+(NUMGRID(2)-1)*GRDSPAC(2)
  WRITE(IU,*)MN,MX
  DO I=1,NUMGRID(1)
  DO J=1,NUMGRID(2)
    IF(UV(I,J,K,T).GT.9000.0)UV(I,J,K,T)=2.E38
  ENDDO
  ENDDO
  DO J=1,NUMGRID(2)
    WRITE(IU,*)(UV(I,J,K,T),I=1,NUMGRID(1))
  ENDDO
  CLOSE(IU)
  RETURN
END SUBROUTINE RADIALWND



!!added by shuyuan 20100729!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine put_laps_3d_multi_R(i4time,EXT,var_3d &  ! in 
             ,units_3d,comment_3d       &            ! in    
             ,array1,array2              &           ! in                                                 
             ,NX_L_1,NY_L_1,NZ_L_1     &             ! in   
             ,NX_L_2,NY_L_2,NZ_L_2      &            ! in                             
             ,N_3D_FIELDS,istatus)              !in out

   real array1(NX_L_1,NY_L_1,NZ_L_1) !variable for output 
   real array2(NX_L_2,NY_L_2,NZ_L_2) !variable for output
       
   character*125 comment_3D(N_3D_FIELDS) !variable comment 
   character*10 units_3D(N_3D_FIELDS)    !variable unit
   character*3 var_3D(N_3D_FIELDS)      !variable name
   character*(*) EXT                    ! file name suffix and file derictory
   integer :: l

   write(6,*)' Subroutine put_laps_3d_multi_R...'

   l = 1
   call put_laps_multi_3d(i4time,EXT,var_3d(l),units_3d(l), &
        comment_3d(l),array1,NX_L_1,NY_L_1,NZ_L_1,1,istatus)       
   if(istatus .ne. 1)return
   if(l .eq. N_3D_FIELDS)return

   l = 2
   call put_laps_multi_3d_append(i4time,EXT,var_3d(l),units_3d(l),   &    
        comment_3d(l),array2,NX_L_2,NY_L_2,NZ_L_2,1,istatus)       
   if(istatus .ne. 1)return
   if(l .eq. N_3D_FIELDS)return

   write(6,*)' Error: N_3D_FIELDS exceeds limit ',N_3D_FIELDS

   return
 end subroutine put_laps_3d_multi_R

!!!!from cloud_deriv_subs!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 subroutine put_laps_multi_3d_append(i4time,EXT,var_2d,units_2d, &
                         comment_2d,field_3d,ni,nj,nk,nf,istatus)

        logical ltest_vertical_grid

        character*150 DIRECTORY
        character*(*) EXT

        character*125 comment_3d(nk*nf),comment_2d(nf)
        character*10 units_3d(nk*nf),units_2d(nf)
        character*3 var_3d(nk*nf),var_2d(nf)
        integer LVL_3d(nk*nf)
        character*4 LVL_COORD_3d(nk*nf)

        real field_3d(ni,nj,nk,nf)

        istatus = 0

        call get_directory(ext,directory,len_dir)

        do l = 1,nf
            write(6,11)directory,ext(1:5),var_2d(l)
11          format(' Writing 3d ',a50,1x,a5,1x,a3)
        enddo ! l

        do l = 1,nf
          do k = 1,nk

            iscript_3d = (l-1) * nk + k

            units_3d(iscript_3d)   = units_2d(l)
            comment_3d(iscript_3d) = comment_2d(l)
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(iscript_3d) = zcoord_of_level(k)/10
                lvl_coord_3d(iscript_3d) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then

                lvl_3d(iscript_3d) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(iscript_3d) = 'HPA'
            else
                write(6,*)' Error, vertical grid not supported,'  &
                         ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(iscript_3d) = var_2d(l)

          enddo ! k
        enddo ! l

        CALL WRITE_LAPS_MULTI(I4TIME,DIRECTORY,EXT,ni,nj, &
        nk*nf,nk*nf,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,  &
                           COMMENT_3D,field_3d,ISTATUS)

        if(istatus .ne. 1)return

        istatus = 1

        return
        end subroutine put_laps_multi_3d_append
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE OUTPUT_ANALS
