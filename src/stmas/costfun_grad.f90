MODULE COSTFUN_GRAD
!*************************************************
! CALCULATE COST FUNCTION AND THE RELEVENT GRADIENTS
! HISTORY: DEPARTED FROM STMAS4D_CORE MODULE BY ZHONGJIE HE, AUGUST 2007.
!*************************************************

  USE PRMTRS_STMAS
  USE WCOMPT_GRADT
  USE SMCOSTF_GRAD
  USE GSBCOST_GRAD
  USE HYDCOST_GRAD,ONLY:HYDROCOST,HYDROGRAD,HYDROCOST_XIE,HYDROGRAD_XIE,HYDROCOST_SHUYUAN,HYDROGRAD_SHUYUAN

  PUBLIC    COSTFUNCT, COSTGRADT, COSTFUNCT1, COSTGRADT1, COSTFUNCT2, COSTGRADT2, CG_VALUE, CG_GRAD

!***************************************************
!!COMMENT:
!   THIS MODULE IS USED BY THE MODULE OF STMAS4D_CORE, THE FUNCTION IS TO CALCULATE THE VALUE OF COSTFUNCTION AND GRADIENTS.
!   SUBROUTINES:
!      COSTFUNCT: CALCULATE THE VALUE OF COSTFUNCTION, WHILE FOR THE RADIAL WIND DATA, THE VERTICAL VELOCITY IS NOT CONSIDERED.
!      COSTGRADT: CALCULATE THE GRADIENT OF COSTFUNCTIONS CORRESPOND TO COSTFUNCT.
!      COSTFUNCT1: CALCULATE THE VALUE OF COSTFUNCTION, FOR THE CASE THAT VERTICAL VELOCITYS ARE STATE VARIABLES.
!      COSTGRADT1: CALCULATE THE GRADIENT OF COSTFUNCTIONS CORRESPOND TO COSTFUNCT1.
!      COSTFUNCT2: JUST LIKE THE SUBROUTINE OF COSTFUNCT, WHILE THE VERTICAL VELOCITY IS CONSIDERED, BUT NOT STATE VARIABLES.
!      COSTGRADT2: CALCULATE THE GRADIENT OF COSTFUNCTIONS CORRESPOND TO COSTFUNCT2.
!***************************************************
CONTAINS


SUBROUTINE COSTFUNCT
!*************************************************
! CALCULATE COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV
  INTEGER  :: NP(MAXDIMS)
  REAL  :: HT,OI,HU,HV,DG
  REAL  :: HT0,HU0,HV0
  REAL  :: CS(NUMSTAT),CM,CB
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
! INITIALIZATION
  COSTFUN=0.0
  DO S=1,NUMSTAT
    CS(S)=0.0
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))/(NOBSTAT(S)*1.0)
      CS(S)=CS(S)+0.5*(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))/(NOBSTAT(S)*1.0)
    ENDDO
  ENDDO
  S=NUMSTAT+1                        !  FOR RADAR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    HT=HU*SIND(DG)+HV*COSD(DG)
    COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))*OBSRADAR
  ENDDO
  S=NUMSTAT+2                       ! FOR SFMR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HU0=0.0
    HV0=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0))
    HT0=SQRT(HU0*HU0+HV0*HV0)
    COSTFUN=COSTFUN+(HT-HT0-OBSVALUE(O))*OI*(HT-HT0-OBSVALUE(O))*OBS_SFMR
  ENDDO

  COSTFUN=0.5D0*COSTFUN

! SMOOTH TERM
  CM=COSTFUN
  CALL SMOTHCOST
  CM=COSTFUN-CM

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    CB=COSTFUN
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN      ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNCOST_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                    ! FOR PRESSURE COORDINATE
      CALL GSBLNCOST_P_HS
    ENDIF
    CB=COSTFUN-CB
  ENDIF

!  WRITE(100,*)(CS(S),S=1,NUMSTAT),CM,CB
  RETURN
END SUBROUTINE COSTFUNCT

SUBROUTINE COSTGRADT
!*************************************************
! CALCULATE GRADIENT OF COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV
  INTEGER  :: NP(MAXDIMS)
  REAL     :: HT,OI,HU,HV,DG
  REAL     :: HT0,HU0,HV0
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
! INITIALIZATION
  DO S=1,NUMSTAT
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      GRADINT(I,J,K,T,S)=0.0
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        GRADINT(I,J,K,T,S)=GRADINT(I,J,K,T,S)  &
        +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)/(NOBSTAT(S)*1.0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDDO
  S=NUMSTAT+1                                    !  FOR RADAR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    HT=HU*SIND(DG)+HV*COSD(DG)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
      +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*SIND(DG)*OBSRADAR
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
      +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*COSD(DG)*OBSRADAR
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  S=NUMSTAT+2                       !  FOR SFMR OBSERVATION
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)*OBSERROR(O)
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HU0=0.0
    HV0=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0))
    HT0=SQRT(HU0*HU0+HV0*HV0)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
      +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HU+HU0)/HT*OBS_SFMR
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
      +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HV+HV0)/HT*OBS_SFMR
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

! SMOOTH TERM
  CALL SMOTHGRAD

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNGRAD_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                     ! FOR PRESSURE COORDINATE
      CALL GSBLNGRAD_P_HS
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE COSTGRADT

SUBROUTINE COSTFUNCT1
!*************************************************
! CALCULATE COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV,WW
  INTEGER  :: NP(MAXDIMS)
  REAL  :: HT,OI,HU,HV,HW,DG,DV
  REAL  :: HT0,HU0,HV0,HW0
  REAL  :: CS(NUMSTAT),CM,CB,CH
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
  WW=W_CMPNNT
! INITIALIZATION
  COSTFUN=0.0
  DO S=1,NUMSTAT
    CS(S)=0.0
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))/(NOBSTAT(S)*1.0)
      CS(S)=CS(S)+0.5*(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))/(NOBSTAT(S)*1.0)
    ENDDO
  ENDDO
  S=NUMSTAT+1                                  !   FOR RADAR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    DV=OBSEINF2(NO)
    HT=HU*SIND(DG)*COSD(DV)+HV*COSD(DG)*COSD(DV)+HW*SIND(DV)
    COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))*OBSRADAR
  ENDDO
  S=NUMSTAT+2                        ! FOR SRMR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    HU0=0.0
    HV0=0.0
    HW0=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,WW)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
      HW0=HW0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0)+(HW+HW0)*(HW+HW0))
    HT0=SQRT(HU0*HU0+HV0*HV0+HW0*HW0)
    COSTFUN=COSTFUN+(HT-HT0-OBSVALUE(O))*OI*(HT-HT0-OBSVALUE(O))*OBS_SFMR
  ENDDO

  COSTFUN=0.5D0*COSTFUN

! SMOOTH TERM
  CM=COSTFUN
  CALL SMOTHCOST
  CM=COSTFUN-CM

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    CB=COSTFUN
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN           ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNCOST_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                         ! FOR PRESSURE COORDINATE
!    CALL GSBLNCOST_P_HS
      CALL CNTNSCOST
    ENDIF
    CB=COSTFUN-CB
  ENDIF

!  WRITE(100,*)(CS(S),S=1,NUMSTAT),CM,CB

! HYDROSTATIC CONDITION TERM                 ! ADDED BY ZHONGJIE HE
  IF(GRDLEVL .LE. ENDHYLV) THEN
    CH=COSTFUN
    CALL HYDROCOST
    CH=COSTFUN-CH
  ENDIF

  RETURN
END SUBROUTINE COSTFUNCT1

SUBROUTINE COSTGRADT1
!*************************************************
! CALCULATE GRADIENT OF COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV,WW
  INTEGER  :: NP(MAXDIMS)
  REAL     :: HT,OI,HU,HV,HW,DG,DV
  REAL     :: HT0,HU0,HV0,HW0
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
  WW=W_CMPNNT
! INITIALIZATION
  DO S=1,NUMSTAT
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      GRADINT(I,J,K,T,S)=0.0
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        GRADINT(I,J,K,T,S)=GRADINT(I,J,K,T,S)  &
        +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)/(NOBSTAT(S)*1.0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDDO
  S=NUMSTAT+1                         !  FOR RADAR DATA
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    DV=OBSEINF2(NO)
    HT=HU*SIND(DG)*COSD(DV)+HV*COSD(DG)*COSD(DV)+HW*SIND(DV)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
      +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*SIND(DG)*COSD(DV)*OBSRADAR
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
      +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*COSD(DG)*COSD(DV)*OBSRADAR
      GRADINT(I,J,K,T,WW)=GRADINT(I,J,K,T,WW)  &
      +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*SIND(DV)*OBSRADAR
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  S=NUMSTAT+2                !  FOR SFMR OBSERVATION
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)*OBSERROR(O)
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    HU0=0.0
    HV0=0.0
    HW0=0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,WW)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
      HW0=HW0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,WW)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0)+(HW+HW0)*(HW+HW0))
    HT0=SQRT(HU0*HU0+HV0*HV0+HW0*HW0)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
      +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HU+HU0)/HT*OBS_SFMR
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
      +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HV+HV0)/HT*OBS_SFMR
      GRADINT(I,J,K,T,WW)=GRADINT(I,J,K,T,WW)  &
      +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HW+HW0)/HT*OBS_SFMR
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

! SMOOTH TERM
  CALL SMOTHGRAD

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNGRAD_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                     ! FOR PRESSURE COORDINATE
!    CALL GSBLNGRAD_P_HS
      CALL CNTNSGRAD
    ENDIF
  ENDIF

! HYDROSTATIC CONDITION TERM                 ! ADDED BY ZHONGJIE HE
  IF(GRDLEVL .LE. ENDHYLV) THEN
    CALL HYDROGRAD
  ENDIF


  RETURN
END SUBROUTINE COSTGRADT1


SUBROUTINE COSTFUNCT2
!*************************************************
! CALCULATE COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
! HISTORY: MODIFIED BY ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV
  INTEGER  :: NP(MAXDIMS)
  REAL  :: HT,OI,HU,HV,HW,DG,DV
  REAL  :: HT0,HU0,HV0
  REAL  :: CS(NUMSTAT+3),CM,CB,CH,rdr,ref
! added by shuyuan 20100907  for calculate ref, qr,qs..
  REAL  ::rc     !   desity*rain water mixing ration    
  REAL  ::sc     !   desity*snow water mixing ration 
  REAL  :: Segma     !variance  for ref_obs
  REAL  ::Temp_ref,tempr,temps
  INTEGER  :: RR,RS  

  UU=U_CMPNNT
  VV=V_CMPNNT
!added by shuyuan 20100728
  RR=ROUR_CMPNNT
  RS=ROUS_CMPNNT

!  CALL WCOMPGERNL        !  MODIFIED BY ZHONGJIE HE.
! INITIALIZATION
  COSTFUN=0.0
  DO S=1,NUMSTAT+3
    CS(S)=0.0
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      COSTFUN=COSTFUN+(HT-OBSVALUE(O))*(HT-OBSVALUE(O))!/(NOBSTAT(S)*1.0)!!!modified by shuyuan   20101028
      CS(S)=CS(S)+0.5*(HT-OBSVALUE(O))*(HT-OBSVALUE(O))!/(NOBSTAT(S)*1.0)!!!modified by shuyuan   20101028
    ENDDO
  ENDDO
  S=NUMSTAT+1                        ! FOR RADIAL WIND VELOCITY OBSERVATIONS
rdr = COSTFUN
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*WWW(I,J,K,T)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    DV=OBSEINF2(NO)
    HT=HU*SIND(DG)*COSD(DV)+HV*COSD(DG)*COSD(DV)+HW*SIND(DV)
    COSTFUN=COSTFUN+(HT-OBSVALUE(O))*(HT-OBSVALUE(O))!*OBSRADAR!!!modified by shuyuan   20101028
  ENDDO

print*,'Radar cost: ',COSTFUN-rdr
  S=NUMSTAT+2                        ! FOR SFMR WIND VELOCITY OBSERVATIONS
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI 
    HU=0.0
    HV=0.0
    HW=0.0
    HU0=0.0
    HV0=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*WWW(I,J,K,T)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0)+HW*HW)
    HT0=SQRT(HU0*HU0+HV0*HV0)
    COSTFUN=COSTFUN+(HT-HT0-OBSVALUE(O))*(HT-HT0-OBSVALUE(O))!*OBS_SFMR!!!modified by shuyuan   20101028
  ENDDO


ref = COSTFUN
! --------------------
! added by shuyuan 20100721    
  Segma=1.
   S=NUMSTAT+3
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      Temp_ref=0.0
      M=0
      rc=0.
      sc=0.
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
            
        M=M+1        
        rc=rc+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,RR)
      ENDDO
      ENDDO
      ENDDO
      ENDDO     
       !changed 20100907  shuyuan
  
     ! Temp_ref=(OBSVALUE(O)-43.1)/17.5
     ! Temp_ref=(10.**Temp_ref)              !rc  unit is g/m3 
      Temp_ref=OBSVALUE(O)
      Temp_ref=(rc-Temp_ref)*(rc-Temp_ref)

      !COSTFUN=COSTFUN+Temp_ref*OI/(Segma*Segma)/(NOBSTAT(S)*1.0) !!!!20101025
     COSTFUN=COSTFUN+Temp_ref!!!!20101025

     CS(S)=CS(S)+0.5*(rc-Temp_ref)*(rc-Temp_ref)/(Segma*Segma)
   
      ENDDO
  print*,'Radar cost: ',COSTFUN, OBSVALUE(O)



  COSTFUN=0.5D0*COSTFUN

! SMOOTH TERM
  CM=COSTFUN
  CALL SMOTHCOST
  CM=COSTFUN-CM

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    CB=COSTFUN
    IF(IFPCDNT.EQ. 0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNCOST_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                     ! FOR PRESSURE COORDINATE
      CALL GSBLNCOST_P_HS
    ENDIF
    CB=COSTFUN-CB
  ENDIF

! HYDROSTATIC CONDITION TERM                 ! ADDED BY ZHONGJIE HE
  IF(GRDLEVL .LE. ENDHYLV) THEN
    CH=COSTFUN
    ! CALL HYDROCOST
    !CALL HYDROCOST_XIE
    call HYDROCOST_SHUYUAN
    CH=COSTFUN-CH
  ENDIF

   OPEN(907,file='costfun.txt',position='append')!!2010823 liu for test
   WRITE(907,1)(CS(S),S=1,NUMSTAT+3),CM,CB,CH
   CLOSE(907)
  WRITE(*,1)(CS(S),S=1,NUMSTAT+3),CM,CB,CH
1 FORMAT('Each state: ',10e12.4,' Smoothing: ',e12.4,' Geostropic: ',e12.4,' Hydrostatic: ',e12.4)
  

  RETURN
END SUBROUTINE COSTFUNCT2

SUBROUTINE COSTGRADT2
!*************************************************
! CALCULATE GRADIENT OF COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
! HISTORY: JANUARY 2008, MODIFIED BY ZHONGJIE HE
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,I1,I2,J1,J2,K1,S,O,NO,M,N,UU,VV,ZZ,RR,RS
  INTEGER  :: NP(MAXDIMS)
  REAL     :: HT,OI,HU,HV,HW,DG,DV,Z1,Z2
  REAL     :: HT0,HU0,HV0
  REAL     :: CC(NGPTOBS,NALLOBS)
!-------------------------------------------
! added by shuyuan 20100713  for calculate ref, qr,qs..
  REAL  ::RouR ,rc     !   desity*rain water mixing ration,rc=ROUR  
  REAL  ::RouS  ,sc    !   desity*snow water mixing ration,sc=ROUS
  REAL  :: C1,C2    ! coeff  for rour and rous
  REAL  :: Segma     !variance  for ref_obs
  REAL  ::Temp_ref,temp,temp1
   

!  DECLARATION :
!                'CC' IS THE COEFFICENT USED TO CALCULATE GRADIENT OF W TO CONTROL VARIABLE
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT
  ZZ=PRESSURE
!added by shuyuan 20100728
  RR=ROUR_CMPNNT
  RS=ROUS_CMPNNT
!  CALL WCOMPGERNL        ! MODIFIED BY ZHONGJIE HE
! INITIALIZATION

  DO S=1,NUMSTAT+1
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      GRADINT(I,J,K,T,S)=0.0
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)*OBSERROR(O)
      OI=1.0/OI
      HT=0.0
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        GRADINT(I,J,K,T,S)=GRADINT(I,J,K,T,S)  &
        +(HT-OBSVALUE(O))*OBSCOEFF(M,O)!!!modified by shuyuan   20101028
        ! +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)/(NOBSTAT(S)*1.0)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDDO
  S=NUMSTAT+1                           ! FOR RADIAL WIND OBSERVATIONS
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)*OBSERROR(O)
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*WWW(I,J,K,T)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
     
    DG=OBSEINF1(NO)
    DV=OBSEINF2(NO)
    HT=HU*SIND(DG)*COSD(DV)+HV*COSD(DG)*COSD(DV)+HW*SIND(DV)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
      +(HT-OBSVALUE(O))*OBSCOEFF(M,O)*SIND(DG)*COSD(DV)
     ! +(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*SIND(DG)*COSD(DV)*OBSRADAR!!!modified by shuyuan   20101028
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
       +(HT-OBSVALUE(O))*OBSCOEFF(M,O)*COSD(DG)*COSD(DV)
      !+(HT-OBSVALUE(O))*OI*OBSCOEFF(M,O)*COSD(DG)*COSD(DV)*OBSRADAR!!!modified by shuyuan   20101028
      CC(M,O)=(HT-OBSVALUE(O))*SIND(DV)*OBSCOEFF(M,O)!*OI*OBSRADAR!!!modified by shuyuan   20101028
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

  CALL WWGRADIENT(CC,S)  ! ADDED BY ZHONGJIE HE TO CALCULATE THE GRADIENT OF WWW TO U, V AND Z

  S=NUMSTAT+2                !  FOR SFMR OBSERVATION
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)*OBSERROR(O)
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    HU0=0.0
    HV0=0.0
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*WWW(I,J,K,T)
      HU0=HU0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,UU)
      HV0=HV0+OBSCOEFF(M,O)*GRDBKGND(I,J,K,T,VV)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    HT=SQRT((HU+HU0)*(HU+HU0)+(HV+HV0)*(HV+HV0)+HW*HW)
    HT0=SQRT(HU0*HU0+HV0*HV0)
    M=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      GRADINT(I,J,K,T,UU)=GRADINT(I,J,K,T,UU)  &
       +(HT-HT0-OBSVALUE(O))*OBSCOEFF(M,O)*(HU+HU0)/HT
     ! +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HU+HU0)/HT*OBS_SFMR!!!modified by shuyuan   20101028
      GRADINT(I,J,K,T,VV)=GRADINT(I,J,K,T,VV)  &
      +(HT-HT0-OBSVALUE(O))*OBSCOEFF(M,O)*(HV+HV0)/HT
  !    +(HT-HT0-OBSVALUE(O))*OI*OBSCOEFF(M,O)*(HV+HV0)/HT*OBS_SFMR!!!modified by shuyuan   20101028
      CC(M,O)=(HT-HT0-OBSVALUE(O))*OBSCOEFF(M,O)*HW/HT!*OI*OBS_SFMR!!!modified by shuyuan   20101028
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

  CALL WWGRADIENT(CC,S)
! for radar  reflectivity !  shuyuan 20100721
! --20100907------------------
   
  Segma=1.
   S=NUMSTAT+3
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)*OBSERROR(O)
      OI=1.0/OI
      Temp_ref=0.0
      M=0
      rc=0.
      sc=0.
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
       
        M=M+1
        rc=rc+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,RR)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      M=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1     
      !changed 20100907 shuyuan  rc unit is  g/m3   
      !Temp_ref=(OBSVALUE(O)-43.1)/17.5
      !Temp_ref=10.**Temp_ref
      Temp_ref=OBSVALUE(O)
      temp=(rc-Temp_ref)/(Segma*Segma)
      temp1=temp*OBSCOEFF(M,O)!!*OI/(NOBSTAT(S)*1.0)  !!!modified by shuyuan   20101028
      GRADINT(I,J,K,T,RR)=GRADINT(I,J,K,T,RR)+temp1  
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO


! SMOOTH TERM
  CALL SMOTHGRAD

! GEOSTROPHIC BALANCE TERM
  IF(GRDLEVL .LE. ENDGSLV) THEN
    IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN     ! FOR SIGMA AND HEIGHT COORDINATE
      CALL GSBLNGRAD_NON_UNIFORM_Z
    ELSEIF(IFPCDNT.EQ.1)THEN                   ! FOR PRESSURE COORDINATE
      CALL GSBLNGRAD_P_HS
    ENDIF
  ENDIF

! HYDROSTATIC CONDITION TERM                 ! ADDED BY ZHONGJIE HE
  IF(GRDLEVL .LE. ENDHYLV) THEN
    ! CALL HYDROGRAD
    !CALL HYDROGRAD_XIE
     call HYDROGRAD_SHUYUAN
  ENDIF

  RETURN
END SUBROUTINE COSTGRADT2

SUBROUTINE TEST_FUN2
!*************************************************
! SHOW THE DIFFERENCE BETWEEN THE ANALYSIS AND THE OBSERVATIONS
! HISTORY: JANUARY 2008, CODED by ZHONGJIE HE.

!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S,O,NO,M,N,UU,VV
  INTEGER  :: NP(MAXDIMS)
  REAL  :: HT,OI,HU,HV,HW,DG,DV
  REAL  :: CS(NUMSTAT),CM,CB
  
  REAL  :: OBS_X,OBS_Y,OBS_Z
! --------------------
  UU=U_CMPNNT
  VV=V_CMPNNT

  CALL WCOMPGERNL        !  MODIFIED BY ZHONGJIE HE.

! INITIALIZATION
  COSTFUN=0.0
  DO S=1,NUMSTAT
    CS(S)=0.0
  ENDDO
  IF(NALLOBS.EQ.0)RETURN
  O=0
  DO S=1,NUMSTAT
    DO NO=1,NOBSTAT(S)
      O=O+1
      DO N=1,MAXDIMS
        NP(N)=OBSIDXPC(N,O)
      ENDDO
      OI=OBSERROR(O)**2
      OI=1.0/OI
      HT=0.0
      M=0
      OBS_X=0
      OBS_Y=0
      OBS_Z=0
      DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
      DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
      DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
      DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
        M=M+1
        HT=HT+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,S)
        OBS_X=OBS_X+OBSCOEFF(M,O)*I
        OBS_Y=OBS_Y+OBSCOEFF(M,O)*J
        OBS_Z=OBS_Z+OBSCOEFF(M,O)*K
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))
      CS(S)=CS(S)+0.5*(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))
    ENDDO
  ENDDO
  S=NUMSTAT+1
  DO NO=1,NOBSTAT(S)
    O=O+1
    DO N=1,MAXDIMS
      NP(N)=OBSIDXPC(N,O)
    ENDDO
    OI=OBSERROR(O)**2
    OI=1.0/OI
    HU=0.0
    HV=0.0
    HW=0.0
    M=0
    OBS_X=0
    OBS_Y=0
    OBS_Z=0
    DO T=NP(4),MIN0(NP(4)+1,NUMGRID(4))
    DO K=NP(3),MIN0(NP(3)+1,NUMGRID(3))
    DO J=NP(2),MIN0(NP(2)+1,NUMGRID(2))
    DO I=NP(1),MIN0(NP(1)+1,NUMGRID(1))
      M=M+1
      HU=HU+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,UU)
      HV=HV+OBSCOEFF(M,O)*GRDANALS(I,J,K,T,VV)
      HW=HW+OBSCOEFF(M,O)*WWW(I,J,K,T)
      OBS_X=OBS_X+OBSCOEFF(M,O)*I
      OBS_Y=OBS_Y+OBSCOEFF(M,O)*J
      OBS_Z=OBS_Z+OBSCOEFF(M,O)*K
    ENDDO
    ENDDO
    ENDDO
    ENDDO
    DG=OBSEINF1(NO)
    DV=OBSEINF2(NO)
    HT=HU*SIND(DG)*COSD(DV)+HV*COSD(DG)*COSD(DV)+HW*SIND(DV)
    COSTFUN=COSTFUN+(HT-OBSVALUE(O))*OI*(HT-OBSVALUE(O))*OBSRADAR
  ENDDO
  COSTFUN=0.5D0*COSTFUN

! SMOOTH TERM
  CM=COSTFUN
  CALL SMOTHCOST
  CM=COSTFUN-CM

! GEOSTROPHIC BALANCE TERM
  CB=COSTFUN
  IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN      ! FOR SIGMA AND HEIGHT COORDINATE
    CALL GSBLNCOST_NON_UNIFORM_Z
  ELSEIF(IFPCDNT.EQ.1)THEN                    ! FOR PRESSURE COORDINATE
    CALL GSBLNCOST_P_HS
  ENDIF
  CB=COSTFUN-CB
!  WRITE(100,*)(CS(S),S=1,NUMSTAT),CM,CB

  RETURN
END SUBROUTINE TEST_FUN2

SUBROUTINE CG_VALUE(F,X,NM)
!*************************************************
! CALCULATE THE COSTFUNCTION, USED BY MINIMIZER_CG
! HISTORY: APRIL 2008, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER    :: I,J,K,T,S,N,NM
!  DOUBLE PRECISION :: X(NM),F
  REAL :: X(NM),F
  N=0
  DO S=1,NUMSTAT
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      N=N+1
      GRDANALS(I,J,K,T,S)=X(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

  IF(W_CMPNNT.NE.0) THEN
    CALL COSTFUNCT1
  ELSE
    CALL WCOMPGERNL
    CALL COSTFUNCT2
  ENDIF

  F=COSTFUN

END SUBROUTINE CG_VALUE

SUBROUTINE CG_GRAD(G,X,NM)
!*************************************************
! CALCULATE THE GRADIENT, USED BY MINIMIZER_CG
! HISTORY: APRIL 2008, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER    :: I,J,K,T,S,N,NM
!  DOUBLE PRECISION :: X(NM),G(NM)
  REAL :: X(NM),G(NM)

  N=0
  DO S=1,NUMSTAT
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      N=N+1
      GRDANALS(I,J,K,T,S)=X(N)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO
  IF(W_CMPNNT.NE.0) THEN
    CALL COSTGRADT1
  ELSE
    CALL WCOMPGERNL
    CALL COSTGRADT2
  ENDIF
  N=0
  DO S=1,NUMSTAT
    DO T=1,NUMGRID(4)
    DO K=1,NUMGRID(3)
    DO J=1,NUMGRID(2)
    DO I=1,NUMGRID(1)
      N=N+1
      G(N)=GRADINT(I,J,K,T,S)
    ENDDO
    ENDDO
    ENDDO
    ENDDO
  ENDDO

END SUBROUTINE CG_GRAD

END MODULE COSTFUN_GRAD
