MODULE SMCOSTF_GRAD
! CALCULATE THE PENALTY COST FUNCTION OF SMOOTH AND THE RELERTIVE GRADIENT
! HISTORY: JANUARY 2008, SEPARATED FROME THE MODULE 'STMAS4D_CORE' BY ZHONGJIE HE

  USE PRMTRS_STMAS
  USE GENERALTOOLS, ONLY : G2ORDERIT

  PUBLIC       SMOTHCOST, SMOTHGRAD

!**************************************************
!COMMENT:
!   THIS MODULE IS USED BY THE MODULE OF COSTFUN_GRAD TO CALCULATE THE PENALTY COST FUNCTION OF SMOOTH AND THE RELERTIVE GRADIENT
!   SUBROUTINES:
!      SMOTHCOST: CALCULATE SMOOTH PENALTY TERM OF COST FUNCTION.
!      SMOTHGRAD: CALCULATE GRADIENTS OF SMOOTH TERM.
!**************************************************
CONTAINS

SUBROUTINE SMOTHCOST
!*************************************************
! SMOOTH PENALTY TERM OF COST FUNCTION
! HISTORY: AUGUST 2007, CODED by WEI LI.
!
!          MAY 2015 MODIFIED by YUANFU XIE
!          add a penalty parameter check
!          and laplacian check: skip if small enough
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S
  REAL     :: Z1,Z2,Z3,AZ,BZ,CZ,sm
  double precision :: smx,smy,smz,smt
! --------------------
! PENALTY TERM

  PRINT*,'SMOOTHING X: ',PENAL_X(1:NUMSTAT)
  PRINT*,'SMOOTHING Y: ',PENAL_Y(1:NUMSTAT)
  PRINT*,'SMOOTHING Z: ',PENAL_Z(1:NUMSTAT)
  PRINT*,'SMOOTHING T: ',PENAL_T(1:NUMSTAT)
  smx = 0.0d0
  IF(NUMGRID(1).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_X(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO K=1,NUMGRID(3)
      DO J=1,NUMGRID(2)
      DO I=2,NUMGRID(1)-1
        ! COSTFUN=COSTFUN+PENAL_X(S)* &
        sm = &
          ((GRDANALS(I-1,J,K,T,S)-2*GRDANALS(I,J,K,T,S)+GRDANALS(I+1,J,K,T,S))/ &
          (GRDSPAC(1)*GRDSPAC(1)))**2
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_X(S)*0.1) smx=smx+PENAL_X(S)*sm
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF

  smy = 0.0d0
  IF(NUMGRID(2).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_Y(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO K=1,NUMGRID(3)
      DO J=2,NUMGRID(2)-1
      DO I=1,NUMGRID(1)
        ! COSTFUN=COSTFUN+PENAL_Y(S)* &
        sm = &
          ((GRDANALS(I,J-1,K,T,S)-2*GRDANALS(I,J,K,T,S)+GRDANALS(I,J+1,K,T,S))/ &
          (GRDSPAC(2)*GRDSPAC(2)))**2
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_Y(S)*0.1) smy=smy+PENAL_Y(S)*sm
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF

  smz = 0.0d0
  IF(NUMGRID(3).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_Z(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO K=2,NUMGRID(3)-1
      DO J=1,NUMGRID(2)
      DO I=1,NUMGRID(1)
        IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
          Z1=ZZZ(I,J,K-1,T)
          Z2=ZZZ(I,J,K  ,T)
          Z3=ZZZ(I,J,K+1,T)
        ELSEIF(IFPCDNT.EQ.1)THEN                     ! FOR PRESSURE COORDINATE
          Z1=PPP(K-1)
          Z2=PPP(K  )
          Z3=PPP(K+1)
        ENDIF

        CALL G2ORDERIT(Z1,Z2,Z3,AZ,BZ,CZ)

        AZ=AZ*(Z2-Z1)*(Z3-Z2)
        BZ=BZ*(Z2-Z1)*(Z3-Z2)
        CZ=CZ*(Z2-Z1)*(Z3-Z2)
        ! COSTFUN=COSTFUN+PENAL_Z(S)* &
        sm = &
          (AZ*GRDANALS(I,J,K-1,T,S)+BZ*GRDANALS(I,J,K,T,S)+CZ*GRDANALS(I,J,K+1,T,S))**2
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_Z(S)*0.1) smz=smz+PENAL_Z(S)*sm
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF

  smt = 0.0d0
  IF(NUMGRID(4).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_T(S) .GT. 0.0) THEN
      DO T=2,NUMGRID(4)-1
      DO K=1,NUMGRID(3)
      DO J=1,NUMGRID(2)
      DO I=1,NUMGRID(1)
        ! COSTFUN=COSTFUN+PENAL_T(S)* &
        sm = &
          ((GRDANALS(I,J,K,T-1,S)-2*GRDANALS(I,J,K,T,S)+GRDANALS(I,J,K,T+1,S))/ &
          (GRDSPAC(4)*GRDSPAC(4)))**2
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_T(S)*0.1) smt=smt+PENAL_T(S)*sm
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF

  COSTFUN = COSTFUN+smx+smy+smz+smt

  RETURN
END SUBROUTINE SMOTHCOST

SUBROUTINE SMOTHGRAD
!*************************************************
! SMOOTH PENALTY TERM OF GRADIENT
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,T,S
  REAL     :: Z1,Z2,Z3,AZ,BZ,CZ
  REAL     :: SM
! --------------------
! PENALTY TERM
  IF(NUMGRID(1).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_X(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO K=1,NUMGRID(3)
      DO J=1,NUMGRID(2)
      DO I=2,NUMGRID(1)-1
        SM=(GRDANALS(I-1,J,K,T,S)-2.0*GRDANALS(I,J,K,T,S)+GRDANALS(I+1,J,K,T,S))/ &
          (GRDSPAC(1)*GRDSPAC(1))
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_X(S)*0.1) THEN
          GRADINT(I-1,J,K,T,S)=GRADINT(I-1,J,K,T,S)+2.0*PENAL_X(S)/(GRDSPAC(1)*GRDSPAC(1))*SM
          GRADINT(I  ,J,K,T,S)=GRADINT(I  ,J,K,T,S)+2.0*PENAL_X(S)/(GRDSPAC(1)*GRDSPAC(1))*(-2.0)*SM
          GRADINT(I+1,J,K,T,S)=GRADINT(I+1,J,K,T,S)+2.0*PENAL_X(S)/(GRDSPAC(1)*GRDSPAC(1))*SM
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF(NUMGRID(2).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_Y(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO K=1,NUMGRID(3)
      DO I=1,NUMGRID(1)
      DO J=2,NUMGRID(2)-1
        SM=(GRDANALS(I,J-1,K,T,S)-2.0*GRDANALS(I,J,K,T,S)+GRDANALS(I,J+1,K,T,S))/ &
          (GRDSPAC(2)*GRDSPAC(2))
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_Y(S)*0.1) THEN
          GRADINT(I,J-1,K,T,S)=GRADINT(I,J-1,K,T,S)+2.0*PENAL_Y(S)/(GRDSPAC(2)*GRDSPAC(2))*SM
          GRADINT(I,J  ,K,T,S)=GRADINT(I,J  ,K,T,S)+2.0*PENAL_Y(S)/(GRDSPAC(2)*GRDSPAC(2))*(-2.0)*SM
          GRADINT(I,J+1,K,T,S)=GRADINT(I,J+1,K,T,S)+2.0*PENAL_Y(S)/(GRDSPAC(2)*GRDSPAC(2))*SM
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF(NUMGRID(3).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_Z(S) .GT. 0.0) THEN
      DO T=1,NUMGRID(4)
      DO J=1,NUMGRID(2)
      DO I=1,NUMGRID(1)
      DO K=2,NUMGRID(3)-1
        IF(IFPCDNT.EQ.0 .OR. IFPCDNT.EQ.2)THEN       ! FOR SIGMA AND HEIGHT COORDINATE
          Z1=ZZZ(I,J,K-1,T)
          Z2=ZZZ(I,J,K  ,T)
          Z3=ZZZ(I,J,K+1,T)
        ELSEIF(IFPCDNT.EQ.1)THEN                     ! FOR PRESSUR COORDINATE
          Z1=PPP(K-1)
          Z2=PPP(K  )
          Z3=PPP(K+1)
        ENDIF
        CALL G2ORDERIT(Z1,Z2,Z3,AZ,BZ,CZ)
        AZ=AZ*(Z2-Z1)*(Z3-Z2)
        BZ=BZ*(Z2-Z1)*(Z3-Z2)
        CZ=CZ*(Z2-Z1)*(Z3-Z2)
        SM=AZ*GRDANALS(I,J,K-1,T,S)+BZ*GRDANALS(I,J,K,T,S)+CZ*GRDANALS(I,J,K+1,T,S)
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_Z(S)*0.1) THEN
          GRADINT(I,J,K-1,T,S)=GRADINT(I,J,K-1,T,S)+2.0*PENAL_Z(S)*AZ*SM
          GRADINT(I,J,K  ,T,S)=GRADINT(I,J,K  ,T,S)+2.0*PENAL_Z(S)*BZ*SM
          GRADINT(I,J,K+1,T,S)=GRADINT(I,J,K+1,T,S)+2.0*PENAL_Z(S)*CZ*SM
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF
  IF(NUMGRID(4).GE.3)THEN
    DO S=1,NUMSTAT
      IF (PENAL_T(S) .GT. 0.0) THEN
      DO K=1,NUMGRID(3)
      DO J=1,NUMGRID(2)
      DO I=1,NUMGRID(1)
      DO T=2,NUMGRID(4)-1
        SM=(GRDANALS(I,J,K,T-1,S)-2.0*GRDANALS(I,J,K,T,S)+GRDANALS(I,J,K,T+1,S))/ &
          (GRDSPAC(4)*GRDSPAC(4))
        ! YUANFU added a check to skip to small Laplacian: May 2015
        IF (sm .GT. PENAL_T(S)*0.1) THEN
          GRADINT(I,J,K,T-1,S)=GRADINT(I,J,K,T-1,S)+2.0*PENAL_T(S)/(GRDSPAC(4)*GRDSPAC(4))*SM
          GRADINT(I,J,K,T  ,S)=GRADINT(I,J,K,T  ,S)+2.0*PENAL_T(S)/(GRDSPAC(4)*GRDSPAC(4))*(-2.0)*SM
          GRADINT(I,J,K,T+1,S)=GRADINT(I,J,K,T+1,S)+2.0*PENAL_T(S)/(GRDSPAC(4)*GRDSPAC(4))*SM
        ENDIF
      ENDDO
      ENDDO
      ENDDO
      ENDDO
      ENDIF
    ENDDO
  ENDIF
  RETURN
END SUBROUTINE SMOTHGRAD

END MODULE SMCOSTF_GRAD
