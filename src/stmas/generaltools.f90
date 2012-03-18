MODULE GENERALTOOLS

CONTAINS

SUBROUTINE INTERPLTN(ND,NG,CO,AC,OC)
!*************************************************
! CALCULATE THE INTERPOLATING COEFFICENT (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
! 2     , 3     , 5     , 9
! 2**0+1, 2**1+1, 2**2+1, 2**3+1
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J
  INTEGER , INTENT(IN)  :: ND,NG
  REAL    , INTENT(IN)  :: AC(ND,NG),OC(ND)
  REAL    , INTENT(OUT) :: CO(NG)
  REAL                  :: DD(ND,NG),DC(ND)
! --------------------
! DECLARE:
!        'ND' IS THE NUMBER OF DIMENSIONS OF THE FIELD
!        'NG' IS THE NUMBER OF GRIDPOINTS NEIGHBOR TO THE POSITION, WHICH IS TO BE INTERPOLATED, NG IS EQUAL TO 2**ND
!        'CO' IS THE INTERPOLATING COEFFICIENTS.
!        'OC' IS THE POSITION WHERE TO BE INTERPOLATED 
!        'AC' IS THE INDEX OF THE GRIDPOINTS, WHICH IS NEIGHBOR TO THE POSITIOIN
! --------------------
  DO I=1,ND
    IF((OC(I)-AC(I,1))*(OC(I)-AC(I,2**(I-1)+1)).GT.0.0) THEN
      PRINT*,'INTERPLTN: OC IS NOT IN THE DOMAIN ',I,OC(I),AC(I,1),AC(I,2**(I-1)+1)
      STOP
    ENDIF
  ENDDO
  DO I=1,ND
    DO J=1,NG
      DD(I,J)=AC(I,J)-OC(I)
    ENDDO
    DC(I)=ABS(AC(I,2**(I-1)+1)-AC(I,1))
!    IF(DC(I).LE.0.0)STOP 'GRID ARRANGE WRONG!!!'            ! omitted by zhongjie he
  ENDDO
  DO I=1,ND
    IF(DC(I).NE.0) THEN 
      DO J=1,NG
        DD(I,J)=DD(I,J)/DC(I)
        IF(DD(I,J).GT.0.0)THEN
          DD(I,J)=DD(I,J)-1.0
        ELSE
          DD(I,J)=DD(I,J)+1.0
        ENDIF
        DD(I,J)=ABS(DD(I,J))
      ENDDO
    ELSE                                     !     added by zhongjie he 
      DO J=1,NG
        DD(I,J)=1.
      ENDDO
    ENDIF
  ENDDO
  DO J=1,NG
    CO(J)=1.0
    DO I=1,ND
      CO(J)=CO(J)*DD(I,J)
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE INTERPLTN

SUBROUTINE INTERPLTN_XIE(ND,NG,CO,AC,OC,IP,NP,PV)
!*************************************************
! CALCULATE THE INTERPOLATING COEFFICENT (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!
!          MARCH 31, 2009 BY YUANFU XIE FOR LOG(P)
!          INTERPOLATION.
!*************************************************
! 2     , 3     , 5     , 9
! 2**0+1, 2**1+1, 2**2+1, 2**3+1
  IMPLICIT NONE
! --------------------
  INTEGER  :: I,J,K,K1
  INTEGER , INTENT(IN)  :: ND,NG,IP,NP
  REAL    , INTENT(IN)  :: AC(ND,NG),OC(ND),PV(NP)
  REAL    , INTENT(OUT) :: CO(NG)
  REAL                  :: DD(ND,NG),DC(ND),P
! --------------------
! DECLARE:
!        'ND' IS THE NUMBER OF DIMENSIONS OF THE FIELD
!        'NG' IS THE NUMBER OF GRIDPOINTS NEIGHBOR TO THE POSITION, WHICH IS TO BE INTERPOLATED, NG IS EQUAL TO 2**ND
!        'CO' IS THE INTERPOLATING COEFFICIENTS.
!        'OC' IS THE POSITION WHERE TO BE INTERPOLATED 
!        'AC' IS THE INDEX OF THE GRIDPOINTS, WHICH IS NEIGHBOR TO THE POSITIOIN
! --------------------
  DO I=1,ND
    IF((OC(I)-AC(I,1))*(OC(I)-AC(I,2**(I-1)+1)).GT.0.0) THEN
      PRINT*,'INTERPLTN: OC IS NOT IN THE DOMAIN ',I,OC(I),AC(I,1),AC(I,2**(I-1)+1)
      STOP
    ENDIF
  ENDDO
  DO I=1,ND
    DO J=1,NG
      DD(I,J)=AC(I,J)-OC(I)
    ENDDO
    DC(I)=ABS(AC(I,2**(I-1)+1)-AC(I,1))
  ENDDO
  DO I=1,ND
    IF(DC(I).NE.0) THEN 
      DO J=1,NG
        DD(I,J)=DD(I,J)/DC(I)
        IF(DD(I,J).GT.0.0)THEN
          ! LOG(P):
          IF (I .EQ. IP) THEN
            K1 = AC(I,J)
            IF (K1 .LE. 1) THEN
              DD(I,J) = 0.0
            ELSE
              K = K1-1
              P = PV(K)*DD(I,J)+(1.0-DD(I,J))*PV(K+1)
              DD(I,J) = (LOG(P)-LOG(PV(k)))/(LOG(PV(K1))-LOG(PV(K)))
            ENDIF
          ELSE
            DD(I,J)=DD(I,J)-1.0
          ENDIF
        ELSE
          IF (I .EQ. IP) THEN
            K = AC(I,J)
            IF (K .GE. NP) THEN
              DD(I,J) = 1.0
            ELSE
              K1 = K+1
              P = (1.0+DD(I,J))*PV(K)-DD(I,J)*PV(K1)
              DD(I,J) = (LOG(PV(K1))-LOG(P))/(LOG(PV(K1))-LOG(PV(K)))
            ENDIF
          ELSE
            DD(I,J)=DD(I,J)+1.0
          ENDIF
        ENDIF
        DD(I,J)=ABS(DD(I,J))
      ENDDO
    ELSE                                     !     added by zhongjie he 
      DO J=1,NG
        DD(I,J)=1.
      ENDDO
    ENDIF
  ENDDO

  DO J=1,NG
    CO(J)=1.0
    DO I=1,ND
      CO(J)=CO(J)*DD(I,J)
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE INTERPLTN_XIE

SUBROUTINE GDERIVEIT(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL DERIVATIVE OF INTERIOR (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL  :: X,Y
  REAL , INTENT(IN)  :: Z1,Z2,Z3
  REAL , INTENT(OUT) :: A,B,C
! --------------------
  X=Z2-Z1
  Y=Z3-Z2
  A=-Y/(X*X+X*Y)
  B=(-X+Y)/(X*Y)
  C=X/(X*Y+Y*Y)
  RETURN
END SUBROUTINE GDERIVEIT

SUBROUTINE GDERIVELB(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL DERIVATIVE OF LEFT BOUNDARY (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL  :: X,Y
  REAL , INTENT(IN)  :: Z1,Z2,Z3
  REAL , INTENT(OUT) :: A,B,C
! --------------------
  X=Z2-Z1
  Y=Z3-Z2
  A=(-2.0*X-Y)/(X*X+X*Y)
  B=(X+Y)/(X*Y)
  C=-X/(X*Y+Y*Y)
  RETURN
END SUBROUTINE GDERIVELB

SUBROUTINE GDERIVERB(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL DERIVATIVE OF RIGHT BOUNDARY (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL  :: X,Y
  REAL , INTENT(IN)  :: Z1,Z2,Z3
  REAL , INTENT(OUT) :: A,B,C
! --------------------
  X=Z2-Z1
  Y=Z3-Z2
  A=Y/(X*X+X*Y)
  B=(-X-Y)/(X*Y)
  C=(X+2.0*Y)/(X*Y+Y*Y)
  RETURN
END SUBROUTINE GDERIVERB

SUBROUTINE G2ORDERIT(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL 2-ORDER DERIVATIVE OF INTERIOR (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL  :: X,Y
  REAL , INTENT(IN)  :: Z1,Z2,Z3
  REAL , INTENT(OUT) :: A,B,C
! --------------------
  X=Z2-Z1
  Y=Z3-Z2
  A=2.0/(X*X+X*Y)
  B=-2.0/(X*Y)
  C=2.0/(X*Y+Y*Y)
  RETURN
END SUBROUTINE G2ORDERIT

SUBROUTINE PSTN2NUMB(NN,NP,NM,NC)
!*************************************************
! TRANSFORM POSITION TO NUMBER (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER               :: PT,M,N
  INTEGER , INTENT(IN)  :: NN
  INTEGER , INTENT(IN)  :: NP(NN),NM(NN)
  INTEGER , INTENT(OUT) :: NC
! --------------------
! DECLARE:
!        'NN' IS THE NUMBER OF DIMENSIONS OF THE FIELD
!        'NP' IS THE INDEX OF THE GRIDPOINT AT EACH DIMENSION
!        'NM' IS THE TOTAL GRID NUMBER OF EACH DIMENSIOIN 
!        'NC' IS THE INDEX OF THE GRIDPOINT, WHEN THE FIELD IS TRANSLATED TO A 1 DIMENSIONAL FILED
! --------------------
  NC=NP(1)
  IF(NN.EQ.1)RETURN
  DO N=2,NN
    PT=NP(N)-1
    DO M=1,N-1
      PT=PT*NM(M)
    ENDDO
    NC=NC+PT
  ENDDO
  RETURN
END SUBROUTINE PSTN2NUMB

SUBROUTINE NUMB2PSTN(NN,NP,NM,NC)
!*************************************************
! TRANSFORM NUMBER TO POSITION (GENERAL) (NOT USED)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER               :: N,N0
  INTEGER , INTENT(IN)  :: NN,NC
  INTEGER , INTENT(IN)  :: NM(NN)
  INTEGER , INTENT(OUT) :: NP(NN)
! --------------------
! DECLARE:
!        'NN' IS THE NUMBER OF DIMENSIONS OF THE FIELD
!        'NP' IS THE INDEX OF THE GRIDPOINT AT EACH DIMENSION, WHEN THE FIELD IS TRANSLATED TO A NN DIMENSIONAL FILED
!        'NM' IS THE TOTAL GRID NUMBER OF EACH DIMENSIOIN
!        'NC' IS THE INDEX OF THE GRIDPOINT IN THE 1 DIMENSIONAL FIELD
! --------------------
  N0=NC
  NP(1)=MOD(N0,NM(1))
  IF(NP(1).EQ.0)NP(1)=NM(1)
  IF(NN.EQ.1)RETURN
  DO N=2,NN
    N0=(N0-NP(N-1))/NM(N-1)+1
    NP(N)=MOD(N0,NM(N))
    IF(NM(N).EQ.1)NP(N)=1
    IF(NP(N).EQ.0)NP(N)=NM(N)
  ENDDO
  RETURN
END SUBROUTINE NUMB2PSTN

SUBROUTINE VRTCLPSTN8(KM,LM,PP,P,PS,IS)
!*************************************************
! GET VERTICAL POSITION (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: K
  INTEGER ,INTENT(IN)  :: KM
  REAL ,   INTENT(IN)  :: PP(KM)
  REAL ,   INTENT(IN)  :: P,LM
  REAL ,   INTENT(OUT) :: PS
  INTEGER ,INTENT(OUT) :: IS
! --------------------
! DECLARE:
!        'KM' IS THE NUMBER OF LEVELS OF THE FIELD
!        'PP' IS THE PRESURE OF EACH LEVEL
!        'P'  IS THE PRESURE OF THE POSITION, WHICH WE WANT TO GET ITS INDEX IN THE VERTICAL COORDINATE
!        'PS' IS THE INDEX OF THE POINT IN THE VERTICAL COORDINATE
!        'IS' IS THE STATE OF THE RESULT OF THIS SUBROUTINE, 0 MEANS WRONG, 1 MEANS OK
!        'LM' IS A LIMITATION TO DECIDE WHETHER THE OBSERVATION IS AVAILABLE (JUST USED FOR THE CASE: KM=1)
! --------------------
  IS=0
  IF(KM.GE.2) THEN
    DO K=1,KM-1
      IF((P-PP(K))*(P-PP(K+1)).GT.0.0)CYCLE
      PS=K+ABS(P-PP(K))/ABS(PP(K+1)-PP(K))
      IS=1
    ENDDO
  ELSEIF(KM.EQ.1) THEN
    IF(ABS(P-PP(1)).LE.LM) THEN
      PS=1
      IS=1
    ENDIF
  ENDIF
  RETURN
END SUBROUTINE VRTCLPSTN8

SUBROUTINE VRTCLPSTN(KM,LM,PP,P,PS,IS)
!*************************************************
! GET VERTICAL POSITION (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: K
  INTEGER ,INTENT(IN)  :: KM
  REAL ,   INTENT(IN)  :: PP(KM),P,LM
  REAL ,   INTENT(OUT) :: PS
  INTEGER ,INTENT(OUT) :: IS
! --------------------
! DECLARE:
!        'KM' IS THE NUMBER OF LEVELS OF THE FIELD
!        'PP' IS THE PRESURE OF EACH LEVEL
!        'P'  IS THE PRESURE OF THE POSITION, WHICH WE WANT TO GET ITS INDEX IN THE VERTICAL COORDINATE
!        'PS' IS THE INDEX OF THE POINT IN THE VERTICAL COORDINATE
!        'IS' IS THE STATE OF THE RESULT OF THIS SUBROUTINE, 0 MEANS WRONG, 1 MEANS OK
!        'LM' IS A LIMITATION OF PREESURE TO DECIDE WHETHER THE OBSERVATION IS AVAILABLE (JUST USED FOR THE CASE: KM=1)
! --------------------
  IS=0
  IF(KM.GE.2) THEN
    DO K=1,KM-1
      IF((P-PP(K))*(P-PP(K+1)).GT.0.0)CYCLE
      PS=K+ABS(P-PP(K))/ABS(PP(K+1)-PP(K))
      IS=1
    ENDDO
  ELSEIF(KM.EQ.1) THEN
    IF(ABS(P-PP(1)).LE.LM) THEN
      PS=1
      IS=1
    ENDIF
  ENDIF

  RETURN
END SUBROUTINE VRTCLPSTN

SUBROUTINE ZPCONVERT(KM,ZZ,PP,Z,P,IS)
!*************************************************
! FIND PRESSURE FROM HEIGHT (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER              :: K
  INTEGER ,INTENT(IN)  :: KM
  REAL ,   INTENT(IN)  :: ZZ(KM),Z
  REAL ,   INTENT(IN)  :: PP(KM)
  REAL ,   INTENT(OUT) :: P
  INTEGER ,INTENT(OUT) :: IS
! --------------------
! DECLARE:
!        'KM' IS THE NUMBER OF LEVELS OF THE FIELD
!        'ZZ' IS THE HEIGHT OF EACH LEVEL
!        'PP' IS THE PRESURE OF EACH LEVEL
!        'Z'  IS THE HEIGHT OF POINT, WHICH WE WANT TO GET ITS PRESURE
!        'P'  IS THE PRESURE CALCULATED FROM THE ZZ, PP, AND Z
!        'IS' IS THE STATE OF THE RESULT OF THIS SUBROUTINE, 0 MEANS WRONG, 1 MEANS OK 
! --------------------
  IS=0
  DO K=1,KM-1
    IF((Z-ZZ(K))*(Z-ZZ(K+1)).GT.0.0)CYCLE
    P=PP(K)*ABS(ZZ(K+1)-Z)/ABS(ZZ(K+1)-ZZ(K))+PP(K+1)*ABS(Z-ZZ(K))/ABS(ZZ(K+1)-ZZ(K))
    IS=1
  ENDDO
  RETURN
END SUBROUTINE ZPCONVERT

SUBROUTINE FCST2BKGD(ND,NG,NS,NF,XF,YF,ZF,TF,NB,XB,YB,ZB,TB,FC,BK)
!*************************************************
! INTERPOLATE MODEL FORECAST TO BACKGROUND USED IN ANALYZING (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER , PARAMETER  :: MD=4
  INTEGER , INTENT(IN) :: ND,NG,NS
  INTEGER , INTENT(IN) :: NF(MD),NB(MD)
  REAL    , INTENT(IN) :: XF(NF(1)),YF(NF(2)),ZF(NF(3)),TF(NF(4))
  REAL    , INTENT(IN) :: XB(NB(1)),YB(NB(2)),ZB(NB(3)),TB(NB(4))
  REAL    , INTENT(IN) :: FC(NF(1),NF(2),NF(3),NF(4),NS)
  INTEGER              :: I,J,K,L,I0,J0,K0,L0,I1,J1,K1,L1,S,M,N,NN
  REAL                 :: PP(MD)
  REAL                 :: AC(ND,NG),OC(ND),CO(NG)
  REAL    , INTENT(OUT):: BK(NB(1),NB(2),NB(3),NB(4),NS)
  INTEGER              :: LM(ND)
! --------------------
! DECLARE:
!        'ND' IS THE NUMBER OF DIMENSIONS OF THE FIELD
!        'NG' IS THE NUMBER OF GRIDPOINTS NEIGHBOR TO THE POSITION, WHICH IS TO BE INTERPOLATED, NG IS EQUAL TO 2**ND
!        'NS' IS THE NUMBER OF VARIABLES THAT TO BE INTERPOLATED
!        'NF' THE GRID NUMBER OF THE ORIGINAL FIELD
!        'NB' THE GRID NUMBER OF THE TARGET FIELD
!        'XF', 'YF', 'ZF' AND 'TF' ARE THE SCALES OF THE ORIGINAL FILED IN THE X, Y, Z AND T COORDINATES RESPECTIVELY
!        'XB', 'YB', 'ZB' AND 'TB' ARE THE SCALES OF THE TARGET FIELD N THE X, Y, Z AND T COORDINATES
!        'FC' THE ORIGINAL FIELD
!        'BK' THE TARGET FIELD
!        'CO' IS THE INTERPOLATING COEFFICIENTS.
!        'OC' IS THE POSITION WHERE TO BE INTERPOLATED
!        'AC' IS THE INDEX OF THE GRIDPOINTS, WHICH IS NEIGHBOR TO THE POSITIOIN
! --------------------
  DO L=1,NB(4)
  DO K=1,NB(3)
  DO J=1,NB(2)
  DO I=1,NB(1)
    CALL CHECKPSTN(NF(1),XF,XB(I),I0)
    CALL CHECKPSTN(NF(2),YF,YB(J),J0)
    CALL CHECKPSTN(NF(3),ZF,ZB(K),K0)
    CALL CHECKPSTN(NF(4),TF,TB(L),L0)
    IF(I0.EQ.0.OR.J0.EQ.0.OR.K0.EQ.0.OR.L0.EQ.0)STOP 'BACKGOUND CAN NOT BE GENERATED'
    M=0
!===============================================
!    DO L1=L0,MIN0(L0+1,NF(4))
!    DO K1=K0,MIN0(K0+1,NF(3))
!    DO J1=J0,MIN0(J0+1,NF(2))
!    DO I1=I0,MIN0(I0+1,NF(1))
!      M=M+1
!      PP(1)=XF(I1)
!      PP(2)=YF(J1)
!      PP(3)=ZF(K1)
!      PP(4)=TF(L1)
!      NN=0
!      DO N=1,MD
!        IF(NF(N).GE.2)THEN
!          NN=NN+1
!          AC(NN,M)=PP(N)
!        ENDIF
!      ENDDO
!    ENDDO
!    ENDDO
!    ENDDO
!    ENDDO
!================ modified by zhongjie he =====
    LM(1)=I0+1
    LM(2)=J0+1
    LM(3)=K0+1
    LM(4)=L0+1
    DO N=1,ND
      IF(ND.LT.N) LM(N)=LM(N)-1
    ENDDO
    DO L1=L0,LM(4)
    DO K1=K0,LM(3)
    DO J1=J0,LM(2)
    DO I1=I0,LM(1)
      M=M+1
      PP(1)=XF(MIN0(I1,NF(1)))
      PP(2)=YF(MIN0(J1,NF(3)))
      PP(3)=ZF(MIN0(K1,NF(2)))
      PP(4)=TF(MIN0(L1,NF(1)))
      NN=0
      DO N=1,MD
        IF(NF(N).GE.2)THEN
          NN=NN+1
          AC(NN,M)=PP(N)
        ENDIF
      ENDDO
    ENDDO
    ENDDO
    ENDDO
    ENDDO
!===============================================
    PP(1)=XB(I)
    PP(2)=YB(J)
    PP(3)=ZB(K)
    PP(4)=TB(L)
    NN=0
    DO N=1,MD
      IF(NF(N).GE.2)THEN
        NN=NN+1
        OC(NN)=PP(N)
      ENDIF
    ENDDO
    CALL INTERPLTN(ND,NG,CO,AC,OC)
    DO S=1,NS
      M=0
      BK(I,J,K,L,S)=0.0
      DO L1=L0,MIN0(L0+1,NF(4))
      DO K1=K0,MIN0(K0+1,NF(3))
      DO J1=J0,MIN0(J0+1,NF(2))
      DO I1=I0,MIN0(I0+1,NF(1))
        M=M+1
        BK(I,J,K,L,S)=BK(I,J,K,L,S)+CO(M)*FC(I1,J1,K1,L1,S)
      ENDDO
      ENDDO
      ENDDO
      ENDDO
    ENDDO
  ENDDO
  ENDDO
  ENDDO
  ENDDO
  RETURN
END SUBROUTINE FCST2BKGD

SUBROUTINE CHECKPSTN(KM,ZZ,Z,KZ)
!*************************************************
! FIND THE POSITION (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: K
  INTEGER , INTENT(IN)  :: KM
  REAL    , INTENT(IN)  :: ZZ(KM),Z
  INTEGER , INTENT(OUT) :: KZ
! --------------------
  KZ=0
  IF(KM.EQ.1.AND.ZZ(1).EQ.Z)THEN
    KZ=1
    RETURN
  ENDIF
  DO K=1,KM-1
    IF((Z-ZZ(K))*(Z-ZZ(K+1)).GT.0.0)CYCLE
    KZ=K
    EXIT
  ENDDO
  RETURN
END SUBROUTINE CHECKPSTN

SUBROUTINE GRIDTRANS(NV,N1,X1,Y1,Z1,T1,N2,X2,Y2,Z2,T2,G1,G2)
!*************************************************
!  THIS ROUTINE MAPS GRID FUNCTIONS TO OTHER GRID.
!
!  HISTORY: APR. 2008, CODED BY YUANFU XIE.
!*************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NV	! NUMBER OF VAR TO MAP
  INTEGER, INTENT(IN) :: N1(4)	! NUMBERS OF INPUT GRIDPTS
  INTEGER, INTENT(IN) :: N2(4)	! NUMBERS OF OUTPUT GRIDPTS

  REAL,    INTENT(IN) :: X1(N1(1)),Y1(N1(2)),Z1(N1(3)),T1(N1(4))
				! BKG GRID POSITIONS
  REAL,    INTENT(IN) :: X2(N2(1)),Y2(N2(2)),Z2(N2(3)),T2(N2(4))
				! FINE GRID POSITIONS
  REAL,    INTENT(IN) :: G1(N1(1),N1(2),N1(3),N1(4),NV)
				! BACKGROUND VARIABLES
  REAL,    INTENT(OUT) :: G2(N2(1),N2(2),N2(3),N2(4),NV)

  ! LOCAL VARIABLES:
  INTEGER :: I,J,K,L
  INTEGER :: IDX(2,N2(1)),IDY(2,N2(2)),IDZ(2,N2(3)),IDT(2,N2(4))
  REAL    :: COX(2,N2(1)),COY(2,N2(2)),COZ(2,N2(3)),COT(2,N2(4)),S

  ! COEFFICIENTS AND INDICES:
  ! LEFT:
  IDX(1:2,1) = 1
  IDY(1:2,1) = 1
  IDZ(1:2,1) = 1
  IDT(1:2,1) = 1
  COX(1,1) = 1.0
  COY(1,1) = 1.0
  COZ(1,1) = 1.0
  COT(1,1) = 1.0
  COX(2,1) = 0.0
  COY(2,1) = 0.0
  COZ(2,1) = 0.0
  COT(2,1) = 0.0
  ! RIGHT:
  IDX(1:2,N2(1)) = N1(1)
  IDY(1:2,N2(2)) = N1(2)
  IDZ(1:2,N2(3)) = N1(3)
  IDT(1:2,N2(4)) = N1(4)
  COX(1,N2(1)) = 1.0
  COY(1,N2(2)) = 1.0
  COZ(1,N2(3)) = 1.0
  COT(1,N2(4)) = 1.0
  COX(2,N2(1)) = 0.0
  COY(2,N2(2)) = 0.0
  COZ(2,N2(3)) = 0.0
  COT(2,N2(4)) = 0.0

  ! X: UNIFORM GRID
  S = X1(2)-X1(1)
  DO I=2,N2(1)-1
    IDX(1,I) = INT( (X2(I)-X1(1))/S ) +1
    IDX(2,I) = IDX(1,I)+1
    COX(2,I) = (X2(I)-X1(IDX(1,I)))/S
    COX(1,I) = 1.0-COX(2,I)
  ENDDO

  ! Y: UNIFORM GRID
  S = Y1(2)-Y1(1)
  DO J=2,N2(2)-1
    IDY(1,J) = INT( (Y2(J)-Y1(1))/S ) +1
    IDY(2,J) = IDY(1,J)+1
    COY(2,J) = (Y2(J)-Y1(IDY(1,J)))/S
    COY(1,J) = 1.0-COY(2,J)
  ENDDO

  ! Z: NON-UNIFORM GRID
  DO K=2,N2(3)-1
    IDZ(1,K) = IDZ(1,K-1)
    IDZ(2,K) = IDZ(1,K)+1
    DO WHILE ( (IDZ(1,K) .LT. N1(3)) .AND. ( &
      ( (Z1(1) .LT. Z1(2)) .AND. (Z2(K) .GE. Z1(IDZ(2,K))) ) .OR. &
      ( (Z1(1) .GT. Z1(2)) .AND. (Z2(K) .LE. Z1(IDZ(2,K))) ) ) )
      IDZ(1:2,K) = IDZ(1:2,K)+1
    ENDDO
  ENDDO

  ! T: UNIFORM GRID:
  IF (N1(4) .GT. 1) S = T1(2)-T1(1)
  DO L=2,N2(4)-1
    IDT(1,L) = INT( (T2(L)-T1(1))/S ) +1
    IDT(2,L) = IDT(1,L)+1
    COT(2,L) = (T2(L)-T1(IDT(1,L)))/S
    COT(1,L) = 1.0-COT(2,L)
  ENDDO

  ! INTERPOLATIONS:
  DO L=1,N2(4)
  DO K=1,N2(3)
  DO J=1,N2(2)
  DO I=1,N2(1)
    G2(I,J,K,L,1:NV) = &
      G1(IDX(1,I),IDY(1,J),IDZ(1,K),IDT(1,L),1:NV)*COX(1,I)*COY(1,J)*COZ(1,K)*COT(1,L)+ &

      G1(IDX(2,I),IDY(1,J),IDZ(1,K),IDT(1,L),1:NV)*COX(2,I)*COY(1,J)*COZ(1,K)*COT(1,L)+ &
      G1(IDX(1,I),IDY(2,J),IDZ(1,K),IDT(1,L),1:NV)*COX(1,I)*COY(2,J)*COZ(1,K)*COT(1,L)+ &
      G1(IDX(1,I),IDY(1,J),IDZ(2,K),IDT(1,L),1:NV)*COX(1,I)*COY(1,J)*COZ(2,K)*COT(1,L)+ &
      G1(IDX(1,I),IDY(1,J),IDZ(1,K),IDT(2,L),1:NV)*COX(1,I)*COY(1,J)*COZ(1,K)*COT(2,L)+ &

      G1(IDX(2,I),IDY(2,J),IDZ(1,K),IDT(1,L),1:NV)*COX(2,I)*COY(2,J)*COZ(1,K)*COT(1,L)+ &
      G1(IDX(2,I),IDY(1,J),IDZ(2,K),IDT(1,L),1:NV)*COX(2,I)*COY(1,J)*COZ(2,K)*COT(1,L)+ &
      G1(IDX(2,I),IDY(1,J),IDZ(1,K),IDT(2,L),1:NV)*COX(2,I)*COY(1,J)*COZ(1,K)*COT(2,L)+ &
      G1(IDX(1,I),IDY(2,J),IDZ(2,K),IDT(1,L),1:NV)*COX(1,I)*COY(2,J)*COZ(2,K)*COT(1,L)+ &
      G1(IDX(1,I),IDY(2,J),IDZ(1,K),IDT(2,L),1:NV)*COX(1,I)*COY(2,J)*COZ(1,K)*COT(2,L)+ &
      G1(IDX(1,I),IDY(1,J),IDZ(2,K),IDT(2,L),1:NV)*COX(1,I)*COY(1,J)*COZ(2,K)*COT(2,L)+ &

      G1(IDX(2,I),IDY(2,J),IDZ(2,K),IDT(1,L),1:NV)*COX(2,I)*COY(2,J)*COZ(2,K)*COT(1,L)+ &
      G1(IDX(2,I),IDY(2,J),IDZ(1,K),IDT(2,L),1:NV)*COX(2,I)*COY(2,J)*COZ(1,K)*COT(2,L)+ &
      G1(IDX(2,I),IDY(1,J),IDZ(2,K),IDT(2,L),1:NV)*COX(2,I)*COY(1,J)*COZ(2,K)*COT(2,L)+ &
      G1(IDX(1,I),IDY(2,J),IDZ(2,K),IDT(2,L),1:NV)*COX(1,I)*COY(2,J)*COZ(2,K)*COT(2,L)+ &

      G1(IDX(2,I),IDY(2,J),IDZ(2,K),IDT(2,L),1:NV)*COX(2,I)*COY(2,J)*COZ(2,K)*COT(2,L)
  ENDDO
  ENDDO
  ENDDO
  ENDDO

END SUBROUTINE GRIDTRANS

SUBROUTINE BKGTOFINE(NS,NB,XB,YB,ZB,TB,NF,XF,YF,ZF,TF,BK,FG)
!*************************************************
!  THIS ROUTINE MAPS MODEL BACKGROUND VARIABLES TO 
!  THE FINEST GRID OF STMAS MULTIGRID.
!
!  HISTORY: APR. 2008, CODED BY YUANFU XIE.
!*************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: NS	! NUMBER OF VAR TO MAP
  INTEGER, INTENT(IN) :: NB(4)	! NUMBERS OF BK GRIDPOINT
  INTEGER, INTENT(IN) :: NF(4)	! NUMBERS OF FINE GRIDPTS

  REAL,    INTENT(IN) :: XB(NB(1)),YB(NB(2)),ZB(NB(3)),TB(NB(4))
				! BKG GRID POSITIONS
  REAL,    INTENT(IN) :: XF(NF(1)),YF(NF(2)),ZF(NF(3)),TF(NF(4))
				! FINE GRID POSITIONS
  REAL,    INTENT(IN) :: BK(NB(1),NB(2),NB(3),NB(4),NS)
				! BACKGROUND VARIABLES
  REAL,    INTENT(OUT) :: FG(NF(1),NF(2),NF(3),NF(4),NS)

  ! LOCAL VARIABLES:
  INTEGER :: I,J,K,L,ID(2,4)
  REAL    :: CO(2,4)

  ! FOR ALL GRID POINTS OVER THE FINEST GRID:
  ID(1,4) = 1
  ID(2,4) = MIN0(2,NB(4))
  DO L=1,NF(4)

    ! FINE GRID POINTS PASSES BK GRID POINTS:
    IF (NB(4) .EQ. 1) THEN	! ONLY ONE GRIDPOINT CASE
      ID(2,4) = ID(1,4)+1
    ELSE
      DO WHILE ( (ID(1,4) .LT. NB(4)) .AND. (TF(L) .GE. TB(ID(2,4))) )
        ID(1,4) = ID(1,4)+1
        ID(2,4) = ID(1,4)+1
      ENDDO
    ENDIF

    ! COEFFICIENTS:
    IF (ID(2,4) .GT. NB(4)) THEN
      ID(1:2,4) = NB(4)
      CO(2,4) = 0.0
    ELSE
      CO(2,4) = (TF(L)-TB(ID(1,4)))/(TB(ID(2,4))-TB(ID(1,4)))
    ENDIF
    CO(1,4) = 1.0-CO(2,4)

    ! INITIALIZE GRID LEVEL 3:
    ID(1,3) = 1
    ID(2,3) = MIN0(2,NB(3))
    DO K=1,NF(3)

      ! FINE GRID POINTS PASSES BK GRID POINTS:
      IF (NB(3) .EQ. 1) THEN		! ONLY ONE GRIDPOINT CASE
        ID(2,3) = ID(1,3)+1
      ELSE
        DO WHILE ( (ID(1,3) .LT. NB(3)) .AND. ( &
	   ((ZB(1) .GT. ZB(2)) .AND. (ZF(K) .LE. ZB(ID(2,3)))) .OR. &
           ((ZB(1) .LT. ZB(2)) .AND. (ZF(K) .GE. ZB(ID(2,3)))) ) )
        ID(1,3) = ID(1,3)+1
        ID(2,3) = ID(1,3)+1
        ENDDO
      ENDIF

      ! COEFFICIENTS:
      IF (ID(2,3) .GT. NB(3)) THEN
        ID(1:2,3) = NB(3)
        CO(2,3) = 0.0
      ELSE
        CO(2,3) = (ZF(K)-ZB(ID(1,3)))/(ZB(ID(2,3))-ZB(ID(1,3)))
      ENDIF
      CO(1,3) = 1.0-CO(2,3)

      ! INITIALIZE GRID LEVEL 2:
      ID(1,2) = 1
      ID(2,2) = MIN0(2,NB(2))
      DO J=1,NF(2)

        ! FINE GRID POINTS PASSES BK GRID POINTS:
        IF (NB(2) .EQ. 1) THEN	! ONLY ONE GRIDPOINT CASE
          ID(2,2) = ID(1,2)+1
        ELSE
          DO WHILE ( (ID(1,2) .LT. NB(2)) .AND. (YF(J) .GE. YB(ID(2,2))) )
            ID(1,2) = ID(1,2)+1
            ID(2,2) = ID(1,2)+1
          ENDDO
        ENDIF

        ! COEFFICIENTS:
        IF (ID(2,2) .GT. NB(2)) THEN
          ID(1:2,2) = NB(2)
          CO(2,2) = 0.0
        ELSE
          CO(2,2) = (YF(J)-YB(ID(1,2)))/(YB(ID(2,2))-YB(ID(1,2)))
        ENDIF
        CO(1,2) = 1.0-CO(2,2)

        ! INITIALIZE GRID LEVEL 2:
        ID(1,1) = 1
        ID(2,1) = MIN0(2,NB(1))
        DO I=1,NF(1)

          ! FINE GRID POINTS PASSES BK GRID POINTS:
          IF (NB(1) .EQ. 1) THEN	! ONLY ONE GRIDPOINT CASE
            ID(2,1) = ID(1,1)+1
          ELSE
            DO WHILE ( (ID(1,1) .LT. NB(1)) .AND. (XF(I) .GE. XB(ID(2,1))) )
              ID(1,1) = ID(1,1)+1
              ID(2,1) = ID(1,1)+1
            ENDDO
          ENDIF

          ! COEFFICIENTS:
          IF (ID(2,1) .GT. NB(1)) THEN
            ID(1:2,1) = NB(1)
            CO(2,1) = 0.0
          ELSE
            CO(2,1) = (XF(I)-XB(ID(1,1)))/(XB(ID(2,1))-XB(ID(1,1)))
          ENDIF
          CO(1,1) = 1.0-CO(2,1)

	  FG(I,J,K,L,1:NS) = &
	    BK(ID(1,1),ID(1,2),ID(1,3),ID(1,4),1:NS)*CO(1,1)*CO(1,2)*CO(1,3)*CO(1,4)+ &

            BK(ID(2,1),ID(1,2),ID(1,3),ID(1,4),1:NS)*CO(2,1)*CO(1,2)*CO(1,3)*CO(1,4)+ &
            BK(ID(1,1),ID(2,2),ID(1,3),ID(1,4),1:NS)*CO(1,1)*CO(2,2)*CO(1,3)*CO(1,4)+ &
            BK(ID(1,1),ID(1,2),ID(2,3),ID(1,4),1:NS)*CO(1,1)*CO(1,2)*CO(2,3)*CO(1,4)+ &
            BK(ID(1,1),ID(1,2),ID(1,3),ID(2,4),1:NS)*CO(1,1)*CO(1,2)*CO(1,3)*CO(2,4)+ &

            BK(ID(2,1),ID(2,2),ID(1,3),ID(1,4),1:NS)*CO(2,1)*CO(2,2)*CO(1,3)*CO(1,4)+ &
            BK(ID(2,1),ID(1,2),ID(2,3),ID(1,4),1:NS)*CO(2,1)*CO(1,2)*CO(2,3)*CO(1,4)+ &
            BK(ID(2,1),ID(1,2),ID(1,3),ID(2,4),1:NS)*CO(2,1)*CO(1,2)*CO(1,3)*CO(2,4)+ &
            BK(ID(1,1),ID(2,2),ID(2,3),ID(1,4),1:NS)*CO(1,1)*CO(2,2)*CO(2,3)*CO(1,4)+ &
            BK(ID(1,1),ID(2,2),ID(1,3),ID(2,4),1:NS)*CO(1,1)*CO(2,2)*CO(1,3)*CO(2,4)+ &
            BK(ID(1,1),ID(1,2),ID(2,3),ID(2,4),1:NS)*CO(1,1)*CO(1,2)*CO(2,3)*CO(2,4)+ &

            BK(ID(2,1),ID(2,2),ID(2,3),ID(1,4),1:NS)*CO(2,1)*CO(2,2)*CO(2,3)*CO(1,4)+ &
            BK(ID(2,1),ID(2,2),ID(1,3),ID(2,4),1:NS)*CO(2,1)*CO(2,2)*CO(1,3)*CO(2,4)+ &
            BK(ID(2,1),ID(1,2),ID(2,3),ID(2,4),1:NS)*CO(2,1)*CO(1,2)*CO(2,3)*CO(2,4)+ &
            BK(ID(1,1),ID(2,2),ID(2,3),ID(2,4),1:NS)*CO(1,1)*CO(2,2)*CO(2,3)*CO(2,4)+ &

            BK(ID(2,1),ID(2,2),ID(2,3),ID(2,4),1:NS)*CO(2,1)*CO(2,2)*CO(2,3)*CO(2,4)
       ENDDO
      ENDDO
    ENDDO
  ENDDO
END SUBROUTINE BKGTOFINE

SUBROUTINE DIRECTION(U,V,DG)
!*************************************************
! CALCULATE DIRECTION (GENERAL) (NOT USED)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL  :: U,V,DG,RA
! --------------------
  RA=ATAN(1.0D0)/45.0D0
  IF(U.EQ.0.0D0)THEN
    IF(V.GT.0.0D0)THEN
      DG=90.0D0
    ELSE
      DG=270.0D0
    ENDIF
  ELSE
    DG=ATAN2(V,U)/RA
    IF(DG.LT.0.0D0)DG=DG+360.0D0
  ENDIF
  RETURN
END SUBROUTINE DIRECTION

SUBROUTINE GETOBDATE(AD,DH,BD)
!*************************************************
! CALCULATE OBSERVATION DATE, MODIFIED FROM 'raob2dwl.f' (GENERAL) (NOT USED)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!          MARCH 2008, MODIFIED BY ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER  :: MN(12)
  INTEGER  :: DT,IY,IM,ID
  REAL  :: AD,BD,HR,DH
  LOGICAL      :: FG
  DATA MN/31,28,31,30,31,30,31,31,30,31,30,31/
!  REAL ,EXTERNAL :: RELDAY
! --------------------
  DT = NINT(AD)
  IY = MOD(DT/1000000,100 )
  IM = MOD(DT/10000  ,100 )
  ID = MOD(DT/100    ,100 )
  HR = AD - DT/100*100 + DH
  IF(MOD(IY,4).EQ.0) MN(2) = 29
  FG=.TRUE.
  DO WHILE (FG)
    IF(HR.LT.0) THEN
      HR = HR+24
      ID = ID-1
      IF(ID.EQ.0) THEN
        IM = IM-1
        IF(IM.EQ.0) THEN
          IM = 12
          IY = IY-1
          IF(IY.LT.0) IY = 99
        ENDIF
        ID = MN(IM)
      ENDIF
    ELSEIF(HR.GE.24) THEN
      HR = HR-24
      ID = ID+1
      IF(ID.GT.MN(IM)) THEN
        ID = 1
        IM = IM+1
        IF(IM.GT.12) THEN
          IM = 1
          IY = MOD(IY+1,100)
        ENDIF
      ENDIF
    ELSE
      FG=.FALSE.
    ENDIF
  ENDDO
!  BD = IY*1000000 + IM*10000 + ID*100 + HR
!  IF(ABS(BD-1.0D0*IDNINT(BD)).LE..01)BD=1.0D0*IDNINT(BD)
!  IY=IY+2000
!  BD=RELDAY(0,0,0,ID,IM,IY)-RELDAY(0,0,0,1,1,2000)
!  BD=BD+HR/24.0
  RETURN
END SUBROUTINE GETOBDATE

SUBROUTINE OBSTOGRID(Y0,X0,Y00,X00,IM,JM,X,Y,IS)
!*************************************************
! CALCULATE THE POSITION OF OBSERVATION AT THE MODEL GRID.
! HISTORY: FEBURARY 2008, CODED by ZHONGJIE HE.
!*************************************************
  IMPLICIT NONE
  INTEGER  :: I,J
  INTEGER , INTENT(IN)  :: IM,JM
  REAL    , INTENT(IN)  :: X0,Y0
  REAL    , INTENT(IN)  :: X00(IM,JM),Y00(IM,JM)
  REAL    , INTENT(OUT) :: X,Y
  INTEGER , INTENT(OUT) :: IS
!--------------------
  IS=0
  X=0
  Y=0

  IF(X0.LE.X00(1,1)) THEN
    IF(ABS(X00(1,1)-X00(2,1)).LE.1E-10) RETURN 
    X=1-(X0-X00(1,1))/(X00(1,1)-X00(2,1))
  ENDIF
  IF(X0.GE.X00(IM,JM)) THEN
    IF(ABS(X00(IM,JM)-X00(IM-1,JM)).LE.1E-10) RETURN
    X=IM+(X0-X00(IM,JM))/(X00(IM,JM)-X00(IM-1,JM))
  ENDIF
  DO I=2,IM
    IF(X0.GE.X00(I-1,1) .AND. X0.LT.X00(I,1)) THEN
      IF(ABS(X00(I-1,1)-X00(I,1)).LE.1E-10) RETURN
      X=I-1+(X0-X00(I-1,1))/(X00(I,1)-X00(I-1,1))
    ENDIF
  ENDDO  

  IF(Y0.LE.Y00(1,1)) THEN
    IF(ABS(Y00(1,1)-Y00(1,2)).LE.1E-10) RETURN
    Y=1-(Y0-Y00(1,1))/(Y00(1,1)-Y00(1,2))
  ENDIF
  IF(Y0.GE.Y00(IM,JM)) THEN
    IF(ABS(Y00(IM,JM)-Y00(IM,JM-1)).LE.1E-10) RETURN
    Y=JM+(Y0-Y00(IM,JM))/(Y00(IM,JM)-Y00(IM,JM-1))
  ENDIF
  DO J=2,JM
    IF(Y0.GE.Y00(1,J-1) .AND. Y0.LT.Y00(1,J)) THEN
      IF(ABS(Y00(1,J-1)-Y00(1,J)).LE.1E-10) RETURN
      Y=J-1+(Y0-Y00(1,J-1))/(Y00(1,J)-Y00(1,J-1))
    ENDIF
  ENDDO
  
  IF(X.GE.1 .AND. X.LE.IM .AND. Y.GE.1 .AND. Y.LE.JM) IS=1

  RETURN
END SUBROUTINE OBSTOGRID

END MODULE GENERALTOOLS
