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
  INTEGER(KIND=4) :: I,J,ND,NG
  REAL(KIND=8)    :: CO(NG),AC(ND,NG),OC(ND),DD(ND,NG),DC(ND)
! --------------------
  DO I=1,ND
    IF((OC(I)-AC(I,1))*(OC(I)-AC(I,2**(I-1)+1)).GT.0.0) &
    STOP 'IS NOT IN THE DOMAIN'
  ENDDO
  DO I=1,ND
    DO J=1,NG
      DD(I,J)=AC(I,J)-OC(I)
    ENDDO
    DC(I)=DABS(AC(I,2**(I-1)+1)-AC(I,1))
    IF(DC(I).LE.0.0)STOP 'GRID ARRANGE WRONG!!!'
  ENDDO
  DO I=1,ND
    DO J=1,NG
      DD(I,J)=DD(I,J)/DC(I)
      IF(DD(I,J).GT.0.0)THEN
        DD(I,J)=DD(I,J)-1.0
      ELSE
        DD(I,J)=DD(I,J)+1.0
      ENDIF
      DD(I,J)=DABS(DD(I,J))
    ENDDO
  ENDDO
  DO J=1,NG
    CO(J)=1.0
    DO I=1,ND
      CO(J)=CO(J)*DD(I,J)
    ENDDO
  ENDDO
  RETURN
END SUBROUTINE INTERPLTN

SUBROUTINE GDERIVEIT(Z1,Z2,Z3,A,B,C)
!*************************************************
! GENERAL DERIVATIVE OF INTERIOR (GENERAL)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL(KIND=8) :: Z1,Z2,Z3,X,Y,A,B,C
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
  REAL(KIND=8) :: Z1,Z2,Z3,X,Y,A,B,C
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
  REAL(KIND=8) :: Z1,Z2,Z3,X,Y,A,B,C
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
  REAL(KIND=8) :: Z1,Z2,Z3,X,Y,A,B,C
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
  INTEGER(KIND=4) :: NN,NC,PT,M,N
  INTEGER(KIND=4) :: NP(NN),NM(NN)
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
  INTEGER(KIND=4) :: NN,NC,N,N0
  INTEGER(KIND=4) :: NP(NN),NM(NN)
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

SUBROUTINE VRTCLPSTN8(KM,PP,P,PS,IS)
!*************************************************
! GET VERTICAL POSITION (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: KM,K,IS
  REAL(KIND=8)    :: PP(KM)
  REAL(KIND=4)    :: P,PS
! --------------------
  IS=0
  DO K=1,KM-1
    IF((P-PP(K))*(P-PP(K+1)).GT.0.0)CYCLE
    PS=K+ABS(P-PP(K))/ABS(PP(K+1)-PP(K))
    IS=1
  ENDDO
  RETURN
END SUBROUTINE VRTCLPSTN8

SUBROUTINE VRTCLPSTN(KM,PP,P,PS,IS)
!*************************************************
! GET VERTICAL POSITION (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: KM,K,IS
  REAL(KIND=4)    :: PP(KM),P,PS
! --------------------
  IS=0
  DO K=1,KM-1
    IF((P-PP(K))*(P-PP(K+1)).GT.0.0)CYCLE
    PS=K+ABS(P-PP(K))/ABS(PP(K+1)-PP(K))
    IS=1
  ENDDO
  RETURN
END SUBROUTINE VRTCLPSTN

SUBROUTINE ZPCONVERT(KM,ZZ,PP,Z,P,IS)
!*************************************************
! FIND PRESSURE FROM HEIGHT (GENERAL)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: K,KM,IS
  REAL(KIND=8)    :: ZZ(KM),Z,P
  REAL(KIND=4)    :: PP(KM)
! --------------------
  IS=0
  DO K=1,KM-1
    IF((Z-ZZ(K))*(Z-ZZ(K+1)).GT.0.0)CYCLE
    P=PP(K)*DABS(ZZ(K+1)-Z)/DABS(ZZ(K+1)-ZZ(K))+PP(K+1)*DABS(Z-ZZ(K))/DABS(ZZ(K+1)-ZZ(K))
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
  INTEGER(KIND=4),PARAMETER :: MD=4
  INTEGER(KIND=4) :: ND,NG,NS
  INTEGER(KIND=4) :: NF(MD),NB(MD)
  REAL(KIND=8)    :: XF(NF(1)),YF(NF(2)),ZF(NF(3)),TF(NF(4))
  REAL(KIND=8)    :: XB(NB(1)),YB(NB(2)),ZB(NB(3)),TB(NB(4)),PP(MD)
  REAL(KIND=8)    :: FC(NF(1),NF(2),NF(3),NF(4),NS)
  REAL(KIND=8)    :: BK(NB(1),NB(2),NB(3),NB(4),NS)
  REAL(KIND=8)    :: AC(ND,NG),OC(ND),CO(NG)
  INTEGER(KIND=4) :: I,J,K,L,I0,J0,K0,L0,I1,J1,K1,L1,S,M,N,NN
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
    DO L1=L0,MIN0(L0+1,NF(4))
    DO K1=K0,MIN0(K0+1,NF(3))
    DO J1=J0,MIN0(J0+1,NF(2))
    DO I1=I0,MIN0(I0+1,NF(1))
      M=M+1
      PP(1)=XF(I1)
      PP(2)=YF(J1)
      PP(3)=ZF(K1)
      PP(4)=TF(L1)
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
  INTEGER(KIND=4) :: K,KM,KZ
  REAL(KIND=8)    :: ZZ(KM),Z
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

SUBROUTINE CONVERTUV(X1,Y1,X2,Y2,U1,V1,U2,V2)
!*************************************************
! CONVERT VECTOR FROM ONE COORDINATE TO ANOTHER (GENERAL) (NOT USED)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL(KIND=8) :: X1,Y1,X2,Y2,U1,V1,U,V,DG,U2,V2
! --------------------
  U=X2-X1
  V=Y2-Y1
  CALL DIRECTION(U,V,DG)
  U2= U1*DCOSD(DG)+V1*DSIND(DG)
  V2=-U1*DSIND(DG)+V1*DCOSD(DG)
  RETURN
END SUBROUTINE CONVERTUV

SUBROUTINE DIRECTION(U,V,DG)
!*************************************************
! CALCULATE DIRECTION (GENERAL) (NOT USED)
! HISTORY: AUGUST 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  REAL(KIND=8) :: U,V,DG,RA
! --------------------
  RA=DATAN(1.0D0)/45.0D0
  IF(U.EQ.0.0D0)THEN
    IF(V.GT.0.0D0)THEN
      DG=90.0D0
    ELSE
      DG=270.0D0
    ENDIF
  ELSE
    DG=DATAN2(V,U)/RA
    IF(DG.LT.0.0D0)DG=DG+360.0D0
  ENDIF
  RETURN
END SUBROUTINE DIRECTION

SUBROUTINE GETOBDATE(AD,DH,BD)
!*************************************************
! CALCULATE OBSERVATION DATE, MODIFIED FROM 'raob2dwl.f' (GENERAL) (NOT USED)
! HISTORY: SEPTEMBER 2007, CODED by WEI LI.
!*************************************************
  IMPLICIT NONE
! --------------------
  INTEGER(KIND=4) :: MN(12)
  INTEGER(KIND=4) :: DT,IY,IM,ID
  REAL(KIND=8) :: AD,BD,HR,DH
  DATA MN/31,28,31,30,31,30,31,31,30,31,30,31/
!  REAL(KIND=8),EXTERNAL :: RELDAY
! --------------------
  DT = IDNINT(AD)
  IY = MOD(DT/1000000,100 )
  IM = MOD(DT/10000  ,100 )
  ID = MOD(DT/100    ,100 )
  HR = AD - DT/100*100 + DH
  IF(MOD(IY,4).EQ.0) MN(2) = 29
1 IF(HR.LT.0) THEN
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
    GOTO 1
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
    GOTO 1
  ENDIF
!  BD = IY*1000000 + IM*10000 + ID*100 + HR
!  IF(DABS(BD-1.0D0*IDNINT(BD)).LE..01)BD=1.0D0*IDNINT(BD)
!  IY=IY+2000
!  BD=RELDAY(0,0,0,ID,IM,IY)-RELDAY(0,0,0,1,1,2000)
!  BD=BD+HR/24.0
  RETURN
END SUBROUTINE GETOBDATE

END MODULE GENERALTOOLS
