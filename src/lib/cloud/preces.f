
      SUBROUTINE PRECES(TI,TF,X,Y,Z,MODE)
      IMPLICIT REAL*8(A,B,C,D,E,F,G,H,O,P,Q,R,S,T,U,V,W,X,Y,Z)

      T0=(TI-2415020.313417D0)/36524.21988D0
      TT=(TF-TI)/36524.21988D0
      ZT=(.011171319D0+6.767999D-6*T0)*TT+1.4641D-6*TT*TT+ 8.68D-8*TT**3
      ZZ=                              ZT+3.8349D-6*TT*TT
      TH=(.009718973D0-4.135460D-6*T0)*TT-2.0653D-6*TT*TT-20.36D-8*TT**3

      IF(MODE.EQ.1)then
          X0=X
          Y0=Y
          Z0=Z

      else ! MODE .EQ. 2
          X0=DCOS(Y)*DCOS(X)
          Y0=DSIN(Y)*DCOS(X)
          Z0=DSIN(X)

      endif

      SINZT=SIN(ZT)
C     PRINT*,'X0,Y0,Z0',X0,Y0,Z0
      COSZT=COS(ZT)
      SINTH=SIN(TH)
      COSTH=COS(TH)
      SINZZ=SIN(ZZ)
      COSZZ=COS(ZZ)
      XX=-SINZT*SINZZ+COSZT*COSTH*COSZZ
      XY=-COSZT*SINZZ-SINZT*COSTH*COSZZ
      XZ=-SINTH*COSZZ
      YX= SINZT*COSZZ+COSZT*COSTH*SINZZ
      YY= COSZT*COSZZ-SINZT*COSTH*SINZZ
      YZ=-SINTH*SINZZ
      ZX= SINTH*COSZT
      ZY=-SINTH*SINZT
      ZZ= COSTH
      X   =XX*X0+XY*Y0+XZ*Z0
      Y   =YX*X0+YY*Y0+YZ*Z0
      ZDUM=ZX*X0+ZY*Y0+ZZ*Z0
C     PRINT*,'X,YDUM,ZDUM',X,YDUM,ZDUM

      IF(MODE.EQ.1)THEN
          Z=ZDUM
C         PRINT*,'PRECES - X,Y,Z',X,Y,Z

      ELSE ! MODE .EQ. 2
          Y=ATAN3(Y,X)
          X=ASIN(ZDUM)
C         PRINT*,'PRECES RA=',Y,' DEC = ',X

      ENDIF

      RETURN
      END
