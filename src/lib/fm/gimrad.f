      FUNCTION GIMRAD(TAU,TEMP,TSKIN,KCHAN,LSFC,PSFC,EMISS)
c  12 dec 94  correct error in sfc tau extrapolation  
C $ FUNCTION GIMRAD(TAU,TEMP,TSKIN,KCHAN,LSFC,PSFC,EMISS)   (HMW)
C $ GIMRAD - compute GOES radiance from temperature profile
C $ Input:
C $     TAU = (R) array of GOES transmittances
C $     TEMP = (R) array of atmospheric temperatures
C $     TSKIN = (R) temperature of surface
C $     KCHAN = (I) channel number
C $     LSFC = (I) subscript of level below or at surface
C $     PSFC = (I) pressure of surface ( or cloud)
C $     EMISS = (R) surface emissivity
C $ Output description:
C $     GOES radiance
C $$ GIMRAD = SOUNDER, GOES, CONVERT
      DIMENSION TAU(*),TEMP(*)
      COMMON/ATMOS/P(40),T(40),W(40),OZO(40)
      COMMON /ATMRAD/ATRD,TAUS,REFL
C     PREPARE SURFACE CORRECTION
      DP = ALOG(PSFC/P(LSFC))/ALOG(P(LSFC)/P(LSFC-1))
C     DP IS NEGATIVE TO INTERPOLATE, POSITIVE TO EXRAPOLATE (PSFC>NL)
      DTAU = TAU(LSFC-1)-TAU(LSFC)
      TAUS = TAU(LSFC)-DTAU*DP
      REFL = 0.
      T1 = TEMP(1)
      B1 = PLANGO(T1,KCHAN)
      TAU1 = TAU(1)
      RAD = 0.
C     INTEGRATE DOWN TO LEVEL AT OR BELOW SURFACE
      DO 110 I = 2,LSFC
      T2 = TEMP(I)
      B2 = PLANGO(T2,KCHAN)
      TAU2 = TAU(I)
      DTAU = TAU1-TAU2
      DR =.5*(B1+B2)*DTAU
      IF (TAUS .GT. 0.1.AND.EMISS .LT. 1.00)THEN
C     DO NOT ADD REFLECTED FOR LAST LEVEL UNLESS PSFC .GT. 1000
      IF (I .EQ. LSFC.AND.PSFC .LE. 1000.)GO TO 105
      TAUB = 0.5*(TAU1+TAU2)
      TAUFAC = TAUS/TAUB
      REFL = REFL+TAUFAC*TAUFAC*DR
      ENDIF
  105 RAD = RAD+DR
      B1 = B2
  110 TAU1 = TAU2
C     ADD (SUBTRACT) INCREMENT OF ATMOS RADIANCE TO REACH SURFACE
C     DP WILL BE NEGATIVE IF PSFC  < 1000 MB
C     DR FALLS OUT AS THE DELTA RADIANCE OF LAYER
      RAD = RAD+DR*DP
C     ADD INCREMENT OF REFLECTED RADIANCE FOR LAYER DOWN TO SURFACE
      IF (TAUS .GT. 0.1.AND.EMISS .LT. 1.00)THEN
         IF (PSFC .LT. 1000.)THEN
         TAUB = 0.5*(TAU(LSFC-1)+TAUS)
C     CHANGE DP TO INCREMENT RATHER THAN DECREMENT
         DP = 1.+DP
         ELSE
         TAUB = 0.5*(TAU(LSFC)+TAUS)
         ENDIF
      TAUFAC = TAUS/TAUB
      REFL = REFL+TAUFAC*TAUFAC*DR*DP
      ENDIF
      ATRD = RAD
      RAD = RAD+(1.-EMISS)*REFL
      BS = PLANGO(TSKIN,KCHAN)
      RAD = RAD+EMISS*BS*TAUS
      RAD = AMAX1(RAD,.001)
      GIMRAD = RAD
      RETURN
      END
