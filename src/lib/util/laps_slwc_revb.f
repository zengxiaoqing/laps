cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

      Subroutine LAPS_SLWC_REVb(CB_PA,CB_K,GT_PA,GT_K,CT_K,
     1           ADIABATIC_LWC,ADJUSTED_LWC,ADJUSTED_SLWC,
     2           I_STATUS1,I_STATUS2)
C
C.......................HISTORY.............................
C
C        WRITTEN: CA. 1982 BY W. A. COOPER IN HP FORTRAN 4
C
C....... CALCULATES TEMPERATURE T AND LIQUID WATER CONTENT FROM
C..      CLOUD BASE PRESSURE P0 AND TEMPERATURE T0, FOR ADIABATIC
C..      ASCENT TO THE PRESSURE P.
C..     ->  INPUT:  CLOUD BASE PRESSURE P0 AND TEMPERATURE T0
C..                 PRESSURE AT OBSERVATION LEVEL P
C..     ->  OUTPUT: "ADIABATIC" TEMPERATURE T AND LIQUID WATER CONTENT
C
C        MODIFIED: November 1989 by Paul Lawson for LAPS/WISP.  Routine
c                  now calculates adiabatic liquid water content
c                  (ADIABATIC_LWC) using cloud base pressure and grid-top
c                  temperature and pressure.  Also calculated are ADJUSTED_LWC,
c                  which adjusts ADIABATIC_LWC using an empirical cloud
c                  water depletion algorithm, and ADJUSTED_SLWC, which is
c                  ADIABATIC_LWC in regions where T < 0 C adjusted
c                  using an empirical algorithm by Marcia Politovich.
c
c                  Subroutine is now hardwired for stratiform cloud only.
c                  Can be modified to include Cu with input from LAPS main.
c
c                  revb: ca 12/89 Calculate adiabatic lwc by going from cloud
c                        base to LAPS grid level instead to cloud top, thus
c                        helping to better calculate in layer clouds.
c                        Add TG (grid temperature) to calcualtion.
c
c                  revc: 2/27/90 Correct error in code.  Zero-out slwc when grid
c                        temperature (GT) > 0.
c
c                  revd: 1/27/99 Correct apparent error in test for I_STATUS1
c                        (Steve Albers)
c
c
c        OUTPUTS:  ADIABATIC_LWC
c                  ADJUSTED_LWC
c                  ADJUSTED_SLWC
c                  I_STATUS1 - 1 when -20 < cld_top_temp < 0 for Stratus
c                              0 Otherwise
c                  I_STATUS2 - 1 when valid input data provided from main
c
      DATA EPS/0.622/,CPD/1.0042E3/,CW/4.218E3/,RD/287.05/,ALHV/2.501E6/
      Integer CTY
      Integer I_STATUS1, I_STATUS2
      I_STATUS1=1
      I_STATUS2=1
c   2 Print *,'ENTER: P-BASE(mb), T-BASE(C), P-TOP, T-TOP, CLD TYPE'
c     READ(5,*) P0, T0, P, CTT, CTY
c     If(CTY.ne.0.and.CTY.ne.1) Go to 2
c
c     Hardwire cloud type (CTY) for stratus for now
c
      CTY=0
c
c.....Convert Pa to mb and Kelvin to Celcius
c
      P0 = CB_PA/100.
      P  = GT_PA/100.
      T0 = CB_K - 273.15
      TG = GT_K - 273.15
      CTT= CT_K - 273.15
c     Print *, 'CTT in Sub = ', CTT
c
c     Check for valid input data...
c
        If(P0.gt.1013..or.P0.lt.50.) Then
          I_STATUS2=0
          RETURN
        Else
        Endif
c
c
        If(T0.gt.50..or.T0.lt.-70.) Then
          I_STATUS2=0
          RETURN
        Else
        Endif
c
c
        If(P.gt.1013..or.P.lt.50.) Then
          I_STATUS2=0
          RETURN
        Else
        Endif
c
c     Set I_STATUS1 = F if 0 < cld top < -20 C (for stratus).
c
      If(CTT.GE.0..OR.CTT.LE.-20.) I_STATUS1=0
c
      TK=T0+273.15
      E=VAPOR(T0)
      R=EPS*E/(P0-E)
      CPT=CPD+R*CW
      THETAQ=TK*(1000./(P0-E))**(RD/CPT)*EXP(ALHV*R/(CPT*TK))
C 1ST APPROX
      T1=TK
      E=VAPOR(T1-273.15)
      RV=EPS*E/(P-E)
      T1=THETAQ/((1000./(P-E))**(RD/CPT)*EXP(ALHV*RV/(CPT*T1)))
C SUCCESSIVE APPROXIMATIONS
      DO 1 I=1,10
      E=VAPOR(T1-273.15)
      RV=EPS*E/(P-E)
      T1=(THETAQ/((1000./(P-E))**(RD/CPT)*EXP(ALHV*RV/(CPT*T1)))
     $   +T1)/2.
      T=T1-273.15
C     Print *, P0,T0,P,T,E,RV,THETAQ
1     CONTINUE
C GET LWC
      E=VAPOR(T)
      RV=EPS*E/(P-E)
      TW=R-RV
      ADIABATIC_LWC=TW*P*28.9644/(8.314E7*T1)*1.E9
      If(ADIABATIC_LWC.lt.0.) ADIABATIC_LWC=0.
c     Print *, 'Adiabtic LWC = ', ADIABATIC_LWC
      If(TG.ge.0.) Then
c
      ADJUSTED_Slwc=0.                                          ! Added 2/27/90
c

         If(CTY.eq.0.) Then
           If(CTT.lt.-20.) Then
             ADJUSTED_LWC=0.
           Elseif(CTT.lt.-15..and.CTT.ge.-20.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/8.
           Elseif(CTT.lt.-10..and.CTT.ge.-15.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/4.
           Else
             ADJUSTED_LWC=ADIABATIC_LWC/2.
           Endif
         Else
           If(CTT.lt.-25.) Then
             ADJUSTED_LWC=0.
           Elseif(CTT.lt.-15..and.CTT.ge.-25.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/8.
           Elseif(CTT.lt.-10..and.CTT.ge.-15.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/4.
           Else
             ADJUSTED_LWC=ADIABATIC_LWC/2.
           Endif
         Endif
      Else
         If(CTY.eq.0.) Then
           If(CTT.lt.-20.) Then
             ADJUSTED_LWC=0.
             ADJUSTED_SLWC=0.
           Elseif(CTT.lt.-15..and.CTT.ge.-20.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/8.
             ADJUSTED_SLWC=ADIABATIC_LWC/8.
           Elseif(CTT.lt.-10..and.CTT.ge.-15.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/4.
             ADJUSTED_SLWC=ADIABATIC_LWC/4.
           Else
             ADJUSTED_LWC=ADIABATIC_LWC/2.
             ADJUSTED_SLWC=ADIABATIC_LWC/2.
           Endif
         Else
           If(CTT.lt.-25.) Then
             ADJUSTED_LWC=0.
             ADJUSTED_SLWC=0.
           Elseif(CTT.lt.-15..and.CTT.ge.-25.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/8.
             ADJUSTED_SLWC=ADIABATIC_LWC/8.
           Elseif(CTT.lt.-10..and.CTT.ge.-15.) Then
             ADJUSTED_LWC=ADIABATIC_LWC/4.
             ADJUSTED_SLWC=ADIABATIC_LWC/4.
           Else
             ADJUSTED_LWC=ADIABATIC_LWC/2.
             ADJUSTED_SLWC=ADIABATIC_LWC/2.
           Endif
         Endif
      Endif
c     Print *,'Adjusted LWC = ', ADJUSTED_LWC
c     Print *,'Adjusted SLWC = ', ADJUSTED_SLWC
      END
C  FUNCTION TO CALCULATE VAPOR PRESSURE:
C
      FUNCTION VAPOR(TFP)
C INPUT IS IN DEGREES C.  IF GT 0, ASSUMED TO BE DEW POINT.  IF
C LESS THAN 0, ASSUMED TO BE FROST POINT.
C ROUTINE CODES GOFF-GRATCH FORMULA
      TVAP=273.16+TFP
      IF(TFP.GT.0.) GO TO 1
C THIS IS ICE SATURATION VAPOR PRESSURE
      IF(TVAP.LE.0) TVAP=1E-20
      E=-9.09718*(273.16/TVAP-1.)-3.56654*ALOG10(273.16/TVAP)
     $  +0.876793*(1.-TVAP/273.16)
      VAPOR=6.1071*10.**E
      RETURN
 1    CONTINUE
C THIS IS WATER SATURATION VAPOR PRESSURE
      IF(TVAP.LE.0) TVAP=1E-20
      E=-7.90298*(373.16/TVAP-1.)+5.02808*ALOG10(373.16/TVAP)
     $  -1.3816E-7*(10.**(11.344*(1.-TVAP/373.16))-1.)
     $  +8.1328E-3*(10.**(3.49149*(1-373.16/TVAP))-1)
      VAPOR=1013.246*10.**E
      RETURN
      END

