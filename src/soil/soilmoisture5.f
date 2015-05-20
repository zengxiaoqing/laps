cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
      Subroutine Soil_Moisture(imax,jmax,
     &           Laps_u,Laps_v,Laps_T,Laps_TD,
     &           Laps_Rain,Laps_sc,Laps_IN,Laps_WFZ,
     &           Laps_MWF,Laps_MWF_pre,Laps_Wx,SoilType,    !Inputs down to here
     &           GridDry,Laps_Evap,Laps_SMC_3D,istatus)

C     LAPS SoilMoisture Subroutine. Driver is program LSM (LAPS Soil Moisture).
C     Original version implemented for FSL Demonstration PC Workstation by
C     Chandran Subramaniam
C     2/8/93
      
C     John Smart 12/1/93: Adapt the software to run in real time on the
C     			  UNIX platform.  Set up LAPS standard I/O.
c     John Smart 9/22/97: Dynamic array mods

      integer*4 imax,jmax
      Include 'soilm.inc'
C
C**** MODEL12 is a Soil Moisture Content model developed in June l986, and
C     based upon the Remote Sensing Based Watershed Model developed at
C     the University of Maryland from 1982 to 1985.  This is Version 12.4,
C     as submitted to Water Resources Bulletin in January 1989.
C     Created by Groves.


      DIMENSION R(100)
      REAL 	KSAT,LAMDA,IN
      real*4    rmiss
      DATA 	DAY,SUMR,IN,OLDWEA/4*0./
      INTEGER   Istatus
      INTEGER   icycle_time
      LOGICAL 	Dry, GridDry
      INTEGER   SoilType(Imax,Jmax)
      REAL      Laps_u(Imax,Jmax)     !u-component
      REAL      Laps_v(Imax,Jmax)     !v-component
      REAL      Laps_T(Imax,Jmax)     !Temperature
      REAL      Laps_TD(Imax,Jmax)    !Dewpoint Temperature
      REAL      Laps_SC(Imax,Jmax)    !Snow Cover - Estimate - (0.0 to 1.0)
      REAL      Laps_IN(Imax,Jmax)    !Infiltration
      REAL      Laps_Wx(Imax,Jmax)    !Weather data
      REAL      Laps_WFZ(Imax,Jmax)   !Wetting Front Depth, z
      REAL      Laps_MWF(Imax,Jmax)   !Wetting Front Moisture Content
      REAL      Laps_MWF_pre(Imax,Jmax)!Previous Wetting Front Moisture Content
      REAL 	Laps_Rain(Imax, Jmax) !Radar-estimated Liq precip
      REAL      Laps_Evap(Imax,Jmax)  !Amount of Evaporation (Calc within)
      REAL      Laps_SMC_3D(Imax,Jmax,3)!Three layer Soil Moisture Content

      call get_laps_cycle_time(icycle_time,istat)
      if(istat.ne.1)then
         write(6,*)'Error getting laps cycle time'
         write(6,*)'WARNING: Setting DELT = 1.0'
         write(6,*)'Assuming cycle time = 3600 sec'
         icycle_time = 3600
      endif
      DELT = float(icycle_time)/3600.

      call get_r_missing_data(rmiss,istat)
      if(istat.ne.1)goto 1000

      Write(6,*)'Calculating Pan Evaporation'
      Call Calc_Evap(imax,jmax,Laps_u,Laps_v,
     &               Laps_T,Laps_TD,    !Input to here.
     &               Laps_Evap,IStatus)
      If (IStatus.eq.1) Then
 	 Write(6,*)'Completed Calculating Evaporation'
      Else
         Write(6,*)'Error computing evaporation, returning'
         istatus = -1
         return
      Endif

      Write(6,*)' Begin Soil Moisture Calculation for grid'

      Do J = 1, Jmax
      Do I = 1, Imax

	ISOIL = SoilType(I,J)
        CALL SOILS(ISOIL,KSAT,THS,THR,PSIF,PSIAE,LAMDA)
        CALL AMC(THFC,ISOIL,THS)

        OLDWEA = Laps_Wx(I,J)	     
        IN = Laps_IN(I,J)
        Z = Laps_WFZ(I,J)
	THI = Laps_MWF(I,J)
        XTHI = Laps_MWF_pre(I,J)

        if(IN .ne.rmiss.and.  Z .ne.rmiss.and.
     &     THI.ne.rmiss.and.XTHI.ne.rmiss)then

        If((Laps_sc(i,j).lt.snowthres).or.(laps_sc(i,j).gt.1.e30))then

          Dry = .FALSE.
          IF((GridDry).OR.(Laps_Rain(I,J).LT.RainThres))dry=.true.
          IF (Dry) THEN
	    PAN = Laps_Evap(I,J)
C	    IF(OLDWEA.LE.0.)XTHI=THI     ! I'm Not Sure  about this statement
            CALL MOIST(XTHI,IN,XDAY,KSAT,THS,THR,PSIAE,LAMDA,THI,
     1                 DEPTH,PAN,CUMD,CUMET,BAL,Z,RZST,OLDWEA)
	    CALL SPLIT(Z,THI,XTHI,THI1,THI2,THI3,HOR1,HOR2,HOR3)
            OLDWEA = 1.0/24.
            IN = 0.0
          ELSE        ! Rain 
            IF(Z.LT.DCM)XTHI=(THI*Z+XTHI*(DCM-Z))/DCM
            R(1) = Laps_Rain(I,J)
	    CALL ENTRY(CUMQ,IN,XTHI,KSAT,THS,PSIF,DELT,R,DEPTH,N)
	    Z=IN/(THS-XTHI)
	    CALL SPLIT(Z,THS,XTHI,THI1,THI2,THI3,HOR1,HOR2,HOR3)
            OLDWEA = -1./24
            THI = THS
          ENDIF
          Laps_MWF(I,J) = THI
          Laps_MWF_pre(i,j) = XTHI
          Laps_WFZ(i,j)= Z
  	  Laps_IN(I,J) = IN
	  Laps_Wx(I,J) = OLDWEA
          if(THS.eq.0.0)then
             Laps_SMC_3D(I,J,1) = 0.0
             Laps_SMC_3D(I,J,2) = 0.0
             Laps_SMC_3D(I,J,3) = 0.0
          else
	     THSPER = THS / 100.0 
	     Laps_SMC_3D(I,J,1) = THI1 ! / THSPER
	     Laps_SMC_3D(I,J,2) = THI2 ! / THSPER
	     Laps_SMC_3D(I,J,3) = THI3 ! / THSPER

          endif

        else

          Laps_MWF(i,j)=ths
          Laps_MWF_pre(i,j)=xthi
	  Laps_WFZ(I,J) = hor1*2.54
          Laps_IN(i,j)=0.1  !Not sure about this
	  Laps_Wx(I,J) = 1./24. 

          if(THS.eq.0.0)then
             Laps_SMC_3D(I,J,1) = 0.0
             Laps_SMC_3D(I,J,2) = 0.0
             Laps_SMC_3D(I,J,3) = 0.0
          else
	     THSPER = ths / 100.0 
	     Laps_SMC_3D(I,J,1) = ths ! / THSPER
	     Laps_SMC_3D(I,J,2) = ths ! / THSPER
	     Laps_SMC_3D(I,J,3) = ths ! / THSPER
          endif

        endif

        endif

      Enddo
      Enddo	
      Write(6,*)' Completed SM Calculation for grid'
      istatus = 1
C
100   WRITE(6,*)'End of Simulation'
      goto 1000
C
1000  Return
      END

C ====================================================================

      SUBROUTINE SOILS(ISOIL,KSAT,THS,THR,PSIF,PSIAE,LAMDA)
C**** Subroutine provides default values for soil hydraulic parameters
      DIMENSION SATCON(6),RESPOR(6),EFFPOR(6),BUBPR(6)
      REAL KSAT,LAM(6),LAMDA
      DATA SATCON/10.,3.5,.65,1.3,.08,.03/
      DATA RESPOR/.04,.05,.08,.10,.08,.11/
      DATA EFFPOR/.39,.4,.39,.32,.4,.38/
      DATA BUBPR/10.,15.,27.,15.,80.,130./
      DATA LAM/.43,.38,.31,.23,.23,.20/
      THE=EFFPOR(ISOIL)
      KSAT=SATCON(ISOIL)/60.
      THR=RESPOR(ISOIL)
      PSIBUB=BUBPR(ISOIL)
      LAMDA=LAM(ISOIL)
      THS=THR+THE
      PSIAE=PSIBUB/2.
      PSIF=PSIAE*(2.+3.*LAMDA)/(1.+3.*LAMDA)
      RETURN
      END
C
      SUBROUTINE AMC(THI,ISOIL,THS)
C**** Subroutine assigns initial soil moisture value
      DIMENSION THFC(6),THWILT(6)
      DATA THWILT/.06,.075,.13,.16,.22,.27/
      DATA THFC/.125,.175,.25,.23,.37,.425/
      THI=THFC(ISOIL)
      IF(THI.LT.THWILT(ISOIL))XTHI=THWILT(ISOIL)
      IF(THI.GT.THS)XTHI=THS
      RETURN
      END
C
      SUBROUTINE ENTRY(RO,IN,THI,KSAT,THS,PSIF,DELT,R,IDEPTH,N)
C**** Subroutine calculates infiltration volume for rainfall period
      DIMENSION CUMF(100),R(100),CUMQ(100),Q(100)
      REAL KSAT,IN
      INTEGER IDEPTH
      IF(THI.GT.THS)THI=THS
      SMAX=(THS-THI)*float(IDEPTH)*2.54
      CUMF(1)=IN
C      CUMF(1)=0.
      CUMQ(1)=0.
      FC=KSAT*DELT*60.
      DO 10  IT=2,N+1
	IF(CUMF(IT-1).GT.SMAX) GO TO 2
	CALL INFIL(THI,DELT,KSAT,THS,PSIF,CUMF(IT-1),DELI)
	GO TO 4
2	DELI=FC
4	IF(DELI.GT.R(IT-1)) DELI=R(IT-1)
	CUMF(IT)=CUMF(IT-1)+DELI
	Q(IT)=R(IT-1)-DELI
	CUMQ(IT)=CUMQ(IT-1)+Q(IT)
	RO=CUMQ(IT)
10	IN=CUMF(IT)
      RETURN
      END
C
      SUBROUTINE INFIL(THI,DELT,KSAT,THS,PSIF,IN,DELI)
C**** Subroutine calculates incremental infiltration volume
      REAL KSAT,IN,ITX
      DT=DELT*60.
      BB=PSIF*(THS-THI)
      IF(IN.LE.0.)IN=0.01
      XE=EXP(-IN/BB)
      TI=(IN-BB+BB*XE)/KSAT
      TX=TI+DT
      ITX=IN
      DO 10 M=1,5
	If(ITX .GT. 5)GOTO 20
	XE=EXP(-ITX/BB)
10	ITX=ITX-(ITX-BB+BB*XE-KSAT*TX)/(1.-XE)
20    DELI=ITX-IN
      IF(IN.EQ..01)IN=0.
      RETURN
      END
C
      SUBROUTINE MOIST(XTHI,XIN,XDAY,KSAT,TS,TR,PS,LA,
     1THETA,IRD,PR,CUMD,CUMTR,BAL,Z,RZST,OLDWEA)
C**** Subroutine simulates redistribution of soil moisture between storms
      COMMON/COM1/PTR,DTHETA,THETAA,THETAB
      COMMON/COM2/KS,THS,THR,PSIAE,LAMDA,RZD,T
      REAL KSAT,KS,LAMDA,LA
      DATA PSIA,PSIB/15500.,5166./
      KS=KSAT/60.
      THS=TS
      THR=TR
      PSIAE=PS
      LAMDA=LA
      RZD=IRD*2.54
c      IF(RZD.LT.50.)RZD=50.
      PTR=PR*2.54/86400.
C**** Change time units to seconds
      TSTART=0.
      T=TSTART*86400.
      TOUT=XDAY*86400.
C**** Initialize THETA, Z, storages
      IF(OLDWEA.LT.0.)THEN
	 THETA=THS
	 Z=XIN/(THETA-XTHI)
	 OLDRZS=THETA*RZD
	 IF(Z.LT.RZD)OLDRZS=Z*THETA+(RZD-Z)*XTHI
      ELSE
	 OLDRZS=RZST
      ENDIF
      CUMD=0.
      CUMTR=0.
C**** Determine ET thresholds
      PSIBUB=PSIAE*2.
      THE=THS-THR
      THETAA=THE*(PSIBUB/PSIA)**LAMDA+THR
c ******************************************
c  Hard wire minimum theta to be lower here. This sets it
c    to the wilting point.
c  Try 20% of THS first
      THETAA=.2*THS
c ******************************************
      THETAB=THE*(PSIBUB/PSIB)**LAMDA+THR
C**** Advance soil moisture from TSTART to XDAY
      CALL INTGRL(THETA,XTHI,Z,CUMD,CUMTR,TOUT)
C**** Calculate storages
      RZST=THETA*RZD
      IF(Z.LT.RZD)RZST=Z*THETA+(RZD-Z)*XTHI
      BAL=-CUMTR-CUMD-(RZST-OLDRZS)
C     IF(Z.EQ.0.)Z=RZD
      RETURN
      END
C
      SUBROUTINE INTGRL(THETA,XTHI,Z,CUMD,CUMTR,TOUT)
C**** Subroutine integrates soil moisture variables
      COMMON/COM1/PTR,DTHETA,THETAA,THETAB
      COMMON/COM2/KS,THS,THR,PSIAE,LAMDA,RZD,T
      DIMENSION X(5),DX(5),DX1(5),X1(5)
      REAL KS,LAMDA
      X(1)=THETA
      X(2)=XTHI
      X(3)=Z
      X(4)=CUMD
      X(5)=CUMTR
      DTHETA=1.
      DX1(1)=1.
      GO TO 10
1     CONTINUE
C**** Check for small change in THETA at wetting front
      IF(Z.EQ.0.)GO TO 10
      IF(DTHETA.GT..03.AND.Z.LT.(.9*RZD))GO TO 10
      IF(Z.LT.RZD)X(1)=X(2)+X(3)*(X(1)-X(2))/RZD
      X(2)=X(1)
      X(3)=0.
C**** Calculate DT
10    DT=AMIN1(-.004/DX1(1),TOUT-T)
      IF(DT.LT.1.)DT=1.
      DT2=DT/2.
C**** Integration routine
      CALL RATE (X(1),X(2),X(3),DX1(1),DX1(2),DX1(3),DX1(4),DX1(5))
      DO 20 I=1,5
20    X1(I)=X(I)+DT2*DX1(I)
      T=T+DT
      CALL RATE (X(1),X(2),X(3),DX(1),DX(2),DX(3),DX(4),DX(5))
      DO 30 I=1,5
30    X(I)=X(I)+DT2*(DX1(I)+DX(I))
      THETA=X(1)
      XTHI=X(2)
      Z=X(3)
      CUMD=X(4)
      CUMTR=X(5)
      IF(T.LT.TOUT)GO TO 1
      RETURN
      END
C
      SUBROUTINE RATE(THETA,XTHI,Z,DTHDT,DTHIDT,DZDT,DRN,TR)
C**** Subroutine calculates rates of change of soil moisture redistribution
      COMMON/COM1/PTR,DTHETA,THETAA,THETAB
      COMMON/COM2/KS,THS,THR,PSIAE,LAMDA,RZD,T
      REAL KS,LAMDA,LAM23
      THE=THS-THR
      C1=PSIAE/(1.+3.*LAMDA)
      LAM23=3.+2./LAMDA
      TR1=0.
      IF(Z.EQ.0.)GO TO 10
C**** Determine rates above the wetting front
      DTHETA=THETA-XTHI
      If (Dtheta .EQ. 0) Then
        Z = 0
	Go To 10
      Endif
      QZ=KS*(C1*((THETA-THR)/THE)**(3.+1./LAMDA)/Z+((THETA-THR)
     1/THE)**LAM23)
      DZDT=QZ/DTHETA
      F1=1.
      IF(THETA.LT.THETAB)F1=(THETA-THETAA)/(THETAB-THETAA)
      IF(THETA.LT.THETAA)F1=0.
      TR1=F1*PTR*Z/RZD
      DTHDT=-((QZ)+TR1)/Z
C**** Determine ET below the wetting front
10    F1=1.
      IF(XTHI.LT.THETAB)F1=(XTHI-THETAA)/(THETAB-THETAA)
      IF(XTHI.LT.THETAA)F1=0.
      TR2=F1*PTR*(RZD-Z)/RZD
C**** Determine moisture loss rates
      TR=TR1+TR2
      DRN=KS*((XTHI-THR)/THE)**LAM23
      DTHIDT=-(DRN+TR2)/(RZD-Z)
      IF(Z.NE.0.)RETURN
      DZDT=0.
      DTHDT=DTHIDT
      RETURN
      END
C
      SUBROUTINE SPLIT(Z,THI,XTHI,THI1,THI2,THI3,HOR1,HOR2,HOR3)
C**** Subroutine calculates moisture content for three soil horizons
      H1=HOR1*2.54
      H2=HOR2*2.54
      H3=HOR3*2.54
      IF(Z.LE.H1)THEN
	THI1=(Z*THI+(H1-Z)*XTHI)/H1
	THI2=XTHI
	THI3=XTHI
	RETURN
      ELSE IF (Z.LE.H2)THEN
	THI1=THI
	THI2=((Z-H1)*THI+(H2-Z)*XTHI)/(H2-H1)
	THI3=XTHI
	RETURN
      ELSE IF(Z.LE.H3)THEN
	THI1=THI
	THI2=THI
	THI3=((Z-H2)*THI+(H3-Z)*XTHI)/(H3-H2)
	RETURN
      ELSE
	THI1=THI
	THI2=THI
	THI3=THI
      ENDIF
      RETURN
      END	
