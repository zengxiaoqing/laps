      SUBROUTINE ETALL(im,jm,tph0d,tlm0d,dlmd,dphd,HLAT,HLON,VLAT,VLON)
C$$$  SUBPROGRAM DOCUMENTATION BLOCK
C                .      .    .                                       .
C SUBPROGRAM: ETALL         COMPUTE EARTH LATITUDE & LONIGTUDE OF
C                           ETA GRID POINTS
C   PRGMMR: ROGERS          ORG: W/NP22     DATE: 90-06-13
C
C ABSTRACT: COMPUTES THE EARTH LATITUDE AND LONGITUDE OF ETA GRID
C   POINTS (BOTH H AND V POINTS)
C
C PROGRAM HISTORY LOG:
C   90-06-13  E.ROGERS
C   98-06-09  M.BALDWIN - CONVERT TO 2-D CODE
C   01-01-03  T BLACK   - MODIFIED FOR MPI
C
C USAGE:    CALL ETALL(HLAT,HLON,VLAT,VLON)
C   INPUT ARGUMENT LIST:
C     NONE
C
C   OUTPUT ARGUMENT LIST:
C     HLAT     - LATITUDE OF H GRID POINTS IN RADIANS (NEG=S)
C     HLON     - LONGITUDE OF H GRID POINTS IN RADIANS (E)
C     VLAT     - LATITUDE OF V GRID POINTS IN RADIANS (NEG=S)
C     VLON     - LONGITUDE OF V GRID POINTS IN RADIANS (E)
C
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM RS/6000 SP
C-----------------------------------------------------------------------
C-----------------------------------------------------------------------
                             D I M E N S I O N
     & IDAT  (3)
C
     &,GLATH (IM,JM),GLONH (IM,JM),GLATV (IM,JM),GLONV (IM,JM)
C
        REAL HLAT(IM,JM),HLON(IM,JM),VLAT(IM,JM),VLON(IM,JM)
C-----------------------------------------------------------------------
C--------------LOGICAL FILE NAMES---------------------------------------
C-----------------------------------------------------------------------
                             D A T A
     &  LIST  /06/
C-----------------------------------------------------------------------
C--------------UNIVERSAL CONSTANTS--------------------------------------
C-----------------------------------------------------------------------
                             D A T A
     & PI/3.141592654/
C-----------------------------------------------------------------------
 1000 FORMAT(' I J=',2I4,' GLAT=',E12.5,' GLON=',E12.5,' ELAT=',E12.5
     &,' ELON=',E12.5,' WLON = ',E12.5)
C----------------------------------------------------------------------
C--------------DERIVED GEOMETRICAL CONSTANTS----------------------------
C----------------------------------------------------------------------
CC
	WBD=-(IM-1)*DLMD
	SBD=-(JM-1)/2*DPHD
      DTR = PI / 180.0
      TPH0 = TPH0D * DTR
      WB = WBD * DTR
      SB = SBD * DTR
      DLM = DLMD * DTR
      DPH = DPHD * DTR

	write(6,*) 'TPH0,WB,SB,DLM,DPH,DTR: ', TPH0,WB,SB,DLM,DPH,DTR

      TDLM = DLM + DLM
      TDPH = DPH + DPH
C
      STPH0 = SIN(TPH0)
      CTPH0 = COS(TPH0)
C
C-----------------------------------------------------------------------
C---COMPUTE GEOGRAPHIC LAT AND LONG OF ETA GRID POINTS (H & V POINTS)---
C-----------------------------------------------------------------------
      DO 200 J = 1,JM
C
         TLMH = WB - TDLM + MOD(J+1,2) * DLM
         TPHH = SB+(J-1)*DPH
         TLMV = WB - TDLM + MOD(J,2) * DLM
         TPHV = TPHH
         STPH = SIN(TPHH)
         CTPH = COS(TPHH)
         STPV = SIN(TPHV)
         CTPV = COS(TPHV)
C----------------------------------------------------------------------
C---------- COMPUTE EARTH LATITUDE/LONGITUDE OF H POINTS --------------
C----------------------------------------------------------------------
         DO 201 I = 1,IM
           TLMH = TLMH + TDLM
           SPHH = CTPH0 * STPH + STPH0 * CTPH * COS(TLMH)
           GLATH(I,J) = ASIN(SPHH)
           CLMH = CTPH * COS(TLMH) / (COS(GLATH(I,J)) * CTPH0)
     1               - TAN(GLATH(I,J)) * TAN(TPH0)
           IF(CLMH .GT. 1.) CLMH = 1.0
           IF(CLMH .LT. -1.) CLMH = -1.0
           FACTH = 1.
           IF(TLMH .GT. 0.) FACTH = -1.
C          WRITE(6,88888) I,J, CLMH
C8888      FORMAT(2X,2I6,1X,E12.5)
           GLONH(I,J) = -TLM0D * DTR + FACTH * ACOS(CLMH)
C          IF(I .EQ. 1) THEN
C           WRITE(LIST,99995) I,J,GLATH(I,J),GLONH(I,J)
C9995       FORMAT(2X,2(I6,1X),2(E12.5,1X))
C          END IF
C
C
           HLAT(I,J) = GLATH(I,J) / DTR
C           HLON(I,J) = 360.0 - GLONH(I,J) / DTR
	    HLON(I,J)= -GLONH(I,J)/DTR
	   
           IF(HLON(I,J) .GT. 180.) HLON(I,J) = HLON(I,J) - 360.
           IF(HLON(I,J) .LT. -180.) HLON(I,J) = HLON(I,J) + 360.
  201    CONTINUE
C----------------------------------------------------------------------
C---------- COMPUTE EARTH LATITUDE/LONGITUDE OF V POINTS --------------
C----------------------------------------------------------------------
         DO 202 I = 1,IM
           TLMV = TLMV + TDLM
           SPHV = CTPH0 * STPV + STPH0 * CTPV * COS(TLMV)
           GLATV(I,J) = ASIN(SPHV)
           CLMV = CTPV * COS(TLMV) / (COS(GLATV(I,J)) * CTPH0)
     1          - TAN(GLATV(I,J)) * TAN(TPH0)
           IF(CLMV .GT. 1.) CLMV = 1.
           IF(CLMV .LT. -1.) CLMV = -1.
           FACTV = 1.
           IF(TLMV .GT. 0.) FACTV = -1.
           GLONV(I,J) = -TLM0D * DTR + FACTV * ACOS(CLMV)
C          IF(I.EQ.1) THEN
C           WRITE(LIST,99995) I,J,GLATV(I,J),GLONV(I,J)
C          END IF
C
C    CONVERT INTO DEGREES AND EAST LONGITUDE
C
           VLAT(I,J) = GLATV(I,J) / DTR
C           VLON(I,J) = 360.0 - GLONV(I,J) / DTR
           VLON(I,J) = -GLONV(I,J) / DTR
           IF(VLON(I,J) .GT. 180.) VLON(I,J) = VLON(I,J) - 360.
           IF(VLON(I,J) .LT. -180.) VLON(I,J) = VLON(I,J) + 360.
  202    CONTINUE
  200 CONTINUE
C
C     DO 210 J = 1,JM
C         WRITE(LIST,88888) J,HLAT(1,J),HLON(1,J),VLAT(1,J),VLON(1,J)
C8888    FORMAT(2X,I5,1X,4(E12.5,1X))
C      END IF
C 210 CONTINUE
      RETURN
      END

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC


       subroutine vecrot_rotlat(IM,JM,TPH0D,TLM0D,VLAT,VLON,DRAY,CRAY)
                                                                                         
C  SUBPROGRAM: ROTLLE        ROTATE WINDS ON LAT/LONG GRID TO E-GRID
C   PRGMMR: T.BLACK         ORG: W/NP22     DATE: ??-??-??
C
C ABSTRACT: ROTATES WINDS ON THE LAT/LONG GRID TO THE ETA MODEL GRID
C               (same angles can be used to do reverse function)
C
C PROGRAM HISTORY LOG:
C   ??-??-??  T.BLACK
C   98-06-08  M.BALDWIN - CONVERT TO 2-D CODE
C   01-01-03  T BLACK - MODIFIED FOR MPI
C
C   INPUT ARGUMENT LIST:
C     VLAT     - LATITUDE OF E-GRID V POINTS (DEGREES)
C     VLON     - LONGITUDE OF E-GRID V POINTS (DEGREES)
C
C   OUTPUT ARGUMENT LIST:
C     DRAY     - ROTATION COSINE
C     CRAY     - ROTATION SINE
C
C ATTRIBUTES:
C   LANGUAGE: FORTRAN 90
C   MACHINE:  IBM RS/6000 SP
C
C***
C*** ROTATE THE LAT-LON WINDS TO/FROM THE E-GRID
C***
C
C    N O T E : INPUT LAT/LONG MUST BE CONVERTED TO RADIANS !!!
C
C----------------------------------------------------------------------
        parameter(D2RAD=1.745329E-2)
                                                                                         
                         D I M E N S I O N
     1 VLAT(IM,JM),VLON(IM,JM)
     2,CRAY(IM,JM),DRAY(IM,JM)
C----------------------------------------------------------------------
                                                                                         
        ERPHI0=TPH0D*D2RAD

        if (TLM0D .lt. 0) then
        ERLAM0=(TLM0D+360.)*D2RAD
        else
        ERLAM0=TLM0D*D2RAD
        endif
                                                                                         
      SPHI0 = SIN(ERPHI0)
      CPHI0 = COS(ERPHI0)
C
	
	write(6,*) 'TPH0D, TLM0D: ', TPH0D, TLM0D
	write(6,*) 'IM,JM: ', IM,JM

      DO J = 1, JM
      DO I = 1, IM

!mp DANGEROUS TEST (which had no impact)
!	if (VLON(I,J) .lt. 0) VLON(I,J)=VLON(I,J)+360.

        TLAT = VLAT(I,J) * D2RAD
        TLON = VLON(I,J) * D2RAD
        RELM = TLON - ERLAM0
        SRLM = SIN(RELM)
        CRLM = COS(RELM)
        SPH = SIN(TLAT)
        CPH = COS(TLAT)
        CC = CPH * CRLM
        TPH = ASIN(CPHI0 * SPH - SPHI0 * CC)
        RCTPH = 1.0 / COS(TPH)
        CRAY(I,J) = SPHI0 * SRLM * RCTPH
        DRAY(I,J) = (CPHI0 * CPH + SPHI0 * SPH * CRLM) * RCTPH

	if (J .eq. JM) then
	write(6,*) 'I,RELM,CRAY,DRAY: ',I,RELM,CRAY(I,J),DRAY(I,J)
	endif

      ENDDO
      ENDDO
	  end subroutine vecrot_rotlat

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCc

	subroutine vecrot_rotlat_new(IM,JM,TPH0D,TLM0D,VLAT,VLON,
     &                     COSALP,SINALP)

          
       REAL:: VLAT(IM,JM),VLON(IM,JM),COSALP(IM,JM),SINALP(IM,JM)


      PI =3.141592654
      D2R=PI/180.
      R2D=1./D2R
C-----------------------------------------------------------------------
C***
C***  READY TO COMPUTE ROTATION ANGLES.
C***
C***  FORMULAS USE GEODETIC LONGITUDES POSITIVE WEST
C***  SO NEGATE LONGITUDES FROM ETALL (VLON AND VLONI)
C***  AND THE CENTRAL LONGITUDES OF THE GRIDS
C***  SINCE THEY ARRIVE HERE AS NEGATIVE WEST.

!mp	is this sufficiently general?  assumption made that all will be west?

C***
      ERLAM0=-TLM0D*D2R
      ERPHI0=TPH0D*D2R
      ERL0_OUT=ERLAM0/D2R
      CPHI0_OUT=COS(ERPHI0)
      SPHI0_OUT=SIN(ERPHI0)
C***
C***  COMPUTE EARTH ROTATION ANGLES FOR WINDS RELATIVE TO BOTH
C***  THE INPUT AND OUTPUT GRIDS ON THE OUTPUT GRID POINTS.
C***
      DO J=1,JM
      DO I=1,IM
        X=CPHI0_OUT*COS(VLAT(I,J)*D2R)*COS((-VLON(I,J)-ERL0_OUT)*D2R)
     1   +SPHI0_OUT*SIN(VLAT(I,J)*D2R)
        Y=-COS(VLAT(I,J)*D2R)*SIN((-VLON(I,J)-ERL0_OUT)*D2R)
        TVLON=ATAN(Y/X)
        IF(X.LT.0.)TVLON=TVLON+PI
        ARG=SPHI0_OUT*SIN(TVLON)/COS(VLAT(I,J)*D2R)
        ARG=AMIN1(ARG,1.)
        ARG=AMAX1(ARG,-1.)
	ALPHA=ASIN(ARG)
	COSALP(I,J)=COS(ALPHA)
	SINALP(I,J)=SIN(ALPHA)
C
	if (J .eq. JM) then
	write(6,*) 'I,COSALP,SINALP: ',I,COSALP(I,J),SINALP(I,J)
	endif
      ENDDO
      ENDDO

	  end subroutine vecrot_rotlat_new
