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

      subroutine get_smf_1d(nk,cbase_m,ctop_m,itype_in
     1  ,height_laps_1d,pres_laps_1d,temp_laps_1d,
     1                          rlwc_laps,prob_laps,mode)

!     1991 Steve Albers (Smith-Feddes Code adapted to LAPS)

!     The array pres_laps_1d begins at 900mb and decreases at 50mb steps.

!     real rlwc_laps(nk),prob_laps(nk)
!     real height_laps_1d(nk),pres_laps_1d(nk),temp_laps_1d(nk)

!     INPUTS
      integer nk              ! Number of LAPS vertical levels
      real cbase_m            ! Cloud base (Meters MSL)
      real ctop_m             ! Cloud top (Meters MSL)
      integer itype_in        ! Cloud type (Now Hardwired to 1) (Stratus)
      real height_laps_1d(nk) ! Vert array, Heights (Meters MSL)
      real pres_laps_1d(nk)   ! Vert array, LAPS Pressure Levels (MB)
      real temp_laps_1d(nk)   ! Vert array, Ambient Temperature (Deg K)
      integer mode            ! Not used (Can be removed)


!     OUTPUTS
      real rlwc_laps(nk)      ! Vert array, LWC (G/M**3)
                                ! Should be assigned for all LAPS levels within
                                ! The cloud layer.
      real prob_laps(nk)      ! Vert array, Not used right now.

      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /IODEV/K5,K6,K7,K8
C
C     DIMENSION WATR(4,2,260),XMVD(4),ICLD(8),ICTP(4),ICP(4)
      DIMENSION WATR(4,2,260),XMVD(4),ICTP(4),ICP(4)
      CHARACTER*70 TITLE
      CHARACTER*2 CLDTYP(4)
C     CHARACTER*1 ANS
      CHARACTER*11 AINDEX(4,260),HINDEX(4,260)

      COMMON /PROF/ NOBHGT,HEIGHT(200),PRESSP(200),TEMPRO(200)
C
C-------Define the Input/Output Logical Devices
C
      K5 = 5
      K6 = 6
      K7 = 7
      K8 = 8

      bot(1) = cbase_m    ! Bottom Height m
      top(1) = ctop_m    ! Top Height m

      if(mode .eq. 1)then ! LAPS data is passed in through the call
          nobhgt = nk
          do i = 1,nobhgt
              height(i) = height_laps_1d(i)
              pressp(i) = pres_laps_1d(i)
              tempro(i) = temp_laps_1d(i)
          enddo ! i

      elseif(mode .eq. 2)then ! Test Mode
!         Read in the test sounding
          call RAOB(NXFLAG,NVFLAG,MAXOBSR) ! This routine was modified for PROFS

      endif

      ndecks = 1

!     0        CLDTYP(I) = '  '
!     1        CLDTYP(I) = 'St'
!     2        CLDTYP(I) = 'Sc'
!     3        CLDTYP(I) = 'Cu'
!     4        CLDTYP(I) = 'Ns'
!     5        CLDTYP(I) = 'Ac'
!     6        CLDTYP(I) = 'As'
!     7        CLDTYP(I) = 'Cs'
!     8        CLDTYP(I) = 'Ci'
!     9        CLDTYP(I) = 'Cc'
!    10        CLDTYP(I) = 'Cb'

      ictp(1) = itype_in     ! Type of layer 1
      icp(1)  = 100   ! Percent Cloud Cover of layer 1

      if(top(1) .le. bot(1))return

      call INTERP_smf(NDECKS)
      call LWC(NDECKS,ICP,ICTP,WATR,XMVD) ! Add this to link into .NEW routines

!     call a modified version of the output routine
c     call OUTPUT(CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE)

      call OUT_LAPS
     1  (nk,CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE,
     1         rlwc_laps,prob_laps)

!     write(6,*)' height  LWC  %Liquid'
      i = 1
      n1 = watr(i,1,1)
        iscript = 0
        do k = 3,iscript ! 260
            r_height = (k+n1-2)*100.
            if(watr(i,1,k) .gt. 0.)
     1  write(6,101)r_height,(watr(i,j,k),j=1,2)
101         format(1x,f8.1,'m ',2f7.3)
        enddo ! k

      return

      end
      SUBROUTINE RAOB(NXFLAG,NVFLAG,MAXOBSR)
C***********************************************************************
C                          SUBROUTINE RAOB
C***********************************************************************
C<Begin>
C<Identification>          Name:  RAOB
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  USRINT
C=======================================================================
C<Description>
C    This subroutine reads in the vertical profile or upperair data
C    needed to run the Smith-Feddes model.
C=======================================================================
C<Called routines>
C    HSORT - (subroutine) eliminates redundant height levels and sorts
C            the RAOB data in ascending order of height
C    INKEY - (subroutine) gets a character from the keyboard without
C            echo.  The routine is not available in source code form.
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL RAOB(NXFLAG,NVFLAG,MAXOBSR)
C    Input:
C       MAXOBSR  - (integer) the maximum number of upperair levels
C                  allowed
C    Output:
C       NXFLAG   - (integer) a flag
C                  0 - the option to compute the icing index was not
C                      selected
C                  1 - the option to compute the icing index was
C                      selected
C       NVFLAG   - (integer) a flag
C                  0 - the vertical profile data is not yet available
C                      for use
C                  1 - the vertical profile data is available for use
C    Output (in Common):
C       NOBHGT   - (integer) the number of levels in the temperature-
C                  pressure-height vertical profile
C       HEIGHT   - (real) the height of each level in the profile
C       PRESSP   - (real) the pressure profile
C       TEMPRO   - (real) the temperature profile
C    Common:
C       K5       - (integer) the keyboard input unit
C       K6       - (integer) the screen output unit
C       K7       - (integer) the disk file unit
C=======================================================================
C<History>
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Original code.
C    09/14/89  ASL    (505) 678-1570    Elton P. Avara
C              Created subroutine HSORT from part of this code
C    12/21/89  Cut out screen i/o for PROFS version, all this does now
C              is essentially reads the sounding from a file.
C=======================================================================
C<End>
C***********************************************************************
C
C-------Set the minimums and the maximums of the different input
C       parameters so there are no magic numbers floating around.
C
      real TEMPMIN,TEMPMAX,PRESMIN,PRESMAX,HEITMIN,HEITMAX
      PARAMETER (TEMPMIN  =   150.000)
      PARAMETER (TEMPMAX  =   349.999)
      PARAMETER (PRESMIN  =    50.000)
      PARAMETER (PRESMAX  =  1099.999)
      PARAMETER (HEITMIN  =    50.000)
      PARAMETER (HEITMAX  = 24999.999)
C
      COMMON /PROF/ NOBHGT,HEIGHT(200),PRESSP(200),TEMPRO(200)
      COMMON /IODEV/K5,K6,K7,K8
C
      CHARACTER*39 DATAFILE
C     CHARACTER*10 SCRCLR,ANSWER
C     CHARACTER*8 CURPOS
C     CHARACTER*4 SCRBLK,HUNITS(2),TUNITS(3),PUNITS(2)
      CHARACTER*4 HUNITS(2),TUNITS(3),PUNITS(2)
C     CHARACTER*3 LERASE
C     CHARACTER*1 NUM(10),ANS,ANS2
      CHARACTER*1 NUM(10)
C
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/
      DATA HUNITS/'(m) ','(ft)'/
      DATA TUNITS/'(K) ','(C) ','(F) '/
      DATA PUNITS/'(mb)','(in)'/
      DATA NHEIT,NTEMP,NPRES/3*1/
C
        DATAFILE='sounding.dat'
        OPEN(K7,FILE=DATAFILE,ERR=19,STATUS='OLD')
        READ(K7,*) NOBHGT
        DO 21 J=1,NOBHGT
   21     READ(K7,*) HEIGHT(J),TEMPRO(J),PRESSP(J)
        CLOSE(K7)
C
C-------Get the desired units for height, temperature, and pressure.
C
      NHEIT=1
      NTEMP=1
      NPRES=1
C
C-------If vertical profile data are available from a file, convert the
C       data to the desired units.
C
   50 HEITMN=HEITMIN
      HEITMX=HEITMAX
C
      TEMPMN=TEMPMIN
      TEMPMX=TEMPMAX
C
      PRESMN=PRESMIN
      PRESMX=PRESMAX
C
C-------Initialize the cursor position and the Icing Index computation
C       flag
C
      ITEM=0
      NXFLAG=0
C
C-------Set the vertical profile flag and return.
C
      NVFLAG=1
C
19    RETURN
      END
      SUBROUTINE INTERP_smf(NDECKS)
C***********************************************************************
C                          SUBROUTINE INTERP
C***********************************************************************
C<Begin>
C<Identification>          Name:  INTERP
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  INPUT
C=======================================================================
C<Description>
C    This subroutine takes the temperature (TEMPRO) and pressure
C    (PRESSP) profiles and interpolates between them to yield the
C    temperatures at the geometric layer bottoms (BTEM), 100 m height
C    levels (TEMP), and geometric layer tops (TTEM); and the pressures
C    at the geometric layer bottoms (PRB), 100 m height levels (PRES),
C    and geometric layer tops (PRT).
C=======================================================================
C<Called routines>
C    None
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL INTERP_smf(NDECKS)
C    Input:
C       NDECKS   - (integer) the number of cloud decks (1-4)
C    Input (in Common):
C       NOBHGT   - (integer) the number of levels in the temperature-
C                  pressure-height vertical profile
C       HEIGHT   - (real) the height of each level in the profile
C       PRESSP   - (real) the pressure profile
C       TEMPRO   - (real) the temperature profile
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C    Output (in Common):
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C    Common:
C       K6       - (integer) the screen output unit
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C    05/03/89  ASL    (505) 678-1570    Elton P. Avara
C              Minor modifications to the code and added the variable
C              PRES computations.
C    09/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definitions of PRES and TEMP to give values
C              every 100 m.  Also added NLVL, PRT, TTEM, and HEITL.
C              Major modifications.
C=======================================================================
C<End>
C***********************************************************************
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /PROF/ NOBHGT,HEIGHT(200),PRESSP(200),TEMPRO(200)
      COMMON /IODEV/K5,K6,K7,K8
C
C-------Initialize the height indices
C
      NLVL=0
      L=2
C
C-------Find the greatest height level required
C
      XMAXHT=0
      DO 5 J=1,NDECKS
    5 IF(XMAXHT.LT.TOP(J)) XMAXHT=TOP(J)
C
C-------The first height level is ground level
C
      NLVL=NLVL+1
      HEITL(NLVL)=HEIGHT(1)
      PRES(NLVL)=PRESSP(1)
      TEMP(NLVL)=TEMPRO(1)
C
C-------Find the next height above ground level which is a multiple of
C       100 m
C
      HEITZ=100*INT(0.01*HEIGHT(1))
      IF(HEITZ.LT.HEIGHT(1)) HEITZ=HEITZ+100.0
C
C-------Increment the height by 100 meters and get the temperature
C       and pressure
C
   10 NLVL=NLVL+1
      HEITZ=HEITZ+100.0
      HEITL(NLVL)=HEITZ
C
C-------Check the height level against the maximum height required
C
      IF(HEITL(NLVL)-XMAXHT.GE.100.0) THEN
        NLVL=NLVL-1
        GOTO 30
      ENDIF
C
C-------Loop over the temperature and pressure heights from the
C       radiosonde
C
      DO 20 J=L,NOBHGT
C
C-------Now as soon as the observation height reaches the desired
C       height level, interpolate the temperature to that height.
C
      IF (HEIGHT(J).GE.HEITL(NLVL)) THEN
        SLOPE = (TEMPRO(J)-TEMPRO(J-1))/(HEIGHT(J)-HEIGHT(J-1))
        TRCPT = -SLOPE * HEIGHT(J) + TEMPRO(J)
        TEMP(NLVL) = SLOPE * HEITL(NLVL) + TRCPT
C
C-------Now interpolate to find the pressure
C
C       IF(ABS(SLOPE).GT.1.0E-4) THEN
C
C--------------Linear trend to temperature in vertical.  The number
C              0.034109527 is the acceleration due to gravity divided
C              by the gas constant
C
C         RATIO = ((TEMPRO(J-1) + SLOPE * (HEITL(NLVL)-HEIGHT(J-1))) /
C    1    TEMPRO(J-1))**(-0.034109527/SLOPE)
C       ELSE
C
C--------------Use average layer temperature-AVTEMP-to calculate pressure
C
C
          AVTEMP = (TEMPRO(J-1) + TEMPRO(J))/2.
          RATIO =EXP(-0.034109527*(HEITL(NLVL)-HEIGHT(J-1))/AVTEMP)
C       ENDIF
        PRES(NLVL) = PRESSP(J-1) * RATIO
C
C-------Go get another height level
C
        L=J
        GOTO 10
C
C-------Have we run out of upperair data?
C
      ELSEIF (J.EQ.NOBHGT) THEN
        WRITE(K6,*)'Temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        WRITE(K6,*)'The model needs this information to perform the calc
     1ulations.'
C
        STOP
      ENDIF
   20 CONTINUE
C
C-------Loop over the geometric layers and get the layer bottom values
C       of temperature and pressure.
C
   30 L=2
      DO 40 I=1,NDECKS
C
C-------Loop over the temperature and pressure heights from the
C       radiosonde
C
      DO 50 J=L,NOBHGT
C
C-------Now as soon as the observation height reaches the desired
C       height level, interpolate the temperature to that height.
C
      IF (HEIGHT(J).GE.BOT(I)) THEN
        SLOPE = (TEMPRO(J)-TEMPRO(J-1))/(HEIGHT(J)-HEIGHT(J-1))
        TRCPT = -SLOPE * HEIGHT(J) + TEMPRO(J)
        BTEM(I) = SLOPE * BOT(I) + TRCPT
C
C-------Now interpolate to find the pressure
C
C       IF(ABS(SLOPE).GT.1.0E-4) THEN
C
C--------------Linear trend to temperature in vertical.  The number
C              0.034109527 is the acceleration due to gravity divided
C              by the gas constant
C
C         RATIO = ((TEMPRO(J-1) + SLOPE * (BOT(I)-HEIGHT(J-1))) /
C    1    TEMPRO(J-1))**(-0.034109527/SLOPE)
C       ELSE
C
C--------------Use average layer temperature-AVTEMP-to calculate pressure
C
          AVTEMP = (TEMPRO(J-1) + TEMPRO(J))/2.
          RATIO =EXP(-0.034109527*(BOT(I)-HEIGHT(J-1))/AVTEMP)
C       ENDIF
        PRB(I) = PRESSP(J-1) * RATIO
C
C-------Go get another geometric layer
C
        L=J
        GOTO 40
C
C-------Have we run out of upperair data?
C
      ELSEIF (J.EQ.NOBHGT) THEN
        WRITE(K6,*)'Temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        WRITE(K6,*)'The model needs this information to perform the calc
     1ulations.'
C
        STOP
      ENDIF
   50 CONTINUE
C
   40 CONTINUE
C
C-------Loop over the geometric layers and get the layer top values
C       of temperature and pressure.
C
      L=2
      DO 60 I=1,NDECKS
C
C-------Loop over the temperature and pressure heights from the
C       radiosonde
C
      DO 70 J=L,NOBHGT
C
C-------Now as soon as the observation height reaches the desired
C       height level, interpolate the temperature to that height.
C
      IF (HEIGHT(J).GE.TOP(I)) THEN
        SLOPE = (TEMPRO(J)-TEMPRO(J-1))/(HEIGHT(J)-HEIGHT(J-1))
        TRCPT = -SLOPE * HEIGHT(J) + TEMPRO(J)
        TTEM(I) = SLOPE * TOP(I) + TRCPT
C
C-------Now interpolate to find the pressure
C
C       IF(ABS(SLOPE).GT.1.0E-4) THEN
C
C--------------Linear trend to temperature in vertical.  The number
C              0.034109527 is the acceleration due to gravity divided
C              by the gas constant
C
C         RATIO = ((TEMPRO(J-1) + SLOPE * (TOP(I)-HEIGHT(J-1))) /
C    1    TEMPRO(J-1))**(-0.034109527/SLOPE)
C       ELSE
C
C--------------Use average temperature in layer-AVTEMP-to calculate pressure
C
          AVTEMP = (TEMPRO(J-1) + TEMPRO(J))/2.
          RATIO =EXP(-0.034109527*(TOP(I)-HEIGHT(J-1))/AVTEMP)
C       ENDIF
        PRT(I) = PRESSP(J-1) * RATIO
C
C-------Go get another geometric layer
C
        L=J
        GOTO 60
C
C-------Have we run out of upperair data?
C
      ELSEIF (J.EQ.NOBHGT) THEN
        WRITE(K6,*)'Temperature and pressure profiles do not extend 100
     1m above highest cloud top.'
        WRITE(K6,*)'The model needs this information to perform the calc
     1ulations.'
C
        STOP
      ENDIF
   70 CONTINUE
C
   60 CONTINUE
C
      RETURN
      END


      SUBROUTINE AL(LVLCT,LAYR,TI,N1,N2,CALW)
C***********************************************************************
C                          SUBROUTINE AL
C***********************************************************************
C<Begin>
C<Identification>          Name:  AL
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  LW
C=======================================================================
C<Description>
C    This subroutine computes moist adiabatic Liquid Water Content of
C    non-cirriform cloud types
C=======================================================================
C<Called routines>
C    TZ     - (subroutine) calculates moist adiabatic lapse rate of
C             temperature and saturation vapor pressure
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL AL(LVLCT,LAYR,TI,N1,N2,CALW)
C    Input:
C       LVLCT    - (integer) cloud type of the cloud deck
C       LAYR     - (integer) cloud layer
C       TI       - (real) temperature at bottom of cloud layer (C)
C       N1       - (integer) the lowest 100 m height level within the
C                  cloud layer
C       N2       - (integer) the highest 100 m height level within the
C                  cloud layer
C    Output:
C       CALW     - (real) array of adiabatic Liquid Water Content at
C                  each 100 m within the cloud layer
C    Common:
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code.
C    05/03/89  ASL    (505) 678-1570    Elton P. Avara
C              Minor modifications to the code and added the variable
C              PRES for calculating DPZ.
C    09/12/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definition of PRES to give values every
C              100 m.  Also added NLVL, TOP, BTEM, TTEM, PRT, TEMP,
C              and HEITL.  Major modifications.
C=======================================================================
C<End>
C***********************************************************************
C
C  Author - C. William Rogers
C
C  PROCEDURE:
C       Compute moist adiabatic temperature lapse rate in S/R TZ.  Then
C  compute water condensed out in air parcel which rises 10M and cools
C  at lapse rate from S/R TZ.  LWC is difference between saturated
C  conditions at bottom and top of 10M thick layer.  This procedure
C  is repeated in 10M steps or fraction thereof until LWC at each 100 m
C  in cloud is obtained.  All sublayer variables are for 10M increments
C  except at top of layer or bottom of next higher layer.  For
C  cumuliform cloud reduce LWC to account for entrainment.
C
C  REFERENCE:
C       Rogers, C.W., J.T. Hanley and E.J. Mack, 1985: "Updating the
C       Smith-Feddes Model", Final Report Contract No. N00228-84-C-3157
C       Calspan Report No. 7330-1.  Calspan Corp., P.O. Box 400 Buffalo,
C       New York 14225. (See Section 3.5)
C
C  NOTE:  All computations are carried out for 10M thick layers except
C         for shallower layers ending at bottom of next higher layer.
C
C  LOCAL GLOSSARY:
C       T      = bottom temperature of 10M layer (C) as input to S/R TZ
C                and in computing TN for ES2 computation.  It will be
C                (K) when returned by S/R TZ and in computation of RHO.
C       ALW    = adiabatic Liquid Water Content in cloud layer
C                bottom to bottom of next layer and at 10M intervals in
C                between
C       DP     = pressure change of 10M or less interval
C       DTZ    = moist adiabatic lapse rate of temperature
C       DZ     = 10M interval or less, the latter at BOT(LAYR) and
C                TOP(LAYR)
C       ES     = saturation vapor pressure at bottom of 10M layer
C       ES2    = saturation vapor pressure at top of 10M or less layer
C       HT     = height of height level in cloud layer relative to cloud
C                deck base
C       IGO    = control variable
C                if IGO = 1 - have reached top of layer
C                       = 0 - have not reached top of layer
C       P      = pressure at bottom of 10M layer
C       RHO    = absolute humidity at bottom of 10M layer
C       TN     = top temperature of 10M layer (C) in ES2 computation
C                and (K) in RHO2 computations
C       RHO2   = absolute humidity at top of 10M or less layer
C       RS     = universal gas constant
C       WR     = gas constant for water vapor
C       Y      = fraction of adiabatic LWC produced by entrainment
C       Z      = height parameter within cloud layer used in S/R AL
C
C-----------------------------------------------------------------------
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
C
      DIMENSION CALW(260)
      LOGICAL NSKIP
C
      DATA A1,B1,C1,A2,B2,C2/8.4897, -13.2191, 4.7295, 10.357,
     1-28.2416, 8.8846 /
      DATA RS/8.313E07/
      DATA WR/2.16528E-7/
C
C-------Initialize some variables
C
      ALW = 0.0
      IGO=0
      INC=0
      T=TI
C
      CALW(1)=ALW
      NSKIP=.FALSE.
      IF(N1.EQ.0.OR.N2.EQ.0) NSKIP=.TRUE.
C
C-------Begin computing the LWC each 100 m within the cloud.
C
      deltaz = 100.
      NZ=N1-1
      DO 100 I=1,2600
C
C-------Compute the LWC in 10 m (or less) height layers
C
      DZ=deltaz
C
C-------Compute the RAOB pressure and temperature gradients and height
C       increment for each 100 m height level.
C
      IF(INC.EQ.0) THEN
        NZ=NZ+1
        LVL=NZ-N1+2
        IF(NSKIP) THEN
          P=PRB(LAYR)
          Z=BOT(LAYR)
CC        T=BTEM(LAYR)-273.16
          DELZ=TOP(LAYR)-BOT(LAYR)
          DZ=DELZ-deltaz*AINT(DELZ/deltaz)
          DPZ=(PRT(LAYR)-PRB(LAYR))/DELZ
          DTZ=(TTEM(LAYR)-BTEM(LAYR))/DELZ
          INC=1
        ELSEIF(NZ.EQ.N1) THEN
          P=PRB(LAYR)
          Z=BOT(LAYR)
CC        T=BTEM(LAYR)-273.16
          DELZ=HEITL(N1)-BOT(LAYR)
          DZ=DELZ-deltaz*AINT(DELZ/deltaz)
          DPZ=(PRES(N1)-PRB(LAYR))/DELZ
          DTZ=(TEMP(N1)-BTEM(LAYR))/DELZ
          INC=1
        ELSEIF(NZ.GT.N1.AND.NZ.LE.N2) THEN
          P=PRES(NZ-1)
          Z=HEITL(NZ-1)
CC        T=TEMP(NZ-1)-273.16
          DPZ=(PRES(NZ)-PRES(NZ-1))/(HEITL(NZ)-HEITL(NZ-1))
          DTZ=(TEMP(NZ)-TEMP(NZ-1))/(HEITL(NZ)-HEITL(NZ-1))
          INC=1
        ELSEIF(NZ.GT.N2) THEN
          P=PRES(N2-1)
          Z=HEITL(N2-1)
CC        T=TEMP(N2-1)-273.16

          DPZ=(PRT(LAYR)-PRES(N2))/(TOP(LAYR)-HEITL(N2))
          DTZ=(TTEM(LAYR)-TEMP(N2))/(TOP(LAYR)-HEITL(N2))
          INC=1
        ENDIF
      ENDIF
C
C-------If at the top of the cloud layer, change the index to show it.
C
      IF(TOP(LAYR)-Z.LT.deltaz) THEN
        DZ=TOP(LAYR)-Z
        IGO=1
      ENDIF
C
C-------Compute the moist adiabatic temperature lapse rate and
C       saturation vapor pressure
C
!     Code for TZ has been pulled in. Note that DTZ was changed to DMTZ.
!     The value of T is also modified in that subroutine.
!     CALL TZ(P,T,DMTZ,ES)
!     SUBROUTINE TZ(P,T,DTZ,ES)

      G=980.
      CP=.239
      R=.06855
C
C  COMPUTE SATURATION VAPOR PRESSURE AT BOTTOM OF 10M SUBLAYER
C
      ES=6.112*EXP(17.67*T/(T+243.5))
      XL=595.-0.5*T
C
C  COMPUTE DIFFERENTIAL OF ES WITH CELSIUS TEMPERATURE
C
      DES=ES*((17.67/(T+243.5))-(17.67*T/(T+243.5)**2))
      TOLD = T
      T=T+273.16
C
C  COMPUTE MOIST ADIABATIC LAPSE RATE OF TEMPERATURE
C
      DMTZ=-G*((1.0+.621*ES*XL/(P*R*T))/(CP+.621*XL*DES/P))
      DMTZ=DMTZ*2.39E-6


C
C-------Compute the mixing ratio and absolute humidity at the bottom of
C       the layer
C
      AMR1=0.622*ES/(P-ES)
      RHO=1000.*P/(2.8704*T)
C
C-------Get the height, pressure, and moist adiabatic temperature (C) of
C       the top of the 10 m (or less) layer
C
      Z=Z+DZ
      P=P+DPZ*DZ
      TN=(T-273.16)+DMTZ*DZ
CC    TN=(T-273.16)+DTZ*DZ
C
C-------Compute the and saturation vapor pressure, mixing ratio, and
C       absolute humidity at the top of the layer
C
C-------Try an approximate form for the saturation vapor pressure
C       at the higher level
      DTEMP = TN - TOLD
      ES2 = ES + DES*DTEMP
C     ES2=6.112*EXP(17.67*TN/(TN+243.5))
      TN=TN+273.16
      AMR2=0.622*ES2/(P-ES2)
      RHO2=1000.*P/(2.8704*TN)
C
C-------Compute the increase in the adiabatic liquid water over the 10 m
C       ascent and add to the total liquid water content thus far.
C
      ALW=ALW+(AMR1-AMR2)*((RHO+RHO2)*0.5)
      T=TN-273.16
      CALW(LVL)=ALW
C
C-------If at the top of a 100 m height level, set the index to show it.
C
      IF(IGO.NE.1) THEN
        IF(NSKIP) THEN
          GOTO 100
        ELSE
          IF(ABS(Z-HEITL(NZ)).GT.1.0) GOTO 100
        ENDIF
      ENDIF
      INC=0
C                                           !MODIFICATIONS NOV/87
C-------Time to reduce the CALW due to entrainment.
C       There are two methods for reducing the CALW.  The first is via
C       Skatskii's curve which is used on the cumuliform clouds
C       (Cb, Sc, Cu, and Ac) and on the unknown cloud type.
C       The second is via Warner's curve which is used on the stratiform
C       clouds (St, As, and Ns).
C
C
C-------Calculate the thickness of the cloud so far in kilometers.
C
      HT=(Z-BOT(LAYR))*.001
C
      IF((LVLCT.EQ.10).OR.(LVLCT.EQ.2).OR.(LVLCT.EQ.3).OR.(LVLCT.EQ.5)
     1.OR.(LVLCT.EQ.25)) THEN
C
C-------Reduce CALW due to entrainment via Skatskii's curve
C                                        !MODIFICATIONS OCT/86
        IF(HT.LT.0.3) THEN
          Y = -1.667 * (HT - 0.6)
        ELSEIF (HT.LT.1.0) THEN
          ARG1 = B1 * B1 - 4.0 * A1 * (C1 - HT)
          Y = (-B1 - SQRT(ARG1)) / (2.0 * A1)
        ELSEIF (HT.LT.2.9) THEN
          ARG2 = B2 * B2 - 4.0 * A2 * (C2 - HT)
          Y = (-B2 - SQRT(ARG2)) / (2.0 * A2)
        ELSE
          Y = 0.26
        ENDIF

C                                        !END OF MODIFICATIONS OCT/86
      ELSEIF ((LVLCT.EQ.1).OR.(LVLCT.EQ.4).OR.(LVLCT.EQ.6)) THEN
C
C-------Reduce CALW due to entrainment via Warner's curve.
C
        IF(HT.LE.0.032) THEN
          Y = -11.0 * HT + 1.0
        ELSEIF (HT.LE.0.177) THEN
          Y =  -1.4 * HT + 0.6915
        ELSEIF (HT.LE.0.726) THEN
          Y =  -0.356 * HT + 0.505
        ELSEIF (HT.LE.1.5) THEN
          Y =  -0.0608 * HT + 0.2912
        ELSE
          Y =   0.20
        ENDIF
C
      ELSE
C
C-------Do nothing because cirriform clouds or clear and then we
C       shouldn't even be in this subroutine.
C
        Y = 1.0
C
      ENDIF
C                                        !END OF MODIFICATIONS NOV/87
C
C-------Now calculate the reduced CALW.
C
      CALW(LVL) = CALW(LVL) * Y
C
C-------If at the top of the cloud layer, exit the loop.
C
      IF(IGO .EQ. 1) GOTO 150
C
  100 CONTINUE
C
  150 RETURN
      END
      SUBROUTINE TZ(P,T,DTZ,ES)
C***********************************************************************
C                          SUBROUTINE TZ
C***********************************************************************
C<Begin>
C<Identification>          Name:  TZ
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  AL
C=======================================================================
C<Description>
C    Compute moist adiabatic lapse rate of temperature and saturation
C    vapor pressure at bottom of 10M layer.
C=======================================================================
C<Called routines>
C    None
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL TZ(P,T,DTZ,ES)
C    Input:
C       P        - (real) pressure at bottom of 10M layer
C       T        - (real) bottom temperature of 10M layer (C)
C    Output:
C       DTZ      - (real) moist adiabatic lapse rate based on variables
C                  at bottom of 10M layer
C       ES       - (real) vapor pressure at bottom of 10M layer
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C=======================================================================
C<End>
C***********************************************************************
C
C  Author - C. William Rogers
C
C  LOCAL GLOSSARY:
C       CP      = specific heat of air
C       G       = acceleration of gravity
C       XL      = latent heat of vaporization
C       R       = gas constant
C       DES     = differential of ES with temperature (C)
C
C  REFERENCE:
C       Rogers, C.W., J.T. Hanley and E.J. Mack, 1985: "Updating the
C       Smith Feddes Model", Final Report Contract No. N00228-84-C-3157
C       Calspan Report No. 7330-1.  Calspan Corp., P.O. Box 400 Buffalo,
C       New York 14225. (See Section 3.5)
C
C-----------------------------------------------------------------------
C
      G=980.
      CP=.239
      R=.06855
C
C  COMPUTE SATURATION VAPOR PRESSURE AT BOTTOM OF 10M SUBLAYER
C
      ES=6.112*EXP(17.67*T/(T+243.5))
      XL=595.-0.5*T
C
C  COMPUTE DIFFERENTIAL OF ES WITH CELSIUS TEMPERATURE
C
      DES=ES*((17.67/(T+243.5))-(17.67*T/(T+243.5)**2))
      T=T+273.16
C
C  COMPUTE MOIST ADIABATIC LAPSE RATE OF TEMPERATURE
C
      DTZ=-G*((1.0+.621*ES*XL/(P*R*T))/(CP+.621*XL*DES/P))
      DTZ=DTZ*2.39E-6
C
      RETURN
      END
      SUBROUTINE MVD(NDECKS,ICTP,XMVD)
C***********************************************************************
C                          SUBROUTINE MVD
C***********************************************************************
C<Begin>
C<Identification>          Name:  MVD
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  SMFD
C=======================================================================
C<Description>
C    Determines the Mean Volume Diameter of the droplets within the
C    cloud decks
C=======================================================================
C<Called routines>
C    None
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL MVD(NDECKS,ICTP,XMVD)
C    Input:
C       NDECKS   - (integer) the number of cloud decks
C       ICTP     - (integer) code for the cloud type (up to 4 cloud
C                  decks)
C    Output:
C       XMVD     - (real) Mean Volume Diameter for each cloud
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code as part of main SMFD routine
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C    09/14/89  ASL    (505) 678-1570    Elton P. Avara
C              Created a separate subroutine for XMVD determination
C=======================================================================
C<End>
C***********************************************************************
C
      DIMENSION ICTP(4),XMVD(4)
C
C-------The MVDs will be defined to a default representative value for
C       each of the different cloud types
C
      DO 10 I = 1,NDECKS
C                                                 ! no cloud
      IF(ICTP(I).EQ.0) THEN
        XMVD(I) = 0.0
C                                                 ! St
      ELSEIF (ICTP(I).EQ. 1) THEN
        XMVD(I) = 12.0
C                                                 ! Sc
      ELSEIF (ICTP(I).EQ. 2) THEN
        XMVD(I) = 10.0
C                                                 ! Cu
      ELSEIF (ICTP(I).EQ. 3) THEN
        XMVD(I) = 18.0
C                                                 ! Ns
      ELSEIF (ICTP(I).EQ. 4) THEN
        XMVD(I) = 12.0
C                                                 ! Ac
      ELSEIF (ICTP(I).EQ. 5) THEN
        XMVD(I) = 18.0
C                                                 ! As
      ELSEIF (ICTP(I).EQ. 6) THEN
        XMVD(I) = 10.0
C                                                 ! Cb
      ELSEIF (ICTP(I).EQ.10) THEN
        XMVD(I) = 25.0
C                                                 ! cirriform
      ELSE
        XMVD(I) = 10.0
C
      ENDIF
   10 CONTINUE
C
      RETURN
      END
      SUBROUTINE IINDEX(XMVD,WATR,NDECKS,AINDEX,HINDEX)
C***********************************************************************
C                          SUBROUTINE IINDEX
C***********************************************************************
C<Begin>
C<Identification>          Name:  IINDEX
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  SMFD
C=======================================================================
C<Description>
C    This subroutine will determine the icing severity index for both
C    fixed-wing aircraft (AINDEX) and rotary-wing helicopters (HINDEX),
C    given the XMVD, WATR, and height level temperatures TEMP, BTEM, and
C    TTEM.
C=======================================================================
C<Called routines>
C    None
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL IINDEX(XMVD,WATR,NDECKS,AINDEX,HINDEX)
C    Input:
C       XMVD     - (real) Mean Volume Diameter for each cloud
C       WATR     - (real) output array containing LWC every 100 m
C       NDECKS   - (integer) the number of cloud decks
C    Output:
C       AINDEX   - (character) Icing Severity Index for fixed wing
C                  aircraft every 100 m within the cloud
C       HINDEX   - (character) Icing Severity Index for helicopters
C                  every 100 m within the cloud
C    Common:
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C       K6       - (integer) the screen output unit
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C    05/03/89  ASL    (505) 678-1570    Elton P. Avara
C              Corrected a bug in determining severe LWC and minor
C              modifications to other code.
C    09/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definition of TEMP, WATR, AINDEX and HINDEX
C              to give values every 100 m.  Also added NLVL.
C=======================================================================
C<End>
C***********************************************************************
C
      INTEGER NDIAMS
      PARAMETER   (NDIAMS = 6)

      REAL FREEZE,TOOCLD
      PARAMETER   (FREEZE = 273.15)
      PARAMETER   (TOOCLD = 243.15)
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /IODEV/K5,K6,K7,K8
C
      CHARACTER*11 AINDEX(4,260),HINDEX(4,260)
      DIMENSION XMVD(4),WATR(4,2,260),DIAM(NDIAMS),XLWCD(NDIAMS)
C
      DATA DIAM/10.0,15.0,20.0,30.0,40.0,50.0/
      DATA XLWCD/2.50,1.30,0.85,0.65,0.55,0.50/
C
C-------Loop over the number of geometric layers which have
C       clouds
C
      DO 300 I = 1,NDECKS
C
C-------First thing we need to do is to determine the threshold of the
C       severe liquid water content (XLWCS)
C
      IF(DIAM(1).LE.XMVD(I).AND.XMVD(I).LE.DIAM(NDIAMS)) THEN
C
        DO 10 K=2,NDIAMS
C
        IF(XMVD(I).GE.DIAM(K-1).AND.XMVD(I).LE.DIAM(K)) THEN
          SLOPE = (XLWCD(K) - XLWCD(K-1))/(DIAM(K) - DIAM(K-1))
          TRCPT = -SLOPE * DIAM(K) + XLWCD(K)
          XLWCS = SLOPE * XMVD(I) + TRCPT
        ENDIF
C
   10   CONTINUE
C
      ELSE
C
        WRITE(K6,1111)'The calculated MVD = ',XMVD(I),' microns.'
! Hongli Jiang: change f3.1 to F4.1. 11/27/2013
 1111   FORMAT(1X,A,F4.1,A)
        WRITE(K6,*)'This is not within the  possible values of 10-50 mic
     1rons.  Thus, the'
        WRITE(K6,*)'icing severity index is not available for this case.
     1'
C
        RETURN
C
      ENDIF
C
C-------Now we can determine the thresholds for the moderate and light
C       liquid water contents, XLWCM and XLWCL.
C
      XLWCM = 0.5 * XLWCS
      XLWCL = 0.1 * XLWCS
C
C-------Compute the Icing Severity Index for each height level within
C       the cloud deck
C
      N1=WATR(I,1,1)+0.001
      N2=WATR(I,1,2)+0.001
C
      IF(N1.LT.0.OR.N2.LT.0) THEN
        AINDEX(I,1) = '    missing'
        HINDEX(I,1) = '    missing'
        GOTO 300
      ELSEIF(N1.EQ.0.OR.N2.EQ.0) THEN
        MAXCNT=2
      ELSE
        MAXCNT=N2-N1+3
      ENDIF
C
      DO 200 L=1,MAXCNT
C
C-------Get the temperature for the height level
C
      IF(L.EQ.1) THEN
        T=BTEM(I)
      ELSEIF(L.EQ.MAXCNT) THEN
        T=TTEM(I)
      ELSE
        T=TEMP(N1+L-2)
      ENDIF
C
C-------Now to determine the actual icing severity index for both the
C       fixed wing aircraft (AINDEX) and the rotary wing helicopters
C       (HINDEX), if the temperature is below freezing.
C
      IF((TOOCLD.LT.T).AND.(T.LE.FREEZE)) THEN
C
C-------First, the fixed wing (AINDEX)
C
        IF(WATR(I,1,L+2).LT.XLWCL) THEN
          AINDEX(I,L) = '  A1: TRACE'
        ELSEIF (WATR(I,1,L+2).LT.XLWCM) THEN
          AINDEX(I,L) = '  A2: LIGHT'
        ELSEIF (WATR(I,1,L+2).LT.XLWCS) THEN
C
          IF(T.LE.FREEZE - 5.0) THEN
            AINDEX(I,L) = 'A3:MODERATE'
          ELSE
            AINDEX(I,L) = '  A2: LIGHT'
          ENDIF
C
        ELSE
C
          IF(T.LE.FREEZE - 5.0) THEN
            AINDEX(I,L) = ' A4: SEVERE'
          ELSEIF (T.LE.FREEZE - 3.0) THEN
            AINDEX(I,L) = 'A3:MODERATE'
          ELSE
            AINDEX(I,L) = '  A2: LIGHT'
          ENDIF
C
        ENDIF
C
C-------Now for the rotary wing (HINDEX)
C
        IF(WATR(I,1,L+2).LT.XLWCL) THEN
          HINDEX(I,L) = '  H1: TRACE'
        ELSEIF (WATR(I,1,L+2).LT.XLWCM) THEN
C
          IF(T.LE.FREEZE - 5.0) THEN
            HINDEX(I,L) = '  H2: LIGHT'
          ELSE
            HINDEX(I,L) = '  H1: TRACE'
          ENDIF
C
        ELSEIF (WATR(I,1,L+2).LT.XLWCS) THEN
C
          IF(T.LE.FREEZE - 10.0) THEN
            HINDEX(I,L) = 'H3:MODERATE'
          ELSEIF (T.LE.FREEZE - 5.0) THEN
            HINDEX(I,L) = '  H2: LIGHT'
          ELSE
            HINDEX(I,L) = '  H1: TRACE'
          ENDIF
C
        ELSE
C
          IF(T.LE.FREEZE - 10.0) THEN
            HINDEX(I,L) = ' H4: SEVERE'
          ELSEIF (T.LE.FREEZE - 5.0) THEN
            HINDEX(I,L) = 'H3:MODERATE'
          ELSE
            HINDEX(I,L) = '  H2: LIGHT'
          ENDIF
C
        ENDIF
C
      ELSE
C
C-------The temperature is too warm or too cold to produce freezing
C
        AINDEX(I,L) = 'A0:NO ICING'
        HINDEX(I,L) = 'H0:NO ICING'
C
      ENDIF
C
  200 CONTINUE
C
  300 CONTINUE
C
      RETURN
      END
      SUBROUTINE CLTYPE(ICTP,NDECKS,CLDTYP)
C***********************************************************************
C                          SUBROUTINE CLTYPE
C***********************************************************************
C<Begin>
C<Identification>          Name:  CLTYPE
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  SMFD
C=======================================================================
C<Description>
C    This subroutine will determine the two character abbreviation for
C    the cloud type from the 3DNEPH integer code for the cloud type.
C=======================================================================
C<Called routines>
C    None
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL CLTYPE(ICTP,NDECKS,CLDTYP)
C    Input:
C       ICTP     - (integer) code for the cloud type (up to 4 cloud
C                  decks)
C       NDECKS   - (integer) the number of cloud decks
C    Output:
C       CLDTYP   - (character) cloud type in two ASCII characters
C    Common:
C       K6       - (integer) the screen output unit
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C=======================================================================
C<End>
C***********************************************************************
C
      DIMENSION ICTP(4)
      CHARACTER*2 CLDTYP(4)
C
C-------Loop over the number of cloud decks
C
      DO 100 I = 1,NDECKS
C
      IF(ICTP(I).EQ.0) THEN
        CLDTYP(I) = '  '
      ELSEIF (ICTP(I).EQ. 1) THEN
        CLDTYP(I) = 'St'
      ELSEIF (ICTP(I).EQ. 2) THEN
        CLDTYP(I) = 'Sc'
      ELSEIF (ICTP(I).EQ. 3) THEN
        CLDTYP(I) = 'Cu'
      ELSEIF (ICTP(I).EQ. 4) THEN
        CLDTYP(I) = 'Ns'
      ELSEIF (ICTP(I).EQ. 5) THEN
        CLDTYP(I) = 'Ac'
      ELSEIF (ICTP(I).EQ. 6) THEN
        CLDTYP(I) = 'As'
      ELSEIF (ICTP(I).EQ. 7) THEN
        CLDTYP(I) = 'Cs'
      ELSEIF (ICTP(I).EQ. 8) THEN
        CLDTYP(I) = 'Ci'
      ELSEIF (ICTP(I).EQ. 9) THEN
        CLDTYP(I) = 'Cc'
      ELSEIF (ICTP(I).EQ.10) THEN
        CLDTYP(I) = 'Cb'
      ELSEIF (ICTP(I).EQ.25) THEN
        CLDTYP(I) = '??'
      ELSE
        WRITE(6,*)'SOMETHING IS WRONG WITH THE CLOUD TYPE ARRAY,'
        WRITE(6,1001) I,ICTP(I)
 1001   FORMAT(1X,'ICTP ARRAY.  ICTP(',I1,') = ',I2)
C
        STOP
      ENDIF
C
  100 CONTINUE
C
      RETURN
      END
      SUBROUTINE OUTPUT(CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE)
C***********************************************************************
C                          SUBROUTINE OUTPUT
C***********************************************************************
C<Begin>
C<Identification>          Name:  OUTPUT
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  SMFD
C=======================================================================
C<Description>
C    This subroutine takes the output parameters and writes them out to
C    the screen or printer.
C=======================================================================
C<Called routines>
C    INKEY - (subroutine) gets a character from the keyboard without
C            echo.  The routine is not available in source code form.
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL OUTPUT(CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE)
C    Input:
C       CLDTYP   - (character) cloud type in two ASCII characters
C       WATR     - (real) output array containing LWC and percent liquid
C                  water every 100 m in the cloud
C       XMVD     - (real) Mean Volume Diameter for each cloud
C       AINDEX   - (character) Icing Severity Index for fixed wing
C                  aircraft every 100 m
C       HINDEX   - (character) Icing Severity Index for helicopters
C                  every 100 m
C       ICP      - (integer) percent (1-100%) cloud coverage for up to 4
C                  cloud decks
C       NDECKS   - (integer) the number of cloud decks
C       TITLE    - (character) the title to be displayed in the output
C    Output:
C       None
C    Common:
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C       K5       - (integer) the keyboard input unit
C       K6       - (integer) the screen output unit
C       K7       - (integer) the disk file unit
C       K8       - (integer) the printer output unit
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Major modifications to the code.  Completely rewritten.
C    05/03/89  ASL    (505) 678-1570    Elton P. Avara
C              Minor modifications to the code.
C    09/14/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definition of TEMP, WATR, AINDEX and HINDEX
C              to give values every 100 m.  Also added NLVL, HEITL,
C              PRES, PRT, and TTEM.  Major modifications.
C    12/21/89  Steve Albers - Modified for PROFS
C=======================================================================
C<End>
C***********************************************************************
C
      REAL       FREEZE
      PARAMETER   (FREEZE = 273.15)
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /IODEV/K5,K6,K7,K8
C
      DIMENSION ICP(4),WATR(4,2,260),XMVD(4)
C
      CHARACTER*70 TITLE
C     CHARACTER*14 FILENM
      CHARACTER*11 AINDEX(4,260),HINDEX(4,260)
C     CHARACTER*10 SCRCLR
C     CHARACTER*8 CURPOS
C     CHARACTER*4 SCRBLK
C     CHARACTER*3 LERASE
      CHARACTER*2 CLDTYP(4)
C     CHARACTER*1 NUM(10),ANS
      CHARACTER*1 NUM(10)
C
C     LOGICAL NPRT,NFIL
C
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/

    1 FORMAT(A,$)
C
C-------Now to write the cloud data starting with the lowest cloud
C       to the highest cloud.
C
C-------First write the cloud deck data that is not height conditional
C       to the screen (and printer/disk file).
C
      WRITE(K6,26) (INT(TOP(J)),J=1,NDECKS)
   26 FORMAT(' Cloud Top (m):',16X,4I12)
C
      WRITE(K6,24) (INT(BOT(J)),J=1,NDECKS)
   24 FORMAT(' Cloud Base (m):',15X,4I12)
C
      WRITE(K6,30) (CLDTYP(J),J=1,NDECKS)
   30 FORMAT(' Cloud Type:',19X,4(10X,A2))
C
      WRITE(K6,32) (XMVD(J),J=1,NDECKS)
   32 FORMAT(' Drop MVD (microns):',11X,4F12.1)
C
      WRITE(K6,34) (ICP(J),J=1,NDECKS)
   34 FORMAT(' Prob. of Encountering Cloud:',2X,4I12)
C
C-------Now write the height conditional cloud data to the screen (and
C       printer/disk file).
C
      DO 100 I = 1,NDECKS
C
      N1=WATR(I,1,1)+0.001
      N2=WATR(I,1,2)+0.001
C
      IF(N1.LT.0.OR.N2.LT.0) THEN
        GOTO 100
      ELSEIF(N1.EQ.0.OR.N2.EQ.0) THEN
        MAXCNT=2
      ELSE
        MAXCNT=N2-N1+3
      ENDIF
C
C-------Write the cloud layer heading
C
      WRITE(K6,41) I
   41 FORMAT(/35X,'Layer - ',I1//' Height Pressure  Temp    LWC    Prob
     1Cloud    Fix-Wng Icing    Rot-Wng Icing'/'   (m)    (mb)     (C) (
     2gm/cu m) All Liquid   Severity Index   Severity Index'/)
C
      DO 200 L=1,MAXCNT
C
C-------Get the temperature for the height level and convert to Celsius.
C
      IF(L.EQ.1) THEN
        Z=BOT(I)
        P=PRB(I)
        T=BTEM(I)-FREEZE
      ELSEIF(L.EQ.MAXCNT) THEN
        Z=TOP(I)
        P=PRT(I)
        T=TTEM(I)-FREEZE
      ELSE
        NZ=N1+L-2
        Z=HEITL(NZ)
        P=PRES(NZ)
        T=TEMP(NZ)-FREEZE
      ENDIF
C
      ALW=WATR(I,1,L+2)
C
C-------Convert WATR(2,I) from decimal to percent format
C
      PROB=100.0*WATR(I,2,L+2)
C
C-------Write the data to the screen (and printer/disk file).
C
      WRITE(K6,42) INT(Z),P,T,ALW,PROB,AINDEX(I,L),HINDEX(I,L)
 42   FORMAT(I7,F9.1,F7.1,F8.3,F10.1,8X,A11,6X,A11)

C
  200 CONTINUE
C
  100 CONTINUE
C
   99 RETURN
      END
      SUBROUTINE OUT_LAPS
     1  (nk,CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE,
     1         rlwc_laps,prob_laps)
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL OUTPUT(CLDTYP,WATR,XMVD,AINDEX,HINDEX,ICP,NDECKS,TITLE)
C    Input:
C       CLDTYP   - (character) cloud type in two ASCII characters
C       WATR     - (real) output array containing LWC and percent liquid
C                  water every 100 m in the cloud
C       XMVD     - (real) Mean Volume Diameter for each cloud
C       AINDEX   - (character) Icing Severity Index for fixed wing
C                  aircraft every 100 m
C       HINDEX   - (character) Icing Severity Index for helicopters
C                  every 100 m
C       ICP      - (integer) percent (1-100%) cloud coverage for up to 4
C                  cloud decks
C       NDECKS   - (integer) the number of cloud decks
C       TITLE    - (character) the title to be displayed in the output
C    Output:
C       None
C    Common:
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C       K5       - (integer) the keyboard input unit
C       K6       - (integer) the screen output unit
C       K7       - (integer) the disk file unit
C       K8       - (integer) the printer output unit
C=======================================================================
C    12/21/89  Steve Albers - Modified for PROFS Map LWC output to
C                             LAPS grid
C=======================================================================
C<End>
C***********************************************************************
C
      REAL       FREEZE
      PARAMETER   (FREEZE = 273.15)
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /IODEV/K5,K6,K7,K8
C
      DIMENSION ICP(4),WATR(4,2,260),XMVD(4)
C
      real rlwc_laps(nk),prob_laps(nk)

      CHARACTER*70 TITLE
C     CHARACTER*14 FILENM
      CHARACTER*11 AINDEX(4,260),HINDEX(4,260)
      CHARACTER*2 CLDTYP(4)
C     CHARACTER*1 NUM(10),ANS
      CHARACTER*1 NUM(10)
C
C     LOGICAL NPRT,NFIL
C
      DATA NUM/'0','1','2','3','4','5','6','7','8','9'/

    1 FORMAT(A,$)

!     Initialize laps arrays
!     do k = 1,nk
!         rlwc_laps(k) = 0.
!         prob_laps(k) = 0.
!     enddo

C
C-------Now write the height conditional cloud data to the screen (and
C       printer/disk file).
C
      DO 100 I = 1,NDECKS
C
      N1=WATR(I,1,1)+0.001
      N2=WATR(I,1,2)+0.001
C
      IF(N1.LT.0.OR.N2.LT.0) THEN
        GOTO 100
      ELSEIF(N1.EQ.0.OR.N2.EQ.0) THEN
        MAXCNT=2
      ELSE
        MAXCNT=N2-N1+3
      ENDIF

      z_new = -999.
      DO 200 L=1,MAXCNT
C
C-------Get the temperature for the height level and convert to Celsius.
C
      IF(L.EQ.1) THEN
        Z=BOT(I)
        P=PRB(I)
        T=BTEM(I)-FREEZE
      ELSEIF(L.EQ.MAXCNT) THEN
        Z=TOP(I)
        P=PRT(I)
        T=TTEM(I)-FREEZE
      ELSE
        NZ=N1+L-2
        Z=HEITL(NZ)
        P=PRES(NZ)
        T=TEMP(NZ)-FREEZE
      ENDIF
C
      ALW=WATR(I,1,L+2)
C
C-------Convert WATR(2,I) from decimal to percent format
C
      PROB=100.0*WATR(I,2,L+2)
C
C-------Write the data to the screen (and printer/disk file).
C

      alw_old = alw
      prob_old = prob
      z_old = z_new

      z_new = zcoord_of_pressure(p*100.)

      if(int(z_old) .ne. int(z_new) .and. z_old .ne. -999.)then
!         Interpolate to LAPS grid point
          iz = int(z_new)
          frac = (z_new - float(iz)) / (z_new - z_old)
          rlwc_laps(iz) = alw_old * frac + alw * (1.0 - frac)
          prob_laps(iz) = prob_old * frac + prob * (1.0 - frac)
      endif


!      WRITE(K6,42) INT(Z),P,T,ALW,PROB,AINDEX(I,L),HINDEX(I,L)
! 42   FORMAT(I7,F9.1,F7.1,F8.3,F10.1,8X,A11,6X,A11)

C
  200 CONTINUE
C
  100 CONTINUE
C
   99 RETURN
      END

      SUBROUTINE LWC(NDECKS,ICP,ICTP,WATR,XMVD)
C***********************************************************************
C                          SUBROUTINE LWC
C***********************************************************************
C<Begin>
C<Identification>          Name:  LWC
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  SMFD
C=======================================================================
C<Description>
C    This routine controls the computation of liquid water content (LWC)
C    and the probability the cloud is all liquid.  It calls the
C    microphysics portion of the code by RTNEPH cloud deck.
C=======================================================================
C<Called routines>
C    LW     - (subroutine) calculates LWC and the probability the cloud
C             is all liquid for super-cooled clouds
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL LWC(NDECKS,ICP,ICTP,WATR)
C    Input:
C       NDECKS   - (integer) number of cloud decks present
C       ICP      - (integer) percent cloud amount in layers
C       ICTP     - (integer) code for the cloud type (up to 4 cloud
C                  decks)
C    Output:
C       WATR     - (real) output array containing LWC and percent liquid
C                  water for every 100 m in the cloud layers
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C    09/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definition of WATR to give values every
C              100 m.  Minor modifications.
C=======================================================================
C<End>
C***********************************************************************
C
C  DESCRIPTION:
C       Loop 10 controls the stepping through the observed RTNEPH cloud
C  decks.  Subroutine LW is called to calculate the Liquid Water Content
C  (LWC) and Thermodynamic Phase (TDP).
C
C  LOCAL GLOSSARY:
C       IRAIN     - NO RAIN (1) Flag
C       LAY       - RTNEPH cloud deck index
C       ICTP(4)   - cloud type by 3DNEPH cloud deck
C
C-----------------------------------------------------------------------
C
      DIMENSION WATR(4,2,260),ICTP(4),ICP(4),XMVD(4)
C
C ------------------- START OF EXECUTABLE CODE -------------------------
C
C-------Set the rain flag to NO RAIN permenantly
C
      IRAIN = 1
C
      DO 10 LAY=1,NDECKS
C
C-------Set the first two values of WATR to indicate no clouds
C
      WATR(LAY,1,1)=-2.0
      WATR(LAY,1,2)=-2.0
C
C-------If the cloud type is 'clear', skip it.
C
      IF(ICTP(LAY).NE.0) CALL LW(LAY,IRAIN,ICTP(LAY),WATR,ICP(LAY)
     1                                          ,XMVD(LAY))
C
  10  CONTINUE
C
      RETURN
      END
      SUBROUTINE LW(LAY,IRAIN,LVLCT,WATR,NEPH,RMVD)
C***********************************************************************
C                          SUBROUTINE LW
C***********************************************************************
C<Begin>
C<Identification>          Name:  LW
C                          Type:  FORTRAN-77 Subroutine
C                      Filename:  ICE.FOR
C                        Parent:  LWC
C=======================================================================
C<Description>
C    This subroutine calculates the LWC and probability the cloud water
C    is all liquid for super-cooled clouds, and calls S/R AL for each
C    height level within a cloud deck.
C=======================================================================
C<Called routines>
C    AL     - (subroutine) calculates moist adiabatic LWC of non-
C             cirriform cloud types
C=======================================================================
C<Parameters>
C    Formal declaration:
C       CALL LW(LAY,IRAIN,LVLCT,WATR,NEPH)
C    Input:
C       LAY      - (integer) cloud layer index (1-4)
C       IRAIN    - (integer) RAIN (2) or NO RAIN (1)
C       LVLCT    - (integer) cloud type of the cloud deck
C       NEPH     - (integer) percent cloud amount in cloud deck
C    Output:
C       WATR     - (real) output array containing LWC and percent liquid
C                  water for every 100 m in the cloud layers
C    Common:
C       BOT      - (real) the height of the (up to 4) cloud layer
C                  bottoms
C       TOP      - (real) the height of the (up to 4) cloud layer tops
C       PRB      - (real) the pressures at the cloud layer bottoms
C       PRT      - (real) the pressures at the cloud layer tops
C       BTEM     - (real) the temperatures at the cloud layer bottoms
C       TTEM     - (real) the temperatures at the cloud layer tops
C       NLVL     - (integer) number of 100 m height levels
C       HEITL    - (real) array of 100 m heights
C       PRES     - (real) the pressure every 100 m
C       TEMP     - (real) the temperature every 100 m
C       K6       - (integer) the screen output unit
C=======================================================================
C<History>
C    11/15/87  UDRI   (513) 229-3921    James K. Luers
C              Delivered basic source code
C    02/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Cleaned up the code
C    05/03/89  ASL    (505) 678-1570    Elton P. Avara
C              Minor modifications to the code and added the variables
C              PRES and TOP.
C    09/13/89  ASL    (505) 678-1570    Elton P. Avara
C              Modified the definitions of PRES and WATR to give
C              values every 100 m.  Also added NLVL, TTEM, TEMP, LAY,
C              and HEITL.  Major modifications.
C=======================================================================
C<End>
C***********************************************************************
C
C  DESCRIPTION:
C       Authors - M. D. Dykton and C. W. Rogers.  Loop 200 finds ILOC
C                 and ITMP and puts the calculated LWC in WATR.
C
C  LOCAL GLOSSARY:
C       ILOC       = (integer) location of point within cloud deck
C       ITMP       = (integer) temperature converted to Table Index
C       LAY        = (integer) cloud layer
C       PERMAX(5,10,2)
C                  = (real) table of the percent of maximum LWC for a
C                    given cloud type (10) at intervals above the cloud
C                    base (5) for raining or non-raining clouds (2).
C                    See TN 74-4, Figures 4, 5, & 6 on PP 6-8 and
C                    paragraph 4A, P 4.  ILOC, LVLCT, and IRAIN are the
C                    indices for PERMAX.
C       T          = (real) temperature in degrees Kelvin
C       THEMAX(10,10)
C                  = (real) table of maximum LWC in G/M**3 that can
C                    occur for a given cloud type (10) and for a given 5
C                    degree temperature interval (10).  See Calspan
C                    Final Report, Table 5, P 26.
C
C  REFERENCES:
C       1) USAFETAC Technical Note 74-4, "A Synoptic-Scale Model for
C          Simulating Condensed Atmospheric Moisture", April 1974, by
C          CAPT Robert G. Feddes.
C
C       2) Rogers, C.W., Hanley J.T. and E.J. Mack, 1985: "Updating the
C          Smith-Feddes Model", Final Report Contract No. N00228-84-C-
C          3157 Calspan Report No. 7330-1. Calspan Corp., P.O. Box 400
C          Buffalo, New York 14225.
C
C-----------------------------------------------------------------------
C
      COMMON /PANT/ BOT(4),TOP(4),PRB(4),PRT(4),BTEM(4),TTEM(4)
      COMMON /QANT/ NLVL,HEITL(260),PRES(260),TEMP(260)
      COMMON /IODEV/K5,K6,K7,K8
C
      DIMENSION CALW(260),WATR(4,2,260),THEMAX(10,10),PERMAX(5,10,2)
C
      DATA THEMAX/.2,.25,.30,.35,.4,.45,.50,.55,.6,.6,
     1            .35,.4,.45,.50,.55,.60,.65,3*.7,10*3.,.35,.4,
     2            .45,.5,.6,.6,.75,3*.90,.25,.30,.35,.40,.40,
     3            .45,.6,3*.7,.20,.25,.25,.3,.35,.40,.45,3*.5,
     4            3*.15,3*.2,4*.25,4*.1,3*.15,3*.2,
     5            4*.05,3*.1,3*.15,10*6.5/
C
      DATA PERMAX/5*.4,.38,.62,.74,.58,.37,.40,.60,.80,
     1           .95,.74,.38,.62,.74,.58,.37,.38,.62,.74,.58,.37,
     2            20*.4,.37,.57,.76,.90,.82,.93,.77,.62,.47,.32,
     3           .96,.88,.74,.58,.37,1.9,1.6,1.,.5,.4,.96,.88,
     4           .74,.58,.37,.96,.88,.74,.58,.37,.93,.77,.62,
     5           .47,.32,.93,.77,.62,.47,.32,.93,.77,.62,.47,
     6           .32,.93,.77,.62,.47,.32,2.75,1.95,1.,.63,.47/
C
C ******************* START OF EXECUTABLE CODE *************************
C
C  IF CLOUD PERCENTAGE IS ZERO, SKIP MICROPHYSICS.
C
      IF(NEPH.EQ.0) GOTO 100
C
C  SET S/R AL PARAMETER
C
      TI=BTEM(LAY)-273.16
C
C-------Determine which 100 m height levels are within this cloud layer.
C
      N1=1
      N2=NLVL
      DO 400 J=1,NLVL
      IF(BOT(LAY).GE.HEITL(J)) N1=J+1
      IF(TOP(LAY).GT.HEITL(J)) N2=J
  400 CONTINUE
      NMAX=N2-N1+3
C
      IF(N1.GT.N2) THEN
        N1=0
        N2=0
        NMAX=2
      ENDIF
C
      DIF=TOP(LAY)-BOT(LAY)
C
C-------Store the indices which indicate which 100 m height levels are
C       within the cloud in the first two elements of WATR.
C
      WATR(LAY,1,1)=N1
      WATR(LAY,1,2)=N2
C
C-------If the cloud type is not cirriform, calculate the adiabatic
C       liquid water content within the cloud layer
C
      IF(LVLCT.LT.7.OR.LVLCT.GT.9) CALL AL(LVLCT,LAY,TI,N1,N2,CALW)
C
C-------Get the liquid water content at all height levels within the
C       cloud layer
C
      DO 200 L=1,NMAX
C
C-------Get the height and temperature for each height level within the
C       cloud layer
C
      IF(L.EQ.1) THEN
        Z=BOT(LAY)
        T=BTEM(LAY)
      ELSEIF(L.LT.NMAX) THEN
        Z=HEITL(N1+L-2)
        T=TEMP(N1+L-2)
      ELSE
        Z=TOP(LAY)
        T=TTEM(LAY)
      ENDIF
C
C  SCALE TEMPERATURE INTO AN INDEX (ITMP) FOR THE ARRAY THEMAX
C
      ITMP = (T - 238.0) / 5.0
      IF(ITMP.GT.10) ITMP = 10
      IF(ITMP.LT. 1) ITMP = 1
C
C  CALCULATE PERCENT HEIGHT OF POINT ABOVE BASE WITHIN CLOUD
C  DECK AND SCALE IT INTO AN INDEX FOR THE TABLE PERMAX.
C
      HGTPCT = (Z - BOT(LAY)) / DIF
      ILOC = HGTPCT * 5. + 1.
      IF(ILOC.LT.1) ILOC = 1
      IF(ILOC.GT.5) ILOC = 5
C
C  IF NOT CIRRIFORM CLOUD, MODIFY ADIABATIC CMC.
C
      IF(LVLCT .LT. 7 .OR. LVLCT .GT. 9) THEN
C
C  NO REDUCTION IN CMC NEAR CLOUD TOP IF STRATUS TYPE CLOUDS
C  (ST, AS, AND NS).  FOR CUMULIFORM CLOUDS REDUCE CMC NEAR CLOUD TOP
C  IF HEIGHT GREATER THAN 80% OF CLOUD DEPTH
C
        IF(LVLCT.NE.1.AND.LVLCT.NE.4.AND.LVLCT.NE.6) THEN
          IF(HGTPCT.GE.0.8) THEN
            LL=LVLCT
            IF(LL.EQ.25) LL=3
            DELTOP=5.0*(HGTPCT-0.8)
            CALW(L)=CALW(L)*(1.0-DELTOP*(1.0-PERMAX(ILOC,LL,1)))
          ENDIF
        ENDIF
C
C  IF PRECIPITATION IS PRESENT, MODIFY THE LIQUID WATER CONTENT
C
        IF(IRAIN.EQ.2) THEN
          PERDIF=PERMAX(ILOC,LVLCT,2)-PERMAX(ILOC,LVLCT,1)
          IF(PERDIF .LE. 0.0) THEN
C
C  DECREASE CMC IN PRECIPITATION
C
            CALW(L)=CALW(L)*PERMAX(ILOC,LVLCT,2)/PERMAX(ILOC,LVLCT,1)
          ELSE
C
C  INCREASE CMC IN PRECIPITATION
C
            PDIF=PERDIF/PERMAX(ILOC,LVLCT,2)
            CALW(L)=CALW(L)/(1.0-PDIF)
          ENDIF
        ENDIF
C
        WATR(LAY,1,L+2)=CALW(L)
C
C  CALCULATE FRACTIONAL PROBABILITY CLOUD WATER IS ALL LIQUID
C
        IF( T .GT. 273. ) THEN
C
C  ALL CLOUD WATER IS LIQUID IF TEMPERATURE > 273.
C
          WATR(LAY,2,L+2) = 1.
        ELSEIF( T .LT. 233. ) THEN
C
C  ALL CLOUD WATER IS ICE IF TEMPERATURE < 233.
C
          WATR(LAY,2,L+2) = 0.
        ELSE ! T is between 233 and 273
C
C  CHECK CLOUD TYPE IF (233. <= TEMPERATURE <= 273.)
C  THERE IS SUPER-COOLED WATER
C
          IF(LVLCT.EQ.3.OR.LVLCT.EQ.10) THEN ! Cu type cloud, 233 < T < 273
C
C  SUPER-COOLED UNSTABLE - CUMULUS OR CUMULONIMBUS CLOUDS
C  -.03927 = (PI/2.) / -40.
C
            WATR(LAY,2,L+2)=COS(-.03927*(T-273.))
          ELSE ! not CU type cloud (Stratiform cloud), 233 < T < 273
C
C  SUPER-COOLED STABLE - STRATIFORM CLOUDS (0C to -30C)
C  CLOUD WATER IS DEPLETED BY A TEMPERATURE DEPENDENT HYPERBOLIC
C  TANGENT FUNCTION-NO CLOUD WATER AT -40C AND ALL WATER AT 0C
C  THE "WIDTH" OF THE DROP-OFF BETWEEN ALL WATER AND NO WATER
C  IN THE CLOUD IS CONTROLLED BY THE DENOMINATOR OF ARGUMENT.
C  THE MIDDLE OF THE DROPOFF IS CURRENTLY AT -20C(253K) SINCE
C  THIS VALUE ENSURES MATHEMATICAL CONSISTENTCY
C
C
C           ARG=(T-253.)/7.5
C           FAC=(TANH(ARG)+1.0)/2.
C
C    11/6/90 CHANGE-USE A LINEAR FACTOR FOR LWC DEPLETION BY GLACIATION
C    INSTEAD OF THE HYPERBOLIC PROFILE ABOVE (Ramps from -10C to -30C)
C
            IF(T.GE.263.) THEN
              FAC=1.0
            ELSEIF(T.LE.243.) THEN
              FAC=0.0
            ELSE
              FAC=(T-243.)/20.
            ENDIF

            WATR(LAY,1,L+2)=WATR(LAY,1,L+2)*FAC
            WATR(LAY,2,L+2) = EXP(.0909091*(T-273.)) - .0263
          ENDIF ! Cu type Cloud
        ENDIF ! T

        IF(T .LE. 233.)WATR(LAY,1,L+2)=0. ! Case for Cu or Stratiform cloud < 40C

      ELSE ! Cirroform cloud
C
C  IF CIRRIFROM CLOUD DO NOT USE ADIABATIC LWC, USE TABLE.
C
        WATR(LAY,1,L+2) = THEMAX(ITMP,LVLCT) * PERMAX(ILOC,LVLCT,IRAIN)
C
C  CLOUD IS ALL ICE
C
        WATR(LAY,2,L+2) = 0.0
      ENDIF
C
C-------Check to make sure that the LWC calculated by this program
C       is non-negative.  If it is negative then tell the user to
C       check their input values for errors.
C
      IF(WATR(LAY,1,L+2).LT.0.0) THEN
        WRITE(K6 ,*)'The model has calculated a negative LWC.'
        WRITE(K6 ,*)'Please check your input data for errors.'
        STOP
      ENDIF
C
  200 CONTINUE
      NDIF=N2-N1
      NMID=(NDIF)/2+2
      NMID1=NMID+1
      IF((NDIF/2)*2.EQ.NDIF) THEN
       TMID=TEMP(NMID)
      ELSE
       TMID=(TEMP(NMID)+TEMP(NMID1))/2.
      ENDIF
!     CALL DSD(LAY,LVLCT,IRAIN,ILOC,PERMAX,WATR,TMID,N1,N2,RMVD)
C
  100 RETURN
      END

      SUBROUTINE DSD
     G              ( LAY, LVLCT, IRAIN, ILOC, PERMAX,
     B                WATR,TEMPK,N1,N2,MVD )
C
C    PURPOSE
C         THIS SUBROUTINE CALCULATES DROP SIZE DISTRIBUTION
C    PARAMETER LIST
C         GIVEN
C              LVLCT (I) = CLOUD TYPE (1-10)
C              IRAIN  (I) = NO RAIN (1) OR RAIN (2)
C              ILOC   (I) = LOCATION OF A CALCULATION POINT
C                           WITHIN A CLOUD DECK
C              LVL    (I) = CLOUD LEVEL
C              PERMAX(5,10,2) (R) = PERCENT OF MAXIMUM LWC THAT CAN
C                                   A FUNCTION OF CLOUD TYPE AND
C                                   RELATIVE CLOUD POSITION
C         YIELDED
C              RNDRPS(4,40) (R)    = DROP SIZE DITRIBUTION AT
C                                    THE MID-POINT OF EACH CLOUD
C                                    DECK
C                        MVD = MEAN VOLUME DIAMETER PER LAYER
C
C    LOCAL GLOSSARY
C              DRPCLD (R) = RADII AT WHICH CLOUD DSD IS EVALUATED
C              CLDRPS (R) = NUMBER OF DROPS/CM**3/MICRON RADIUS INTERVAL
C                           CENTERED AT A PARTICULAR RADIUS
C              ILAYR  (R) = LAYER NUMBER (SAME AS LAY)
C              LWCC   (R) = CLOUD CONDENSED MOISTURE
C              LWCCL  (R) = CLOUD LIQUID CONTENT
C              LWCCI  (R) = CLOUD ICE CONTENT
C              LWCR   (R) = RAIN CONDENSED MOISTURE
C              LWCRL  (R) = RAIN LIQUID CONTENT
C              LWCRI  (R) = RAIN ICE CONTENT
C LWC=TOTAL CONDENSED MOISTURE, PRECIPITATION OR NO PRECIPITATION, LIQUI
C  CLOUD.  USED IN BOTH NO PRECIPITATION AND PRECIPITATION CASES.
C  TDP=TOTAL LIQUID CONDENSED MOISTURE, PRECIPITATION OR NO PRECIPITATIO
C  TDP=LWC FOR LIQUID CLOUD
C  TDP=0.0 FOR ICE CLOUD
C  USED ONLY IF ALL MOISTURE IS ONLY IN CLOUD
C              TEMPK  (R) = LAYER MIDPOINT TEMPERATURE (DEGREES KELVIN)
C
C.......................................................................
C
C
      DIMENSION  VOLI(40),CLDRPS(4,40),RNDRPS(10,4,2),
     1           PERMAX(5,10,2), DRPCLD(40),WATR(4,2,260)
      REAL       LWC, LWCC, LWCR, LWCCL, LWCCI, LWCRL, LWCRI
      REAL MVD,MVDI(40)

      DATA DRPCLD/1.0,3.0,5.0,7.0,9.0,11.0,13.0,15.0,17.0,19.0,
     1 21.0,23.0,25.0,27.0,29.0,31.0,33.0,35.0,37.0,39.0,41.0,
     2 43.0,45.0,47.0,49.0,51.0,53.0,55.0,57.0,59.0,61.0,63.0,
     3 65.0,67.0,69.0,71.0,73.0,75.0,77.0,79.0/
C
C ******************** START OF EXECUTABLE CODE ************************
C
C     PRESETTING LWC VARIABLES
C
      LWC   = -99.0
      LWCC  = -99.0
      LWCR  = -99.0
      LWCCL = -99.0
      LWCCI = -99.0
      LWCRL = -99.0
      LWCRI = -99.0
C
C   FOR THE TIME BEING WE WILL CALCULATE A DROP SIZE DISTRIBUTION ONLY
      ILAYR=LAY

C   AT THE MID-POINT OF A CLOUD DECK
      NDIF=N2-N1
      NMID=NDIF/2+1
      NMID1=NMID+1
      IF((NDIF/2)*2.EQ.NDIF) THEN
       LWC=WATR(LAY,1,NMID)
      ELSE
       LWC=(WATR(LAY,1,NMID)+WATR(LAY,1,NMID1))/2.
      ENDIF
      TDP=LWC
C  IF NO PRECIPITATION AT THIS LEVEL ALL MOISITURE IS CLOUD WATER
      IF ( .NOT. ( IRAIN .EQ. 2 ) ) GO TO 20
C
C         IT'S PRECIPITATING; CALCULATE THE LWC PERCENT DIFFERENCE
C         BETWEEN PRECIPITATING & NON-PRECIPITATING CLOUD TYPES
 10   DIF  = PERMAX(ILOC,LVLCT,2) - PERMAX(ILOC,LVLCT,1)
      PDIF = DIF / PERMAX(ILOC,LVLCT,2)
      IF (.NOT.( DIF .LE. 0. )) GO TO 60
C
C     ALL LIQUID IS IN CLOUD FORM - NO PRECIPITATION
C
   20 LWCCL = TDP
      LWCCI = LWC - LWCCL
C          CLOUD LIQUID CONTENT
C          CLOUD ICE CONTENT
C
C     FOLLOWING COMPUTES AND STORES DROP SIZE DISTRIBUTIONS
C     CORRESPONDING TO LWCCL AND LWCCI
C
C
      IF(LWCCL .GT. 0.)then
          CALL NDROPS
     G           ( LVLCT, LWCCL, DRPCLD, 40,LAY,
     Y             CLDRPS    )
C
C-----CALCULATE THE VOLUMES OF EACH OF THE RADII CATEGORIES
C
      TOTVOL=0.0
      PI=ACOS(-1.0)
      DO 25 K=1,40
         VOLI(K)=(4.0/3.0)*PI*(DRPCLD(K)*1E-4)**3
         IWATR=INT(CLDRPS(LAY,K))
         TOTVOL=TOTVOL+IWATR*VOLI(K)
  25  CONTINUE
C
C-----CALCULATE THE MEAN VOLUME DIAMETER FOR THIS LAYER
C
      DO 35 K=1,40
         L=41-K
!        WRITE(3,*) CLDRPS(LAY,K)
         IWATR=INT(CLDRPS(LAY,L))
C         WRITE(3,*)IWATR,VOLI(L),TOTVOL
         MVDI(L)=IWATR*VOLI(L)/TOTVOL
C         WRITE(3,*)MVDI(L)
   35 CONTINUE
      TOTMVD=0.0

      MVD=0.0
      DO 40 K=1,40
         L=41-K
         TOTMVD=TOTMVD+MVDI(L)
         IF (TOTMVD.GE.0.5) THEN
             MVD=DRPCLD(L)
C            WRITE(3,*)MVD
             GO TO 45
         ENDIF
   40 CONTINUE
      ENDIF ! IF NO LWC
C
   45 IF(LWCCI.EQ.0.)GO TO 50
C
C     CALL NDROPS
C    G           (8, LWCCI, DRPCLD, 8,
C    Y             WATR(13,ILAYR)     )
C
C
   50 CONTINUE
      RETURN
C
C     FOLLOWING TREATS CASES WHERE SOME LIQUID IS IN RAIN FORM
C
 60   CONTINUE
C SPLIT LWC INTO RAIN AND CLOUD WATER
      LWCR = PDIF * LWC
      LWCC = LWC - LWCR
C
C     TEMPK IS TEMPERATURE AT THE LAYER MIDPOINT
      IF((LVLCT .GE. 7) .AND. (LVLCT .LE. 9)) GO TO 400
C          IF TEMPERATURE > 273. THEN ALL MOISTURE IS LIQUID
      IF ( TEMPK .GT. 273. ) LWCCL = 1.
C          IF TEMPERATURE < 233. THEN ALL MOISTURE IS ICE
      IF ( TEMPK .LT. 233. ) LWCCL = 0.
C          SKIP IF NO SUPER-COOLED LIQUID IS PRESENT
      IF ( TEMPK .GT. 273. .OR. TEMPK .LT. 233. ) GO TO 210
      IF ( LVLCT .EQ. 3 .OR. LVLCT .EQ. 10 ) GO TO 200
C
C     SUPER-COOLED STABLE -- STRATIFORM CLOUDS
      LWCCL = EXP(.0909091*(TEMPK-273.))-.0263
      GO TO 210
C
 200  CONTINUE
C     SUPER-COOLED UNSTABLE -- CUMULUS OR CUMULONIMBUS CLOUDS
C     -.03927 = ( PI / 2. ) / -40.
      LWCCL = COS (-.03927 * (TEMPK - 273.))
C
 210  CONTINUE
C  FRACTIONAL PROBABILITY CLOUD IS ALL WATER
      WATER=LWCCL
C SET LWCCL FOR LIQUID CLOUD PROCESSING
C  LWCCL=1.0. INSURES LWCCL CONVERTS TO LWCC IN 2ND STATEMENT HENCE.
      LWCCL=1.0
      IF(WATER .EQ. 0.0) LWCCL=0.0
C  LIQUID CLOUD VARIABLES SET FOR S/R NDROPS AND S/R WRTOUT.
      LWCCL = LWCCL * LWCC
      LWCCI = LWCC - LWCCL
      WATL = LWCCL
      WATI = LWCCI
      GO TO 500
C  ICE CLOUD VARIABLE SET FOR S/R NDROPS AND S/R WRTOUT
C  IF CIRRIFORM CLOUD THEN NO CLOUD LIQUID CONTENT AND PROBABILITY OF CL
C  CONDENSED MOISTURE BEING ALL LIQUID IS ZERO
  400 WATL=0.0
      WATI=LWCC
      LWCCI=LWCC
      LWCCL=0.0
C
C     FOLLOWING COMPUTES AND STORES DROP SIZE DISTRIBUTIONS
C     CORRESPONDING TO LWCCL AND LWCCI
  500 CONTINUE
C
C
      IF(LWCCL .GT. 0.)then
          CALL NDROPS
     G           ( LVLCT, LWCCL, DRPCLD, 40,LAY,
     :             CLDRPS  )
C
C-----CALCULATE THE VOLUMES OF EACH OF THE RADII CATEGORIES
C
      TOTVOL=0.0
      PI=ACOS(-1.0)
      DO 525 K=1,40
         VOLI(K)=(4.0/3.0)*PI*(DRPCLD(K)*1E-4)**3
         IWATR=INT(CLDRPS(LAY,K))
         TOTVOL=TOTVOL+IWATR*VOLI(K)
C         TWATR=WATR(J,ILAYR)
C         TOTVOL=TOTVOL+TWATR*VOLI(K)
 525  CONTINUE
C
C-----CALCULATE THE MEAN VOLUME DIAMETER FOR THIS LAYER
C
      DO 535 K=1,40
         L=41-K
C         TWATR=WATR(J,ILAYR)
         IWATR=INT(CLDRPS(LAY,L))
C        WRITE(3,*)IWATR,VOLI(L),TOTVOL
C         MVDI(L)=TWATR*VOLI(L)/TOTVOL
         MVDI(L)=IWATR*VOLI(L)/TOTVOL
C         WRITE(3,*)MVDI(L)
  535 CONTINUE
      TOTMVD=0.0
      MVD=0.0
      DO 540 K=1,40
         L=41-K
         TOTMVD=TOTMVD+MVDI(L)
         IF (TOTMVD.GE.0.5) THEN
             MVD=DRPCLD(L)
C            WRITE(3,*)MVD
             GO TO 45
         ENDIF
  540 CONTINUE
      ENDIF ! IF NO LWC
C
C
      IF(LWCCI.EQ.0.)GO TO 90
C
C      CALL NDROPS
C    G           ( 8, LWCCI, DRPCLD,40,
C    Y             WATR(13,ILAYR)  )
C
C
   90 CONTINUE
C     THE FACTOR QZ IS INCORPORATED INTO THE MODEL TO ACCOUNT FOR
C     THE FACT THAT "THE DIFFERENTIAL FALL VELOCITIES FOR LIQUID
C     AND SOLID WATER WILL DECREASE THE TOTAL CONDENSED MOISTURE
C     BELOW THE FREEZING LEVEL."  (TN 74-4, P 11, PARA. A)
C     IN THOSE LAYERS WHOSE HEIGHTS ARE BELOW THE FREEZING LEVEL,
C     THE CONDENSED MOISTURE IS HALVED, FOR CONVECTIVE CLOUD TYPES.
      QZ = 1.
C
C     FOLLOW SAME COMPUTATIONAL PROCEDURE FOR RAIN AS FOR CLOUDS (ABOVE)
      IF ( TEMPK .GT. 273. ) LWCRL = 1.
      IF ( TEMPK .LT. 233. ) LWCRL = 0.
      IF ( TEMPK .GT. 273. .AND. LVLCT .EQ. 10 ) GO TO 230
      IF ( TEMPK .GT. 273. .AND. LVLCT .EQ. 3 ) GO TO 230
      IF ( TEMPK .GT. 273. .OR.  TEMPK .LT. 233. ) GO TO 240
      IF(LVLCT .EQ. 3 .OR. LVLCT .EQ. 10 ) GO TO 220
C
C     SUPER-COOLED STABLE -- NON-CONVECTIVE SO FACTOR DOESN'T APPLY
      LWCRL = EXP(.0909091 * (TEMPK-273.))-.0263
C  PRECIPITATION AT TEMPERATURES BELOW ZERO IS FROZEN. STRATIFORM CLOUDS
      LWCRL=0.0
      GO TO 240
C
 220  CONTINUE
C     SUPER-COOLED UNSTABLE
      LWCRL =  COS(-.03927  *  (TEMPK - 273.))
C  PRECIPTATION AT TEMPERATURES BELOW ZERO IS FROZEN.  CUMULIFORM CLOUDS
      LWCRL=0.0
C
 230  CONTINUE
C     BRANCH HERE IF CONVECTIVE CLOUDS; TEST TO SEE IF TEMPERATURE
C     ABOVE FREEZING (SO THAT LAYER IS BELOW FREEZING POINT).
      IF ( TEMPK .GT. 273.) QZ = .5
C
 240  CONTINUE
C
      LWCRL = (LWCRL * QZ) * LWCR
      LWCRI = (LWCR * QZ) - LWCRL
C
C     LWC = LWCCL + LWCCI + LWCRL + LWCRI
C
C
C     COMPUTE DROP SIZE DISTRIBUTION (NUMBER/M**3/MICRON) AT 300
C     MICRON INTERVALS FROM 150 TO 2850 MICRONS, CORRESPONDING TO ...
C
C     LWCRL (RAIN LIQUID):
      IF ( .NOT. ( LWCRL .NE. 0 ) ) GO TO 105
      R = -150.
      DO 100 IA = 1,10
           R = R+300.
  100      RNDRPS(IA,LAY,1)= 20.*EXP((-.004555*R)/LWCRL**.25)
C          SEE EQUATION, P 8, TN 74-4 (WHICH IS FROM P 17, TN 74-1)
C     ENDDO
C
  105 CONTINUE
C
C     LWCRI (RAIN ICE):
      IF (.NOT.( LWCRI .NE. 0 )) GO TO 115
      R = -150.
      DO 110 IA = 1,10
           R = R+300.
  110      RNDRPS(IA,LAY,2)= 20. * EXP((-.004555*R)/LWCRI**.25)
C          SEE EQUATION, P 8, TN 74-4
C     ENDDO
C
C
  115 RETURN
      END
C
C

      SUBROUTINE NDROPS(LVLCT,LWC,DROPSZ,N,LAY,
     1 CLDRPS)
CC
CC   THE CLOUD DROP SIZE DISTRIBUTION CALCULATED HEREIN IS BASED
CC   ON BERRY AND REINHARDT(JOURNAL OF ATMOSPHERIC SCIENCE, 1974)
CC   THE TOTAL DROP CONCENTRATION MUST BE SPECIFIED.  FOR NOW,
CC   THIS IS SPECIFIED FOR ALL CLOUD TYPES.
CC   LATER I INTEND TO MAKE THIS A FUNCTION OF CLOUD BASE VERTICAL
CC   VELOCITY AND TEMPERATURE.
CC
CC   PARAMETER LIST
CC
CC    GIVEN
CC     LVLCT(I) = CLOUD TYPE
CC     LWC(R)    = LIQUID WATER CONTENT(GM/CM**3)
CC     DROPSZ(R) = NUMBER OF DROPS PER 2 MICRON DIAMETER INTERVAL
CC     N(I)      = DIMENSION OF DROPSZ
CC     LAY(I)    = CLOUD DECK INDEX
CC    YIELDED
CC     CLDRPS(R) = NUMBER OF DROPS/CM**3/2 MICRON INTERVAL
CC
CC   -------------------------------------------------------------
      REAL PI,DENSITY
      PARAMETER(PI=3.14159,DENSITY=1.0)
      REAL LWC,CLDRPS(4,40),DROPSZ(40),TND(10)
      INTEGER GMAFCT
      DATA GAMMA/2./
      DATA TND/2.5E2,3.0E2,5.0E2,2.5E2,2.5E2,3.0E2,2.5E2,
     1 3.0E2,2.5E2,5.0E2/
CC
CC
      LWC=LWC*1E-6
      GMAFCT=1
      DO 100 I=1,GAMMA
       GMAFCT=GMAFCT*I
 100  CONTINUE
CC
      G=(1+GAMMA)**(1+GAMMA)/FLOAT(GMAFCT)
CC
      DDIAM=2.0E-4
      DO 200 I=1,40
       DIA=DROPSZ(I)/1E4
       DVOL=PI*DIA*DIA*DDIAM/2.
       VOL=PI*DIA**3/6.
       S=TND(LVLCT)*VOL*DENSITY/LWC
       ARG=(1+GAMMA)*S
       IF(ARG.GT.75) ARG=75
       F=TND(LVLCT)**2*G*S**GAMMA*EXP(-ARG)/LWC
       CLDRPS(LAY,I)=F*DVOL*DENSITY
 200  CONTINUE
      RETURN
      END




