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
From benjamin@apache.fsl.noaa.gov Fri Jun  9 19:04:57 1995
Return-Path: <benjamin@apache.fsl.noaa.gov>
Received: from fslg8.fsl.noaa.gov by peaks.fsl.noaa.gov (4.1/SMI-4.1)
        id AA00472; Fri, 9 Jun 95 19:04:55 GMT
Received: by fslg8.fsl.noaa.gov (5.57/Ultrix3.0-C)
        id AA28278; Fri, 9 Jun 95 13:03:42 -0600
Received: by apache.fsl.noaa.gov (4.1/SMI-4.1)
        id AA03698; Fri, 9 Jun 95 13:04:47 DST
From: benjamin@apache.fsl.noaa.gov (Stan Benjamin)
Message-Id: <9506091904.AA03698@apache.fsl.noaa.gov>
Subject: color NCAR GKS
To: albers@apache.fsl.noaa.gov
Date: Fri, 9 Jun 1995 13:04:47 -0600 (DST)
X-Mailer: ELM
Content-Type: text
Content-Length: 21881
Status: RO

Steve,

Sorry it took me so long to get back to you.
Here are 2 routines.  The first is for changing
the colors of contours, for instance
for different fields overlaid.  The second is for solid
color for a single field.

Note that you need to call GSCR to set up the
color identifiers.

A lot of this stuff, of course, is my code, and you
will want to extract the essential nuggets.  There
are a few routines that need to be called before
and after in your driver.  I assume you already kno
w about those for GKS, but let me know if you don't

Stan

Routine 1


      SUBROUTINE PLOT_CONTOUR(NF, MX, NX, NY, F, LAB1, FINT, FMIN, FMAX,
     + U,V,IVAR,TIMELABEL,NC,ISM,IHL)
C
C PURPOSE       Plot contours of F(NX,NY).
C
C ARGUMENTS     NF --- I  flash buffer to be used
C               MX --- I  maximum size of first dimension of F
C               NX --- I  actual size of first dimension of F
C               NY --- I  actual size of second dimension of F
C               F ---- I  two dimensional array to be contoured
C               LAB1 - I  user-supplied plot label
C               FINT - I  contour interval (if 0., .1*(FMAX-FMIN) is used)
C               FMIN -  O minimum grid value
C               FMAX -  O maximum grid value
C               ISM  - I  flag as whether to smooth on 2x grid
C               IHL  - I  flag as whether to plot high/low labels
C
C      INCLUDE 'CONFIG1'
C      INCLUDE 'CONFIG3'
C


      PARAMETER (N2X=400)
      INTEGER NF, MX, NX, NY, IVAR
      REAL F1(N2X,N2X),HOLD(N2X,2)
      REAL F(MX,NY), FINT, FMIN, FMAX, machine_epsilon_p
     +   ,U(MX,NY),V(MX,NY)

      CHARACTER*48 LAB3
      CHARACTER*48 TIMELABEL
      CHARACTER*48 LAB1

        PARAMETER (LRWK=5500, LIWK=4000)

        REAL RWRK(LRWK)
        INTEGER IWRK(LIWK)

        REAL SCALE
        DATA SCALE /1.1/

      parameter (machine_epsilon_p = 1.19e-07)  ! from iftran.im - bj
      COMMON /SAVMAP/
     .  MX1, MX2, MY1, MY2,
     .  U1,  U2,  V1,  V2
      COMMON /FXFY1/ XA, YA, UA, VA, U1A, V1A, DUDI, DVDJ

        COMMON /NXNY/ NX1,NY1

c ... compare,spv def'n from iftran.im/b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      spv(a) = compare(a,SPVAL_P,0.01)

        SPVAL_P=99999.

        FMIN0 = FMIN
        FMAX0 = FMAX

        NX1 = NX
        NY1 = NY
        NX2X = NX*2-1
        NY2X = NY*2-1

        MAXCHR = 'X'
        MINCHR = 'N'

C
C     CALC MAX AND MIN OF F
C
      FMAX = -1.E8
      FMIN =  1.E8
      DO 55 J=1,NY
        DO 60 I=1,NX
          IF (SPV(F(I,J)) .EQ. NO_P)
     .    THEN
            IF (F(I,J).GT.FMAX) FMAX=F(I,J)
            IF (F(I,J).LT.FMIN) FMIN=F(I,J)
          END IF
60      CONTINUE
55    CONTINUE
C
C     TEST FOR ZERO FINT
C
      IF (FINT .LE. MACHINE_EPSILON_P) THEN
        DF = .1*(FMAX-FMIN)
      ELSE
        DF = FINT
      END IF

      WRITE(6,5) FMIN,FMAX,DF
5     FORMAT(/,' MIN = ',G12.5,5X,'MAX = ',
     .  G12.5,5X,' PLOT INTERVAL = ',G12.5)
C
C    PLOT IF THERE IS ADEQUATE FIELD VARIATION
C

        write(LAB3,101)  FMAX,FMIN,DF
  101   FORMAT(2X,'MAX=',F10.3,2X,'MIN=',F10.3,2X,'INT=',F8.3,2X)

        call gflas3(NF)

        CALL SET (0.,1.,0.,1.,0.,1.,0.,1.,1)

        IF (NC.EQ.1) CALL PLCHHQ (0.0,0.88,TIMELABEL,0.018,0.,-1.)
        CALL PCSETI ('CC - CHARACTER COLOR', 7+NC)
        CALL PLCHHQ (0.45,0.14-(NC-1)*0.02, LAB3   , .01,
     1       0., -1.)
        CALL PLCHHQ (0.05,0.14-(NC-1)*0.02, LAB1   , .01,
     1       0., -1.)

        CALL SET(MX1,MX2,MY1,MY2,U1,U2,
     .    V1,V2,1)                          ! SUPMAP SET
        FMIN = CEILING(FMIN,DF,MACHINE_EPSILON_P)   ! Minimum contour
        IOFFP=1
        SPVALU=SPVAL_P
        IOFFM=1            ! CONREC parameters
        CALL FXFY(NF)                               ! Initialize FX and FY fcns


        NULBLL = 1         ! number of unlabelled lines between labelled ones
        NHI = -1
        IF (IVAR.EQ.8) NHI=0

       IF (ISM.EQ.1) THEN
C --- CALL SET FOR 2X HIGH-RESOLUTION GRID

        CALL SET(MX1,MX2,MY1,MY2,
     .    U1,U1+(NX2X-1)*DUDI,
     .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
C --- Interpolate to 2x resolution grid.
        CALL SUBGRD  (F,MX,NX,NY,F1,N2X,NX2X,NY2X)
C --- Light smoother on high-resolution grid
        CALL SMOOTH2 (F1,HOLD,N2X,NX2X,NY2X, 0.5)
        CALL SMOOTH2 (F1,HOLD,N2X,NX2X,NY2X,-0.5)
       END IF
C --- Do contouring

        call cpseti ('CLS - contour level selection flag',20)
        call cpsetr ('CIS',DF)
        call cpseti ('LIS - label interval specifier',2)
c       call cpsetr ('CMN',FMIN0)
c       call cpsetr ('CMX',FMAX0)

       if (ihl.ge.1) then
        call CPSETR('HLS  - HIGH/LOW LABEL SIZE',.025)
        call CPSETC('ILT',' ')
        call cpseti('NSD',2)
        call cpseti('NLS',2)
        CALL CPSETI ('HLC - HIGH/LOW LABEL COLOR INDEX', 7+NC)
        if (ihl.eq.2) call cpsetc('LOT',' ')
        if (ihl.eq.3) call cpsetc('HIT',' ')
       end if

       IF (ISM.EQ.1) THEN
        CALL CPRECT (F1,N2X,NX2X,NY2X,RWRK,LRWK,IWRK,LIWK)
        CALL CPPKCL (F1,RWRK,IWRK)
       if (ihl.ge.1) then
        call cplbdr (f1,rwrk,lwrk)
       end if
       ELSE
        CALL CPRECT (F,NX,NX,NY,RWRK,LRWK,IWRK,LIWK)
        CALL CPPKCL (F,RWRK,IWRK)
       if (ihl.ge.1) then
        call cplbdr (f,rwrk,lwrk)
       end if
       END IF


        CALL CPGETI ('NCL - NUMBER OF CONTOUR LEVELS', NCON)
        DO 111 I=1,NCON
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX', I)
          call cpgetr ('CLV - contour level values',cval)
          if (cval.lt.0.) then
            call cpseti ('CLD',61166)
          else
            call cpseti ('CLD',65535)
          end if
          CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', 7+NC)
111     CONTINUE


       IF (ISM.EQ.1) THEN
c       CALL SET(MX1,MX2,MY1,MY2,
c    .    U1,U1+(NX2X-1)*DUDI,
c    .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
c       CALL CONREC (F1,N2X,NX2X,NY2X,FMIN,FMAX,DF,1,NHI,-'1470'O) ! Plot contour
        CALL CPCLDR (F1,RWRK,IWRK)
       ELSE
c       CALL SET(MX1,MX2,MY1,MY2,
c    .    U1,U1+(NX2X-1)*DUDI,
c    .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
c       CALL CONREC (F1,N2X,NX2X,NY2X,FMIN,FMAX,DF,1,NHI,-'1470'O) ! Plot contour
        CALL CPCLDR (F,RWRK,IWRK)
       END IF

        IF (IVAR.EQ.5) THEN
        print *,' enter scale (suggested value = 1.5)'
        accept *,scale
        print *,' enter increment for plotting analyzed winds (suggested
     1 = 2)'
        accept *,iskip
        print *,' ENTER MIN SPEED FOR PLOTTING WINDS'
        ACCEPT *,SPDMIN
      CALL GETSET (MX1,MX2,MY1,MY2,
     1  DUMMY,DUMMY,DUMMY,DUMMY,IDUMMY)
      CALL SET (MX1,MX2,MY1,MY2,1.
     1  ,FLOAT(NX),1.,FLOAT(NY),1)

          DO J = 2,NY,iskip
          DO I = 2,NX,iskip
            CALL BARBS_R(FLOAT(I),FLOAT(J),U(I,J),V(I,J),SCALE,
     1  SPDMIN)
          END DO
          END DO
        END IF
        isize=7
        CALL PWRIT (MX1,MY2,' ',1,0,0,0)
c       CALL FRAME                                  ! Flush frame buffer
        WRITE(6,155)
155     FORMAT(' THE FIELD HAS BEEN PLOTTED')
      RETURN
      END

      SUBROUTINE FXFY(IBUF)
C
C PURPOSE       Initialize parameters for FX and FY functions
C
C      INCLUDE 'CONFIG1'
c     INCLUDE 'CONFIG3'
C
      COMMON /NXNY/ NX,NY
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
      COMMON /SAVMAP/
     .  MX1, MX2, MY1, MY2,
     .  UU1, UU2, VV1, VV2
        SPVAL_P=99999.
      X = SPVAL_P
      Y = SPVAL_P
      U = SPVAL_P
      V = SPVAL_P
      U1 = UU1
      DUDI = (UU2-UU1)/(NX-1.)
      V1 = VV1
      DVDJ = (VV2-VV1)/(NY-1.)
      RETURN
      END
      FUNCTION FX(XIN,YIN)
C
C PURPOSE       Transform grid indices (XIN,YIN) to an x-coordinate in the
C               SUPMAP (U,V) coordinate space.  This routine is used by
C               CONREC and VELVCT.
C
C REMARKS       Subroutine FXFY must be called prior to the first call to this
C               function.
C
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
c ... compare, machine_epsilon from iftran.im / b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      IF (COMPARE(XIN,X,1.19e-07).EQ.NO_P) THEN
        X = XIN
        U = U1 + (XIN-1.)*DUDI
      END IF
      FX = U
      RETURN
      END
      FUNCTION FY(XIN,YIN)
C
C PURPOSE       Transform grid indices (XIN,YIN) to a y-coordinate in the
C               SUPMAP (U,V) coordinate space.  This routine is used by
C               CONREC and VELVCT.
C
C REMARKS       Subroutine FXFY must be called prior to the first call to this
C               function.
C
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
c ... compare, machine_epsilon from iftran.im / b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      IF (COMPARE(YIN,Y,1.19e-07).EQ.NO_P) THEN
        Y = YIN
        V = V1 + (YIN-1.)*DVDJ
      END IF
      FY = V
      RETURN
      END
cccccccccccccccccccccccccccccc  ceiling  cccccccccccccccccccccccccccc
c
c  ceiling - ceiling of x, over interval y, to accuracy z
c
        real function ceiling(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0) then                ! ceiling(x,y,z) = trunc(x+y-z,y)
          ceiling=trunc(x+y-z,y)
        else if (x.lt.0) then           ! ceiling(x,y,z) = -trunc(-x,y)
          ceiling= -trunc(-x,y)
        else                            ! ceiling(x,y,z) = 0.
          ceiling = 0.
        endif
c
        return
        end
c
ccccccccccccccccccccccccccccc  floor  cccccccccccccccccccccccccccc
c
c  floor - floor of x, over interval y, to accuracy z
c
        real function floor(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0.) then               ! floor(x,y,z) = trunc(x,y)
          floor = trunc(x,y)
        else if (x.lt.0.) then          ! floor(x,y,z) = -trunc(-x+y-z,y)
          floor = -trunc(-x+y-z,y)
        else                            ! floor(x,y,z) = 0.
          floor = 0.
        endif
c
        return
        end



Routine 2


      sUBROUTINE PLOT_CONTOUR(NF, MX, NX, NY, F, LAB1, FINT, FMIN, FMAX,
     + U,V,IVAR,TIMELABEL,NC,ISM,IHL)
C
C PURPOSE       Plot contours of F(NX,NY).
C
C ARGUMENTS     NF --- I  flash buffer to be used
C               MX --- I  maximum size of first dimension of F
C               NX --- I  actual size of first dimension of F
C               NY --- I  actual size of second dimension of F
C               F ---- I  two dimensional array to be contoured
C               LAB1 - I  user-supplied plot label
C               FINT - I  contour interval (if 0., .1*(FMAX-FMIN) is used)
C               FMIN -  O minimum grid value
C               FMAX -  O maximum grid value
C               ISM  - I  flag as whether to smooth on 2x grid
C               IHL  - I  flag as whether to plot high/low labels
C
C      INCLUDE 'CONFIG1'
C      INCLUDE 'CONFIG3'
C


      PARAMETER (N2X=200)
      INTEGER NF, MX, NX, NY, IVAR
      REAL F1(N2X,N2X),HOLD(N2X,2)
      REAL F(MX,NY), FINT, FMIN, FMAX, machine_epsilon_p
     +   ,U(MX,NY),V(MX,NY)

      CHARACTER*48 LAB3
      CHARACTER*48 TIMELABEL
      CHARACTER*48 LAB1

        PARAMETER (LRWK=50000,LIWK=50000,
     $             LMAP=500000,NWRK=50000,NOGRPS=5)
        REAL RWRK(LRWK), XWRK(NWRK), YWRK(NWRK)
        INTEGER IWRK(LIWK)
        INTEGER MAP(LMAP),IAREA(NOGRPS),IGRP(NOGRPS)



        REAL RWRK(LRWK)
        INTEGER IWRK(LIWK),map(liwk)

        REAL SCALE
        DATA SCALE /1.1/

      parameter (machine_epsilon_p = 1.19e-07)  ! from iftran.im - bj
      COMMON /SAVMAP/
     .  MX1, MX2, MY1, MY2,
     .  U1,  U2,  V1,  V2
      COMMON /FXFY1/ XA, YA, UA, VA, U1A, V1A, DUDI, DVDJ

        COMMON /NXNY/ NX1,NY1
      external fill
      external color
c ... compare,spv def'n from iftran.im/b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      spv(a) = compare(a,SPVAL_P,0.01)


        SPVAL_P=99999.

        FMIN0 = FMIN
        FMAX0 = FMAX

        NX1 = NX
        NY1 = NY
        NX2X = NX*2-1
        NY2X = NY*2-1

        MAXCHR = 'X'
        MINCHR = 'N'

C
C     CALC MAX AND MIN OF F
C
      FMAX = -1.E8
      FMIN =  1.E8
      DO 55 J=1,NY
        DO 60 I=1,NX
          IF (SPV(F(I,J)) .EQ. NO_P)
     .    THEN
            IF (F(I,J).GT.FMAX) FMAX=F(I,J)
            IF (F(I,J).LT.FMIN) FMIN=F(I,J)
          END IF
60      CONTINUE
55    CONTINUE
C
C     TEST FOR ZERO FINT
C
      IF (FINT .LE. MACHINE_EPSILON_P) THEN
        DF = .1*(FMAX-FMIN)
      ELSE
        DF = FINT
      END IF

      WRITE(6,5) FMIN,FMAX,DF
5     FORMAT(/,' MIN = ',G12.5,5X,'MAX = ',
     .  G12.5,5X,' PLOT INTERVAL = ',G12.5)
C
C    PLOT IF THERE IS ADEQUATE FIELD VARIATION
C

        write(LAB3,101)  FMAX,FMIN,DF
  101   FORMAT(2X,'MAX=',F10.3,2X,'MIN=',F10.3,2X,'INT=',F8.3,2X)

C Set up color table
c       CALL COLOR
C Initialize Areas
        CALL ARINAM(MAP,LMAP)



        CALL SET (0.,1.,0.,1.,0.,1.,0.,1.,1)

        CALL PCSETI ('CC - CHARACTER COLOR', 27+NC)
        IF (NC.EQ.1) CALL PLCHHQ (0.0,0.88,TIMELABEL,0.018,0.,-1.)
        CALL PLCHHQ (0.45,0.14-(NC-1)*0.02, LAB3   , .01,
     1       0., -1.)
        CALL PLCHHQ (0.05,0.14-(NC-1)*0.02, LAB1   , .01,
     1       0., -1.)

        CALL SET(MX1,MX2,MY1,MY2,U1,U2,
     .    V1,V2,1)                          ! SUPMAP SET
        FMIN = CEILING(FMIN,DF,MACHINE_EPSILON_P)   ! Minimum contour
        IOFFP=1
        SPVALU=SPVAL_P
        IOFFM=1            ! CONREC parameters
        CALL FXFY(NF)                               ! Initialize FX and FY fcns


        NULBLL = 1         ! number of unlabelled lines between labelled ones
        NHI = -1
        IF (IVAR.EQ.8) NHI=0

       IF (ISM.EQ.1) THEN
C --- CALL SET FOR 2X HIGH-RESOLUTION GRID

        CALL SET(MX1,MX2,MY1,MY2,
     .    U1,U1+(NX2X-1)*DUDI,
     .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
C --- Interpolate to 2x resolution grid.
        CALL SUBGRD  (F,MX,NX,NY,F1,N2X,NX2X,NY2X)
C --- Light smoother on high-resolution grid
        CALL SMOOTH2 (F1,HOLD,N2X,NX2X,NY2X, 0.5)
        CALL SMOOTH2 (F1,HOLD,N2X,NX2X,NY2X,-0.5)
       END IF
C --- Do contouring

        call cpseti ('CLS - contour level selection flag',20)
        call cpsetr ('CIS',DF)
        call cpseti ('LIS - label interval specifier',2)
c       call cpsetr ('CMN',FMIN0)
c       call cpsetr ('CMX',FMAX0)

       if (ihl.ge.1) then
        call CPSETR('HLS  - HIGH/LOW LABEL SIZE',.020)
        call CPSETC('ILT',' ')

         print *,' How many significant digits do you want in label?'
         accept *, isd
        call cpseti('NSD',isd)
        call cpseti('NLS',isd)
        CALL CPSETI ('HLC - HIGH/LOW LABEL COLOR INDEX', 27+NC)
        if (ihl.eq.2) call cpsetc('LOT',' ')
        if (ihl.eq.3) call cpsetc('HIT',' ')
       end if

       IF (ISM.EQ.1) THEN
        CALL CPRECT (F1,N2X,NX2X,NY2X,RWRK,LRWK,IWRK,LIWK)
cx       CALL CPPKCL (F1,RWRK,IWRK)

       call cplbam (f1,rwrk,iwrk,map)
       if (ihl.ge.1) then
        call cplbdr (f1,rwrk,lwrk)
       end if
        call cpclam (f1,rwrk,iwrk,map)
C Set fill style to solid, and fill contours
        CALL GSFAIS(1)
        call arscam (map,xwrk,ywrk,nwrk,iarea,igrp,nogrps,fill)
cx      call cpcldm (f1,rwrk,iwrk,map,cpdrpl)
        call cpback (f1,rwrk,iwrk)
       ELSE
        CALL CPRECT (F,NX,NX,NY,RWRK,LRWK,IWRK,LIWK)
        call cplbam (f,rwrk,iwrk,map)
        CALL CPPKCL (F,RWRK,IWRK)
       if (ihl.ge.1) then
        call cplbdr (f,rwrk,lwrk)
       end if
       END IF


        CALL CPGETI ('NCL - NUMBER OF CONTOUR LEVELS', NCON)
        DO 111 I=1,NCON
          CALL CPSETI ('PAI - PARAMETER ARRAY INDEX', I)
          call cpgetr ('CLV - contour level values',cval)
          if (cval.lt.0.) then
            call cpseti ('CLD',61166)
          else
            call cpseti ('CLD',65535)
          end if
cx          CALL CPSETI ('CLC - CONTOUR LINE COLOR INDEX', 27+NC)
111     CONTINUE

       if (ihl.ge.1) then
        call CPSETR('HLS  - HIGH/LOW LABEL SIZE',.025)
        call CPSETC('ILT',' ')
        call cpseti('NSD',2)
        call cpseti('NLS',2)
        CALL CPSETI ('HLC - HIGH/LOW LABEL COLOR INDEX', 27+NC)
        if (ihl.eq.2) call cpsetc('LOT',' ')
        if (ihl.eq.3) call cpsetc('HIT',' ')
       end if


       IF (ISM.EQ.1) THEN
c       CALL SET(MX1,MX2,MY1,MY2,
c    .    U1,U1+(NX2X-1)*DUDI,
c    .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
c       CALL CONREC (F1,N2X,NX2X,NY2X,FMIN,FMAX,DF,1,NHI,-'1470'O) ! Plot contour
        CALL CPCLDR (F1,RWRK,IWRK)
       ELSE
c       CALL SET(MX1,MX2,MY1,MY2,
c    .    U1,U1+(NX2X-1)*DUDI,
c    .    V1,V1+(NY2X-1)*DVDJ,1)                          ! SUPMAP SET
c       CALL CONREC (F1,N2X,NX2X,NY2X,FMIN,FMAX,DF,1,NHI,-'1470'O) ! Plot contour
        CALL CPCLDR (F,RWRK,IWRK)
       END IF


        IF (IVAR.EQ.5) THEN
        print *,' enter scale (suggested value = 1.5)'
        accept *,scale
        print *,' enter increment for plotting analyzed winds (suggested
     1 = 2)'
        accept *,iskip
        print *,' ENTER MIN SPEED FOR PLOTTING WINDS'
        ACCEPT *,SPDMIN
      CALL GETSET (MX1,MX2,MY1,MY2,
     1  DUMMY,DUMMY,DUMMY,DUMMY,IDUMMY)
      CALL SET (MX1,MX2,MY1,MY2,1.
     1  ,FLOAT(NX),1.,FLOAT(NY),1)

          DO J = 2,NY,iskip
          DO I = 2,NX,iskip
            CALL BARBS_R(FLOAT(I),FLOAT(J),U(I,J),V(I,J),SCALE,
     1  SPDMIN)
          END DO
          END DO
        END IF

        call gflas3(NF)

cisize=7
cCALL PWRIT (MX1,MY2,' ',1,0,0,0)
c       CALL FRAME                                  ! Flush frame buffer
        WRITE(6,155)
155     FORMAT(' THE FIELD HAS BEEN PLOTTED')
      RETURN
      END

      SUBROUTINE FXFY(IBUF)
C
C PURPOSE       Initialize parameters for FX and FY functions
C
C      INCLUDE 'CONFIG1'
c     INCLUDE 'CONFIG3'
C
      COMMON /NXNY/ NX,NY
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
      COMMON /SAVMAP/
     .  MX1, MX2, MY1, MY2,
     .  UU1, UU2, VV1, VV2
        SPVAL_P=99999.
      X = SPVAL_P
      Y = SPVAL_P
      U = SPVAL_P
      V = SPVAL_P
      U1 = UU1
      DUDI = (UU2-UU1)/(NX-1.)
      V1 = VV1
      DVDJ = (VV2-VV1)/(NY-1.)
      RETURN
      END
      FUNCTION FX(XIN,YIN)
C
C PURPOSE       Transform grid indices (XIN,YIN) to an x-coordinate in the
C               SUPMAP (U,V) coordinate space.  This routine is used by
C               CONREC and VELVCT.
C
C REMARKS       Subroutine FXFY must be called prior to the first call to this
C               function.
C
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
c ... compare, machine_epsilon from iftran.im / b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      IF (COMPARE(XIN,X,1.19e-07).EQ.NO_P) THEN
        X = XIN
        U = U1 + (XIN-1.)*DUDI
      END IF
      FX = U
      RETURN
      END
      FUNCTION FY(XIN,YIN)
C
C PURPOSE       Transform grid indices (XIN,YIN) to a y-coordinate in the
C               SUPMAP (U,V) coordinate space.  This routine is used by
C               CONREC and VELVCT.
C
C REMARKS       Subroutine FXFY must be called prior to the first call to this
C               function.
C
      COMMON /FXFY1/ X, Y, U, V, U1, V1, DUDI, DVDJ
c ... compare, machine_epsilon from iftran.im / b. jewett
      compare(a,b,c)=((int(sign(1.,c-abs(a-b)))+1)/2)
      IF (COMPARE(YIN,Y,1.19e-07).EQ.NO_P) THEN
        Y = YIN
        V = V1 + (YIN-1.)*DVDJ
      END IF
      FY = V
      RETURN
      END
cccccccccccccccccccccccccccccc  ceiling  cccccccccccccccccccccccccccc
c
c  ceiling - ceiling of x, over interval y, to accuracy z
c
        real function ceiling(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0) then                ! ceiling(x,y,z) = trunc(x+y-z,y)
          ceiling=trunc(x+y-z,y)
        else if (x.lt.0) then           ! ceiling(x,y,z) = -trunc(-x,y)
          ceiling= -trunc(-x,y)
        else                            ! ceiling(x,y,z) = 0.
          ceiling = 0.
        endif
c
        return
        end
c
ccccccccccccccccccccccccccccc  floor  cccccccccccccccccccccccccccc
c
c  floor - floor of x, over interval y, to accuracy z
c
        real function floor(x,y,z)
        implicit none
        real x,y,z,trunc,a,b
        trunc(a,b) = (b*int(a/b))
c
        if (x.gt.0.) then               ! floor(x,y,z) = trunc(x,y)
          floor = trunc(x,y)
        else if (x.lt.0.) then          ! floor(x,y,z) = -trunc(-x+y-z,y)
          floor = -trunc(-x+y-z,y)
        else                            ! floor(x,y,z) = 0.
          floor = 0.
        endif
c
        return
        end

        SUBROUTINE FILL (XWRK,YWRK,NWRK,IAREA,IGRP,NGRPS)
C
        DIMENSION XWRK(*),YWRK(*),IAREA(*),IGRP(*)

        DO 10, I=1,NGRPS
          IF (IGRP(I).EQ.3) IAREA3=IAREA(I)
 10     CONTINUE

        IF (IAREA3 .GT. 0) THEN
C If the area is defined by 3 or more points, fill it
           CALL GSFACI(IAREA3+1)
           CALL GFA(NWRK,XWRK,YWRK)
        ENDIF
      return
      end

      SUBROUTINE COLOR
      print *, ' White or Black background. Enter 1/2'
      accept 5, iwhite
5     format (i)

C    BACKGROUND COLOR
C  The background is white here for better visibility on paper

      if (iwhite.eq.1) then
        CALL GSCR (1,0,1.,1.,1.)
        CALL GSCR (1,1,0.8,0.8,1.)
      else
        CALL GSCR (1,0,0.,0.,0.)
        CALL GSCR (1,1,0.,0.,0.)
      end if

C
C     BACKGROUND COLOR
C     BLACK
c     CALL GSCR(1,0,0.,0.,0.)
C
C     FORGROUND COLORS
C White
      CALL GSCR(1,  1, 1.0, 1.0, 1.0)
        CALL GSCR (1,1,0.8,0.8,1.)
C Orchid
c     CALL GSCR(1,  2, 0.85, 0.45, 0.8)
C white
      CALL GSCR(1,  2, 1.0, 1.0, 1.0)
c     CALL GSCR(1,  2, 0.85, 0.45, 0.8)
C Red
      CALL GSCR(1,  3, 0.9, 0.25, 0.0)
C OrangeRed
      CALL GSCR(1,  4, 1.0, 0.0, 0.2)
C Orange
      CALL GSCR(1,  5, 1.0, 0.65, 0.0)
C Gold
      CALL GSCR(1,  6, 1.0, 0.85, 0.0)
C Yellow
      CALL GSCR(1,  7, 1.0, 1.0, 0.0)
C GreenYellow
      CALL GSCR(1,  8, 0.7, 1.0, 0.2)
C Chartreuse
      CALL GSCR(1,  9, 0.5, 1.0, 0.0)
C Celeste
      CALL GSCR(1, 10, 0.2, 1.0, 0.5)
C Green
      CALL GSCR(1, 11, 0.2, 0.8, 0.2)
C Aqua
      CALL GSCR(1, 12, 0.0, 0.9, 1.0)
C DeepSkyBlue
      CALL GSCR(1, 13, 0.0, 0.75, 1.0)
C RoyalBlue
      CALL GSCR(1, 14, 0.25, 0.45, 0.95)
C SlateBlue
      CALL GSCR(1, 15, 0.4, 0.35, 0.8)
C DarkViolet
      CALL GSCR(1, 16, 0.6, 0.0, 0.8)
C Lavender
      CALL GSCR(1, 17, 0.8, 0.8, 1.0)
      CALL GSCR (1,21,0.,0.,0.)
      CALL GSCR (1,22,0.,.7,0.)
      CALL GSCR (1,23,.0,.0,.7)
      CALL GSCR (1,24,.0,.7,.0)
      CALL GSCR (1,25,.3,.3,.5)
      CALL GSCR (1,26,.3,.3,.5)
      CALL GSCR (1,27,.3,.3,.5)
      CALL GSCR (1,28,.0,.0,.95)
      CALL GSCR (1,29,.95,.0,.0)
      CALL GSCR (1,30,.0,.95,.0)
      CALL GSCR (1,31,.2,.7,.7)

C Done.
C
        RETURN
C
      END
