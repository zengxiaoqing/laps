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
      SUBROUTINE SUPMAP_local 
     +                (JPROJ,POLAT,POLON,RROT,PL1,PL2,PL3,PL4,JJLTS,
     +                   JGRID,JUS,JDOT,IER)
C   SupMap is used to plot map outlines according to one of nine projections.
C The origin and orientation of the projection are selected by the user.
C Points on the Earth defined by latitude and longitude are transformed to
C points in the u,v plane, the plane of projection.  The u and v axes are
C respectively parallel to the x and y axes of the plotter.  A rectangular
C frame parallel to the u and v axes is chosen and only material within the
C frame is plotted.
C   Continental and U.S. state and county outlines are available for plotting.
C In addition, a western U.S. subset of the counties has been extracted, and
C up to three user files may be read.

C            Name          Date                 Description
C       --------------  -----------     ---------------------------------------
C       ?? (NCAR)          JAN 69       "Revised"
C                          MAY 71       "Revised"
C                          OCT 73       "Standardized"
C                          JUL 74       "Revised"
C                          AUG 76       "Revised"
C       M. Cairns       27 Jan 81       Added MapCol common for intensities
C       J. Wakefield    19 Jun 81       Added MapDas common for dash patterns
C                       29 Jun 81       Removed trailing blanks and modified
C                                       to read unformatted files.
C                       28 Sep 81       Added common block SupMpB.
C                       23 Mar 82       Added negative JDot function, to allow
C                                       suppression of SupMap call printing if
C                                       JUS=0.  See comments below.
C                       16 Feb 83       Added SysDisk to data file names and
C                                       added error code for open failure.
C                       10 Apr 84       Added user files option.
C                       14 Jan 85       Fixed bug in Part < 0 option.
C       J. Ramer        13 Nov 85       Added common blocks necessary for using
C                                       MapInv_GC with all projections and for
C                                       recreating a map from the contents of
C                                       the common blocks.  Added subroutine
C                                       RESUP which recreates the map.
C       J. Wakefield    19 Nov 86       Chgd SysDisk to Lib_Dev for files.
c
c       P. McDonald     17 Aug 93       Fixed things up to run on Stardent
c                                       and (hopefully) other UNIX machines.

C       Usage
C                       Call SupMap (JProj,PoLat,PoLong,RRot,PL1,PL2,PL3,PL4,
C                                    JJLTS,JGrid,JUS,JDot,IER)

C       Dimension of    PL1(2),PL2(2),PL3(2),PL4(2)
C       arguments

C       On input        JProj
C       for SupMap        |JProj| defines the projection type
C                         according to the following code:
C                               1  Stereographic
C                               2  Orthographic
C                               3  Lambert Conformal Conic with two standard
C                                  parallels
C                               4  Lambert Equal Area
C                               5  Gnomonic
C                               6  Azimuthal Equidistant
C                               7  Dummy -- this code is not used
C                               8  Cylindrical Equidistant
C                               9  Mercator
C                              10  Mollweide type
C                         If JProj < 0, the map and grid lines are omitted.

C                       PoLat,PoLong,RRot
C                         If (|JProj|.ne.3)
C                         . PoLat and PoLong define in degrees the latitude and
C                           longitude of the point on the globe which is to
C                           transform to the origin of the u,v plane.
C                               -90 .le. PoLat .le. 90
C                              -180 .le. PoLong .le. 180
C                           degrees of latitude north of the Equator and
C                           degrees of longitude east of the Greenwich Meridian
C                           are positive.  If the origin is at the North Pole,
C                           "north" is considered to be in the direction of
C                           (PoLong+180.). If the origin is at the South Pole,
C                           "north" is in the direction of PoLong.
C                         . RRot is the angle between the v axis and north at
C                           the origin.  It is measured in degrees and is taken
C                           to be positive if the angular movement from north
C                           to the v axis is counter-clockwise.  For the
C                           cylindrical projections (8,9,10), the axis of the
C                           projection is parallel to the v axis.
C                         If (|JProj|.eq.3) (Lambert Conformal Conic)
C                         . PoLong = central meridian of projection in degrees.
C                         . PoLat,RRot are the two standard parallels in deg.

C                       JJLTS,PL1,PL2,PL3,PL4
C                         |JJLTS| can take the values 1 through 5 and specifies
C                         one of five options on the way in which the limits of
C                         the rectangular map are defined by the parameters
C                         PL1, PL2, PL3, and PL4.

C                         |JJLTS| = 1
C                           The maximum useful area produced by the projection
C                           is plotted.  PL1, PL2, PL3, and PL4 are not used
C                           and may be set to zero.

C                         |JJLTS| = 2
C                           In this case (PL1,PL2) and (PL3,PL4) are the
C                           latitudes and longitudes in degrees of two points
C                           which are to be at opposite corners of the map
C                           (upper right and lower left, respectively).
C                           Care must be taken when using cylindrical
C                           projections and this option.

C                         |JJLTS| = 3
C                           The minimum and maximum values of u and v are
C                           specified by PL1 through PL4.  PL1 = UMin,
C                           PL2 = UMax, PL3 = VMin, PL4 = VMax.  Knowledge of
C                           the transformation equations is necessary for this
C                           option to be used (see below).

C                         |JJLTS| = 4
C                           Here PL1 = AUMin, PL2 = AUMax, PL3 = AVMin,
C                           PL4 = AVMax, where
C                              AUMin = Angular distance from origin to left
C                                      frame of map.
C                              AUMax = Angular distance from origin to right
C                                      frame of map.
C                              AVMin = Angular distance from origin to lower
C                                      frame.
C                              AVMax = Angular distance from origin to upper
C                                      frame.
C                           AUMin, AUMax, AVMin, AVMax must be positive and the
C                           origin must be within the rectangular limits of the
C                           map.  This option is useful for polar projections.
C                           It is not appropriate for the Lambert Conformal
C                           with two standard parallels.  An error message is
C                           printed if an attempt is made to use JJLTS = 4 when
C                           JProj = 3, (see below).

C                         |JJLTS| = 5
C                           PL1 through PL4 are two element arrays giving the
C                           latitudes and longitudes of four points which are
C                           to be on the four sides of the rectangular frame.
C                           PL1(1), PL1(2) are respectively the latitude and
C                           longitude of a point on the left frame.  Similarly,
C                           PL2 lies on the right frame, PL3 lies on the lower
C                           frame and PL4 lies on the upper frame.  Note that
C                           in the calling program PL1 through PL4 will be
C                           dimensioned:
C                               Real  PL1(2),PL2(2),PL3(2),PL4(2)

C                          .If JJLTS is positive, the SupMap call is plotted
C                           below the map.  This is omitted if JJLTS is < 0.

C                       JGrid
C                         |JGrid| gives in degrees the interval at which lines
C                         of latitude and longitude are to be plotted.  A value
C                         in the range 1 through 10 will usually be appropriate
C                         but higher values are acceptable.
C                         If JGrid < 0 the border around the map is omitted.
C                         If JGrid = 0 no grid lines are plotted.

C                       JUS
C                         |JUS|
C                         1  World outlines
C                         2  U.S. state outlines
C                         4  U.S. counties and states
C                         8  PROFS-region counties and states (western U.S.)
C                  16,32,64  SupMap user file(s) 1,2,3
C                         To access a combination of datasets, use the sum of
C                         their values (e.g. 10 = U.S. + PROFS area counties).
C                         .If JUS is positive, the SupMap call and values of
C                          of UMin, UMax, VMin, VMax are printed as an aid to
C                          debugging.  This is omitted if JUS is negative.

C                       JDot
C                         |JDot|
C                         0  for continuous outlines.
C                         1  for dotted outlines.
C                         .If JDot is negative, the SupMap call is neither
C                          printed nor plotted.

C       On output       All arguments except IER are unchanged.
C       for SupMap

C                       IER
C                         Error flag with the following meanings:
C                         If IER =
C                          0  Map successfully plotted.
C                         29  Error opening data file.
C                         33  Attempt to use non-existent projection.
C                         34  Map limits inappropriate.
C                         35  Angular limits too great.
C                         36  Map has zero area.

C       Entry points    MaPlot, SupCon, SupFst, SupMap, SupTrp, SupVec, QCon,
C                       QVec, VecPlt, ReSup

C                       MaPlot
C                         Actually draws the map.

C                       SupCon
C                         Once the transformation has been set up by an initial
C                         call to SupMap, the subroutine SupCon may be called
C                         to transform a point, (latitude, longitude) to the
C                         corresponding point, (u, v) on the plane.  Contours
C                         may thus be readily drawn against the map background.
C                         (See SupFst and SupVec below).

C                            Call SupCon(RLat,RLon,U,V)

C                         On input:
C                           RLat,RLon are the latitude and longitude of a point
C                           to be transformed to the u,v plane.
C                            -90. .le. RLat .le.  90.
C                           -180. .le. RLon .le. 180.

C                         On output:
C                           RLat,RLon are unchanged.
C                           U,V are the transformed coordinates of the point.

C                       QCon
C                         Actually performs the above mentioned transformation.

C                       SupFst
C                       SupVec
C                         To facilitate drawing lines on the map these routines
C                         which act like the plotting routines FrstPt and
C                         Vector are included.  They are subject to the same
C                         restrictions as SupCon above.

C                            Call SupFst (RLat,RLon)
C                            Call SupVec (RLat,RLon)

C                       QVec
C                         Decides what lines are to be drawn and where.

C                       SupTrp
C                         Performs interpolation to the edges of the frame.

C                       VecPlt
C                         Called by QVec to draw (dot) lines on the plotter.

C       Common blocks    Name   Length
C                       SupMp1  9 FP
C                       SupMp2  1 Int + 204 FP
C                       SupMp3  2 Int + 5 FP
C                       SupMp4  5 Int + 2 FP
C                       SupMp5  1 Int + 5 FP
C                       SupMp6  6 FP
C                       SupMp7  3 Int + 2 FP
C                       SupMp8  6 FP
C                       SupMp9  3 FP
C                       SupMpA  1 Int
C                       SupMpB  4 FP
C                       SupMpC  4 Int + 1 FP
C                       SupMpD  9 FP
C                       MapCol  4 Int
C                       MapDas  4 Int

C       I/O             Map plotted.  Outline data is read from any of several
C                       disk files.  SupMap call printed.

C       Precision       Single

C       Language        FORTRAN

C       Algorithm       The latitudes and longitudes of successive outline
C                       points are transformed to coordinates in the plane of
C                       projection and joined by a vector.

C       References      Hershey, A. V., The plotting of maps on a CRT printer.
C                         NWL Report No. 1844, 1963.
C                       Lee, Tso-Hwa, Students Summary Reports, Work-study
C                         program in scientific computing.  NCAR 1968.
C                       Parker, R. L., UCSD SuperMap:  World Plotting Package.
C                       Steers, J.A., An Introduction to the Study of Map
C                         Projections.  Univ. of London Press, 1962.

C       Accuracy        The definition of the map produced is limited by the
C                       fact that the resolution of the virtual plotter space
C                       is 1024 units in the x and y directions.

C       Plotting routines       PWRT, FrstPt, Vector, Point, DashLn, Perim, Set
C       used

C       Required resident       ATan, Tan, Sin, Cos, ALog, SqRt, ATan2, ACos
C       routines

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMP9/DS,DI,DSRDI
        COMMON/SUPMPA/IIER
        COMMON/SUPMPB/X1,Y1,X2,Y2
        COMMON/SUPMPC/LPROJ,ROT,JLTS,LGRID,IUS
        COMMON/SUPMPD/AAA,BBB,CCC,DDD,EEE,FFF,GGG,HHH,III

      DIMENSION PL1(2),PL2(2),PL3(2),PL4(2),LABA(20),LABB(18)

      REAL LAT1,LAT2,AAA,BBB,CCC,DDD,EEE,FFF,GGG,HHH,III

      EQUIVALENCE (PLA1,AUMIN),(PLA2,AUMAX),(PLA3,AVMIN),(PLA4,AVMAX),
     1            (PHIA,LAT1),(ROT,LAT2)

      DATA PLTRES,RESLIM/ 1024.,10./

!     EXTERNAL SUPMBD

      SQU(X) = X*X
c
c       SUPMBD subroutine changed to BLOCK DATA
c     IGDFLT = %LOC(SUPMBD)             !FORCE LOAD OF BLOCK DATA from library

      write(6,*)' Subroutine supmap...',jproj

      ROT = RROT
      ILTS = IABS(JJLTS)
      JLTS = JJLTS
      LGRID = JGRID
      IGRID = IABS(LGRID)
c     JGR = ISIGN(1,LGRID)
      JGR = LGRID
      IUS = JUS
      IUSGN = ISIGN(1,IUS)
      PLA1 = PL1(1)
      PLA2 = PL2(1)
      PLA3 = PL3(1)
      PLA4 = PL4(1)
      LABL = 77
      IDOT = IABS(JDOT)

C INITIALIZATION
      LPROJ = JPROJ
      IPROJ = IABS(LPROJ)
      JPR = ISIGN(1,LPROJ)
C     JGR=JPR
      IOUT = ABS(IUS)           !USE IOUT FOR CHOOSING DATA FILE(S)
      PHIA = POLAT
      POLONG = POLON
      PHIO = POLON
      PHIOC = 540.-PHIO
      ICROSS = 0
      IIER = 0
      ILF = 0

C COMPUTE CONSTANTS APPROPRIATE TO EACH PROJECTION
      IF (IPROJ .NE. 3) GO TO 30

C LAMBERT CONFORMAL CONIC
      SGN = SIGN(1.,0.5*(LAT1+LAT2))
      CHI1 = (90.-SGN*LAT1)*DTR
      IF (LAT1 .EQ. LAT2) GO TO 20
      CHI2 = (90.-SGN*LAT2)*DTR
      CONE = ALOG(SIN(CHI1)/SIN(CHI2))/ALOG(TAN(0.5*CHI1)/TAN(0.5*CHI2))
      GO TO 60
   20 CONE = COS(CHI1)
      GO TO 60

C THE OTHERS
   30 XT = ROT*DTR
      SINR = SIN(XT)
      COSR = COS(XT)
      XT = PHIA*DTR
      SINO = SIN(XT)
      COSO = COS(XT)

C  CALCULATE COEFFICIENTS NECESSARY TO CONVERT DISPLACEMENT VECTOR ON EARTH'S
C  SURFACE FROM CARTESIAN SYSTEM WITH ORIGIN AT CENTER OF EARTH, UNIT VECTOR K
C  K AT THE PROJECTION CENTER, AND UNIT VECTOR J IN THE POSITIVE U DIRECTION TO
C  A SYSTEM WITH UNIT VECTOR K AT NORTH POLE AND UNIT VECTOR I AT INTERSECTION
C  OF EQUATOR AND GREENWICH MERIDIAN.  NEEDED LATER FOR MAPINV_GC.
      AAA=COS(ROT*DTR)*COS(POLONG*DTR)*SIN(PHIA*DTR)-
     &    SIN(ROT*DTR)*SIN(POLONG*DTR)
      BBB=COS(ROT*DTR)*SIN(POLONG*DTR)*SIN(PHIA*DTR)+
     &    SIN(ROT*DTR)*COS(POLONG*DTR)
      CCC=-COS(ROT*DTR)*COS(PHIA*DTR)
      DDD=-SIN(ROT*DTR)*COS(POLONG*DTR)*SIN(PHIA*DTR)-
     &    COS(ROT*DTR)*SIN(POLONG*DTR)
      EEE=-SIN(ROT*DTR)*SIN(POLONG*DTR)*SIN(PHIA*DTR)+
     &    COS(ROT*DTR)*COS(POLONG*DTR)
      FFF=SIN(ROT*DTR)*COS(PHIA*DTR)
      GGG=COS(POLONG*DTR)*COS(PHIA*DTR)
      HHH=SIN(POLONG*DTR)*COS(PHIA*DTR)
      III=SIN(PHIA*DTR)

      IF (IPROJ-7) 60,67,40

C CYLINDRICAL PROJECTIONS               [8,9,10]
   40 IF (PHIA .NE. 0.0) GO TO 42
      IF (ROT .EQ. 0.0) GO TO 45
      IF (ABS(ROT) .EQ. 180.) GO TO 50
   42 SINO1 = COSO*COSR
      COSO1 = SQRT(CON1-SINO1*SINO1)
      OVC1 = 1./COSO1
      PHIO = PHIO-ATAN2(SINR*OVC1,-COSR*SINO*OVC1)*RTD
      PHIOC = 540.-PHIO
      SINR = SINR*COSO*OVC1
      COSR = -SINO*OVC1
      SINO = SINO1
      COSO = COSO1
      GO TO 60

C USE SIMPLE TRANSFORMS FOR CYLINDRICAL PROJECTIONS IF ROT = POLAT = .0
C   I.E. IPROJ = 11, 12, 13
   45 SINO = 1.0
      IPROJ = IPROJ+3
      GO TO 55

   50 SINO = -1.0
      PHIO = PHIO+180.
      PHIOC = PHIOC+180.
   55 COSO = 0.0
      SINR = 0.0
      COSR = 1.0
      ILF = 1

C ILTS = 1         THE MAXIMUM USEFUL AREA IS PLOTTED.
C ---------
   60 GO TO (61 ,62 ,62 ,61 ,61 ,66 ,67 ,68 ,66 ,70 ,
     1       68 ,66 ,70 ),IPROJ

C STEREOGRAPHIC                         [1]
C LAMBERT EQUAL AREA                    [4]
C GNOMONIC                              [5]
   61 UMIN = -2.0
      VMIN = -2.0
      UMAX = 2.0
      VMAX = 2.0
      GO TO 80

C ORTHOGRAPHIC                          [2]
C LAMBERT CONFORMAL CONIC               [3]
   62 UMIN = -1.0
      VMIN = -1.0
      UMAX = 1.0
      VMAX = 1.0
      GO TO 80

C AZIMUTHAL EQUIDISTANT                 [6]
C MERCATOR WITH ARBITRARY POLE          [9]
C MERCATOR                              [12]
   66 UMAX = PI
      VMAX = PI
      UMIN = -PI
      VMIN = -PI
      GO TO 80

C DUMMY  --  ERROR EXIT                 [7]
   67 IIER = 33
      CALL ULIBER2 (IIER
     1           ,' SUPMAP-ATTEMPT TO USE NON-EXISTENT PROJECTION')     
      GO TO 700

C CYLINDRICAL EQUIDISTANT               [8,11]
   68 UMAX = 180.
      UMIN = -180.
      VMAX = 90.
      VMIN = -90.
      GO TO 80

C MOLLWEIDE TYPE                        [10,13]
   70 UMAX = 2.0
      UMIN = -2.0
      VMAX = 1.0
      VMIN = -1.0

   80 UEPS = 0.5*(UMAX-UMIN)
      VEPS = 0.5*(VMAX-VMIN)
      IF (IPROJ .EQ. 3) UEPS = 180.

C COMPUTE THE APPROPRIATE MAP BOUNDARIES.
      GO TO (600,200,300,400,500),ILTS

C ILTS = 2         POINT (PL1,PL2) IN UPPER RIGHT CORNER , (PL3,PL4) IN
C ---------        LOWER LEFT CORNER OF PLOT.
  200 RLAT = PLA1
      RLON = PLA2
      CALL QCON
      U1 = U
      V1 = V
      RLAT = PLA3
      RLON = PLA4
      CALL QCON
      UMAX = AMAX1(U1,U)
      UMIN = AMIN1(U1,U)
      VMAX = AMAX1(V1,V)
      VMIN = AMIN1(V1,V)
      GO TO 600

C ILTS = 3         SET PLOT LIMITS DIRECTLY.
C ----------
  300 UMAX = PLA2
      UMIN = PLA1
      VMAX = PLA4
      VMIN = PLA3
      GO TO 600

C ILTS = 4         USE ANGULAR DISTANCES TO SET PLOT LIMITS.
C ----------
  400 COSUMI = COS(AUMIN*DTR)
      SINUMI = SQRT(CON1-COSUMI*COSUMI)
      COSUMA = COS(AUMAX*DTR)
      SINUMA = SQRT(CON1-COSUMA*COSUMA)
      COSVMI = COS(AVMIN*DTR)
      SINVMI = SQRT(CON1-COSVMI*COSVMI)
      COSVMA = COS(AVMAX*DTR)
      SINVMA = SQRT(CON1-COSVMA*COSVMA)

      GO TO (401,402,403,404,405,406,407,408,409,410,
     1       408,409,410),IPROJ

C STEREOGRAPHIC                         [1]
  401 UMAX = (1.-COSUMA)/SINUMA
      UMIN = -(1.-COSUMI)/SINUMI
      VMAX = (1.-COSVMA)/SINVMA
      VMIN = -(1.-COSVMI)/SINVMI
      GO TO 600

C ORTHOGRAPHIC                          [2]
  402 IF (AMAX1(AUMIN,AUMAX,AVMIN,AVMAX) .GT. 90.) GO TO 900
      UMAX = SINUMA
      UMIN = -SINUMI
      VMAX = SINVMA
      VMIN = -SINVMI
      GO TO 600

C LAMBERT CONFORMAL CONIC               [3]
  403 IIER = 34
      CALL ULIBER2 (IIER,' SUPMAP-MAP LIMITS INAPPROPRIATE')
      GO TO 700

C LAMBERT EQUAL AREA                    [4]
  404 UMAX = (1.+COSUMA)/SINUMA
      UMIN = (1.+COSUMI)/SINUMI
      VMAX = (1.+COSVMA)/SINVMA
      VMIN = (1.+COSVMI)/SINVMI
      UMAX = 2./SQRT(1.+UMAX*UMAX)
      UMIN = -2./SQRT(1.+UMIN*UMIN)
      VMAX = 2./SQRT(1.+VMAX*VMAX)
      VMIN = -2./SQRT(1.+VMIN*VMIN)
      GO TO 600

C GNOMONIC                              [5]
  405 IF (AMAX1(AUMIN,AUMAX,AVMIN,AVMAX) .GE. 90.) GO TO 900
      UMAX = SINUMA/COSUMA
      UMIN = -SINUMI/COSUMI
      VMAX = SINVMA/COSVMA
      VMIN = -SINVMI/COSVMI
      GO TO 600

C AZIMUTHAL EQUIDISTANT                 [6]
  406 UMAX = AUMAX*DTR
      UMIN = -AUMIN*DTR
      VMAX = AVMAX*DTR
      VMIN = -AVMIN*DTR
      GO TO 600

C DUMMY  --  ERROR EXIT                 [7]
  407 GO TO 67

C CYLINDRICAL EQUIDISTANT               [8,11]
  408 UMAX = AUMAX
      UMIN = -AUMIN
      VMAX = AVMAX
      VMIN = -AVMIN
      GO TO 600

C MERCATOR                              [9,12]
  409 IF (AMAX1(AVMIN,AVMAX) .GE. 90.) GO TO 900
      UMAX = AUMAX*DTR
      UMIN = -AUMIN*DTR
      VMAX = ALOG((1.+SINVMA)/COSVMA)
      VMIN = -ALOG((1.+SINVMI)/COSVMI)
      GO TO 600

C MOLLWEIDE TYPE                        [10,13]
  410 UMAX = AUMAX*OV90
      UMIN = -AUMIN*OV90
      VMAX = SINVMA
      VMIN = -SINVMI
      GO TO 600

C ILTS = 5         USE FOUR EDGE POINTS TO SET LIMITS.
C ----------
  500 PLB1 = PL1(2)
      RLAT = PLA1
      RLON = PLB1+EPS
      CALL QCON
      UMIN = U
      PLB2 = PL2(2)
      RLAT = PLA2
      RLON = PLB2-EPS
      CALL QCON
      UMAX = U
      PLB3 = PL3(2)
      RLAT = PLA3
      RLON = PLB3
      CALL QCON
      VMIN = V
      PLB4 = PL4(2)
      RLAT = PLA4
      RLON = PLB4
      CALL QCON
      VMAX = V

C COMPUTE MAP LIMITS FOR PLOT
  600 DU = UMAX-UMIN
      DV = VMAX-VMIN

C ERROR IF MAP HAS ZERO AREA
      IF (DU.EQ.0.0 .OR. DV.EQ.0.0) GO TO 905
        IF(PART.LE..0)GOTO 620          !USER-SUPPLIED X1,Y1,X2,Y2
      IF (DU .GT. DV) GO TO 610
      Y1 = 0.5*(1.-PART)
      Y2 = 1.-Y1
      X1 = 0.5*(1.-PART*DU/DV)
      X2 = 1.-X1
      GO TO 620

  610 X1 = 0.5*(1.-PART)
      X2 = 1.-X1
      Y1 = 0.5*(1.-PART*DV/DU)
      Y2 = 1.-Y1

C ERROR IF MAP HAS ESSENTIALLY ZERO AREA
  620 IF (AMIN1(X2-X1,Y2-Y1)*PLTRES .LT. RESLIM) GO TO 905
      CALL SET (X1,X2,Y1,Y2,UMIN,UMAX,VMIN,VMAX,1)
      call line(umin,vmin,umin,vmax)
      call line(umin,vmax,umax,vmax)
      call line(umax,vmax,umax,vmin)
      call line(umax,vmin,umin,vmin)
      DS = SQU(((X2-X1)*PLTRES)/DU)
      DSRDI = SQRT(DS/DI)
      write(6,*)' DSRDI = ',DSRDI

C DO WE WRITE ANYTHING?
      IF (JLTS.LT.0 .AND. IUSGN.LT.0) GO TO 640
        IF(JDOT.LT.0)GOTO 640

C CREATE THE LABEL
c     ENCODE (LABL,7000,LABA(1)) LPROJ,PHIA,POLONG,ROT,PLA1,PLA2,PLA3,
c    1                           PLA4,JLTS,LGRID,IUS,IDOT
c7000 FORMAT (8H SUPMAP(,I3,7(1H,,F6.1),4(1H,,I3),1H))

c     IF (ILTS .EQ. 5) ENCODE (61,7010,LABB(1)) PLB1,PLB2,PLB3,PLB4
c7010 FORMAT (7X,1H(,24X,4(1H,,F6.1),1H))

      IF (JLTS .LT. 0) GO TO 630

C WRITE SUPMAP CALL BENEATH THE MAP
      CALL PWRT (240,17,LABA(1),LABL,0,0)
      IF (ILTS .EQ. 5) CALL PWRT (240,1,LABB(1),61,0,0)

  630 IF (IUSGN .LT. 0) GO TO 640

C PRINT OUT THE CALL ET AL.
      K = (LABL+3)/4
      WRITE (6,6000) (LABA(I),I=1,K)
 6000 FORMAT (30A4)
      IF (ILTS .EQ. 5) WRITE (6,6000) (LABB(I),I=1,18)
      WRITE (6,6010) UMIN,UMAX,VMIN,VMAX
 6010 FORMAT (8H UMIN = ,F11.6,9H  UMAX = ,F11.6,9H  VMIN = ,
     +        F11.6,9H  VMAX = ,F11.6)

C DRAW THE MAP
  640 IF (IOUT.NE.0 .OR. IGRID.NE.0 .OR. JGR.GE.0) CALL MAPLOT_local
      IDOT = 0

C RETURN IER
  700 IER = IIER
      RETURN

C ERROR RETURNS
  900 IIER = 35
      CALL ULIBER2 (IIER,' SUPMAP-ANGULAR LIMITS TOO GREAT')
      GO TO 700
  905 IIER = 36
      CALL ULIBER2 (IIER,' SUPMAP-MAP HAS ZERO AREA')
      GO TO 700

      END
C-------------------------------------------------------------------------------
      SUBROUTINE MAPLOT_local
C THIS SUBROUTINE PLOTS THE CONTINENTAL AND U.S. STATE OUTLINES,
C MERIDIANS, PARALLELS, LIMBS WHERE APPROPRIATE. IT LABELS KEY MERIDIANS
C AND POLES, AND IT DRAWS A BORDER.

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP2/NPTS,MAXLAT,MINLAT,MAXLON,MINLON,PTS(1200000)
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMPA/IIER
        COMMON/MAPCOL/MPCOL1,MPCOL2,MPCOL3,MPCOL4
        COMMON/MAPDAS/LDASH1,LDASH2,LDASH3,LDASH4
c       DATA                    !Default line intensities and dash patterns
c    1          MPCOL1,LDash1   /255,'1777'O/,  !Map lines
c    2          MPCOL2,LDash2   /128,'1756'O/,  !Grid lines
c    3          MPCOL3,LDash3   /192,'1777'O/,  !Limb lines
c    4          MPCOL4,LDash4   /255,'1777'O/   !Perimeter

        character       supmap_dir*150
        integer       lsdir

        data iwrite /0/

        character*1000 c_line

      DIMENSION SPLAT(2)
      REAL MAXLAT,MINLAT,MAXLON,MINLON,MIDLAT,MIDLON
      CHARACTER*180 NAMFIL
      DATA   SINLMB,COSLMB /0.017452406, 0.99984765/
      DATA FLOORC / 10000. /

      FLOOR(X) = AINT(X+FLOORC)-FLOORC
      CLING(X) = FLOOR(X)+1.


      call get_directory('static',supmap_dir,lsdir)

      supmap_dir = supmap_dir(1:lsdir)//'/ncarg/'
c        supmap_dir = '../static/ncarg/'
      lsdir = index (supmap_dir, ' ') - 1

      CALL DASHLN(LDASH1)
        CALL OPTN(2HIN,MPCOL1)
      GO TO (10,10,10,5,10,5,905,5,5,5,5,5,5),IPROJ
    5 ICF=1
   10 IF(IOUT.EQ.0)GOTO 100

C***Select appropriate file(s) according to bits set in IOut
        IBS=1                           !START WITH 1ST BIT
        MASK=1
   11   DO 12 IBIT=IBS,16               !CHECK LOWER 16 BITS
         IBS=IBS+1
c        IF((IOUT.AND.MASK).NE.0)GOTO 13        !EXIT LOOP
         IF(iand(IOUT,MASK).NE.0)GOTO 13        !EXIT LOOP
         MASK=MASK*2
   12   CONTINUE
        GO TO 100                       !FELL OUT OF LOOP SO DONE
   13   MASK=MASK*2

        if(ibit .eq. 1)goto14
        if(ibit .eq. 2)goto15
        if(ibit .eq. 3)goto16
        if(ibit .eq. 4)goto17
        if(ibit .eq. 5)goto18
        if(ibit .eq. 6)goto19
        if(ibit .eq. 7)goto20
        GO TO 100                       !OUT OF RANGE -- RETURN

c  14   IF((IOUT.AND.2).NE.0)THEN       ! STATES WITH CONTINENTS?
   14   IF(iand(IOUT,2).NE.0)THEN       ! STATES WITH CONTINENTS?
         IBS=IBS+1      !SKIP BIT LATER
         MASK=MASK*2
c        NAMFIL='Lib_Dev:[GUDAT]CONANDSTA.DAT'  !3 (SPECIAL CASE)
c       namfil = '/home/star1/b/mcdonald/data/supmap/conandsta.dat'
        namfil = supmap_dir (1:lsdir) // 'conandsta.dat'
        ELSE
c        NAMFIL='Lib_Dev:[GUDAT]CONTINENT.DAT'  !1
c       namfil = '/home/star1/b/mcdonald/data/supmap/continent.dat'
        namfil = supmap_dir (1:lsdir) // 'continent_minus_us.dat'
        ENDIF
        GO TO 25
c  15   NAMFIL='Lib_Dev:[GUDAT]STATE.DAT'       !2
c  15 namfil = '/home/star1/b/mcdonald/data/supmap/state.dat'
   15   namfil = supmap_dir (1:lsdir) // 'state.dat'
        GO TO 25
c  16   NAMFIL='Lib_Dev:[GUDAT]USCOUNTY.DAT'    !4
c  16 namfil = '/home/star1/b/mcdonald/data/supmap/uscounty.dat'
   16   namfil = supmap_dir (1:lsdir) // 'uscounty.dat'
        GO TO 25
c  17   NAMFIL='Lib_Dev:[GUDAT]COUNTY.DAT'      !8
c  17 namfil = '/home/star1/b/mcdonald/data/supmap/county.dat'
   17   namfil = supmap_dir (1:lsdir) // 'state_from_counties.dat'
        GO TO 25
c  18   NamFil='SupMap_UserFile1'               !16
   18   namfil = supmap_dir (1:lsdir) // 'userfile1.dat'
        GoTo 25
c  19   NamFil='SupMap_UserFile2'               !32
   19   namfil = supmap_dir (1:lsdir) // 'userfile2.dat'
        GoTo 25
c  20   NamFil='SupMap_UserFile3'               !64
   20   namfil = supmap_dir (1:lsdir) // 'userfile3.dat'
        GoTo 25

C***Open file
c  25   Open(3,Name=NamFil,Type='Old',Form='UnFormatted',ReadOnly,Err=26)
   25   continue
        write(6,*)Namfil(1:60)
        Open(3,file=NamFil,status='Old',Form='UnFormatted',err=26)
        GoTo 30

C***Error opening file -- return
   26   IIER=29
        write(6,*)'Supmap - Error opening: ',NamFil
        Return

C***Read next line
   30   continue

        if(.true.)then
            if(iwrite .eq. 0)
     1        write(6,*)' Using simple read to read binary map info...'
            Read(3,End=99)NPts,MaxLat,MinLat,MaxLon,MinLon
     1                   ,(Pts(M),M=1,NPts)

        else ! This may be needed for DEC Alpha but will not work on LINUX
            if(iwrite .eq. 0)
     1        write(6,*)' Using cio.c to read binary map info...'
            Read(3,End=99,err=41)NPts,MaxLat,MinLat,MaxLon,MinLon
     1                   ,(Pts(M),M=1,200)
 41         continue

!           Convert between Bigendian and Littleendian (under construction)
            call in_to_im(4,4,NPts,1)
            call in_to_im(4,4,Pts,NPts)
            call in_to_im(4,4,MaxLat,1)
            call in_to_im(4,4,MinLat,1)
            call in_to_im(4,4,MaxLon,1)
            call in_to_im(4,4,MinLon,1)

        endif

        iwrite = iwrite + 1

   99   NPts=NPts/2
      IF (NPTS .EQ. 0)THEN
        CLOSE(3)
        GOTO 11                 !CHECK NEXT BIT
      ENDIF
      IF (NPTS .LE. 16) GO TO 70
      IF (ICF .NE. 0) GO TO 70

C DOES THIS LINE INTERSECT THE SCREEN?
C       1  2  3
C       4     5
C       6  7  8
      MIDLAT = (MAXLAT+MINLAT)*0.5
      MIDLON = (MAXLON+MINLON)*0.5
      RLAT = MAXLAT
      RLON = MAXLON
      CALL QCON
      X3 = U
      Y3 = V
      RLON = MIDLON
      CALL QCON
      X2 = U
      Y2 = V
      RLON = MINLON
      CALL QCON
      X1 = U
      Y1 = V
      RLAT = MIDLAT
      CALL QCON
      X4 = U
      Y4 = V
      RLON = MAXLON
      CALL QCON
      X5 = U
      Y5 = V
      RLAT = MINLAT
      CALL QCON
      X8 = U
      Y8 = V
      RLON = MIDLON
      CALL QCON
      X7 = U
      Y7 = V
      RLON = MINLON
      CALL QCON
      X6 = U
      Y6 = V
      XMN = AMIN1(X1,X2,X3,X4,X5,X6,X7,X8)
      XMX = AMAX1(X1,X2,X3,X4,X5,X6,X7,X8)
      YMN = AMIN1(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)
      YMX = AMAX1(Y1,Y2,Y3,Y4,Y5,Y6,Y7,Y8)

      DX = AMIN1(XMX-XMN,180.)
      DY = AMIN1(YMX-YMN,180.)
      XMX = XMX+0.01*DX
      XMN = XMN-0.01*DX
      YMX = YMX+0.01*DY
      YMN = YMN-0.01*DY

      IF (XMN.GT.UMAX .OR. XMX.LT.UMIN .OR. YMN.GT.VMAX .OR.
     1    YMX.LT.VMIN) GO TO 30

   70 RLAT = PTS(1)
      RLON = PTS(2)
      IFST = 1
      IGO = 0
      CALL QVEC
      DO 75 J=2,NPTS
         RLAT = PTS(2*J-1)
         RLON = PTS(2*J)
         CALL QVEC
   75 CONTINUE

      GO TO 30

C***DONE PLOTTING -- DRAW AND LABEL GRID, IF DESIRED
  100 SPLAT(2) = 90.
      SPLAT(1) = -90.
      IF (IGRID .EQ. 0) GO TO 300
      IDOTS = IDOT
      IDOT = 0

c     if (idot .eq. 0) go to 200

C LETTER KEY MERIDIANS AND POLES

C       SOUTH POLE
      ISPF = 0
      RLAT = -90.
      RLON = 0.0
      CALL QCON
      IF ((U .GT. UMAX) .OR. (U .LT. UMIN) .OR. (V .GT. VMAX) .OR.
     1    (V .LT. VMIN)) GO TO 110
      USP = U
      VSP = V
      ISPF = 1
      IPF = 1
      If (phia .LT. 0) CALL PWRT (U,V,2HSP,2,1,0)

C       NORTH POLE
110   INPF = 0
      IPF = 0
      RLAT = 90.
      CALL QCON
      IF ((U .GT. UMAX) .OR. (U .LT. UMIN) .OR. (V .GT. VMAX) .OR.
     1    (V .LT. VMIN)) GO TO 120
      UNP = U
      VNP = V
      INPF = 1
      IPF = 1
      If (phia .GT. 0) CALL PWRT (U,V,2HNP,2,1,0)

C       EQUATOR
  120 RLON = PHIO-10.
      RLAT = 0.0
      DO 125 I=1,36
         RLON = RLON+10.
         CALL QCON
         IF (U.LE.UMAX .AND. U.GE.UMIN .AND. V.LE.VMAX .AND. V.GE.VMIN)
     1       GO TO 130
  125 CONTINUE
      GO TO 140

  130 CALL PWRT (U,V,2HEQ,2,1,0)

C       GREENWICH MERIDIAN
  140 RLAT = 85.
      RLON = 0.0
      DO 145 I=1,16
         RLAT = RLAT-10.
         CALL QCON
         IF (U.LE.UMAX .AND. U.GE.UMIN .AND. V.LE.VMAX .AND. V.GE.VMIN)
     1       GO TO 150
  145 CONTINUE
      GO TO 160

  150 CALL PWRT (U,V,2HGM,2,1,0)

C       DATE LINE
  160 RLAT = 85.
      RLON = 180.
      DO 165 I=1,16
         RLAT = RLAT-10.
         CALL QCON
         IF (U.LE.UMAX .AND. U.GE.UMIN .AND. V.LE.VMAX .AND. V.GE.VMIN)
     1       GO TO 170
  165 CONTINUE
      GO TO 200

  170 CALL PWRT (U,V,1HI,1,1,0)

  200 RGRID = IGRID
c      CALL DASHLN (LDASH2)
c       CALL OPTN(2HIN,MPCOL2)

C SHOULD WE BOTHER LIMITING GRID POINTS TRANSFORMED?
      IF (ICF .NE. 0) GO TO 270
      IF (IPROJ.GE.8 .AND. IPROJ.LE.10) GO TO 270

C SET UP TO FIND EXTREMA
      DLON = RGRID
      STLON = FLOOR(POLONG/RGRID)*RGRID
      IF (ISPF.NE.0 .AND. INPF.EQ.0) STLON = STLON+180.
      RLON = STLON-DLON
      SPLON = STLON+360.
      J = 0
      PSIGN = 1.

C CHECK FOR SOUTH POLE
      IF (ISPF .NE. 0) PSIGN = -1.

C DO WE GRID POLES SPECIALLY?
      SPLAT(2) = 90.*PSIGN
      SPLAT(1) = SPLAT(2)

C IF BOTH POLES WITHIN FRAME JUMP.
      IF (INPF.NE.0 .AND. ISPF.NE.0) GO TO 270

C IF EITHER IN FRAME USE AS BASE
      IF (INPF.NE.0 .OR. ISPF.NE.0) GO TO 230

C NO POLE IS CLOSE TO THE WINDOW
      J = -1
      SPLAT(2) = FLOOR(PHIA/RGRID)*RGRID
      IF (ABS(SPLAT(2)) .EQ. 90.) SPLAT(2) = 0.0

C SEARCH FOR FIRST POINT WITHIN FRAME.
  210 RLON = RLON+DLON
      DLAT = RGRID
      RLAT = SPLAT(2)-DLAT
  215 RLAT = RLAT+DLAT
      CALL QCON
      IF ((U .LE. UMAX) .AND. (U .GE. UMIN) .AND. (V .LE. VMAX) .AND.
     1    (V .GE. VMIN)) GO TO 225
      IF (ABS(RLAT) .LT. 90.) GO TO 215
      IF (DLAT .LT. 0.0) GO TO 220

C REVERSE LATITUDE SEARCH DIRECTION
      RLAT = SPLAT(2)+DLAT
      DLAT = -DLAT
      GO TO 215

C UPDATE LONGITUDE ! QUIT.
  220 J = 0
      IF (RLON-SPLON) 210,300,300

C SET UP FOR LIMIT SEARCH
  225 J = J+1
      STLON = RLON
      RLON = STLON-DLON
      IF (RLAT .EQ. 0.0) RLAT = SIGN(RLAT,-PSIGN)
      SPLAT(2) = RLAT
      SPLAT(1) = SPLAT(2)

      splat(1) = 90.0
      splat(2) = -90.0
      stlon = stlon - dlon
      splon = splon + dlon
      rlon  = stlon
  226 if (rlon .le. splon) then
          rlat2 = 80.0
          if (amod(rlon,90.0) .eq. 0.0) rlat2 = 90.0
          rlat  = -rlat2
          ifst  = 1
          igo   = 0
          call qvec
  227     rlat = rlat + 1.0
          if (rlat .le. rlat2) then
              if (iproj .eq. 1) rlat = rlat2
              call qvec
              if (igo .ne. 0) then
                  if (rlat .lt. splat(1)) splat(1) = rlat
                  if (rlat .gt. splat(2)) splat(2) = rlat
              endif
              go to 227
          endif
          rlon  = rlon + dlon
          go to 226
      endif
      if (dlon .ne. 0.0) go to 285

C LONGITUDE LOOP
C       IGF     FLAG TO SIGNAL NO POINTS WITHIN WINDOW.
C       IPF     FLAG SIGNALS WHETHER A POLE LIES WITHIN THE FRAME.
C       ILF     FLAG SIGNALS WHETHER TO PLOT COMPLETE LONGITUDES
C               (I.E. TO POLE FOR ALL LATITUDES.)
  230 RLON = RLON+DLON
      IF (RLON.GE.SPLON .OR. RLON.LT.STLON) GO TO 285
      I1 = IPF
      I2 = MOD(I1+1,2)
      TSA = PSIGN
      DLAT = -PSIGN
      DX = AMOD(90.,RGRID)
      IF (DX .EQ. 0.0) DX = RGRID
      XLAT = 90.-DX
      IF (ILF.NE.0 .OR. AMOD(RLON,90.).EQ.0.0) XLAT = 90.
      OLAT = SIGN(AMIN1(ABS(SPLAT(I2+1)),XLAT),SPLAT(I2+1))
      IGF = 0
  235 IFST = 1
      IGO = 0
      RLAT = OLAT
      CALL QVEC

C LATITUDE LOOP.
  240 RLAT = RLAT+DLAT
      IGF = MAX0(IGO,IGF)
      CALL QVEC
      IF (IGO .NE. 0) GO TO 245

C THIS POINT OUTSIDE THE FRAME
      IF (RLAT*TSA .LE. SPLAT(I1+1)*TSA) GO TO 250
  245 IF (ABS(RLAT) .LT. XLAT) GO TO 240
      RLAT = SIGN(AMAX1(ABS(SPLAT(I1+1)),XLAT),SPLAT(I1+1))

C POSSIBLE NEW LATITUDE EXTREME.
  250 SPLAT(I1+1) = RLAT

C REVERSE LATITUDE SEARCH DIRECTION
      I1 = I2
      I2 = MOD(I1+1,2)
      TSA = -PSIGN
      DLAT = PSIGN
      IF (I1 .NE. 0) GO TO 235

C LATITUDE LOOP FINISHED.
      IF (ABS(SPLAT(I2+1)) .LT. 90.) GO TO 255
      IPF = 1
      PSIGN = SIGN(1.,SPLAT(I2+1))
      SPLAT(I2+1) = SPLAT(I1+1)
      SPLAT(I1+1) = 90.*PSIGN
  255 IF (IGF .NE. 0) GO TO 230

C LONGITUDE EXTREME REACHED.
      IF (J .NE. 0) GO TO 260

C CHANGE LONGITUDE DIRECTION.
      J = 1
      SPLON = RLON
      RLON = STLON
      DLON = -DLON
      STLON = SPLON-360.
      GO TO 230

C SET UP LAST LONGITUDE EXTREME
  260 IF (DLON .LT. 0.0) GO TO 265
      SPLON = RLON
      GO TO 285
  265 STLON = RLON
      GO TO 285

C DRAW ALL MERIDIANS.
  270 DLON = RGRID
      STLON = 0.0
      SPLON = 360.
      RLON = 0.0
      SPLAT(2) = 90.
      SPLAT(1) = -90.
      DX = AMOD(90.,RGRID)
      IF (DX .EQ. 0.0) DX = RGRID
      OLAT = 90.-DX

  275 RLON = RLON+DLON
      IGO = 0
      IFST = 1
      XLAT = OLAT

      IF (ILF.NE.0 .OR. AMOD(RLON,90.).EQ.0.0) then
         If (phia.LT.0) then
            rlat = 0.0
            xlat = 90.0

         Else if (phia.GE.0) then
            rlat = 90.0
            xlat = 0.0

         End if

      Else
         rlat = xlat

      End if

      CALL QVEC
  280 RLAT = RLAT-1.
      CALL QVEC
      IF (RLAT .GT. -XLAT) GO TO 280
      IF (RLON .LT. SPLON) GO TO 275

C DRAW PARALLELS
  285 DLAT = RGRID
      RLAT = AMIN1(SPLAT(2),SPLAT(1))
      OLAT = AMAX1(SPLAT(2),SPLAT(1))
      SPLAT(2) = FLOOR(RLAT/RGRID)*RGRID
      SPLAT(1) = AMIN1(CLING(OLAT/RGRID)*RGRID,90.)
      RLAT = AMAX1(DLAT-90.,SPLAT(2))-DLAT
      OLAT = AMIN1(90.-DLAT,SPLAT(1))
c      if (dlon .gt. 0.0) then
c          stlon = stlon - dlon
c      else
c          stlon = stlon + dlon
c      endif
      DLON = 1.
      IF (ILF .NE. 0) IPF = 0

  290 RLAT = RLAT+DLAT
      IF (IPF .NE. 0) DLON = 1./COS(DTR*RLAT)
      IGO = 0
      IFST = 1
      RLON = STLON
      CALL QVEC

  295 RLON = RLON+DLON
      CALL QVEC
      IF (RLON .LE. SPLON) GO TO 295
      IF (RLAT .LT. OLAT) GO TO 290

      IDOT = IDOTS
c      CALL DASHLN (LDASH3)
c       CALL OPTN(2HIN,MPCOL3)

C DRAW LIMB LINES
  300 IDOTS = IDOT
      IDOT = 0
      GO TO (400,330,305,335,400,340,400,400,400,345,
     1       400,400,345),IPROJ

C LAMBERT CONFORMAL CONIC           [3]
  305 continue
      go to 400 ! test to eliminate spurious line
      DLAT = 1.
      RLON = PHIO+CON2
      OLAT = AMAX1(-90.,SPLAT(2)-DLAT)
      K = CLING(SPLAT(1)-SPLAT(2))
      DO 320 I=1,2
         IGO = 0
         IFST = 1
         RLAT = OLAT
         CALL QVEC
         DO 310 J=1,K
            RLAT = RLAT+DLAT
            CALL QVEC
  310    CONTINUE
         RLON = PHIO-CON2
  320 CONTINUE
      GO TO 400

C ORTHOGRAPHIC                  [2]
  330 RADIUS = 1.
      AXIS = 1.
      GO TO 350

C LAMBERT EQUAL AREA            [4]
  335 RADIUS = 2.
      AXIS = 1.
      GO TO 350

C AZIMUTHAL EQUDISTANT          [6]
  340 RADIUS = PI
      AXIS = 1.
      GO TO 350

C MOLLWEIDE                     [10,13]
  345 RADIUS = 2.
      AXIS = 0.5

  350 U = RADIUS
      V = 0.0
      W = 0.0
      IGO = 0
      IFST = 1
      DO 370 I=1,361
         V = AXIS*V
         IF (U.LE.UMAX .AND. U.GE.UMIN .AND. V.LE.VMAX .AND. V.GE.VMIN)
     1       GO TO 355
         IGO = 0
         GO TO 365
  355    IF (IGO .NE. 0) GO TO 360
         CALL FRSTPT (U,V)
         IGO = 1
         GO TO 365

  360    CALL VECTOR (U,V)
  365    V = U*SINLMB+W*COSLMB
         U = U*COSLMB-W*SINLMB
         W = V
  370 CONTINUE

C DRAW BORDER
  400 IF (JGR .GT. 0)THEN
c       CALL DASHLN(LDASH4)
c       CALL OPTN(2HIN,MPCOL4)
c       CALL PERIM(1,1,1,1)
c        call frstpt (umin,vmin)
c        call vector (umax,vmin)
c        call vector (umax,vmax)
c        call vector (umin,vmax)
c        call vector (umin,vmin)
      ENDIF
      IDOT = IDOTS
      RETURN

  905 RETURN
      END

C-------------------------------------------------------------------------------
      SUBROUTINE QCON
C THIS SUBROUTINE TRANSFORMS THE POINT (RLAT,RLON), IN DEGREES,
C TO (U,V) ON THE MAP PLANE DEPENDENT UPON THE PROJECTION, IPROJ.

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMPA/IIER
c
c
      common /supmp7/   phio,phia,igrid,idot,ilts
c
      data      resl    / 2.9765625 /
c
c
      DATA OLDU,OLDV / 0.,0./

      U = AMOD(RLON+PHIOC,360.)-180.

      GO TO (50 ,50 ,130,50 ,50 ,50 ,170,50 ,50 ,50 ,
     1       210,220,230),IPROJ

   50 T1 = U*DTR
      T2 = RLAT*DTR
      SINPH = SIN(T1)
      SINLA = SIN(T2)
      COSPH = COS(T1)
      COSLA = SQRT(1.-SINLA*SINLA)
      TCOS = COSLA*COSPH
      COSA = SINLA*SINO+TCOS*COSO
      SINA = SQRT(CON1-COSA*COSA)

C PATCH TO AVOID DIVIDE BY ZERO
      OVSINA = 1.E12
      IF(SINA.NE.0.0)OVSINA = 1./SINA
C END PATCH

      SINB = COSLA*SINPH*OVSINA
      COSB = (SINLA*COSO-TCOS*SINO)*OVSINA

C PERFORM TRANSFORMATION APPROPRIATE TO THE PROJECTION
      GO TO (110,120,130,140,150,160,170,180,190,200),IPROJ

C STEREOGRAPHIC                         [1]
  110 R = (1.-COSA)*OVSINA
        a = atan2 (sina, cosa) * 0.5
        r = tan (a)
      go to 300
c 110 ihemi = 1
c     if (phia .lt. 0.0) ihemi = -1
c     call maproj_poster (rlat, rlon, u, v, resl, polong, ihemi, 1)
c     v     = 8193 - v
c     go to 305

C ORTHOGRAPHIC                          [2]
  120 R = SINA
      IF (COSA) 320,320,300

C LAMBERT CONFORMAL CONIC               [3]
  130 UDIF = ABS(U-OLDU)
      OLDU = U
      CHI = 90.-SGN*RLAT
      IF (CHI .GE. CON2) GO TO 320
      R = TAN(0.5*DTR*CHI)**CONE
      U = U*CONE*DTR
      V = -R*SGN*COS(U)
      U = R*SIN(U)
      GO TO 310

C LAMBERT EQUAL AREA                    [4]
  140 IF (ABS(COSA+1.) .LT. 1.E-6) GO TO 320
      R = (1.+COSA)*OVSINA
      R = 2./SQRT(1.+R*R)
      GO TO 300

C GNOMONIC                              [5]
  150 IF (COSA .LE. 0.0) GO TO 320
      R = SINA/COSA
      GO TO 300

C AZIMUTHAL EQUIDIDSANT                 [6]
  160 IF (ABS(COSA+1.) .LT. 1.E-6) GO TO 320
      R = ACOS(COSA)
      GO TO 300

C DUMMY   --  ERROR                     [7]
  170 IIER = 33
      CALL ULIBER2 (IIER,
     1             ' SUPMAP-ATTEMPT TO USE NON-EXISTENT PROJECTION')
      GO TO 320

C CYLINDRICAL EQUIDISTANT,  ARBITRARY POLE AND ORIENTATION.
  180 IF (ABS(1.-COSA*COSA) .LT. 1.E-4) GO TO 320
      U = ATAN2(SINB*COSR+COSB*SINR,SINB*SINR-COSB*COSR)*RTD
      V = 90.-ACOS(COSA)*RTD
      GO TO 305

C MERCATOR, ARBITRARY POLE AND ORIENTATION.
  190 IF ((1.-COSA*COSA) .LT. 2.E-6) GO TO 320
      U = ATAN2(SINB*COSR+COSB*SINR,SINB*SINR-COSB*COSR)
      V = ALOG((1.+COSA)*OVSINA)
      GO TO 305

C MOLLWEIDE, ARBITRARY POLE AND ORIENTATION.
  200 IF (ABS(1.-COSA*COSA) .LT. 2.E-6) GO TO 320
      U = ATAN2(SINB*COSR+COSB*SINR,SINB*SINR-COSB*COSR)*TOVPI
      UDIF = ABS(U-OLDU)
      OLDU = U
      V = COSA
C     U = U*SQRT(1.-V*V)
      U = U*SQRT(ABS(1.-V*V))
      GO TO 310

C CYLINDRICAL EQUIDISTANT FOR POLAT = ROT = 0.  [11]
  210 V = RLAT
      GO TO 305

C MERCATOR                              [12]
  220 U = U*DTR
      V = ALOG(TAN(0.00872664*(RLAT+90.0001)))
      GO TO 305

C MOLLWEIDE                             [13]
  230 U = U*OV90
      V = SIN(RLAT*DTR)
      UDIF = ABS(U-OLDU)
      OLDU = U
      U = U*SQRT(1.-V*V)
      GO TO 310

C TERMINAL PHASE    [1,2,4,5,6]
  300 U = R*(SINB*COSR+COSB*SINR)
      V = R*(COSB*COSR-SINB*SINR)

C CHECK FOR CROSSOVER
  305 UDIF = ABS(U-OLDU)
      OLDU = U
  310 VDIF = ABS(V-OLDV)
      OLDV = V
      ICROSS = 0
      IF (UDIF.GT.UEPS .OR. VDIF.GT.VEPS) ICROSS = 1
      RETURN

C DISPENSE WITH UNDEFINED POINTS
  320 U = 1.E12
      ICROSS = 0
      IF (ABS(U-OLDU) .GT. UEPS) ICROSS = 1
      OLDU = U
      RETURN

      END
C-------------------------------------------------------------------------------
      SUBROUTINE QVEC
C THIS SUBROUTINE TRANSFORMS AND PLOTS LINE SEGMENTS FOR SUPMAP AND OTHERS

C INPUTS (PASSED THROUGH COMMON.)

C  (RLAT,RLON)  NEXT POINT TO BE PLOTTED
C  IFST         A FLAG USED TO SIGNAL THE FIRST POINT OF A LINE SEGMENT
C               = 0  -  START A NEW LINE
C               = 1  -  CONTINUATION OF A LINE

C OTHER VARIABLES

C  (U,V)        NEXT POINT TRANSFORMED TO THE VIRTUAL SCREEN BY SUPCONQ
C  ICROSS       A FLAG RETURNED BY SUPCONQ FOR CYLINDRICAL PROJECTIONS
C  IGO          = 0  -  LAST POINT NOT PLOTTED
C               = 1  -  LAST POINT WAS PLOTTED.
C  (U1,V1),(U2,V2)  PARAMETERS PASSED TO SUPTRP.

        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMP9/DS,DI,DSRDI

c     SQU(X) = (X)*(X)

C TRANSFORM THE POINT
      CALL QCON

C HAVE WE FLIPPED TO OTHER SIDE OF FRAME?
      IF (ICROSS .NE. 0) IGO = 0

C ARE WE WITHIN THE FRAME?
      IF (U.GT.UMAX .OR. U.LT.UMIN .OR. V.GT.VMAX .OR. V.LT.VMIN)
     1    GO TO 20
      IF (IGO .EQ. 0) GO TO 30

C CONTINUE LINE
C CHECK PROXIMITY TO PREVIOUS POINT.
c   5 IF ((SQU(U-UOLD)+SQU(V-VOLD))*DS .LE. DI) RETURN
    5 continue
      CALL VECPLT
   10 UOLD = U
      VOLD = V
      IGOLD = IGO
      RETURN

C THIS POINT LIES OUTSIDE THE FRAME
   20 IGO = 0
      IF (IFST .NE. 0) GO TO 65

      IF (IGOLD .EQ. 0) GO TO 10

C IT WAS INSIDE - INTERPOLATE TO EDGE OF FRAME
C STATUS OF LAST POINT.   IF NOT INSIDE FRAME, GO ON

C IF UNINTERPOLATABLE
      IF (ICROSS .NE. 0) GO TO 70

      U1 = UOLD
      V1 = VOLD
      U2 = U
      V2 = V
      CALL SUPTRP

C CHECK PROXIMITY TO PREVIOUS POINT.
c     IF ((SQU(U-UOLD)+SQU(V-VOLD))*DS .LE. DI) GO TO 25
      CALL VECPLT
   25 UOLD = U2
      VOLD = V2
      IGOLD = 0
      RETURN

C THIS POINT IS WITHIN THE FRAME

C IS IT THE FIRST POINT OF A LINE?
   30 IF (IFST .NE. 0) GO TO 60
      IF (IGOLD .EQ. 0) GO TO 50

C THE PREVIOUS POINT WAS INSIDE THE FRAME ON THE OTHER SIDE.
C START A NEW LINE
   40 CALL FRSTPT (U,V)
      IGO = 1
      GO TO 10

C LAST POINT NOT IN FRAME - THIS ONE IS
   50 IF (ICROSS .NE. 0) GO TO 40

C INTERPOLATE BACK TO EDGE
      U1 = U
      V1 = V
      U2 = UOLD
      V2 = VOLD
      CALL SUPTRP
      CALL FRSTPT (U,V)
      IGO = 1
      IGOLD = 1
      UOLD = U
      VOLD = V
      U = U1
      V = V1
      GO TO 5

C FIRST POINT ON LINE SEGMENT  -  CHECK FOR DUPLICATION OF END POINT
   60 IF (U.NE.UOLD .OR. V.NE.VOLD) CALL FRSTPT (U,V)
      IGO = 1
   65 IFST = 0
      GO TO 10

C IGNORE UNDEFINED POINT
   70 IFST = 1
      GO TO 10

      END
C-------------------------------------------------------------------------------
      SUBROUTINE SUPTRP
C THE INTERPOLATION ROUTINE.  FINDS (U,V) ON THE EDGE OF THE FRAME NEAREST
C (U1,V1).  (U1,V1) MUST LIE WITHIN THE FRAME, (U2,V2) WITHOUT.

        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2

      F(V) = (V-V2)*DU/DV+U2
      G(U) = (U-U2)*DV/DU+V2

C FIND INDEX TO (U2,V2)
C       5 | 4 | 6
C       --+---+--
C       2 | 1 | 3
C       --+---+--
C       8 | 7 | 9
      I = 1
      DU = U1-U2
      DV = V1-V2
      A = U2-UMIN
      B = U2-UMAX
      C = V2-VMIN
      D = V2-VMAX
      IF (A) 110,140,120
  110 I = I+1
      GO TO 140
  120 IF (B) 140,140,130
  130 I = I+2
  140 IF (C) 150,200,160
  150 I = I+6
      GO TO 200
  160 IF (D) 200,200,170
  170 I = I+3

  200 GO TO (900,210,220,230,240,250,260,270,280),I

  210 U = UMIN
      GO TO 300
  220 U = UMAX
      GO TO 300
  230 V = VMAX
      GO TO 350
  240 IF (F(VMAX)-UMIN) 210,230,230
  250 IF (F(VMAX)-UMAX) 230,230,220
  260 V = VMIN
      GO TO 350
  270 IF (F(VMIN)-UMIN) 210,260,260
  280 IF (F(VMIN)-UMAX) 260,260,220

C INTERPOLATE
  300 V = G(U)
      RETURN
  350 U = F(V)
      RETURN

C ERROR EXIT
  900 U = U2
      V = V2
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE VECPLT
C PLOTS THE LINE SEGMENT FROM (UOLD,VOLD) TO (U,V)
C INPUTS (PASSED THROUGH COMMON)
C  (UOLD,VOLD)  THE LAST POINT PLOTTED
C  (U,V)        THE NEXT POINT
C  IDOT         CONTROL FLAG  [ DOT VS PLOT ]

        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMP9/DS,DI,DSRDI

C DO WE DOT OR PLOT?
      IF (IDOT .NE. 0) GO TO 10

C PLOT
      CALL VECTOR (U,V)
      RETURN

C DOT
   10 DU = U-UOLD
      DV = V-VOLD

      call line(UOLD,VOLD,U,V)
      RETURN

      I = (ABS(DU)+ABS(DV))*DSRDI
      IF (I .LE. 1) GO TO 30
      A = 1./FLOAT(I)
      I = I-1
      DU = DU*A
      DV = DV*A
      UO = U
      VO = V
      U = UOLD
      V = VOLD
      DO 20 K=1,I
         U = U+DU
         V = V+DV
         CALL POINT (U,V)
   20 CONTINUE
      U = UO
      V = VO
   30 CALL POINT (U,V)
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE SUPVEC (XLAT,XLON)
C THIS SUBROUTINE ALLOWS THE USER TO DRAW LINES ON THE VIRTUAL SCREEN SET UP BY
C SUPMAP, UNENCUMBERED BY DECISIONS AS TO WHETHER IT WILL BE VISIBLE THROUGH
C THE WINDOW.
C USE SUPFST AND SUPVEC IN EXACTLY THE SAME MANNER AS FRSTPT AND VECTOR.

        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN

      RLAT = XLAT
      RLON = XLON
      CALL QVEC
      RETURN
      END
C-------------------------------------------------------------------------------
      SUBROUTINE SUPFST (XLAT,XLON)

        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD

      IGO = 0
      IFST = 1
      RLAT = XLAT
      RLON = XLON
      CALL QVEC
      RETURN

      END
C-------------------------------------------------------------------------------
      SUBROUTINE SUPCON (XLAT,XLON,XU,XV)
C THIS SUBROUTINE IS PROVIDED TO RETAIN COMPATIBILITY WITH USER'S PROGRAMS.

        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP8/U,V,U1,V1,U2,V2

      RLAT = XLAT
      RLON = XLON
      CALL QCON
      XU = U
      XV = V

      RETURN
      END
C-------------------------------------------------------------------------------
!       Moved to LAPS library (Steve Albers 1997)
!       BLOCK DATA
c     SUBROUTINE SUPMBD
C FORCE-LOAD BLOCK DATA
!       COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
!       COMMON/SUPMP2/NPTS,MAXLAT,MINLAT,MAXLON,MINLON,PTS(1000)
!       COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
!       COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
!       COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
!       COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
!       COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
!       COMMON/SUPMP8/U,V,U1,V1,U2,V2
!       COMMON/SUPMPA/IIER
!       COMMON/MAPCOL/MPCOL1,MPCOL2,MPCOL3,MPCOL4
!       COMMON/MAPDAS/LDASH1,LDASH2,LDASH3,LDASH4
!       DATA                    !Default line intensities and dash patterns
!    1          MPCOL1,LDash1   /255,'1777'O/,  !Map lines
!    2          MPCOL2,LDash2   /128,'1756'O/,  !Grid lines
!    3          MPCOL3,LDash3   /192,'1777'O/,  !Limb lines
!    4          MPCOL4,LDash4   /255,'1777'O/   !Perimeter


c       COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
c       COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
!       COMMON/SUPMP9/DS,DI,DSRDI

!     DATA   CON1 / 1.00001/
!     DATA   CON2 / 179.99999/
!     DATA     DI / 16./
!     DATA    DTR / 1.7453292519943E-2/
!     DATA    EPS / 1.E-6/
!     DATA   OV90 / 1.11111111111111E-2/
!     DATA     PI / 3.1415926535898/
!     DATA    RTD / 57.295779513082/
!     DATA  TOVPI / 0.63661977236758/
!     DATA   UOLD / 0.0 /
!     DATA   VOLD / 0.0 /
!     DATA PART/1.0/          !SIZE OF PICTURE (90% OF SCREEN)

c       RETURN
!     END
C-------------------------------------------------------------------------------
        SUBROUTINE RESUP(IER)
C RECALL SUPMAP USING CONTENTS OF COMMON BLOCKS.

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMPB/X1,Y1,X2,Y2
        COMMON/SUPMPC/LPROJ,ROT,JLTS,LGRID,IUS

        LPRJ=LPROJ
        PHA=PHIA
        PLONG=POLONG
        RT=ROT
        UMN=UMIN
        UMX=UMAX
        VMN=VMIN
        VMX=VMAX
        LGRD=LGRID
        IU=IUS
        IDT=IDOT
        PART=AMAX1(X2-X1,Y2-Y1)

        CALL SUPMAP(LPRJ,PHA,PLONG,RT,UMN,UMX,VMN,VMX,-3,LGRD,IU,IDT,IER
     1)

        RETURN

        END
c
c
c
        subroutine supset (jproj, polat, polon, rot,
     +                     pl1, pl2, pl3, pl4, jlts,
     +                     jgrid, jout, jdot, jerr)
c
c
        call supmap (jproj, polat, polon, rot,
     +               pl1, pl2, pl3, pl4, jlts,
     +               0, 0, 0, jerr)
c
c
        return
        end



c!!!!!!
c*******************************************************************************
c
c       The following routines are replacements for NCAR Graphics routines
c       that are called by SUPMAP.  They are included to supercede the
c       original routines.
c
c*******************************************************************************




        subroutine ichar_to_str (chars, nchars, str)

        integer       chars(*), nchars
        character       str*(*)


        integer       ic, ib, iw

        integer       iword, ichr
        character       ibyte(4)
        equivalence     (ibyte(1), iword)


        iw = 0
        ib = 4
        ic = 0
        do while (ic .lt. nchars)

            ib = ib + 1
            if (ib .gt. 4) then
                iw    = iw + 1
                iword = chars(iw)
                ib    = 1
            end if

            ic       = ic + 1
            ichr     = byte_to_i4(ibyte(ib))
            str(ic:) = char (ichr)

        end do

        return
        end




        subroutine dashln (idash_pat)

        integer       idash_pat


        return
        end




        subroutine optn (iopt, ival)

        integer       iopt, ival


        return
        end




        subroutine pwrt (u, v, chars, nchars, isize, iorient)

        real          u, v
        integer       chars(*), nchars, isize, iorient


        character       str*256


        call ichar_to_str (chars, nchars, str)

        write (*,'(2f10.5,3i5,5x,a)') u, v, isize, iorient, nchars, str

        return
        end




c       subroutine set (x1, x2, y1, y2, u1, u2, v1, v2, itype)
c
c       real          x1, x2, y1, y2, u1, u2, v1, v2
c       integer       itype
c
c
c       return
c       end




        subroutine uliber2 (ier, chars)

        integer       ier

        character       chars*(*)

        write(6,*)' Supmap error: ',ier,chars

        return
        end


        subroutine supmap_block_data()

c     routine supmap.f

        COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
        COMMON/SUPMP2/NPTS,MAXLAT,MINLAT,MAXLON,MINLON,PTS(1200000)
        COMMON/SUPMP3/POLONG,CONE,RLAT,RLON,JGR,ILF,SGN
        COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP5/PHIOC,SINO,COSO,SINR,COSR,IPROJ
        COMMON/SUPMP6/UMIN,UMAX,VMIN,VMAX,UEPS,VEPS
        COMMON/SUPMP7/PHIO,PHIA,IGRID,IDOT,ILTS
        COMMON/SUPMP8/U,V,U1,V1,U2,V2
        COMMON/SUPMPA/IIER
        COMMON/MAPCOL/MPCOL1,MPCOL2,MPCOL3,MPCOL4
        COMMON/MAPDAS/LDASH1,LDASH2,LDASH3,LDASH4

!       DATA                    !Default line intensities and dash patterns
!    1          MPCOL1,LDash1   /255,1023/,  !Map lines
!    2          MPCOL2,LDash2   /128,1006/,  !Grid lines
!    3          MPCOL3,LDash3   /192,1023/,  !Limb lines
!    4          MPCOL4,LDash4   /255,1023/   !Perimeter


c       COMMON/SUPMP1/PI,TOVPI,DTR,RTD,EPS,OV90,CON1,CON2,PART
c       COMMON/SUPMP4/IFST,IGO,IGOLD,ICROSS,IOUT,UOLD,VOLD
        COMMON/SUPMP9/DS,DI,DSRDI

        MPCOL1 = 255
        MPCOL2 = 128
        MPCOL3 = 192
        MPCOL4 = 255

        LDash1 = 1023 
        LDash2 = 1006
        LDash3 = 1023
        LDash4 = 1023

        CON1 = 1.00001
        CON2 = 179.99999
        DI = 16.
        DTR = 1.7453292519943E-2
        EPS = 1.E-6
        OV90 = 1.11111111111111E-2
        PI = 3.1415926535898
        RTD = 57.295779513082
        TOVPI = 0.63661977236758
        UOLD = 0.0 
        VOLD = 0.0 
        PART   = 1.0           !SIZE OF PICTURE (90% OF SCREEN)
  

        return
        end
