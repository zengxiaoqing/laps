!dis   
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis    
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis    
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis    
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis    
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis    
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis    
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis    
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis   
!dis 

MODULE constants

  ! Constants used by mm5post

      REAL,       PARAMETER        :: CM2INCH  =   1.0 / 2.54
      REAL,       PARAMETER        :: CP       =   1005.7
      REAL,       PARAMETER        :: GRAV     =      9.81
      INTEGER,    PARAMETER        :: IMISS    = -99999
      INTEGER,    PARAMETER        :: ISMTH    =      0
      INTEGER,    PARAMETER        :: JSMTH    =      0
      REAL,       PARAMETER        :: LAPSE    =      6.5E-03
      REAL,       PARAMETER        :: LV       =      2.5E+06
      REAL,       PARAMETER        :: M2FEET   =      3.281
      REAL,       PARAMETER        :: MPS2KNTS =      1.944
      REAL,       PARAMETER        :: MPS2MPH  =      2.237
      REAL,       PARAMETER        :: PI       =      3.1415927
      REAL,       PARAMETER        :: P0       =   100000.0
      REAL,       PARAMETER        :: R        =    287.04
      REAL,       PARAMETER        :: RV       =    461.5
      REAL,       PARAMETER        :: T0       =    273.15
      REAL,       PARAMETER        :: XMISS    = -99999.9
      REAL,       PARAMETER        :: CPOG     = CP / GRAV
      REAL,       PARAMETER        :: CPOR     = CP / R
      REAL,       PARAMETER        :: DEG2RAD  = PI / 180.0
      REAL,       PARAMETER        :: E        = R / RV
      REAL,       PARAMETER        :: GOR      = GRAV / R
      REAL,       PARAMETER        :: KAPPA    = R / CP
      REAL,       PARAMETER        :: RAD2DEG  = 180.0 / PI
      REAL,       PARAMETER        :: ROG      = R / GRAV
      REAL,       PARAMETER        :: RVOLV    = RV / LV

END MODULE constants


