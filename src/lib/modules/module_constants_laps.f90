MODULE constants_laps

  ! Constants used by laps     

      REAL,       PARAMETER        :: CM2INCH  =   1.0 / 2.54
      REAL,       PARAMETER        :: M2INCH   =   100. / 2.54
      REAL,       PARAMETER        :: INCH2M   =   .0254
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
      REAL,       PARAMETER        :: R_GAS    =    287.04
      REAL,       PARAMETER        :: R        =    R_GAS
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
      REAL,       PARAMETER        :: rmissing = 1.e9
      INTEGER,    PARAMETER        :: imissing = NINT(rmissing)
  
END MODULE constants_laps


