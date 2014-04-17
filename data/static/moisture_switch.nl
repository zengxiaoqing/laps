 &moisture_switch_nl
 COVAR_SWITCH = 0,
 PRINT_SWITCH = 0,
 RAOB_SWITCH = 1,
 RADIOMETER_SWITCH = 2000,
 RAOB_LOOKBACK = 4500,
 ENDIAN = 1,
 RAOB_RADIUS = 45000.0,
 GOES_SWITCH = 0,
 CLOUD_SWITCH = 1,
 CLOUD_D = 0,
 MAX_CDELRH_NL= 1.0,
 CF_SET_NL = 0.2,
 CLOUD_WEIGHT_NL = 10.0,
 RADIO_WT_NL = 1.0,
 SOUNDER_SWITCH = 0,
 TIROS_SWITCH = 0,
 SAT_SKIP = 1,
 GVAP_SWITCH = 1,
 IHOP_FLAG = 1,
 TIME_DIFF = 3900,
 GPS_SWITCH = 1,
 SFC_MIX = 1,
 MOD_4DDA_1 = 0,
 MOD_4DDA_FACTOR = 0.02,
 T_REF = -10.0,
 PATH_TO_GVAP12 = '/public/data/sat/nesdis/goes12/tpw/sfov_ihop/ascii/',
 PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
 PATH_TO_GPS = '/public/data/gpsmet/netcdf/',
 PATH2COVAR = ' '
/
c   
c# PATH_TO_GVAP8 = '/public/data/sat/nesdis/goes8/tpw/sfov_ihop/ascii/',
c   This is a switch for the module  lq3.  The data are delivered in 
c   default form.  It is up to the user to define these for their 
c   particular applicaton.  Here are the definitions of the current 
c   values used.
c  PRINT_SWITCH = 0, disables long printouts, now default, use for debugging
c  RAOB_SWITCH = 0,    raob on/off
c  RADIOMETER_SWITCH = 2000,  meters from ground level to accept data in vert.
c  RAOB_LOOKBACK = 3600, raob latency seconds
c  ENDIAN = 1 (big endian machines [default], =0 little endian machines)
c  GOES_SWITCH = 12,  goes switch on/off assign sat
c  CLOUD_SWITCH = 1, cloud usage on/off
c  CLOUD_D = 1, cloud field presence is mandated to produce output
c  MAX_CDELRH_NL = 0.1,  Steve's factor value for the cloud saturature (cloud_sat.f) module function
c  CF_SET_NL = =.3,  Adjustable cutoff (cloud fraction) above which the cloud_sat.f functionality is invoked.
C  cloud_weight_nl = 0.5, adjustable weight for cloud analysis weight.  nominal 0.5 (higher = more weight to term in variational solution)
c  radio_wt_nl = 1.0, radiometer weight control in variational step
c  SOUNDER_SWITCH = 1,  sounder switch  on/off, [ONLY ACTIVE IF GOES_SWITCH IS NON-ZERO]
c  TIROS_SWITCH = 0,   use of tiros on/off (inactive)
c  SAT_SKIP = 2,      process every 2nd, 3rd, 10th... whatever gridpoint, if 0 or neg, will force fail the call to powell.
c  GVAP_SWITCH = 1,    gvap on/off
c  IHOP_flag = 1,     (off = 0), (normally on = 1) on for nesdis off for cimss
c  TIME_DIFF = 3900,  gvap AND gps latency seconds
c  GPS_SWITCH = 1,   (off = 0, 1 = public, default, 2 = MADIS)
c  SFC_MIX = 0,        surface mixing on/off
c  MOD_4DDA_1 = 0,     4dda airforce compounded drying factor on/off
c  MOD_4DDA_FACTOR = 0.02,  drying value (+ = drying)
c  T_REF = -132.0,     reference temp in C where all is ice
c  PATH_TO_GVAP8 = '/public/data/sat/ssec/goes8/ascii/',
c  PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
c  PATH_TO_GPS = '/public/data/gpsmet/netcdf/'  madis path different if (2) above selected
