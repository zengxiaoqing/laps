 &moisture_switch_nl
 RAOB_SWITCH = 0,
 RAOB_LOOKBACK = 3600,
 RAOB_RADIUS = 45000.0,
 GOES_SWITCH = 8,
 CLOUD_SWITCH = 1,
 CLOUD_D = 0,
 SOUNDER_SWITCH = 1,
 TIROS_SWITCH = 0,
 SAT_SKIP = 2,
 GVAP_SWITCH = 1,
 IHOP_flag = 0,
 TIME_DIFF = 3600, 
 GPS_SWITCH = 1,
 SFC_MIX = 1,
 MOD_4DDA_1 = 0,
 MOD_4DDA_FACTOR = 0.02,
 T_REF = -132.0,
 PATH_TO_GVAP8 = '/public/data/sat/ssec/goes8/ascii/',
 PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
 PATH_TO_GPS = '/public/data/gpsmet/netcdf/'
 

# PATH_TO_GVAP8 = '/public/data/sat/nesdis/goes8/tpw/sfov_ihop/ascii/',
/
c   
c   This is a switch for the module  lq3.  The data are delivered in 
c   default form.  It is up to the user to define these for their 
c   particular applicaton.  Here are the definitions of the current 
c   values used.
c  RAOB_SWITCH = 0,    raob on/off
c  RAOB_LOOKBACK = 3600, raob latency seconds
c  GOES_SWITCH = 8,  goes switch on/off assign sat
c  CLOUD_SWITCH = 1, cloud usage on/off
c  CLOUD_D = 1, cloud field presence is mandated to produce output
c  SOUNDER_SWITCH = 1,  sounder switch  on/off
c  TIROS_SWITCH = 0,   use of tiros on/off (inactive)
c  SAT_SKIP = 2,      process every 2nd, 3rd, 10th... whatever gridpoint
c  GVAP_SWITCH = 1,    gvap on/off
c  IHOP_flag = 0,     (normally off = 0), (on = 1) initiates IHOP special data
c  TIME_DIFF = 70800,  gvap latency seconds
c  SFC_MIX = 0,        surface mixing on/off
c  MOD_4DDA_1 = 0,     4dda airforce compounded drying factor on/off
c  MOD_4DDA_FACTOR = 0.02,  drying value (+ = drying)
c  T_REF = -132.0,     reference temp in C where all is ice
c  PATH_TO_GVAP8 = '/public/data/sat/ssec/goes8/ascii/',
c  PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
c  PATH_TO_GPS = '/public/data/gpsmet/netcdf/'
c  RAOB_SWITCH = 0,    raob on/off
c  GOES_SWITCH = 8,  goes switch on/off assign sat
c  CLOUD_SWITCH = 1, cloud usage on/off
c  SOUNDER_SWITCH = 1,  sounder switch  on/off
c  TIROS_SWITCH = 0,   use of tiros on/off (inactive)
c  GVAP_SWITCH = 1,    gvap on/off
c  SFC_MIX = 0,        surface mixing on/off
c  MOD_4DDA_1 = 0,     4dda airforce compounded drying factor on/off
c  MOD_4DDA_FACTOR = 0.02,  drying value (+ = drying)
c  PATH_TO_GVAP8 = '/public/data/sat/ssec/goes8/ascii/',
c  PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
c  PATH_TO_GPS = '/public/data/gpsmet/netcdf/'
