 &moisture_switch
 RAOB_SWITCH = 0,
 RAOB_LOOKBACK = 3600,
 GOES_SWITCH = 8,
 CLOUD_SWITCH = 1,
 SOUNDER_SWITCH = 1,
 TIROS_SWITCH = 0,
 SAT_SKIP = 2,
 GVAP_SWITCH = 0,
 SFC_MIX = 0,
 MOD_4DDA_1 = 0,
 MOD_4DDA_FACTOR = 0.02,
 T_REF = -132.0,
 PATH_TO_GVAP8 = '/public/data/sat/ssec/goes8/ascii/',
 PATH_TO_GVAP10 = '/public/data/sat/ssec/goes10/ascii/',
 
 /
c   
c   This is a switch for the module  lq3.  The data are delivered in 
c   default form.  It is up to the user to define these for their 
c   particular applicaton.  Here are the definitions of the current 
c   values used.
c   
c    RAOB_SWITCH = 0,  (0=off, 1=on, controls use of snd files)
c    RAOB_LOOKBACK = 3600, (seconds) [inactive]
c    GOES_SWITCH = 8, (8,10, GOES sat number= on, 0=off)
c    CLOUD_SWITCH = 1, (0=off, 1=on, controls use of lc3)
c    SOUNDER_SWITCH = 0, (0=off [implies imager on], 1=on)
c    TIROS_SWITCH = 0, (0=off, 1=on) [currently inactive]
c    SAT_SKIP = 3, (3= every third gridpoint, 1=all)
c    GVAP_SWITCH = 0, (0=off, 1=on) [currently inactive]
c    SFC_MIX = 1, (1=on, 0=off) [controls interaction with sfc Q, lsx]
c    MOD_4DDA_1 = 0, (0=off,1=on) switch for using MOD_4DDA_FACTOR
c    MOD_4DDA_FACTOR = 0.02, (vertical compounded drying constant)
c    T_REF = -47.0, (centigrade temperature for ice phase ref. in rhl)
