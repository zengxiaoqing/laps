 &moisture_switch
 RAOB_SWITCH = 0,
 RAOB_LOOKBACK = 3600,
 GOES_SWITCH = 8,
 CLOUD_SWITCH = 1,
 SOUNDER_SWITCH = 0,
 TIROS_SWITCH = 0,
 SAT_SKIP = 3,
 GVAP_SWITCH = 0,
 SFC_MIX = 1,
 MOD_4DDA_1 = 0,
 MOD_4DDA_FACTOR = 0.02,
 T_REF = -47.0,
 
 /
c   This is a switch for the module  lq3.  If the value of the first record is 0,
c   raob data (i.e., .snd files) will be ignored even if they are present ... (but
c   only by the moisture analysis).  If the switch is 1, then raob data will be
c   used if present.  If .snd files are not present, they will not be used
c   (obviously).
c
c   The second parameter is for the raob_file lookback time.  If this is 3600
c   it will accept data that is up to one hour old.
c   
c   The NEXT parameter is for using goes sbn data.   If this is desired, it
c   should be set to 8 for goes 8 or 9 for goes 9, if you wish to turn off using
c   satellite data this should be = 0.  Using sbn data in laps requires the
c   forward model code be present.
c   
c   The NEXT parameter is used to activate (1) or deactivate (0) using cloud data
c   in saturating the air in cloudy areas.  This is especially useful when the
c   cloud analysis is not working correctly... Just turn it off.
c   
c   The NEXT parameter is the sounder switch.  (0) use imager data for moisture
c   adjustment, (1) use sounder data.  Note that if the GOES switch is off (second
c   parameter) this switch is ignored.
c   
c   
c   The NEXT parameter is used to indicate tiros data.  If non-zero it must
c   reflect the satellite number (i.e. 12).
c
c   The sat_skip parameter will not skip=1, or skip the number of gridpoints
c   each application
C
C   GVAP_SWITCH enables use of derived goes precipitable water
