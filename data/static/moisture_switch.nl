 &moisture_switch
 raob_switch = 0,
 goes_switch = 9,
 cloud_switch = 1,
 sounder_switch = 1,
 tiros_switch = 0,
 /

c   This is a switch for the module  lq3.  If the value of the first record is 0,
c   raob data (i.e., .snd files) will be ignored even if they are present ... (but
c   only by the moisture analysis).  If the switch is 1, then raob data will be
c   used if present.  If .snd files are not present, they will not be used
c   (obviously).
c   
c   The second parameter is for using goes sbn data.   If this is desired, it
c   should be set to 8 for goes 8 or 9 for goes 9, if you wish to turn off using
c   satellite data this should be = 0.  Using sbn data in laps requires the
c   forward model code be present.
c   
c   The third parameter is used to activate (1) or deactivate (0) using cloud data
c   in saturating the air in cloudy areas.  This is especially useful when the
c   cloud analysis is not working correctly... Just turn it off.
c   
c   The forth parameter is the sounder switch.  (0) use imager data for moisture
c   adjustment, (1) use sounder data.  Note that if the GOES switch is off (second
c   parameter) this switch is ignored.
c   
c   
c   The fifth parameter is used to indicate tiros data.  If non-zero it must
c   reflect the satellite number (i.e. 12).
