c
c.....      LAPS Surface analysis include file.
c.....	      Last update:  12 Dec 1999  PAS
c............................................
c.....  This file contains variables that are common to many
c.....  parts of the LAPS Surface analysis and that will need 
c.....  to be changed for porting LAPS to different domains.
c.....  Most, if not all, porting changes will be done automatically
c.....  when 'lapsparms.for' is changed.
c
      parameter(badflag = -99.9)


!       Namelist data structure
	type sfc_nl_type
          sequence

          real rms_wind
          real rms_temp
          real rms_dewpoint

          real bad_tgd_land
          real bad_tgd_water

	end type sfc_nl_type

        type (sfc_nl_type) sfc_nl_parms

