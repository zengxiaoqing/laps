 &obs_driver_nl
 path_to_metar='/public/data/metar/netcdf/',
 path_to_local_data='/public/data/ldad/mesonet/netcdf/',
 path_to_buoy_data='/public/data/maritime/netcdf/',
 path_to_gps_data='/public/data/gpsmet/netcdf/',
 metar_format='default',
 minutes_to_wait_for_metars=10,
 ick_metar_time=0,
 itime_before=900,
 itime_after=900,
 maxobs=6000,
 l_allow_empty_lso=.false.,
 /

c Obs Driver ingest (obs_driver.exe)
c
c 'path_to_metar' - Directory for metar/synop data
c
c 'path_to_local_data' - Directory for local/mesonet or LDAD data
c
c 'path_to_buoy_data' - Directory for buoy/ship (maritime) data
c
c 'path_to_gps_data' - Directory for gps/profiler surface data in the event
c                      this is stored separately from LDAD
c
c 'metar_format' - A value of 'default' means that we are using the 
c                  'c8_project' in 'nest7grid.parms' to specify the METAR
c                  format. An override to this can be specified as follows:
c
c                  'NIMBUS' denotes FSL NetCDF format and NIMBUS file timestamp
c                  'WFO' denotes FSL NetCDF format and WFO file timestamp
c                  'CWB' is Central Weather Bureau in Taiwan
c                  'AFWA' is Air Force Weather Agency
c
c 'minutes_to_wait_for_metars' - Especially helpful on WFO
c
c 'ick_metar_time' - 0=don't check, 1=check
c.....      Do we check the METARs for time?  If yes, then off hour LSO
c.....      or LS2 files will not include top of the hour data.  This will
c.....      probably be ok if there are a lot of off hour mesonet data, but
c.....      not so good otherwise.
c
c 'itime_before'
c 'itime_after'
c.....      Time window.  Set these variables for time filtering the data.
c.....      Units are seconds.  For example, if you want data with an 
c.....      observation time from 15 min before to 5 min after the analysis 
c.....      time to be included in the LSO file, use 900 and 300 for 
c.....      time_beforeand time_after, respectively.
c
c 'maxobs' - max number of surface stations in raw data
c
c 'l_allow_empty_lso' - .true. means we will write out an empty LSO file
c                       even if there are no obs. .false. means no LSO file
c                       will be written when there are no obs.
