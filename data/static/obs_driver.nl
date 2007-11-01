 &obs_driver_nl
 path_to_metar='/public/data/metar/netcdf/',
 path_to_local_data='/public/data/ldad/mesonet/netcdf/',
 path_to_buoy_data='/public/data/maritime/netcdf/',
 path_to_gps_data='/null/public/data/gpsmet/netcdf/',
 path_to_tower_data='/data/fxa/LDAD/mesonet/met-tower/netCDF/',
 metar_format='default',
 minutes_to_wait_for_metars=10,
 ick_metar_time=0,
 itime_before=900,
 itime_after=900,
 maxobs=160000,
 i4wait_local_obs_max=0,
 local_obs_thresh=0,
 l_allow_empty_lso=.false.,
 /

c Obs Driver ingest (obs_driver.exe)
c
c 'path_to_metar' - Directory for metar/synop data
c
c 'path_to_local_data' - Directory for local/mesonet or LDAD mesonet data.
c                        This can include ASOS data or MADIS LDAD mesonet.
c                        If the data are non-MADIS the full path should 
c                        generally be given. If it's MADIS data, then give the 
c                        parent directory just above the level where the 
c                        'mesonet/netCDF' and 'urbanet/netCDF' data are being
c                        stored. For MADIS data the 'path_to_local_data'
c                        variable must include the string 'madis'.
c
c 'path_to_buoy_data' - Directory for buoy/ship (maritime) data
c
c 'path_to_gps_data' - Directory for gps/profiler surface data in the event
c                      this is stored separately from LDAD
c
c 'path_to_tower_data' - Path to Met-tower data for soil moisture (when
c                        c8_project='RSA')
c
c 'metar_format' - A value of 'default' means that we are using the 
c                  'c8_project' in 'nest7grid.parms' to specify the METAR
c                  format. An override to this can be specified as follows:
c
c                  'NIMBUS' denotes GSD NetCDF format and NIMBUS file timestamp
c                  'WFO' denotes AWIPS NetCDF format and WFO file timestamp,
c                        as well as MADIS
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
c
c 'i4wait_local_obs_max' - number of seconds (since systime) to wait for a 
c                          sufficient number of LDAD local obs
c
c 'local_obs_thresh' - threshold number of LDAD local obs to trigger wait loop
