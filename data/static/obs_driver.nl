 &obs_driver_nl
 path_to_metar='/public/data/metar/netcdf/',
 path_to_local_data='/public/data/madis/LDAD/',
 path_to_buoy_data='/public/data/maritime/netcdf/',
 path_to_synop_data='/null/public/data/synop/netcdf/',
 path_to_gps_data='/null/public/data/gpsmet/netcdf/',
 path_to_tower_data='/data/fxa/LDAD/mesonet/met-tower/netCDF/',
 metar_format='default',
 minutes_to_wait_for_metars=10,
 ick_metar_time=0,
 itime_before=900,
 itime_after=900,
 madis_dirs='mesonet','urbanet','hfmetar','nepp','crn','hcn','hydro',
 maxobs=320000,
 i4wait_local_obs_max=0,
 local_obs_thresh=0,
 l_allow_empty_lso=.false.,
 l_multiple_reports=.false.,
 l_dupe_names=.true.,
 itest_madis_qc=15,
 n_cycles=1,
 nominal_latency=-1,
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
c 'path_to_synop_data' - Directory for SYNOP data
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
c                  'MADIS' denotes MADIS NetCDF format and MADIS file timestamp
c                  'WFO' denotes AWIPS NetCDF format and WFO file timestamp
c                  'CWB' is Central Weather Bureau in Taiwan
c                  'AFWA' is Air Force Weather Agency
c                  'AOML' is this NOAA agency's CSV format data
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
c 'madis_dirs' - first entry should be 'mesonet'
c
c 'maxobs' - max number of surface stations in raw data
c
c 'i4wait_local_obs_max' - number of seconds (since systime) to wait for a 
c                          sufficient number of LDAD local obs
c
c 'local_obs_thresh' - threshold number of LDAD local obs to trigger wait loop.
c                      A value of zero will disable the threshold. With MADIS
c                      data the threshold is applied just to the 'mesonet'
c                      and isn't applied to the 'urbanet'.
c
c 'l_allow_empty_lso' - .true. means we will write out an empty LSO file
c                       even if there are no obs. .false. means no LSO file
c                       will be written when there are no obs.
c
c 'l_multiple_reports' - .false. means that for local mesonet data we pick
c                         only the observation within the time window that is
c                         closest to the current 'systime' for each station
c                      - .true. means we will select all reports from a given
c                         station falling within the time window
c
c 'l_dupe_names' - .true. means subroutine 'check_for_dupes' will check and
c                  set to missing any stations that duplicate the name of
c                  a prior station in the list. This parameter is used only
c                  when 'l_multiple_reports' is set to .false.
c
c 'itest_madis_qc' -  0 means we do not test MADIS QC flags
c                     1 means level 1 test for MADIS QC flags with local data
c                       along with MADIS subjective black list
c                     2 means level 1,2,3 tests are all applied along with
c                       the MADIS subjective black list
c                    15 means level 1 test for MADIS QC flags and the spatial
c                       statistical QC check along with MADIS subjective black 
c                       list
c
c 'n_cycles' - number of time cycles to process LSO data for, looking back
c              in time from the current 'systime'
c
c 'nominal_latency' - number of seconds after 'systime' we want the 'obs_driver'
c                     to run with the most recent (current) cycle. This is
c                     disabled by using a negative number
                      
