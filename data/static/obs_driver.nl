 &obs_driver_nl
 path_to_metar='/public/data/metar/netcdf/',
 path_to_local_data='/public/data/ldad/mesonet/netcdf/',
 path_to_buoy_data='/data/fxa/point/maritime/netcdf/',
 metar_format='FSL',
 minutes_to_wait_for_metars=10,
 /

c Obs Driver ingest (obs_driver.exe)
c
c 'path_to_metar' - Directory for metars
c
c 'path_to_local_data' - Directory for local data
c
c 'path_to_buoy data' - Directory for buoy data
c
c 'metar_format' - Format of metar data. Valid values are 'FSL' and 'CWB'
c
c 'minutes_to_wait_for_metars' - Especially helpful on WFO
