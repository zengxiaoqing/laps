 &remap_nl
 n_radars_remap=2,
 path_to_radar_a=
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'',
'',
'',
'',
'',
'',
'',
'',
 c4_radarname_a=
'KFTG',
'KFTG',
'',
'',
'',
'',
'',
'',
'',
'',
 laps_radar_ext_a=
'v02',
'vrc',
'',
'',
'',
'',
'',
'',
'',
'',
 /
c
c 'n_radars_remap'   - Number of radars to loop through and process. For each
c                      radar, one radar time will be processed. This time 
c                      represents the most recent time that we have new input 
c                      data for. The cron should therefore run at least as often
c                      as the most frequently updated radar. Note that there
c                      is an upper limit of 10 to this parameter; this is
c                      also the number of array elements for each subsequent
c                      parameter.
c
c 'path_to_radar_a'  - Full path to each directory containing a set of radar
c                      tilts. 
c 
c 'c4_radarname_a'   - Name for each radar (4 characters).
c
c 'laps_radar_ext_a' - Extension name representing the type of output that
c                      will be generated. For full volume data (e.g. wideband)
c                      set this to 'v01', or 'v02', etc. For single tilt data
c                      or if just a few low-level tilts are available
c                      (e.g. WFO narrowband data), set this to 'vrc'.
