 &remap_nl
 n_radars_remap=3,
 path_to_vrc_nl='rdr',
 path_to_radar_a=
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'/data/fxa/laps/lapsprd/rdr/001/raw',
'',
'',
'',
'',
'',
'',
'',
 laps_radar_ext_a=
'v01',
'vrc',
'vrc',
'',
'',
'',
'',
'',
'',
'',
 /

c
c 'n_radars_remap'   - Number of radars (and/or radar types) to loop through 
c                      and process. For each radar, one radar time will be 
c                      processed. This time represents the most recent time 
c                      that we have new input data for. The cron should 
c                      therefore run at least as often as the most frequently 
c                      updated radar. Note that there is an upper limit of 10
c                      to this parameter; this is also the number of array 
c                      elements for each subsequent parameter.
c
c 'path_to_radar_a'  - Full path to each directory containing a set of input
c                      radar tilts/volumes. Max of 10 elements. 
c
c 'path_to_vrc_nl'   - Path used if we have 'vrc' output. Can be either 'rdr' 
c                      or 'lapsprd'. If we're using the mosaicing program 
c                      (needed if we have multiple radars), we can set this 
c                      parameter to the desired value of 'rdr'. If we're not 
c                      using the mosaicing program we set the parameter to 
c                      'lapsprd' and the use of only one radar is implied. 
c                      Additional radars would be overwriting in the same 
c                      directory which is somewhat undesirable. 
c 
c 'laps_radar_ext_a' - Extension name representing the type of output that
c                      will be generated. For full volume data (e.g. wideband)
c                      set this to 'v01', or 'v02', etc. For single tilt data
c                      or if just a few low-level tilts are available
c                      (e.g. WFO narrowband data), set this to 'vrc'.
c                      Max of 10 elements.
