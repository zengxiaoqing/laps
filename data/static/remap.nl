 &remap_nl
 n_radars_remap=7,
 path_to_vrc_nl='rdr',
 ref_min=0.0,
 min_ref_samples=4,
 min_vel_samples=4,
 dgr=1.1,
 path_to_radar_a=
'/public/data/radar/wsr88d/wideband/kama/netcdf',
'/public/data/radar/wsr88d/wideband/kcys/netcdf',
'/public/data/radar/wsr88d/wideband/kddc/netcdf',
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'/public/data/radar/wsr88d/wideband/kgld/netcdf',
'/public/data/radar/wsr88d/wideband/kama/netcdf',
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
 laps_radar_ext_a=
'v01',
'v02',
'v03',
'v04',
'v05',
'vrc',
'vrc',
 /

c
c 'n_radars_remap'   - Number of radars (and/or radar types) to loop through 
c                      and process. If this parameter is set to '-1' it        
c                      automatically switches over to the value of 
c                      'max_radars_cmn'. Using '-1' is potentially recommended 
c                      only if one type of radar data is being used.
c
c                      For each radar, up to 3 radar times will be processed. 
c                      These are chosen among the times that we have new input 
c                      data for.  
c
c 'path_to_radar_a'  - Full path to each directory containing a set of input
c                      radar tilts/volumes. 
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
c
c 'ref_min'          - Minimum threshold applied to reflectivity gates during
c                      remapping step.
c
c 'min_ref_samples'  - Minimum number of valid reflectivity gates needed
c                      within a LAPS grid volume to do the averaging
c
c 'min_vel_samples'  - Minimum number of valid velocity gates needed
c                      within a LAPS grid volume to do the averaging
c
c 'dgr'              - Maximum gap in degrees of azimuth upon which the 
c                      horizontal reflectivity averaging is performed
