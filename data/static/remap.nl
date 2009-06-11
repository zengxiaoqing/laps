 &remap_nl
 n_radars_remap=15,
 max_times=5,
 path_to_vrc_nl='rdr',
 ref_min=0.0,
 min_ref_samples=-1,
 min_vel_samples=-1,
 dgr=1.1,
 abs_vel_min=2.0,
 l_line_ref_qc=.false.,
 l_hybrid_first_gate=.false.,
 l_unfold=.false.,
 l_ppi_mode=.false.,
 path_to_radar_a=
'/public/data/radar/wsr88d/wideband/kama/netcdf',
'/public/data/radar/wsr88d/wideband/kcys/netcdf',
'/public/data/radar/wsr88d/wideband/kddc/netcdf',
'/public/data/radar/wsr88d/wideband/kftg/netcdf',
'/public/data/radar/wsr88d/wideband/kgjx/netcdf',
'/public/data/radar/wsr88d/wideband/kgld/netcdf',
'/public/data/radar/wsr88d/wideband/klnx/netcdf',
'/public/data/radar/wsr88d/wideband/kmtx/netcdf',
'/public/data/radar/wsr88d/wideband/kpux/netcdf',
'/public/data/radar/wsr88d/wideband/kriw/netcdf',
'/public/data/radar/wsr88d/wideband/ksfx/netcdf',
'/public/data/radar/wsr88d/wideband/kudx/netcdf',
'/public/data/radar/wsr88d/wideband/kuex/netcdf',
 laps_radar_ext_a=
'v01',
'v02',
'v03',
'v04',
'v05',
'v06',
'v07',
'v08',
'v09',
'v10',
'v11',
'v12',
'v13',
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
c 'max_times'        - Maxiumum number of volume times to process for each radar
c                      This should be set to a large number (e.g. 999) for 
c                      archive runs and a much smaller number for real-time to
c                      help with load balancing.
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
c 'laps_radar_ext_a' - Extension names representing the type of output that
c                      will be generated for the set of radars. For each 
c                      radar having full volume data (e.g. wideband / level2) 
c                      set this to 'v01', or 'v02', etc. Otherwise for each 
c                      radar having single tilt data or if just a few low-level 
c                      tilts are available (e.g. WFO Level 3 narrowband data), 
c                      set this to 'vrc'.
c
c 'ref_min'          - Minimum threshold applied to reflectivity gates during
c                      remapping step.
c
c 'min_ref_samples'  - Minimum number of valid reflectivity gates needed
c                      within a LAPS grid volume to do the averaging
c                      -1 defaults to letting the software make the choice
c
c 'min_vel_samples'  - Minimum number of valid velocity gates needed
c                      within a LAPS grid volume to do the averaging
c                      -1 defaults to letting the software make the choice
c
c 'dgr'              - Maximum gap in degrees of azimuth upon which the 
c                      horizontal reflectivity averaging is performed
c
c 'abs_vel_min'      - Minimum absolute value of velocity needed for gate
c                      to be used
c
c 'l_line_ref_qc'    - Flag to call Dr. Deng's (CWB) QC routine for continuous
c                      line echo deleting
c                      
c 'l_hybrid_first_gate - Flag to mask out gates near the radar for low tilts
c                        to emulate what happens with the hybrid reflectivity
c                        scan
c                      
c 'l_unfold'         - Perform unfolding on the polar NetCDF data right after
c                      it is read in. The Nyquist velocity will then be set
c                      to 'r_missing_data' to prevent further unfolding of the
c                      Cartesian data.
c                      
c 'l_ppi_mode'       - True will map one radar tilt in each LAPS level (for 
c                      testing only). Set to False when running for normal 
c                      operations.

