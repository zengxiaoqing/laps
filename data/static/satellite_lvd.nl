  &satellite_lvd_nl
 NSATS = 1,
 NTYPES = 1,
 NCHANNEL = 4,
 CSATID = 'goes11',
 CSATTYPES = 'gvr',
 CCHANNELTYPES = 'vis','4u ','wv ','11u',
 CPATH2SAT = '/public/data/sat/fsl-gs/goes-west/raw/image/fsl-pacus/netcdf/',
 L_CELL_AFWA = .FALSE.,
 L_NATIONAL = .FALSE.,
 ISKIP_BILIN = 1,
 I_DELTA_SAT_T_SEC = 900,
 I_MSNG_SAT_FLAG =
0,   0,   0,   0,
 N_IMAGES = 1,
/
c
c Further details are in 'src/ingest/satellite/lvd/README'
c
c 'NTYPES' - listed for each satellite
c
c 'NCHANNEL' - listed for each satellite
c
c 'CSATID' - listed for each satellite
c          - values can be 'goes11', 'goes12', 'goes16', 'coms', 'him8', etc.
c          - if you are using GOES 13 data please specify 'goes12' until we update the software.
c
c 'CSATTYPES' - values listed for each satellite can be as follows:
c               'gvr' for raw gvar
c               'rll' for gvar NetCDF with lat/lon arrays
c               'cdf' for remapped fsl-conus [lambert]
c               'wfo' for wfo SBN (also netcdf and lambert projection)
c               'gwc' for AirForce Global Weather Center
c
c 'CCHANNELTYPES' - values can be 'vis', '4u', 'wv', '11u' and
c                   others for 6.7 micron and 12 micron
c
c 'CPATH2SAT' - listed for each satellite
