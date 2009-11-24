  &satellite_lvd_nl
 NSATS = 1
 NTYPES = 1
 NCHANNEL = 4
 CSATID = 'goes11'
 CSATTYPES = 'gvr'
 CCHANNELTYPES = 'vis','4u ','wv ','11u'
 CPATH2SAT = '/public/data/sat/fsl-gs/goes11/raw/image/fsl-pacus/netcdf/'
 L_CELL_AFWA = .FALSE.,
 L_NATIONAL = .FALSE.,
 ISKIP_BILIN = 1,
 I_DELTA_SAT_T_SEC = 900,
 I_MSNG_SAT_FLAG =
0,   0,   0,   0,
 N_IMAGES = 1,
/
c
c 'CSATTYPES' - values can be as follows:
c               'gvr' for raw gvar
c               'cdf' for remapped fsl-conus [lambert]
c               'wfo' for wfo SBN (also netcdf and lambert projection)
c               'gwc' for AirForce Global Weather Center
c
c 'CSATID' - values can be 'goes08', 'goes09', goes10', 'goes11', 'goes12', etc.
c
c 'CCHANNELTYPES' - values can be 'vis', '4u', 'wv', '11u' and
c                   others for 6.7 micron and 12 micron
