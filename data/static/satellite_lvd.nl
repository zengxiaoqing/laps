  &satellite_lvd_nl
 NSATS = 2,
 NTYPES = 1,1,
 NCHANNEL = 4,4,
 CSATID = 'goes11','goes12',
 CSATTYPES = 'gvr','gvr',
 CCHANNELTYPES = 'vis','4u ','wv ','11u',
'vis','4u ','wv ','11u',
 CPATH2SAT = '/public/data/sat/fsl-gs/goes11/raw/image/fsl-pacus/netcdf/',
             '/public/data/sat/fsl-gs/goes13/raw/image/fsl-conus/netcdf/',
 L_CELL_AFWA = .FALSE.,
 L_NATIONAL = .FALSE.,
 ISKIP_BILIN = 1,
 I_DELTA_SAT_T_SEC = 900,
 SOLALT_THR_VIS = 15.,
 I_MSNG_SAT_FLAG =
0,   0,   0,   0,
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
c          - if you are using GOES 13 data please specify 'goes12' until we update the software.
c
c 'CCHANNELTYPES' - values can be 'vis', '4u', 'wv', '11u' and
c                   others for 6.7 micron and 12 micron
c
c 'SOLALT_THR_VIS' - solar elevation angle should be at least this high to
c                    allow processing of visibile satellite data in the
c                    normalization (albedo or reflectance) steps
