&lapsprep_nl
  hotstart = .true.,
  balance  = .true., 
  use_sfc_bal = .false., 
  hydrometeor_scale_factor_pcp = 1.0,
  hydrometeor_scale_factor_cld = 1.0,
  output_format = 'wps',
  snow_thresh = 0.3,
  lwc2vapor_thresh = 0.,
  ice2vapor_thresh = 0.,
  make_sfc_uv = .false.,
  use_laps_skintemp = .true.,	
  use_laps_vv = .false.,	
  lapsprep_always_write = .false.,
/
 
c
c  hotstart:
c    Logical flag, set to true to pull in the five hydrometeor species
c    (from the 'lwc' file) into the output files.  
c
c  balance:
c    Logical flag, set to true to use the balanced wind, temp, height, and
c    potentially sfc fields.  Normally set to .true. if 'hotstart' is true. If 
c    this is set to .true. then 'lrunbal' in 'balance.nl' should also be set to 
c    .true.. 
c
c  use_sfc_bal:
c    Logical flag, set this to use balanced surface fields (active only when 
c    when 'balance' is set to .true.)
c
c  hydrometeor_scale_factor:	
c    A factor which scales the hydrometeor concentrations for a grid
c    spacing. (hydrometeor_scale = -hydrometeor_scale_factor/dx)
c    Note that dx is in kilometers. For negative values, the absolute value
c    is the scale that LAPS input values are multiplied by to account for 
c    sub-grid scale effects. For positive values the scale factor is applied
c    without the grid spacing (dx) term.
c    
c
c  output_format:
c    List of character strings, one specifying each output format to
c    be made per run of lapsprep.  Valid values:
c      'mm5':  Makes files suitable for ingest into regridder
c      'wrf':  Makes files for hinterp ingest (e.g. for WRF/WRFSI version 2)
c      'wps':  Suggested for WRF/WPS version 3
c      'rams':  Makes RALPH2 format
c      'cdf':  Generic netCDF format used by FSL RAMS for hot start.
c
c  snow_thresh:
c    Real value, controls the setting of the snow cover flag in the output.
c
c    WPS format:
c    Threshold LAPS snow cover fraction used to set the SNOWH variable to either
c    0.0 or 1.0 in the WPS output. To disable the writing of the field, 
c    set > 1.0
c 
c    Gribprep format:
c    Any value of snow cover fraction from the LAPS analysis (0.->1.0) 
c    exceeding this threshold will cause the snow cover flag to be set
c    in the output.  To prevent any snow cover flags from being set, set
c    this value > 1.0.
c
c  lwc2vapor_thresh:
c    Real value, controls the saturation of cloud condensate.  Set
c    to 0 to disable.  If enabled, typical values are going to be around
c    1.0.  If set to 1.0, locations with mainly cloud liquid will be given an 
c    RH of 100%. This is done with respect to liquid saturation. Values greater 
c    than 1.0 allow for supersaturation (e.g., 1.1 allows 110% RH).
c
c  ice2vapor_thresh:
c    Real value, controls the saturation of cloud condensate.  Set
c    to 0 to disable.  If enabled, typical values are going to be around
c    1.0.  If set to 1.0, locations with mainly cloud ice will be given an 
c    RH of 100%. This is done with respect to ice saturation. Values greater 
c    than 1.0 allow for supersaturation (e.g., 1.1 allows 110% RH).
c
c  make_sfc_uv:
c    Logical flag. If set to true, then the surface u/v fields from lsx
c    will be replaced with winds interpolated from the 3D isobaric 
c    u/v fields.
c
c  use_laps_skintemp:
c    Logical flag.  If set to true, then the LAPS tsk field will be output
c    as SKINTEMP (only applies to "wps" format)
c
c  use_laps_vv
c    Logical flag.  If set to true, then the LAPS omega field will be output
c    as WW (only applies to "wps" format)
c
c  lapsprep_always_write
c    Logical flag. If set to tre, then output will be written for every analysis
c    cycle, instead of just those cycles having a model initialization.