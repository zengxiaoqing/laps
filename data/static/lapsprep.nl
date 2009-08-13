&lapsprep_nl
  hotstart = .true.,
  balance  = .true., 
  hydrometeor_scale_factor_pcp = 5.0,
  hydrometeor_scale_factor_cld = 0.5,
  output_format = 'wps',
  snow_thresh = 1.1,
  lwc2vapor_thresh = 1.01,
  make_sfc_uv = .false.,
	
/
 
c
c  hotstart:
c    Logical flag, set to true to pull in the five hydrometeor species
c    into the output files.  
c
c  balance:
c    Logical flag, set to true to use the balanced wind, temp, height 
c    fields.  Normally set to true if hotstart is true. If this is set
c    to true then lrunbal in balance.nl should also be set to true. 
c
c  hydrometeor_scale_factor:	
c    A factor which scales the hydrometeor concentrations for a grid
c    spacing. (hydrometeor_scale = hydrometeor_scale_factor/dx)
c    Note that dx is in kilometers. This is the scale that LAPS input values
c    are multiplied by to account for sub-grid scale effects.
c    For setting the internal hydrometeor scale to 1.0, the scale factor should
c    be equal to the grid resolution in km.
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
c    Any value of snow cover fraction from the LAPS analysis (0.->1.0) 
c    exceeding this threshold will cause the snow cover flag to be set
c    in the output.  To prevent any snow cover flags from being set, set
c    this value > 1.0.
c
c  lwc2vapor_thresh:
c    Real value, controls the conversion of cloud liquid to vapor.  Set
c    to 0 to disable.  If enabled, typical values are going to be around
c    1.0 (default value is 1.1).  If set to 1.0, cloud water will be converted
c    to vapor up until the RH for that point reaches 100%.  Any remaining 
c    cloud water will be left in place.  Values greater than 1.0 allow
c    for supersaturation (e.g., 1.1 allows 110% max RH).
c
c  make_sfc_uv:
c    Logical flag. If set to true, then the surface u/v fields from lsx
c    will be replaced with winds interpolated from the 3D isobaric 
c    u/v fields.
