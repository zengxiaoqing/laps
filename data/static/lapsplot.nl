 &lapsplot_nl
 latlon_int=0,
 continent_line_width=1.0,
 c3_time_zone='UTC',
 c_institution='NOAA/ESRL',
 time_zone=0.,
 c_vnt_units='KT-FT',
 c_units_type='english',
 c_pbl_depth_units='metric',
 l_discrete=.false.,
 l_sphere=.false.,
 l_low_fill=.true.,
 l_high_fill=.true.,
 mode_supmap=3,
 iraster=0,
 icol_barbs=0,
 dist_plot_ua=0.,
 /

c LAPSPLOT PARAMETERS
c
c latlon_int - interval for plotting lat/lon lines in integer degrees
c              (0) means no latlon lines plotted
c
c continent_line_width - Width of continental boundaries is adjusted by this
c                        parameter. County line width is unaffected by this 
c                        parameter and will have a line width of 1.
c
c                        If this is set to 0 then all boundaries (including
c                        continental, state, and county) will be suppressed.
c
c c3_time_zone - initials representing time zone (up to 3 characters)
c
c c_institution - institution used in label (up to 9 characters)
c
c time_zone - local time minus UTC in hours (real number)
c
c c_vnt_units - units for use in Ventilation Index (valid values are 
c               'KT-FT' or 'default')
c
c c_units_type - default type of units ('english' or 'metric')
c
c c_pbl_depth_units - PBL depth units ('english' or 'metric')
c
c l_discrete - default type of image colortable, '.false.' means more
c              continuous colors, '.true.' means fewer more discrete colors
c
c l_sphere - apply adjustements and compensating distortions to plots so they
c            will appear better when viewed on a spherical projection device
c            such as Science On a Sphere (TM)
c
c l_low_fill - do vertical extrapolation of individual 3-D radar reflectivity
c              plots (vxx files)
c
c l_high_fill - do vertical interpolation of individual 3-D radar reflectivity
c               plots (vxx files)
c
c mode_supmap - (1) use map database in 'data/static/ncarg/*'
c               (3) use ezmap database in 'NCARG_ROOT/lib/ncarg/database'
c               (4) use rangs database linked into 'data/static/ncarg/rangs/*'
c
c iraster - (0)  use default settings for raster image plots (vs contour fill)
c           (+1) use faster raster image plots whenever possible    
c           (-1) use contour fill image plots exclusively for better quality
c
c icol_barbs - (0) use default situation dependent settings for wind barb color
c              (1) always have colored wind barbs (instead of black)
c
c dist_plot_ua - used to thin out the plotted upper air (ACARS & Cloud Drift 
c                Wind) observations and represents the threshold distance in 
c                grid points between obs (e.g. 1.0). A value of zero will 
c                disable the thinning to see all the obs.
