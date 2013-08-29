 &lapsplot_nl
 latlon_int=0,
 continent_line_width=1.0,
 country_line_width=1.0,
 state_line_width=1.0,
 county_line_width=1.0,
 c3_time_zone='UTC',
 c_institution='NOAA/ESRL LAPS',
 time_zone=0.,
 c_vnt_units='KT-FT',
 c_units_type='english',
 c_pbl_depth_units='metric',
 chigh_sfcwind=50.,
 chigh_3dwind=200.,
 chigh_cape=7000.,
 chigh_tpw=7.,
 power_tpw=0.7,
 c_ob_color='default',
 i_background_color=2,
 l_discrete=.false.,
 l_sphere=.false.,
 l_low_fill=.true.,
 l_high_fill=.true.,
 i_pcp_sto_colorbar=3,
 i_sno_sto_colorbar=4,
 mode_supmap=3,
 iraster=0,
 icol_barbs=0,
 icol_continent=7,
 dist_plot_ua=0.,
 montage=0,
 /

c LAPSPLOT PARAMETERS
c
c latlon_int - interval for plotting lat/lon lines in integer degrees
c              (0) means no latlon lines plotted
c
c continent_line_width - Width of continental boundaries is adjusted by this
c                        parameter. If this is set to 0 then the boundaries
c                        will be suppressed.
c
c country_line_width - Width of country boundaries is adjusted by this
c                      parameter. If this is set to 0 then the boundaries
c                      will be suppressed.
c
c state_line_width - Width of state boundaries is adjusted by this
c                    parameter. If this is set to 0 then the boundaries
c                    will be suppressed. 
c
c county_line_width - Width of county boundaries is adjusted by this
c                     parameter. If this is set to 0 then the boundaries
c                     will be suppressed.
c
c c3_time_zone - initials representing time zone (up to 3 characters)
c
c c_institution - institution used in label (up to 30 characters, though they
c                 may be truncated for cross-sections down to 14)
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
c chigh_sfcwind - maximum of colorbar range for surface wind speed plots
c                 (recommended values are 50., 100., 200.)
c
c chigh_3dwind - maximum of colorbar range for 3-D wind plots
c                (recommended values are 100., 200.)
c
c chigh_cape - maximum of colorbar range for CAPE plots in J/kg
c              (recommended values are 7000., 3500.)
c
c chigh_tpw - maximum of colorbar range for CAPE plots in centimeters
c              (recommended values are 7., or smaller positive integer values)
c
c power_tpw - lower numbers skew the colors more toward the left side of the
c             colorbar
c
c c_ob_color - color to plot surface obs 'default', or 'white'
c
c i_background_color - 1 is white background, 2 is black background
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
c i_pcp_sto_colorbar - set range of storm total precip colorbar
c               (1) 0-10 inches
c               (3) 0-25 inches
c               (2) 0-40 inches
c
c i_sno_sto_colorbar - set range of storm total snow colorbar
c               (2) 0-40  inches
c               (4) 0-100 inches
c
c mode_supmap - (1) use map database in 'data/static/ncarg/*'
c               (3) use ezmap database in 'NCARG_ROOT/lib/ncarg/database'
c               (4) use rangs database with a default location of 
c                   'NCARG_ROOT/lib/ncarg/database/RANGS_GSHHS'
c                   see http://ncl.ucar.edu/Document/Graphics/rangs.shtml
c
c iraster - (0)  use default settings for raster image plots (vs contour fill)
c           (+1) use faster raster image plots whenever possible    
c           (-1) use contour fill image plots exclusively for better quality
c
c icol_barbs - (0) use default situation dependent settings for wind barb color
c              (1) always have colored wind barbs (instead of black)
c
c icol_continent - (7) use yellow color for continental outlines
c                  (1) use white color for continental outlines
c
c dist_plot_ua - used to thin out the plotted upper air (ACARS & Cloud Drift 
c                Wind) observations and represents the threshold distance in 
c                grid points between obs (e.g. 1.0). A value of zero will 
c                disable the thinning to see all the obs.
c
c montage - (0)  use default settings on whether to do montage or animation
c                with the "on-the-fly" page
c           (+1) preferentially chooses montages over animations
c           (-1) preferentially chooses animations over montages