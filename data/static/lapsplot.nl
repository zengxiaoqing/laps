 &lapsplot_nl
 latlon_int=0,
 continent_line_width=1.0,
 c3_time_zone='UTC',
 c_institution='NOAA/FSL',
 time_zone=0.,
 c_vnt_units='KT-FT',
 c_units_type='english',
 c_pbl_depth_units='metric',
 l_discrete=.false.,
 mode_supmap=3,
 /

c LAPSPLOT PARAMETERS
c
c latlon_int - interval for plotting lat/lon lines in integer degrees
c              (0) means no latlon lines plotted
c
c continent_line_width - width of continental boundaries (outside the 
c                        contiguous U.S.)
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
c mode_supmap - (1) use map database in 'data/static/ncarg/*'
c               (3) use ezmap database in 'NCARG_ROOT/lib/ncarg/database'
c

