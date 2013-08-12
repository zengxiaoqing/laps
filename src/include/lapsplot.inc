

!       Namelist data structure
	type namelist_type
          sequence

          logical l_discrete
          logical l_sphere
          logical l_low_fill
          logical l_high_fill

          real continent_line_width
          real country_line_width
          real state_line_width
          real county_line_width
          real time_zone
          real dist_plot_ua
          real dist_plot_sfc
          real chigh_sfcwind
          real chigh_3dwind
          real chigh_cape
          real chigh_tpw
          real power_tpw

          integer latlon_int
          integer mode_supmap
          integer icol_barbs
          integer icol_continent
          integer icol_country   
          integer icol_state     
          integer icol_county    
          integer iraster
          integer i_background_color
          integer i_pcp_sto_colorbar
          integer i_sno_sto_colorbar

          character*3 c3_time_zone
          character*30 c_institution
          character*6 c_vnt_units
          character*7 c_tpw_units
          character*7 c_units_type
          character*7 c_pbl_depth_units
          character*10 c_ob_color
          character*20 btemp_colortable

	end type namelist_type

        type (namelist_type) namelist_parms


!       Other plot parameters
	type plot_type
          sequence

          logical l_discrete
          logical l_hinterp_zoom

          real rimage_intensity
          real zoom
          real xcen
          real ycen
          real zoom_wdw
          real obs_size
          real color_power

          integer contour_line_width
          integer ncols
          integer icol_barbs
          integer iraster

	end type plot_type

        type (plot_type) plot_parms
