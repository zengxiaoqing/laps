
       subroutine get_lapsplot_parms(namelist_parms,istatus)

       include 'lapsplot.inc'
 
       character*150 static_dir,filename
       character*3 c3_time_zone
       character*30 c_institution
       character*6 c_vnt_units
       character*7 c_tpw_units
       character*7 c_units_type
       character*7 c_pbl_depth_units
       character*10 c_ob_color
       character*20 btemp_colortable
       logical l_discrete,l_sphere,l_low_fill,l_high_fill       
       real time_zone

       namelist /lapsplot_nl/ latlon_int,continent_line_width
     1                       ,country_line_width
     1                       ,state_line_width
     1                       ,county_line_width
     1                       ,c3_time_zone,time_zone
     1                       ,c_institution,c_vnt_units,c_tpw_units
     1                       ,c_units_type,c_pbl_depth_units
     1                       ,chigh_sfcwind,chigh_3dwind,chigh_cape
     1                       ,chigh_tpw,power_tpw,scale_omega
     1                       ,l_discrete, l_sphere
     1                       ,l_low_fill, l_high_fill       
     1                       ,mode_supmap, iraster, icol_barbs
     1                       ,icol_continent,icol_country
     1                       ,icol_state,icol_county
     1                       ,dist_plot_ua, dist_plot_sfc
     1                       ,c_ob_color, i_background_color
     1                       ,btemp_colortable
     1                       ,i_pcp_sto_colorbar
     1                       ,i_sno_sto_colorbar
     1                       ,montage

!      Set defaults
       latlon_int = 0
       continent_line_width = 1.0
       country_line_width = 1.0
       state_line_width = 1.0
       county_line_width = 1.0
       c3_time_zone = 'UTC'
       time_zone = 0.0
       c_institution = 'NOAA/FSL LAPS'
       c_vnt_units = 'M**2/S'
       c_tpw_units = 'CM'
       mode_supmap = 3
       iraster = 0
       l_sphere = .false.
       icol_barbs = 0
       icol_continent = 7 ! yellow
       icol_country = 7   ! yellow
       icol_state = 7     ! yellow
       icol_county = 7    ! yellow
       chigh_cape = 7000.
       chigh_tpw = 7.
       power_tpw = 0.7
       scale_omega = 100. ! relative to Pa/S
       c_ob_color = 'default'
       i_background_color = 2
       btemp_colortable = 'linear'
       i_pcp_sto_colorbar = 3
       i_sno_sto_colorbar = 4 
       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/lapsplot.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,lapsplot_nl,err=901)
       close(1)

       print*,'success reading lapsplot_nl in ',filename
       write(*,lapsplot_nl)

!      Set namelist structure
       namelist_parms%latlon_int = latlon_int
       namelist_parms%continent_line_width = continent_line_width
       namelist_parms%country_line_width   = country_line_width
       namelist_parms%state_line_width     = state_line_width
       namelist_parms%county_line_width    = county_line_width
       namelist_parms%c3_time_zone = c3_time_zone
       namelist_parms%c_institution = c_institution
       namelist_parms%time_zone = time_zone
       namelist_parms%c_vnt_units = c_vnt_units
       namelist_parms%c_tpw_units = c_tpw_units
       namelist_parms%c_units_type = c_units_type
       namelist_parms%c_pbl_depth_units = c_pbl_depth_units
       namelist_parms%chigh_sfcwind = chigh_sfcwind
       namelist_parms%chigh_3dwind = chigh_3dwind
       namelist_parms%chigh_cape = chigh_cape
       namelist_parms%chigh_tpw = chigh_tpw
       namelist_parms%power_tpw = power_tpw
       namelist_parms%scale_omega = scale_omega
       namelist_parms%c_ob_color = c_ob_color
       namelist_parms%btemp_colortable = btemp_colortable
       namelist_parms%i_background_color = i_background_color
       namelist_parms%l_discrete = l_discrete
       namelist_parms%l_sphere = l_sphere
       namelist_parms%l_low_fill = l_low_fill
       namelist_parms%l_high_fill = l_high_fill
       namelist_parms%mode_supmap = mode_supmap
       namelist_parms%iraster = iraster
       namelist_parms%icol_barbs = icol_barbs
       namelist_parms%icol_continent = icol_continent
       namelist_parms%icol_country = icol_country
       namelist_parms%icol_state = icol_state
       namelist_parms%icol_county = icol_county
       namelist_parms%dist_plot_ua = dist_plot_ua
       namelist_parms%dist_plot_sfc = dist_plot_sfc
       namelist_parms%i_pcp_sto_colorbar = i_pcp_sto_colorbar
       namelist_parms%i_sno_sto_colorbar = i_sno_sto_colorbar

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading lapsplot_nl in ',filename
       write(*,lapsplot_nl)
       istatus = 0
       return

       end
