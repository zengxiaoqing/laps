      subroutine get_wrfsi_config(istatus)

c
c routine returns a common set of variables that describe a domain projection.
c

        implicit none
c
c all these include files are needed to marry and namelist besides nest7grid.parms
c to the lapsparms.cmn.  At the moment we only deal with one additional namelist wrfsi.nl 
c
c       include 'wrf_dims.inc'

        include 'lapsparms.cmn'

        integer      nx_l,ny_l,nz_l
        character*8  c_vcoordinate
        real*4       PRESSURE_BOTTOM
        real*4       PRESSURE_INTERVAL
        integer      laps_cycle_time
        integer      i_perimeter
        integer      i2_missing_data
        logical*1    l_highres_laps,lpad1,lpad2,lpad3
        real*4       r_missing_data
        integer      MAX_RADARS
        real*4       ref_base
        real*4       ref_base_useable
        integer      maxstns
        integer      N_PIREP
        integer      vert_rad_meso
        integer      vert_rad_sao
        integer      vert_rad_pirep
        integer      vert_rad_prof     
        character*3  craddat_type
        character*50 c50_lowres_dir
        character*8  radarext_3d
        character*8  radarext_3d_accum
        real*4       silavwt_parm
        real*4       toptwvl_parm
        character*8  c8_project
        character*9  fdda_model_source(10)  !models_max)

        include 'wrf_horzgrid.cmn'
        include 'wrf_vgridspec.cmn'
        include 'wrf_sfcfiles.cmn'
        include 'wrf_projname.cmn'
        include 'wrf_rawdatapaths.cmn'

c       include 'wrf_laps_analysis.cmn'

 
        integer       i,nest,lvc
        integer       istatus

        real*4        grid_spacing_wrf_m
        real*4        grid_spacing_m
        character*6   wrftolaps_c6_maprojname
        character*8   c_analysis_type

        call read_wrfsi_hgridspec (istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_hgridspec'
           return
        endif

        call read_wrfsi_sfcfiles (istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_sfcfiles'
           return
        endif

        call read_wrfsi_project_id (istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_grid_fname'
           return
        endif

        call read_wrfsi_rawdatapaths (istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_rawdatapaths'
           return
        endif

c       call read_wrfsi_vgridspec (istatus)
c       if(istatus.ne.1)then
c          print*,'error reading wrfsi_vgridspec'
c          return
c       endif

        call read_wrfsi_laps_control (
     +  nx_l,ny_l,nz_l,c_vcoordinate,grid_spacing_m
     + ,pressure_bottom,pressure_interval
     + ,laps_cycle_time
     + ,l_highres_laps,lpad1,lpad2,lpad3
     + ,i_perimeter,c50_lowres_dir
     + ,craddat_type
     + ,radarext_3d,radarext_3d_accum
     + ,i2_missing_data,r_missing_data
     + ,max_radars,ref_base,ref_base_useable
     + ,maxstns,n_pirep
     + ,vert_rad_meso,vert_rad_sao,vert_rad_pirep
     + ,vert_rad_prof,silavwt_parm,toptwvl_parm
     + ,c8_project,fdda_model_source,istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_laps_control'
           return
        endif
c
c somehow we need to determine which nest we are processing
c
        nest = 1
        standard_latitude  =moad_stand_lats(1)
        standard_latitude2 =moad_stand_lats(2)
        standard_longitude =moad_stand_lons(1)
        grid_cen_lat_cmn   =moad_known_lat
        grid_cen_lon_cmn   =moad_known_lon
        path_to_topt10m=topo_10m
        path_to_topt30s=topo_30s
        path_to_pctl10m=pctland_10m
        path_to_soil2m = '/null'

        c6_maproj=wrftolaps_c6_maprojname(map_proj_name)

        c80_description = simulation_name_cmn

c hardwires for the time being since these are not in wrfsi.nl
        silavwt_parm_cmn = 0.
        toptwvl_parm_cmn = 2.

        path_to_raw_pirep_cmn = path_to_raw_pirep
        path_to_raw_rass_cmn = path_to_raw_rass
        path_to_raw_profiler_cmn = path_to_raw_profiler
        path_to_raw_blprass_cmn = path_to_raw_blprass
        path_to_raw_blpprofiler_cmn = path_to_raw_blpprofiler
        path_to_wsi_2d_radar_cmn = path_to_wsi_2d_radar
        path_to_wsi_3d_radar_cmn = path_to_wsi_3d_radar
        path_to_qc_acars_cmn = path_to_qc_acars

        c_analysis_type = 'laps'
        call s_len(c_vcoordinate,lvc)

        if(c_analysis_type .eq. 'laps')then   !this means that laps will be configured
!                                              using laps_analysis_control
           nx_l_cmn = nx_l
           ny_l_cmn = ny_l
           nk_laps = nz_l
           grid_spacing_m_cmn = grid_spacing_m
           vertical_grid = c_vcoordinate(1:lvc)

c       else                                  !this means that laps will be configured
!                                              with the wrf h/v gridspecs.

c we also need to get the level info (p_bot, p_int, etc)
c          nx_l_cmn = xdim(nest)
c          ny_l_cmn = ydim(nest)
c          nk_laps = moad_nz
c          vertical_grid = verticalcoord
c          grid_spacing_m_cmn = grid_spacing_wrf_m(nest)

        endif

        pressure_bottom_l = pressure_bottom
        pressure_interval_l = pressure_interval
        laps_cycle_time_cmn = laps_cycle_time
        l_highres = l_highres_laps
        i_perimeter_cmn = i_perimeter
        c50_lowres_directory = c50_lowres_dir
        c_raddat_type = craddat_type
        radarext_3d_cmn = radarext_3d
        radarext_3d_accum_cmn = radarext_3d_accum
        i2_missing_data_cmn = i2_missing_data
        r_missing_data_cmn = r_missing_data
        max_radars = max_radars_cmn
        ref_base_cmn = ref_base
        ref_base_useable_cmn = ref_base_useable
        maxstns_cmn = maxstns
        n_pirep_cmn = n_pirep
        vert_rad_meso_cmn = vert_rad_meso
        vert_rad_sao_cmn = vert_rad_sao
        vert_rad_pirep_cmn = vert_rad_pirep
        vert_rad_prof_cmn = vert_rad_prof
        c8_project_common = c8_project

        do i=1,10
           fdda_model_source_cmn(i) = fdda_model_source(i)
        enddo

        PRESSURE_0_L = PRESSURE_BOTTOM_L + PRESSURE_INTERVAL_L

        iflag_lapsparms_cmn = 1

        return
        end
c
c---------------------------------------------------------

        function grid_spacing_wrf_m(nest)

        include 'wrf_horzgrid.cmn'
        grid_spacing_wrf_m = moad_delta_x
        i = nest
 
        do while (i.ne.1)
         grid_spacing_wrf_m = grid_spacing_wrf_m / ratio_to_parent(i)
         i = parent_id(i)
        end do
 
        return
        end  
