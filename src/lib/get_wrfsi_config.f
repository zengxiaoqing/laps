      subroutine get_wrfsi_config(nest,istatus)

c
c routine returns a common set of variables that describe a domain projection.
c
        implicit none
c
c all these include files are needed to marry namelist besides nest7grid.parms
c to the lapsparms.cmn.  At the moment we only deal with one additional namelist wrfsi.nl 
c
c       include 'wrf_dims.inc'

        include 'lapsparms.cmn'

        integer      nx_l,ny_l,nz_l
        character*8  c_vcoordinate
        real       PRESSURE_BOTTOM
        real       PRESSURE_INTERVAL
        integer      laps_cycle_time
        integer      i_perimeter
c       integer      i2_missing_data
c       real       r_missing_data
        integer      MAX_RADARS
        real       ref_base
        real       ref_base_useable
        real       silavwt_parm,toptwvl_parm
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
        character*8  c8_project
        character*9  fdda_model_source(10) !maxbgmodels = 9 in src/include/bgdata.inc

        include 'wrf_horzgrid.cmn'
        include 'wrf_vgridspec.cmn'
        include 'wrf_sfcfiles.cmn'
        include 'wrf_projname.cmn'
        include 'wrf_rawdatapaths.cmn'

c       include 'wrf_laps_analysis.cmn'

        integer       num_staggers
        integer       i,nest,lvc,ltyp
        integer       istatus

        real        grid_spacing_wrf_m
        real        grid_spacing_m
        character*6   wrftolaps_c6_maprojname
        character*10  c_analysis_type

        integer       iflag_config_wrfsi
        data          iflag_config_wrfsi/0/
        save          iflag_config_wrfsi


        if(iflag_config_wrfsi.eq.1)then

           istatus = 1
           grid_spacing_m_cmn =grid_spacing_wrf_m(nest)
           lli_orig_cmn       =domain_origin_lli(nest)
           llj_orig_cmn       =domain_origin_llj(nest)
           uri_orig_cmn       =domain_origin_uri(nest)
           urj_orig_cmn       =domain_origin_urj(nest)
           parent_id_cmn      =parent_id(nest)
           ratio_to_parent_cmn=ratio_to_parent(nest)
           if(nest.eq.1)then
              nx_l_cmn           =xdim
              ny_l_cmn           =ydim
           else
              nx_l_cmn=(domain_origin_uri(nest)-
     +                  domain_origin_lli(nest))*
     +                  ratio_to_parent(nest)+1
              ny_l_cmn=(domain_origin_urj(nest)-
     +                  domain_origin_llj(nest))*
     +                  ratio_to_parent(nest)+1
           endif

           return

        endif

        call read_wrfsi_hgridspec (istatus)
        if(istatus.ne.1)then
           print*,'error reading wrfsi_hgridspec'
           return
        endif

	write(6,*) 'in config, past hgridspec call'

        r_missing_data=+1e37
        i2_missing_data=-99

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

c        call read_wrfsi_rawdatapaths (istatus)
c        if(istatus.ne.1)then
c           print*,'error reading wrfsi_rawdatapaths'
c           return
c        endif

c       call read_wrfsi_vgridspec (istatus)
c       if(istatus.ne.1)then
c          print*,'error reading wrfsi_vgridspec'
c          return
c       endif

c        call read_analysis_control(c_analysis_type,istatus)
c        if(istatus.ne.1)then
c           print*,'error reading wrfsi c_analysis_control'
c           return
c        endif

c        call read_wrfsi_laps_control (
c     +  nx_l,ny_l,nz_l,c_vcoordinate,grid_spacing_m
c     + ,pressure_bottom,pressure_interval
c     + ,laps_cycle_time
c     + ,l_highres_laps,lpad1,lpad2,lpad3
c     + ,i_perimeter,c50_lowres_dir
c     + ,craddat_type
c     + ,radarext_3d,radarext_3d_accum
c     + ,i2_missing_data,r_missing_data
c     + ,max_radars,ref_base,ref_base_useable
c     + ,maxstns,n_pirep
c     + ,vert_rad_meso,vert_rad_sao,vert_rad_pirep
c     + ,vert_rad_prof,silavwt_parm,toptwvl_parm
c     + ,c8_project,fdda_model_source,istatus)
c        if(istatus.ne.1)then
c           print*,'error reading wrfsi_laps_control'
c           return
c        endif
c
        if(nest.eq.0)then
           nest=1    !this would be the MOAD
        endif

        standard_latitude  =moad_stand_lats(1)
        standard_latitude2 =moad_stand_lats(2)
        standard_longitude =moad_stand_lons(1)
        grid_cen_lat_cmn   =moad_known_lat(nest)
        grid_cen_lon_cmn   =moad_known_lon(nest)
        num_domains_cmn    =num_domains
        lli_orig_cmn       =domain_origin_lli(nest)
        llj_orig_cmn       =domain_origin_llj(nest)
        uri_orig_cmn       =domain_origin_uri(nest)
        urj_orig_cmn       =domain_origin_urj(nest)
        parent_id_cmn      =parent_id(nest)
        ratio_to_parent_cmn=ratio_to_parent(nest)

        path_to_topt10m=topo_10m
        path_to_topt30s=topo_30s
        path_to_pctl10m=pctland_10m
        path_to_soiltype_top30s = soiltype_top_30s
        path_to_soiltype_bot30s = soiltype_bot_30s
        path_to_landuse30s = landuse_30s
        path_to_greenfrac = greenfrac
        path_to_soiltemp1deg = soiltemp_1deg
        path_to_albedo = albedo_ncep
        path_to_maxsnoalb = maxsnowalb
        path_to_islope = islope
        path_to_sst  = sstemp

        c6_maproj=wrftolaps_c6_maprojname(map_proj_name)
        c80_description = simulation_name_cmn

        path_to_raw_pirep_cmn = path_to_raw_pirep
        path_to_raw_rass_cmn = path_to_raw_rass
        path_to_raw_profiler_cmn = path_to_raw_profiler
        path_to_raw_blprass_cmn = path_to_raw_blprass
        path_to_raw_blpprofiler_cmn = path_to_raw_blpprofiler
        path_to_wsi_2d_radar_cmn = path_to_wsi_2d_radar
        path_to_wsi_3d_radar_cmn = path_to_wsi_3d_radar
        path_to_qc_acars_cmn = path_to_qc_acars

        c_analysis_type ='wrfsi'
        c_vcoordinate='pressure'

        call s_len(c_vcoordinate,lvc)
        call s_len(c_analysis_type,ltyp)
        vertical_grid = c_vcoordinate(1:lvc)

        num_staggers=num_staggers_wrf

        nx_l_cmn = xdim
        ny_l_cmn = ydim
        grid_spacing_m_cmn = grid_spacing_wrf_m(nest)
        silavwt_parm_cmn = silavwt_parm_wrf
        toptwvl_parm_cmn = toptwvl_parm_wrf

        i2_missing_data_cmn = i2_missing_data
        r_missing_data_cmn = r_missing_data

        iflag_lapsparms_cmn = 1
        iflag_config_wrfsi = 1

        return
        end
c
c---------------------------------------------------------

        function grid_spacing_wrf_m(nest)

        include 'wrf_horzgrid.cmn'
        grid_spacing_wrf_m = moad_delta_x
        i = nest
 
        do while (i.ne.1)
         if(i.eq.0)then
            print*,'Error:  parent_id = 0. Nest = ',nest
            print*,'Error:  Check wrfsi.nl variable parent_id.'
            print*,'***** Terminating *****'
            stop
         endif
         grid_spacing_wrf_m = grid_spacing_wrf_m / ratio_to_parent(i)
         i = parent_id(i)
        end do
 
        return
        end  

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCCC

	subroutine get_nmm_grid_spacing(dlmd,dphd,res)

        include 'wrf_horzgrid.cmn'
	
	if (moad_delta_x .gt. 2. .or. moad_delta_y .gt. 2.) then
	write(6,*) 'unexpectedly large deltax/deltay values'
	write(6,*) moad_delta_x, moad_delta_y
	write(6,*) 'these seem inappropriate for nmm grid'
	STOP
	endif

	dlmd=moad_delta_x
	dphd=moad_delta_y

	res=111200.*(dlmd**2.+dphd**2.)**0.5

	return
	end subroutine get_nmm_grid_spacing

CCCCCCCCCCCCCCCCCCCCCCCCCCCCCC
