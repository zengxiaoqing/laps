!-------------------------------
      subroutine open_namelist (lun)
      implicit none
      integer lun
      character*80 filename

      lun = 10
      filename = 'wrfsi.nl'
      open(lun,file=filename,status='old',err=900)
      return

 900  print*,'error opening namelist file ', filename
      stop

      end

!-----------------------------
      subroutine read_wrfsi_project_id (istatus)

      implicit none

      include 'wrf_projname.cmn'
      include 'grid_fname.cmn'
      integer lenr
      integer istatus

      integer iflag_simname_cmn
      data iflag_simname_cmn/0/
      save iflag_simname_cmn

      character*80 filename
      character*80 simulation_name
      character*80 user_desc

      namelist /project_id/ simulation_name,user_desc

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_simname_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,project_id,err=901)
         close(93)
         iflag_simname_cmn=1
         simulation_name_cmn=simulation_name
      endif
      istatus = 1
      return
 
 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading', filename(1:8),' for project_id'
      write(*,project_id)
      return

      end

!-------------------------------
      subroutine read_wrfsi_filetimespec (lun,
     +            fname, ftype, order,
     +            start_year, start_month, start_day, start_hour,
     +            end_year, end_month, end_day, end_hour, interval)

      implicit none
      integer lun

      character*80 fname
      character*80 ftype
      character*80 order
      integer start_year, start_month, start_day, start_hour
      integer end_year, end_month, end_day, end_hour, interval

      namelist /filetimespec/ 
     +            fname,
     +            start_year, start_month, start_day, start_hour,
     +            end_year, end_month, end_day, end_hour, interval,
     +            ftype, order

      rewind(lun)
      read(lun,filetimespec,err=901)
      return

 901  print*,'error reading namelist file for filetimespec'
      write(*,filetimespec)
      stop
      end

!-------------------------------
      subroutine read_wrfsi_hgridspec (istatus)



      implicit none

      include 'wrf_horzgrid.cmn'
      include 'grid_fname.cmn'

      character*80 filename

      integer      istatus
      integer      iflag_horzgrid_cmn
      data         iflag_horzgrid_cmn/0/
      save         iflag_horzgrid_cmn

      integer      lenr

      namelist /hgridspec/
     +          num_domains,xdim,ydim,parent_id
     +         ,ratio_to_parent
     +         ,domain_origin_parent_x
     +         ,domain_origin_parent_y
     +         ,stagger_type
     +         ,map_proj_name,latlon_grid 
     +         ,moad_known_lat, moad_known_lon, moad_known_loc
     +         ,moad_stand_lats, moad_stand_lons 
     +         ,moad_delta_x, moad_delta_y
     +         ,silavwt_parm_wrf, toptwvl_parm_wrf


      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_horzgrid_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,hgridspec,err=901)
         close(93)
         iflag_horzgrid_cmn=1
      endif
      istatus = 1
      return

 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for hgridspec'
      write(*,hgridspec)
      return
      end

!-----------------------------
      subroutine read_wrfsi_vgridspec (istatus)

      implicit none
      include 'wrf_vgridspec.cmn'
      include 'grid_fname.cmn'

      namelist /vgridspec/
     + verticalcoord, moad_nz, vertical_increment
     +,vertical_stretch, max_vertical_inc, levels
     +,vstagger_type

      character*80 filename

      integer      istatus
      integer      iflag_vgridspec_cmn
      data         iflag_vgridspec_cmn/0/
      save         iflag_vgridspec_cmn

      integer      lenr

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_vgridspec_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,vgridspec,err=901)
         close(93)
         iflag_vgridspec_cmn=1
      endif
      istatus = 1
      return

 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for vgridspec'
      write(*,vgridspec)
      return
      end
!---------------------------------------------------------------
      subroutine read_wrfsi_sfcfiles (istatus)

      implicit none
      include 'wrf_sfcfiles.cmn'
      include 'grid_fname.cmn'

      namelist /sfcfiles/ topo_30s, topo_10m, pctland_10m
     &,landuse_30s, soiltype_top_30s, soiltype_bot_30s

      character*80 filename

      integer      istatus
      integer      iflag_sfcfiles_cmn
      data         iflag_sfcfiles_cmn/0/
      save         iflag_sfcfiles_cmn

      integer      lenr

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_sfcfiles_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,sfcfiles,err=901)
         close(93)
         iflag_sfcfiles_cmn=1
      endif
      istatus = 1
      return

 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for sfcfiles'
      write(*,sfcfiles)
      return
      end
!----------------------------
      subroutine read_wrfsi_rawdatapaths (istatus)

      implicit none
      include 'wrf_rawdatapaths.cmn'
      include 'grid_fname.cmn'

      character*80 filename

      integer      istatus
      integer      iflag_rawdatapaths_cmn
      data         iflag_rawdatapaths_cmn/0/
      save         iflag_rawdatapaths_cmn

      namelist /paths_to_raw_data/
     +          path_to_raw_pirep
     +         ,path_to_raw_rass, path_to_raw_profiler
     +         ,path_to_raw_blprass
     +         ,path_to_raw_blpprofiler,path_to_wsi_2d_radar
     +         ,path_to_wsi_3d_radar, path_to_qc_acars

      integer      lenr

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_rawdatapaths_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,paths_to_raw_data,err=901)
         close(93)
         iflag_rawdatapaths_cmn=1
      endif
      istatus = 1
      return

 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for paths_to_raw_data'
      write(*,paths_to_raw_data)
      return
      end
c
c -------------------------------------------------------------------
c
      subroutine read_analysis_control(c_analysis_type,istatus)

      implicit none
      include 'grid_fname.cmn'

      character*80 filename
      character*10 c_analysis_type

      integer      istatus
      integer      iflag_anl_control_cmn
      data         iflag_anl_control_cmn/0/
      save         iflag_anl_control_cmn

      namelist /analysis_control/
     +          c_analysis_type

      integer      lenr

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_anl_control_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,analysis_control,err=901)
         close(93)
         iflag_anl_control_cmn=1
      endif
      istatus = 1
      return
 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for paths_to_raw_data'
      write(*,analysis_control)
      return
      end
!
!----------------------------
!
      subroutine read_wrfsi_laps_control(
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

      implicit none

c     include 'wrf_dims.inc'
c     include 'wrf_laps_analysis.cmn'

      integer      nx_l,ny_l,nz_l
      character*8  c_vcoordinate
      real*4       grid_spacing_m
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

      include 'grid_fname.cmn'

      character*80 filename

      integer      istatus
      integer      iflag_wrf_laps_cmn
      data         iflag_wrf_laps_cmn/0/
      save         iflag_wrf_laps_cmn

      namelist /laps_analysis_control/
     +          nx_l,ny_l,nz_l
     +         ,c_vcoordinate
     +         ,grid_spacing_m
     +         ,pressure_bottom
     +         ,pressure_interval
     +         ,laps_cycle_time
     +         ,l_highres_laps
     +         ,i_perimeter
     +         ,c50_lowres_dir
     +         ,craddat_type
     +         ,radarext_3d
     +         ,radarext_3d_accum
     +         ,i2_missing_data
     +         ,r_missing_data
     +         ,max_radars
     +         ,ref_base
     +         ,ref_base_useable
     +         ,maxstns
     +         ,n_pirep
     +         ,vert_rad_meso
     +         ,vert_rad_sao
     +         ,vert_rad_pirep
     +         ,vert_rad_prof
     +         ,silavwt_parm
     +         ,toptwvl_parm
     +         ,c8_project
     +         ,fdda_model_source

      integer      lenr

      istatus = 0
      call s_len(generic_data_root,lenr)
      filename = generic_data_root(1:lenr)//'/static/wrfsi.nl'
      if(iflag_wrf_laps_cmn.ne.1)then
         open(93,file=filename,status='old',err=900)
         rewind(93)
         read(93,laps_analysis_control,err=901)
         close(93)
         iflag_wrf_laps_cmn=1
      endif
      istatus = 1
      return

 900  print*,'error opening namelist file ', filename
      return

 901  print*,'error reading namelist file for laps_analysis_control'
      write(*,laps_analysis_control)
      return
      end
c
c --------------------------------------------------------
      subroutine get_grid_dim_xy_wrf(nest,nx,ny,istatus)

      implicit  none

      include 'wrf_horzgrid.cmn'

      integer  istatus
      integer  nest
      integer  nx,ny

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_stand_lon_wrf'
         return
      endif

      istatus = 1

      nx=xdim(nest)
      ny=ydim(nest)

      return
      end
!----------------------------

      subroutine get_stand_lon_wrf(std_lon,istatus)

      implicit  none

      include 'wrf_horzgrid.cmn'

      integer  istatus
      real     std_lon

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_stand_lon_wrf'
         return
      endif

      istatus = 1

      std_lon=moad_stand_lons(1)

      return
      end
!-----------------------------

      subroutine get_stand_lats_wrf(std_lat,std_lat2
     +                              ,istatus)

      implicit  none

      include 'wrf_horzgrid.cmn'
      integer  istatus

      real std_lat,std_lat2

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_stand_lon_wrf'  
         return
      endif
 
      istatus = 1

c Note: this is not designed for nesting atm.
      std_lat = moad_stand_lats(1)
      std_lat2= moad_stand_lats(2)

      return
      end
!-----------------------------

      subroutine get_c6_maproj_wrf(c6_maproj
     +                              ,istatus)

      implicit  none

      include 'wrf_horzgrid.cmn'
      character c6_maproj*6
      character wrftolaps_c6_maprojname*6
      integer  istatus
      integer  nest

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_c6_maproj_wrf'
         return
      endif

      istatus = 1

c Note: this is not designed for nesting atm.
      c6_maproj=wrftolaps_c6_maprojname(map_proj_name)

      return
      end
!-----------------------------

      subroutine get_grid_center_lon_wrf(grid_center_lon
     +                              ,istatus)

      implicit  none

      include 'wrf_horzgrid.cmn'
      real     grid_center_lon
      integer  istatus
      integer  nest

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_grid_center_lon_wrf'
         return
      endif

      istatus = 1

      grid_center_lon=moad_known_lon

      return
      end

!-----------------------------

      subroutine get_grid_spacing_wrf(grid_spacing_m,istatus)
      implicit  none

      include 'wrf_horzgrid.cmn'
      real     grid_spacing_m
      integer  istatus
      integer  nest

      call read_wrfsi_hgridspec (istatus)
      if(istatus.ne.1)then
         print*,'error in get_grid_spacing_wrf'
         return
      endif

      istatus = 1

      grid_spacing_m=moad_delta_x

      return
      end
