      subroutine config_satellite_lvd(istatus)
c
cdoc  Reads static/satellite_lvd.nl file.

      character nest7grid*150
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'
      include 'grid_fname.cmn'

      include 'satdata_lvd.for'

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_lvd_nl,err=901)
      close(1)
      istatus = 1
      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading satellite_nl in ',nest7grid
      write(*,satellite_lvd_nl)
      return
      end
c
c =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine get_sat_sounder_info(n_sndr,
     +c_sndr_id,n_sndr_channels,path_to_sat_sounder,n_elems,
     +n_lines,channel_wavelength_u,imsng_sndr_pix,pct_req_lsr,
     +istatus)
c
cdoc  Reads static/sat_sounder.nl file.

      include 'lsr_dims.inc'
      include 'grid_fname.cmn'       !grid_fnam_common

      integer len_dir
      integer n_sndr
      integer n_sndr_channels
      integer n_elems(max_sat)
      integer n_lines(max_sat)
      integer imsng_sndr_pix
      integer istatus
      character*6   c_sndr_id(max_sat)
      character*150 nest7grid
      character*200 path_to_sat_sounder(max_sat)

      real*4 channel_wavelength_u(max_ch,max_sat)
      real*4 pct_req_lsr

      namelist /satellite_sounder_nl/ n_sndr,c_sndr_id,path_to_sat_sound
     +er,n_elems,n_lines,n_sndr_channels,channel_wavelength_u,imsng_sndr
     +_pix,pct_req_lsr

      call get_directory(grid_fnam_common,nest7grid,len_dir)

      nest7grid = nest7grid(1:len_dir)//'/sat_sounder.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_sounder_nl,err=901)
      close(1)

      istatus = 1
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading satellite_sounder_nl in ',nest7grid
      write(*,satellite_sounder_nl)
      stop
      end
c
c-----------------------------------------------------------
c
      subroutine get_balance_nl(lrunbal,gamo,delo,istatus)
c
cdoc  Reads static/balance.nl file.

      implicit none

      integer    istatus
      integer    len_dir
      logical    lrunbal
      real*4     gamo,delo
      character  nest7grid*150

      include   'grid_fname.cmn'       !grid_fnam_common

      namelist /balance_nl/lrunbal,gamo,delo

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'balance.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,balance_nl,err=901)
      close(1)
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading balance.nl in ',nest7grid
      write(*,balance_nl)
      stop
      end
c
c ---------------------------------------------------------
c
      subroutine mosaic_radar_nl(c_radar_mosaic_type,n_radars,
     & c_radar_ext,i_window,imosaic_3d,istatus)
c
cdoc  Reads static/radar_mosaic.nl file

      implicit none

      include    'radar_mosaic_dim.inc'
      include    'grid_fname.cmn'              !grid_fnam_common

      integer    istatus
      integer    len_dir
      integer    n_radars
      integer    i_window
      integer    imosaic_3d
      character  c_radar_mosaic_type*3
      character  c_radar_ext(max_radars_mosaic)*3
      character  nest7grid*150

      namelist /radar_mosaic_nl/c_radar_mosaic_type,n_radars,
     & c_radar_ext,i_window,imosaic_3d

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'radar_mosaic.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,radar_mosaic_nl,err=901)
      close(1)
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading radar_mosaic.nl in ',nest7grid
      write(*,radar_mosaic_nl)
      stop
      end
c
c -------------------------------------------------------------
c
      subroutine get_background_info(bgpaths,bgmodels
     +,oldest_forecast,max_forecast_delta,use_analysis,cmodel
     +,itime_inc)

cdoc reads static/background.nl

      implicit none
      include 'bgdata.inc'
      include 'grid_fname.cmn'             !grid_fnam_common

      character*150 nest7grid
      character*256 bgpaths(maxbgmodels)
      character*132 cmodel(maxbgmodels)
      integer bgmodels(maxbgmodels), len_dir
      integer oldest_forecast, max_forecast_delta
      integer itime_inc
      logical use_analysis
      namelist /background_nl/bgpaths,bgmodels
     +,oldest_forecast,max_forecast_delta,use_analysis,cmodel
     +,itime_inc

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif
      nest7grid = nest7grid(1:len_dir)//'background.nl'

      open(1,file=nest7grid(1:len_dir+13),status='old',err=900)
      read(1,background_nl,err=901)
      close(1)
      return
 900  print*,'error opening file ',nest7grid
      stop
 901  print*,'error reading background_nl in ',nest7grid
      write(*,background_nl)
      stop
      end
c
c --- OSSE namelist reader ---
c
      subroutine get_osse_information(path_to_model,
     1                                c_obs_types,
     1                                n_sim_obs,
     1                                istatus)
c
cdoc  Reads static/osse.nl file.

      implicit none

      include 'osse.inc'
      include 'grid_fname.cmn'                 !grid_fnam_common

      integer        istatus
      integer        len_dir
      integer        i,n_sim_obs

      character      nest7grid*150
      character      path_to_model*150
      character      c_obs_types(max_ob_types)*15

      namelist /osse_nl/path_to_model,c_obs_types

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif


      do i = 1,max_ob_types
         c_obs_types(i) = ' '
      enddo

      nest7grid = nest7grid(1:len_dir)//'osse.nl'
      print*,'nest7grid = ',nest7grid(1:len_dir+7)

      open(1,file=nest7grid,status='old',err=900)
      read(1,osse_nl,err=901)
      close(1)
      do i = 1,max_ob_types
         if(c_obs_types(i).ne.' ')then
            n_sim_obs = n_sim_obs + 1
         endif
      enddo
      print*,'n_sim_obs types from namelist= ',n_sim_obs
      istatus = 1
     
      return
       
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading osse.nl in ',nest7grid
      write(*,osse_nl)
      stop
      end
c
c ---------------------------------------------------------------
c
       subroutine get_wind_parms(l_use_raob,l_use_cdw,l_use_radial_vel
     1                          ,weight_bkg_const_wind
     1                          ,rms_thresh_wind
     1                          ,max_pr,max_pr_levels,i_3d
     1                          ,istatus)

       include 'grid_fname.cmn'                          !grid_fnam_common

       logical l_use_raob, l_use_cdw, l_use_radial_vel
       namelist /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel
     1                   ,weight_bkg_const_wind
     1                   ,rms_thresh_wind
     1                   ,max_pr,max_pr_levels,i_3d
 
       character*150 static_dir,filename
 
       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/wind.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,wind_nl,err=901)
       close(1)

       if(max_pr .le. 0)then
           write(6,*)' ERROR: invalid or uninitialized value of '
     1              ,'max_pr in wind.nl ',max_pr
           istatus = 0
           return
       endif

       if(max_pr_levels .le. 0)then
           write(6,*)' ERROR: invalid or uninitialized value of '
     1              ,'max_pr_levels in wind.nl ',max_pr_levels
           istatus = 0
           return
       endif

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading wind_nl in ',filename
       write(*,wind_nl)
       istatus = 0
       return

       end
