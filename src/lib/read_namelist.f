      subroutine config_satellite_lvd(istatus)
c
cdoc  Reads static/satellite_lvd.nl file.

      character nest7grid*150
      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satellite_namelist_lvd.cmn'
      include 'grid_fname.cmn'

      integer num_sat
      integer num_types(maxsat)
      integer num_channels(maxtype,maxsat)
      integer i_qc_sat_flag(maxchannel*maxtype*maxsat)

      character*6   csatid(maxsat)
      character*3   csattypes(maxtype*maxsat)
      character*3   cchanneltypes(maxchannel*maxtype*maxsat)
      character*200 cpath2sat(maxtype*maxsat)
      integer   max_files
      parameter (max_files=20000)

      character*200 c_filenames(max_files)

      integer istatus

      include 'satdata_lvd.for'

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd.nl'
      open(1,file=nest7grid,status='old',err=900)
      read(1,satellite_lvd_nl,err=901)
      close(1)
      istatus = 1

      call conf_isjtkc(maxsat,maxtype,maxchannel,
     &csatid,csattypes,cchanneltypes,cpath2sat,
     &i_qc_sat_flag,istatus)

      if(.false.)then

      is=0
      kk=0
      do i=1,maxsat
       if(isats(i).eq.1)is=is+1
      enddo
      print*
      print*,'Number of satellites to process: ',is
      do i=1,maxsat
       if(isats(i).eq.1)then
        print*
        print*,'Sat ',i,' ID: ',c_sat_id(i)
        print*,'================================='
       endif
       do j=1,maxtype
        if(itypes(j,i).eq.1)then
         print*,'  Type: ',j, ' = ',c_sat_types(j,i)
         print*,'  ***************************************'
         print*,'  Attributes of this satellite and type: '
         print*,'  ---------------------------------------'
         print*,'   Latin/Lap/Lov: ',r_latin(j,i),r_lap(j,i),r_lov(j,i)
         print*,'   vis/ir(x/y) res (km): ', r_resolution_x_vis(j,i),
     &r_resolution_y_vis(j,i),r_resolution_x_ir(j,i),r_resolution_y_ir
     &(j,i)
         print*,'  ---------------------------------------'
         do k=1,maxchannel
          if(ichannels(k,j,i).eq.1)then
           print*,'   Chn ',k,' = ',c_channel_types(k,j,i)(1:3)
           print*,'   Path      = ',TRIM(path_to_raw_sat(k,j,i))
           if(kk.lt.1)then
            call get_file_names(path_to_raw_sat(k,j,i)
     1,numoffiles,c_filenames, max_files,istatus)
c           print*,'filenames: ',c_filenames(1:numoffiles)
            if(numoffiles.le.0)then
             print*
             print*,'!!!! Error: No data in given path'
             print*,'Path: ',TRIM(path_to_raw_sat(k,j,i))
             print*
             return
            endif
            kk=10
           endif
          endif
         enddo
        endif
       enddo
      enddo

      endif

      return

 900  print*,'error opening file ',nest7grid
      return
 901  print*,'error reading satellite_nl in ',nest7grid
      write(*,satellite_lvd_nl)

c----- don't worry about the rep namelist atm.

      return

      print*,'**************************************************'
      print*,'Now using repository version: satellite_lvd_rep.nl'
      print*,'**************************************************'

      nest7grid = nest7grid(1:len_dir)//'/satellite_lvd_rep.nl'
      open(1,file=nest7grid,status='old',err=902)
      read(1,satellite_lvd_nl,err=903)
      close(1)
      istatus = 1
      return

 902  print*,'error opening file ',nest7grid
      return
 903  print*,'error reading satellite_lvd_rep.nl in ',nest7grid
      write(*,satellite_lvd_nl)
      return
      end
c
c =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
c
      subroutine get_sat_sounder_info(n_sndr,
     +c_sndr_id,n_sndr_channels,path_to_sat_sounder,
     +channel_wavelength_u,imsng_sndr_pix,pct_req_lsr,
     +istatus)
c
cdoc  Reads static/sat_sounder.nl file.

      include 'lsr_dims.inc'
      include 'grid_fname.cmn'       !grid_fnam_common

      integer len_dir
      integer n_sndr
      integer n_sndr_channels
      integer imsng_sndr_pix
      integer istatus
      character*6   c_sndr_id(max_sat)
      character*150 nest7grid
      character*200 path_to_sat_sounder(max_sat)

      real*4 channel_wavelength_u(max_ch,max_sat)
      real*4 pct_req_lsr

      namelist /satellite_sounder_nl/ n_sndr,c_sndr_id,path_to_sat_sound
     +er,n_sndr_channels,channel_wavelength_u,imsng_sndr
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
      subroutine get_balance_nl(lrunbal,adv_anal_by_t_min
     .           ,cpads_type,incl_clom,setdelo0,istatus)
c
cdoc  Reads static/balance.nl file.

      implicit none

      integer    istatus
      integer    len_dir
      integer    adv_anal_by_t_min
      logical    lrunbal
      logical    incl_clom 
      logical    setdelo0
      character  nest7grid*150
      character  cpads_type*3

      include   'grid_fname.cmn'       !grid_fnam_common

      namelist /balance_nl/lrunbal,adv_anal_by_t_min,cpads_type
     1        ,incl_clom,setdelo0

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'balance.nl'

      incl_clom=.true.

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
c-----------------------------------------------------------
c
      subroutine get_sfcqc_nl(lrunqc,istatus)
c
cdoc  Reads static/balance.nl file.

      implicit none

      integer    istatus
      integer    len_dir
      logical    lrunqc
      character  nest7grid*150

      include   'grid_fname.cmn'       !grid_fnam_common

      namelist /sfc_qc_nl/lrunqc

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif

      nest7grid = nest7grid(1:len_dir)//'sfc_qc.nl'

      open(1,file=nest7grid,status='old',err=900)
      read(1,sfc_qc_nl,err=901)
      close(1)
      istatus = 1
      return
900   print*,'error opening file ',nest7grid
      stop
901   print*,'error reading sfc_qc.nl in ',nest7grid
      write(*,sfc_qc_nl)
      stop
      end
c
c ---------------------------------------------------------
c
      subroutine mosaic_radar_nl(c_radar_mosaic_type,n_radars,
     & c_radar_ext,i_window,mosaic_cycle_time,imosaic_3d,
     & n_radars_wideband,n_radars_narrowband,istatus)
c
cdoc  Reads static/radar_mosaic.nl file

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
     & c_radar_ext,i_window,mosaic_cycle_time,imosaic_3d,
     & n_radars_wideband,n_radars_narrowband

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
     +,forecast_length
     +,use_analysis,cmodel,itime_inc,smooth_fields,luse_sfc_bkgd
     +,lgb_only)

cdoc reads static/background.nl

      implicit none
      include 'bgdata.inc'
      include 'grid_fname.cmn'             !grid_fnam_common

      character*150 nest7grid
      character*256 bgpaths(maxbgmodels)
      character*132 cmodel(maxbgmodels)
      integer bgmodels(maxbgmodels), len_dir
      integer forecast_length
      integer itime_inc
      logical luse_sfc_bkgd
      logical use_analysis
      logical smooth_fields
      logical lgb_only
      namelist /background_nl/bgpaths,bgmodels
     +,forecast_length
     +,use_analysis,cmodel,itime_inc,smooth_fields,luse_sfc_bkgd
     +,lgb_only

      smooth_fields = .false.
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
     1                                cmodel,
     1                                c_obs_types,
     1                                n_sim_obs,
     1                                c_simob_fnames,
     1                                a9_time_init,
     1                                a4_time_fcst,
     1                                ifcst_intrvl,
     1                                isim_time_hr,
     1                                istatus)
c
cdoc  Reads static/osse.nl file.

      implicit none

      include 'grid_fname.cmn'                 !grid_fnam_common
      include 'osse.inc'

      integer        istatus
      integer        len_dir
      integer        i,n_sim_obs
      integer        ifcst_intrvl
      integer        isim_time_hr
      integer        n_simfiles

      integer        nsimobs
      data           nsimobs/0/
      save           nsimobs

      character      nest7grid*150
      character      path_to_model*150
      character      c_obs_types(maxobtype)*15
      character      c_simob_fnames(nsimfiles)*50
      character      a9_time_init*9
      character      a4_time_fcst*4
      character      cmodel*10

      namelist /osse_nl/path_to_model,cmodel,c_obs_types
     1,c_simob_fnames,a9_time_init,a4_time_fcst,ifcst_intrvl
     1,isim_time_hr

      istatus = 0

      call get_directory(grid_fnam_common,nest7grid,len_dir)
      if(nest7grid(len_dir:len_dir).ne.'/') then
        len_dir=len_dir+1
        nest7grid(len_dir:len_dir)='/'
      endif


      do i = 1,maxobtype
         c_obs_types(i) = ' '
      enddo
      do i = 1,nsimfiles
         c_simob_fnames(i) = ' '
      enddo


      nest7grid = nest7grid(1:len_dir)//'osse.nl'
      print*,'nest7grid = ',nest7grid(1:len_dir+7)

      open(1,file=nest7grid,status='old',err=900)
      read(1,osse_nl,err=901)
      close(1)

      if(nsimobs.eq.0)then
         do i = 1,maxobtype
            if(c_obs_types(i).ne.' ')then
               n_sim_obs = n_sim_obs + 1
            endif
         enddo
         nsimobs=n_sim_obs
      endif
      n_sim_obs=nsimobs
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
     1                          ,thresh_2_radarobs_lvl_unfltrd
     1                          ,thresh_4_radarobs_lvl_unfltrd
     1                          ,thresh_9_radarobs_lvl_unfltrd
     1                          ,weight_bkg_const_wind
     1                          ,weight_radar
     1                          ,rms_thresh_wind
     1                          ,max_pr,max_pr_levels,i_3d
     1                          ,istatus)

       include 'grid_fname.cmn'                          !grid_fnam_common

       logical l_use_raob, l_use_cdw, l_use_radial_vel
       integer*4 thresh_2_radarobs_lvl_unfltrd
     1          ,thresh_4_radarobs_lvl_unfltrd
     1          ,thresh_9_radarobs_lvl_unfltrd

       namelist /wind_nl/ l_use_raob, l_use_cdw, l_use_radial_vel
     1                   ,thresh_2_radarobs_lvl_unfltrd
     1                   ,thresh_4_radarobs_lvl_unfltrd
     1                   ,thresh_9_radarobs_lvl_unfltrd
     1                   ,weight_bkg_const_wind
     1                   ,weight_radar
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
c
c----------------------------------------
c
       subroutine get_gridnl(mode)
       implicit none
       integer mode
       integer len
       integer istatus
       character*256 directory
       character*256 fname
       character*200 cdataroot
       character*10  c10_grid_fname
       namelist /grid_nl/ mode

       mode = 0
       call find_domain_name(cdataroot,c10_grid_fname,istatus)
       if(istatus.ne.1)then
          print*,'Error returned from find_domain_name'
          return
       endif
       call get_directory(c10_grid_fname,directory,len)
       fname = directory(1:len)//'grid.nl'
       open(3,file=fname,status='old',err=101)
       read(3,grid_nl,err=101,end=101)
       close(3)

101    return
       end
c
c --------------------------------------------------------------
c
       subroutine read_verif_nl(type_obs,path_to_raw_profiler,
     1path_to_raw_sounding, raob_process_lag,  raob_process_lag_bal,
     1                          max_verif, verif_output_dir,
     1                          verif_missing_data, n_verif, istatus)

       implicit none

       include 'grid_fname.cmn'                          !grid_fnam_common

       integer          max_verif
       character*1	type_obs
       character*150	path_to_raw_profiler
       character*150	path_to_raw_sounding
       integer*4	raob_process_lag_Bal
       integer*4	raob_process_lag
       integer          n_verif
       integer          i,len
       character*150	verif_output_dir(4)
       real*4		verif_missing_data
       integer		istatus
       
       namelist /verif_nl/ type_obs,path_to_raw_profiler,
     1path_to_raw_sounding, raob_process_lag, raob_process_lag_bal,
     1                  verif_output_dir,
     1                  verif_missing_data

       character*150    static_dir,filename
       integer		len_dir

       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/verif.nl'

       open(1,file=filename,status='old',err=900)
       read(1,verif_nl,err=901)
       close(1)

       n_verif=0
       do i=1,max_verif
          call s_len(verif_output_dir(i),len)
          if(len.gt.0)then
             n_verif=n_verif+1
          endif
       enddo
       if(n_verif.le.0)then
          print*,'Error! Check namelist verif.nl variable'
          print*,'verif_output_dir. It appears to be empty'
          return
       endif
c      print*,'Number of verification types from verif.nl: ',n_verif
c      call s_len(verif_output_dir(1),len)
c      do i=1,n_verif
c         print*,i,' ',verif_output_dir(i)(1:len)
c      enddo

       istatus = 1
       return

  900  print*,'error opening file ',filename
       return

  901  print*,'error reading verif_nl in ',filename
       write(*,verif_nl)
       return

       end
c
c --------------------------------------------------------------
c
       subroutine read_sfc_nl(use_lso_qc,skip_internal_qc 
     1                       ,itheta, redp_lvl, del, gam, ak
     1                       ,l_require_lso
     1                       ,bad_t,bad_td,bad_u,bad_v,bad_p
     1                       ,bad_mp,bad_th,bad_the
     1                       ,bad_vis,bad_tb8
     1                       ,thresh_t,thresh_td,thresh_mslp
     1                       ,sfc_nl_parms,istatus)

       implicit none

       real    badflag
       include 'laps_sfc.inc'

       include 'grid_fname.cmn'                          !grid_fnam_common
       integer use_lso_qc, skip_internal_qc, itheta
       logical l_require_lso
       real    redp_lvl,del,gam,ak
       real    bad_t,bad_td,bad_u,bad_v,bad_p
       real    bad_mp,bad_th,bad_the
       real    bad_vis,bad_tb8
       real    bad_tgd_land,bad_tgd_water
       real    thresh_t,thresh_td,thresh_mslp
       real    rms_wind, rms_temp, rms_dewpoint

       integer istatus
       
       namelist /surface_analysis/  use_lso_qc,skip_internal_qc,
     1                              itheta, redp_lvl, del, gam, ak,       
     1                              l_require_lso,
     1                              bad_t,bad_td,bad_u,bad_v,bad_p,
     1                              bad_mp,bad_th,bad_the,
     1                              bad_tgd_land,bad_tgd_water,
     1                              bad_vis,bad_tb8,
     1                              thresh_t,thresh_td,thresh_mslp,
     1                              rms_wind, rms_temp, rms_dewpoint

       character*150    static_dir,filename
       integer len_dir

       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/surface_analysis.nl'

       open(1,file=filename,status='old',err=900)
       read(1,surface_analysis,err=901)
       close(1)

       sfc_nl_parms%rms_wind = rms_wind
       sfc_nl_parms%rms_temp = rms_temp
       sfc_nl_parms%rms_dewpoint = rms_dewpoint
       sfc_nl_parms%bad_tgd_land  = bad_tgd_land
       sfc_nl_parms%bad_tgd_water = bad_tgd_water


       istatus = 1
       return

  900  print*,'error opening file ',filename
       return

  901  print*,'error reading verif_nl in ',filename
       write(*,surface_analysis)
       return

       end
c
c --------------------------------------------------------------
c
       subroutine get_laps_redp(redp_lvl,istatus)

       implicit none
       real    badflag
       include 'laps_sfc.inc'
       include 'grid_fname.cmn'                          !grid_fnam_common
       integer use_lso_qc, skip_internal_qc, itheta
       logical l_require_lso
       real    redp_lvl,del,gam,ak
       real    bad_t,bad_td,bad_u,bad_v,bad_p
       real    bad_mp,bad_th,bad_the
       real    bad_vis,bad_tb8
       real    thresh_t,thresh_td,thresh_mslp
       real    rms_wind, rms_temp, rms_dewpoint
       integer istatus
       namelist /surface_analysis/  use_lso_qc,skip_internal_qc,
     1                              itheta, redp_lvl, del, gam, ak,
     1                              l_require_lso,
     1                              bad_t,bad_td,bad_u,bad_v,bad_p,
     1                              bad_mp,bad_th,bad_the,
     1                              bad_vis,bad_tb8,
     1                              thresh_t,thresh_td,thresh_mslp,
     1                              rms_wind, rms_temp, rms_dewpoint
       character*150    static_dir,filename
       integer len_dir
       istatus = 0
       call get_directory(grid_fnam_common,static_dir,len_dir)
       filename = static_dir(1:len_dir)//'/surface_analysis.nl'
       open(1,file=filename,status='old',err=903)
       read(1,surface_analysis,err=904)
       close(1)
       sfc_nl_parms%rms_wind = rms_wind
       sfc_nl_parms%rms_temp = rms_temp
       sfc_nl_parms%rms_dewpoint = rms_dewpoint
       istatus = 1
       return
  903  print*,'error opening file ',filename
       return
  904  print*,'error reading sfc_nl in ',filename
       write(*,surface_analysis)
       return

       end
c ------------------------------------------------------------------
c New routine that configures the satellite information needed in
c software using a limited satellite lvd namelist file
c J. Smart 5/07
c
      subroutine conf_isjtkc(nsdim,ntdim,ncdim,
     &csatid,csattypes,cchanneltypes,cpath2sat,
     &i_qc_sat_flag,istatus)

      include 'satellite_dims_lvd.inc'
      include 'satellite_common_lvd.inc'
      include 'satdata_lvd_table.for'

      integer       i_qc_sat_flag(ncdim*ntdim*nsdim)
      character*6   csatid(nsdim)
      character*3   csattypes(ntdim*nsdim)
      character*3   cchanneltypes(ncdim*ntdim*nsdim)
      character*200 cpath2sat(ntdim*nsdim)

      ISATS=0
      do i=1,nsats
       if(csatid(i).eq.'goes08')then
         ISATS(1)=1
         c_sat_id(1)=csatid(i)
       elseif(csatid(i).eq.'meteos')then
         ISATS(2)=1
         c_sat_id(2)=csatid(i)
       elseif(csatid(i).eq.'goes10')then
         ISATS(3)=1
         c_sat_id(3)=csatid(i)
       elseif(csatid(i).eq.'gmssat')then
         ISATS(4)=1
         c_sat_id(4)=csatid(i)
       elseif(csatid(i).eq.'goes12')then
         ISATS(5)=1
         c_sat_id(5)=csatid(i)
       elseif(csatid(i).eq.'goes09')then
         ISATS(6)=1
         c_sat_id(6)=csatid(i)
       elseif(csatid(i).eq.'goes11')then
         c_sat_id(7)=csatid(i)
         ISATS(7)=1
       elseif(csatid(i).eq.'noaapo')then
         c_sat_id(8)=csatid(i)
         ISATS(8)=1
       endif
      enddo

      ITYPES=0
      jj=0
      do i=1,nsats

c first satellite (goes08)
       if(csatid(i).eq.'goes08')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'gvr')then
          ITYPES(1,1)=1
          c_sat_types(1,1)=csattypes(jj)
         elseif(csattypes(jj).eq.'wfo')then
          ITYPES(2,1)=1
          c_sat_types(2,1)=csattypes(jj)
         elseif(csattypes(jj).eq.'cdf')then
          ITYPES(3,1)=1
          c_sat_types(3,1)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,1)=1
          c_sat_types(4,1)=csattypes(jj)
         endif
        enddo

c second satellite (meteosat)
       elseif(csatid(i).eq.'meteos')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'rll')then
          ITYPES(4,2)=1
          c_sat_types(4,2)=csattypes(jj)
         elseif(csattypes(jj).eq.'cdf')then
          ITYPES(3,2)=1
          c_sat_types(3,2)=csattypes(jj)
         endif
        enddo

c third satellite (goes10)
       elseif(csatid(i).eq.'goes10')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'gvr')then
          ITYPES(1,3)=1
          c_sat_types(1,3)=csattypes(jj)
         elseif(csattypes(jj).eq.'wfo')then
          ITYPES(2,3)=1
          c_sat_types(2,3)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,3)=1
          c_sat_types(4,3)=csattypes(jj)
         endif
        enddo

c forth satellite (gms)
       elseif(csatid(i).eq.'gmssat')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'hko')then
          ITYPES(2,4)=1
          c_sat_types(2,4)=csattypes(jj)
         elseif(csattypes(jj).eq.'twn')then
          ITYPES(3,4)=1
          c_sat_types(3,4)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,4)=1
          c_sat_types(4,4)=csattypes(jj)
         endif
        enddo

c fifth satellite (goes12)
       elseif(csatid(i).eq.'goes12')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'gvr')then
          ITYPES(1,5)=1
          c_sat_types(1,5)=csattypes(jj)
         elseif(csattypes(jj).eq.'wfo')then
          ITYPES(2,5)=1
          c_sat_types(2,5)=csattypes(jj)
         elseif(csattypes(jj).eq.'cdf')then
          ITYPES(3,5)=1
          c_sat_types(3,5)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,5)=1
          c_sat_types(4,5)=csattypes(jj)
         endif
        enddo

c sixth satllite (goes09)
       elseif(csatid(i).eq.'goes09')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'gvr')then
          ITYPES(1,6)=1
          c_sat_types(1,6)=csattypes(jj)
         elseif(csattypes(jj).eq.'cdf')then
          ITYPES(3,6)=1
          c_sat_types(3,6)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,6)=1
          c_sat_types(4,6)=csattypes(jj)
         endif
        enddo
c seventh satellite (goes11)
       elseif(csatid(i).eq.'goes11')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'gvr')then
          ITYPES(1,7)=1
          c_sat_types(1,7)=csattypes(jj)
         elseif(csattypes(jj).eq.'wfo')then
          ITYPES(2,7)=1
          c_sat_types(2,7)=csattypes(jj)
         elseif(csattypes(jj).eq.'cdf')then
          ITYPES(3,7)=1
          c_sat_types(3,7)=csattypes(jj)
         elseif(csattypes(jj).eq.'rll')then
          ITYPES(4,7)=1
          c_sat_types(4,7)=csattypes(jj)
         endif
        enddo
c eigth satellite (noaa polar orbiter)
       elseif(csatid(i).eq.'noaapo')then
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq.'rll')then      !raw lat-lon > could be (likely is) netcdf format
          ITYPES(4,8)=1
          c_sat_types(4,8)=csattypes(jj)
         elseif(csattypes(jj).eq.'ncp')then  !earth-projected polar stereographic netcdf format (FMI)
          ITYPES(1,8)=1
          c_sat_types(1,8)=csattypes(jj)
         endif
        enddo
       endif
      enddo

c ----
c goes08 (first satellite type)

      ICHANNELS = 0
      jj=0
      kk=0
      do i=1,nsats
       if(csatid(i).eq.'goes08')then
        do j=1,ntypes(i)
         jj=jj+1
c format type 1 (gvr)
         if(csattypes(jj).eq.'gvr')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,1)=1
            c_channel_types(1,1,1)=cchanneltypes(kk)
            i_msng_sat_flag(1,1,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,1,1)=1
            c_channel_types(2,1,1)=cchanneltypes(kk)
            i_msng_sat_flag(2,1,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,1,1)=1
            c_channel_types(3,1,1)=cchanneltypes(kk)
            i_msng_sat_flag(3,1,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,1,1)=1
            c_channel_types(4,1,1)=cchanneltypes(kk)
            i_msng_sat_flag(4,1,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,1,1)=1
            c_channel_types(5,1,1)=cchanneltypes(kk)
            i_msng_sat_flag(5,1,1)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,1,1)=cpath2sat(jj)

c format type 2 (wfo)
         elseif(csattypes(jj).eq.'wfo')then
          call s_len(cpath2sat(jj),n)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,2,1)=1
            c_channel_types(1,2,1)=cchanneltypes(kk)
            path_to_raw_sat(1,2,1)=cpath2sat(jj)(1:n)//'vis/regclip/'
            i_msng_sat_flag(1,2,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,2,1)=1
            c_channel_types(2,2,1)=cchanneltypes(kk)
            path_to_raw_sat(2,2,1)=cpath2sat(jj)(1:n)//'i39/regclip/'
            i_msng_sat_flag(2,2,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'iwv')then
            ICHANNELS(3,2,1)=1
            c_channel_types(3,2,1)=cchanneltypes(kk)
            path_to_raw_sat(3,2,1)=cpath2sat(jj)(1:n)//'iwv/regclip/'
            i_msng_sat_flag(3,2,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,2,1)=1
            c_channel_types(4,2,1)=cchanneltypes(kk)
            path_to_raw_sat(4,2,1)=cpath2sat(jj)(1:n)//'i11/regclip/'
            i_msng_sat_flag(4,2,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i12')then
            ICHANNELS(5,2,1)=1
            c_channel_types(5,2,1)=cchanneltypes(kk)
            path_to_raw_sat(5,2,1)=cpath2sat(jj)(1:n)//'i12/regclip/'
            i_msng_sat_flag(5,2,1)=i_qc_sat_flag(kk)
           endif
          enddo

c format type 3 (netcdf)
         elseif(csattypes(jj).eq.'cdf')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,1)=1
            c_channel_types(1,3,1)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,3,1)=1
            c_channel_types(2,3,1)=cchanneltypes(kk)
            i_msng_sat_flag(2,3,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'iwv')then
            ICHANNELS(3,3,1)=1
            c_channel_types(3,3,1)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,3,1)=1
            c_channel_types(4,3,1)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i12')then
            ICHANNELS(5,3,1)=1
            c_channel_types(5,3,1)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,1)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,3,1)=cpath2sat(jj)

c format type 4 (raw lat lon)
         elseif(csattypes(jj).eq.'rll')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,1)=1
            c_channel_types(1,4,1)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,1)=1
            c_channel_types(2,4,1)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,1)=1
            c_channel_types(3,4,1)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,1)=1
            c_channel_types(4,4,1)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,1)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,1)=1
            c_channel_types(5,4,1)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,1)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,1)=cpath2sat(jj)
         endif
        enddo !type for this sat -> end of types for goes08

       elseif(csatid(i).eq.'meteos')then
c
c meteos    #goes09 --- 5-12-99 J. Smart changed to AFWA METEOSAT
         
        do j=1,ntypes(i)
         jj=jj+1
c format type 3 (raw lat lon)
         if(csattypes(jj).eq.'cdf')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,2)=1
            c_channel_types(1,3,2)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,3,2)=1
            c_channel_types(2,3,2)=cchanneltypes(kk)
            i_msng_sat_flag(2,3,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,3,2)=1
            c_channel_types(3,3,2)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,3,2)=1
            c_channel_types(4,3,2)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,3,2)=1
            c_channel_types(5,3,2)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,2)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,2)=cpath2sat(jj)
c format type 4 (raw lat lon)
         elseif(csattypes(jj).eq.'rll')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,2)=1
            c_channel_types(1,4,2)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,2)=1
            c_channel_types(2,4,2)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,2)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,2)=1
            c_channel_types(3,4,2)=cchanneltypes(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,2)=1
            c_channel_types(4,4,2)=cchanneltypes(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,2)=1
            c_channel_types(5,4,2)=cchanneltypes(kk)
           endif
          enddo
         endif
        enddo
c end of types for meteos
c
c goes10 (third satellite type)
       elseif(csatid(i).eq.'goes10')then
c format type 1 (gvr)
        do j=1,ntypes(i)
         jj=jj+1
         if(csattypes(jj).eq. 'gvr')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,3)=1
            c_channel_types(1,1,3)=cchanneltypes(kk)
            i_msng_sat_flag(1,1,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,1,3)=1
            c_channel_types(2,1,3)=cchanneltypes(kk)
            i_msng_sat_flag(2,1,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,1,3)=1
            c_channel_types(3,1,3)=cchanneltypes(kk)
            i_msng_sat_flag(3,1,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,1,3)=1
            c_channel_types(4,1,3)=cchanneltypes(kk)
            i_msng_sat_flag(4,1,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,1,3)=1
            c_channel_types(5,1,3)=cchanneltypes(kk)
            i_msng_sat_flag(5,1,3)=i_qc_sat_flag(kk)
           endif
          enddo
c
c     format type 2 (wfo)
         elseif(csattypes(jj).eq. 'wfo')then
          call s_len(cpath2sat,n)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,2,3)=1
            c_channel_types(1,2,3)=cchanneltypes(kk)
            path_to_raw_sat(1,2,3)=cpath2sat(jj)(1:n)//'vis/regclip/'
            i_msng_sat_flag(1,2,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,2,3)=1
            c_channel_types(2,2,3)=cchanneltypes(kk)
            path_to_raw_sat(2,2,3)=cpath2sat(jj)(1:n)//'i39/regclip/'
            i_msng_sat_flag(2,2,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'iwv')then
            ICHANNELS(3,2,3)=1
            c_channel_types(3,2,3)=cchanneltypes(kk)
            path_to_raw_sat(3,2,3)=cpath2sat(jj)(1:n)//'iwv/regclip/'
            i_msng_sat_flag(3,2,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,2,3)=1
            c_channel_types(4,2,3)=cchanneltypes(kk)
            path_to_raw_sat(4,2,3)=cpath2sat(jj)(1:n)//'i11/regclip/'
            i_msng_sat_flag(4,2,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i12')then
            ICHANNELS(5,2,3)=1
            c_channel_types(5,2,3)=cchanneltypes(kk)
            path_to_raw_sat(5,2,3)=cpath2sat(jj)(1:n)//'i12/regclip/'
            i_msng_sat_flag(5,2,3)=i_qc_sat_flag(kk)
           endif
          enddo

c     format type 4 (raw lat lon)
         elseif(csattypes(jj).eq. 'rll')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,3)=1
            c_channel_types(1,4,3)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,3)=1
            c_channel_types(2,4,3)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,3)=1
            c_channel_types(3,4,3)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,3)=1
            c_channel_types(4,4,3)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,3)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,3)=1
            c_channel_types(5,4,3)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,3)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,3)=cpath2sat(jj)

         endif  ! switch for types for goes10 

        enddo   ! all types for goes10
c gmssat (fourth satellite type)
       elseif(csatid(i).eq.'gmssat')then
        do j=1,ntypes(i)
c     format type 4 (rll)
         jj=jj+1
         if(csattypes(jj).eq.'rll')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,4)=1
            c_channel_types(1,4,4)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wvp')then
            ICHANNELS(3,4,4)=1
            c_channel_types(3,4,4)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,4)=1
            c_channel_types(4,4,4)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,4)=1
            c_channel_types(5,4,4)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,4)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,4)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'hko')then

c     format type 2 (hko):  for HKO (JS  and PW Li 3-20-03)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,2,4)=1
            c_channel_types(1,2,4)=cchanneltypes(kk)
            i_msng_sat_flag(1,2,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv ')then
            ICHANNELS(3,2,4)=1
            c_channel_types(3,2,4)=cchanneltypes(kk)
            i_msng_sat_flag(3,2,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'ir1')then
            ICHANNELS(4,2,4)=1
            c_channel_types(4,2,4)=cchanneltypes(kk)
            i_msng_sat_flag(4,2,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,2,4)=1
            c_channel_types(5,2,4)=cchanneltypes(kk)
            i_msng_sat_flag(5,2,4)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,2,4)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'twn')then

c     format type 3 (twn):  for taiwan (JS  and BS Wang 6-7-01)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,4)=1
            c_channel_types(1,3,4)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv ')then
            ICHANNELS(3,3,4)=1 
            c_channel_types(3,3,4)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,3,4)=1  !ir1
            c_channel_types(4,3,4)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,4)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,3,4)=1  !<-- end 4th type
            c_channel_types(5,3,4)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,4)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,3,4)=cpath2sat(jj)
         endif
        enddo  !all types for gmssat

c goes12 (fifth satellite type)
       elseif(csatid(i).eq.'goes12')then
        do j=1,ntypes(i)
         jj=jj+1
c format type 1 (gvr)
         if(csattypes(jj).eq.'gvr')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,5)=1
            c_channel_types(1,1,5)=cchanneltypes(kk)
            i_msng_sat_flag(1,1,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,1,5)=1
            c_channel_types(2,1,5)=cchanneltypes(kk)
            i_msng_sat_flag(2,1,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,1,5)=1
            c_channel_types(3,1,5)=cchanneltypes(kk)
            i_msng_sat_flag(3,1,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,1,5)=1
            c_channel_types(4,1,5)=cchanneltypes(kk)
            i_msng_sat_flag(4,1,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'13u')then
            ICHANNELS(6,1,5)=1
            c_channel_types(6,1,5)=cchanneltypes(kk)
            i_msng_sat_flag(6,1,5)=i_qc_sat_flag(kk)
c      C_CHANNEL_TYPES(5,1,5)='   '  !doesn't appear to be 12u for goes12
           endif
          enddo
          path_to_raw_sat(1:6,1,5)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'wfo')then
          call s_len(cpath2sat(jj),n)
c format type 2 (wfo)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,2,5)=1
            c_channel_types(1,2,5)=cchanneltypes(kk)
            path_to_raw_sat(1,2,5)=cpath2sat(jj)(1:n)//'vis/regclip/'
            i_msng_sat_flag(1,2,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,2,5)=1
            c_channel_types(2,2,5)=cchanneltypes(kk)
            path_to_raw_sat(2,2,5)=cpath2sat(jj)(1:n)//'i39/regclip/'
            i_msng_sat_flag(2,2,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'iwv')then
            ICHANNELS(3,2,5)=1
            c_channel_types(3,2,5)=cchanneltypes(kk)
            path_to_raw_sat(3,2,5)=cpath2sat(jj)(1:n)//'iwv/regclip/'
            i_msng_sat_flag(3,2,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,2,5)=1
            c_channel_types(4,2,5)=cchanneltypes(kk)
            path_to_raw_sat(4,2,5)=cpath2sat(jj)(1:n)//'i11/regclip/'
            i_msng_sat_flag(4,2,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i12')then
            ICHANNELS(5,2,5)=1
            c_channel_types(5,2,5)=cchanneltypes(kk)
            path_to_raw_sat(5,2,5)=cpath2sat(jj)(1:n)//'i12/regclip/'
            i_msng_sat_flag(5,2,5)=i_qc_sat_flag(kk)
           endif
          enddo

         elseif(csattypes(jj).eq.'cdf')then
c format type 3 (cdf)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,5)=1
            c_channel_types(1,3,5)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,3,5)=1
            c_channel_types(2,3,5)=cchanneltypes(kk)
            i_msng_sat_flag(2,3,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,3,5)=1
            c_channel_types(3,3,5)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,3,5)=1
            c_channel_types(4,3,5)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,3,5)=1
            c_channel_types(5,3,5)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,5)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,3,5)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'rll')then
c format type 3 (raw lat lon)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,5)=1
            c_channel_types(1,4,5)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,5)=1
            c_channel_types(2,4,5)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,5)=1
            c_channel_types(3,4,5)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,5)=1
            c_channel_types(4,4,5)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,5)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,5)=1
            c_channel_types(5,4,5)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,5)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,5)=cpath2sat(jj)
         endif
        enddo

       elseif(csatid(i).eq.'goes09')then
c goes09 (sixth satellite type)
        do j=1,ntypes(i)
         jj=jj+1
c format type 1 (gvr)
         if(csattypes(jj).eq.'gvr')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,6)=1
            c_channel_types(1,1,6)=cchanneltypes(kk)
            i_msng_sat_flag(1,1,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,1,6)=1
            c_channel_types(2,1,6)=cchanneltypes(kk)
            i_msng_sat_flag(2,1,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,1,6)=1
            c_channel_types(3,1,6)=cchanneltypes(kk)
            i_msng_sat_flag(3,1,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,1,6)=1
            c_channel_types(4,1,6)=cchanneltypes(kk)
            i_msng_sat_flag(4,1,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,1,6)=1
            c_channel_types(5,1,6)=cchanneltypes(kk)
            i_msng_sat_flag(5,1,6)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,1,6)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'cdf')then
c format type 3 (cdf)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,6)=1
            c_channel_types(1,3,6)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,3,6)=1
            c_channel_types(2,3,6)=cchanneltypes(kk)
            i_msng_sat_flag(2,3,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,3,6)=1
            c_channel_types(3,3,6)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,3,6)=1
            c_channel_types(4,3,6)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,3,6)=1
            c_channel_types(5,3,6)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,6)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,3,6)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'rll')then
c format type 4 (rll)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,6)=1
            c_channel_types(1,4,6)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,6)=1
            c_channel_types(2,4,6)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,6)=1
            c_channel_types(3,4,6)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,6)=1
            c_channel_types(4,4,6)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,6)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,6)=1
            c_channel_types(5,4,6)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,6)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,6)=cpath2sat(jj)

         endif
        enddo

c satellite = goes11
       elseif(csatid(i).eq.'goes11')then
        do j=1,ntypes(i)
         jj=jj+1
c format type 1 (gvr)
         if(csattypes(jj).eq.'gvr')then
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,7)=1
            c_channel_types(1,1,7)=cchanneltypes(kk)
            i_msng_sat_flag(1,1,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,1,7)=1
            c_channel_types(2,1,7)=cchanneltypes(kk)
            i_msng_sat_flag(2,1,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,1,7)=1
            c_channel_types(3,1,7)=cchanneltypes(kk)
            i_msng_sat_flag(3,1,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,1,7)=1
            c_channel_types(4,1,7)=cchanneltypes(kk)
            i_msng_sat_flag(4,1,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'13u')then
            ICHANNELS(6,1,7)=1
            c_channel_types(6,1,7)=cchanneltypes(kk)
            i_msng_sat_flag(6,1,7)=i_qc_sat_flag(kk)
c      C_CHANNEL_TYPES(5,1,7)='   '  !doesn't appear to be 12u for goes11 -> ? maybe there is?
           endif
          enddo
          path_to_raw_sat(1:6,1,7)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'wfo')then
          call s_len(cpath2sat(jj),n)
c format type 2 (wfo)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,2,7)=1
            c_channel_types(1,2,7)=cchanneltypes(kk)
            path_to_raw_sat(1,2,7)=cpath2sat(jj)(1:n)//'vis/regclip/'
            i_msng_sat_flag(1,2,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,2,7)=1
            c_channel_types(2,2,7)=cchanneltypes(kk)
            path_to_raw_sat(2,2,7)=cpath2sat(jj)(1:n)//'i39/regclip/'
            i_msng_sat_flag(2,2,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'iwv')then
            ICHANNELS(3,2,7)=1
            c_channel_types(3,2,7)=cchanneltypes(kk)
            path_to_raw_sat(3,2,7)=cpath2sat(jj)(1:n)//'iwv/regclip/'
            i_msng_sat_flag(3,2,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,2,7)=1
            c_channel_types(4,2,7)=cchanneltypes(kk)
            path_to_raw_sat(4,2,7)=cpath2sat(jj)(1:n)//'i11/regclip/'
            i_msng_sat_flag(4,2,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i12')then
            ICHANNELS(5,2,7)=1
            c_channel_types(5,2,7)=cchanneltypes(kk)
            path_to_raw_sat(5,2,7)=cpath2sat(jj)(1:n)//'i12/regclip/'
            i_msng_sat_flag(5,2,7)=i_qc_sat_flag(kk)
           endif
          enddo

         elseif(csattypes(jj).eq.'cdf')then
c format type 3 (cdf)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,3,7)=1
            c_channel_types(1,3,7)=cchanneltypes(kk)
            i_msng_sat_flag(1,3,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,3,7)=1
            c_channel_types(2,3,7)=cchanneltypes(kk)
            i_msng_sat_flag(2,3,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,3,7)=1
            c_channel_types(3,3,7)=cchanneltypes(kk)
            i_msng_sat_flag(3,3,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,3,7)=1
            c_channel_types(4,3,7)=cchanneltypes(kk)
            i_msng_sat_flag(4,3,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,3,7)=1
            c_channel_types(5,3,7)=cchanneltypes(kk)
            i_msng_sat_flag(5,3,7)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,3,7)=cpath2sat(jj)

         elseif(csattypes(jj).eq.'rll')then
c format type 4 (raw lat lon)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,7)=1
            c_channel_types(1,4,7)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u')then
            ICHANNELS(2,4,7)=1
            c_channel_types(2,4,7)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'wv')then
            ICHANNELS(3,4,7)=1
            c_channel_types(3,4,7)=cchanneltypes(kk)
            i_msng_sat_flag(3,4,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,7)=1
            c_channel_types(4,4,7)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,7)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'12u')then
            ICHANNELS(5,4,7)=1
            c_channel_types(5,4,7)=cchanneltypes(kk)
            i_msng_sat_flag(5,4,7)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,7)=cpath2sat(jj)

         endif
        enddo

c satellite = NOAA POLAR ORBITER NETCDF
       elseif(csatid(i).eq.'noaapo')then
        do j=1,ntypes(i)
         jj=jj+1
c format type 1 (ncp): netcdf polar (FMI's data type). Stored like wfo
         if(csattypes(jj).eq.'ncp')then
          call s_len(cpath2sat(jj),n)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,1,8)=1
            c_channel_types(1,1,8)=cchanneltypes(kk)
            path_to_raw_sat(1,1,8)=cpath2sat(kk)(1:n)//'/vis/'
            i_msng_sat_flag(1,1,8)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i39')then
            ICHANNELS(2,1,8)=1
            c_channel_types(2,1,8)=cchanneltypes(kk)
            path_to_raw_sat(2,1,8)=cpath2sat(kk)(1:n)//'/i39/'
            i_msng_sat_flag(2,1,8)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'i11')then
            ICHANNELS(4,1,8)=1
            c_channel_types(4,1,8)=cchanneltypes(kk)
            path_to_raw_sat(4,1,8)=cpath2sat(kk)(1:n)//'/i11/'
            i_msng_sat_flag(4,1,8)=i_qc_sat_flag(kk)
           endif
          enddo

         elseif(csattypes(jj).eq.'rll')then
c format type 4 (rll)
          do k=1,nchannel(jj)
           kk=kk+1
           if(cchanneltypes(kk).eq.'vis')then
            ICHANNELS(1,4,8)=1
            c_channel_types(1,4,8)=cchanneltypes(kk)
            i_msng_sat_flag(1,4,8)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'4u ')then
            ICHANNELS(2,4,8)=1
            c_channel_types(2,4,8)=cchanneltypes(kk)
            i_msng_sat_flag(2,4,8)=i_qc_sat_flag(kk)
           elseif(cchanneltypes(kk).eq.'11u')then
            ICHANNELS(4,4,8)=1
            c_channel_types(4,4,8)=cchanneltypes(kk)
            i_msng_sat_flag(4,4,8)=i_qc_sat_flag(kk)
           endif
          enddo
          path_to_raw_sat(1:6,4,8)=cpath2sat(jj)
         endif
        enddo
       endif
      enddo
      print*,'Done in conf_isjtkc. Returning to conf_satellite'
      return
      end
