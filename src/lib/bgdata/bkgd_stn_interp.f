cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

      subroutine bkgd_stn_interp(nstns,slat,slon,selev,stime
     &,stn_ht,stn_tp,stn_td,stn_uw,stn_vw,stn_pr,stn_mslp
     &,bkg_status)

      implicit none
      include 'bgdata.inc'

c eventually these parameters should go into bgdata.inc
c     integer mxvars,mxlvls
c     parameter(mxvars=10,mxlvls=10)
c     character*10 cvars(mxvars)
c     integer      levels(mxlvls,mxvars)

c
c
c------------------> BACKGROUND MODEL DESIGNATION <-----------------------------
c
c *** The following variable designates which model to use as the background.
c *** Read from standard input.
c
c        bgmodel = 1 ---> RUC (60 km native grid)
c        bgmodel = 2 ---> ETA (48 km conus-c grid)
c        bgmodel = 3 ---> CWB FA (5km) model
c        bgmodel = 4 ---> SBN Conus211 (Eta or RUC)
c        bgmodel = 5 ---> RUC (40 km native grid)
c        bgmodel = 6 ---> AVN (360 x 181 lat-lon grid)
c        bgmodel = 7 ---> ETA (48 km from grib file)
c        bgmodel = 8 ---> NOGAPS (1.0 deg)
c        bgmodel = 9 ---> NWS Conus (RUC, ETA, NGM, AVN)
c
c
      integer bgmodel
      integer reject_cnt
      integer accepted_files
      integer forecast_length
c
c------------------> GRID DIMENSION SPECIFICATION <-----------------------------
c
c *** The following variables specify the various model grid dimensions.
c *** Read from standard input.
c *** The first line specifies the LAPS grid dimensions and p level specs.
c *** The second line specifies the LAPS cycle time.
c *** The third line specifies the LAPS root path directory.
c *** The fourth line specifies the LAPS domain designation.
c *** Each following line specifies the grid dimensions for each corresponding
c        bgmodel defined above.
c
      integer   nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .          laps_cycle_time,             !LAPS cycle time
     .          nx_bg  ,ny_bg  ,nz_bg,       !Background model grid dimensions
     .          bkg_status                   !status returned from lga_driver
c                                            ! 1=good 0=bad
      real*4    prbot,delpr                  !LAPS bottom and delta pressures

      integer  nstns
      integer  stime(nstns)
      real*4   slat(nstns)
      real*4   slon(nstns)
      real*4   selev(nstns)

      real*4    stn_ht(nstns)
     .         ,stn_tp(nstns)
     .         ,stn_td(nstns)
     .         ,stn_uw(nstns)
     .         ,stn_vw(nstns)
     .         ,stn_pr(nstns)
     .         ,stn_mslp(nstns)
c
c
c------------------> BACKGROUND GRID DATA PATH <--------------------------------
c
c *** The following variables specify the various model data directory paths.
c *** Read from standard input.
c *** Each line contains the background data directory path for each
c        corresponding bgmodel defined above.
c
      character*256 bgpath
      character*256 cfilespec

c
c  This is the max number of paths allowed in nest7grid.parms
c  and should match the value in lib/lapsgrid.f 
c 
      character*256 bgpaths(maxbgmodels)
      integer bgmodels(maxbgmodels), bglen
      integer istatus

      character*13 cvt_i4time_wfo_fname13
      character*13 wfo_fname13_to_fname9
      character*9 a9

c     integer  plvls(mxvars,mxlvls)

      integer i4time_now,i4time_latest
      integer i4time_now_lga,i4time_now_gg
      integer max_files,bg_files
      integer itime_inc
      parameter (max_files=600)
      character*256 names(max_files)
      character*256 reject_names(max_files)
      integer nbgmodels
      integer lbgp

      integer nvars
      integer nx,ny

c     integer idum,jdum,kdum
c     integer idims(mxvars,mxlvls)
c     real centralLat,centralLon,rlat00,rlon00
c    +,latNxNy,lonNxNy,latdxdy,londxdy,dx,dy
c
c-------------------------------------------------------------------------------
c
c *** Comments used in writing netcdf files only.
c
      integer istat, i,l, no_infinite_loops
c
c cmodel is really only 12 chars but the SBN netcdf carrys 132
c
      character*132 cmodels(maxbgmodels)
      character*132 cmodel
      integer oldest_forecast, max_forecast_delta
      logical use_analysis, use_systime

c_______________________________________________________________________________
c
c Read information from static/background.nl
c
      call get_laps_cycle_time(laps_cycle_time,istatus)
      if(istatus.ne.1)then
         print*,'Error getting analysis cycle time'
         return
      endif
      call get_background_info(bgpaths,bgmodels,oldest_forecast
     +,max_forecast_delta,forecast_length,use_analysis,cmodels
     +,itime_inc)

      nbgmodels=0
      do i=1,maxbgmodels
         if(index(bgpaths(i),' ').ne.0)then
            nbgmodels=nbgmodels+1
         endif
      enddo

c
c *** Initialize esat table.
c
      call es_ini
c
c *** Get current time from systime.dat
c
      use_systime=.true.
      if(use_systime)then
         call get_systime(i4time_now,a9,bkg_status)
         print*,'Analysis systime: ',a9,' ',i4time_now
      else
         i4time_now = i4time_now_gg()
         print*,'Using i4time now'
      endif

      print*,'Time adjustment calc from namelist itime_inc',
     +' (itime_inc) = ',itime_inc
      i4time_now_lga=i4time_now+(itime_inc*laps_cycle_time)
      call make_fnam_lp(i4time_now_lga,a9,istatus)
      print*,'processing background data for ',a9
      print*

c*****************************************
      bg_files=0
      i=1
      bgmodel = bgmodels(i)
      bgpath =  bgpaths(i)
      call s_len(bgpath,lbgp)
      cmodel =  cmodels(i)
      bkg_status = -99 

      reject_cnt=0
      bg_files=0

      no_infinite_loops=0
      do while((bkg_status.le.0 .and. i.le.nbgmodels)
     +         .and. (no_infinite_loops.le.nbgmodels))


         print *,'HERE:',bkg_status, i, bg_files,reject_cnt

         if(reject_cnt.gt.bg_files )then    !.and. bg_files.gt.0) then
            i=i+1
            bgmodel = bgmodels(i)
            bgpath =  bgpaths(i)
            cmodel =  cmodels(i)
            reject_cnt = 0
         endif

         if (bgmodel .lt. 1 .or. bgmodel .gt. maxbgmodels) then
            print*
            print*,' Cannot proceed with model specification in LGA'
            print*,' Check bgpaths in static/background.nl'
            print*,' LGA process aborted...'
            stop
         endif

c
c this commented subroutine is under development to replace
c get_acceptable_files.  At the moment get_acceptable_files
c is still doing the job even though it is difficult software
c to work with. get_acceptable_files is also used in dprep.f 
c
c        call get_bkgd_files(i4time_now_lga,bgpath,bgmodel


         print*,'Calling get_acceptable_files '
         print*,'bgpath:  ',bgpath(1:lbgp)
         print*,'bgmodel: ',bgmodel
         print*

         call get_acceptable_files(i4time_now_lga,bgpath,bgmodel
     +        ,names,max_files,oldest_forecast,max_forecast_delta
     +        ,use_analysis,bg_files,accepted_files,forecast_length
     +        ,cmodel,nx_bg,ny_bg,nz_bg,reject_names,reject_cnt)


        if(accepted_files.eq.0.and.bg_files.eq.0) then

           print*,'No Acceptable files found for background model:'
           print*,'bgpath =  ',bgpath(1:lbgp)
           print*,'bgmodel = ',bgmodel 

           no_infinite_loops=no_infinite_loops+1
           reject_cnt=reject_cnt+1
           if(i.eq.nbgmodels)bkg_status = 0 

c             if(bg_files.eq.0)then
c                bkg_status = 1
c             endif

        elseif(accepted_files.eq.0.and.reject_cnt.eq.bg_files)then

              reject_cnt=reject_cnt+1

        else

           print *

c          print *, 'Input Parameters'
c          print *, '----------------'
c          print *
c          print *, ' Analysis setup: '
c          print *, laps_cycle_time

           print *



c************************************
c     bg_files=0
c     i=0
c     bkg_status = 0 

c     no_infinite_loops=0
c     do while(bkg_status.le.0 .and. i.le.maxbgmodels
c    +     .and. no_infinite_loops.lt.30)
c        no_infinite_loops=no_infinite_loops+1
c        print *,'HERE:',bkg_status, i, bg_files,reject_cnt
c        if(bg_files.le.reject_cnt) then
c           i=i+1 
c           bgmodel = bgmodels(i)
c           if(bgmodel .eq. 0) goto 965
c           bgpath =  bgpaths(i)
c           call s_len(bgpath,bglen)
c           if(bgmodel .eq. 4)then
c              cfilespec=bgpath(1:bglen)//'/*'
c              call get_latest_file_time(cfilespec,i4time_latest)
c              cfname = cvt_i4time_wfo_fname13(i4time_latest)
c              cfilespec=bgpath(1:bglen)//cfname

c              call get_sbn_dims(cfilespec,cmodel,mxvars,
c    +                  mxlvls,idims,cvars,plvls,ntbg,istatus)
c              if(istatus .ne. 1)then
c                 print*,'Error - get_sbn_dims'
c                 return
c              endif
c              call get_attribute_sbn(cfilespec,centralLat
c    +,centralLon,rlat00,rlon00,latNxNy,lonNxNy,latdxdy,londxdy
c    +,dx,dy,nx,ny,istatus)
c              if(istatus .ne. 1)then
c                 print*,'Error - get_attribute_sbn'
c                 return
c              endif

c           endif
c        endif
c        if (bgmodel .lt. 1 .or. bgmodel .gt. maxbgmodels) then
c           print*,'Bad model specification in LGA, bgmodel =',bgmodel
c           print*,'   LGA process aborted...'
c           stop
c        endif

c        call get_acceptable_files(i4time_now_lga,bgpath,bgmodel
c    +        ,names,max_files,oldest_forecast,max_forecast_delta
c    +        ,use_analysis,bg_files,0,cmodel(i)
c    +        ,nx_bg,ny_bg,nz_bg,reject_files,reject_cnt)


c       if(bg_files.eq.0) then
c          print*,'No Acceptable files found for model: ',bgpath,
c    +          bgmodel 
c          bkg_status = 0
c       else

c          print *, 'bgmodel is set to ', bgmodel
c          print *, ' '
c          print *, 'input parameters'
c          print *,  ' '
c          print *, nx_bg,ny_bg,nz_bg
c          print *, 'bgpath ', bgpath(1:bglen)
c          print *, 'cmodel ',cmodel(i)
c970       continue

c          call get_bkgd_mdl_info(bgmodel,cmodel,cfname
c    &,mxvars,mxlvls,nvars,cvars,levels,gproj,nx_bg,ny_bg,nz_bg
c    &,dlat,dlon,Lat1,Lat2,Lon0,sw,ne,istatus)
c          if(istatus.ne.1)then
c             print*,'Error getting background model information'
c             return
c          endif
c
c *** Call bkgd driver.
c
           call  bkgd_stn_interp_sub(
     + nstns,slat,slon,selev,stime
     +,bgmodel,bgpath,names,cmodel,max_files,bg_files
     +,laps_cycle_time,reject_names,reject_cnt,accepted_files
     +,nx_bg, ny_bg, i4time_now_lga, stn_ht,stn_tp
     +,stn_td,stn_uw,stn_vw,stn_pr,stn_mslp, bkg_status)

c          if(bkg_status.lt.0) then
c             do l = 1,bg_files
c                reject_cnt=reject_cnt+1
c                reject_files(reject_cnt)=names(l)       ! (-bkg_status), name was names(-bkg_status)
c             enddo
c             bg_files = 0                               ! reset the counter
c          endif

ccc - to here
        endif
        
      enddo

      print*,'station interp done'
      print*
      do i=1,nstns
         print*,'istn/ht/t/td/pr/mslp ',i,stn_ht(i)
     &,stn_tp(i),stn_td(i),stn_pr(i),stn_mslp(i)
      enddo
      print*

      if(no_infinite_loops.gt.oldest_forecast) then
         print*,'ERROR: Infinite loop condition found'
      endif
c
      if(bkg_status.eq.0) goto 965
      print *,'bkgd...Normal completion.'
      stop
c
c *** Error trapping.
c
900   continue
      print *,'Error reading background model number.'
      stop
c
910   continue
      print *,'Error reading LAPS grid information.'
      stop
c
920   continue
      print *,'Error reading LAPS cycle time.'
      stop
c
930   continue
      print *,'Error reading LAPS root directory path.'
      stop
c
940   continue
      print *,'Error reading LAPS domain specification.'
      stop
c
950   continue
      print *,'Error reading background model grid dimensions.'
      stop
c
960   continue
      print *,'Error reading background data directory paths.'
      stop
 965  continue
      print *,'No acceptable background model found.'
      
c
      end
c
c===============================================================================
c
      subroutine bkgd_stn_interp_sub(
     + nstns,slat,slon,selev,stime
     +,bgmodel,bgpath,bg_names,cmodel,max_files,bg_files
     +,laps_cycle_time,reject_names,reject_cnt,nbg       !<-- nbg is same as accepted_files
     +,nx_bg,ny_bg,i4time_now,stn_ht,stn_tp,stn_td
     +,stn_uw,stn_vw,stn_pr,stn_mslp, bkg_status)

c
      implicit none

c     include 'bgdata.inc'
c
      integer   nbg
      integer   bg_files
      integer   nzbg_ht,nzbg_tp,nzbg_sh
     +         ,nzbg_uv,nzbg_ww
     +         ,nz_laps

      integer   nstns,laps_cycle_time
     +         ,nx_bg  ,ny_bg           !Background model horiz grid dimensions
     +         ,max_files, bkg_status
     +         ,bgmodel

      character*256 bgpath
     +             ,bg_names(max_files)
     +             ,reject_names(max_files)

      character*132 cmodel

      character*256 fname_bg(max_files)
      real*4    diff_mn,diff_mx
      integer   icnt
      integer   i_mx,j_mx
      integer   i_mn,j_mn

      logical   l_done_this
      logical   allocated_for_vinterp
c
c *** Background model grid data.
c
      real, allocatable  :: prbg_sfc(:,:)
      real, allocatable  :: uwbg_sfc(:,:)
      real, allocatable  :: vwbg_sfc(:,:)
      real, allocatable  :: shbg_sfc(:,:)
      real, allocatable  :: tpbg_sfc(:,:)
      real, allocatable  :: htbg_sfc(:,:)
      real, allocatable  :: mslpbg(:,:)
c
c *** 3D background arrays.
c
      real, allocatable  :: prbght(:,:,:)
      real, allocatable  :: prbgsh(:,:,:)
      real, allocatable  :: prbguv(:,:,:)
      real, allocatable  :: prbgww(:,:,:)
      real, allocatable  :: htbg(:,:,:)
      real, allocatable  :: tpbg(:,:,:)
      real, allocatable  :: shbg(:,:,:)
      real, allocatable  :: uwbg(:,:,:)
      real, allocatable  :: vwbg(:,:,:)
      real, allocatable  :: wwbg(:,:,:)

      real, allocatable  :: hgt1dz(:)
     .                     ,tmp1dz(:)
     .                     ,ssh1dz(:)
     .                     ,uuw1dz(:)
     .                     ,vvw1dz(:)
     .                     ,www1dz(:)

      real, allocatable  :: pr1dz_ht(:)
     .                     ,pr1dz_sh(:)
     .                     ,pr1dz_uv(:)
     .                     ,pr1dz_ww(:)

      real, allocatable  :: hgt1dzl(:)
     .                     ,tmp1dzl(:)
     .                     ,ssh1dzl(:)
     .                     ,uuw1dzl(:)
     .                     ,vvw1dzl(:)
     .                     ,www1dzl(:)


      real, allocatable  ::  pr_laps(:)

c these are the station interpolated data at the background times
      real*4    tp_sfc(nstns,nbg)
     .         ,sh_sfc(nstns,nbg)
     .         ,uw_sfc(nstns,nbg)
     .         ,vw_sfc(nstns,nbg)
     .         ,pr_sfc(nstns,nbg)
     .         ,mslp(nstns,nbg)
     .         ,alt_sfc(nstns,nbg)

c these are the time interpolated (from the background times) stn info
      real*4    stn_ht(nstns)
     .         ,stn_tp(nstns)
     .         ,stn_td(nstns)
     .         ,stn_uw(nstns)
     .         ,stn_vw(nstns)
     .         ,stn_pr(nstns)
     .         ,stn_mslp(nstns)
     .         ,stn_alt(nstns)

      real*4    weight
      real*4    rmissingflag
      real*4    dlat,dlon,Lat1,Lat2,Lon0,Lat0
      real*4    sw(2),ne(2)
      real*4    dx,dy,cenlon,cenlat

c input station information
      real*4    slat(nstns)
     .         ,slon(nstns)
     .         ,selev(nstns)
     .         ,grx(nstns)     !calculated by init_hinterp
     .         ,gry(nstns)     !      "
      integer   stime(nstns)   !preferrably the i4 time of stn ob

      integer   ct,n,newfcst
     +         ,ihour, itstatus
     +         ,lga_files,lga_times(max_files)
     +         ,lga_valid
     +         ,i4time_lga_valid(max_files),i4time_now
     +         ,bg_times(max_files)
     +         ,i4time_bg_valid(max_files)
     +         ,bg_valid(max_files)
     +         ,valid_bg(max_files),time_bg(max_files)
     +         ,bgvalid, lenp
     +         ,i,ic,ii,j,jj,k,kk,l,ldl,lf
     +         ,istatus,istatus_prep
     +         ,init_timer,reject_cnt

      integer   nf,bglen
c
      character*256 names(max_files)
      character*13  fname13,fname9_to_wfo_fname13
      character*9   wfo_fname13_to_fname9
      character*256 cfname
      character*2   gproj
      character*200 fullname
      character*256 outdir
      character*4   af_bg(max_files)
      character*1   cgrddef

      interface

       subroutine read_bgdata(nx_bg,ny_bg
     +,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,ctype
     +,bgpath,fname_bg,af_bg,fullname,cmodel,bgmodel
     +,prbght,prbgsh,prbguv,prbgww
     +,htbg,tpbg,uwbg,vwbg,shbg,wwbg
     +,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +,uwbg_sfc,vwbg_sfc,mslpbg,istatus)
c
         real  :: prbg_sfc(nx_bg,ny_bg)
         real  :: uwbg_sfc(nx_bg,ny_bg)
         real  :: vwbg_sfc(nx_bg,ny_bg)
         real  :: shbg_sfc(nx_bg,ny_bg)
         real  :: tpbg_sfc(nx_bg,ny_bg)
         real  :: htbg_sfc(nx_bg,ny_bg)
         real  :: mslpbg(nx_bg,ny_bg)
c
         real  :: prbght(nx_bg,ny_bg,nzbg_ht)
         real  :: prbgsh(nx_bg,ny_bg,nzbg_sh)
         real  :: prbguv(nx_bg,ny_bg,nzbg_uv)
         real  :: prbgww(nx_bg,ny_bg,nzbg_ww)
         real  :: htbg(nx_bg,ny_bg,nzbg_ht)
         real  :: tpbg(nx_bg,ny_bg,nzbg_tp)
         real  :: shbg(nx_bg,ny_bg,nzbg_sh)
         real  :: uwbg(nx_bg,ny_bg,nzbg_uv)
         real  :: vwbg(nx_bg,ny_bg,nzbg_uv)
         real  :: wwbg(nx_bg,ny_bg,nzbg_ww)

         character*4   af_bg
         character*5   ctype
         character*132 cmodel
         character*256 bgpath
         character*256 fname_bg
         character*200 fullname
         integer       bgmodel
         integer       nx_bg
         integer       ny_bg
         integer       nzbg_ht
         integer       nzbg_tp
         integer       nzbg_sh
         integer       nzbg_uv
         integer       nzbg_ww
         integer       istatus
       end subroutine

       subroutine vinterp(nz_laps,nx,ny,
     .  nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww,
     .  prlaps, prbght,prbgsh,prbguv,prbgww,
     .  htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .  htvi,tpvi,shvi,uwvi,vwvi,wwvi)

         integer nx,ny
         integer nzbg_ht
         integer nzbg_tp
         integer nzbg_sh
         integer nzbg_uv
         integer nzbg_ww

         integer nz_laps

         real*4  ::  prbght(nx,ny,nzbg_ht)
         real*4  ::  prbgsh(nx,ny,nzbg_sh)
         real*4  ::  prbguv(nx,ny,nzbg_uv)
         real*4  ::  prbgww(nx,ny,nzbg_ww)
         real*4  ::  tpbg(nx,ny,nzbg_tp)
         real*4  ::  htbg(nx,ny,nzbg_ht)
         real*4  ::  shbg(nx,ny,nzbg_sh)
         real*4  ::  uwbg(nx,ny,nzbg_uv)
         real*4  ::  vwbg(nx,ny,nzbg_uv)
         real*4  ::  wwbg(nx,ny,nzbg_ww)


         real*4  ::  tpvi(nx,ny,nz_laps)
         real*4  ::  htvi(nx,ny,nz_laps)
         real*4  ::  shvi(nx,ny,nz_laps)
         real*4  ::  uwvi(nx,ny,nz_laps)
         real*4  ::  vwvi(nx,ny,nz_laps)
         real*4  ::  wwvi(nx,ny,nz_laps)
c
         real*4  ::  prlaps(nz_laps)
 
       end subroutine

       subroutine get_bkgd_mdl_info(bgmodel
     &,cmodel,fullname,nx,ny,nzbg_ht,nzbg_tp
     &,nzbg_sh,nzbg_uv,nzbg_ww
     &,gproj,dlat,dlon,centrallat,centrallon,dx,dy
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

        character*200 fullname
        character*132 cmodel
        character*2   gproj
        character*1   cgrddef
        integer       istatus
        integer       bgmodel
        integer       nx,ny
        integer       nzbg_ht
        integer       nzbg_tp
        integer       nzbg_sh
        integer       nzbg_uv
        integer       nzbg_ww
        real          Lat1
        real          Lat0,Lon0
        real          sw(2),ne(2)
        real          centrallat
        real          centrallon
        real          dx,dy

       end subroutine

       subroutine init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .     lat,lon,grx,gry,bgmodel,cmodel)

         character  gproj*2
         character  cmodel*132
         integer nx_bg,ny_bg,nx_laps,ny_laps,bgmodel
         real*4 lat(nx_laps,ny_laps)
     .         ,lon(nx_laps,ny_laps)
     .         ,grx(nx_laps,ny_laps)
     .         ,gry(nx_laps,ny_laps)

       end subroutine

       subroutine hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz
     .,grx,gry,fvi,flaps,bgmodel)

         integer nx_bg,ny_bg,nx_laps,ny_laps,nz,bgmodel
         real*4 fvi(nx_bg,ny_bg,nz)
     .         ,flaps(nx_laps,ny_laps,nz)
     .         ,grx(nx_laps,ny_laps)
     .         ,gry(nx_laps,ny_laps)

       end subroutine

      end interface

c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c

      print *,'in bkgd_stn_interp sub'
      print*,'Interpolating for ',nstns,' stations'
      bkg_status=0

      call get_r_missing_data(rmissingflag,istatus)
      if(istatus.ne.1)then
          print*,'Error getting r_missing_data'
          return
      endif

      call s_len(bgpath,lenp)
      if(bgpath(lenp:lenp).ne.'/')then
         lenp=lenp+1
         bgpath(lenp:lenp)='/'
      endif
      fullname=bgpath(1:lenp)//bg_names(1)

      call get_bkgd_mdl_info(bgmodel,cmodel,fullname
     &,nx_bg,ny_bg,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv
     &,nzbg_ww ,gproj,dlat,dlon,cenlat,cenlon,dx,dy
     &,Lat0,Lat1,Lon0,sw,ne,cgrddef,istatus)

      if(istatus.ne.1)then
         print*,'Error getting background model information'
         return
      endif
c
c allocate model arrays
c
      allocate (htbg(nx_bg,ny_bg,nzbg_ht))
      allocate (tpbg(nx_bg,ny_bg,nzbg_tp))
      allocate (shbg(nx_bg,ny_bg,nzbg_sh))
      allocate (uwbg(nx_bg,ny_bg,nzbg_uv))
      allocate (vwbg(nx_bg,ny_bg,nzbg_uv))
      allocate (wwbg(nx_bg,ny_bg,nzbg_ww))
      allocate (prbght(nx_bg,ny_bg,nzbg_ht))
      allocate (prbguv(nx_bg,ny_bg,nzbg_uv))
      allocate (prbgsh(nx_bg,ny_bg,nzbg_sh))
      allocate (prbgww(nx_bg,ny_bg,nzbg_ww))

      allocate (htbg_sfc(nx_bg,ny_bg))
      allocate (prbg_sfc(nx_bg,ny_bg))
      allocate (shbg_sfc(nx_bg,ny_bg))
      allocate (uwbg_sfc(nx_bg,ny_bg))
      allocate (vwbg_sfc(nx_bg,ny_bg))
      allocate (tpbg_sfc(nx_bg,ny_bg))
      allocate (mslpbg(nx_bg,ny_bg))

      allocate (hgt1dz(nzbg_ht))
      allocate (tmp1dz(nzbg_tp))
      allocate (ssh1dz(nzbg_sh))
      allocate (uuw1dz(nzbg_uv))
      allocate (vvw1dz(nzbg_uv))
      allocate (www1dz(nzbg_ww))


c
c initialize projection common block and return need info
c
      call init_gridconv_cmn(gproj,nx_bg,ny_bg,nzbg_ht
     &,dlat,dlon,cenlat,cenlon,Lat0,Lat1,Lon0
     &,sw(1),sw(2),ne(1),ne(2),cgrddef,istatus)


      do j=1,max_files
         names(j)=bg_names(j)
      enddo

c add in here the determination of which two model files are
c needed such that we can time interpolate
      do j=1,nbg
         call i4time_fname_lp(names(j)(1:9),bg_times(j),istatus)
         read(bg_names(j)(12:13),'(i2)')ihour
         bg_valid(j)=ihour*3600
         i4time_bg_valid(j)=bg_times(j)+bg_valid(j)
         af_bg(j)=names(j)(10:13)
      enddo

      allocate (pr1dz_ht(nzbg_ht))
      allocate (pr1dz_sh(nzbg_tp))
      allocate (pr1dz_uv(nzbg_sh))
      allocate (pr1dz_ww(nzbg_uv))

      allocated_for_vinterp=.false.
c
c ****** Read background model data for nbg times (=2 or more).
c
      do nf = 1,nbg

        fullname = bgpath(1:lenp)//bg_names(nf)

        call read_bgdata(nx_bg,ny_bg
     +    ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +    ,'lapsb',bgpath,fname_bg(nf),af_bg(nf)
     +    ,fullname,cmodel,bgmodel
     +    ,prbght,prbgsh,prbguv,prbgww
     +    ,htbg,tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg,istatus_prep)

        if (istatus_prep .ne. 0) then

          print *,'Background model data not returned from ',
     +'read_bgdata: ',bgpath(1:lenp)//fname_bg(nf)
          print *,'Process aborted for this file.'

c         convert to wfo if necessary

          print*
          print*,'Updating reject information: '

          if(bgmodel.eq.4)then
             print*,'WFO Update'
             fname13=wfo_fname13_to_fname9(fname_bg(nf))
             fname13=fname13(1:9)//af_bg(nf)
             do ii=1,nbg
                if(fname13.eq.bg_names(ii))then
                   reject_cnt=reject_cnt+1
                   reject_names(reject_cnt)=bg_names(ii)
                   print*,'reject_cnt/reject_names'
                   print*,'cnt/time: ',reject_cnt
     +,reject_names(reject_cnt)

                endif
             enddo
          else
             print*
             do ii=1,nbg
                if(fname_bg(nf).eq.bg_names(ii))then
                   reject_cnt=reject_cnt+1
                   reject_names(reject_cnt)=bg_names(ii)

                   print*,'reject_cnt/reject_names'
                   print*,'cnt/time: ',reject_cnt
     +,reject_names(reject_cnt)
                endif
             enddo
          endif
          bkg_status= -nf

          deallocate(htbg, tpbg, shbg, uwbg, vwbg, wwbg
     +,prbght, prbguv, prbgsh, prbgww, htbg_sfc, prbg_sfc
     +,shbg_sfc, uwbg_sfc, vwbg_sfc, tpbg_sfc, mslpbg)

 
        else   !processing the file

c        endif
c
         itstatus=init_timer()
c
c determine grid-x and grid-y (ri and rj) for station locations
c
         l_done_this = .false.
         if(.not.l_done_this)then
             do i = 1,nstns
                call init_hinterp(nx_bg,ny_bg,1,1,gproj
     +,slat(i),slon(i),grx(i),gry(i),bgmodel,cmodel)
             enddo
             l_done_this = .true.
         endif


c
c ****** Horizontally interpolate background data to station location (lat/lon).
c
         do i = 1,nstns


c ****** Horizontally interpolate background surface data to station location. 
c
         if(bgmodel.ne.1.and.bgmodel.ne.9)then

c need sfc T for subroutine sfcbkgd as first guess estimate.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     +        grx(i),gry(i),tpbg_sfc,tp_sfc(i,nf),bgmodel)
c simple interpolation to get the mslp.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     +        grx(i),gry(i),mslpbg,mslp(i,nf),bgmodel)
c need sfc moisture variable for subroutine sfcbkgd as first guess estimate.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     +        grx(i),gry(i),shbg_sfc,sh_sfc(i,nf),bgmodel)

c we need the 3D data interpolated to the station location.
c we'll need one for each pr field depending on the model!!!!

            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_ht,
     +        grx(i),gry(i),prbght(1,1,1),pr1dz_ht,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_sh,
     +        grx(i),gry(i),prbgsh(1,1,1),pr1dz_sh,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_uv,
     +        grx(i),gry(i),prbguv(1,1,1),pr1dz_uv,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_ww,
     +        grx(i),gry(i),prbgww(1,1,1),pr1dz_ww,bgmodel)


            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_ht,
     +        grx(i),gry(i),htbg(1,1,1),hgt1dz,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_tp,
     +        grx(i),gry(i),tpbg(1,1,1),tmp1dz,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_sh,
     +        grx(i),gry(i),shbg(1,1,1),ssh1dz,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_uv,
     +        grx(i),gry(i),uwbg(1,1,1),uuw1dz,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_uv,
     +        grx(i),gry(i),vwbg(1,1,1),vvw1dz,bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nzbg_ww,
     +        grx(i),gry(i),wwbg(1,1,1),www1dz,bgmodel)

c
c first interpolate vertically to LAPS pressure levels
c this is needed for the call to sfcbkgd which requires all vars
c to be on the same set of levels.

            if(.not.allocated_for_vinterp)then
               call get_laps_dimensions(nz_laps,istatus)
               allocate (hgt1dzl(nz_laps)
     +               ,tmp1dzl(nz_laps)
     +               ,ssh1dzl(nz_laps)
     +               ,pr_laps(nz_laps)
     +               ,uuw1dzl(nz_laps)
     +               ,vvw1dzl(nz_laps)
     +               ,www1dzl(nz_laps))
               allocated_for_vinterp=.true.

               call get_pres_1d(i4time_now,nz_laps,pr_laps,istatus)
               do k = 1,nz_laps
                  pr_laps(k)=pr_laps(k)/100.
               enddo
            endif

c questionable whether we need this... if interpolating upper
c air obs we should use data on model surfaces directly?!

            print*,'Call vinterp '

            call vinterp(nz_laps,1,1
     +  ,nzbg_ht,nzbg_tp,nzbg_sh,nzbg_uv,nzbg_ww
     +  ,pr_laps,pr1dz_ht,pr1dz_sh,pr1dz_uv,pr1dz_ww
     +  ,hgt1dz,tmp1dz,ssh1dz,uuw1dz,vvw1dz,www1dz
     +  ,hgt1dzl,tmp1dzl,ssh1dzl,uuw1dzl,vvw1dzl,www1dzl)
c
c
c *****
c all calculations from here on are for sfc station variables.
c
c compute sfc values from model using station elev and model info.
c note1: sh_sfc can have q, rh, or Td depending on which model background
c       is used.  sfcbkgd deals with this and always returns
c       Td (deg K) (in array sh_sfc)
c note2: this subroutine may only work for variables on LAPS p surfaces
c        we should be interpolating directly from model levels; however,
c        some models do not provide sh on the same number of vertical
c        levels as T and hgt.

            call sfcbkgd(bgmodel, tmp1dzl,ssh1dzl,hgt1dzl
     +  ,tp_sfc(i,nf), sh_sfc(i,nf), selev(i)
     +  ,pr_laps, 1,1, nz_laps, pr_sfc(i,nf))
 
c check for T > Td before sfc p computation. Due to large scale
c interpolation we can have slightly larger (fractional) Td than T.
c
            call tdcheck(1,1,sh_sfc(i,nf),tp_sfc(i,nf),
     &icnt,i_mx,j_mx,diff_mx,diff_mn)

            if(icnt .gt. 0) then
               print *,' Dewpt greater than temp at station'
               print*,'Max diff of ',diff_mx
               print*,'Min diff of ',diff_mn
            endif
c
c..... Do the winds.
c      Note: arrays are on model vertical grid, not laps pressure.
c            this is ok if u/v are on the same model levels as hgt.
c
            if(selev(i).gt.hgt1dz(1))then
               call interp_to_sfc(selev(i),uuw1dz,hgt1dz,
     +  1,1,nzbg_uv,rmissingflag,uw_sfc(i,nf))
               call interp_to_sfc(selev(i),vvw1dz,hgt1dz,
     +  1,1,nzbg_uv,rmissingflag,vw_sfc(i,nf))
            else
c
c we probably can do better here if elev is under model terrain
c but for now we will set wind equal to lowest level
               uw_sfc(i,nf)=uuw1dz(1)
               vw_sfc(i,nf)=vvw1dz(1)
            endif
c
c now compute alt setting and other variables needed at ea stn loc.
c
         endif
         enddo
        endif
      enddo

      bkg_status = 1
c
c time interpolate between  model time if necessary
c--------------------------------------------------
c
      if(nbg.gt.1)then

         newfcst=bg_times(nbg-1)-(i4time_bg_valid(nbg)-i4time_now)
         weight=float(i4time_bg_valid(nbg)-i4time_now)/
     &float(i4time_bg_valid(nbg)-i4time_bg_valid(nbg-1))

      endif

      print*,'Time interp weight = ',weight

      do i=1,nbg
         call s_len(bg_names(i),j)
         print*,'bg_name: ',i,bg_names(i)(1:j)
      enddo
c
c *** Determine if interpolated stn data needs
c     linear time interpolation to satisfy analysis time.

      if(nbg.gt.1)then

         print*,i,bg_times(nbg),bg_times(nbg-1),
     +bg_valid(nbg),bg_valid(nbg-1),laps_cycle_time,
     +i4time_bg_valid(nbg),i4time_bg_valid(nbg-1)

         do i=1,nstns
            stn_tp(i)= (1.- weight)*tp_sfc(i,nbg) +
     +                      weight *tp_sfc(i,nbg-1)
            stn_td(i)= (1.- weight)*sh_sfc(i,nbg) +
     +                      weight *sh_sfc(i,nbg-1)
            stn_uw(i)= (1.- weight)*uw_sfc(i,nbg) +
     +                      weight *uw_sfc(i,nbg-1)
            stn_vw(i)= (1.- weight)*vw_sfc(i,nbg) +
     +                      weight *vw_sfc(i,nbg-1)
            stn_pr(i)= (1.- weight)*pr_sfc(i,nbg) +
     +                      weight *pr_sfc(i,nbg-1)
            stn_mslp(i)= (1.- weight)*mslp(i,nbg) +
     +                        weight *mslp(i,nbg-1)
            stn_alt(i)= (1.- weight)*alt_sfc(i,nbg) +
     +                       weight *alt_sfc(i,nbg-1)
         enddo
      else
         print*,'No time interp when nbg < 2'
         print*
         do i=1,nstns
            stn_tp(i)  = tp_sfc(i,nbg)
            stn_td(i)  = sh_sfc(i,nbg)
            stn_uw(i)  = uw_sfc(i,nbg)
            stn_vw(i)  = vw_sfc(i,nbg)
            stn_pr(i)  = pr_sfc(i,nbg)
            stn_mslp(i)= mslp(i,nbg)
            stn_alt(i) = alt_sfc(i,nbg)
         enddo
      endif

      bkg_status=1

      return
      end
c
c---------------------------------------------------------------------------
c
      subroutine checknan_3d(x,ni,nj,nk,nan_flag)
c
c     Routine to check a real 3-d array for NaN's.
c
      integer ni,nj,nk
      real*4 x(ni,nj,nk)
c
      nan_flag = 1
c
      do k=1,nk
      do j=1,nj
      do i=1,ni
         if( nan( x(i,j,k) ) .eq. 1) then
            print *,' ** ERROR. Found a NaN at ', i, j, k
            nan_flag = -1
            return
         endif
      enddo !i
      enddo !j
      enddo !k
c
      return
      end
c
c
      subroutine tdcheck(nx_laps,ny_laps,sh_sfc,tp_sfc,
     &icnt,i_mx,j_mx,dmax,dmin)

      integer icnt,i_mx,j_mx

      real sh_sfc(nx_laps,ny_laps)
      real tp_sfc(nx_laps,ny_laps)
      real dmax,dmin

      icnt=0
      diff_mx = -1.e30
      diff_mn =  1.e30
      i_mx = 0
      j_mx = 0
      i_mn = 0
      j_mn = 0
      do j=1,ny_laps
      do i=1,nx_laps
         if(sh_sfc(i,j).gt.tp_sfc(i,j))then
            diff = sh_sfc(i,j) - tp_sfc(i,j)
            sh_sfc(i,j)=tp_sfc(i,j)
            icnt=icnt+1
            if(diff .gt. diff_mx) then
               diff_mx = diff
               i_mx = i
               j_mx = j
            endif
            if(diff .lt. diff_mn) then
               diff_mn = diff
               i_mn = i
               j_mn = j
            endif
         endif
      enddo
      enddo

      return
      end
