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
      program lga
      implicit none
      include 'bgdata.inc'
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
      integer nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .          laps_cycle_time,             !LAPS cycle time
     .          nx_bg  ,ny_bg  ,nz_bg,       !Background model grid dimensions
     .          lga_status                   !status returned from lga_driver
c                                            ! 1=good 0=bad
      integer   np
c
      character*255 lapsroot                 !LAPS root path
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
c
      character*9 a9
      integer i4time_now,i4time_latest
      integer i4time_now_lga,i4time_now_gg
      integer max_files,bg_files
      integer itime_inc
      parameter (max_files=900)
      character*256 names(max_files)
      character*256 reject_files(max_files)
      integer reject_cnt
      data reject_cnt/0/
      integer ntbg
      integer idum,jdum,kdum
c
c-------------------------------------------------------------------------------
c
c *** Comments used in writing netcdf files only.
c
      integer istat, i,l, no_infinite_loops
c
c cmodel is really only 12 chars but the SBN netcdf carrys 132
c
      character*31  laps_domain_file
      character*132 cmodel(maxbgmodels)
      character*255 cfname
      integer oldest_forecast, max_forecast_delta
      logical use_analysis, use_systime
      logical time_plus_one,time_minus_one
c_______________________________________________________________________________
c
c Read information from static/nest7grid.parms

      call get_grid_dim_xy(nx_laps,ny_laps,istatus)
      call get_laps_dimensions(nz_laps,istatus)
      call get_laps_cycle_time(laps_cycle_time,istatus)
c
c Read information from static/background.nl
c
      call get_background_info(bgpaths,bgmodels,oldest_forecast
     +,max_forecast_delta,use_analysis,cmodel,itime_inc)

c
c *** Initialize esat table.
c
      call es_ini
c
c *** Get current time from systime.dat
c
      use_systime=.true.
      if(use_systime)then
         call get_systime(i4time_now,a9,lga_status)
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

      time_plus_one=.true.
      time_minus_one=.true.


      bg_files=0
      i=0
      lga_status = 0 

      no_infinite_loops=0
      do while(lga_status.le.0 .and. i.le.maxbgmodels
     +     .and. no_infinite_loops.lt.30)
         no_infinite_loops=no_infinite_loops+1
         print *,'HERE:',lga_status, i, bg_files,reject_cnt
         if(bg_files.le.reject_cnt) then
            i=i+1 
            bgmodel = bgmodels(i)
            if(bgmodel .eq. 0) goto 965
            bgpath =  bgpaths(i)
            call s_len(bgpath,bglen)
            if(bgmodel .eq. 4)then
               cfilespec=bgpath(1:bglen)//'/*'
               call get_latest_file_time(cfilespec,i4time_latest)
               cfname = cvt_i4time_wfo_fname13(i4time_latest)
               call get_sbn_dims(bgpath,cfname,idum,jdum,kdum
     .,ntbg)
            endif
         endif
         if (bgmodel .lt. 1 .or. bgmodel .gt. maxbgmodels) then
            print*,'Bad model specification in LGA, bgmodel =',bgmodel
            print*,'   LGA process aborted...'
            stop
         endif

         call get_acceptable_files(i4time_now_lga,bgpath,bgmodel
     +        ,names,max_files,oldest_forecast,max_forecast_delta
     +        ,use_analysis,bg_files,0,cmodel(i),ntbg
     +        ,nx_bg,ny_bg,nz_bg,reject_files,reject_cnt)


        if(bg_files.eq.0) then
           print*,'No Acceptable files found for model: ',bgpath,
     +          bgmodel 
           lga_status = 0
        else

           print *, 'bgmodel is set to ', bgmodel
           print *, ' '
           print *, 'Input Parameters'
           print *, '----------------'
           print *, ' Analysis dimensions: '
           print *, nx_laps,ny_laps,nz_laps
           print *, laps_cycle_time
           print *
           print *, ' Background dimensions: '
           print *, nx_bg,ny_bg,nz_bg
           print *, 'bgpath ', bgpath(1:bglen)
           print *, 'cmodel ',cmodel(i)
           print *
 970       continue
c
c *** Call lga driver.
c
           call lga_driver(nx_laps,ny_laps,nz_laps,
     .          laps_cycle_time,lapsroot,laps_domain_file,
     .          bgmodel,bgpath,names,cmodel(i),max_files,bg_files,
     .          nx_bg, ny_bg, nz_bg, i4time_now_lga, lga_status)

           if(lga_status.lt.0) then
              do l = 1,bg_files
                 reject_cnt=reject_cnt+1
                 reject_files(reject_cnt)=names(l)       ! (-lga_status), name was names(-lga_status)
              enddo
              bg_files = 0                               ! reset the counter
           endif
c
c these constructs force t-1 and t+1 cycle time background generation.
c these should be removed when lga runs outside sched.pl
c
           if(lga_status.eq.1.and.time_minus_one)then
              lga_status = 0
              i4time_now_lga = i4time_now-laps_cycle_time
              time_minus_one = .false.
              print*
              print*,'Start time-minus-one cycle time interpolation'
              print*
           endif

           if(lga_status.eq.1.and.time_plus_one)then
              lga_status = 0
              i4time_now_lga = i4time_now + laps_cycle_time
              time_plus_one = .false.
              print*
              print*,'Start time-plus-one cycle time interpolation'
              print*
           endif

        endif
        
      enddo
      if(no_infinite_loops.gt.oldest_forecast) then
         print*,'ERROR: LGA infinite loop condition found'
      endif
c
      if(lga_status.eq.0) goto 965
      print *,'LGA...Normal completion.'
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
      subroutine lga_driver(nx_laps,ny_laps,nz_laps,
     .     laps_cycle_time,lapsroot,laps_domain_file,
     .     bgmodel,bgpath,bg_names,cmodel,max_files,bg_files,
     .     nx_bg,ny_bg,nz_bg,i4time_now, lga_status)

c
      implicit none
      include 'bgdata.inc'
c
      integer nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .          nx_bg  ,ny_bg  ,nz_bg,       !Background model grid dimensions
     .          max_files, lga_status,
     .          laps_cycle_time,
     .          bgmodel,nbg
      integer icnt
      integer i_mx, i_mn, j_mx, j_mn, nan_flag
      real diff, diff_mx, diff_mn
      real lon0,lat1,lat2

      character*(*) lapsroot
     +             ,laps_domain_file
     +             ,bgpath
     +             ,bg_names(max_files)
     +             ,cmodel
      
      integer warncnt
      
c
c *** Background model grid data.
c
      real*4    prbg(nx_bg,ny_bg,nz_bg),     !Pressure (mb)
     .          htbg(nx_bg,ny_bg,nz_bg),     !Height (m)
     .          tpbg(nx_bg,ny_bg,nz_bg),     !Temperature (K)
     .          shbg(nx_bg,ny_bg,nz_bg),     !Specific humidity (kg/kg)
     .          uwbg(nx_bg,ny_bg,nz_bg),     !U-wind (m/s)
     .          vwbg(nx_bg,ny_bg,nz_bg),     !V-wind (m/s)
     .          mslpbg(nx_bg,ny_bg),         !mslp  (mb)
     .          wwbg(nx_bg,ny_bg,nz_bg),     !W-wind (pa/s)
     .          htbg_sfc(nx_bg,ny_bg),
     .          prbg_sfc(nx_bg,ny_bg), 
     .          shbg_sfc(nx_bg,ny_bg), 
     .          uwbg_sfc(nx_bg,ny_bg), 
     .          vwbg_sfc(nx_bg,ny_bg), 
     .          tpbg_sfc(nx_bg,ny_bg)


c
c *** Background data vertically interpolated to LAPS isobaric levels.
c
      real*4    htvi(nx_bg,ny_bg,nz_laps),   !Height (m)
     .          tpvi(nx_bg,ny_bg,nz_laps),   !Temperature (K)
     .          shvi(nx_bg,ny_bg,nz_laps),   !Specific humidity (kg/kg)
     .          uwvi(nx_bg,ny_bg,nz_laps),   !U-wind (m/s)
     .          vwvi(nx_bg,ny_bg,nz_laps),   !V-wind (m/s)
     .          wwvi(nx_bg,ny_bg,nz_laps)    !W-wind (pa/s)
      integer msgpt(nx_bg,ny_bg)
c
c *** Background data interpolated to LAPS grid.
c
      real*4    ht(nx_laps,ny_laps,nz_laps), !Height (m)
     .          tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .          sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .          uw(nx_laps,ny_laps,nz_laps), !!U-wind (m/s)
     .          vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .          ww(nx_laps,ny_laps,nz_laps), !W-wind (pa/s)
     .          pr(nz_laps),                 !LAPS pressures
     .          lat(nx_laps,ny_laps),        !LAPS lat
     .          lon(nx_laps,ny_laps),        !LAPS lon
     .          topo(nx_laps,ny_laps),       !LAPS avg terrain
     .          grx(nx_laps,ny_laps),        !hinterp factor
     .          gry(nx_laps,ny_laps),         !hinterp factor
     .          ht_sfc(nx_laps,ny_laps),
     .          tp_sfc(nx_laps,ny_laps),
     .          Tdsfc(nx_laps,ny_laps),
     .          sh_sfc(nx_laps,ny_laps),
     .          qsfc(nx_laps,ny_laps),
     .          uw_sfc(nx_laps,ny_laps),
     .          vw_sfc(nx_laps,ny_laps),
     .          pr_sfc(nx_laps,ny_laps),
     .          mslp(nx_laps,ny_laps)
c
      real      ssh2,                        !Function name
     .          shsat,cti,
     .          htave,tpave,shave,uwave,vwave,
     .          std_lon,std_lat1,std_lat2
c
      integer   ct,
     .          ihour,imin,
     .          lga_files,lga_times(max_files),
     .          lga_valid,
     .          i4time_lga_valid(max_files),i4time_now,
     .          bg_times(max_files),bg_files,
     .          i4time_bg_valid(max_files),
     .          bg_valid(max_files),
     .          valid_bg(max_files),time_bg(max_files),
     .          bgvalid,
     .          i,ic,ii,j,jj,k,kk,l,ldl,lf,
     .          istatus,istatus_prep

      integer   i4lgatime(max_files)
      integer   i4bgtime
      integer   nlga
c
      character*255 lgapath
      character*256 lga_names(max_files)
      character*256 names(max_files)
      character*256 fname_bg(max_files)
      character*13  fname13,fname9_to_wfo_fname13
      character*6   c6_maproj
      character*2   gproj
      character*256 fullname,outdir
      character*31  ext
c     character*3   var(nz_laps)
c     character*4   lvl_coord(nz_laps)
c     character*10  units(nz_laps)
      character*125 comment(nz_laps)
      character*4   af_bg(max_files)
      character*10  c_domain_name
      character*200 c_dataroot
      integer len_dir,ntime, nf
c
      data ntime/0/
      data ext/'lga'/

      warncnt=0 
c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c
      print *,'in lga_sub'
      lga_status=0

      call get_directory('static',outdir,len_dir)    

      call find_domain_name(c_dataroot,c_domain_name,istatus)
      if(istatus.ne.1)then
          print*,'Error returned: find_domain_name'
          return
      endif
      call get_laps_domain(nx_laps,ny_laps,c_domain_name
     +,lat,lon,topo,istatus)
      if (istatus.lt.1)then
          print *,'Error reading lat, lon, topo data.'
          return
      endif
c
c *** Specify model path, extension for write laps routine.
c
      call get_directory('lga',outdir,len_dir)
      print *,'writing to dir ',outdir
c
c
c *** get LAPS pressure levels.  Using pressures.nl
c
      call get_pres_1d(i4time_now,nz_laps,pr,istatus)
      do k = 1,nz_laps
         pr(k)=pr(k)/100.
      enddo
c
c *** Determine which of the "bg_names" has not already been processed
c
      do j=1,max_files
         names(j)=bg_names(j)
      enddo
      do j=1,bg_files
         call i4time_fname_lp(names(j)(1:9),bg_times(j),istatus)
         read(bg_names(j)(12:13),'(i2)')ihour
         bg_valid(j)=ihour*3600
         i4time_bg_valid(j)=bg_times(j)+bg_valid(j)
      enddo

      call get_directory('lga',lgapath,ldl)
      call get_file_times(lgapath,max_files,lga_names
     1                      ,lga_times,nlga,istatus)
      if(nlga.gt.0)then
         k=1
         do while (k.le.nlga)

            read(lga_names(k)(ldl+10:ldl+11),'(i2)')ihour
            read(lga_names(k)(ldl+12:ldl+13),'(i2)')imin
            lga_valid=ihour*3600+imin*60

            do j=1,bg_files
               if(bg_times(j).eq.lga_times(k).and.
     1            bg_valid(j).eq.lga_valid)then
                  names(j)=' '
               endif
            enddo
            k=k+1
         enddo
      endif
      nbg=0
      do j=1,bg_files
         if(names(j).ne.' ')then
            nbg = nbg+1
            fname_bg(nbg)=bg_names(j)
            af_bg(nbg)=bg_names(j)(10:13)
            time_bg(nbg)=bg_times(j)
            valid_bg(nbg)=bg_valid(j)
         endif
      enddo

      if(nbg.eq.0)print*,'No new model background to process'
c
c ****** Read background model data.
c
      do nf=1,nbg
 
c Removal of this loop causes already existing lga files to be overwritten
c possibly with the same data.  However the error frequency on SBN may warrent
c this extra work.  
         if(.false.) then
            do i=1,lga_files
               if (fname_bg(nf) .eq. lga_names(i)(1:9) .and.
     .              af_bg(nf)(3:4) .eq. lga_names(i)(10:11)) then
                  
                  print *,i,lga_names(i),':',fname_bg(nf),':',af_bg(nf)
                  call get_lga_source(nx_laps,ny_laps,nz_laps
     +                 ,fname_bg(nf),af_bg(nf),comment(1))
                  call s_len(cmodel,ic)

                  if(cmodel(1:ic) .eq. comment(1)(1:ic)) then
                     print *,'LGA file exists, not making one. ' 
     +                    ,lga_names(i)
                     lga_status=1
                     goto 80
                  else
                     print *,'Overwritting LGA file ',lga_names(i)
     +                    ,'from model ',comment(1)(1:i),' with a new '
     +                    ,'one from model ',cmodel(1:ic)
                  endif
               endif
            enddo          
         endif     

         call s_len(bgpath,i)
         fullname = bgpath(1:i)//'/'//fname_bg(nf)
         call s_len(fullname,i)
c
c new subroutine to read appropriate background model file.
c
         print*
         print*,'Reading - ',fullname(1:i)
         print*
         call read_bgdata(nx_bg,ny_bg,nz_bg
     +    ,bgpath,fname_bg(nf),af_bg(nf),fullname,cmodel,bgmodel
     +    ,htbg, prbg,tpbg,uwbg,vwbg,shbg,wwbg
     +    ,htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc
     +    ,uwbg_sfc,vwbg_sfc,mslpbg
     +    ,gproj,lon0,lat1,lat2,istatus_prep)

         if (istatus_prep .ne. 0) then

            call s_len(fname_bg(nf),lf)
            call s_len(bgpath,l)
            if (bgmodel .gt. 1 .and. bgmodel .le. 3) then
               fname13=fname_bg(nf)(1:lf)//af_bg(nf)
            elseif (bgmodel .eq. 4) then
               fname13=fname9_to_wfo_fname13(fname_bg(nf))
            endif
            print *,'Error reading background model data for: ',
     .         bgpath(1:l)//'/'//fname13
            print *,'Process aborted for this file.'
            lga_status= -nf
            return
 
         endif
c
c ****** Vertically interpolate background data to LAPS isobaric levels.
c
         call vinterp(nx_bg,ny_bg,nz_bg,nz_laps,pr,
     .                prbg,htbg,tpbg,shbg,uwbg,vwbg,wwbg,
     .                htvi,tpvi,shvi,uwvi,vwvi,wwvi)
c
c ****** Run 2dx filter on vertically interpolated fields.
c ****** Only run filter on sh up to 300 mb, since the filter may
c           create small neg values when the field is small to begin with.
c ****** If bgmodel=4, there exists missing data.  Don't use these points
c           in the filter.
c
         if (bgmodel .eq. 4 .or. bgmodel .eq. 9) then
            do j=1,ny_bg
            do i=1,nx_bg
               if (htvi(i,j,1) .eq. missingflag) then
                  msgpt(i,j)=0
               else
                  msgpt(i,j)=1
               endif
            enddo
            enddo
c
            do k=1,nz_laps
            do j=2,ny_bg-1
            do i=2,nx_bg-1
               if (max(msgpt(i-1,j-1),
     .                 msgpt(i-1,j  ),
     .                 msgpt(i-1,j+1),
     .                 msgpt(i  ,j-1),
     .                 msgpt(i  ,j  ),
     .                 msgpt(i  ,j+1),
     .                 msgpt(i+1,j-1),
     .                 msgpt(i+1,j  ),
     .                 msgpt(i+1,j+1)) .eq. 1) then
                  ct=0
                  htave=0.
                  tpave=0.
                  shave=0.
                  uwave=0.
                  vwave=0.
                  do jj=j-1,j+1
                  do ii=i-1,i+1
                     if (msgpt(ii,jj) .eq. 1) then
                        ct=ct+1
                        htave=htave+htvi(ii,jj,k)
                        tpave=tpave+tpvi(ii,jj,k)
                        shave=shave+shvi(ii,jj,k)
                        uwave=uwave+uwvi(ii,jj,k)
                        vwave=vwave+vwvi(ii,jj,k)
                     endif
                  enddo
                  enddo
                  cti=1./float(ct)
                  htvi(i,j,k)=htave*cti
                  tpvi(i,j,k)=tpave*cti
                  shvi(i,j,k)=shave*cti
                  uwvi(i,j,k)=uwave*cti
                  vwvi(i,j,k)=vwave*cti
               endif
            enddo
            enddo
            enddo
         endif
c
         do kk=nz_laps,1,-1
            if (pr(kk) .ge. 300.) goto 20
         enddo
20       continue
         call filter_2dx(htvi,nx_bg,ny_bg,nz_laps, 0.5)
         call filter_2dx(htvi,nx_bg,ny_bg,nz_laps,-0.5)
         call filter_2dx(tpvi,nx_bg,ny_bg,nz_laps, 0.5)
         call filter_2dx(tpvi,nx_bg,ny_bg,nz_laps,-0.5)
         call filter_2dx(shvi,nx_bg,ny_bg,kk     , 0.5)
         call filter_2dx(shvi,nx_bg,ny_bg,kk     ,-0.5)
         call filter_2dx(uwvi,nx_bg,ny_bg,nz_laps, 0.5)
         call filter_2dx(uwvi,nx_bg,ny_bg,nz_laps,-0.5)
         call filter_2dx(vwvi,nx_bg,ny_bg,nz_laps, 0.5)
         call filter_2dx(vwvi,nx_bg,ny_bg,nz_laps,-0.5)
         call filter_2dx(wwvi,nx_bg,ny_bg,nz_laps, 0.5)
         call filter_2dx(wwvi,nx_bg,ny_bg,nz_laps,-0.5)
c
         if (bgmodel .eq. 4) then
            do j=1,ny_bg
            do i=1,nx_bg
               if (msgpt(i,j) .eq. 0) then
                  do k=1,nz_laps
                     htvi(i,j,k)=missingflag
                     tpvi(i,j,k)=missingflag
                     shvi(i,j,k)=missingflag
                     uwvi(i,j,k)=missingflag
                     vwvi(i,j,k)=missingflag
                  enddo
               endif
            enddo
            enddo
         endif
c
c ****** Horizontally interpolate background data to LAPS grid points.
c
         call init_hinterp(nx_bg,ny_bg,nx_laps,ny_laps,gproj,
     .        lat,lon,grx,gry,bgmodel)

         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,htvi,ht,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,shvi,sh,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,uwvi,uw,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,vwvi,vw,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,tpvi,tp,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,
     .        grx,gry,wwvi,ww,bgmodel)

c
c ****** Check for missing value flag in any of the fields.
c ****** Check for NaN's in any of the fields.
c
         do k=1,nz_laps
            do j=1,ny_laps
               do i=1,nx_laps
                  if((abs(ht(i,j,k)) .gt. 100000.) .or.
     +                 (tp(i,j,k) .gt. 1000.) .or.
     +                 (tp(i,j,k) .le. 0.) .or.
     +                 (abs(sh(i,j,k)) .gt. 1.) .or.
     +                 (abs(uw(i,j,k)) .gt. 150.) .or.
     +                 (abs(vw(i,j,k)) .gt. 150.))then
c
c ww may be missing from some models! Don't stop just because of that.
C    +             .or.(abs(ww(i,j,k)) .gt. 100.) )then
c                  if (max(ht(i,j,k),tp(i,j,k),sh(i,j,k),
c     .                 uw(i,j,k),vw(i,j,k)) .ge. missingflag) then

               print*,'ERROR: Missing or bad value detected: ',i,j,k
               print*,'ht/tp/sh/uw/vw/ww: ',ht(i,j,k),tp(i,j,k)
     +                  ,sh(i,j,k), uw(i,j,k),vw(i,j,k),ww(i,j,k)
                       lga_status = -nf
                       return
                  endif
c-----------------------------------
cc            if (ht(i,j,k) .ne. ht(i,j,k) .or. 
cc     .           tp(i,j,k) .ne. tp(i,j,k) .or.
cc     .           sh(i,j,k) .ne. sh(i,j,k) .or. 
cc     .           uw(i,j,k) .ne. uw(i,j,k) .or.
cc     .           vw(i,j,k) .ne. vw(i,j,k)) then
cc               print *,'ERROR: NaN detected:',i,j,k
cc               lga_status = -nf
cc               return
cc            endif
c-----------------------------------
               enddo
            enddo
         enddo
c
      call checknan_3d(ht,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array ht'
         lga_status = -nf
         return
      endif
c
      call checknan_3d(tp,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array tp'
         lga_status = -nf
         return
      endif
c
      call checknan_3d(sh,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array sh'
         lga_status = -nf
         return
      endif
c
      call checknan_3d(uw,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array uw'
         lga_status = -nf
         return
      endif
c
      call checknan_3d(vw,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array vw'
         lga_status = -nf
         return
      endif
      call checknan_3d(ww,nx_laps,ny_laps,nz_laps,nan_flag)
      if(nan_flag .ne. 1) then
         print *,' ERROR: NaN found in array vw'
         lga_status = -nf
         return
      endif
c
c ****** Horizontally interpolate background surface data to LAPS grid points.
c
      if(bgmodel.ne.1.and.bgmodel.ne.9)then

         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,htbg_sfc,ht_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,tpbg_sfc,tp_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,shbg_sfc,sh_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,uwbg_sfc,uw_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,vwbg_sfc,vw_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,prbg_sfc,pr_sfc,bgmodel)
         call hinterp_field(nx_bg,ny_bg,nx_laps,ny_laps,1,
     .        grx,gry,mslpbg,mslp,bgmodel)
c
c check for T > Td before sfc p computation. Due to large scale
c interpolation we can have slightly larger (fractional) Td than T.
c
         call tdcheck(nx_laps,ny_laps,sh_sfc,tp_sfc,
     &icnt,i_mx,j_mx,diff_mx,diff_mn)

         print *,' Dewpoint check (before call sfcbkgd):'
         print *,'     Dewpt greater than temp at ',icnt,' points.'
         if(icnt .gt. 0) then

            print*,'Max diff of ',diff_mx,' at ',i_mx,',',j_mx
            print*,'Min diff of ',diff_mn,' at ',i_mn,',',j_mn

         endif
      endif
c
c... Do the temp, moisture (sh_sfc returns with Td), and pressures
c
      call sfcbkgd(bgmodel, tp, sh, ht, tp_sfc, sh_sfc, topo, pr,
     .            nx_laps, ny_laps, nz_laps, pr_sfc)

      call tdcheck(nx_laps,ny_laps,sh_sfc,tp_sfc,
     &icnt,i_mx,j_mx,diff_mx,diff_mn)

      print *,' Dewpoint check (after call sfcbkgd):'
      print *,'     Dewpt greater than temp at ',icnt,' points.'
      if(icnt .gt. 0) then

         print*,'Max diff of ',diff_mx,' at ',i_mx,',',j_mx
         print*,'Min diff of ',diff_mn,' at ',i_mn,',',j_mn

      endif
c
c..... Do the winds
c
      call interp_to_sfc(topo,uw,ht,nx_laps,ny_laps,
     &                         nz_laps,missingflag,uw_sfc)
      call interp_to_sfc(topo,vw,ht,nx_laps,ny_laps,
     &                         nz_laps,missingflag,vw_sfc)
c
c ****** Eliminate any supersaturations or negative sh generated 
c           through interpolation (set min sh to 1.e-6).
c
      do k=1,nz_laps
         do j=1,ny_laps
         do i=1,nx_laps
            shsat=ssh2(pr(k),tp(i,j,k)-273.15,
     .             tp(i,j,k)-273.15,-47.0)*0.001
            sh(i,j,k)=max(1.0e-6,min(sh(i,j,k),shsat))
         enddo
         enddo
      enddo
c
      do j=1,ny_laps
      do i=1,nx_laps
         if(sh_sfc(i,j).gt.tp_sfc(i,j))then
            sh_sfc(i,j)=tp_sfc(i,j)
         endif
      enddo
      enddo
c
c the wind components are still on the native grid projection;
c rotate them to the LAPS (output) domain as necessary.

      print*

      if(.true.)then
         call rotate_background_uv(nx_laps,ny_laps,nz_laps,lon
     +,gproj,lon0,lat1,lat2,uw,vw,uw_sfc,vw_sfc,istatus)
      else
         print*,'Not rotating u/v'
      endif

      if(istatus.ne.1)then
         print*,'Error in rotate_background_uv '
         return
      endif

c
c Write LGA
c ---------
      bgvalid=time_bg(nf)+valid_bg(nf)
      call write_lga(nx_laps,ny_laps,nz_laps,time_bg(nf),
     .bgvalid,cmodel,missingflag,pr,ht,tp,sh,uw,vw,ww,istatus)
      if(istatus.ne.1)then
         print*,'Error writing lga - returning to main'
         return
      endif
c         
c Write LGB
c ---------

      if(bgmodel.ne.7)then
 
         do j=1,ny_laps
         do i=1,nx_laps
            if(pr_sfc(i,j) .lt. missingflag) then
               qsfc(i,j)=ssh2(pr_sfc(i,j)*0.01,
     +                   tp_sfc(i,j)-273.15,
     +                   sh_sfc(i,j)-273.15,-47.)*0.001
c              sfcgrid(i,j,kk+4)=qsfc(i,j)
            else
               qsfc(i,j) = missingflag
c              sfcgrid(i,j,kk+4)=missingflag
            endif
         enddo
         enddo

         call write_lgb(nx_laps,ny_laps,time_bg(nf),bgvalid
     .,cmodel,missingflag,uw_sfc,vw_sfc,tp_sfc,qsfc,pr_sfc,mslp
     .,sh_sfc,istatus)
         if(istatus.ne.1)then
            print*,'Error writing lgb - returning to main'
            return
         endif

      endif

      lga_status = 1
c
80    continue

      enddo  !Main loop through two model backgrounds

c time interpolate between existing lga (bg) files.
c-------------------------------------------------------------------------------
c
      do i=1,nbg
         call s_len(bg_names(i),j)
         print*,'bg_name: ',i,bg_names(i)(1:j)
      enddo
c
c *** Determine if new file needs to (can) be created and perform
c        linear time interpolation.
      if(bg_files.gt.1)then
         i=bg_files
         print*,i,bg_files,bg_times(i),bg_times(i-1),
     +     bg_valid(i),bg_valid(i-1),laps_cycle_time,
     +     i4time_bg_valid(i),i4time_bg_valid(i-1)

         if(i4time_bg_valid(i-1) .gt.i4time_now   .and.
     +      i4time_bg_valid(i)   .lt.i4time_now  )then
            ext = 'lga'
            call get_directory(ext,outdir,len_dir) 
            print*,outdir,ext,nz_laps

c interp 3D fields
            call time_interp(outdir,ext,
     +           nx_laps,ny_laps,nz_laps,6,pr,
     +           i4time_bg_valid(i),i4time_bg_valid(i-1),
     +           i4time_now,bg_times(i-1),bg_valid(i-1),
     +           bg_times(i  ),bg_valid(i  ))

            if(bgmodel.ne.1.or.bgmodel.ne.7)then
               ext = 'lgb'
               call get_directory(ext,outdir,len_dir) 
               print*,outdir,ext

c interp 2D fields
               call time_interp(outdir,ext,
     +           nx_laps,ny_laps,1,7,pr(1),
     +           i4time_bg_valid(i),i4time_bg_valid(i-1),
     +           i4time_now,bg_times(i-1),bg_valid(i-1),
     +           bg_times(i  ),bg_valid(i  ))
            endif
            lga_status = 1
         else
            print*,'Time Interpolation Not Necessary!'
            lga_status = 1
         endif

      else
         print*,'No time interp when bg_files < 2'
         print*
         lga_status=1
      endif

      return
      end

      subroutine get_lga_source(nx,ny,nz,fname,af,source)
      character*9 fname
      character*4 af
      character*10 dumb1
      character*125 comment
      character*12 source
      integer ihour, bgtime, nx,ny,nz
      real dumb2(nx,ny,nz)

      call i4time_fname_lp(fname,bgtime,istatus)
      read(af,'(i4)') ihour
      call get_lapsdata_3d(bgtime,bgtime+ihour*3600,nx,ny,
     +        nz,'lga ','HT ',dumb1,comment,dumb2,istatus)
      if(istatus.lt.1) then
         stop 'error returned from get_lapsdata_3d'
      endif


      source = comment(1:12)
      return
      end
c
c
      subroutine checknan_2d(x,ni,nj,nan_flag)
c
c     Routine to check a real 2-d array for NaN's.
c
      integer ni,nj
      real*4 x(ni,nj)
c
      nan_flag = 1
c
      do j=1,nj
      do i=1,ni
         if( nan( x(i,j) ) .eq. 1) then
            print *,' ** ERROR. Found a NaN at ', i, j
            nan_flag = -1
            return
         endif
      enddo !i
      enddo !j
c
      return
      end
c
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
      dmax = -1.e30
      dmin =  1.e30
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
            if(diff .gt. dmax) then
               dmax = diff
               i_mx = i
               j_mx = j
            endif
            if(diff .lt. dmin) then
               dmin = diff
               i_mn = i
               j_mn = j
            endif
         endif
      enddo
      enddo

      return
      end
