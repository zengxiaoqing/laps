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
c        bgmodel = 3 ---> NOGAPS (2.5 deg)
c        bgmodel = 4 ---> SBN Conus211 (Eta or RUC)
c        bgmodel = 5 ---> RUC (40 km native grid)
c        bgmodel = 6 ---> AVN (360 x 181 lat-lon grid)
c        bgmodel = 7 ---> ETA (48 km from grib file)
c        bgmodel = 8 ---> NOGAPS (1.0 deg)
c        bgmodel = 9 ---> NWS Conus (RUC, ETA, NGM, AVN)
c
      integer nbgmodel
      parameter (nbgmodel=9)
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
     .          lga_status                      !status returned from lga_driver
c                                            ! 1=good 0=bad
      real*4    prbot,delpr                  !LAPS bottom and delta pressures
c
      character*255 lapsroot                 !LAPS root path
      character*31  laps_domain_file         !LAPS domain desinator
c
c------------------> BACKGROUND GRID DATA PATH <--------------------------------
c
c *** The following variables specify the various model data directory paths.
c *** Read from standard input.
c *** Each line contains the background data directory path for each
c        corresponding bgmodel defined above.
c
      character*150 bgpath

c
c  This is the max number of paths allowed in nest7grid.parms
c  and should match the value in lib/lapsgrid.f 
c 
      integer maxbgmodels
      parameter (maxbgmodels=10)
      character*150 bgpaths(maxbgmodels)
      integer bgmodels(maxbgmodels), bglen
c
      character*9 a9
      integer i4time_now
      integer*4 max_files,bg_files
      parameter (max_files=200)
      character*100 names(max_files)
      character*100 reject_files(max_files)
      integer reject_cnt
      data reject_cnt/0/

c
c-------------------------------------------------------------------------------
c
c *** Comments used in writing netcdf files only.
c
      integer istat, i
c
c cmodel is really only 12 chars but the SBN netcdf carrys 132
c
      character*132 cmodel
      include 'lapsparms.cmn'
      integer oldest_forecast, max_forecast_delta
      logical use_analysis
c_______________________________________________________________________________
c
c Read background model from nest7grid.parms
      laps_domain_file = 'nest7grid'
c      istat = index(laps_domain_file,' ')-1

      call s_len(laps_domain_file,istat)
      call get_laps_config(laps_domain_file(1:istat),istat)
      nx_laps = NX_L_CMN
      ny_laps = NY_L_CMN
      nz_laps = nk_laps
      prbot = PRESSURE_BOTTOM_L/100.
      delpr = PRESSURE_INTERVAL_L/100.
      laps_cycle_time = laps_cycle_time_cmn

      i=1

      call get_background_info(150,bgpaths,bgmodels
     +     ,oldest_forecast,max_forecast_delta,use_analysis) 
      lga_status = 0

c
c *** Initialize esat table.
c
      call es_ini
      bg_files=0
      i=0
      do while(lga_status.le.0 .and. i.le.nbgmodel)
         if(bg_files.le.reject_cnt) then
            i=i+1 
            bgmodel = bgmodels(i)
            if(bgmodel .eq. 0) goto 965
            call s_len(bgpaths(i),bglen)
            bgpath =  bgpaths(i)
         endif
         if (bgmodel .lt. 1 .or. bgmodel .gt. nbgmodel) then
            print*,'Bad model specification in LGA, bgmodel =',bgmodel
            print*,'   LGA process aborted...'
            stop
         endif

c *** Get current time from systime.dat
c
         call get_systime(i4time_now,a9,lga_status)
         call get_acceptable_files(i4time_now,bgpath,bgmodel
     +        ,names,max_files,oldest_forecast,max_forecast_delta
     +        ,use_analysis,bg_files,0,cmodel,nx_bg,ny_bg,nz_bg
     +        ,reject_files,reject_cnt)

        if(bg_files.le.1) then
           print*,'No Acceptable files found for model: ',bgpaths(i),
     +          bgmodels(i) 
           lga_status = 0
        else


           print *, 'bgmodel is set to ', bgmodel

           print *, ' '
           print *, 'input parameters'
           print *,  ' '
           print *, nx_laps,ny_laps,nz_laps,prbot,delpr
           print *, laps_cycle_time
ccc   print *, lapsroot
           print *, laps_domain_file
           print *, nx_bg,ny_bg,nz_bg
           print *, 'bgpath ', bgpath(1:bglen)
           print *, 'cmodel ',cmodel
 970       continue
c
c *** Call lga driver.
c

           call lga_driver(nx_laps,ny_laps,nz_laps,prbot,delpr,
     .          laps_cycle_time,lapsroot,laps_domain_file,
     .          bgmodel,bgpath,names,cmodel,
     .          nx_bg,ny_bg,nz_bg, 5*nz_laps, lga_status)

           if(lga_status.lt.0) then
              reject_cnt=reject_cnt+1
              reject_files(reject_cnt)=names(-lga_status)
           endif



        endif
      enddo

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
      subroutine lga_driver(nx_laps,ny_laps,nz_laps,prbot,delpr,
     .     laps_cycle_time,lapsroot,laps_domain_file,
     .     bgmodel,bgpath,bg_names,cmodel,
     .     nx_bg,ny_bg,nz_bg, kdim, lga_status)

c
      implicit none
      include 'bgdata.inc'
c
      integer nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .          nx_bg  ,ny_bg  ,nz_bg,       !Background model grid dimensions
     .          max_files, lga_status, laps_cycle_time,
     .          bgmodel
      real prbot, delpr
      character*(*) lapsroot, laps_domain_file, bgpath, bg_names(2)
     +     , cmodel
      
c
      parameter (max_files=2000)
      
c
c *** Background model grid data.
c
      real*4    prbg(nx_bg,ny_bg,nz_bg),     !Pressure (mb)
     .          htbg(nx_bg,ny_bg,nz_bg),     !Height (m)
     .          tpbg(nx_bg,ny_bg,nz_bg),     !Temperature (K)
     .          shbg(nx_bg,ny_bg,nz_bg),     !Specific humidity (kg/kg)
     .          uwbg(nx_bg,ny_bg,nz_bg),     !U-wind (m/s)
     .          vwbg(nx_bg,ny_bg,nz_bg),      !V-wind (m/s)
     .          mslpbg(nx_bg,ny_bg),          !mslp  (mb)
     .          wwbg(nx_bg,ny_bg,nz_bg),      !W-wind (m/s)
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
     .          vwvi(nx_bg,ny_bg,nz_laps)    !V-wind (m/s)
      integer msgpt(nx_bg,ny_bg)
c
c *** Background data interpolated to LAPS grid.
c
      integer kdim
      real*4    ht(nx_laps,ny_laps,nz_laps), !Height (m)
     .          tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .          sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .          uw(nx_laps,ny_laps,nz_laps), !!U-wind (m/s)
     .          vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .          grid(nx_laps,ny_laps,kdim), !Full LAPS array for write_laps
     .          pr(nz_laps),     !LAPS pressures
     .          lat(nx_laps,ny_laps),        !LAPS lat
     .          lon(nx_laps,ny_laps),         !LAPS lon
     .          grx(nx_laps,ny_laps),         !hinterp factor
     .          gry(nx_laps,ny_laps),         !hinterp factor
     .          ht_sfc(nx_laps,ny_laps),
     .          tp_sfc(nx_laps,ny_laps),
     .          sh_sfc(nx_laps,ny_laps),
     .          uw_sfc(nx_laps,ny_laps),
     .          vw_sfc(nx_laps,ny_laps),
     .          pr_sfc(nx_laps,ny_laps),
     .          mslp(nx_laps,ny_laps)
c
      real      ssh2,                        !Function name
     .          shsat,cti,
     .          htave,tpave,shave,uwave,vwave
c
      integer   ct,
     .          ihour,imin,
     .          lga_files,lga_times(max_files),lga_valid(max_files),
     .          bgtime,bgvalid,
     .          ip(kdim),
     .          i,ic,ii,j,jj,k,kk,l,
     .          istatus
c
      character*255 lgapath
      character*100 lga_names(max_files)
      character*13  fname13,fname9_to_wfo_fname13
      character*9   fname
      character*2   gproj
      character*150  outdir, fullname
      character*31  ext
      character*3   var(kdim)
      character*4   af
      character*4   lvl_coord(kdim)
      character*10  units(kdim)
      character*125 comment(kdim)
      integer len_dir, ntime, nf
      integer nxbg, nybg, nzbg(5),ntbg

c

      data ntime/0/
      data ext/'lga'/
c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c


      lga_status=0
      call get_directory('static',outdir,len_dir)    

      call get_laps_lat_lon(outdir(1:len_dir),laps_domain_file,
     .                      nx_laps,ny_laps,lat,lon,istatus)
      if (istatus .ne. 1) print *,'Error reading LAPS lat, lon data.'
c
c *** Specify model path, extension for write laps routine.
c
      call get_directory('lga',outdir,len_dir)
      print *,'writing to dir ',outdir
c
c
c *** Compute LAPS pressure levels.
c
      do k=1,nz_laps
         pr(k)=prbot-delpr*float(k-1)
      enddo
c
c *** Get .lga file names.
c
      lgapath=outdir(1:len_dir)//'/*.lga'
      call get_file_names(lgapath,lga_files,lga_names,max_files,istatus)
      do i=1,lga_files
c         j=index(lga_names(i),' ')-18
         call s_len(lga_names(i),j)
         j=j-17
 
         if (j .ge. 0) lga_names(i)=lga_names(i)(j+1:j+13)
ccc         print *,i,lga_names(i)
      enddo
c

c
c ****** Read background model data.
c
      do nf=1,2
         call s_len(bg_names(nf),j)
         j=j-13
         fname = bg_names(nf)(j+1:j+9)
         af = bg_names(nf)(j+10:j+13)
         call i4time_fname_lp(fname,bgtime,istatus)

         do i=1,lga_files
            if (fname .eq. lga_names(i)(1:9) .and.
     .          af(3:4) .eq. lga_names(i)(10:11)) then

               print *,i,lga_names(i),':',fname,':',af
               call get_lga_source(nx_laps,ny_laps,nz_laps
     +              ,fname,af,comment(1))
               call s_len(cmodel,ic)
               if(cmodel(1:ic) .eq. comment(1)(1:ic)) then
                  print *,'LGA file exists, not making one. ' 
     +                   ,lga_names(i)
                  lga_status=1
                  goto 80
               else
                  print *,'Overwritting LGA file ',lga_names(i)
     +            ,'from model ',comment(1)(1:i),' with a new one'
     +            ,' from model ',cmodel(1:ic)
               endif
            endif
         enddo               


         call s_len(bgpath,i)
         fullname = bgpath(1:i)//'/'//bg_names(nf)
         if (bgmodel .eq. 1) then     ! Process 60 km RUC data
            call read_ruc60_native(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                             prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                             gproj,istatus)
 
         elseif (bgmodel .eq. 2) then ! Process 48 km ETA conus-c grid data
            call read_eta_conusC(fullname,nx_bg,ny_bg,nz_bg,
     .                          htbg, prbg,tpbg,uwbg,vwbg,shbg,
     .                          htbg_sfc,prbg_sfc,shbg_sfc,tpbg_sfc,
     .                          uwbg_sfc,vwbg_sfc,mslpbg,
     .                          istatus)
            if(istatus.gt.0) then
               call lprep_eta_conusc(nx_bg,ny_bg,nz_bg,htbg
     +              ,prbg,tpbg,uwbg,vwbg,shbg,gproj,istatus)
            endif
         elseif (bgmodel .eq. 3 .or. 
     .           bgmodel .eq. 8) then ! Process NOGAPS data
            call read_nogaps(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                       prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                       gproj,istatus)
 
         elseif (bgmodel .eq. 4) then ! Process SBN Conus 211 data (Eta or RUC)
            ntbg=10 
            
            call get_sbn_dims(bgpath,fname,nxbg,nybg,nzbg,ntbg)

            call read_conus_211(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .           nxbg,nybg,nzbg,ntbg,
     .           prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .           prbg_sfc,uwbg_sfc,vwbg_sfc,shbg_sfc,tpbg_sfc,
     .           mslpbg,gproj,1,istatus)
c
         elseif (bgmodel .eq. 5) then ! Process 40 km RUC data
            call read_ruc2_hybb(fullname,nx_bg,ny_bg,nz_bg
     +                     ,mslpbg,htbg,prbg,shbg,uwbg,vwbg,tpbg,wwbg
     +                    ,istatus)
            if(istatus.gt.0) then

               print*,'Read complete: entering prep'
               call lprep_ruc2_pub(nx_bg,ny_bg,nz_bg
     +              ,htbg,prbg,shbg,uwbg,vwbg,tpbg,gproj)

               print*,'Data prep complete'
            endif

c     
         elseif (bgmodel .eq. 6 .or.
     .           bgmodel .eq. 7) then ! Process AVN or ETA grib data
            call read_dgprep(bgmodel,bgpath,fname,af,nx_bg,ny_bg,nz_bg
     .                      ,prbg,htbg,tpbg,shbg,uwbg,vwbg
     .                      ,gproj,istatus)
         elseif (bgmodel .eq. 9) then ! Process NWS Conus data (RUC,ETA,NGM,AVN)
            call read_conus_nws(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                          prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                          gproj,istatus)
c
         endif
         
         if (istatus .ne. 1) then
c            l=index(bgpath,' ')-1

            call s_len(bgpath,l)
            if (bgmodel .gt. 1 .and. bgmodel .lt. 3) then
               fname13=fname//af
            elseif (bgmodel .eq. 4) then
               fname13=fname9_to_wfo_fname13(fname)
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
     .                prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                htvi,tpvi,shvi,uwvi,vwvi)
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


c
c ****** Horizontally interpolate background surface data to LAPS grid points.
c
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
     +                 (abs(uw(i,j,k)) .gt. 100.) .or.
     +                 (abs(vw(i,j,k)) .gt. 100.) ) then
c                  if (max(ht(i,j,k),tp(i,j,k),sh(i,j,k),
c     .                 uw(i,j,k),vw(i,j,k)) .ge. missingflag) then
               print*,'ERROR: Missing or bad value detected: ',i,j,k
               print*,ht(i,j,k),tp(i,j,k),sh(i,j,k), uw(i,j,k),vw(i,j,k)
               lga_status = -nf
               return
            endif
            if (ht(i,j,k) .ne. ht(i,j,k) .or. 
     .           tp(i,j,k) .ne. tp(i,j,k) .or.
     .           sh(i,j,k) .ne. sh(i,j,k) .or. 
     .           uw(i,j,k) .ne. uw(i,j,k) .or.
     .           vw(i,j,k) .ne. vw(i,j,k)) then
               print *,'ERROR: NaN detected:',i,j,k
               lga_status = -nf
               return
            endif
         enddo
      enddo
      enddo
c
c ****** Eliminate any supersaturations or negative sh generated 
c           through interpolation (set min sh to 1.e-6).
c
         do k=1,nz_laps
         do j=1,ny_laps
         do i=1,nx_laps
            shsat=ssh2(pr(k),tp(i,j,k)-273.15,
     .             tp(i,j,k)-273.15,0.0)*0.001
            sh(i,j,k)=max(1.0e-6,min(sh(i,j,k),shsat))
         enddo
         enddo
         enddo
c
c ****** Fill grid for LAPS write routine.
c
c         i=index(cmodel,' ')-1

         call s_len(cmodel,ic)
         do k=1,nz_laps                 ! Height
            do j=1,ny_laps
            do i=1,nx_laps
               grid(i,j,k)=ht(i,j,k)
            enddo
            enddo
            ip(k)=int(pr(k))
            var(k)='HT '
            lvl_coord(k)='mb  '
            units(k)='Meters'
            comment(k)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         enddo
c
         do k=1,nz_laps                 ! Temperature
            kk=k+nz_laps
            do j=1,ny_laps
            do i=1,nx_laps
               grid(i,j,kk)=tp(i,j,k)
            enddo
            enddo
            ip(kk)=int(pr(k))
            var(kk)='T3 '
            lvl_coord(kk)='mb  '
            units(kk)='Kelvin'
            comment(kk)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         enddo
c
         do k=1,nz_laps                 ! Specific humidity
            kk=k+2*nz_laps
            do j=1,ny_laps
            do i=1,nx_laps
               grid(i,j,kk)=sh(i,j,k)
            enddo
            enddo
            ip(kk)=int(pr(k))
            var(kk)='SH '
            lvl_coord(kk)='mb  '
            units(kk)='kg/kg'
            comment(kk)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         enddo
c
         do k=1,nz_laps                 ! u-component wind
            kk=k+3*nz_laps
            do j=1,ny_laps
            do i=1,nx_laps
               grid(i,j,kk)=uw(i,j,k)
            enddo
            enddo
            ip(kk)=int(pr(k))
            var(kk)='U3 '
            lvl_coord(kk)='mb  '
            units(kk)='m/s'
            comment(kk)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         enddo
c
         do k=1,nz_laps                 ! v-component wind
            kk=k+4*nz_laps
            do j=1,ny_laps
            do i=1,nx_laps
               grid(i,j,kk)=vw(i,j,k)
            enddo
            enddo
            ip(kk)=int(pr(k))
            var(kk)='V3 '
            lvl_coord(kk)='mb  '
            units(kk)='m/s'
            comment(kk)=cmodel(1:ic)//' interpolated to LAPS isobaric.'
         enddo


c
         read(af,'(i4)') ihour
         bgvalid=bgtime+ihour*3600
c
         ext = 'lga'
         call get_directory(ext,outdir,len_dir) 
         print *,'Writing - ',fname//af(3:4),'00.',ext

         call write_laps(bgtime,bgvalid,outdir,ext,
     .                   nx_laps,ny_laps,nz_laps,kdim,var,
     .                   ip,lvl_coord,units,comment,grid,istatus)
         

         if (istatus .eq. 1) then
            lga_status = 1
         else
           print *,'Error writing interpolated data to LAPS database.'
         endif

         if(bgmodel.eq.2.or.bgmodel.eq.4) then

c
c Write the 2d fields to lgb
c
            kk=0
            var(kk+1)='USF'
            units(kk+1)='m/s'
            var(kk+2)='VSF'
            units(kk+2)='m/s'
            var(kk+3)='TSF'
            units(kk+3)='K'
            var(kk+4)='RSF'
            units(kk+4)='none'
            var(kk+5)='PSF'
            units(kk+5)='PA'
            var(kk+6)='SLP'
            units(kk+6)='PA'
            do j=1,ny_laps
               do i=1,nx_laps
                  grid(i,j,kk+1) = uw_sfc(i,j)
                  grid(i,j,kk+2) = vw_sfc(i,j)
                  grid(i,j,kk+3) = tp_sfc(i,j)
                  grid(i,j,kk+4) = sh_sfc(i,j)
                  grid(i,j,kk+5) = pr_sfc(i,j)
                  grid(i,j,kk+6) = mslp(i,j)

               enddo
            enddo
            do kk=1,6
               ip(kk)=0
            enddo
         

            ext = 'lgb'
            call get_directory(ext,outdir,len_dir)
            print *,'writing to dir ',outdir
            print *,'Writing - ',fname//af(3:4),'00.',ext
            call write_laps(bgtime,bgvalid,outdir,ext,
     .           nx_laps,ny_laps,6,6,var,
     .           ip,lvl_coord,units,comment,grid,istatus)
         

            if (istatus .eq. 1) then
               lga_status = 1
            else
               print*,'Error writing interpolated data to LAPS lgb'
            endif
         endif

c
 80      continue
      enddo

      if(lga_status.le.0) return

c
c-------------------------------------------------------------------------------
c *** Backfill any missing file times using linear interpolation
c        between existing files.
c-------------------------------------------------------------------------------
c
c *** Get existing .lga file names.
c
      call get_file_names(lgapath,lga_files,lga_names,max_files,istatus)
c
c *** Convert file names to i4times and sort into one array.
c
      do i=1,lga_files
c         j=index(lga_names(i),' ')-18

         call s_len(lga_names(i),j)
         j=j-17
         if (j .ge. 0) then
            call i4time_fname_lp(lga_names(i)(j+1:j+9),
     .                           lga_times(i),istatus)
            read(lga_names(i)(j+10:j+11),'(i2)') ihour
            read(lga_names(i)(j+12:j+13),'(i2)') imin
            lga_valid(i)=ihour*3600+imin*60
         endif
      enddo
c
c *** Determine if new file needs to be created and perform
c        linear interpolation.
c
      do i=lga_files,2,-1
         if (lga_times(i) .eq. lga_times(i-1) .and.
     .       lga_valid(i)-lga_valid(i-1) .gt. laps_cycle_time) then
            ext = 'lga'
            call get_directory(ext,outdir,len_dir) 
            print*,outdir,ext,nz_laps
            call time_interp(outdir,ext,
     .           nx_laps,ny_laps,nz_laps,nz_laps*5,
     .           pr,laps_cycle_time,
     .           lga_times(i-1),lga_valid(i-1),
     .           lga_times(i  ),lga_valid(i  ))
            if(bgmodel.eq.2.or.bgmodel.eq.4) then
               ext = 'lgb'
               call get_directory(ext,outdir,len_dir) 
               print*,outdir,ext
               call time_interp(outdir,ext,
     .              nx_laps,ny_laps,1,6,
     .              pr,laps_cycle_time,
     .              lga_times(i-1),lga_valid(i-1),
     .              lga_times(i  ),lga_valid(i  ))
            endif
         endif
      enddo
c
      return
c
c *** Error traps.
c
900   print *,'Could not find systime.dat.'
      stop
c
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
      if(istatus.ne.1) then
         stop 'error returned from get_lapsdata_3d'
      endif


      source = comment(1:12)
      return
      end
