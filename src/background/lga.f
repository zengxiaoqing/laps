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
      parameter (maxbgmodels=4)
      character*150 bgpaths(maxbgmodels)
      integer bgmodels(maxbgmodels), len
c
      integer*4 max_files,bg_files
      parameter (max_files=2000)
      character*100 names(max_files)
c
c-------------------------------------------------------------------------------
c
c *** Comments used in writing netcdf files only.
c
      integer istat, i
      character*12 cmodel
      include 'lapsparms.cmn'
      integer oldest_forecast
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

      call get_background_info(150,bgpaths,bgmodels,oldest_forecast
     +                      ,use_analysis) 
      lga_status = 0
      do while(lga_status.eq.0 .and. i.le.nbgmodel)
         bgmodel = bgmodels(i)
         if(bgmodel .eq. 0) goto 965
c         len = index(bgpaths(i),' ')

         call s_len(bgpaths(i),len)
         bgpath =  bgpaths(i)
          
         if(bgmodel.eq.1) then
            nx_bg = 81
            ny_bg = 62
            nz_bg = 25        
            cmodel = 'RUC60_NATIVE'   
         else if(bgmodel.eq.2) then
            nx_bg = 93
            ny_bg = 65
            nz_bg = 19        
            cmodel = 'ETA48_CONUS'        
         else if(bgmodel.eq.3) then
            nx_bg = 144
            ny_bg = 73
            nz_bg = 16        
            cmodel = 'NOGAPS (2.5)'            
         else if(bgmodel.eq.4) then
            nx_bg = 93
            ny_bg = 65
            nz_bg = 19        
            cmodel = 'SBN CONUS211'            
         else if(bgmodel.eq.5) then
            nx_bg = 151
            ny_bg = 113
            nz_bg = 40       
            cmodel = 'RUC40_NATIVE'            
         else if(bgmodel.eq.6) then
            nx_bg = 360
            ny_bg = 181
            nz_bg = 16
            cmodel = 'AVN_LL_GRIB'
         else if(bgmodel.eq.7) then
            nx_bg = 185
            ny_bg = 129
            nz_bg = 39
            cmodel = 'ETA48_GRIB'
         else if(bgmodel.eq.8) then
            nx_bg = 360
            ny_bg = 181
            nz_bg = 16        
            cmodel = 'NOGAPS (1.0)'            
         else if(bgmodel.eq.9) then
            call get_file_names(bgpath,bg_files,names,max_files,istat)
            if (istat .ne. 1) print *,'Error in get_file_names.'
            call get_conus_dims(names(1),nx_bg,ny_bg,nz_bg)
            print *,'nws:',nx_bg,ny_bg,nz_bg
            cmodel='NWS_CONUS'
         endif
         

         if (bgmodel .lt. 1 .or. bgmodel .gt. nbgmodel) then
            print *,'Bad model specification in LGA, bgmodel =',bgmodel
            print *,'   LGA process aborted...'
            stop
         endif

         print *, 'bgmodel is set to ', bgmodel

         print *, ' '
         print *, 'input parameters'
         print *,  ' '
         print *, nx_laps,ny_laps,nz_laps,prbot,delpr
         print *, laps_cycle_time
ccc         print *, lapsroot
         print *, laps_domain_file
         print *, nx_bg,ny_bg,nz_bg

c
c *** Read background model data directory paths from standard input.
c
         print *, 'bgpath ', bgpath(1:len)
c
c *** Read background model name for netcdf comment output only.
c
         print *, 'cmodel ',cmodel
 970     continue
c
c *** Initialize esat table.
c
         call esat_init
c
c *** Call lga driver.
c

         call lga_driver(nx_laps,ny_laps,nz_laps,prbot,delpr,
     .        laps_cycle_time,lapsroot,laps_domain_file,
     .        bgmodel,bgpath,oldest_forecast,use_analysis,cmodel,
     .        nx_bg,ny_bg,nz_bg, lga_status)
         i=i+1 
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
     .       laps_cycle_time,lapsroot,laps_domain_file,
     .       bgmodel,bgpath,oldest_forecast,use_analysis,cmodel,
     .       nx_bg,ny_bg,nz_bg, lga_status)

c
      implicit none
      include 'netcdf.inc'
c
      integer nx_laps,ny_laps,nz_laps,     !LAPS grid dimensions
     .          nx_bg  ,ny_bg  ,nz_bg,       !Background model grid dimensions
     .          max_files, lga_status, ncid
      integer oldest_forecast
      logical use_analysis
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
     .          vwbg(nx_bg,ny_bg,nz_bg)      !V-wind (m/s)
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
      real*4    ht(nx_laps,ny_laps,nz_laps), !Height (m)
     .          tp(nx_laps,ny_laps,nz_laps), !Temperature (K)
     .          sh(nx_laps,ny_laps,nz_laps), !Specific humidity (kg/kg)
     .          uw(nx_laps,ny_laps,nz_laps), !!U-wind (m/s)
     .          vw(nx_laps,ny_laps,nz_laps), !V-wind (m/s)
     .          grid(nx_laps,ny_laps,nz_laps*5), !Full LAPS array for write_laps
c
     .          pr(nz_laps),prbot,delpr,     !LAPS pressures
     .          lat(nx_laps,ny_laps),        !LAPS lat
     .          lon(nx_laps,ny_laps)         !LAPS lon
c
      real      ssh2,                        !Function name
     .          shsat,cti,
     .          htave,tpave,shave,uwave,vwave,
     .          msgflg                       !Missing data value
c
      integer bgmodel,laps_cycle_time,
     .          i4time_now,ct4,ct,
     .          ihour,imin,
     .          lga_files,lga_times(max_files),lga_valid(max_files),
     .          bg_files,bgtime,bgvalid,
     .          ip(5*nz_laps),
     .          i,ii,j,jj,k,kk,l,ll,n,
     .          istatus
c
      character*(*) bgpath,lapsroot,laps_domain_file
      character*255 lgapath
      character*100 bg_names(max_files)
      character*100 lga_names(max_files)
      character*100 names(max_files)
      character*13  fname13,fname9_to_wfo_fname13,file_list(2)
      character*9   a9,fname,wfo_fname13_to_fname9
      character*2   gproj
      character*150  outdir
      character*31  ext
      character*3   var(5*nz_laps)
      character*4   af
      character*4   lvl_coord(5*nz_laps)
      character*10  units(5*nz_laps)
      character*125 comment(5*nz_laps)
      character*12  cmodel
      integer len_dir, ntime, last_time,next_time, nf
      integer nxbg, nybg, nzbg(5),ntbg, ivaltimes(20)
      integer bigint
      parameter(bigint=2000000000)
c
      data msgflg/1.e30/
      data ntime/0/
      data last_time/0/
      data next_time/bigint/
      data ext/'lga'/
c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c
      lga_status=0
      call get_directory('static',outdir,len_dir)    

      call get_laps_lat_lon(outdir(1:len_dir),laps_domain_file,
     .                      nx_laps,ny_laps,lat,lon,istatus)
c      call get_laps_lat_lon(lapsroot(1:l)//'/'//laps_domain_file(1:ll)//
c     .                      '/static',laps_domain_file,
c     .                      nx_laps,ny_laps,lat,lon,istatus)
      if (istatus .ne. 1) print *,'Error reading LAPS lat, lon data.'
c
c *** Specify model path, extension for write laps routine.
c

c      l=index(lapsroot,' ')-1
c      ll=index(laps_domain_file,' ')-1
c      outdir=lapsroot(1:l)//'/'//laps_domain_file(1:ll)//'/lapsprd/lga/'
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
c *** Process background model file names.
c
      ct4=0
      call get_file_names(bgpath,bg_files,names,max_files,istatus)
      if (istatus .ne. 1) print *,'Error in get_file_names.'

      print *,'Number of bg_files found:',bg_files


      do i=1,bg_files
c         j=index(names(i),' ')-14

         call s_len(names(i),j)
         j=j-13
         if (j .ge. 0) then
            if (names(i)(j+1:j+1) .eq. '1' .or. 
     .          names(i)(j+1:j+1) .eq. '9') then
              if (bgmodel .eq. 4) then
                   fname=wfo_fname13_to_fname9(names(i)(j+1:j+13))

                   call open_sbn_netcdf(bgpath,fname,ncid,ntbg
     +                  ,istatus)
                   if(istatus.eq.0) then
                      print*,'Not enough records in file ',fname
                   else
                      call get_sbn_dims(ncid,nxbg,nybg,nzbg,ntbg
     +                                ,ivaltimes)
                      if(nxbg.lt.nx_bg.and.nybg.lt.ny_bg)then
                         cmodel='RUC 60 SBN'
                      else
                         cmodel='ETA SBN'
                      endif

                      istatus= NF_CLOSE(ncid)
                      do k=1,ntbg
                         write(af,'(i4.4)') ivaltimes(k)/3600
                         ct4=ct4+1
                         bg_names(ct4)=fname//af
                      enddo
                   endif
               else
                 ct4=ct4+1
                 bg_names(ct4)=names(i)(j+1:j+13)
               endif
            endif
         endif
      enddo
      bg_files=ct4

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
c *** Get current time from systime.dat
c
      call get_systime(i4time_now,a9,istatus)
c      call get_directory('time',lapsroot,l)
c      open(1,file=lapsroot(1:l)//'systime.dat',status='old',err=900)
c      read(1,*) i4time_now
c      close(1)
c
c *** Create new lga file if it does not already exist.
c
      do n=bg_files,1,-1
c         do i=1,lga_files
c            if (bg_names(n)(1:9) .eq. lga_names(i)(1:9) .and.
c     .          bg_names(n)(12:13) .eq. lga_names(i)(10:11)) then
c                  print *,'LGA file exists, not making one.'
c                  goto 40
c            endif
c         enddo
         fname=bg_names(n)(1:9)
         call i4time_fname_lp(fname,bgtime,istatus)
         af=bg_names(n)(10:13)
         read(af,'(i4)') ihour

c
c *** File names are returned sorted from newest to oldest ***
c     process only the newest which is at least as old as the laps time
c
      
         if(bgtime.gt.i4time_now .and. .not. use_analysis) then
           print*,
     +      'Background model newer than requested LAPS - skipping'
           print*,'This behaivior can be changed in the namelist '
           print*,'parameter use_analysis '
           goto 40
         endif  
         if(ntime.gt.bgtime) then
c           print*,'Background older than latest available - skipping'
           goto 40
         endif
            
c
c ****** Do NOT process model fcst if fcst is greater than 18 hours.
c
         if (ihour .gt. oldest_forecast) then
            print *,'IHOUR > ',oldest_forecast,' no lga file created.'
            print*,'oldest_forecast is a namelist parameter'
            goto 40
         endif
c
c ****** Do NOT process model fcst file if model init is more than 6 hours old.
c ****** Current time already read from systime.dat
c
c         if (bgtime .lt. i4time_now-21600) then
c            print *,'BGTIME < I4TIME_NOW-21600, no lga file created.'
c            goto 40   
c         endif
c
c ****** Do NOT process model fcst file if model init time is not a LAPS time.
c
         if (mod(bgtime,laps_cycle_time) .ne. 0) THEN
            print *,'Model init time is not a LAPS time',
     .              ' no lga file created.'
            goto 40
         endif
c
c ****** Do NOT process model fcst if the fcst valid time is not a LAPS time.
c
         read(af,'(i4)') ihour
         if (laps_cycle_time .ge. 3600 .and.
     .       mod(laps_cycle_time,3600) .eq. 0 .and.
     .       mod(ihour,laps_cycle_time/3600) .ne. 0) then
            print *,'Model fcst valid time is not a LAPS time.'
            goto 40
         endif
         if(ntime.eq.0) then
            ntime = bgtime
         endif   


         if(bgtime+ihour*3600.gt.i4time_now.and.
     +        bgtime+ihour*3600.lt.next_time) then
            next_time = bgtime+ihour*3600
            file_list(1) = bg_names(n)
          endif
         if(bgtime+ihour*3600.le.i4time_now.and.
     +        bgtime+ihour*3600.gt.last_time) then
            if(next_time.lt.bigint) then
               last_time = bgtime+ihour*3600
               file_list(2) = bg_names(n)
            else
c
c otherwise the latest model is incomplete and doesn't have a current forecast 
c look for an older one
c
               ntime=0
            endif
          endif
40       continue

      enddo
c
c ****** Read background model data.
c
      if (last_time.eq.0 .or. next_time.eq.bigint) then
        print*, 'not enough files:',last_time,next_time 
        lga_status = 0
        return        
      endif

      do nf=1,2
         fname = file_list(nf)(1:9)
         af = file_list(nf)(10:13)
         call i4time_fname_lp(fname,bgtime,istatus)

         do i=1,lga_files
            if (fname .eq. lga_names(i)(1:9) .and.
     .          af(3:4) .eq. lga_names(i)(10:11)) then

               print *,i,lga_names(i),':',fname,':',af
               call get_lga_source(nx_laps,ny_laps,nz_laps
     +              ,fname,af,comment(1))
               call s_len(cmodel,i)
               if(cmodel(1:i) .eq. comment(1)(1:i)) then
                  print *,'LGA file exists, not making one. ' 
     +                   ,lga_names(i)
                  lga_status=1
                  goto 80
               else
                  print *,'Overwritting LGA file ',lga_names(i)
     +            ,'from model ',comment(1)(1:i),' with a new one'
     +            ,' from model ',cmodel(1:i)
               endif
            endif
               
         enddo


         if (bgmodel .eq. 1) then     ! Process 60 km RUC data
            call read_ruc60_native(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                             prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                             gproj,istatus)
 
         elseif (bgmodel .eq. 2) then ! Process 48 km ETA conus-c grid data
           call read_eta_conusC(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                          prbg, htbg,tpbg,shbg,uwbg,vwbg,
     .                          gproj,istatus)

         elseif (bgmodel .eq. 3 .or. 
     .           bgmodel .eq. 8) then ! Process NOGAPS data
            call read_nogaps(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                       prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                       gproj,istatus)
 
         elseif (bgmodel .eq. 4) then ! Process SBN Conus 211 data (Eta or RUC)
            call read_conus_211(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                            nxbg,nybg,nzbg,ntbg,
     .                            prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                            gproj,istatus)
c
         elseif (bgmodel .eq. 5) then ! Process 40 km RUC data
            call read_ruc40_native(bgpath,fname,af,nx_bg,ny_bg,nz_bg,
     .                             prbg,htbg,tpbg,shbg,uwbg,vwbg,
     .                             gproj,istatus)
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
            goto 80
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
               if (htvi(i,j,1) .eq. msgflg) then
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
                     htvi(i,j,k)=msgflg
                     tpvi(i,j,k)=msgflg
                     shvi(i,j,k)=msgflg
                     uwvi(i,j,k)=msgflg
                     vwvi(i,j,k)=msgflg
                  enddo
               endif
            enddo
            enddo
         endif
         

c
c ****** Horizontally interpolate background data to LAPS grid points.
c
         call hinterp(nx_bg,ny_bg,nx_laps,ny_laps,nz_laps,gproj,
     .        lat,lon,
     .        htvi,tpvi,shvi,uwvi,vwvi,
     .        ht,tp,sh,uw,vw,
     .        bgmodel)
c
c ****** Check for missing value flag in any of the fields.
c ****** Check for NaN's in any of the fields.
c
         do k=1,nz_laps
         do j=1,ny_laps
         do i=1,nx_laps
            if (max(ht(i,j,k),tp(i,j,k),sh(i,j,k),
     .              uw(i,j,k),vw(i,j,k)) .ge. msgflg) then
               print *,'Missing value flag detected...Abort...',i,j,k
               print*,ht(i,j,k),tp(i,j,k),sh(i,j,k), uw(i,j,k),vw(i,j,k)
               goto 80
            endif
            if (ht(i,j,k) .ne. ht(i,j,k) .or. 
     .          tp(i,j,k) .ne. tp(i,j,k) .or.
     .          sh(i,j,k) .ne. sh(i,j,k) .or. 
     .          uw(i,j,k) .ne. uw(i,j,k) .or.
     .          vw(i,j,k) .ne. vw(i,j,k)) then
               print *,'NaN detected...Abort...',i,j,k
               goto 80
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

         call s_len(cmodel,i)
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
            comment(k)=cmodel(1:i)//' interpolated to LAPS isobaric.'
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
            comment(kk)=cmodel(1:i)//' interpolated to LAPS isobaric.'
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
            comment(kk)=cmodel(1:i)//' interpolated to LAPS isobaric.'
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
            comment(kk)=cmodel(1:i)//' interpolated to LAPS isobaric.'
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
            comment(kk)=cmodel(1:i)//' interpolated to LAPS isobaric.'
         enddo
c
         read(af,'(i4)') ihour
         bgvalid=bgtime+ihour*3600
c
         ext = 'lga'
         print *,'Writing - ',fname//af(3:4),'00.',ext

         call write_laps(bgtime,bgvalid,outdir,ext,
     .                   nx_laps,ny_laps,nz_laps*5,nz_laps*5,var,
     .                   ip,lvl_coord,units,comment,grid,istatus)
         
         if (istatus .eq. 1) then
            lga_status = 1
         else
           print *,'Error writing interpolated data to LAPS database.'
         endif
c
 80      continue
      enddo

      if(lga_status.eq.0) return

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
            print*,outdir,ext
            call time_interp(outdir,ext,
     .           nx_laps,ny_laps,nz_laps,
     .           pr,laps_cycle_time,
     .           lga_times(i-1),lga_valid(i-1),
     .           lga_times(i  ),lga_valid(i  ))
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
