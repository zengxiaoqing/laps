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

c     use laps_static
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
     .          nx_bg,ny_bg,nz_bg,
     .          lga_status                   !status returned from lga_driver
c                                            ! 1=good 0=bad
      integer   np
c
c
c------------------> BACKGROUND GRID DATA PATH <--------------------------------
c
c *** The following variables specify the various model data directory paths.
c *** Read from standard input.
c *** Each line contains the background data directory path for each
c        corresponding bgmodel defined above.
c
      character*256 cfilespec

      character*2   gproj
      real          dlat,dlon
      real          Lat0,Lat1,Lon0
      real          sw(2),ne(2)

c
c  This is the max number of paths allowed in nest7grid.parms
c  and should match the value in lib/lapsgrid.f 
c 
      character*256 bgpaths(maxbgmodels)
      character*256 bgpath
      integer bgmodels(maxbgmodels)
      integer bgmodel
      integer istatus
 
      character*13 cvt_i4time_wfo_fname13
c
      character*9 a9
      integer i4time_now,i4time_latest
      integer i4time_now_lga,i4time_now_gg
      integer max_files,bg_files
      integer itime_inc
      integer itime
c     parameter (max_files=600)
      parameter (max_files=5000)
      character*256 names(max_files)
      character*256 reject_names(max_files)
      integer reject_cnt
      integer accepted_files
      integer lbgp
      integer nbgmodels
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
      character*132 cmodels(maxbgmodels)
      character*132 cmodel
      integer oldest_forecast, max_forecast_delta
      integer forecast_length
      logical use_analysis, use_systime
      logical ltime(-1:1)
      logical smooth_fields
c_______________________________________________________________________________
c
c Read information from static/nest7grid.parms

      print *,'istatus ',istatus
      call get_grid_dim_xy(nx_laps,ny_laps,istatus)
      call get_laps_dimensions(nz_laps,istatus)
      call get_laps_cycle_time(laps_cycle_time,istatus)
c
c Read information from static/background.nl
c
      call get_background_info(bgpaths,bgmodels,oldest_forecast
     +,max_forecast_delta,forecast_length,use_analysis,cmodels
     +,itime,smooth_fields)

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
         call get_systime(i4time_now,a9,lga_status)
         print*,'Systime: ',a9,' ',i4time_now
      else
         i4time_now = i4time_now_gg()
         print*,'Using i4time now'
      endif

      print*,'Time adjustment calc from namelist itime_inc',
     +' (itime_inc) = ',itime_inc

      ltime=.true.

      do itime_inc=-1,+1

       print*
       print*,'----------------------------------------------'
       if(itime_inc.lt.0)then
        print*,'Start time-minus-one cycle time interpolation'
       elseif(itime_inc.eq.0)then
        print*,'Start cycle time interpolation'
       else
        print*,'Start time-plus-one cycle time interpolation'
       endif
       print*,'----------------------------------------------'
       print*

       i4time_now_lga=i4time_now+(itime_inc*laps_cycle_time)
       call make_fnam_lp(i4time_now_lga,a9,istatus)
       print*,'processing background data for ',a9
       print*

       bg_files=0
       i=1

      bgmodel = bgmodels(i)
      bgpath =  bgpaths(i)
      call s_len(bgpath,lbgp)
      cmodel =  cmodels(i)
      lga_status = -99 

      reject_cnt=0
      bg_files=0
      no_infinite_loops=0

      do while((lga_status.le.0 .and. i.le.nbgmodels)
     +         .and. (no_infinite_loops.le.nbgmodels))


         print *,'HERE:',lga_status, i, bg_files,reject_cnt

         if(reject_cnt.gt.bg_files )then    !.and. bg_files.gt.0) then
            i=i+1
            bgmodel = bgmodels(i)
            bgpath =  bgpaths(i)
            cmodel =  cmodels(i)
            reject_cnt = 0
         endif

         if (bgmodel .lt. 0 .or. bgmodel .gt. maxbgmodels) then
            print*
            print*,' Cannot proceed with model specification in LGA'
            print*,' Check bgpaths in static/background.nl'
            print*,'bgmodel=',bgmodel,' maxbgmodels=,',maxbgmodels
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
           if(i.eq.nbgmodels)lga_status = 0 

        elseif(accepted_files.eq.0.and.reject_cnt.eq.bg_files)then

              reject_cnt=reject_cnt+1

        else

           print *, ' '
           print *, 'Input Parameters'
           print *, '----------------'
           print *
           print *, ' Analysis setup: '
           print *, nx_laps,ny_laps,nz_laps
           print *, laps_cycle_time
           print *
c
c *** Call lga driver and, if necessary, interpolate acceptable files.
c
           call lga_driver(nx_laps,ny_laps,nz_laps,
     .          laps_cycle_time,bgmodel,bgpath,cmodel,reject_cnt,
     .          reject_names,names,max_files,accepted_files,
     .          i4time_now_lga, smooth_fields,lga_status)

c
c these constructs force t-1 and t+1 cycle time background generation.
c these should be removed when lga runs outside sched.pl
c
           print*
           if(lga_status.eq.1.and.ltime(itime_inc))then
c             lga_status = -99 
              ltime(itime_inc) = .false.
              print*,'----------------------------------------------'
              if(itime_inc.lt.0)then
                 print*,'Completed time-minus-one cycle time interp'
              elseif(itime_inc.eq.0)then
                 print*,'Completed cycle time interpolation'
              else
                 print*,'Completed time-plus-one cycle time interp'
              endif
              print*,'----------------------------------------------'
           endif
           print*

           if(lga_status.le.0)no_infinite_loops=no_infinite_loops+1

        endif
        
      enddo !do while

      enddo !itime_inc = -1,+1

      if(no_infinite_loops.gt.nbgmodels) then
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
965   continue
      print *,'No acceptable background model found.'
      
c
      end
c ===============================================================
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
     &icnt,i_mx,j_mx,i_mn,j_mn,dmax,dmin)

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
