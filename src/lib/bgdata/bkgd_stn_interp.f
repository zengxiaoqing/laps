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

      subroutine bkgd_stn_interp(nstns,slat,slon,selev
     &,stn_ht,stn_tp,stn_td,stn_uw,stn_vw,stn_pr,stn_mslp
     &,bkg_status)

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
     .          bkg_status                   !status returned from lga_driver
c                                            ! 1=good 0=bad
      real*4    prbot,delpr                  !LAPS bottom and delta pressures

      integer   nstns
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
      character*132 cmodel(maxbgmodels)
      character*255 cfname
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

      bg_files=0
      i=0
      bkg_status = 0 

      no_infinite_loops=0
      do while(bkg_status.le.0 .and. i.le.maxbgmodels
     +     .and. no_infinite_loops.lt.30)
         no_infinite_loops=no_infinite_loops+1
         print *,'HERE:',bkg_status, i, bg_files,reject_cnt
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
           bkg_status = 0
        else

           print *, 'bgmodel is set to ', bgmodel
           print *, ' '
           print *, 'input parameters'
           print *,  ' '
           print *, nx_bg,ny_bg,nz_bg
           print *, 'bgpath ', bgpath(1:bglen)
           print *, 'cmodel ',cmodel(i)
 970       continue
c
c *** Call bkgd driver.
c
           call  bkgd_stn_interp_sub(nstns,slat,slon,selev,
     .  bgmodel,bgpath,names,cmodel(i),max_files,bg_files,
     .  laps_cycle_time,
     .  nx_bg, ny_bg, nz_bg, i4time_now_lga, stn_ht,stn_tp,
     .  stn_td,stn_uw,stn_vw,stn_pr,stn_mslp, bkg_status)

           if(bkg_status.lt.0) then
              do l = 1,bg_files
                 reject_cnt=reject_cnt+1
                 reject_files(reject_cnt)=names(l)       ! (-bkg_status), name was names(-bkg_status)
              enddo
              bg_files = 0                               ! reset the counter
           endif
c
        endif
        
      enddo
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
      subroutine bkgd_stn_interp_sub(nstns,slat,slon,selev,
     .     bgmodel,bgpath,bg_names,cmodel,max_files,nbg,
     .     laps_cycle_time,
     .     nx_bg,ny_bg,nz_bg,i4time_now,stn_ht,stn_tp,stn_td,
     .     stn_uw,stn_vw,stn_pr,stn_mslp, bkg_status)

c
      implicit none

c     include 'bgdata.inc'
c
      integer   nbg

      integer   nstns,laps_cycle_time
     .         ,nx_bg  ,ny_bg  ,nz_bg       !Background model grid dimensions
     .         ,max_files, bkg_status
     .         ,bgmodel

      character*(*) bgpath
     +             ,bg_names(max_files)
     +             ,cmodel

      real*4    diff_mn,diff_mx
      integer   icnt
      integer i_mx,j_mx
      integer i_mn,j_mn

      logical   l_done_this
c
c *** Background model grid data.
c
      real*4    htbg(nx_bg,ny_bg,nz_bg)      !Height (m)
     .         ,prbg(nx_bg,ny_bg,nz_bg)      !Pressure (mb)
     .         ,tpbg(nx_bg,ny_bg,nz_bg)      !Temperature (K)
     .         ,shbg(nx_bg,ny_bg,nz_bg)      !Specific humidity (kg/kg)
     .         ,uwbg(nx_bg,ny_bg,nz_bg)      !U-wind (m/s)
     .         ,vwbg(nx_bg,ny_bg,nz_bg)      !V-wind (m/s)
     .         ,wwbg(nx_bg,ny_bg,nz_bg)      !W-wind (pa/s)

      real*4    htbg_sfc(nx_bg,ny_bg)

      real*4    mslpbg(nx_bg,ny_bg)         !mslp  (mb)
     .         ,prbg_sfc(nx_bg,ny_bg) 
     .         ,shbg_sfc(nx_bg,ny_bg) 
     .         ,uwbg_sfc(nx_bg,ny_bg) 
     .         ,vwbg_sfc(nx_bg,ny_bg) 
     .         ,tpbg_sfc(nx_bg,ny_bg)

      real*4    hgt1dz(1,1,nz_bg)
     .         ,tmp1dz(1,1,nz_bg)
     .         ,ssh1dz(1,1,nz_bg)
     .         ,prs1dz(1,1,nz_bg)
     .         ,uuw1dz(1,1,nz_bg)
     .         ,vvw1dz(1,1,nz_bg)

      real*4    pr(nz_bg)

c these are the station interpolated data at the background times
      real*4    tp_sfc(nstns,nbg)
     .         ,sh_sfc(nstns,nbg)
     .         ,uw_sfc(nstns,nbg)
     .         ,vw_sfc(nstns,nbg)
     .         ,pr_sfc(nstns,nbg)
     .         ,mslp(nstns,nbg)
     .         ,alt_sfc(nstns,nbg)

      real*4    slat(nstns)
     .         ,slon(nstns)
     .         ,selev(nstns)
     .         ,grx(nstns)
     .         ,gry(nstns)

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

      integer   ct,n,newfcst,
     .          ihour,
     .          lga_files,lga_times(max_files),
     .          lga_valid,
     .          i4time_lga_valid(max_files),i4time_now,
     .          bg_times(max_files),
     .          i4time_bg_valid(max_files),
     .          bg_valid(max_files),
     .          valid_bg(max_files),time_bg(max_files),
     .          bgvalid,
     .          i,ic,ii,j,jj,k,kk,l,ldl,lf,
     .          istatus,istatus_prep

      integer   nf
c
      character*256 names(max_files)
      character*13  fname13,fname9_to_wfo_fname13
      character*2   gproj
      character*256 fullname,outdir
      character*4   af_bg(max_files)

c_______________________________________________________________________________
c *** Get LAPS lat, lons.
c
      print *,'in bkgd_stn_interp sub'
      bkg_status=0

      call get_r_missing_data(rmissingflag,istatus)
      if(istatus.ne.1)then
          print*,'Error getting r_missing_data'
          return
      endif

      do j=1,max_files
         names(j)=bg_names(j)
      enddo
      do j=1,nbg
         call i4time_fname_lp(names(j)(1:9),bg_times(j),istatus)
         read(bg_names(j)(12:13),'(i2)')ihour
         bg_valid(j)=ihour*3600
         i4time_bg_valid(j)=bg_times(j)+bg_valid(j)
         af_bg(j)=names(j)(10:13)
      enddo

c
c ****** Read background model data for nbg times (=2).
c
      do nf = 1,nbg

         call s_len(bgpath,i)
         fullname = bgpath(1:i)//'/'//bg_names(nf)

         call read_bgdata(nx_bg,ny_bg,nz_bg
     +    ,bgpath,bg_names(nf),af_bg(nf),fullname,cmodel,bgmodel
     +    ,htbg, prbg, tpbg, uwbg, vwbg, shbg, wwbg
     +    ,htbg_sfc, prbg_sfc, shbg_sfc, tpbg_sfc
     +    ,uwbg_sfc, vwbg_sfc, mslpbg, gproj, istatus_prep)

         if (istatus_prep .ne. 0) then

            call s_len(bg_names(nf),lf)
            call s_len(bgpath,l)
            if (bgmodel .gt. 1 .and. bgmodel .le. 3) then
               fname13=bg_names(nf)(1:lf)//af_bg(nf)
            elseif (bgmodel .eq. 4) then
               fname13=fname9_to_wfo_fname13(bg_names(nf))
            endif
            print *,'Error reading background model data for: ',
     .         bgpath(1:l)//'/'//fname13
            print *,'Process aborted for this file.'
            bkg_status= -nf
            return
 
         endif
c
c determine grid-x and grid-y (ri and rj) for station locations
      l_done_this = .false.
      if(.not.l_done_this)then
      do i = 1,nstns
         call init_hinterp(nx_bg,ny_bg,1,1,gproj,
     . slat(i),slon(i),grx(i),gry(i),bgmodel)
      enddo
      l_done_this = .true.
      endif
c
c
c ****** Horizontally interpolate background data to station location (lat/lon).
c
         do i = 1,nstns

c
c ****** Horizontally interpolate background surface data to LAPS grid points.
c
         if(bgmodel.ne.1.and.bgmodel.ne.9)then

c need sfc T for subroutine sfcbkgd as first guess estimate.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     .        grx(i),gry(i),tpbg_sfc,tp_sfc(i,nf),bgmodel)
c simple interpolation to get the mslp.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     .        grx(i),gry(i),mslpbg,mslp(i,nf),bgmodel)
c need sfc moisture variable for subroutine sfcbkgd as first guess estimate.
            call hinterp_field(nx_bg,ny_bg,1,1,1,
     .        grx(i),gry(i),shbg_sfc,sh_sfc(i,nf),bgmodel)

c we need the 3D data interpolated to the station location.
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),prbg(1,1,1),prs1dz(1,1,1),bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),htbg(1,1,1),hgt1dz(1,1,1),bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),tpbg(1,1,1),tmp1dz(1,1,1),bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),shbg(1,1,1),ssh1dz(1,1,1),bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),uwbg(1,1,1),uuw1dz(1,1,1),bgmodel)
            call hinterp_field(nx_bg,ny_bg,1,1,nz_bg,
     .        grx(i),gry(i),vwbg(1,1,1),vvw1dz(1,1,1),bgmodel)

c compute sfc values from model using station elev and model info.
c note: sh_sfc can have q, rh, or Td depending on which model background
c       is used.  sfcbkgd deals with this and always returns
c       Td (deg K) (in array sh_sfc)
c
            call sfcbkgd(bgmodel, tmp1dz(1,1,1)
     .  ,ssh1dz(1,1,1), hgt1dz(1,1,1)
     .  ,tp_sfc(i,nf), sh_sfc(i,nf), selev(i)
     .  ,prs1dz(1,1,1), 1,1, nz_bg, pr_sfc(i,nf))
c
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
c..... Do the winds
c
            if(selev(i).gt.hgt1dz(1,1,1))then
               call interp_to_sfc(selev(i),uuw1dz,hgt1dz,
     .  1,1,nz_bg,rmissingflag,uw_sfc(i,nf))
               call interp_to_sfc(selev(i),vvw1dz,hgt1dz,
     .  1,1,nz_bg,rmissingflag,vw_sfc(i,nf))
            else
               uw_sfc(i,nf)=uuw1dz(1,1,1)
               vw_sfc(i,nf)=vvw1dz(1,1,1)
            endif
c
c now compute alt setting and other variables needed at ea stn loc.
c

         endif
         enddo
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
         print*,'No time interp when bg_files < 2'
         print*
         do i=1,nstns
            stn_tp(i)  = tp_sfc(i,nbg)
            stn_td(i)  = sh_sfc(i,nbg)
            stn_uw(i)  = uw_sfc(i,nbg)
            stn_vw(i)  = vw_sfc(i,nbg)
            stn_pr(i)  = pr_sfc(i,nbg)
            stn_mslp(i)= mslp(i,nbg)
            stn_alt(i)  = alt_sfc(i,nbg)
         enddo
      endif

      bkg_status=1

      return
      end
c
c---------------------------------------------------------------------------

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
