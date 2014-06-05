cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis
c
        subroutine get_hydro_obs(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 itest_madis_qc,l_multiple_reports,
     &                 lat,lon,ni,nj,grid_spacing,
     &                 nn,n_obs_g,n_obs_b,stations,
     &                 reptype,atype,weather,wmoid,
     &                 store_1,!store_2,store_2ea,
!    &                 store_3,store_3ea,store_4,store_4ea,
     &                 store_5,store_5ea,store_6,store_6ea,
!    &                 store_7,store_cldht,store_cldamt,
     &                 provider, laps_cycle_time, 
     &                 local_obs_thresh, i4wait_local_obs_max, jstatus)       

c*****************************************************************************
c
c.....  Input variables/arrays
c
        integer maxsta ! processed stations for LSO file
        character*(*) path_to_local_data, local_format
c
c.....  Local variables/arrays
c
	integer    maxobs
	integer    rtime
        integer  i4time_ob_a(maxobs), before, after
        real    lat(ni,nj), lon(ni,nj)
        real  k_to_f
        character*9 a9time_before, a9time_after, a9time_a(maxobs)
        logical l_reject(maxobs), ltest_madis_qc, ltest_madis_qcb
        logical l_multiple_reports, l_same_stn
        logical l_good_global,l_first_solar
c
	integer  wmoid(maxsta)
c
	character  save_stn(maxobs)*6
	character  timech*9, time*4
	character  stations(maxsta)*20
	character  provider(maxsta)*11
	character  presWeather(maxobs)*25, weather(maxsta)*25
	character  reptype(maxsta)*6, atype(maxsta)*6
	character  store_cldamt(maxsta,5)*4 
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
        character*40 string
c
c.....  Declarations for call to NetCDF reading routine (from gennet)

      include 'netcdf.inc'
      integer ICcheckNum, QCcheckNum, maxStaticIds, nInventoryBins,
     +     recNum,nf_fid, nf_vid, nf_status
      parameter (ICcheckNum=100)       ! Manually added
      parameter (QCcheckNum=100)       ! Manually added
      parameter (maxStaticIds=30000)   ! Manually added
      parameter (nInventoryBins=24)    ! Manually added
      integer filterSetNum, firstInBin(nInventoryBins), firstOverflow,
     +     globalInventory, invTime(maxobs), inventory(maxStaticIds),
     +     isOverflow(maxobs), lastInBin(nInventoryBins),
     +     lastRecord(maxStaticIds), nStaticIds,
     +     numericWMOid(maxobs), precip12hrICA(maxobs),
     +     precip12hrICR(maxobs), precip12hrQCA(maxobs),
     +     precip12hrQCR(maxobs), precip1hrICA(maxobs),
     +     precip1hrICR(maxobs), precip1hrQCA(maxobs),
     +     precip1hrQCR(maxobs), precip24hrICA(maxobs),
     +     precip24hrICR(maxobs), precip24hrQCA(maxobs),
     +     precip24hrQCR(maxobs), precip3hrICA(maxobs),
     +     precip3hrICR(maxobs), precip3hrQCA(maxobs),
     +     precip3hrQCR(maxobs), precip5minICA(maxobs),
     +     precip5minICR(maxobs), precip5minQCA(maxobs),
     +     precip5minQCR(maxobs), precip6hrICA(maxobs),
     +     precip6hrICR(maxobs), precip6hrQCA(maxobs),
     +     precip6hrQCR(maxobs), precipAccumICA(maxobs),
     +     precipAccumICR(maxobs), precipAccumQCA(maxobs),
     +     precipAccumQCR(maxobs), prevRecord(maxobs),
     +     secondsStage1_2(maxobs), secondsStage3(maxobs)
      real elevation(maxobs), latitude(maxobs), longitude(maxobs),
     +     precip12hr(maxobs), precip12hrQCD( QCcheckNum, maxobs),
     +     precip1hr(maxobs), precip1hrQCD( QCcheckNum, maxobs),
     +     precip24hr(maxobs), precip24hrQCD( QCcheckNum, maxobs),
     +     precip3hr(maxobs), precip3hrQCD( QCcheckNum, maxobs),
     +     precip5min(maxobs), precip5minQCD( QCcheckNum, maxobs),
     +     precip6hr(maxobs), precip6hrQCD( QCcheckNum, maxobs),
     +     precipAccum(maxobs), precipAccumQCD( QCcheckNum, maxobs),
     +     riverFlow(maxobs), riverStage(maxobs)
      double precision observationTime(maxobs), receivedTime(maxobs),
     +     reportTime(maxobs), riverReportChangeTime(maxobs)
      character precip5minDD(maxobs)
      character precip12hrDD(maxobs)
      character*72 ICT(ICcheckNum)
      character precip1hrDD(maxobs)
      character*11 dataProvider(maxobs)
      character*60 QCT(QCcheckNum)
      character*11 handbook5Id(maxobs)
      character precip24hrDD(maxobs)
      character precipAccumDD(maxobs)
      character*51 stationName(maxobs)
      character precip3hrDD(maxobs)
      character*11 stationType(maxobs)
      character precip6hrDD(maxobs)
      character*256 rawMessage(maxobs)
      character*24 staticIds(maxStaticIds)
      character*12 providerId(maxobs)
      character*4 homeWFO(maxobs)
      character*11 stationId(maxobs)

      real seaSurfaceTemp(maxobs) ! manually added
c
c.....  Output arrays.
c
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)

        integer ibmask(8)
c
c.....  Start.
c
 
        l_first_solar = .true.

        if(itest_madis_qc .gt. 0)then
            if(itest_madis_qc .eq. 15)then  ! call DD & QCR checking routines
                ltest_madis_qc  = .true.    ! for subjective QC reject list
                ltest_madis_qcb = .true.
                ibmask(1) = 0
                ibmask(2) = 1               ! Validity check applied
                ibmask(3) = 0
                ibmask(4) = 0
                ibmask(5) = 0
                ibmask(6) = 1               ! Statistical Spatial Consistency check
                ibmask(7) = 0
                ibmask(8) = 0
                level_qc = 0                ! Subjective QC (reject list) only
            else                            ! values of 1-2 (DD flag check)
                ltest_madis_qc  = .true.
                ltest_madis_qcb = .false.
                level_qc = itest_madis_qc
            endif
        else                                ! value of 0 (neither check routine)
            ltest_madis_qc  = .false.
            ltest_madis_qcb = .false.
        endif

             
        write(6,*)' Subroutine get_hydro_obs:' 
        write(6,*)
     1      ' itest_madis_qc/ltest_madis_qc/ltest_madis_qcb/level_qc = '   
     1       ,itest_madis_qc,ltest_madis_qc,ltest_madis_qcb,level_qc       

c
c.....	Set jstatus flag for the local data to bad until we find otherwise.
c
	jstatus = -1

        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

        call get_box_size(box_size,istatus)
        if(istatus .ne. 1)return
c
c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing / 1000.) !in km
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north

        nn_in = nn
c
c.....	Zero out the counters.
c
 10     nn = nn_in
        n_obs_g = 0	        ! # of local obs in the laps grid
        n_obs_ng = 0	        ! # of local obs not in the laps grid
        n_obs_b = 0	        ! # of local obs in the box
c
c.....  Get the data from the NetCDF file.  First, open the file.
c.....  If not there, return to obs_driver.
c
        ix = 1
c
c.....  Set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after

!       Ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 3599

        call get_filetime_range(before,after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,3600                                     
     1                         ,i4time_file_b,i4time_file_a)              

        I4_elapsed = ishow_timer()

        do i4time_file = i4time_file_a, i4time_file_b, -3600

            call s_len(path_to_local_data,len_path)
            filename13= cvt_i4time_wfo_fname13(i4time_file)
 	    data_file = path_to_local_data(1:len_path)//filename13

            write(6,*)' mesonet file = ',data_file(1:len_path+13)

!           goto 590 ! debugging test
	    nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)

	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *, data_file
	       go to 590
	    endif
c
c.....  Get the dimension of some of the variables.
c
c.....  "recNum"
c
	    nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *,'dim recNum'
	    endif
	    nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *,'dim recNum'
	    endif

            if(recnum .gt. maxobs-ix+1)then
                write(6,*)
     1              ' ERROR: exceeded maxobs limits in get_local_obs'
     1              ,ix-1,recnum,(ix-1)+recnum,maxobs
                write(6,*)' Try increasing "maxobs" in obs_driver.nl'
                go to 590
            endif

c
c.....  Call the read routine.
c
        call read_madis_hydro_netcdf(nf_fid, ICcheckNum, QCcheckNum, 
     +     maxStaticIds, nInventoryBins, recnum, elevation(ix), 
     +     latitude(ix), longitude(ix), precip12hr(ix), 
     +     precip12hrQCD(1,ix), precip1hr(ix), precip1hrQCD(1,ix), 
     +     precip24hr(ix), precip24hrQCD(1,ix), precip3hr(ix), 
     +     precip3hrQCD(1,ix), precip5min(ix), precip5minQCD(1,ix), 
     +     precip6hr(ix), precip6hrQCD(1,ix), precipAccum(ix), 
     +     precipAccumQCD(1,ix), riverFlow(ix), riverStage(ix), ICT, 
     +     QCT, dataProvider(ix), handbook5Id(ix), homeWFO(ix), 
     +     precip12hrDD(ix), precip1hrDD(ix), precip24hrDD(ix), 
     +     precip3hrDD(ix), precip5minDD(ix), precip6hrDD(ix), 
     +     precipAccumDD(ix), providerId(ix), rawMessage(ix), 
     +     staticIds, stationId(ix), stationName(ix), 
     +     stationType(ix), observationTime(ix), receivedTime(ix), 
     +     reportTime(ix), riverReportChangeTime(ix), filterSetNum, 
     +     firstInBin, firstOverflow, globalInventory, invTime(ix), 
     +     inventory, isOverflow(ix), lastInBin, lastRecord, 
     +     nStaticIds, numericWMOid(ix), precip12hrICA(ix), 
     +     precip12hrICR(ix), precip12hrQCA(ix), precip12hrQCR(ix), 
     +     precip1hrICA(ix), precip1hrICR(ix), precip1hrQCA(ix), 
     +     precip1hrQCR(ix), precip24hrICA(ix), precip24hrICR(ix), 
     +     precip24hrQCA(ix), precip24hrQCR(ix), precip3hrICA(ix), 
     +     precip3hrICR(ix), precip3hrQCA(ix), precip3hrQCR(ix), 
     +     precip5minICA(ix), precip5minICR(ix), precip5minQCA(ix), 
     +     precip5minQCR(ix), precip6hrICA(ix), precip6hrICR(ix), 
     +     precip6hrQCA(ix), precip6hrQCR(ix), precipAccumICA(ix), 
     +     precipAccumICR(ix), precipAccumQCA(ix), 
     +     precipAccumQCR(ix), prevRecord(ix), secondsStage1_2(ix), 
     +     secondsStage3(ix),badflag)

            n_local_file = recNum
            write(6,*)'     n_local_file = ',n_local_file

            ix = ix + n_local_file

            I4_elapsed = ishow_timer()

590     enddo                  ! i4time_file

        n_local_all = ix - 1
        write(6,*)' n_local_all = ',n_local_all
c
c.....  First check the data coming from the NetCDF files.  There can be
c.....  "FloatInf" (used as fill value) in some of the variables.  These
c.....  are not handled the same by different operating systems.  For 
c.....  example, IBM systems make "FloatInf" into "NaN" and store them that
c.....  way in the file, which messes up other LAPS routines.  This code
c.....  checks for "FloatInf" and sets the variable to 'badflag'.  If the
c.....  "FloatInf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
        max_write = 100

        if(n_local_all .gt. maxobs)then
           write(6,*)' Error in get_hydro_obs: n_local_all is ',
     1                                         n_local_all
           write(6,*)' Try increasing obs_driver.nl/maxobs from ',maxobs
           stop
        endif
c
c..................................
c.....	First QC loop over all the obs.
c..................................
c
	do i=1,n_local_all
           l_reject(i) = .false.
c
c........  Toss the ob if lat/lon/elev or observation time are bad by setting 
c........  lat to badflag (-99.9), which causes the bounds check to think that
c........  the ob is outside the LAPS domain.
	   if( nanf( latitude(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( longitude(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( elevation(i) ) .eq. 1 ) l_reject(i) = .true.
	   if( nanf( observationTime(i) ) .eq. 1 ) l_reject(i) = .true.

c
c.....  Bounds check: is station in the box?  Find the ob i,j location
c.....  on the LAPS grid, then check if outside past box boundary.
c
!          Test for invalid latitude
           if(latitude(i) .lt. -90 .or. latitude(i) .gt. +90.)then
               if(.true.)then
                   write(6,81,err=82)i,n_local_all
     1                               ,stationId(i)
     1                               ,latitude(i)
 81                format(2i7,1x,a8,' invalid latitude ',e12.5)
               endif
 82            l_reject(i) = .true.
               go to 105
           endif

!          Test for badflag OR (close to S Pole but not quite at it)
!          Check can also be generalized in 'latlon_to_rlapsgrid' for 'lambert'
           if(latitude(i) .lt. -89.999 .and. latitude(i) .ne. -90.)then
               l_reject(i) = .true.
               go to 105
           endif

           call latlon_to_rlapsgrid(latitude(i),longitude(i),lat,lon,       
     &                              ni,nj,ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir
     1   .or. rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) then
               if(i .le. max_write)then
                   write(6,91,err=92)i,stationId(i)
     1                               ,nint(ri_loc),nint(rj_loc)
 91                format(i6,1x,a8,' out of box ',2i12)
               endif
 92            l_reject(i) = .true.
               go to 105
           endif
c
c.....  Elevation ok?
c
	   if(elevation(i).gt.5200. .or. elevation(i).lt.-400.)then
               l_reject(i) = .true.
               go to 105
           endif

!          End of geographic location check

	   i4time_ob_a(i) = nint(observationTime(i)) + 315619200
	   call make_fnam_lp(i4time_ob_a(i),a9time_a(i),istatus)

           call filter_string(stationId(i))
c
c........  Check to see if its in the desired time window.
c
	   if(i4time_ob_a(i) .lt. before 
     1   .or. i4time_ob_a(i) .gt. after) then
               if(i .le. max_write)then
                   write(6,71,err=72)i,stationId(i)
     1                               ,a9time_a(i),i4time_ob_a(i)
     1                               ,before,after
 71		   format(i6,1x,a8,' out of time ',a11,3i12)
               endif
 72            l_reject(i) = .true.
               go to 105
           endif

!          End of time check

!          Pick closest station if multiple stations are in time window
           if(.not. l_multiple_reports)then
             do k = 1,i-1
               if(stationId(i) .eq. stationId(k))then ! possibly the same stn
                 l_same_stn = .true.
                 if(latitude(i)  .ne. latitude(k) .or.
     1              longitude(i) .ne. longitude(k)      )then
                     l_same_stn = .false.
                 endif

                 if(l_same_stn)then ! added to if block for efficiency
                   if( (.not. l_reject(i)) .and. (.not. l_reject(k)) 
     1                                                             )then       
                     i_diff = abs(i4time_ob_a(i) - i4time_sys)
                     k_diff = abs(i4time_ob_a(k) - i4time_sys)

                     if(i_diff .ge. k_diff)then
                         i_reject = i
                     else
                         i_reject = k
                     endif

                     write(6,51)i,k,stationId(i),a9time_a(i),a9time_a(k)       
     1                         ,i_reject
 51		     format(' Duplicate detected ',2i6,1x,a6,1x,a9,1x,a9
     1                     ,1x,i6)

                     l_reject(i_reject) = .true.
                   endif ! both weren't rejected
                 endif ! same station
               endif ! possibly the same stn
             enddo ! k
           endif ! l_multiple_reports
c
c
!          if( nanf( soilMoisturePercent(i)) .eq. 1 ) 
!    1               soilMoisturePercent(i) = badflag

 105       continue
c
	enddo !i

        write(6,*)' Completed 1st QC loop'

        I4_elapsed = ishow_timer()
c
c..................................
c.....	Second QC loop over all the obs.
c..................................
c
	jfirst = 1
c
	do i=1,n_local_all

           if(l_reject(i))go to 125
c
c.....  Right time, right location...

           timech = a9time_a(i)
	   time = timech(6:9)
	   read(time,*) rtime
c
c.....  Check if station is reported more than once this
c.....  time period.
c
           if(.false.)then ! we may not need this second dupe check
	     if(jfirst .eq. 1) then
	       icount = 1
	       save_stn(1) = stationId(i)
	       jfirst = 0
	       go to 150
	     endif
c
	     do k=1,icount
               if(stationId(i) .eq. save_stn(k)) then
                 write(6,*)' Rejecting duplicate ',i,k,stationId(i)
     1                    ,' ',a9time_a(i),' ',a9time_a(k)
                 go to 125
               endif
	     enddo !k
c
	     icount = icount + 1
	     save_stn(icount) = stationId(i)  ! only one...save for checking
c 
           endif ! second dupe check (set to .false.)

 150	   nn = nn + 1

           if(nn .gt. maxsta)then
              write(6,*)' ERROR in get_hydro_obs: increase maxsta '
     1                 ,nn,maxsta
              stop
           endif
 
           n_obs_b = n_obs_b + 1     !station is in the box

           call s_len2(dataProvider(i),lenp)  
c
c.....  Check if its in the LAPS grid.
c
           call latlon_to_rlapsgrid(latitude(i),longitude(i),lat,lon,       
     &                              ni,nj,ri_loc,rj_loc,istatus)
           if(  (ri_loc.lt.1. .or. ri_loc.gt.float(ni)) .OR.
     1          (rj_loc.lt.1. .or. rj_loc.gt.float(nj))     )then
               n_obs_ng = n_obs_ng + 1 ! outside the grid (inside the box)
           else
               n_obs_g = n_obs_g + 1   ! on grid...count it
           endif
c
c
c.....  Convert units for storage.  For those variables with a "change
c.....  time", check to make sure the variable was observed within the
c.....  last cycle (and that they're not just carrying an old ob for the 
c.....  current time).
c
c
c
c..... Precip 1hr
         pcp1 = badflag
         if(.true.                            .AND.
     1      precip1hr(i) .ge. 0.              .AND.    
     1      precip1hr(i) .lt. 10000.          .AND.    
     1      precip1hr(i) .ne. badflag               )then
             pcp1 = precip1hr(i) / 25.4 ! convert mm to inches
             write(6,*)' Found a hydro 1hr precip ob: '
     1                 ,pcp1,' ',dataProvider(i)(1:lenp),' '
     1                 ,precip1hrDD(i)
         endif
c
c..... Precip 3hr
         pcp3 = badflag
         if(.true.                            .AND.
     1      precip3hr(i) .ge. 0.              .AND.    
     1      precip3hr(i) .lt. 10000.          .AND.    
     1      precip3hr(i) .ne. badflag               )then
             pcp3 = precip3hr(i) / 25.4 ! convert mm to inches
             write(6,*)' Found a hydro 3hr precip ob: '
     1                 ,pcp3,' ',dataProvider(i)(1:lenp),' '
     1                 ,precip3hrDD(i)
         endif
c
c..... Precip 6hr
         pcp6 = badflag
         if(.true.                            .AND.
     1      precip6hr(i) .ge. 0.              .AND.    
     1      precip6hr(i) .lt. 10000.          .AND.    
     1      precip6hr(i) .ne. badflag               )then
             pcp6 = precip6hr(i) / 25.4 ! convert mm to inches
             write(6,*)' Found a hydro 6hr precip ob: '
     1                 ,pcp6,' ',dataProvider(i)(1:lenp),' '
     1                 ,precip6hrDD(i)
         endif
c
c..... Precip 24hr
         pcp24 = badflag
         if(.true.                            .AND.
     1      precip24hr(i) .ge. 0.             .AND.    
     1      precip24hr(i) .lt. 10000.         .AND.    
     1      precip24hr(i) .ne. badflag               )then
             pcp24 = precip24hr(i) / 25.4 ! convert mm to inches
             write(6,*)' Found a hydro 24hr precip ob: '
     1                 ,pcp24,' ',dataProvider(i)(1:lenp),' '
     1                 ,precip24hrDD(i)
         endif
c
c..... Other stuff.  
c
c
	 store_6ea(nn,1) = 0.0             ! precipitation (in)
	 store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c..... Output the data to the storage arrays
c
	 call s_len(stationId(i), len)
         if(len .ne. 0)then
             stations(nn)(1:len) = stationId(i)(1:len) ! station name
         else
             write(6,*)' Warning in get_local_obs: blank station name.'
     1                ,' Assigning name ',i
             write(stations(nn),101)i
 101	     format(i5,15x)
         endif
c
         if(lenp .ne. 0) then
	     provider(nn)(1:lenp) = dataProvider(i)(1:lenp)    ! data provider
         endif
         call filter_string(provider(nn))
c
         call s_len(stationType(i), len)
         if(len .ne. 0) then
            ilen = min(len, 6)
            atype(nn)(1:ilen) = stationType(i)(1:ilen) ! auto stn type
         endif
         call filter_string(atype(nn))
c
	 reptype(nn)(1:6) = 'LDAD  '            ! report type
	 wmoid(nn) = ibadflag                   ! WMO ID
c
	 store_1(nn,1) = latitude(i)            ! station latitude
	 store_1(nn,2) = longitude(i)           ! station longitude
	 store_1(nn,3) = elevation(i)           ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
         store_6(nn,1) = pcp1                   ! 1-h precipitation
         store_6(nn,2) = pcp3                   ! 3-h precipitation
         store_6(nn,3) = pcp6                   ! 6-h precipitation
         store_6(nn,4) = pcp24                  ! 24-h precipitation
         store_6(nn,5) = badflag                ! snow cover
c
 125     continue
       enddo !i

       I4_elapsed = ishow_timer()

       i4t_since_sys = i4time_now_gg() - i4time_sys

       write(6,*)' i4t_since_sys / nobs = ',i4t_since_sys,n_obs_b

       if(n_obs_b       .lt. local_obs_thresh .and. 
     1    i4t_since_sys .lt. i4wait_local_obs_max       )then
           write(6,*)' Waiting 60 sec for more obs'
           call snooze_gg(60.,istatus)
           go to 10
       endif
c
c..... That's it...lets go home.
c
       print *,' Found ',n_obs_b,' hydro obs in the LAPS box'
       print *,' Found ',n_obs_g,' hydro obs in the LAPS grid'
       print *,'       ',n_obs_ng,' hydro obs outside the LAPS grid'
       print *,' '
       jstatus = 1            ! everything's ok...
       return
c
 990   continue               ! no data available
       jstatus = 0
       print *,' WARNING: No data available from GET_HYDRO_OBS'
       return
c
       end

