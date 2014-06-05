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
        subroutine get_local_obs(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 itest_madis_qc,l_multiple_reports,
     &                 lat,lon,ni,nj,grid_spacing,
     &                 nn,n_obs_g,n_obs_b,stations,
     &                 reptype,atype,weather,wmoid,
     &                 store_1,store_2,store_2ea,
     &                 store_3,store_3ea,store_4,store_4ea,
     &                 store_5,store_5ea,store_6,store_6ea,
     &                 store_7,store_cldht,store_cldamt,
     &                 provider, laps_cycle_time, 
     &                 local_obs_thresh, i4wait_local_obs_max, jstatus)       

c
c*****************************************************************************
c
c	Routine to gather data from the LDAD mesonet files for LAPS.   
c
c	Changes:
c		P. Stamus  04-24-98  Original version (from get_metar_obs).
c		           05-01-98  Add soil moisture variables.          
c                          08-28-98  Updated read_local call, other stuff.
c                                        Added laps_cycle_time for time 
c                                        checks of the variables.
c                          06-21-99  Change ob location check to gridpt space.
c                                      Figure box size in gridpoint space from
c                                      user-defined size (deg) and grid_spacing.
c                          10-19-99  Added checks on each variable when doing
c                                      units conversion.
c                          01-11-00  Fixed check on ob time (overall), and 
c                                      check on time for individual variables.
c
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
        integer  recNum
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
      integer maxSensor,nf_fid, nf_vid, nf_status
      parameter (maxSensor=2)          ! Manually added
      parameter (maxPSTEntries=3000)   ! Manually added
      integer code1PST(maxPSTEntries), code2PST(maxPSTEntries), 
     +     code4PST(maxPSTEntries),
     +     altimeterQCR(maxobs), dewpointQCR(maxobs),
     +     firstOverflow, globalInventory, nStaticIds,
     +     numericWMOid(maxobs), precipAccumQCR(maxobs),
     +     precipIntensity( maxSensor, maxobs),
     +     precipRateQCR(maxobs), precipType( maxSensor, maxobs),
     +     pressChange3HourQCR(maxobs), pressChangeChar(maxobs),
     +     relHumidityQCR(maxobs), seaLevelPressureQCR(maxobs),
     +     stationPressureQCR(maxobs), temperatureQCR(maxobs),
     +     visibilityQCR(maxobs), windDirQCR(maxobs),
     +     windSpeedQCR(maxobs)
      real altimeter(maxobs), dewpoint(maxobs), elevation(maxobs),
     +     latitude(maxobs), longitude(maxobs),
     +     meanWeightedTemperature(maxobs), precipAccum(maxobs),
     +     precipRate(maxobs), pressChange3Hour(maxobs),
     +     relHumidity(maxobs),
     +     seaLevelPressure(maxobs), soilMoisturePercent(maxobs),
     +     soilTemperature(maxobs), solarRadiation(maxobs),
     +     stationPressure(maxobs), temperature(maxobs),
     +     visibility(maxobs), windDir(maxobs), windDirMax(maxobs),
     +     windGust(maxobs), windSpeed(maxobs)
      double precision observationTime(maxobs), receivedTime(maxobs),
     +     reportTime(maxobs), rhChangeTime(maxobs),
     +     stationPressChangeTime(maxobs), tempChangeTime(maxobs),
     +     windDirChangeTime(maxobs), windGustChangeTime(maxobs),
     +     windSpeedChangeTime(maxobs)
      character windDirDD(maxobs)
      character*11 stationType(maxobs)
      character windSpeedDD(maxobs)
      character relHumidityDD(maxobs)
      character stationPressureDD(maxobs)
      character altimeterDD(maxobs)
      character pressChange3HourDD(maxobs)
      character precipRateDD(maxobs)
      character*11 dataProvider(maxobs)
      character*11 namePST(maxPSTEntries)
      character*6 stationId(maxobs)
      character dewpointDD(maxobs)
      character seaLevelPressureDD(maxobs)
      character visibilityDD(maxobs)
      character precipAccumDD(maxobs)
      character*51 stationName(maxobs)
      character*12 providerId(maxobs)
      character temperatureDD(maxobs)

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

             
        write(6,*)' Subroutine get_local_obs:' 
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
c Get size of maxPSTEntries
c
            nf_status=NF_INQ_DIMID(nf_fid,'maxPSTEntries',nf_vid)
            if(nf_status.ne.NF_NOERR) then
              print *, NF_STRERROR(nf_status)
              print *,'dim maxPSTEntries'
              numPSTEntries = maxPSTEntries

              nf_status=NF_INQ_DIMLEN(nf_fid,nf_vid,numPSTEntries)
              if(nf_status.ne.NF_NOERR) then
                print *, NF_STRERROR(nf_status)
                print *,'dim maxPSTEntries'
                numPSTEntries = maxPSTEntries
              endif

            else
              numPSTEntries = maxPSTEntries

            endif

            if(numPSTEntries .gt. maxPSTEntries)then
              write(6,*)' ERROR, maxPSTEntries should be increased '
     1                 ,maxPSTEntries,numPSTEntries
              istatus = 0
              return
            else
              write(6,*)' numPSTEntries/maxPSTEntries = ',
     1                    numPSTEntries,maxPSTEntries
            endif

c
c.....  Call the read routine.
c
            call read_ldad_madis_netcdf(nf_fid, maxSensor, recNum, 
     +     maxPSTEntries, code1PST, code2PST, namePST,
     +     altimeterQCR(ix), dewpointQCR(ix), firstOverflow, 
     +     globalInventory, nStaticIds, numericWMOid, 
     +     precipAccumQCR(ix), precipIntensity, 
     +     precipRateQCR(ix), precipType, pressChange3HourQCR(ix), 
     +     pressChangeChar, relHumidityQCR(ix), seaLevelPressureQCR(ix),
     +     stationPressureQCR(ix), temperatureQCR(ix), 
     +     visibilityQCR(ix),       
     +     windDirQCR(ix), windSpeedQCR(ix), altimeter(ix), 
     +     dewpoint(ix), 
     +     elevation(ix), latitude(ix), longitude(ix), 
     +     meanWeightedTemperature(ix), precipAccum(ix), 
     +     precipRate(ix), pressChange3Hour(ix), relHumidity(ix), 
     +     seaLevelPressure(ix), seaSurfaceTemp(ix), 
     +     soilMoisturePercent(ix), soilTemperature(ix), 
     +     solarRadiation(ix), 
     +     stationPressure(ix), temperature(ix), visibility(ix), 
     +     windDir(ix), windDirMax(ix), windGust(ix), windSpeed(ix), 
     +     altimeterDD(ix), dataProvider(ix), dewpointDD(ix), 
     +     precipAccumDD(ix), precipRateDD(ix), presWeather(ix), 
     +     pressChange3HourDD(ix), providerId(ix), relHumidityDD(ix), 
     +     seaLevelPressureDD(ix), stationId(ix), stationName(ix), 
     +     stationPressureDD(ix), stationType(ix), temperatureDD(ix), 
     +     visibilityDD(ix), windDirDD(ix), windSpeedDD(ix), 
     +     observationTime(ix), receivedTime(ix), reportTime(ix), 
     +     rhChangeTime(ix), stationPressChangeTime(ix), 
     +     tempChangeTime(ix), windDirChangeTime(ix), 
     +     windGustChangeTime(ix), windSpeedChangeTime(ix),badflag)

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
      
c
c..................................
c.....	First QC loop over all the obs.
c..................................
c
        if(n_local_all .gt. maxobs)then
           write(6,*)' Error in get_local_obs: n_local_all is ',
     1                                         n_local_all       
           write(6,*)' Try increasing obs_driver.nl/maxobs from ',maxobs
           stop   
        endif

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
              if(.not. l_reject(k))then 
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
 51		     format(' Duplicate detected ',2i7,1x,a6,1x,a9,1x,a9
     1                     ,1x,i6)

                     l_reject(i_reject) = .true.
                   endif ! both weren't rejected
                 endif ! same station
               endif ! possibly the same stn
              endif ! more efficient l_reject(k) test
             enddo ! k
           endif ! l_multiple_reports
c
c
	   if( nanf( rhChangeTime(i)   ) .eq. 1 ) 
     1               rhChangeTime(i)   = ibadflag
	   if( nanf( tempChangeTime(i)    ) .eq. 1 ) 
     1               tempChangeTime(i)    = ibadflag
	   if( nanf( stationPressChangeTime(i)    ) .eq. 1 ) 
     1               stationPressChangeTime(i)    = ibadflag
	   if( nanf( windDirChangeTime(i)   ) .eq. 1 ) 
     1               windDirChangeTime(i)   = ibadflag
	   if( nanf( windSpeedChangeTime(i)   ) .eq. 1 ) 
     1               windSpeedChangeTime(i)   = ibadflag
	   if( nanf( windGustChangeTime(i) ) .eq. 1 ) 
     1               windGustChangeTime(i) = ibadflag
c
	   if( nanf( visibility(i)  ) .eq. 1 ) visibility(i) = badflag       
	   if( nanf( seaLevelPressure(i) ) .eq. 1 ) 
     1               seaLevelPressure(i)  = badflag
	   if( nanf( temperature(i) ) .eq. 1 ) temperature(i) = badflag       
	   if( nanf( dewpoint(i)   ) .eq. 1 ) dewpoint(i) = badflag
	   if( nanf( solarRadiation(i)) .eq. 1 ) 
     1               solarRadiation(i) = badflag
	   if( nanf( seaSurfaceTemp(i)) .eq. 1 ) 
     1               seaSurfaceTemp(i) = badflag
	   if( nanf( soilTemperature(i)) .eq. 1 ) 
     1               soilTemperature(i) = badflag
	   if( nanf( soilMoisturePercent(i)) .eq. 1 ) 
     1               soilMoisturePercent(i) = badflag
	   if( nanf( windDir(i)   ) .eq. 1 ) windDir(i) = badflag       
	   if( nanf( windSpeed(i)   ) .eq. 1 ) windSpeed(i) = badflag
	   if( nanf( windGust(i)  ) .eq. 1 ) windGust(i)   = badflag
	   if( nanf( altimeter(i)  ) .eq. 1 ) altimeter(i)   = badflag

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
              write(6,*)' ERROR in get_local_obs: increase maxsta '
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
c.....	Figure out the cloud data.
c.....     NOTE: Not currently reading cloud data from mesonets.
c
           kkk = 0               ! number of cloud layers
c
c
c.....  Convert units for storage.  For those variables with a "change
c.....  time", check to make sure the variable was observed within the
c.....  last cycle (and that they're not just carrying an old ob for the 
c.....  current time).
c
c.....  Temperature, dewpoint and RH.
c
          write(6,*)' Testing ',stationId(i)

	  temp_k = temperature(i) 
	  if(tempChangeTime(i) .gt. 0.) then ! implies that it is a good value
	     if( abs(observationTime(i) - tempChangeTime(i)) 
     1                          .gt. laps_cycle_time) then
		temp_k = badflag
	     endif
	  endif
	  if(temp_k .le. badflag) then !t bad?
	     temp_f = badflag	!then bag it
	  else
             temp_f = k_to_f(temp_k)
	  endif
          call sfc_climo_qc_r('t_f',temp_f)
          if(ltest_madis_qc)
     1        call madis_qc_r(temp_f,temperatureDD(i),level_qc,badflag)       
          if(ltest_madis_qcb)
     1        call madis_qc_b(temp_f,temperatureQCR(i),ibmask,0,badflag) 
c       
	  dewp_k = dewpoint(i)
          call sfc_climo_qc_r('td_k',dewp_k)
	  if(dewp_k .le. badflag) then !dp bad?
	     dewp_f = badflag	       !then bag it
	  else
	     dewp_f = k_to_f(dewp_k)
	  endif
          if(ltest_madis_qc)
     1        call madis_qc_r(dewp_f,dewpointDD(i),level_qc,badflag)       
          if(ltest_madis_qcb)
     1        call madis_qc_b(dewp_f,dewpointQCR(i),ibmask,0,badflag)       
c
	  rh_p = relHumidity(i) 
	  if(rh_p.lt.0. .or. rh_p.gt.100.) rh_p = badflag
	  if(rhChangeTime(i) .gt. 0.) then
	     if( abs(observationTime(i) - rhChangeTime(i)) 
     1                             .gt. laps_cycle_time) then
		rh_p = badflag
	     endif
	  endif
          if(ltest_madis_qc)
     1        call madis_qc_r(rh_p,relHumidityDD(i),level_qc,badflag)        
          if(ltest_madis_qcb)
     1        call madis_qc_b(rh_p,relHumidityQCR(i),ibmask,0,badflag)        
c
c..... Wind speed and direction
c
	  dir = windDir(i) 
          call sfc_climo_qc_r('dir_deg',dir)
	  spd = windSpeed(i)
          call sfc_climo_qc_r('spd_ms',spd)
	  if(windDirChangeTime(i).gt.0. .and. 
     1       windSpeedChangeTime(i).gt.0.     ) then       
	     if( (abs(observationTime(i) - windDirChangeTime(i)) 
     &                          .gt. laps_cycle_time) .or.
     &           (abs(observationTime(i) - windSpeedChangeTime(i)) 
     &                          .gt. laps_cycle_time)      ) then
		dir = badflag
		spd = badflag
	     endif
	  endif
          if(ltest_madis_qc)
     1        call madis_qc_r(dir,windDirDD(i),level_qc,badflag)        
          if(ltest_madis_qcb)
     1        call madis_qc_b(dir,windDirQCR(i),ibmask,1,badflag)        
          if(ltest_madis_qc)
     1        call madis_qc_r(spd,windSpeedDD(i),level_qc,badflag)        
          if(ltest_madis_qcb)
     1        call madis_qc_b(spd,windSpeedQCR(i),ibmask,1,badflag)        
	  if(spd .ne. badflag) spd = 1.94254 * spd !m/s to kt
c
	  dirgust = windDirMax(i)
          call sfc_climo_qc_r('dir_deg',dirgust)
	  spdgust = windGust(i)
          call sfc_climo_qc_r('spd_ms',spdgust)
	  if(windGustChangeTime(i) .gt. 0.) then
	     if( abs(observationTime(i) - windGustChangeTime(i)) 
     1                                      .gt. laps_cycle_time) then
		dirgust = badflag
		spdgust = badflag
	     endif
	  endif

	  if(spdgust .ne. badflag) spdgust = 1.94254 * spdgust !m/s to kt

	  if(spdgust .ne. badflag .and. spd .ne. badflag)then
             if(spd .gt. spdgust)then
                write(6,*)' Speed exceeds Gust at ',trim(stationId(i))       
     1                                             ,spd,spdgust
		dir = badflag
		spd = badflag
		dirgust = badflag
		spdgust = badflag
             endif
          endif 
c
c..... Pressure...Station pressure, MSL and altimeter
c
	  stn_press = stationPressure(i)
          call sfc_climo_qc_r('stnp_pa',stn_press)
          if(ltest_madis_qc)
     1        call madis_qc_r(stn_press,stationPressureDD(i),level_qc
     1                                                      ,badflag)
          if(ltest_madis_qcb)
     1        call madis_qc_b(stn_press,stationPressureQCR(i),ibmask
     1                                                    ,0,badflag)
	  if(stationPressChangeTime(i) .gt. 0.) then
	     if( abs(observationTime(i) - stationPressChangeTime(i))
     1                               .gt. laps_cycle_time ) then
		stn_press = badflag
	     endif
	  endif

	  if(stn_press .ne. badflag) stn_press = stn_press * 0.01 !Pa to mb
c
          call sfc_climo_qc_r('mslp_pa',seaLevelPressure(i))
          if(ltest_madis_qc)
     1        call madis_qc_r(seaLevelPressure(i),seaLevelPressureDD(i)
     1                                           ,level_qc,badflag)       
          if(ltest_madis_qcb)
     1        call madis_qc_b(seaLevelPressure(i),seaLevelPressureQCR(i)       
     1                                           ,ibmask,0,badflag)       
	  if(seaLevelPressure(i) .ne. badflag) seaLevelPressure(i)   
     1                             = seaLevelPressure(i)   * 0.01 !Pa to mb

          call sfc_climo_qc_r('alt_pa',altimeter(i))
          if(ltest_madis_qc)
     1        call madis_qc_r(altimeter(i),altimeterDD(i),level_qc
     1                                                   ,badflag)
          if(ltest_madis_qcb)
     1        call madis_qc_b(altimeter(i),altimeterQCR(i),ibmask,0
     1                                                   ,badflag)
	  if(altimeter(i) .ne. badflag) 
     1                         altimeter(i) = altimeter(i) * 0.01 !Pa to mb
c
c..... Visibility
c
	 if(visibility(i).lt.0. .or. visibility(i).gt.330000.) then
	   visibility(i) = badflag
	 else
	   visibility(i) = visibility(i) * .001      !m to km
	   visibility(i) = 0.621371 * visibility(i)  !km to miles
	 endif

c
c..... Solar Radiation
c
         solar_rad = solarRadiation(i)                         
!        call filter_string(dataProvider(i))
         if(solar_rad .le. badflag .or. solar_rad .ge. 2000.) then !  bad?
            solar_rad = badflag                                    !  bag
         else ! good solar
            write(6,*)' solar dataProvider = '
     1               ,stationId(i),dataProvider(i)(1:lenp)

            l_good_global = .false.

            do iprov = 1,numPSTEntries                  
               if(l_first_solar)write(6,*)iprov,namePST(iprov)           

               if(namePST(iprov) .eq. dataProvider(i)(1:lenp))then
!                 write(6,*)' match,code2PST = ',code2PST(iprov)

                  l_good_global = .false.
                  if(code2PST(iprov) .eq. 0)then
                     string = ' not defined'
                  elseif(code2PST(iprov) .eq. 1)then
                     string = ' diffuse 15 min'
                  elseif(code2PST(iprov) .eq. 2)then
                     string = ' diffuse 1 hr'  
                  elseif(code2PST(iprov) .eq. 3)then
                     string = ' diffuse 24 hr' 
                  elseif(code2PST(iprov) .eq. 4)then
                     string = ' direct 15 min'  
                  elseif(code2PST(iprov) .eq. 5)then
                     string = ' direct 1 hr'   
                  elseif(code2PST(iprov) .eq. 6)then
                     string = ' direct 24 hr'  
                  elseif(code2PST(iprov) .eq. 7)then
                     string = ' global 15 min'  
                     l_good_global = .true.
                  elseif(code2PST(iprov) .eq. 8)then
                     string = ' global 1 hr'    
                  elseif(code2PST(iprov) .eq. 9)then
                     string = ' global 24 hr'   
                  elseif(code2PST(iprov) .eq. 10)then
                     string = ' diffuse 5 min'  
                  elseif(code2PST(iprov) .eq. 11)then
                     string = ' direct 5 min'  
                  elseif(code2PST(iprov) .eq. 12)then
                     string = ' global 5 min'             
                     l_good_global = .true.
                  elseif(code2PST(iprov) .eq. 13)then
                     string = ' diffuse instantaneous'
                  elseif(code2PST(iprov) .eq. 14)then
                     string = ' direct instantaneous'
                  elseif(code2PST(iprov) .eq. 15)then
                     string = ' global instantaneous'
                     l_good_global = .true.
                  elseif(code2PST(iprov) .eq. 16)then
                     string = ' diffuse 15 min'
                  else
                     string = ' invalid code'
                  endif ! determine code
               endif ! found a provider match

            enddo ! iprov

!           write(6,*)' l_good_global = ',l_good_global
            l_first_solar = .false.

            if(.not. l_good_global)then
               solar_rad = badflag                        !  bag
            else
               call s_len2(string,lens)
               write(6,*)' Good global solar ' 
     1                   ,stationId(i),dataProvider(i)(1:lenp)
     1                   ,string(1:lens)
            endif

         endif ! solar QC check                    
c
c..... Sea Surface Temperature
c
         seatemp_k = seaSurfaceTemp(i)                         
         call sfc_climo_qc_r('sst_k',seatemp_k)
         if(seatemp_k .ne. badflag) then          
            seatemp_f = k_to_f(seatemp_k)
         else
            seatemp_f = badflag
         endif

c
c..... Soil Surface Temperature
c
         soiltemp_k = soilTemperature(i)                         
         call sfc_climo_qc_r('tgd_k',soiltemp_k)
         if(soiltemp_k .ne. badflag) then          
            soiltemp_f = k_to_f(soiltemp_k)
         else
            soiltemp_f = badflag
         endif

c
c..... Soil Moisture
c
         soilmoist_p = soilMoisturePercent(i)                         
         if(soilmoist_p.lt.0. .or. soilmoist_p.gt.100.) 
     1      soilmoist_p = badflag
c
c
c..... Fill the expected accuracy arrays.  Values are based on information
c..... in the 'Federal Meteorological Handbook No. 1' for the METARs, 
c..... Appendix C (http://www.nws.noaa.gov/oso/oso1/oso12/fmh1/fmh1appc.htm)
c..... Here however, we know that the local data has wide variations in 
c..... quality, so for now we double the FMH-1 numbers.  Later, we may be
c..... able to better define these numbers as we gain experience with the
c..... different stations that the providers use.
c
c..... Note also that we convert the units in Appendix C to match what we're 
c..... using here.
c
c..... Temperature (deg F)
c
	 fon = 9. / 5.  !ratio when converting C to F
	 store_2ea(nn,1) = 10.0 * fon        ! start...we don't know what we have
	 if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-62.) .and. temp_f.le.c2f(-50.)) then
	      store_2ea(nn,1) = 2.2 * fon  ! conv to deg F
	   elseif(temp_f.gt.c2f(-50.) .and. temp_f.lt.c2f(50.)) then
	      store_2ea(nn,1) = 1.2 * fon  ! conv to deg F
	   elseif(temp_f.ge.c2f(50.) .and. temp_f.le.c2f(54.)) then
	      store_2ea(nn,1) = 2.2 * fon  ! conv to deg F
	   endif
	 endif
c
c..... Dew point (deg F).  Also estimate a RH accuracy based on the dew point.
c..... Estimates for the RH expected accuracy are from playing around with the
c..... Psychrometric Tables for various T/Td combinations (including their
c..... accuracies from the FMH-1 Appendix C).
c
	 store_2ea(nn,2) = 10.0 * fon       ! start...don't know what we have 
	 if(dewp_f .ne. badflag) then
	    if(dewp_f.ge.c2f(-34.) .and. dewp_f.lt.c2f(-24.)) then
	       store_2ea(nn,2) = 2.2 * fon ! conv to deg F
	    elseif(dewp_f.ge.c2f(-24.) .and. dewp_f.lt.c2f(-1.)) then
	       store_2ea(nn,2) = 1.7 * fon ! conv to deg F
	    elseif(dewp_f.ge.c2f(-1.) .and. dewp_f.le.c2f(30.)) then
	       store_2ea(nn,2) = 1.1 * fon ! conv to deg F
	    endif
	 endif
	 store_2ea(nn,3) = 50.0            ! Relative Humidity %
	 if(rh_p .ne. badflag) then
	    if(rh_p .lt. 30.) then
	       store_2ea(nn,3) = 20.0      ! RH (%) 
	    elseif(rh_p.ge.30. .and. rh_p.lt.80.) then
	       store_2ea(nn,3) = 12.0      ! RH (%) 
	    elseif(rh_p.ge.80.) then
	       store_2ea(nn,3) = 8.0       ! RH (%) 
	    endif
	 endif
c
c..... Wind direction (deg) and speed (kts)
c
	 store_3ea(nn,1) = 15.0    ! deg 
	 store_3ea(nn,2) = 10.0    ! kt
	 if(spd .ne. badflag) then
	    if(spd.ge.1.0 .and. spd.lt.10.0) then
	       store_3ea(nn,2) = 2.0          ! kt
	    elseif((spd .gt. 10.0) .and.(spd .lt. 200.0)) then
	       store_3ea(nn,2) = spd * 0.2  ! 20% of speed (kts)
	    endif
c
	    if(spd .ge. 5.0) then    ! dir check
	       store_3ea(nn,1) = 10.0   ! deg
	    endif
	 endif
c
c..... Pressure and altimeter (mb)
c
	 store_4ea(nn,1) = 2.00            ! pressure (mb)
	 store_4ea(nn,2) = 2.00            ! altimeter (mb)
c
c..... Visibility (miles).  For automated stations use a guess based 
c..... on Table C-2 in Appendix C of FMH-1.  For manual stations, use
c..... a guess based on the range between reportable values (e.g., for
c..... reported visibility between 0 and 3/8th mile, set accuracy to 
c..... 1/16th mile).  This isn't ideal, but its a start.
c
	 store_5ea(nn,1) = 10.00         ! Start with this (miles)
	 if(visibility(i) .ne. badflag) then
	    if(visibility(i) .lt. 2.0) then
	       store_5ea(nn,1) = 0.50 ! miles
	    elseif(visibility(i).ge.2.0 .and. visibility(i).lt.3.0) then       
	       store_5ea(nn,1) = 1.00 ! miles
	    elseif(visibility(i) .gt. 3.0) then
	       store_5ea(nn,1) = 2.00 ! miles
	    endif
	 endif
c
c
c..... Precip
         pcp1 = badflag
         pcp24 = badflag

         if(precipAccum(i) .ge. 0.      .AND.    
     1      precipAccum(i) .ne. badflag       )then

            iprov_pst = 0
            do iprov = 1,numPSTEntries                  
               if(namePST(iprov) .eq. dataProvider(i)(1:lenp))then
!                 write(6,*)' match,code1PST = ',code1PST(iprov)
                  iprov_pst = iprov
               endif ! found a provider match
            enddo ! iprov

            if(dataProvider(i)(1:lenp) .eq. 'HADS' .AND.
     1         precipRateDD(i) .ne. 'Z')then
               pcp1 = precipRate(i) * (3600. / .0254) ! convert m/s to inches
               if(pcp1 .lt. 0. .OR. pcp1 .ge. 1e3)pcp1 = badflag       
               write(6,*)' Found a local 1hr HADS precip rate ob: '
     1                 ,pcp1,' ',dataProvider(i)(1:lenp),' '
     1                 ,stationId(i),' '
     1                 ,precipRateDD(i)                                 
            endif

            if(code1PST(iprov_pst) .eq. -3 .AND. 
     1         dataProvider(i)(1:lenp) .eq. 'HADS')then  
               pcp24 = precipAccum(i) / 25.4 ! convert mm to inches
               if(pcp24 .lt. 0. .OR. pcp24 .ge. 1e3)pcp24 = badflag       
               if(precipAccumDD(i) .eq. 'X')pcp24 = badflag
               write(6,*)' Found a local 24hr HADS precip ob    : '
     1                 ,pcp24,' ',dataProvider(i)(1:lenp),' '
     1                 ,stationId(i),' '
     1                 ,precipAccumDD(i),' ',code1PST(iprov_pst)       
            elseif(code1PST(iprov_pst) .eq. -2 .AND. ! 24hr ob at 00UT
     1             i4time_sys .eq. ((i4time_sys/86400) * 86400) )then                         
               pcp24 = precipAccum(i) / 25.4 ! convert mm to inches
               if(pcp24 .lt. 0. .OR. pcp24 .ge. 1e3)pcp24 = badflag       
               write(6,*)' Found a local 24hr 00UTC precip ob: '
     1                 ,pcp24,' ',dataProvider(i)(1:lenp),' '
     1                 ,stationId(i),' '
     1                 ,precipAccumDD(i),' ',code1PST(iprov_pst)       
            elseif(code1PST(iprov_pst) .ne. 0)then
               write(6,*)' Potential precip ob: '                             
     1                 ,precipAccum(i) / 25.4,' '
     1                 ,dataProvider(i)(1:lenp),' '
     1                 ,stationId(i),' '
     1                 ,precipAccumDD(i),' ',code1PST(iprov_pst)       
            endif ! type of precip ob
         endif ! precip was reported
c
c..... Other stuff.  
c
	 store_5ea(nn,2) = 0.0             ! solar radiation 
	 store_5ea(nn,3) = 1.0 * fon       ! soil/water temperature (F)
	 store_5ea(nn,4) = 0.0             ! soil moisture (rh %)
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
         weather(nn)(1:25) = presWeather(i)(1:25) ! present weather
         call filter_string(weather(nn))

	 reptype(nn)(1:6) = 'LDAD  '            ! report type
	 wmoid(nn) = ibadflag                   ! WMO ID
c
	 store_1(nn,1) = latitude(i)            ! station latitude
	 store_1(nn,2) = longitude(i)           ! station longitude
	 store_1(nn,3) = elevation(i)           ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = dewp_f                 ! dew point (deg f) 
	 store_2(nn,3) = rh_p                   ! Relative Humidity (%)
c
	 store_3(nn,1) = dir                    ! wind dir (deg)
	 store_3(nn,2) = spd                    ! wind speed (kt)
	 store_3(nn,3) = dirgust                ! wind gust dir (deg)
	 store_3(nn,4) = spdgust                ! wind gust speed (kt)
c
         store_4(nn,1) = altimeter(i)           ! altimeter setting (mb)
         store_4(nn,2) = stn_press              ! station pressure (mb)
         store_4(nn,3) = seaLevelPressure(i)    ! MSL pressure (mb)
         store_4(nn,4) = badflag                ! 3-h press change character
         store_4(nn,5) = badflag                ! 3-h press change (mb)
c
         store_5(nn,1) = visibility(i)          ! visibility (miles)
         store_5(nn,2) = solar_rad              ! solar radiation (Watts/M**2)

         if(seatemp_f .ne. badflag)then
	     store_5(nn,3) = seatemp_f          ! soil/water temperature (F)
         else
	     store_5(nn,3) = soiltemp_f         ! soil/water temperature (F)
         endif

         store_5(nn,4) = soilmoist_p            ! soil moisture (%)
c
         store_6(nn,1) = pcp1                   ! 1-h precipitation
         store_6(nn,2) = badflag                ! 3-h precipitation
         store_6(nn,3) = badflag                ! 6-h precipitation
         store_6(nn,4) = pcp24                  ! 24-h precipitation
         store_6(nn,5) = badflag                ! snow cover
c
         store_7(nn,1) = float(kkk)             ! number of cloud layers
         store_7(nn,2) = badflag                ! 24-h max temperature
         store_7(nn,3) = badflag                ! 24-h min temperature
c
c.....	Store cloud info if we have any. 
c
	 if(kkk .gt. 0) then
	   do ii=1,kkk
	     store_cldht(nn,ii) = badflag  !ht(ii,i)
	     store_cldamt(nn,ii)(1:1) = ' '
	     store_cldamt(nn,ii)(2:4) = '   '  !cvr(ii,i)(1:3)
	   enddo !ii
	 endif
c
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
       print *,' Found ',n_obs_b,' local obs in the LAPS box'
       print *,' Found ',n_obs_g,' local obs in the LAPS grid'
       print *,'       ',n_obs_ng,' local obs outside the LAPS grid'
       print *,' '
       jstatus = 1            ! everything's ok...
       return
c
 990   continue               ! no data available
       jstatus = 0
       print *,' WARNING: No data available from GET_LOCAL_OBS'
       return
c
       end

