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
        subroutine get_local_cwb(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 eastg,westg,anorthg,southg,
     &                 lat,lon,ni,nj,grid_spacing,
     &                 nn,n_obs_g,n_obs_b,stations,
     &                 reptype,atype,weather,wmoid,
     &                 store_1,store_2,store_2ea,
     &                 store_3,store_3ea,store_4,store_4ea,
     &                 store_5,store_5ea,store_6,store_6ea,
     &                 store_7,store_cldht,store_cldamt,
     &                 provider, laps_cycle_time, jstatus)
c
	include 'netcdf.inc'
c
c.....  Input variables/arrays
c
        integer maxobs ! raw data file
        integer maxsta ! processed stations for LSO file
        character*(*) path_to_local_data, local_format
c
c.....  Local variables/arrays
c
	real lat(ni,nj), lon(ni,nj)
	real lats(maxobs), lons(maxobs), elev(maxobs)
        real t(maxobs), t24max(maxobs), t24min(maxobs), td(maxobs)
	real rh(maxobs) 
        real pcp1hr(maxobs), pcp3hr(maxobs), pcp6hr(maxobs)
        real pcp24hr(maxobs)
        real dd(maxobs), ff(maxobs), wgdd(maxobs), wgff(maxobs)
        real stnp(maxobs), mslp(maxobs), pc(maxobs)
        real sr(maxobs), st(maxobs)

        integer    wmoid(maxobs), pcc(maxobs)
        integer    rtime
        integer  i4time_ob_a(maxobs), before, after

        character  stname(maxobs)*5
	character  timech*9, time*4
        character  provider(maxobs)*11
        character  weather(maxobs)*25
        character  reptype_l(maxobs)*6, atype_l(maxobs)*6
        character*9 a9time_before, a9time_after, a9time_a(maxobs)
        logical l_dupe(maxobs)
c
c.....  Output arrays
c
        real store_1(maxsta,4), 
     &         store_2(maxsta,3), store_2ea(maxsta,3),
     &         store_3(maxsta,4), store_3ea(maxsta,2),
     &         store_4(maxsta,5), store_4ea(maxsta,2),
     &         store_5(maxsta,4), store_5ea(maxsta,4),
     &         store_6(maxsta,5), store_6ea(maxsta,2),
     &         store_7(maxsta,3),
     &         store_cldht(maxsta,5)

        character  stations(maxsta)*20
        character  reptype(maxsta)*6, atype(maxsta)*6
        character  store_cldamt(maxsta,5)*4
c
c.....  Start.
c
        write(6,*)' Subroutine get_local_cwb'
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

c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing / 1000.) !in km
c
c.....	Zero out the counters.
c
        n_obs_g = 0	        ! # of local obs in the laps grid
        n_obs_b = 0	        ! # of local obs in the box
c
c.....  Get the mesonet data.
c
        ix = 1
c
c.....  Set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after

!       Ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 0 ! 3599
c
c modified by min-ken.hsieh
c change file time interval to 900(15min)
c
c
        call get_filetime_range(before,after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,900                                     
     1                         ,i4time_file_b,i4time_file_a)              

c
c modified by min-ken.hsieh
c change file time interval to 900(15min)
c temporally close this loop
c so only the exact time-matched data file can be read
c meanwhile itime_before/itime_after in obs_driver.nl will not be used here
c
c       do i4time_file = i4time_file_a, i4time_file_b, -900
            i4time_file = i4time_sys
            maxobs_in = maxobs-ix+1

            write(6,*)
     1           ' maxobs/maxobs_in/i4time_file (read_local_cwb call)'       
     1            ,maxobs,maxobs_in,i4time_file

!           Note that 'ix' subscripts can be added in if we want to
!           merge data from various filetimes. We currently read with
!           just one filetime as that should cover the observation
!           time window reasonably well for now.

            call read_local_cwb ( path_to_local_data, maxobs_in,
     ~                      badflag, ibadflag, i4time_file, 
     ~                      reptype_l, atype_l, stname, lats, lons,       
     ~                      elev, t, t24max, t24min, td, rh, 
     ~                      pcp1hr, pcp3hr, pcp6hr, pcp24hr,
     ~                      dd, ff, wgdd, wgff, 
     ~                      stnp, mslp, pcc, pc, sr, st,
     ~                      num, istatus )

c           call read_local_cwb(path_to_local_data,maxobs_in
c    1                 ,badflag,ibadflag,i4time_file                     ! I
c    1                 ,stname(ix)                                       ! O
c    1                 ,lats(ix),lons(ix),elev(ix)                       ! O
c    1                 ,i4time_ob_a(ix),t(ix),td(ix),rh(ix)              ! O
c    1                 ,pcp(ix),stnp(ix),mslp(ix),dd(ix),ff(ix)          ! O
c    1                 ,num,istatus)                                     ! O

	    if(istatus .ne. 1)then
                write(6,*)
     1          '     Warning: bad status return from READ_LOCAL_CWB'       
                n_local_time = num
            else
                write(6,*)' Good status return from READ_LOCAL_CWB'
                n_local_time = num
            endif

            write(6,*)'     n_local_time = ',n_local_time

c           ix = ix + n_local_time

c modified by min-ken,hsieh
c close loop
c       enddo ! i4time_file



!       This might make the lines shorter when reading with 'vi'
c       n_local_all = ix - 1
        n_local_all = num
        write(6,*)' n_local_all = ',n_local_all

        max_write = 100
c
c
c..................................
c.....	First QC loop over all the obs.
c..................................
c
	do i=1,n_local_all
           l_dupe(i) = .false.
c
c........  Toss the ob if lat/lon/elev or observation time are bad by setting 
c........  lat to badflag (-99.9), which causes the bounds check to think that
c........  the ob is outside the LAPS domain.
	   if( nanf( lats(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( lons(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( elev(i) ) .eq. 1 ) lats(i)  = badflag

           i4time_ob_a(i) = i4time_sys

	   call make_fnam_lp(i4time_ob_a(i),a9time_a(i),istatus)

           call filter_string(stname(i))

           do k = 1,i-1
             if(       stname(i) .eq. stname(k) 
     1                          .AND.
     1           ( (.not. l_dupe(i)) .and. (.not. l_dupe(k)) )
     1                                                           )then
                 i_diff = abs(i4time_ob_a(i) - i4time_sys)
                 k_diff = abs(i4time_ob_a(k) - i4time_sys)

                 if(i_diff .ge. k_diff)then
                     i_reject = i
                 else
                     i_reject = k
                 endif

                 if(i .le. 100 .or. i .eq. (i/10)*10)then
                     write(6,51)i,k,stname(i),a9time_a(i),a9time_a(k)     
     1                         ,i_reject
51                   format(' Duplicate detected ',2i6,1x,a6,1x,a9
     1                     ,1x,a9,1x,i6)
                 endif

                 lats(i_reject) = badflag ! test with this for now

                 l_dupe(i_reject) = .true.
             endif
           enddo ! k
c
c
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( td(i)   ) .eq. 1 ) td(i)    = badflag
	   if( nanf( dd(i)   ) .eq. 1 ) dd(i)    = badflag
	   if( nanf( ff(i)   ) .eq. 1 ) ff(i)    = badflag
c
	enddo !i
c
c..................................
c.....	Second QC loop over all the obs.
c..................................
c
	jfirst = 1
	box_low = 1. - float(ibox_points)  !buffer on west/south side
	box_idir = float(ni + ibox_points) !buffer on east
	box_jdir = float(nj + ibox_points) !buffer on north

	do i=1,n_local_all
	   if(abs(lats(i)) .gt. +90.) go to 125	
	   call latlon_to_rlapsgrid(lats(i),lons(i),lat,lon,
     &                              ni,nj,ri_loc,rj_loc,istatus)
	   if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir) go to 125
	   if(rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) go to 125
c
c.....  Elevation ok?
c
	   if(elev(i) .eq. badflag) go to 125
	   if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  If you want to check the valid time, or if there are more than
c.....  one report from this station, put that here.
c
c
c.....  Check to see if its in the desired time window.
c
	   if(i4time_ob_a(i) .lt. before 
     1   .or. i4time_ob_a(i) .gt. after) then
               write(6,91,err=125)i,stname(i)
     1                               ,a9time_a(i),before
     1                               ,after
 91            format(i6,1x,a8,' out of time ',a11,2i12)
               go to 125
           endif
c
c.....  Right time, right location...

           timech = a9time_a(i)
	   time = timech(6:9)
	   read(time,*) rtime
c          	
	   nn = nn + 1

           if(nn .gt. maxsta)then
              write(6,*)' ERROR in get_local_obs: increase maxsta '
     1                 ,nn,maxsta
              stop
           endif

	   n_obs_b = n_obs_b + 1
c
c.....  Check if its in the LAPS grid.
c
           if(ri_loc.lt.1 .or. ri_loc.gt.float(ni)) go to 151  !off grid
           if(rj_loc.lt.1 .or. rj_loc.gt.float(nj)) go to 151  !off grid
           n_obs_g = n_obs_g + 1                           !on grid...count it
151        continue
c
c.....  Fill expected accuracy arrays...see the 'get_metar_obs' routine for details.
c.....  Note that these values are only guesses based on US mesonet stations.
c
           store_2ea(nn,1) = 3.0             ! temperature (deg F)               
           store_2ea(nn,2) = 3.0             ! Dew point (deg F)               
           store_2ea(nn,3) = 30.0            ! Relative Humidity %
c
c..... Wind direction (deg) and speed (kts)
c
           store_3ea(nn,1) = 15.0            ! wind direction (dir)
           store_3ea(nn,2) = 6.0             ! wind speed (kt)

c..... Pressure and altimeter (mb)
c
           store_4ea(nn,1) = 2.00            ! pressure (mb)
           store_4ea(nn,2) = 0.00            ! altimeter (mb) (don't have)
c
c..... Other stuff (don't report these). 
c 
           store_5ea(nn,1) = 0.0             ! Visibility 
           store_5ea(nn,2) = 0.0             ! solar radiation       
           store_5ea(nn,3) = 0.0             ! soil/water temperature
           store_5ea(nn,4) = 0.0             ! soil moisture
c
           store_6ea(nn,1) = 0.0             ! precipitation (in)
           store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c.....  Clouds get set to zero, since don't have cloud info from these mesonets.
c
	   kkk = 0
c
c.....  Output the data to the storage arrays.
c
!    	  call s_len(stname(i), len)
!         stations(nn)(1:len) = stname(i)(1:len)

 	  call s_len(stname(i), len)
          if(len .ne. 0)then
              stations(nn)(1:len) = stname(i)(1:len) ! station name
          else
              write(6,*)' Warning in get_local_cwb: blank station name.'
     1                 ,' Assigning name ',i
              write(stations(nn),101)i
 101	      format(i5,15x)
          endif
c
          atype(nn) = atype_l(i)
c
          reptype(nn) = reptype_l(i)
c
	  weather(nn)(1:25) = 'UNK                     '
	  provider(nn)(1:11) = 'CWB        '
	  wmoid(nn) = ibadflag
c 
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation (m)
	 store_1(nn,4) = rtime                  ! observation time
c	
	 store_2(nn,1) = t(i)                   ! temperature (deg F)
	 store_2(nn,2) = td(i)                  ! dew point (deg F)
	 store_2(nn,3) = rh(i)                  ! Relative Humidity
c
         store_3(nn,1) = dd(i)                  ! wind dir (deg)
         store_3(nn,2) = ff(i)                  ! wind speed (kt)
         store_3(nn,3) = wgdd(i)                ! wind gust dir (deg)
         store_3(nn,4) = wgff(i)                ! wind gust speed (kt)
c
         store_4(nn,1) = badflag                ! altimeter setting (mb)
         store_4(nn,2) = stnp(i)                ! station pressure (mb)
         store_4(nn,3) = mslp(i)                ! MSL pressure (mb)
         store_4(nn,4) = pcc(i)                 ! 3-h press change character
         store_4(nn,5) = pc(i)                  ! 3-h press change (mb)
c
         store_5(nn,1) = badflag                ! visibility (miles)
         store_5(nn,2) = sr(i)                  ! solar radiation 
         store_5(nn,3) = st(i)                  ! soil/water temperature
         store_5(nn,4) = badflag                ! soil moisture 
c
         store_6(nn,1) = pcp1hr(i)              ! 1-h precipitation
         store_6(nn,2) = pcp3hr(i)              ! 3-h precipitation
         store_6(nn,3) = pcp6hr(i)              ! 6-h precipitation
         store_6(nn,4) = pcp24hr(i)             ! 24-h precipitation
         store_6(nn,5) = badflag                ! snow cover
c
         store_7(nn,1) = float(kkk)             ! number of cloud layers
         store_7(nn,2) = t24max(i)              ! 24-h max temperature
         store_7(nn,3) = t24min(i)              ! 24-h min temperature
c
c.....  That's it for this station.
c
 125     continue
       enddo !i
c
c.....  That's it...lets go home.
c
	 print *,' Found ',n_obs_b,' local obs in the LAPS box'
	 print *,' Found ',n_obs_g,' local obs in the LAPS grid'
         print *,' '
         jstatus = 1            ! everything's ok...
         return
c
 990     continue               ! no data available
         jstatus = 0
         print *,' WARNING: No data available from GET_LOCAL_CWB'
         return
c
         end

