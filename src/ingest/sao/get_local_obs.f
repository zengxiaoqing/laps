c
        subroutine get_local_obs(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 eastg,westg,anorthg,southg,
     &                 lat,lon,ni,nj,grid_spacing,
     &                 nn,n_local_g,n_local_b,stations,
     &                 reptype,atype,weather,wmoid,
     &                 store_1,store_2,store_2ea,
     &                 store_3,store_3ea,store_4,store_4ea,
     &                 store_5,store_5ea,store_6,store_6ea,
     &                 store_7,store_cldht,store_cldamt,
     &                 provider, laps_cycle_time, jstatus)

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
	include 'netcdf.inc'
c
c.....  Read arrays.
c
        integer maxobs ! raw data file
        integer maxsta ! processed stations for LSO file
	real*8  timeobs(maxobs), rh_time(maxobs), p_time(maxobs)
	real*8  t_time(maxobs), dd_time(maxobs), gust_time(maxobs)
	real*8  ff_time(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs), rh(maxobs), stnp(maxobs)
	real*4  dd(maxobs), ff(maxobs), ddg(maxobs), ffg(maxobs)
	real*4  mslp(maxobs), alt(maxobs), vis(maxobs)
        real    lat(ni,nj), lon(ni,nj)
        logical l_dupe(maxobs)
c
c.....  Output arrays.
c
	real*4  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
	integer*4  before, after, wmoid(maxsta)
        integer*4  i4time_ob_a(maxobs)
	integer    rtime
	integer    recNum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*6, save_stn(maxobs)*6
	character  timech*9, time*4
	character  stations(maxsta)*20
	character  pro(maxobs)*11, provider(maxsta)*11
	character  wx(maxobs)*25, weather(maxsta)*25
	character  reptype(maxsta)*6, atype(maxsta)*6
	character  store_cldamt(maxsta,5)*4, stn_type(maxobs)*11
        character*(*) path_to_local_data, local_format
        character*9 a9time_before, a9time_after, a9time_a(maxobs)
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
c
c.....  Start.
c
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
c
c.....	Zero out the counters.
c
        n_local_g = 0	        ! # of local obs in the laps grid
        n_local_b = 0	        ! # of local obs in the box
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

        do i4time_file = i4time_file_a, i4time_file_b, -3600

            call s_len(path_to_local_data,len_path)
            filename13= cvt_i4time_wfo_fname13(i4time_file)
 	    data_file = path_to_local_data(1:len_path)//filename13

            write(6,*)' mesonet file = ',data_file(1:len_path+13)

	    nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)

	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *, data_file
	       go to 990
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
c
c.....  Call the read routine.
c
	    call read_local(nf_fid, recNum, alt(ix),
     &         pro(ix), td(ix), elev(ix), lats(ix), lons(ix),
     &         timeobs(ix), wx(ix), rh(ix), rh_time(ix),
     &         mslp(ix), stname(ix), p_time(ix), stnp(ix), 
     &         stn_type(ix), t_time(ix), t(ix), vis(ix),
     &         dd(ix), dd_time(ix), ddg(ix), ffg(ix), gust_time(ix), 
     &         ff(ix), ff_time(ix), badflag, istatus)
c
	    if(istatus .ne. 1)then
                write(6,*)
     1          '     Warning: bad status return from READ_LOCAL'       
                n_local_file = 0

            else
                n_local_file = recNum
                write(6,*)'     n_local_file = ',n_local_file

            endif

            ix = ix + n_local_file

        enddo ! i4time_file

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
c.....	First loop over all the obs.
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
	   if( nanf( timeobs(i) ) .eq. 1 ) lats(i) = badflag

	   i4time_ob_a(i) = nint(timeobs(i)) + 315619200
	   call make_fnam_lp(i4time_ob_a(i),a9time_a(i),istatus)

           call filter_string(stname(i))

           do k = 1,i-1
             if(       stname(i) .eq. stname(k) 
!    1                          .AND.
!    1           (lats(i) .ne. badflag .and. lats(k) .ne. badflag)
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

                 write(6,51)i,k,stname(i),a9time_a(i),a9time_a(k)
     1                     ,i_reject
 51		 format(' Duplicate detected ',2i6,1x,a6,1x,a9,1x,a9
     1                 ,1x,i6)

                 lats(i_reject) = badflag ! test with this for now

                 l_dupe(i_reject) = .true.
             endif
           enddo ! k
c
c
	   if( nanf( rh_time(i)   ) .eq. 1 ) rh_time(i)   = ibadflag
	   if( nanf( t_time(i)    ) .eq. 1 ) t_time(i)    = ibadflag
	   if( nanf( p_time(i)    ) .eq. 1 ) p_time(i)    = ibadflag
	   if( nanf( dd_time(i)   ) .eq. 1 ) dd_time(i)   = ibadflag
	   if( nanf( ff_time(i)   ) .eq. 1 ) ff_time(i)   = ibadflag
	   if( nanf( gust_time(i) ) .eq. 1 ) gust_time(i) = ibadflag
c
	   if( nanf( vis(i)  ) .eq. 1 ) vis(i)   = badflag
	   if( nanf( mslp(i) ) .eq. 1 ) mslp(i)  = badflag
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( td(i)   ) .eq. 1 ) td(i)    = badflag
	   if( nanf( dd(i)   ) .eq. 1 ) dd(i)    = badflag
	   if( nanf( ff(i)   ) .eq. 1 ) ff(i)    = badflag
	   if( nanf( ffg(i)  ) .eq. 1 ) ffg(i)   = badflag
	   if( nanf( alt(i)  ) .eq. 1 ) alt(i)   = badflag
c
	enddo !i
c
c..................................
c.....	Second loop over all the obs.
c..................................
c
	jfirst = 1
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north
c
	do 125 i=1,n_local_all
c
c.....  Bounds check: is station in the box?  Find the ob i,j location
c.....  on the LAPS grid, then check if outside past box boundary.
c
!          Test for badflag OR (close to S Pole but not quite at it)
!          Check can also be generalized in 'latlon_to_rlapsgrid' for 'lambert'
           if(lats(i) .lt. -89.999 .and. lats(i) .ne. -90.) go to 125 
           call latlon_to_rlapsgrid(lats(i),lons(i),lat,lon,ni,nj,       
     &                              ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir
     1   .or. rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) then
               if(i .le. max_write)then
                   write(6,81,err=125)i,wmoid(i),stname(i)
     1                               ,nint(ri_loc),nint(rj_loc)
 81                format(i6,i7,1x,a8,' out of box ',2i12)
               endif
               go to 125
           endif
c
c.....  Elevation ok?
c
	   if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125       
c
c.....  Check to see if its in the desired time window.
c
	   if(i4time_ob_a(i) .lt. before 
     1   .or. i4time_ob_a(i) .gt. after) then
               if(i .le. max_write)then
                   write(6,91,err=125)i,wmoid(i),stname(i)
     1                               ,a9time_a(i),before
     1                               ,after
 91		   format(i6,i7,1x,a8,' out of time ',a11,2i12)
               endif
               go to 125
           endif
c
c.....  Right time, right location...

           timech = a9time_a(i)
	   time = timech(6:9)
	   read(time,*) rtime
c
c.....  Check if station is reported more than once this
c.....  time period.
c
	   if(jfirst .eq. 1) then
	     icount = 1
	     save_stn(1) = stname(i)
	     jfirst = 0
	     go to 150
	  endif
c
	  do k=1,icount
             if(stname(i) .eq. save_stn(k)) then
                 write(6,*)' Rejecting duplicate ',i,k,stname(i)
     1                    ,' ',a9time_a(i),' ',a9time_a(k)
                 go to 125
             endif
	  enddo !k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i)  ! only one...save for checking
c
 150	  nn = nn + 1

          if(nn .gt. maxsta)then
              write(6,*)' ERROR in get_local_obs: increase maxsta '
     1                 ,nn,maxsta
              stop
          endif
 
	  n_local_b = n_local_b + 1     !station is in the box
c
c.....  Check if its in the LAPS grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151 !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151 !off grid
	  n_local_g = n_local_g + 1  !on grid...count it
 151	  continue
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
	  temp_k = t(i) 
	  if(temp_k.lt.190. .or. temp_k.gt.345.) temp_k = badflag
	  if(t_time(i) .ge. 0.) then ! implies that it is not set to ibadflag
	     if( (timeobs(i) - t_time(i)) .gt. laps_cycle_time) then
		temp_k = badflag
	     endif
	  endif
	  if(temp_k .le. badflag) then !t bad?
	     temp_f = badflag	!then bag it
	  else
	     temp_f = ((temp_k - 273.16) * 9./5.) + 32. ! K to F
	  endif
c       
	  dewp_k = td(i)
	  if(dewp_k.lt.210. .or. dewp_k.gt.320.) dewp_k = badflag
	  if(dewp_k .le. badflag) then !dp bad?
	     dewp_f = badflag	       !then bag it
	  else
	     dewp_f = ((dewp_k - 273.16) * 9./5.) + 32.	! K to F
	  endif
c
	  rh_p = rh(i) 
	  if(rh_p.lt.0. .or. rh_p.gt.100.) rh_p = badflag
	  if(rh_time(i) .ge. 0.) then
	     if( (timeobs(i) - rh_time(i)) .gt. laps_cycle_time) then
		rh_p = badflag
	     endif
	  endif
c
c..... Wind speed and direction
c
	  dir = dd(i) 
	  if(dir.lt.0. .or. dir.gt.360.) dir = badflag
	  spd = ff(i)
	  if(spd.lt.0 .or. spd.gt.100.) spd = badflag
	  if(dd_time(i).ge.0. .and. ff_time(i).ge.0.) then
	     if( ((timeobs(i) - dd_time(i)) .gt. laps_cycle_time) .or.
     &           ((timeobs(i) - ff_time(i)) .gt. laps_cycle_time) ) then
		dir = badflag
		spd = badflag
	     endif
	  endif
	  if(spd .ne. badflag) spd = 1.94254 * spd !m/s to kt
c
	  dirgust = ddg(i)
	  if(dirgust.lt.0. .or. dirgust.gt.360.) dirgust = badflag
	  spdgust = ffg(i)
	  if(spdgust.lt.0 .or. spdgust.gt.120.) spdgust = badflag
	  if(gust_time(i) .ne. badflag) then
	     if( (timeobs(i) - gust_time(i)) .gt. laps_cycle_time) then
		dirgust = badflag
		spdgust = badflag
	     endif
	  endif
	  if(spdgust .ne. badflag) spdgust = 1.94254 * spdgust !m/s to kt
c
c..... Pressure...Station pressure, MSL and altimeter
c
	  stn_press = stnp(i)
	  if(stn_press.lt.40000. .or. stn_press.gt.120000.) stn_press = badflag
	  if(p_time(i) .ge. 0.) then
	     if( (timeobs(i) - p_time(i)) .gt. laps_cycle_time ) then
		stn_press = badflag
	     endif
	  endif
	  if(stn_press .ne. badflag) stn_press = stn_press * 0.01 !Pa to mb
c
	  if(mslp(i).lt.85000. .or. mslp(i).gt.120000.) then
	     mslp(i) = badflag
	  else
	     mslp(i) = mslp(i) * 0.01 !Pa to mb
	  endif
	  if(alt(i).lt.85000. .or. alt(i).gt.120000.) then
	     alt(i) = badflag
	  else
	     alt(i) =  alt(i) * 0.01 !Pa to mb
	  endif
c
c..... Visibility
c
	if(vis(i).lt.0. .or. vis(i).gt.330000.) then
	   vis(i) = badflag
	else
	   vis(i) = vis(i) * .001      !m to km
	   vis(i) = 0.621371 * vis(i)  !km to miles
	endif
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
	 if(ff(i) .ne. badflag) then
	    if(ff(i).ge.1.0 .and. ff(i).lt.10.0) then
	       store_3ea(nn,2) = 2.0          ! kt
	    elseif(ff(i) .gt. 10.0) then
	       store_3ea(nn,2) = ff(i) * 0.2  ! 20% of speed (kts)
	    endif
c
	    if(ff(i) .ge. 5.0) then    ! dir check
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
	 if(vis(i) .ne. badflag) then
	    if(vis(i) .lt. 2.0) then
	       store_5ea(nn,1) = 0.50 ! miles
	    elseif(vis(i).ge.2.0 .and. vis(i).lt.3.0) then
	       store_5ea(nn,1) = 1.00 ! miles
	    elseif(vis(i) .gt. 3.0) then
	       store_5ea(nn,1) = 2.00 ! miles
	    endif
	 endif
c
c..... Other stuff.  
c
	 store_5ea(nn,2) = 0.0             ! solar radiation 
	 store_5ea(nn,3) = 0.0             ! soil/water temperature
	 store_5ea(nn,4) = 0.0             ! soil moisture
c
	 store_6ea(nn,1) = 0.0             ! precipitation (in)
	 store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c..... Output the data to the storage arrays
c
	 call s_len(stname(i), len)
         if(len .ne. 0)then
             stations(nn)(1:len) = stname(i)(1:len) ! station name
         else
             write(6,*)' Warning in get_local_obs: blank station name.'
     1                ,' Assigning name ',i
             write(stations(nn),101)i
 101	     format(i5,15x)
         endif
c
	 call s_len(pro(i), len)
         if(len .ne. 0) then
	     provider(nn)(1:len) = pro(i)(1:len)    ! data provider
         endif
c
         call s_len(stn_type(i), len)
         if(len .ne. 0) then
            ilen = min(len, 6)
            atype(nn)(1:ilen) = stn_type(i)(1:ilen) ! auto stn type
         endif
c
         weather(nn)(1:25) = wx(i)(1:25)        ! present weather
         call filter_string(weather(nn))

	 reptype(nn)(1:6) = 'LDAD  '            ! report type
	 wmoid(nn) = ibadflag                   ! WMO ID
c
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = dewp_f                 ! dew point (deg f) 
	 store_2(nn,3) = rh_p                   ! Relative Humidity
c
	 store_3(nn,1) = dir                    ! wind dir (deg)
	 store_3(nn,2) = spd                    ! wind speed (kt)
	 store_3(nn,3) = dirgust                ! wind gust dir (deg)
	 store_3(nn,4) = spdgust                ! wind gust speed (kt)
c
	 store_4(nn,1) = alt(i)                 ! altimeter setting (mb)
	 store_4(nn,2) = stn_press              ! station pressure (mb)
	 store_4(nn,3) = mslp(i)                ! MSL pressure (mb)
	 store_4(nn,4) = badflag                ! 3-h press change character
         store_4(nn,5) = badflag                ! 3-h press change (mb)
c
	 store_5(nn,1) = vis(i)                 ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = badflag                ! soil/water temperature
	 store_5(nn,4) = badflag                ! soil moisture 
c
	 store_6(nn,1) = badflag                ! 1-h precipitation
	 store_6(nn,2) = badflag                ! 3-h precipitation
	 store_6(nn,3) = badflag                ! 6-h precipitation
	 store_6(nn,4) = badflag                ! 24-h precipitation
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
  125	 continue
c
c
c.....  That's it...lets go home.
c
	 print *,' Found ',n_local_b,' local obs in the LAPS box'
	 print *,' Found ',n_local_g,' local obs in the LAPS grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' WARNING: No data available from GET_LOCAL_OBS'
	 return
c
	 end
