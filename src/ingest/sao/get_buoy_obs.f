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
c
c
	subroutine get_buoy_obs(maxobs,maxsta,i4time,data_file,
     &                      buoy_format,
     &                      eastg,westg,anorthg,southg,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_buoy_g,n_buoy_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)

c
c*****************************************************************************
c
c	Routine to gather data from the buoy files for LAPS.   
c
c	Changes:
c		P. Stamus  08-10-98  Original version (from get_metar_obs).
c                          06-21-99  Change ob location check to gridpt space.
c                                      Figure box size in gridpoint space from
c                                      user-defined size (deg) and grid_spacing.
c                          10-19-99  Added checks on each variable when doing
c                                      units conversion.
c                          01-11-00  Fixed check on ob time (overall).
c
c       Notes:
c         1. This routine is not set up to collect cloud data (from ship
c            reports).  It could be done at some point if the LAPS cloud
c            analysis could make use of the cloud ob as its reported (limited
c            detail).
c
c*****************************************************************************
c
	include 'netcdf.inc'
	include 'surface_obs.inc'
c
c.....  Read arrays.
c
        integer maxobs,maxsta
	real*8  timeobs(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs)
	real*4  dd(maxobs), ff(maxobs), ffg(maxobs)
	real*4  mslp(maxobs)
	real*4  vis(maxobs), dp(maxobs)
	real*4  pcp1(maxobs), pcp6(maxobs), pcp24(maxobs)
	real*4  equivspd(maxobs), sea_temp(maxobs), t_wet(maxobs)
        real    lat(ni,nj), lon(ni,nj)
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
	integer*4  itime60, before, after, wmoid(maxobs)
	integer    rtime, dpchar(maxobs), iplat_type(maxobs)
	integer    recNum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*8, save_stn(maxobs)*8
	character  data_file*(*), timech*9, time*4, buoy_format*(*)
	character  stations(maxsta)*20, provider(maxsta)*11
	character  weather(maxobs)*25, wx(maxsta)*25
	character  reptype(maxobs)*6, atype(maxobs)*6
	character  store_cldamt(maxsta,5)*4
c
c
c.....  Start.
c
	ibadflag = int( badflag )
c
c.....	Set jstatus flag for the buoy data to bad until we find otherwise.
c
	jstatus = -1
c
c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing / 1000.)  !in km
c
c.....	Zero out the counters.
c
	n_buoy_g = 0		! # of buoy obs in the laps grid
	n_buoy_b = 0		! # of buoy obs in the box

        call s_len(buoy_format, len_buoy_format)
        if(buoy_format(1:len_buoy_format) .eq. 'FSL')then ! FSL NetCDF format
c
c.....      Get the data from the NetCDF file.  First, open the file.
c.....      If not there, return.
c
	    nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)

	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *, data_file
               go to 990
	    endif
c
c.....      Get the dimension of some of the variables.
c
c.....      "recNum"
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

	    call read_buoy(nf_fid , recNum, iplat_type,
     &         td, elev, equivspd, lats, lons, 
     &         pcp1, pcp24, pcp6,
     &         wx, dp, dpchar,
     &         mslp, sea_temp, stname, t,
     &         timeobs, vis, t_wet, dd, ffg, ff,
     &         badflag, istatus)
c
        else ! Read buoy obs in CWB format
            recNum = 150

	    call read_buoy_cwb(data_file, recNum, iplat_type,
     &         td, elev, equivspd, lats, lons, 
     &         pcp1, pcp24, pcp6,
     &         wx, dp, dpchar,
     &         mslp, sea_temp, stname, t,
     &         timeobs, vis, t_wet, dd, ffg, ff,
     &         badflag, istatus)

        endif

	if(istatus .ne. 1) go to 990
	n_buoy_all = recNum
c
c.....  First check the data coming from the NetCDF file.  There can be
c.....  "FloatInf" (used as fill value) in some of the variables.  These
c.....  are not handled the same by different operating systems.  For 
c.....  example, IBM systems make "FloatInf" into "NaN" and store them that
c.....  way in the file, which messes up other LAPS routines.  This code
c.....  checks for "FloatInf" and sets the variable to 'badflag'.  If the
c.....  "FloatInf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
	do i=1,n_buoy_all
c
c.....  Toss the ob if lat/lon/elev or observation time are bad by setting 
c.....  lat to badflag (-99.9), which causes the bounds check to think that
c.....  the ob is outside the LAPS domain.
c
	   if( nanf( lats(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( lons(i) ) .eq. 1 ) lats(i)  = badflag
	   if( nanf( elev(i) ) .eq. 1 ) lats(i)  = badflag
c
	   if( nanf( timeobs(i) ) .eq. 1 ) lats(i) = badflag
c
	   if( nanf( vis(i)  ) .eq. 1 ) vis(i)   = badflag
	   if( nanf( mslp(i) ) .eq. 1 ) mslp(i)  = badflag
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( td(i)   ) .eq. 1 ) td(i)    = badflag
	   if( nanf( dd(i)   ) .eq. 1 ) dd(i)    = badflag
	   if( nanf( ff(i)   ) .eq. 1 ) ff(i)    = badflag
	   if( nanf( ffg(i)  ) .eq. 1 ) ffg(i)   = badflag
	   if( nanf( pcp1(i) ) .eq. 1 ) pcp1(i)  = badflag
	   if( nanf( pcp6(i) ) .eq. 1 ) pcp6(i)  = badflag
	   if( nanf( pcp24(i)) .eq. 1 ) pcp24(i) = badflag
	   if( nanf( dp(i)   ) .eq. 1 ) dp(i)    = badflag
c
	   if( nanf( sea_temp(i) ) .eq. 1 ) sea_temp(i) = badflag
	   if( nanf( t_wet(i)    ) .eq. 1 ) t_wet(i)    = badflag
	   if( nanf( equivspd(i) ) .eq. 1 ) equivspd(i) = badflag
c
	enddo !i
c
c.....  Set up the time window.
c
	before = i4time - time_before
	after  = i4time + time_after
c
c..................................
c.....	Now loop over all the obs.
c..................................
c
	jfirst = 1
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north
c
	do 125 i=1,n_buoy_all
c
c.....  Bounds check: is station in the box?  Find the ob i,j location
c.....  on the LAPS grid, then check if outside past box boundary.
c
           if(lats(i) .lt. -90.) go to 125   !badflag (-99.9)...from nan ck
           call latlon_to_rlapsgrid(lats(i),lons(i),lat,lon,ni,nj,
     &                              ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir) go to 125
           if(rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) go to 125
c
c.....  Elevation ok?
c
          if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  Check to see if its in the desired time window.
c
	  itime60 = nint(timeobs(i)) + 315619200
	  if(itime60.lt.before .or. itime60.gt.after) go to 125
c
c.....  Right time, right location...

 	  call make_fnam_lp(itime60,timech,istatus)
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
	     if(stname(i) .eq. save_stn(k)) go to 125
	  enddo !k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i)  ! only one...save for checking
c
 150	  nn = nn + 1
	  n_buoy_b = n_buoy_b + 1     !station is in the box
c
c.....  Check if its in the LAPS grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_buoy_g = n_buoy_g + 1  !on grid...count it
 151	  continue
c
c.....	Figure out the cloud data.
c.....  NOTE: Not reading cloud data from ship/buoy file.  The data
c.....        is too ambiguous for LAPS use at this time.
c
	  kkk = 0               ! number of cloud layers
c
c
c.....  Convert units for storage.
c
c.....  Temperature and dewpoint
c
	temp_k = t(i)                         
	if(temp_k.lt.190. .or. temp_k.gt.345.) temp_k = badflag
	if(temp_k .le. badflag) then          !t bad?
	   temp_f = badflag                   !          bag
	else
	   temp_f = ((temp_k - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
	dewp_k = td(i)
	if(dewp_k.lt.210. .or. dewp_k.gt.320.) dewp_k = badflag
	if(dewp_k .le. badflag) then           !dp bad?
	   dewp_f = badflag                    !         bag
	else
	   dewp_f = ((dewp_k - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
c..... Wind speed and direction
c
	ddg = badflag
	if(dd(i).lt.0. .or. dd(i).gt.360.) then
	   dd(i) = badflag
	endif
c
	if(ff(i).lt.0. .or. ff(i).gt.100.) then
	   ff(i) = badflag
	else
	   ff(i)  = 1.94254 * ff(i) !m/s to kt
	endif
	if(ffg(i).lt.0. .or. ffg(i).gt.120.) then
	   ffg(i) = badflag
	else
	   ffg(i) = 1.94254 * ffg(i) !m/s to kt
	   ddg = dd(i)
	endif
c
c..... Pressure...MSL and 3-h pressure change
c
	if(mslp(i).lt.85000. .or. mslp(i).gt.120000.) then
	   mslp(i) = badflag
	else
	   mslp(i) = mslp(i) * 0.01 !Pa to mb
	endif
	if(dp(i)   .ne. badflag)   dp(i) =   dp(i) * 0.01 !Pa to mb
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
c..... Climo-type stuff...Precip, sea surface temps
c
	if(pcp1(i)  .ne. badflag)  pcp1(i) =  pcp1(i) * 39.370079 ! m to in
	if(pcp6(i)  .ne. badflag)  pcp6(i) =  pcp6(i) * 39.370079 ! m to in
	if(pcp24(i) .ne. badflag) pcp24(i) = pcp24(i) * 39.370079 ! m to in
c
	seatemp_k = sea_temp(i)                         
	if(seatemp_k .le. badflag) then          !t bad?
	   seatemp_f = badflag                   !  bag
	else
	   seatemp_f = ((seatemp_k - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
c..... Fill the expected accuracy arrays.  Values are based on information
c..... in the 'Coastal-Marine Automated Network (C-MAN) Users Guide', and 
c..... we assume that they are about the same for buoys, ships, and C-MAN
c..... sites.  Note that we convert the units to match what we're using here.
c
c..... Temperature (deg F)
c
	fon = 9. / 5.  !ratio when converting C to F
	store_2ea(nn,1) = 5.0 * fon        ! start...we don't know what we have
	if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-40.) .and. temp_f.le.c2f(50.)) then
	      store_2ea(nn,1) = 1.0 * fon  ! conv to deg F
	   endif
	endif
c
c..... Dew point (deg F).  Also estimate a RH accuracy based on the dew point.
c..... Estimates for the RH expected accuracy are from playing around with the
c..... Psychrometric Tables for various T/Td combinations.
c
	 store_2ea(nn,2) = 5.0 * fon       ! start...don't know what we have 
	 store_2ea(nn,3) = 50.0            ! Relative Humidity %
	 if(dewp_f .ne. badflag) then
	    if(dewp_f.ge.c2f(-35.) .and. dewp_f.le.c2f(-2.)) then
	       store_2ea(nn,2) = 2.0 * fon ! conv to deg F
	       store_2ea(nn,3) = 20.0      ! RH (%) 
	    elseif(dewp_f.gt.c2f(-2.) .and. dewp_f.le.c2f(30.)) then
	       store_2ea(nn,2) = 1.0 * fon ! conv to deg F
	       store_2ea(nn,3) = 8.0       ! RH (%) 
	    endif
	 endif
c
c..... Wind direction (deg) and speed (kts)
c
	 store_3ea(nn,1) = 15.0    ! deg 
	 store_3ea(nn,2) =  1.0    ! kt
	 if(ff(i) .ne. badflag) then
	    if(ff(i).ge.1.0 .and. ff(i).le.40.0) then
	       store_3ea(nn,2) = 2.0          ! kt
	    elseif(ff(i) .gt. 40.0) then
	       store_3ea(nn,2) = ff(i) * 0.05  ! 5% of speed (kts)
	    endif
c
	    if(ff(i) .ge. 5.0) then    ! dir check
	       store_3ea(nn,1) = 10.0  ! deg
	    endif
	 endif
c
c..... Pressure and altimeter (mb)
c
	 store_4ea(nn,1) = 1.00            ! pressure (mb)
	 store_4ea(nn,2) = 0.00            ! altimeter (mb)
c
c..... Visibility (miles).  Use a guess based on the range between 
c..... reportable values (e.g., for reported visibility between 
c..... 0 and 3/8th mile, set accuracy to 1/16th mile).  This is close
c..... to the stated accuracies in the C-MAN user guide.
c
	 store_5ea(nn,1) = 10.00         ! Start with this (miles)
	 if(vis(i) .ne. badflag) then
	    if(vis(i) .le. 0.375) then
	       store_5ea(nn,1) = 0.0625	! miles
	    elseif(vis(i).gt.0.375 .and. vis(i).le.2.0) then
	       store_5ea(nn,1) = 0.125 ! miles
	    elseif(vis(i).gt.2.0 .and. vis(i).le.3.0) then
	       store_5ea(nn,1) = 0.25 ! miles
	    elseif(vis(i).gt.3.0 .and. vis(i).le.8.0) then
	       store_5ea(nn,1) = 1.00 ! miles
	    elseif(vis(i) .gt. 8.0) then
	       store_5ea(nn,1) = 5.00 ! miles
	    endif
	 endif
c
c..... Other stuff.  Note that precip isn't very good.
c
	 store_5ea(nn,2) = 0.0             ! solar radiation 
	 store_5ea(nn,3) = 1.0 * fon       ! soil/water temperature (F)
	 store_5ea(nn,4) = 0.0             ! soil moisture
c
	 store_6ea(nn,1) = 0.20            ! precipitation (in)
	 store_6ea(nn,2) = 0.0             ! snow cover (in) 
c
c
c..... Output the data to the storage arrays
c
	 call s_len(stname(i), len)
	 stations(nn)(1:len) = stname(i)(1:len) ! station name
	 provider(nn)(1:11) = 'NWS        '     ! data provider (all from NWS)
	 reptype(nn)(1:6) = 'MARTIM'            ! station type
	 atype(nn)(1:6) ='      '               ! used here for moving/fixed stns
	 if(iplat_type(i) .eq. 0) then
	    atype(nn)(1:6) = 'FIX   '
	 elseif(iplat_type(i) .eq. 1) then
	    atype(nn)(1:6) = 'MVG   '
	 endif
	 wmoid(nn) = ibadflag                   ! WMO id...not applicable here
	 weather(nn)(1:25) = wx(i)(1:25)        ! present weather
         call filter_string(weather(nn))
c       
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = dewp_f                 ! dew point (deg f)
	 store_2(nn,3) = badflag                ! Relative Humidity
c
	 store_3(nn,1) = dd(i)                  ! wind dir (deg)
	 store_3(nn,2) = ff(i)                  ! wind speed (kt)
	 store_3(nn,3) = ddg                    ! wind gust dir (deg)
	 store_3(nn,4) = ffg(i)                 ! wind gust speed (kt)
c
	 store_4(nn,1) = badflag                ! altimeter setting (mb)
	 store_4(nn,2) = badflag                ! station pressure (mb)
	 store_4(nn,3) = mslp(i)                ! MSL pressure (mb)

         if(dpchar(i) .ne. ibadflag)then        ! 3-h press change character
  	     store_4(nn,4) = float(dpchar(i))       
         else
  	     store_4(nn,4) = badflag
         endif

         store_4(nn,5) = dp(i)                  ! 3-h press change (mb)
c
	 store_5(nn,1) = vis(i)                 ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = seatemp_f              ! soil/water temperature (F)
	 store_5(nn,4) = badflag                ! soil moisture 
c
	 store_6(nn,1) = pcp1(i)                ! 1-h precipitation (in)
	 store_6(nn,2) = badflag                ! 3-h precipitation
	 store_6(nn,3) = pcp6(i)                ! 6-h precipitation
	 store_6(nn,4) = pcp24(i)               ! 24-h precipitation
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
	 print *,' Found ',n_buoy_b,' buoy obs in the LAPS box'
	 print *,' Found ',n_buoy_g,' buoy obs in the LAPS grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' WARNING.  No data available from READ_BUOY.'
	 return
c
	 end
