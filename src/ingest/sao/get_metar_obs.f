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
	subroutine get_metar_obs(maxobs,maxsta,i4time,data_file,
     &                      metar_format,
     &                      eastg,westg,anorthg,southg,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_sao_g,n_sao_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)

c
c*****************************************************************************
c
c	Routine to gather METAR data for LAPS.   
c
c	Changes:
c		P. Stamus  10-30-96  Original version (from get_sao_obs).
c                          11-13-96  Pass in full path to data.
c                          12-06-96  Put Atype (01 or 02) in obstype var.
c	                   06-09-97  Changes for new METAR CDL.
c                          01-29-98  Upgrades for LS2.
c                          05-01-98  Add soil moisture variables.
c                          06-21-99  Change ob location check to gridpt space.
c                                      Figure box size in gridpoint space from
c                                      user-defined size (deg) and grid_spacing.
c                          07-22-99  Add elev ck for badflag.  Fix wmoid var.
c
c
c*****************************************************************************
c
	include 'netcdf.inc'
	include 'surface_obs.inc'
c
        integer maxobs,maxsta
	real*8  timeobs(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs), tt(maxobs), ttd(maxobs)
	real*4  dd(maxobs), ff(maxobs), ffg(maxobs)
	real*4  mslp(maxobs), alt(maxobs)
	real*4  ht(6,maxobs), vis(maxobs), dp(maxobs)
	real*4  pcp1(maxobs), pcp3(maxobs), pcp6(maxobs), pcp24(maxobs)
	real*4  max24t(maxobs), min24t(maxobs), snowcvr(maxobs)
        real    lat(ni,nj), lon(ni,nj)
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
	integer*4  itime60, wmoid(maxobs), wmoid_in(maxobs)
	integer*4  before, after
	integer    rtime, dpchar(maxobs)
	integer    maxSkyCover, recNum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*5, save_stn(maxobs)*5
	character  data_file*(*), timech*9, time*4, metar_format*(*)
	character  cvr(6,maxobs)*8
	character  stations(maxsta)*20, provider(maxsta)*11
	character  weather(maxobs)*25, wx(maxobs)*25
	character  reptype(maxobs)*6, atype(maxobs)*6
	character  reptype_in(maxobs)*6, atype_in(maxobs)*6
	character  store_cldamt(maxsta,5)*4
c
c
c.....	Set jstatus flag for the sao data to bad until we find otherwise.
c
	jstatus = -1
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
	n_sao_g = 0		! # of saos in the laps grid
	n_sao_b = 0		! # of saos in the box

        call s_len(metar_format,len_metar_format)

        if(metar_format(1:len_metar_format) .eq. 'FSL')then
c
c.....      Get the data from the NetCDF file.  First, open the file.
c
	    nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)

	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *, data_file
	    endif
c
c.....      Get the dimension of some of the variables.
c.....      "maxSkyCover"
c
	    nf_status = NF_INQ_DIMID(nf_fid,'maxSkyCover',nf_vid)
	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *,'dim maxSkyCover'
	    endif
	    nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,maxSkyCover)
	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *,'dim maxSkyCover'
	    endif
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

c
c.....      Call the read routine.
c
	    call read_metar(nf_fid , maxSkyCover, recNum, alt,
     &         atype_in, td, ttd, elev,
     &         lats, lons, max24t, min24t,
     &         pcp1, pcp24, pcp3, pcp6,
     &         wx, dp, dpchar,
     &         reptype_in, mslp, cvr, ht,
     &         snowcvr, stname, tt, t,
     &         timeobs, vis, dd, ffg, ff,
     &         wmoid_in, badflag, istatus)

	    if(istatus .ne. 1) go to 990
	    n_sao_all = recNum

            i4time_offset = 315619200

        else ! Read CWB Metar and Synop Obs
            recNum=150
            maxSkyCover=10
            call read_metar_cwb(data_file , maxSkyCover, recNum, alt,    
     &         atype_in, td, ttd, elev,
     &         lats, lons, max24t, min24t,
     &         pcp1, pcp24, pcp3, pcp6,
     &         wx, dp, dpchar,
     &         reptype_in, mslp, cvr, ht,
     &         snowcvr, stname, tt, t,
     &         timeobs, vis, dd, ffg, ff,
     &         wmoid_in, badflag, n_metar_cwb, istatus)

            ix = n_metar_cwb + 1

            if(.false.)then
               call read_synop_cwb(data_file , maxSkyCover, recNum, alt,       
     &         atype_in(ix), td(ix), ttd(ix), elev(ix),
     &         lats(ix), lons(ix), max24t(ix), min24t(ix),
     &         pcp1(ix), pcp24(ix), pcp3(ix), pcp6(ix),
     &         wx(ix), dp(ix), dpchar(ix),
     &         reptype_in(ix), mslp(ix), cvr(1,ix), ht(1,ix),
     &         snowcvr(ix), stname(ix), tt(ix), t(ix),
     &         timeobs(ix), vis(ix), dd(ix), ffg(ix), ff(ix),
     &         wmoid_in(ix), badflag, n_synop_cwb, istatus)
            endif

            n_sao_all = n_metar_cwb + n_synop_cwb
            if(n_sao_all .le. 0) go to 990

            i4time_offset = 0

        endif
c
	if(n_sao_all .gt. maxobs)then
            write(6,*)' ERROR, n_sao_all > maxobs ',n_sao_all,maxobs
            return
        endif
c
c.....  First check the data coming from the NetCDF file.  There can be
c.....  "FloatInf" (used as fill value) in some of the variables.  These
c.....  are not handled the same by different operating systems.  For 
c.....  example, IBM systems make "FloatInf" into "NaN" and store them that
c.....  way in the LSO, which messes up other LAPS routines.  This code
c.....  checks for "FloatInf" and sets the variable to 'badflag'.  If the
c.....  "FloatInf" is in the lat, lon, elevation, or time of observation,
c.....  we toss the whole ob since we can't be sure where it is.
c
	do i=1,n_sao_all
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
	   if( nanf( tt(i)   ) .eq. 1 ) tt(i)    = badflag
	   if( nanf( ttd(i)  ) .eq. 1 ) ttd(i)   = badflag
	   if( nanf( dd(i)   ) .eq. 1 ) dd(i)    = badflag
	   if( nanf( ff(i)   ) .eq. 1 ) ff(i)    = badflag
	   if( nanf( ffg(i)  ) .eq. 1 ) ffg(i)   = badflag
	   if( nanf( alt(i)  ) .eq. 1 ) alt(i)   = badflag
	   if( nanf( pcp1(i) ) .eq. 1 ) pcp1(i)  = badflag
	   if( nanf( pcp3(i) ) .eq. 1 ) pcp3(i)  = badflag
	   if( nanf( pcp6(i) ) .eq. 1 ) pcp6(i)  = badflag
	   if( nanf( pcp24(i)) .eq. 1 ) pcp24(i) = badflag
	   if( nanf( dp(i)   ) .eq. 1 ) dp(i)    = badflag
c
	   if( nanf( snowcvr(i) ) .eq. 1 ) snowcvr(i) = badflag
	   if( nanf( max24t(i)  ) .eq. 1 ) max24t(i)  = badflag
	   if( nanf( min24t(i)  ) .eq. 1 ) min24t(i)  = badflag
c
	   do j=1,5
	      if( nanf( ht(j,i) ) .eq. 1 ) ht(j,i) = badflag
	   enddo !j
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
	do 125 i=1,n_sao_all
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
          if(elev(i) .eq. badflag) go to 125 !from read_metar (missing elev)
          if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125
c
c.....  Check to see if its in the desired time window (if the flag
c.....  says to check this).
c
	  itime60 = nint(timeobs(i)) + i4time_offset
c
	  if(ick_METAR_time .eq. 1) then
	    if(itime60.lt.before .or. itime60.gt.after) go to 125
	  endif
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
	  enddo	!k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i) ! only one...save for checking
c
 150	  nn = nn + 1
	  n_sao_b = n_sao_b + 1	!station is in the box
c
c.....  Check if its in the LAPS grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_sao_g = n_sao_g + 1  !on grid...count it
 151	  continue
c
c.....	Figure out the cloud data.
c
	  k_layers = 0               ! number of cloud layers
c
	  if(cvr(1,i)(1:1) .eq. ' ') then
	     k_layers = 0
	  else
	     do k=1,5
		if(cvr(k,i)(1:1) .ne. ' ') k_layers = k_layers + 1       
	     enddo !k
	  endif
c
	  if(k_layers .eq. 0) then   ! no cloud data...probably AMOS station
	    go to 126		     ! skip rest of cloud stuff
	  endif

          iihigh = k_layers

	  do ii=1,iihigh
	    if(ht(ii,i) .gt. 25000.0) ht(ii,i) = badflag
c
	    if(cvr(ii,i)(1:3) .eq. 'SKC') then
	       ht(ii,i) = 22500.0    ! Manual Ob
	    elseif(cvr(ii,i)(1:3) .eq. 'CLR') then
	       ht(ii,i) = 3657.4     ! Automatic Ob
            else                     ! Check for bad height
               if(ht(ii,i) .gt. 17000.0) then
                   write(6,*)' WARNING in get_metar_obs: '      
     1                      ,' reject cloud ob, height = '
     1                      ,ht(ii,i),stname(i)
                   ht(ii,i) = badflag
                   k_layers = 0
               endif
	    endif

	  enddo !ii

 126	  continue
c
c.....	Check cloud info for very high heights...set to max if greater.
c.....	Also convert agl cloud heights to msl by adding elevation.
c
	if(k_layers .gt. 0) then
	  do ii=1,k_layers
	    if(ht(ii,i) .ge. 0.) then
	      ht(ii,i) = ht(ii,i) + elev(i)		! conv agl to msl
	      if(ht(ii,i) .gt. 22500.) ht(ii,i) = 22500.
	    endif
	  enddo !ii
	endif
c
c
c.....  Convert units for storage.
c
c.....  Temperature and dewpoint
c
	temp_k = tt(i)                        !set to temp_from_tenths
	if(temp_k .eq. badflag) temp_k = t(i) !no temp_from_tenths, set to t
	if(temp_k .eq. badflag) then          !t bad too?
	   temp_f = badflag                   !          bag
	else
	   temp_f = ((temp_k - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
	dewp_k = ttd(i)                        !set to dp_from_tenths
	if(dewp_k .eq. badflag) dewp_k = td(i) !no dp_from_tenths, set to dp
	if(dewp_k .eq. badflag) then           !dp bad too?
	   dewp_f = badflag                    !         bag
	else
	   dewp_f = ((dewp_k - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
c..... Wind speed and direction
c
	ddg = badflag
	if(ff(i)  .ne. badflag) ff(i)  = 1.94254 * ff(i)   !m/s to kt
	if(ffg(i) .ne. badflag) then
	   ffg(i) = 1.94254 * ffg(i) !m/s to kt
	   ddg = dd(i)
	endif
c
c..... Pressure...MSL and altimeter, 3-h pressure change
c
	if(alt(i)  .ne. badflag)  alt(i) =  alt(i) * 0.01   !Pa to mb
	if(mslp(i) .ne. badflag) mslp(i) = mslp(i) * 0.01   !Pa to mb
	if(dp(i)   .ne. badflag)   dp(i) =   dp(i) * 0.01   !Pa to mb
c
c..... Visibility
c
	if(vis(i) .ne. badflag) then
	   vis(i) = vis(i) * .001      !m to km
	   vis(i) = 0.621371 * vis(i)  !km to miles
	endif
c
c..... Climo-type stuff...Precip and snow cover, Max/min temps
c
	if(pcp1(i)  .ne. badflag)  pcp1(i) =  pcp1(i) * 39.370079 ! m to in
	if(pcp3(i)  .ne. badflag)  pcp3(i) =  pcp3(i) * 39.370079 ! m to in
	if(pcp6(i)  .ne. badflag)  pcp6(i) =  pcp6(i) * 39.370079 ! m to in
	if(pcp24(i) .ne. badflag) pcp24(i) = pcp24(i) * 39.370079 ! m to in
	if(snowcvr(i) .ne. badflag) snowcvr(i) = snowcvr(i) * 39.370079 ! m to in
	if(max24t(i) .ne. badflag) then
	   max24t(i) = ((max24t(i) - 273.16) * 9./5.) + 32. ! K to F
	endif
	if(min24t(i) .ne. badflag) then
	   min24t(i) = ((min24t(i) - 273.16) * 9./5.) + 32.  ! K to F
	endif
c
c
c..... Fill the expected accuracy arrays.  Values are based on information
c..... in the 'Federal Meteorological Handbook No. 1' for the METARs, 
c..... Appendix C (http://www.nws.noaa.gov/oso/oso1/oso12/fmh1/fmh1appc.htm)
c..... Note that we convert the units in Appendix C to match what we're 
c..... using here.
c
c..... Temperature (deg F)
c
	fon = 9. / 5.  !ratio when converting C to F
	store_2ea(nn,1) = 5.0 * fon        ! start...we don't know what we have
	if(temp_f .ne. badflag) then
	   if(temp_f.ge.c2f(-62.) .and. temp_f.le.c2f(-50.)) then
	      store_2ea(nn,1) = 1.1 * fon  ! conv to deg F
	   elseif(temp_f.gt.c2f(-50.) .and. temp_f.lt.c2f(50.)) then
	      store_2ea(nn,1) = 0.6 * fon  ! conv to deg F
	   elseif(temp_f.ge.c2f(50.) .and. temp_f.le.c2f(54.)) then
	      store_2ea(nn,1) = 1.1 * fon  ! conv to deg F
	   endif
	endif
c
c..... Dew point (deg F).  Also estimate a RH accuracy based on the dew point.
c..... Estimates for the RH expected accuracy are from playing around with the
c..... Psychrometric Tables for various T/Td combinations (including their
c..... accuracies from the FMH-1 Appendix C).
c
	 store_2ea(nn,2) = 5.0 * fon       ! start...don't know what we have 
	 store_2ea(nn,3) = 50.0            ! Relative Humidity %
	 if(dewp_f .ne. badflag) then
	    if(dewp_f.ge.c2f(-34.) .and. dewp_f.lt.c2f(-24.)) then
	       store_2ea(nn,2) = 2.2 * fon ! conv to deg F
	       store_2ea(nn,3) = 20.0      ! RH (%) 
	    elseif(dewp_f.ge.c2f(-24.) .and. dewp_f.lt.c2f(-1.)) then
	       store_2ea(nn,2) = 1.7 * fon ! conv to deg F
	       store_2ea(nn,3) = 12.0      ! RH (%) 
	    elseif(dewp_f.ge.c2f(-1.) .and. dewp_f.le.c2f(30.)) then
	       store_2ea(nn,2) = 1.1 * fon ! conv to deg F
	       store_2ea(nn,3) = 8.0       ! RH (%) 
	    endif
	 endif
c
c..... Wind direction (deg) and speed (kts)
c
	 store_3ea(nn,1) = 10.0    ! deg 
	 store_3ea(nn,2) =  1.0    ! kt
	 if(ff(i) .ne. badflag) then
	    if(ff(i).ge.1.0 .and. ff(i).le.10.0) then
	       store_3ea(nn,2) = 1.0          ! kt
	    elseif(ff(i) .gt. 10.0) then
	       store_3ea(nn,2) = ff(i) * 0.1  ! 10% of speed (kts)
	    endif
c
	    if(ff(i) .ge. 5.0) then    ! dir check
	       store_3ea(nn,1) = 5.0   ! deg
	    endif
	 endif
c
c..... Pressure and altimeter (mb)
c
	 store_4ea(nn,1) = 0.68            ! pressure (mb)
	 store_4ea(nn,2) = 0.68            ! altimeter (mb)
c
c..... Visibility (miles).  For automated stations use a guess based 
c..... on Table C-2 in Appendix C of FMH-1.  For manual stations, use
c..... a guess based on the range between reportable values (e.g., for
c..... reported visibility between 0 and 3/8th mile, set accuracy to 
c..... 1/16th mile).  This isn't ideal, but its a start.
c
	 store_5ea(nn,1) = 10.00         ! Start with this (miles)
	 if(vis(i) .ne. badflag) then
	    if(atype_in(i)(1:2) .eq. 'A0') then   ! have an auto station
	       if(vis(i) .lt. 2.0) then
		  store_5ea(nn,1) = 0.25         ! miles
	       elseif(vis(i).ge.2.0 .and. vis(i).lt.3.0) then
		  store_5ea(nn,1) = 0.50         ! miles
	       elseif(vis(i) .gt. 3.0) then
		  store_5ea(nn,1) = 1.00         ! miles
	       endif
	    else		! have a manual station
	       if(vis(i) .le. 0.375) then
		  store_5ea(nn,1) = 0.0625       ! miles
	       elseif(vis(i).gt.0.375 .and. vis(i).le.2.0) then
		  store_5ea(nn,1) = 0.125        ! miles
	       elseif(vis(i).gt.2.0 .and. vis(i).le.3.0) then
		  store_5ea(nn,1) = 0.25         ! miles
	       elseif(vis(i).gt.3.0 .and. vis(i).le.15.0) then
		  store_5ea(nn,1) = 1.00         ! miles
	       elseif(vis(i) .gt. 15.0) then
		  store_5ea(nn,1) = 5.00         ! miles
	       endif
	    endif
	 endif
c
c..... Other stuff.  Don't really know about the precip, but probably
c..... worse that this guess.
c
	 store_5ea(nn,2) = 0.0             ! solar radiation 
	 store_5ea(nn,3) = 0.0             ! soil/water temperature
	 store_5ea(nn,4) = 0.0             ! soil moisture 
c
	 store_6ea(nn,1) = 0.01            ! precipitation (in)
	 store_6ea(nn,2) = 1.0             ! snow cover (in) 
c
c
c..... Output the data to the storage arrays
c
	 call s_len(stname(i), len)
	 stations(nn)(1:len) = stname(i)(1:len) ! station name
c
	 call s_len(atype_in(i), len)
	 if(len .ne. 0) then
	    atype(nn)(1:len) = atype_in(i)(1:len) ! auto stn type
	 endif
c
	 call s_len(reptype_in(i), len)
	 if(len .ne. 0) then
	    reptype(nn)(1:len) = reptype_in(i)(1:len) ! report type
	 endif
c
 	 weather(nn)(1:25) = wx(i)(1:25)        ! present weather
         call filter_string(weather(nn))

	 provider(nn)(1:11) = 'NWS        '     ! data provider (all from NWS)
	 wmoid(nn) = wmoid_in(i)
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
	 store_4(nn,1) = alt(i)                 ! altimeter setting (mb)
	 store_4(nn,2) = badflag                ! station pressure (mb)
	 store_4(nn,3) = mslp(i)                ! MSL pressure (mb)
	 store_4(nn,4) = float(dpchar(i))       ! 3-h press change character
         store_4(nn,5) = dp(i)                  ! 3-h press change (mb)
c
	 store_5(nn,1) = vis(i)                 ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = badflag                ! soil/water temperature
	 store_5(nn,4) = badflag                ! soil moisture
c
	 store_6(nn,1) = pcp1(i)                ! 1-h precipitation
	 store_6(nn,2) = pcp3(i)                ! 3-h precipitation
	 store_6(nn,3) = pcp6(i)                ! 6-h precipitation
	 store_6(nn,4) = pcp24(i)               ! 24-h precipitation
	 store_6(nn,5) = snowcvr(i)             ! snow cover
c
	 store_7(nn,1) = float(k_layers)        ! number of cloud layers
	 store_7(nn,2) = max24t(i)              ! 24-h max temperature
	 store_7(nn,3) = min24t(i)              ! 24-h min temperature
c
c.....	Store cloud info if we have any. 
c
	 if(k_layers .gt. 0) then
	   do ii=1,k_layers
	     store_cldht(nn,ii) = ht(ii,i)
	     store_cldamt(nn,ii)(1:1) = ' '
	     store_cldamt(nn,ii)(2:4) = cvr(ii,i)(1:3)
	   enddo !ii
	 endif
c
c
  125	 continue
c
c
c.....  That's it...lets go home.
c
	 print *,' Found ',n_sao_b,' METARs in the LAPS box'
	 print *,' Found ',n_sao_g,' METARs in the LAPS grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' ERROR.  No data available from READ_METAR.'
	 return
c
	 end
c
c
	function c2f(temp_c)
c
c       Takes a single value in deg C and returns deg F.
c
	c2f = (temp_c * 9./5.) + 32.
c
	return
	end
