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
     &                      eastg,westg,anorthg,southg,nn,
     &                      n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,
     &                      stations,store,wx,obstype,
     &                      store_emv,store_amt,store_hgt,
     &                      num_var,jstatus)
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
c
c*****************************************************************************
c
	include 'surface_obs.inc'
c
        integer maxobs,maxsta
	real*8  timeobs(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs), dd(maxobs), ff(maxobs)
	real*4  tt(maxobs), ttd(maxobs)
	real*4  mslp(maxobs), alt(maxobs), ht(6,maxobs), ffg(maxobs)
	real*4  vis(maxobs), precip1(maxobs)
	real*4	store(maxsta,num_var), store_hgt(maxsta,5)
c
	integer*4  correction(maxobs), itime60
	integer*4  before, after
	integer*2  rtime, dpchar(maxobs), dp(maxobs)
c
	character  stname(maxobs)*5, atype(maxobs)*6, timech*9, time*4
	character  weather(maxobs)*25 
	character  store_amt(maxsta,5)*4, store_emv(maxsta,5)*1
	character  stations(maxsta)*3, wx(maxsta)*8,auto*5,rept*4
	character  obstype(maxsta)*8,reptype(maxobs)*6
	character  data_file*80
	character  save_stn(maxobs)*3
c
	character  stname_in(maxobs)*5, cvr(6,maxobs)*8
c
c.....	Set jstatus flag for the sao data to bad until we find otherwise.
c
	jstatus = -1
c
c.....	Areal outline for the 'box' around the LAPS grid.
c
	east = eastg + box_size
	west = westg - box_size
	anorth = anorthg + box_size
	south  = southg - box_size
c
c.....	Zero out the counters.
c
	n_sao_g = 0		! # of saos in the laps grid
	n_sao_pos_g = 0		! total # of saos possible in laps grid
	n_sao_b = 0		! # of saos in the box
	n_sao_pos_b = 0		! total # of saos possible in the box
c
c.....  Call the routine that reads the METAR data files.
c
        call read_metar(data_file,n_sao_all,maxobs,
     &                   stname_in,
     &                   lats,lons,elev,
     &                   timeobs,reptype,atype,
     &                   cvr,ht,vis,weather,
     &                   mslp,t,tt,td,ttd,dd,ff,ffg,alt,
     &                   dpchar,dp,precip1,istatus)
c
	if(istatus .ne. 1) go to 990
c
c.....  First check the data coming from the NetCDF file.  There can be
c.....  "FloatInf" (used as fill value) in some of the variables.  These
c.....  are not handled the same by different operating systems.  In our
c.....  case, IBM systems make "FloatInf" into "NaN" and store them that
c.....  way in the LSO, which messes up other LAPS routines.  This code
c.....  checks for "FloatInf" and sets the variable to 'badflag'.  If the
c.....  "FloatInf" is in the lat, lon, or elevation, we toss the whole ob
c.....  since we don't know where it is.
c
	do i=1,n_sao_all
c
c.....  Toss the ob if lat/lon/elev bad by setting lat to badflag (-99.9),
c.....  which causes the bounds check to think its outside the LAPS domain.
c
	   if(lats(i)+1. .eq. lats(i)) lats(i) = badflag
	   if(lons(i)+1. .eq. lons(i)) lats(i) = badflag
	   if(elev(i)+1. .eq. elev(i)) lats(i) = badflag
c
	   do j=1,5
	      if(ht(j,i)+1. .eq. ht(j,i)) ht(j,i) = badflag
	   enddo !j
c
	   if(vis(i)+1.  .eq.  vis(i)) vis(i)  = badflag
	   if(mslp(i)+1. .eq. mslp(i)) mslp(i) = badflag
	   if(t(i)+1.    .eq.    t(i)) t(i)    = badflag
	   if(td(i)+1.   .eq.   td(i)) td(i)   = badflag
	   if(tt(i)+1.   .eq.   tt(i)) tt(i)   = badflag
	   if(ttd(i)+1.  .eq.  ttd(i)) ttd(i)  = badflag
	   if(dd(i)+1.   .eq.   dd(i)) dd(i)   = badflag
	   if(ff(i)+1.   .eq.   ff(i)) ff(i)   = badflag
	   if(ffg(i)+1.  .eq.  ffg(i)) ffg(i)  = badflag
	   if(alt(i)+1.  .eq.  alt(i)) alt(i)  = badflag
c
	enddo !i
c
c.....  Fix for metars....until can rig better
c
	do i=1,n_sao_all
	   stname(i) = '     '
cc	   write(6,9999) i, stname_in(i)
	   stname(i)(1:3) = stname_in(i)(2:4)
	enddo !i
 9999	format(1x,i4,2x,a5)
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
	do 125 i=1,n_sao_all
c
c.....  Check if station is in the box.
c
	  if(lons(i).gt.east .or. lons(i).lt.west) go to 125
	  if(lats(i).gt.anorth .or. lats(i).lt.south) go to 125
c
c.....  Check to see if its in the desired time window (if the flag
c.....  says to check this).
c
	  itime60 = nint(timeobs(i)) + 315619200
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
	     save_stn(1) = stname(i)(1:3)
	     jfirst = 0
	     go to 150
	  endif
c
	  do k=1,icount
	     if(stname(i)(1:3) .eq. save_stn(k)) go to 125
	  enddo !k
c
	  icount = icount + 1
	  save_stn(icount) = stname(i)(1:3)   ! only one...save for checking
c
 150	  nn = nn + 1
	  n_sao_b = n_sao_b + 1		!station is in the box
c
c.....  Check if its in the LAPS grid.
c
	 if(lons(i).gt.eastg .or. lons(i).lt.westg) go to 151
	 if(lats(i).gt.anorthg .or. lats(i).lt.southg) go to 151
	 n_sao_g = n_sao_g + 1
 151	 continue
c
c.....	Figure out the report type.
c
cc	  rept = '    '
cc	  auto = '     '
cc	  rept = reptype(i)(1:4)
cc	  auto = atype(i)(1:5)

	  obstype(nn) = '        '
	  obstype(nn)(1:5) = reptype(i)(1:5)
	  if(atype(i)(1:1) .eq. 'A') obstype(nn)(8:8) = 'A'
	  if(atype(i)(3:3) .eq. '1') obstype(nn)(7:7) = '1'
	  if(atype(i)(3:3) .eq. '2') obstype(nn)(7:7) = '2'
c
cc	  if(rept(3:3) .eq. ' ') then
cc	    obstype(nn)(1:4) = rept
cc	    obstype(nn)(4:8) = auto
cc	  else
cc	    obstype(nn)(4:8) = auto
cc	    obstype(nn)(1:4) = rept
cc	  endif
c
c+++++++++++++++++++++++++++++++++++
c
c       Temp fix for missing elev
c
c
	  if(stname(i)(1:3) .eq. 'ITR') elev(i) = 1125.

c+++++++++++++++++++++++++++++++++++
c
c
c.....	Figure out the cloud data.
c
	  kkk = 0               ! number of cloud layers
	  cover = badflag
	  hgt_ceil = 22500.0	! set to no ceiling to start
	  hgt_low = 22500.0	! set to no lowest cloud to start
c
	  if(cvr(1,i)(1:1) .eq. ' ') then
	     kkk = 0
	  else
	     do k=1,5
		if(cvr(k,i)(1:1) .ne. ' ') kkk = kkk + 1
	     enddo !k
	  endif
c
c
cc	kkk = 0
c
c
	  if(kkk .eq. 0) then	! no cloud data...probably AMOS station
	    cover    = badflag
	    hgt_ceil = badflag
	    hgt_low  = badflag
	    go to 126		! skip rest of cloud stuff
	  endif
	  iceiling = 0
	  ilow = 0
	  do ii=1,kkk
	    if(ht(ii,i) .gt. 25000.0) ht(ii,i) = badflag
	    if(cvr(ii,i)(1:3).eq.'CLR' .or.
     &                            cvr(ii,i)(1:3).eq.'SKC') then
	      cover = 0.
	      cvr(ii,i)(1:3) = 'CLR'
	      if(atype(i)(1:1) .eq. ' ') then    ! man station
		 hgt_ceil = 22500.0
		 hgt_low = 22500.0
		 ht(ii,i) = 22500.0
	      else                        ! auto station (12000 ft hgt limit)
		 hgt_ceil = 3657.4
		 hgt_low = 3657.4
		 ht(ii,i) = 3657.4
	      endif
	    elseif(cvr(ii,i)(1:3).eq.'FEW' 
     &                      .or. cvr(ii,i)(1:4).eq.'-FEW') then
	      cover = .1
	      if(ilow .eq. 0) then
	        hgt_low = ht(ii,i)
	        ilow = 1
	      endif
	    elseif(cvr(ii,i)(1:3).eq.'SCT' 
     &                      .or. cvr(ii,i)(1:4).eq.'-SCT') then
	      cover = .3
	      if(ilow .eq. 0) then
	        hgt_low = ht(ii,i)
	        ilow = 1
	      endif
	    elseif(cvr(ii,i)(1:3).eq.'BKN' 
     &                      .or. cvr(ii,i)(1:4).eq.'-BKN') then
	      cover = .7
	      if(iceiling .eq. 0) then
	        hgt_ceil = ht(ii,i)
	        iceiling = 1
	      endif
	      if(ilow .eq. 0) then
	        hgt_low = ht(ii,i)
	        ilow = 1
	      endif
	    elseif(cvr(ii,i)(1:3).eq.'OVC' 
     &                      .or. cvr(ii,i)(1:4).eq.'-OVC' 
     &                      .or. cvr(ii,i)(1:2).eq.'VV') then
	      cover = 1.
	      if(iceiling .eq. 0) hgt_ceil = ht(ii,i)
	      if(ilow .eq. 0) hgt_low = ht(ii,i)
	      if(cvr(ii,i)(1:2) .eq. 'VV') cvr(ii,i)(1:2) = 'X '
	    endif
	  enddo !ii
 126	  continue
c
c.....	check cloud info for very high heights...set to max if greater
c.....	also convert agl cloud heights to msl by adding elevation.
c
	if(hgt_ceil .ge. 0.) then	! if not...badflag
	  hgt_ceil = hgt_ceil + elev(i)	! convert cld hts from agl to msl
	  if(hgt_ceil .gt. 22500.) hgt_ceil = 22500.
	endif
	if(hgt_low .ge. 0.) then	! if not...badflag
	  hgt_low = hgt_low + elev(i)	! convert agl to msl
	  if(hgt_low  .gt. 22500.) hgt_low  = 22500.
	endif
	if(kkk .gt. 0) then
	  do ii=1,kkk
	    if(ht(ii,i) .ge. 0.) then
	      ht(ii,i) = ht(ii,i) + elev(i)		! conv agl to msl
	      if(ht(ii,i) .gt. 22500.) ht(ii,i) = 22500.
	    endif
	  enddo !ii
	endif
c
c
c.....  Convert units, while doing a quick and dirty qc on the data.
c

	if(alt(i).gt.25. .and. alt(i).lt.34.) 
     &                            alt(i) = alt(i) * 33.8624

	alt(i) = alt(i) * 0.01                                   !Pa to mb
	mslp(i) = mslp(i) * 0.01                                 !Pa to mb
	if(alt(i).lt.900. .or. alt(i).gt.1150.) alt(i) = badflag !alt (mb)
	if(mslp(i).lt.900. .or. mslp(i).gt.1150.) mslp(i) = badflag !MSL p (mb)
c
	t(i) = ((t(i) - 273.16) * 9./5.) + 32.
	td(i) = ((td(i) - 273.16) * 9./5.) + 32.
	if(t(i).lt.-50. .or. t(i).gt.130.)  t(i) = badflag   !temp
	if(td(i).lt.-50. .or. td(i).gt.90.) td(i) = badflag  !dewpoint
	tt(i) = ((tt(i) - 273.16) * 9./5.) + 32.
	ttd(i) = ((ttd(i) - 273.16) * 9./5.) + 32.
	if(tt(i).lt.-50. .or. tt(i).gt.130.)  tt(i) = badflag   !temp
	if(ttd(i).lt.-50. .or. ttd(i).gt.90.) ttd(i) = badflag  !dewpoint
c
	temp_f = tt(i)
	dewp_f = ttd(i)
	if(temp_f .eq. badflag) temp_f = t(i)
	if(dewp_f .eq. badflag) dewp_f = td(i)
c
	if(dd(i).lt.0. .or. dd(i).gt.360.)dd(i) = badflag
	ff(i) = 1.94254 * ff(i)                             !m/s to kt
	ffg(i) = 1.94254 * ffg(i)                           !m/s to kt
	if(ff(i).lt.0. .or. ff(i).gt.120.)ff(i) = badflag
	if(ffg(i).lt.0. .or. ffg(i).gt.140.)ffg(i) = badflag
	ddg = dd(i)
	if(ffg(i) .eq. badflag) ddg = badflag
c
	vis(i) = vis(i) * .001				!m to km
	vis(i) = 0.621371 * vis(i)                      !km to miles
	if(vis(i) .gt. 200.) vis(i) = badflag			! visibility
	if(vis(i) .lt.   0.) vis(i) = badflag
c
c.....  Figure out the pressure tendency.  NOTE:  As of 11-94, there
c.....  is a bug in the NetCDF SAO files, such that pressure changes
c.....  over 9.9 mb are not recorded properly.  When this is fixed in
c.....  the NetCDF files, this section of code will have to be changed.
c
	dp(i) = dp(i) * 0.01
	if(dpchar(i).le.-1 .or. dpchar(i).gt.9 .or. 
     &                                 dp(i).lt.-999) then
	   rdp3 = badflag
	else
	   rdp3 = float((dpchar(i) * 100) + dp(i))
	endif
c
c..... Output the data to the storage arrays
c
	 stations(nn) = stname(i)(1:3)          ! station name
	 wx(nn) = '        '
	 wx(nn) = weather(i)(1:8)               ! observed wx
	 store(nn,1) = lats(i)			! station latitude
	 store(nn,2) = lons(i)			! station longitude
	 store(nn,3) = elev(i)			! station elevation
	 store(nn,4) = rtime			! observation time
	 store(nn,5) = temp_f			! t (deg f)
	 store(nn,6) = dewp_f			! td (deg f)
	 store(nn,7) = dd(i)			! dd (deg)
	 store(nn,8) = ff(i)			! ff (kt)
	 store(nn,9) = ddg			! gust dd (deg)
	 store(nn,10) = ffg(i)			! gust ff (kt)
	 store(nn,11) = badflag			! station pressure (mb)
	 store(nn,12) = mslp(i)			! MSL pressure (mb)
	 store(nn,13) = alt(i)			! altimeter setting (mb)
	 store(nn,14) = float(kkk)		! number of cloud layers
	 store(nn,15) = hgt_ceil		! ceiling height (m msl)
	 store(nn,16) = hgt_low			! lowest cloud height (m msl)
	 store(nn,17) = cover			! cld cover (tenths)
	 store(nn,18) = vis(i)			! visibility (miles)
	 store(nn,19) = badflag			! solar radiation 
	 store(nn,20) = rdp3			! 3-h pressure tendency (mb)
c
c.....	Store cloud info if we have any.  New...make sure that the first
c.....  character in the string is either a '-' or a blank (necessary for
c.....  the LAPS cloud analysis).
c
	 if(kkk .gt. 0) then
	   do ii=1,kkk
	     if(cvr(ii,i)(1:2) .eq. '-X') ht(ii,i) = badflag !thin obscd=no ht
	     store_hgt(nn,ii) = ht(ii,i)
	     store_emv(nn,ii) = '   '             !emv(ii,i)
	     if(cvr(ii,i)(1:1) .eq. '-') then
	        store_amt(nn,ii)(1:4) = cvr(ii,i)(1:4)
	     else
		store_amt(nn,ii)(1:1) = ' '
		store_amt(nn,ii)(2:4) = cvr(ii,i)(1:3)
	     endif
	   enddo !ii
	 endif
c
c
  125	 continue
c
c
c.....  That's it...lets go home.
c
	 n_sao_pos_g = n_sao_g
	 n_sao_pos_b = n_sao_b
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
