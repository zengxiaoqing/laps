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
	subroutine get_local_obs(maxobs,maxsta,i4time,data_file,
     &                      eastg,westg,anorthg,southg,nn,
     &                      icount,
     &                      stations,store,wx,obstype,
     &                      store_emv,store_amt,store_hgt,
     &                      num_var,jstatus)
c
c*****************************************************************************
c
c	Routine to gather local mesonet data for LAPS.   
c
c	Changes:
c	  P. Stamus, NOAA/FSL  02-04-98  Original version (from get_cdot_obs).
c                              04-20-98  Add flags for bad update times.

c
c       NOTE: This routine is not as yet generic.  There are several places
c             where the code is hardwired and expects the data to be in a
c             known format.  This will change as the LDAD mesonet CDL improves.
c
c*****************************************************************************
c
	include 'netcdf.inc'
	include 'surface_obs.inc'
        integer maxobs, maxsta
c
	real*8  reptime(maxobs), time_t(maxobs), time_p(maxobs)
	real*8  time_dd(maxobs), time_ff(maxobs), time_ffg(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs)
	real*4  dd(maxobs), ff(maxobs), ffg(maxobs)
	real*4  stnp(maxobs)
	real*4	store(maxsta,num_var), store_hgt(maxsta,5)
c
	integer*4  itime60, before, after
	integer    rtime
	integer    recNum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*5, timech*9, time*4
	character  stationId(maxobs)*6
	character  stations(maxsta)*3
	character  store_amt(maxsta,5)*4, store_emv(maxsta,5)*1
	character  wx(maxsta)*8
	character  obstype(maxsta)*8
	character  save_stn(maxobs)*3
	character  data_file*80
	character  dataProvider(maxobs)*11 
c
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
c.....  Get the data from the NetCDF file.  First open the file.
c
	nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)
	if(nf_status .ne. NF_NOERR) then
	   print *, NF_STRERROR(nf_status)
	   print *, data_file
           return
	endif
c
c.....  Get size of recNum
c
	nf_status = NF_INQ_DIMID(nf_fid,'recNum',nf_vid)
	if(nf_status .ne. NF_NOERR) then
	   print *, NF_STRERROR(nf_status)
	   print *, 'dim recNum'
           return
	endif
	nf_status = NF_INQ_DIMLEN(nf_fid,nf_vid,recNum)
	if(nf_status .ne. NF_NOERR) then
	   print *, NF_STRERROR(nf_status)
	   print *, 'dim recNum'
           return
	endif
c
c.....  Call the read routine.
c
	call read_local(nf_fid , recNum, dataProvider, td,
     &     elev, lats, lons, reptime, 
     &     stationId, time_p, stnp, 
     &     time_t, t, dd, time_dd, ffg, time_ffg,
     &     ff, time_ff, istatus) 
c
	if(istatus .ne. 1) then
	   print *,'  No data available from READ_LOCAL.'
	   go to 990   
	endif
c
	n_local_all = recNum
c
c
c.....  First check the data coming from the NetCDF file.  There can be
c.....  "FloatInf" (used as fill value) in some of the variables.  These
c.....  are not handled the same by different operating systems.  In our
c.....  case, IBM systems make "FloatInf" into "NaN" and store them that
c.....  way in the LSO, which messes up other LAPS routines.  This code
c.....  checks for "FloatInf" and sets the variable to 'badflag'.  If the
c.....  "FloatInf" is in the lat or lon, we toss the whole ob since we 
c.....  don't know where it is.
c
	do i=1,n_local_all
c
c.....  Toss the ob if lat/lon bad by setting lat to badflag (-99.9),
c.....  which causes the bounds check to think its outside the LAPS domain.
c
	   if(lats(i)+1. .eq. lats(i)) lats(i) = badflag
	   if(lons(i)+1. .eq. lons(i)) lats(i) = badflag
	   if(elev(i)+1. .eq. elev(i)) lats(i) = badflag
	   if(elev(i) .gt. 5000.) lats(i) = badflag !bag stn 
	   if(elev(i) .lt. -900.) lats(i) = badflag !bag stn 
c
	   if(reptime(i)+1.  .eq.  reptime(i)) reptime(i)  = badflag
	   if(time_t(i)+1.   .eq.   time_t(i)) time_t(i)   = badflag
	   if(time_p(i)+1.   .eq.   time_p(i)) time_p(i)   = badflag
	   if(time_dd(i)+1.  .eq.  time_dd(i)) time_dd(i)  = badflag
	   if(time_ff(i)+1.  .eq.  time_ff(i)) time_ff(i)  = badflag
	   if(time_ffg(i)+1. .eq. time_ffg(i)) time_ffg(i) = badflag
c
	   if(t(i)+1.    .eq.    t(i)) t(i)    = badflag
	   if(td(i)+1.   .eq.   td(i)) td(i)   = badflag
	   if(ff(i)+1.   .eq.   ff(i)) ff(i)   = badflag
	   if(ffg(i)+1.  .eq.  ffg(i)) ffg(i)  = badflag
	   if(dd(i)+1    .eq.   dd(i)) dd(i)   = badflag
c
	enddo !i
c
c
c.....  Move the station names.
c
	do i=1,n_local_all
	   stname(i) = '     '
	   if(dataProvider(i)(1:4) .eq. 'CDOT') then
	      stname(i)(1:1) = stationId(i)(1:1)
	      stname(i)(2:2) = stationId(i)(3:3)
	      if(ichar(stationId(i)(4:4)) .lt. 32) then
		 stname(i)(3:3) = ' '
	      else
		 stname(i)(3:3) = stationId(i)(4:4)
	      endif
	   else
	      stname(i)(1:4) = stationId(i)(1:4)
	   endif
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
	icount = 0
	jfirst = 1
	do 125 i=1,n_local_all
c
c.....  Check if station is in the box.
c
	  if(lons(i).gt.east .or. lons(i).lt.west) go to 125
	  if(lats(i).gt.anorth .or. lats(i).lt.south) go to 125
c
c.....  Check to see if its in the desired time window.
c
	  itime60 = nint(reptime(i)) + 315619200
c
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
c
c.....  Check if its in the LAPS grid.
c
	 if(lons(i).gt.eastg .or. lons(i).lt.westg) go to 151
	 if(lats(i).gt.anorthg .or. lats(i).lt.southg) go to 151
 151	 continue
c
c.....	Set the report type.
c
	  obstype(nn) = '        '
	  obstype(nn)(1:8) = dataProvider(i)(1:8)
	  do k=8,1,-1
	    if(ichar(obstype(nn)(k:k)) .eq. 0) then
	         obstype(nn)(k:k) = " "
	    else
	         go to 152
	    endif
	  enddo !k
152	  continue
c
c.....  Convert units, while doing a quick and dirty qc on the data.
c
	  if(time_t(i) .gt. 0.) then
	     if(abs(time_t(i) - reptime(i)) .gt. 1200.) then
		t(i) = badflag
	     else
		t(i) = ((t(i) - 273.15) * 9./5.) + 32.     ! K to F
		if(t(i).lt.-50. .or. t(i).gt.130.)  t(i) = badflag !temp
	     endif
	  else
	     t(i) = badflag
	  endif
c
	  td(i) = ((td(i) - 273.15) * 9./5.) + 32. ! K to F
	  if(td(i).lt.-50. .or. td(i).gt.90.) td(i) = badflag !dewpoint
c
	  if(time_p(i) .gt. 0.) then
	     if(abs(time_p(i) - reptime(i)) .gt. 1200.) then
		stnp(i) = badflag
	     else
		stnp(i) = stnp(i) * .01
	     endif
	  else
	     stnp(i) = badflag
	  endif
          if(stnp(i).gt.1000. .or. stnp(i).lt.600.) 
     &                                         stnp(i) = badflag
c
	  if(time_dd(i) .gt. 0.) then
	     if(abs(time_dd(i) - reptime(i)) .gt. 1200.) then
		dd(i) = badflag
	     else
		if(dd(i).lt.0. .or. dd(i).gt.360.) dd(i) = badflag
	     endif
	  else
	     dd(i) = badflag
	  endif
c
	  if(time_ff(i) .gt. 0.) then
	     if(abs(time_ff(i) - reptime(i)) .gt. 1200.) then
		ff(i) = badflag
	     else
		ff(i) = 1.94254 * ff(i)	!m/s to kt
		if(ff(i).lt.0. .or. ff(i).gt.120.) ff(i) = badflag
	     endif
	  else
	     ff(i) = badflag
	  endif
c
	  if(time_ffg(i) .gt. 0.) then
	     if(abs(time_ffg(i) - reptime(i)) .gt. 1200.) then
		ffg(i) = badflag
	     else
		ffg(i) = 1.94254 * ffg(i) !m/s to kt
		if(ffg(i).lt.0. .or. ffg(i).gt.140.) ffg(i) = badflag
	     endif
	  else
	     ffg(i) = badflag
	  endif
c
	  ddg = dd(i)
	  if(ffg(i) .eq. badflag) ddg = badflag
c
c..... Output the data to the storage arrays
c
	 stations(nn) = stname(i)(1:3)          ! station name
	 wx(nn) = '        '
	 store(nn,1) = lats(i)			! station latitude
	 store(nn,2) = lons(i)			! station longitude
	 store(nn,3) = elev(i)	                ! station elevation
	 store(nn,4) = rtime			! observation time
	 store(nn,5) = t(i)             	! t (deg f)
	 store(nn,6) = td(i)			! td (deg f)
	 store(nn,7) = dd(i)                   	! dd (deg)
	 store(nn,8) = ff(i)			! ff (kt)
	 store(nn,9) = ddg			! gust dd (deg)
	 store(nn,10) = ffg(i)			! gust ff (kt)
	 store(nn,11) = stnp(i)                 ! station pressure (mb)
	 store(nn,12) = badflag			! MSL pressure (mb)
	 store(nn,13) = badflag			! altimeter setting (mb)
	 store(nn,14) = 0.                      ! number of cloud layers
	 store(nn,15) = badflag	                ! ceiling height (m msl)
	 store(nn,16) = badflag			! lowest cloud height (m msl)
	 store(nn,17) = badflag			! cld cover (tenths)
	 store(nn,18) = badflag			! visibility (miles)
	 store(nn,19) = badflag			! solar radiation 
	 store(nn,20) = badflag			! 3-h pressure tendency (mb)
c
c
  125	 continue
c
c
c.....  That's it...lets go home.
c
	 print *,' Found ',icount,' Local stations in the LAPS box'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 return
c
	 end
