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
	subroutine get_cdot_obs(maxobs,maxsta,i4time,data_file,
     &                      eastg,westg,anorthg,southg,nn,
     &                      icount,
     &                      stations,store,wx,obstype,
     &                      store_emv,store_amt,store_hgt,
     &                      num_var,jstatus)
c
c*****************************************************************************
c
c	Routine to gather CDOT data (mesonets) for LAPS.   
c
c	Changes:
c		P. Stamus  10-30-96  Original version (from get_metar_obs).
c                          11-13-96  Pass in full path to data.
c                          01-09-97  Elev check for CDOT
c                          04-02-97  Change temps in from C to K...adj cks.
c
c*****************************************************************************
c
	include 'surface_obs.inc'
        integer maxobs, maxsta
c
	real*8  timeobs(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs), rh(maxobs)
	real*4  ff(maxobs), ffg(maxobs)
	real*4  precip1(maxobs)
	real*4	store(maxsta,num_var), store_hgt(maxsta,5)
	real*4  idd(maxobs), iddg(maxobs)
c
	integer*4  itime60, before, after
	integer*2  rtime
c
	character  stname(maxobs)*5, timech*9, time*4
	character  store_amt(maxsta,5)*4, store_emv(maxsta,5)*1
	character  stations(maxsta)*3, wx(maxsta)*8
	character  obstype(maxsta)*8
	character  save_stn(maxobs)*3, stname_in(maxobs)*10
	character  data_file*80
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
c.....  Call the routine that reads the local data files.
c
        call read_cdot(data_file,n_cdot,maxobs,
     &                   stname_in,
     &                   lats,lons,elev,
     &                   timeobs,
     &                   t,td,rh,idd,ff,iddg,ffg,
     &                   precip1,istatus)
	if(istatus .ne. 1) then
	   print *,'  No data available from READ_CDOT.'
	   go to 990   !change this when adding other lcl data
	endif
c
c.....  Other local data reads go here.
c
c
c.....  Now add up the total number of local data.
c
	n_local_all = n_cdot
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
ccc	   if(elev(i)+1. .eq. elev(i)) lats(i) = badflag
	   if(t(i)+1.    .eq.    t(i)) t(i)    = badflag
	   if(td(i)+1.   .eq.   td(i)) td(i)   = badflag
	   if(ff(i)+1.   .eq.   ff(i)) ff(i)   = badflag
	   if(ffg(i)+1.  .eq.  ffg(i)) ffg(i)  = badflag
c
	   if(idd(i)+1   .eq.  idd(i)) idd(i)  = badflag
	   if(iddg(i)+1  .eq. iddg(i)) iddg(i) = badflag
c
	enddo !i
c
c++++++++++++++++++++++++++++++++++++++++++
c
c       Temp. fix to stn elevs until they actually tell us what the
c       CDOT elevs are.....
c
c
	do i=1,n_local_all
	   if(stname_in(i)(1:4) .eq. '1-10') then
	      elev(i) = 1634.
	   elseif(stname_in(i)(1:4) .eq. '1-11') then
	      elev(i) = 1615.
	   elseif(stname_in(i)(1:4) .eq. '1-12') then
	      elev(i) = 1646.
	   elseif(stname_in(i)(1:4) .eq. '1-13') then
	      elev(i) = 1615.
	   elseif(stname_in(i)(1:4) .eq. '1-14') then
	      elev(i) = 1621.
	   elseif(stname_in(i)(1:4) .eq. '1-15') then
	      elev(i) = 1621.
	   elseif(stname_in(i)(1:4) .eq. '1-16') then
	      elev(i) = 1883.
	   elseif(stname_in(i)(1:4) .eq. '1-17') then
	      elev(i) = 1700.
	   elseif(stname_in(i)(1:3) .eq. '1-0') then
	      elev(i) = 1829.
	   elseif(stname_in(i)(1:3) .eq. '1-1') then
	      elev(i) = 1829.
	   elseif(stname_in(i)(1:3) .eq. '1-2') then
	      elev(i) = 1591.
	   elseif(stname_in(i)(1:3) .eq. '1-3') then
	      elev(i) = 1640.
	   elseif(stname_in(i)(1:3) .eq. '1-4') then
	      elev(i) = 1603.
	   elseif(stname_in(i)(1:3) .eq. '1-5') then
	      elev(i) = 1584.
	   elseif(stname_in(i)(1:3) .eq. '1-6') then
	      elev(i) = 1676.
	   elseif(stname_in(i)(1:3) .eq. '1-7') then
	      elev(i) = 1901.
	   elseif(stname_in(i)(1:3) .eq. '1-8') then
	      elev(i) = 1901.
	   elseif(stname_in(i)(1:3) .eq. '1-9') then
	      elev(i) = 1573.
c
	   elseif(stname_in(i)(1:3) .eq. '2-0') then
	      elev(i) = 2207.
	   elseif(stname_in(i)(1:3) .eq. '2-1') then
	      elev(i) = 1932.
	   elseif(stname_in(i)(1:3) .eq. '2-2') then
	      elev(i) = 3230.
	   elseif(stname_in(i)(1:3) .eq. '2-3') then
	      elev(i) = 2834.
	   elseif(stname_in(i)(1:3) .eq. '2-4') then
	      elev(i) = 1889.
	   elseif(stname_in(i)(1:3) .eq. '2-5') then
	      elev(i) = 2347.
	   elseif(stname_in(i)(1:3) .eq. '2-6') then
	      elev(i) = 2865.
c
	   elseif(stname_in(i)(1:4) .eq. '3-10') then
	      elev(i) = 3407.
	   elseif(stname_in(i)(1:4) .eq. '3-11') then
	      elev(i) = 2743.
	   elseif(stname_in(i)(1:4) .eq. '3-12') then
	      elev(i) = 1463.
	   elseif(stname_in(i)(1:4) .eq. '3-13') then
	      elev(i) = 1201.
	   elseif(stname_in(i)(1:4) .eq. '3-14') then
	      elev(i) = 1140.
	   elseif(stname_in(i)(1:3) .eq. '3-0') then
	      elev(i) = 1152.
	   elseif(stname_in(i)(1:3) .eq. '3-1') then
	      elev(i) = 1219.
	   elseif(stname_in(i)(1:3) .eq. '3-2') then
	      elev(i) = 2103.
	   elseif(stname_in(i)(1:3) .eq. '3-3') then
	      elev(i) = 1493.
	   elseif(stname_in(i)(1:3) .eq. '3-4') then
	      elev(i) = 1506.
	   elseif(stname_in(i)(1:3) .eq. '3-5') then
	      elev(i) = 1372.
	   elseif(stname_in(i)(1:3) .eq. '3-6') then
	      elev(i) = 1371.
	   elseif(stname_in(i)(1:3) .eq. '3-7') then
	      elev(i) = 1841.
	   elseif(stname_in(i)(1:3) .eq. '3-8') then
	      elev(i) = 1567.
	   elseif(stname_in(i)(1:3) .eq. '3-9') then
	      elev(i) = 1481.
c
	   elseif(stname_in(i)(1:3) .eq. '4-0') then
	      elev(i) = 2255.
	   elseif(stname_in(i)(1:3) .eq. '4-1') then
	      elev(i) = 1719.
	   elseif(stname_in(i)(1:3) .eq. '4-2') then
	      elev(i) = 1445.
	   elseif(stname_in(i)(1:3) .eq. '4-3') then
	      elev(i) = 1384.
	   elseif(stname_in(i)(1:3) .eq. '4-4') then
	      elev(i) = 1634.
	   elseif(stname_in(i)(1:3) .eq. '4-5') then
	      elev(i) = 3170.
c
	   elseif(stname_in(i)(1:3) .eq. '5-0') then
	      elev(i) = 2006.
	   elseif(stname_in(i)(1:3) .eq. '5-1') then
	      elev(i) = 1835.
	   elseif(stname_in(i)(1:3) .eq. '5-2') then
	      elev(i) = 1944.
	   elseif(stname_in(i)(1:3) .eq. '5-3') then
	      elev(i) = 2133.
	   elseif(stname_in(i)(1:3) .eq. '5-4') then
	      elev(i) = 2188.
	   elseif(stname_in(i)(1:3) .eq. '5-5') then
	      elev(i) = 2340.
	   elseif(stname_in(i)(1:3) .eq. '5-6') then
	      elev(i) = 1932.
	   elseif(stname_in(i)(1:3) .eq. '5-7') then
	      elev(i) = 1889.
	   endif
c
	   if(elev(i) .gt. 5000.) lats(i) = badflag !bag stn 
	   if(elev(i) .lt. -900.) lats(i) = badflag !bag stn 
c
	enddo !i
c
c++++++++++++++++++++++++++++++++++++++++++
c
c.....  Move the station names.
c
	do i=1,n_local_all
	   stname(i) = '     '
	   stname(i)(1:1) = stname_in(i)(1:1)
	   stname(i)(2:2) = stname_in(i)(3:3)
	   if(ichar(stname_in(i)(4:4)) .lt. 32) then
	      stname(i)(3:3) = ' '
	   else
	      stname(i)(3:3) = stname_in(i)(4:4)
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
	  itime60 = nint(timeobs(i)) + 315619200
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
	  obstype(nn) = 'CDOT    '
c
c.....	Figure out the cloud data.
c
	  kkk = 0               ! number of cloud layers
c
	  if(kkk .eq. 0) then	! no cloud data...probably AMOS station
	    cover    = badflag
	    hgt_ceil = badflag
	    hgt_low  = badflag
	  endif
c
c.....  Convert units, while doing a quick and dirty qc on the data.
c
	t(i) = ((t(i) - 273.15) * 9./5.) + 32.     ! K to F
	td(i) = ((td(i) - 273.15) * 9./5.) + 32.   ! K to F
	if(t(i).lt.-50. .or. t(i).gt.130.)  t(i) = badflag   !temp
	if(td(i).lt.-50. .or. td(i).gt.90.) td(i) = badflag  !dewpoint
c
	dd = idd(i)
	ddg = iddg(i)
	if(dd.lt.0. .or. dd.gt.360.)dd = badflag
	if(ddg.lt.0. .or. ddg.gt.360.)ddg = badflag
	ff(i) = 1.94254 * ff(i)                             !m/s to kt
	ffg(i) = 1.94254 * ffg(i)                           !m/s to kt
	if(ff(i).lt.0. .or. ff(i).gt.120.)ff(i) = badflag
	if(ffg(i).lt.0. .or. ffg(i).gt.140.)ffg(i) = badflag
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
	 store(nn,7) = dd			! dd (deg)
	 store(nn,8) = ff(i)			! ff (kt)
	 store(nn,9) = ddg			! gust dd (deg)
	 store(nn,10) = ffg(i)			! gust ff (kt)
	 store(nn,11) = badflag			! station pressure (mb)
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
	 print *,' Found ',icount,' CDOTs in the LAPS box'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 return
c
	 end
