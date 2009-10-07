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
c
	subroutine get_gps_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_gps_data,gps_format,
     &                      itime_before,itime_after,
     &                      eastg,westg,anorthg,southg,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_obs_g,n_obs_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)

c
c*****************************************************************************
c
cdoc	Routine to gather data from the gps files for LAPS.   
c
c*****************************************************************************
c
	include 'netcdf.inc'
c
c.....  Input variables/arrays.
c
        integer maxobs ! raw data file
        integer maxsta ! processed stations for LSO file
        character  path_to_gps_data*(*)
c
c.....  Local variables/arrays
c
        real    lat(ni,nj), lon(ni,nj)
	real*8  timeobs(maxobs)
	real  lats(maxobs), lons(maxobs), elev(maxobs)
	real  t(maxobs), rh(maxobs), stnp(maxobs)

	integer  i4time_ob, before, after, wmoid(maxobs)
	integer    rtime, dpchar(maxobs), iplat_type(maxobs)

        character*9 a9time_file
	character  stname(maxobs)*5, save_stn(maxobs)*8
	character  data_file*255, timech*9, time*4, gps_format*(*)
	character  weather(maxobs)*25, wx(maxsta)*25
	character  reptype(maxobs)*6, atype(maxobs)*6
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
c
	integer    recNum, nf_fid, nf_vid, nf_status
	character  stations(maxsta)*20, provider(maxsta)*11
	character  store_cldamt(maxsta,5)*4
c
c
c.....  Start.
c
	ibadflag = int( badflag )
c
c.....	Set jstatus flag for the gps data to bad until we find otherwise.
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
        ibox_points = box_length / (grid_spacing / 1000.)  !in km
c
c.....	Zero out the counters.
c
	n_obs_g = 0		! # of gps obs in the laps grid
	n_obs_b = 0		! # of gps obs in the box

        call s_len(gps_format, len_gps_format)
        if(gps_format(1:len_gps_format) .eq. 'NIMBUS' .or. ! FSL NetCDF format
     1     gps_format(1:len_gps_format) .eq. 'CWB'   )then ! CWB NetCDF format
            call s_len(path_to_gps_data,len_path)

            call get_file_time(path_to_gps_data(1:len_path)
     1                        ,i4time_sys,i4time_nearest)
            if(abs(i4time_sys - i4time_nearest) .le. 900)then
                i4time_file = i4time_nearest
            else
                write(6,*)'No data available within 900s of systime'
                goto990
            endif
c
!           i4time_file = i4time_sys - 900
            call make_fnam_lp(i4time_file,a9time_file,istatus)

            if(gps_format(1:len_gps_format) .eq. 'NIMBUS') then
               data_file = path_to_gps_data(1:len_path)//a9time_file
     1                                                 //'0030o.nc'
            elseif(gps_format(1:len_gps_format) .eq. 'CWB') then
               data_file = path_to_gps_data(1:len_path)//a9time_file
     1                                                 //'00.cwb.nc'
            endif
c
c.....      Get the data from the NetCDF file.  First, open the file.
c.....      If not there, return.
c
	    nf_status = NF_OPEN(data_file,NF_NOWRITE,nf_fid)

	    if(nf_status.ne.NF_NOERR) then
	       print *, NF_STRERROR(nf_status)
	       print *, data_file
               go to 990
            else
               write(6,*)' File opened successfully'
!              write(6,*)' Returning since GPS code is not tested yet'
!              go to 990
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

	    call read_gps(nf_fid , recNum, 
     &         stnp, elev, lats, lons, 
     &         t, rh, timeobs, stname)

            i4time_offset=315619200
c
        else
            write(6,*)' Unknown GPS format ',gps_format
            go to 990

        endif

	if(istatus .ne. 1) go to 990
	n_gps_raw = recNum
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
	do i=1,n_gps_raw
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
	   if( nanf( stnp(i) ) .eq. 1 ) stnp(i)  = badflag
	   if( nanf( t(i)    ) .eq. 1 ) t(i)     = badflag
	   if( nanf( rh(i)   ) .eq. 1 ) rh(i)    = badflag
c
	enddo !i
c
c.....  Set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after
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
	do 125 i=1,n_gps_raw
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
	  i4time_ob = nint(timeobs(i)) + i4time_offset
	  if(i4time_ob.lt.before .or. i4time_ob.gt.after) go to 125
c
c.....  Right time, right location...

 	  call make_fnam_lp(i4time_ob,timech,istatus)
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
	  n_obs_b = n_obs_b + 1     !station is in the box
c
c.....  Check if its in the LAPS grid.
c
          if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151  !off grid
          if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151  !off grid
	  n_obs_g = n_obs_g + 1  !on grid...count it
 151	  continue
c
c.....	Figure out the cloud data.
c.....  NOTE: Not reading cloud data from ship/gps file.  The data
c.....        is too ambiguous for LAPS use at this time.
c
	  kkk = 0               ! number of cloud layers
c
c
c.....  Convert units for storage.
c
c.....  Temperature and dewpoint
c
	temp_c = t(i)                         
	if(temp_c .lt. -85. .or. temp_c .gt. +70.)temp_c = badflag
	if(temp_c .eq. badflag) then                 ! t bad?
	   temp_f = badflag                          !          bag
	else
	   temp_f = c_to_f(temp_c)                   ! C to F
	endif
c
	dewp_f = badflag
c
c..... Pressure...Station Pressure
c
	if(stnp(i).lt.400. .or. stnp(i).gt.1200.) then
	   stnp(i) = badflag
	endif
c
c..... Fill the expected accuracy arrays.  Values are based on information
c..... in the 'Coastal-Marine Automated Network (C-MAN) Users Guide', and 
c..... we assume that they are about the same for gpss, ships, and C-MAN
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
         store_2ea(nn,2) = 0.0             ! Dew Point not reported
         store_2ea(nn,3) = 10.0            ! Relative Humidity %
c
c..... Wind (Not reported)
c
         store_3ea(nn,1) = 0.00            ! deg
         store_3ea(nn,2) = 0.00            ! kt
c
c..... Pressure and altimeter (mb)
c
         store_4ea(nn,1) = 1.00            ! pressure (mb)
         store_4ea(nn,2) = 0.00            ! altimeter (mb)
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
c..... Output the data to the storage arrays
c
	 call s_len(stname(i), len)
	 stations(nn)(1:len) = stname(i)(1:len) ! station name
	 provider(nn)(1:11) = 'NWS        '     ! data provider (all from NWS)
	 reptype(nn)(1:6) = 'GPS   '            ! station type
	 atype(nn)(1:6) ='UNK   '               ! used here for moving/fixed stns
	 wmoid(nn) = ibadflag                   ! WMO id...not applicable here
	 weather(nn)(1:25) = 
     1              'UNK                      ' ! present weather
c       
	 store_1(nn,1) = lats(i)                ! station latitude
	 store_1(nn,2) = lons(i)                ! station longitude
	 store_1(nn,3) = elev(i)                ! station elevation
	 store_1(nn,4) = rtime                  ! observation time
c
	 store_2(nn,1) = temp_f                 ! temperature (deg f)
	 store_2(nn,2) = badflag                ! dew point (deg f)
	 store_2(nn,3) = rh(i)                  ! Relative Humidity
c
	 store_3(nn,1) = badflag                ! wind dir (deg)
	 store_3(nn,2) = badflag                ! wind speed (kt)
	 store_3(nn,3) = badflag                ! wind gust dir (deg)
	 store_3(nn,4) = badflag                ! wind gust speed (kt)
c
	 store_4(nn,1) = badflag                ! altimeter setting (mb)
	 store_4(nn,2) = stnp(i)                ! station pressure (mb)
	 store_4(nn,3) = badflag                ! MSL pressure (mb)
         store_4(nn,4) = badflag
         store_4(nn,5) = badflag                ! 3-h press change (mb)
c
	 store_5(nn,1) = badflag                ! visibility (miles)
	 store_5(nn,2) = badflag                ! solar radiation 
	 store_5(nn,3) = badflag                ! soil/water temperature (F)
	 store_5(nn,4) = badflag                ! soil moisture 
c
	 store_6(nn,1) = badflag                ! 1-h precipitation (in)
	 store_6(nn,2) = badflag                ! 3-h precipitation
	 store_6(nn,3) = badflag                ! 6-h precipitation
	 store_6(nn,4) = badflag                ! 24-h precipitation
	 store_6(nn,5) = badflag                ! snow cover
c
	 store_7(nn,1) = 0                      ! number of cloud layers
	 store_7(nn,2) = badflag                ! 24-h max temperature
	 store_7(nn,3) = badflag                ! 24-h min temperature
c
  125	 continue
c
c
c.....  That's it...lets go home.
c
	 print *,' Found ',n_obs_b,' gps obs in the LAPS box'
	 print *,' Found ',n_obs_g,' gps obs in the LAPS grid'
	 print *,' '
	 jstatus = 1		! everything's ok...
	 return
c
 990	 continue		! no data available
	 jstatus = 0
	 print *,' WARNING.  No data available from READ_GPS'
	 return
c
	 end
