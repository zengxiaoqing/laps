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

c There will probably have to be some adjustments made to the
c read routines (there's one routine for the data and one for 
c the station metadata), and possibly to the
c 'get_tmeso_obs' routine too.  I've made some
c assumptions there, and have tried to document in the 
c code where they are.  
c
c
        subroutine get_local_cwb(maxobs,maxsta,i4time,
     &                 path_to_local_data,
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
c======================================================================
c
c     Routine to gather the CWB ASCII Mesonet data, and store it for
c     LAPS use.
c     
c     Original:  P. Stamus, NOAA/FSL  08 Sep 1999
c     Changes:
c
c     Notes:
c	1. Code assumes that the station metadata file is in the same
c 	   directory as the raw observations.  If they are different,
c 	   change the variable 'infile' in the call to 
c          'read_tmeso_stntbl'.
c
c       2. Code assumes that the reported precip is a 1-h amount.  
c
c	3. The 'priority' array is not used at this time.
c
c======================================================================
c
	include 'surface_obs.inc'
c
c.....  Input 
c
	real lat(ni,nj), lon(ni,nj)
        character path_to_local_data*(*)
c
c.....  Output arrays
c
        real  store_1(maxsta,4), 
     &        store_2(maxsta,3), store_2ea(maxsta,3),
     &        store_3(maxsta,4), store_3ea(maxsta,2),
     &        store_4(maxsta,5), store_4ea(maxsta,2),
     &        store_5(maxsta,4), store_5ea(maxsta,4),
     &        store_6(maxsta,5), store_6ea(maxsta,2),
     &        store_7(maxsta,3),
     &        store_cldht(maxsta,5)
c
        integer    wmoid(maxobs)
c
        character  stations(maxsta)*20
        character  provider(maxobs)*11
        character  weather(maxobs)*25
        character  reptype(maxobs)*6, atype(maxobs)*6
        character  store_cldamt(maxsta,5)*4

c
c.....  Local arrays
c
        real t_in(maxsta), td_in(maxsta), rh_in(maxsta)
        real dd_in(maxsta), ff_in(maxsta)
        real sfcp_in(maxsta), pcp_in(maxsta)
	real lat_in(maxsta), lon_in(maxsta), elev_in(maxsta)

        real rtime_in(maxsta)
c
        character stn_in(maxsta)*5
c
c.....  Stuff for the mesonet metadata.
c
	real lat_master(maxsta),lon_master(maxsta),elev_master(maxsta)
	real priority(maxsta)
c
	character stn_master(maxsta)*5
c
c.....  Start.
c
	ibadflag = int(badflag)
c
c.....	Set jstatus flag for the local data to bad until we find otherwise.
c
	jstatus = -1
c
c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c
	box_length = box_size * 111.137  !km/deg lat (close enough for lon)
	ibox_points = box_length / (grid_spacing / 1000.)  !in km
c
	n_local_g = 0
	n_local_b = 0
c
c.....  Get the mesonet metadata (station information).
c
        call read_tmeso_stntbl(path_to_local_data,maxsta,badflag,
     &                         stn_master,
     &                         lat_master,lon_master,elev_master,
     &                         priority,num_master,istatus)
	if(istatus .ne. 1) go to 990
c
c.....  Get the mesonet data.
c
        call read_tmeso(path_to_local_data,maxsta,badflag,i4time,
     &                  stn_in,rtime_in,
     &                  t_in,td_in,rh_in,pcp_in,sfcp_in,dd_in,
     &                  ff_in,num_in,istatus)
c
	if(istatus .ne. 1) go to 990
	n_local_all = num_in
c
c.....  Match data with metadata for each station, then store
c.....  the metadata in arrays.
c
        do i=1,n_local_all
c
c.....  For each station, search the master list.  
c
	   do j=1,num_master
c
	      if(stn_in(i)(1:5) .eq. stn_master(j)(1:5)) then
c
c.....  Found one...store the location info.
c
		 lat_in(i) = lat_master(j)
		 lon_in(i) = lon_master(j)
		 elev_in(i) = elev_master(j)
c
	      endif
	   enddo !j
	enddo !i
c
c.....  Now that we have a matched set of data and metadata, 
c.....  process the stations.  Check to see if the station is
c.....  within the LAPS grid; if so, store it.
c
	jfirst = 1
	box_low = 1. - float(ibox_points)  !buffer on west/south side
	box_idir = float(ni + ibox_points) !buffer on east
	box_jdir = float(nj + ibox_points) !buffer on north
c
	do i=1,n_local_all
	   if(lat_in(i) .lt. -90.) go to 125	
	   call latlon_to_rlapsgrid(lat_in(i),lon_in(i),lat,lon,
     &                              ni,nj,ri_loc,rj_loc,istatus)
	   if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir) go to 125
	   if(rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) go to 125
c
c.....  Elevation ok?
c
	   if(elev_in(i) .eq. badflag) go to 125
	   if(elev_in(i).gt.5200. .or. elev_in(i).lt.-400.) go to 125
c
c.....  If you want to check the valid time, or if there are more than
c.....  one report from this station, put that here.
c

c          	
	   nn = nn + 1
	   n_local_b = n_local_b + 1
c
c.....  Within LAPS grid?
c
	   if(ri_loc.lt.1 .or. ri_loc.gt.float(ni)) go to 151  !off grid
	   if(rj_loc.lt.1 .or. rj_loc.gt.float(nj)) go to 151  !off grid
	   n_local_g = n_local_g + 1                           !on grid...count it
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
!    	  call s_len(stn_in(i), len)
!         stations(nn)(1:len) = stn_in(i)(1:len)

 	  call s_len(stn_in(i), len)
          if(len .ne. 0)then
              stations(nn)(1:len) = stn_in(i)(1:len) ! station name
          else
              write(6,*)' Warning in get_local_cwb: blank station name.'
     1                 ,' Assigning name ',i
              write(stations(nn),101)i
 101	      format(i5,15x)
          endif
c
	  atype(nn)(1:6) = 'MESONT'
c
	  reptype(nn)(1:6) = 'UNK   '
c
	  weather(nn)(1:25) = 'UNK                     '
	  provider(nn)(1:11) = 'CWB        '
	  wmoid(nn) = ibadflag
c 
	  store_1(nn,1) = lat_in(i)              ! station latitude
	  store_1(nn,2) = lon_in(i)              ! station longitude
	  store_1(nn,3) = elev_in(i)             ! station elevation (m)
	  store_1(nn,4) = rtime_in(i)            ! observation time
c	
	  store_2(nn,1) = t_in(i)                ! temperature (deg F)
	  store_2(nn,1) = td_in(i)               ! dew point (deg F)
	  store_2(nn,1) = rh_in(i)               ! Relative Humidity
c
          store_3(nn,1) = dd_in(i)               ! wind dir (deg)
          store_3(nn,2) = ff_in(i)               ! wind speed (kt)
          store_3(nn,3) = badflag                ! wind gust dir (deg)
          store_3(nn,4) = badflag                ! wind gust speed (kt)
c
          store_4(nn,1) = badflag                ! altimeter setting (mb)
          store_4(nn,2) = sfcp_in(i)             ! station pressure (mb)
          store_4(nn,3) = badflag                ! MSL pressure (mb)
          store_4(nn,4) = badflag                ! 3-h press change character
          store_4(nn,5) = badflag                ! 3-h press change (mb)
c
          store_5(nn,1) = badflag                ! visibility (miles)
          store_5(nn,2) = badflag                ! solar radiation 
          store_5(nn,3) = badflag                ! soil/water temperature
          store_5(nn,4) = badflag                ! soil moisture
c
          store_6(nn,1) = pcp_in(i)              ! 1-h precipitation
          store_6(nn,2) = badflag                ! 3-h precipitation
          store_6(nn,3) = badflag                ! 6-h precipitation
          store_6(nn,4) = badflag                ! 24-h precipitation
          store_6(nn,5) = badflag                ! snow cover
c
          store_7(nn,1) = float(kkk)             ! number of cloud layers
          store_7(nn,2) = badflag                ! 24-h max temperature
          store_7(nn,3) = badflag                ! 24-h min temperature
c
c.....  That's it for this station.
c
 125      continue
        enddo !i
c
c.....  All done.
c
        print *,' Found ',n_local_b,' mesonet stations in the LAPS box'
        print *,' Found ',n_local_g,' mesonet stations in the LAPS grid'
        print *,' '
        jstatus = 1            ! everything's ok...
      	write(6,*)' GET_LOCAL_CWB: exit'
      	return
c
 990  	continue
c
      	write(6,*)' GET_LOCAL_CWB: exit without reading data'
     	return
c
   	end

c
c
        subroutine read_tmeso(infile,maxsta,badflag,i4time,stn,rtime,       
     &                        t,td,rh,pcp,sfcp,dd,ff,num,istatus)
c
c======================================================================
c
c     Routine to read the CWB ASCII Mesonet files.
c     
c     Original:  P. Stamus, NOAA/FSL  08 Sep 1999
c     Changes:
c
c======================================================================
c
        real t(maxsta), td(maxsta), rh(maxsta)
        real dd(maxsta), ff(maxsta)
        real sfcp(maxsta), pcp(maxsta)
        real rtime(maxsta)
c
        character infile*(*), stn_id*5, stn(maxsta)*5, 
     1            a9_to_a8*8, a9time*9, a8time*8, a6time*6, filename*11       
c
c
c.....  Start here.  Fill the output arrays with something, then open
c.....	the file to read.
c
	istatus = 0
	do i=1,maxsta
	   stn(i)(1:5) = '     '
	   rtime(i) = badflag
	   t(i) = badflag
	   td(i) = badflag
	   rh(i) = badflag
 	   pcp(i) = badflag
 	   sfcp(i) = badflag
	   dd(i) = badflag
	   ff(i) = badflag
	enddo !i
c
        call make_fnam_lp(i4time,a9time,istatus)
        if(istatus .ne. 1)go to 990

        a8time = a9_to_a8(a9time)
        a6time = a8time(3:8)
!       a6time = '050826'

        filename = 'MSO.m'//a6time

        call s_len(infile,len_infile)
        open(11,file=infile(1:len_infile)//filename,status='old'       
     1                                                     ,err=980)
c
        num = 0

!       Skip header comments at the top of the file
        do iread = 1,16
            read(11,*,end=550,err=990)
        enddo
c
c.....  This starts the read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
c     
 500    continue
c
        read(11,900,end=550,err=990) idum,stn_id,iyr,imth,idy,ihr,imin, 
     &                               ispd,idir,irh,isfcp,it,itd,iprecip
 900    format(i2,1x,a5,5(1x,i2),1x,i3,1x,i2,1x,i3,1x,i4,2(1x,i3),1x,i5)
c
c.....  Check for valid date/time...if bad, toss this ob.
c
        if(iyr.lt.0 .or. imth.le.0 .or. idy.le.0 .or. 
     &     ihr.lt.0 .or. imin.lt.0) then 
           print *, ' Bad date/time at station: ', stn_id, idum
           go to 500
        endif
c
c.....  Have good date/time...store ob.  Adjust/scale variables while storing.
c
        num = num + 1   !add to count
c
	stn(num)(1:5) = stn_id(1:5)
c
        itime = (ihr * 100) + imin   !figure out the time
        rtime(num) = float(itime)
c
        if(idir.gt.36 .or. idir.lt.0) then
           dd(num) = badflag
        else
           dd(num) = float(idir * 10)
        endif
c
        if(ispd .lt. 0) then
           ff(num) = badflag
        else
           ff(num) = (float(ispd) * 0.1) * 1.94254 !conv m/s to kt
        endif
c
        if(it .le. -90) then
           t(num) = badflag
        else
           t(num) = (float(it) * 0.1) * 9/5 + 32 !conv C to F
        endif
c
        if(itd .le. -90) then
           td(num) = badflag
        else
           td(num) = (float(itd) * 0.1) * 9/5 + 32 !conv C to F
        endif
c
        if(irh .lt. 0) then
           rh(num) = badflag
        else
           rh(num) = float(irh)
        endif
c
        if(iprecip .lt. 0) then
           pcp(num) = badflag
        else
           pcp(num) = float(iprecip) * 0.1 * 0.03937 !conv mm to inch
        endif
c
        if(isfcp .le. 0) then
           sfcp(num) = badflag
        else
           ps = float(isfcp) * 0.1
           if(ps .lt. 500.) then
              ps = ps + 1000.
           endif
           sfcp(num) = ps
        endif
c
c.....  Go back for the next ob.
c
        go to 500
c
c.....  Hit end of file...that's it.
c
 550    continue
c
        print *,' Found ', num, ' mesonet stations.'
	istatus = 1
        return
c     
 980    continue
c
        write(6,*)' WARNING: could not open mesonet data file ',filename
	istatus = -1
        return
c     
 990    continue
c
        print *,' ** ERROR reading mesonet data.'
	istatus = -1
        return
c     
        end
c
c
        subroutine read_tmeso_stntbl(infile,maxsta,badflag,stn,
     &                               lat,lon,elev,prior,num,istatus)
c
c======================================================================
c
c     Routine to read station information for the CWB ASCII Mesonet 
c	data.
c     
c     Original:  P. Stamus, NOAA/FSL  08 Sep 1999
c     Changes:
c
c======================================================================
c
        real lat(maxsta), lon(maxsta), elev(maxsta)
        real prior(maxsta)
c
	integer stn_id
c
        character infile*(*), stn_name*5, stn(maxsta)*5
c
c
c.....  Start here.  Fill the output with something, then open the 
c.....	file to read.
c
	do i=1,maxsta
	   lat(i) = badflag
	   lon(i) = badflag
	   elev(i) = badflag
	   prior(i) = badflag
	   stn(i)(1:5) = '     '
	enddo !i
c
        call s_len(infile,len_infile)
        open(11,file=infile(1:len_infile)//'stn.table',status='old'
     1                                                ,err=990)
c
        num = 0

!       Skip header comments at the top of the file
        do iread = 1,6
            read(11,*,end=550,err=990)
        enddo
c
c.....  This starts the station read loop.  Since we don't know how many 
c.....  stations we have, read until we hit the end of file.
c     
 500    continue
c
        read(11,900,end=550,err=990) stn_id,stn_name,alat,alon,ielev,
     &                               aprior 
 900    format(i5,1x,a5,2x,f7.4,2x,f8.4,2x,i4,1x,f9.3)
c
c.....  Move station info to arrays for sending to calling routine.
c
	num = num + 1
	stn(num)(1:5) = stn_name(1:5)
	lat(num) = alat
	lon(num) = alon
	elev(num) = float(ielev)
	prior(num) = aprior
c
c.....  Go back for the next ob.
c
        go to 500
c
c.....  Hit end of file...that's it.
c
 550    continue
c
        print *,' Found ', num
     1         , ' mesonet stations in the station table.' 
        istatus = 1
        return
c     
 990    continue
c
        print *,' ** ERROR reading mesonet station table'
        istatus = 0
        return
c     
	end


