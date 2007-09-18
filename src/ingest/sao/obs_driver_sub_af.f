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
c******************************************************************************
c
c	Driver program for the LAPS surface data collection for the
c       Air Force version of LAPS.  This program gets the correct time 
c       and other stuff, then calls the routines that read the different 
c       data files.  Checks are made for stations that have reported but
c       are missing from the current hour.  Estimates are made for those 
c       stations, for the basic variables, using a Barnes analysis of the 
c       trends at nearby stations.  A normal LSO is written that includes 
c       the estimated "obs", and station and trend info is stored in a 
c       master file.
c
c	History:
c	   P. Stamus  02-20-97  Original (from AF version of obs_driver).
c                     10-28-97  Made dynamic - newlaps. Add interactive optn.
c                     12-09-98  Changes for new LSO format.
c                     03-29-99  Remove moving buoy/ship obs from master file.
c                     06-11-99  Pass lat/lon and grid size to get_sao_obs.
c                     06-17-99  Also pass grid_spacing to figure box size.
c
c******************************************************************************
c
c
c
c
	subroutine obs_driver_sub_af(ni,nj,maxsta,path_to_obs,
     &                               laps_cycle_time)
c
c
	parameter (badflag = -99.9)
	parameter (num_keep = 12)      !Num hrs to keep in master file.
c
	real fnorm(0:ni-1,0:nj-1)
c
	real  lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
        integer  wmoid(maxsta), jstatus, grid_spacing 
        integer    narg, iargc
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  weather(maxsta)*25 
        character  reptype(maxsta)*6, atype(maxsta)*6
        character  store_cldamt(maxsta,5)*4
	character  atime*24, outfile*256
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9
c
	real rii(maxsta), rjj(maxsta)
	real u_c(maxsta), v_c(maxsta)
c
	real trend_t(maxsta), trend_td(maxsta), trend_alt(maxsta)
	real trend_u(maxsta), trend_v(maxsta)
	real t_est(maxsta), td_est(maxsta), alt_est(maxsta)
	real u_est(maxsta), v_est(maxsta)
	real trend_t_est(maxsta), trend_td_est(maxsta)
	real trend_alt_est(maxsta)
	real trend_u_est(maxsta), trend_v_est(maxsta)
c
	integer ifound(maxsta)
	integer ii(maxsta), jj(maxsta)
c
	character wx(maxsta)*8, obstype(maxsta)*8
c
	character master_file*256, laps_domain*9, time*4, up_c*3
	character*(*) path_to_obs
c
	real mstn_lat(maxsta), mstn_lon(maxsta), mstn_elev(maxsta)
	real mstn_t(maxsta), mtrend_t(maxsta)
	real mstn_td(maxsta), mtrend_td(maxsta)
	real mstn_u(maxsta), mtrend_u(maxsta)
	real mstn_v(maxsta), mtrend_v(maxsta)
	real mstn_alt(maxsta), mtrend_alt(maxsta)
c
	integer mstn_ii(maxsta), mstn_jj(maxsta)
	integer n_updates(maxsta)
	integer rtime
c
        character mstn_name(maxsta)*5
c
	real mstn_lat_new(maxsta), mstn_lon_new(maxsta)
	real mstn_elev_new(maxsta)
	real mstn_t_new(maxsta), mtrend_t_new(maxsta)
	real mstn_td_new(maxsta), mtrend_td_new(maxsta)
	real mstn_u_new(maxsta), mtrend_u_new(maxsta)
	real mstn_v_new(maxsta), mtrend_v_new(maxsta)
	real mstn_alt_new(maxsta), mtrend_alt_new(maxsta)
c
	integer mstn_ii_new(maxsta), mstn_jj_new(maxsta)
	integer n_updates_new(maxsta)
c
        character mstn_name_new(maxsta)*5
c
c
c.....	Start here.  Check to see if this is an interactive run.
c
	narg = iargc()
c
c.....  Get the time from the scheduler, or from the user if interactive.
c
	if(narg .eq. 0) then
	   call get_systime(i4time, filename9, istatus)
	   call cv_i4tim_asc_lp(i4time, atime, status) ! get the atime
c
	else
c	   
 970	   write(6,973)
 973	   format(' Enter input filename (yydddhhmm): ',$)
	   read(5,972) filename9
 972	   format(a9)
	   call i4time_fname_lp(filename9(1:9),i4time,status)
	   i4time = i4time / laps_cycle_time * laps_cycle_time
	   call cv_i4tim_asc_lp(i4time, atime, status) ! get the atime
	endif
c
	write(6,900) filename9
900	format(' Getting surface data for: ',a9)
c
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c
	call get_directory('static',dir_s,len)
	ext_s = 'nest7grid'
c
	var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat,grid_spacing,istatus)
c
	var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon,grid_spacing,istatus)
c
c.....	Find east/west and north/south sides of grid (max extension of grid)
c
	grid_east = -999.
	grid_west = 999.
	grid_north = -90.
	grid_south = 90.
	do i=1,ni
	  if(lat(i,nj) .gt. grid_north) grid_north = lat(i,nj)
	  if(lat(i,1) .lt. grid_south) grid_south = lat(i,1)
	enddo !i
	do j=1,nj	
	  if(lon(ni,j) .gt. grid_east) grid_east = lon(ni,j)
	  if(lon(1,j) .lt. grid_west) grid_west = lon(1,j)
	enddo !j
c
c.....	Set up the counters, and zero/blank arrays.
c
	nn = 0
	n_obs_g = 0
	n_obs_b = 0
	n_sao_g = 0
	n_sao_b = 0
	n_local_g = 0
	n_local_b = 0
	n_buoy_g = 0
	n_buoy_b = 0
c
	do i=1,maxsta
	   stations(i) = '                    '
	   provider(i) = '           '
	   weather(i)  = '                         '
	   reptype(i)  = '      '
	   atype(i)    = '      '
c
	   do j=1,2
	      store_3ea(i,j) = 0.
	      store_4ea(i,j) = 0.
	      store_6ea(i,j) = 0.
	   enddo !j
c
	   do j=1,3
	      store_2(i,j) = badflag
	      store_7(i,j) = badflag
	      store_2ea(i,j) = 0.
	   enddo !j
c
	   do j=1,4
	      store_1(i,j) = badflag
	      store_3(i,j) = badflag
	      store_5(i,j) = badflag
	      store_5ea(i,j) = 0.
	   enddo !j
c
	   do j=1,5
	      store_4(i,j) = badflag
	      store_6(i,j) = badflag
	      store_cldht(i,j) = badflag
	      store_cldamt(i,j) = '    '
	   enddo !j
	enddo !i
c
c.....  Get the current surface data.
c
	iobs_flag = 1
c
cc	path_to_obs = './'
cc      path_to_obs = '/data/lapb/import/lapsdat/afwa/t5/sfc/'
	
	call get_sao_obs_af(filename9,path_to_obs,maxsta,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      nn,lat,lon,ni,nj,grid_spacing,
     &                      n_sao_g,n_sao_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider,jstatus)
c
	if(jstatus .ne. 1) then
	   iobs_flag = 0
	   go to 250
	endif
c
c.....  Finish the obs counts.
c
	n_obs_b = nn
	n_obs_g = n_sao_g 
c
c.....  Compute u and v components from dd and ff
c
	do i=1,n_obs_b
	   call decomp_wind(store_3(i,1),store_3(i,2),    !dd, ff
     &                       u_c(i),v_c(i),istatus)
	enddo !i
c
c.....  Find out where the current stations are on/around the LAPS grid.
c
        call find_ij(store_1(1,1),store_1(1,2),lat,lon,n_obs_b,
     &               maxsta,ni,nj,ii,jj,rii,rjj)
c
c.....  Fill some arrays.
c
 250	continue
	do i=1,maxsta
	   ifound(i) = -1
	   trend_t(i)   = badflag
	   trend_td(i)  = badflag
	   trend_u(i)   = badflag
	   trend_v(i)   = badflag
	   trend_alt(i) = badflag
	   n_updates(i) = 0
	enddo !i
c
c.....  Get the station and trend info from the master file.
c
	call get_directory('static',master_file,len)
	master_file = master_file(1:len) // 'SFC_master_file'
cc	master_file = './SFC_master_file'
	call read_master_file(master_file,maxsta,n_master,
     &               mstn_name,mstn_lat,mstn_lon,mstn_elev,mstn_ii,
     &               mstn_jj,n_updates,mstn_t,mtrend_t,mstn_td,
     &               mtrend_td,mstn_u,mtrend_u,mstn_v,mtrend_v,
     &               mstn_alt,mtrend_alt,istatus)
c
	if(istatus .ne. 1) then
	   print *,
     &     ' Problem with read_master_file. Using current hr only.'
	   go to 500
	endif
c
c.....  Calculate the trends for the current stations, and tag missing stns.
c
c
	do i=1,n_master
	  ifound(i) = 0
	  do j=1,n_obs_b
	    if(mstn_name(i)(1:5) .eq. stations(j)(1:5)) then
	      ifound(i) = 1
c
c.....  Do each variable
c
	      if(mstn_t(i)  .ne. badflag .and.
     &           store_2(j,1) .ne. badflag          ) then
		 trend_t(i) = store_2(j,1) - mstn_t(i) 
	      else
		 trend_t(i) = mtrend_t(i)
	      endif
c
	      if(mstn_td(i) .ne. badflag .and. 
     &           store_2(j,2) .ne. badflag          ) then
		 trend_td(i) = store_2(j,2) - mstn_td(i)
	      else
		 trend_td(i) = mtrend_td(i)
	      endif
c
	      if(mstn_u(i)  .ne. badflag .and. 
     &           u_c(j)     .ne. badflag             ) then
		 trend_u(i) = u_c(j) - mstn_u(i)
	      else
		 trend_u(i) = mtrend_u(i)
	      endif
	      if(mstn_v(i)  .ne. badflag .and. 
     &           v_c(j)     .ne. badflag             ) then
		 trend_v(i) = v_c(j) - mstn_v(i)
	      else
		 trend_v(i) = mtrend_v(i)
	      endif
c
	      if(mstn_alt(i) .ne. badflag .and. 
     &           store_4(j,1) .ne. badflag        ) then
		 trend_alt(i) = store_4(j,1) - mstn_alt(i)
	      else
		 trend_alt(i) = mtrend_alt(i)
	      endif
c
	    endif
	  enddo !j
	enddo !i
c
c.....  Now calculate trends for missing stations.  Analyze the trends at
c.....  stations we have to the missing ones using a Barnes-type analysis.
c
	rom2 = 0.005
	call dynamic_wts_af(ni,nj,n_obs_b,rom2,d,fnorm)
c
	do k=1,n_master
	   if(ifound(k) .eq. 0) then
	      if(iobs_flag .eq. 0) then
		 trend = mtrend_t(k)
	      else
		 call barnes_af(trend,ni,nj,mstn_ii(k),mstn_jj(k),
     &               ii,jj,trend_t,n_obs_b,badflag,fnorm)
	      endif
	      if(trend .eq. badflag .or. 
     &           mstn_t(k) .eq. badflag) then
		 t_est(k) = badflag
		 trend_t_est(k) = badflag
	      else
		 t_est(k) = mstn_t(k) + trend
		 trend_t_est(k) = trend
	      endif
c
	      if(iobs_flag .eq. 0) then
		 trend = mtrend_td(k)
	      else
		 call barnes_af(trend,ni,nj,mstn_ii(k),mstn_jj(k),
     &               ii,jj,trend_td,n_obs_b,badflag,fnorm)
	      endif
	      if(trend .eq. badflag .or.
     &           mstn_td(k) .eq. badflag) then
		 td_est(k) = badflag
		 trend_td_est(k) = badflag
	      else
		 td_est(k) = mstn_td(k) + trend
		 trend_td_est(k) = trend
	      endif
c
	      if(iobs_flag .eq. 0) then
		 trend = mtrend_u(k)
	      else
		 call barnes_af(trend,ni,nj,mstn_ii(k),mstn_jj(k),
     &               ii,jj,trend_u,n_obs_b,badflag,fnorm)
	      endif
	      if(trend .eq. badflag .or. 
     &           mstn_u(k) .eq. badflag) then
		 u_est(k) = badflag
		 trend_u_est(k) = badflag
	      else
		 u_est(k) = mstn_u(k) + trend
		 trend_u_est(k) = trend
	      endif
c
	      if(iobs_flag .eq. 0) then
		 trend = mtrend_v(k)
	      else
		 call barnes_af(trend,ni,nj,mstn_ii(k),mstn_jj(k),
     &               ii,jj,trend_v,n_obs_b,badflag,fnorm)
	      endif
	      if(trend .eq. badflag .or.
     &           mstn_v(k) .eq. badflag) then
		 v_est(k) = badflag
		 trend_v_est(k) = badflag
	      else
		 v_est(k) = mstn_v(k) + trend
		 trend_v_est(k) = trend
	      endif
c
	      if(iobs_flag .eq. 0) then
		 trend = mtrend_alt(k)
	      else
		 call barnes_af(trend,ni,nj,mstn_ii(k),mstn_jj(k),
     &               ii,jj,trend_alt,n_obs_b,badflag,fnorm)
	      endif
	      if(trend .eq. badflag .or.
     &           mstn_alt(k) .eq. badflag) then
		 alt_est(k) = badflag
		 trend_alt_est(k) = badflag
	      else
		 alt_est(k) = mstn_alt(k) + trend
		 trend_alt_est(k) = trend
	      endif
c
	   endif
	enddo !k

c
c.....  Update the master file.
c
c.....  This section updates with current obs found in the master file, 
c.....  and with current obs not previously in the master file.
c
 500    continue
	istn_count = 0
	do k=1,n_obs_b
	   if(reptype(k)(1:4) .eq. 'BUOY' .or.
     &        reptype(k)(1:4) .eq. 'SHIP') go to 410
	   istn_count = istn_count + 1  !count fixed stations
	   i = istn_count
	   mstn_name_new(i)(1:5) = stations(k)(1:5)
	   mstn_lat_new(i)       = store_1(k,1)
	   mstn_lon_new(i)       = store_1(k,2)
	   mstn_elev_new(i)      = store_1(k,3)
	   mstn_ii_new(i)        = ii(k)
	   mstn_jj_new(i)        = jj(k)
	   n_updates_new(i)      = 0
	   mstn_t_new(i)         = store_2(k,1)
	   mtrend_t_new(i)       = 0.
	   mstn_td_new(i)        = store_2(k,2)
	   mtrend_td_new(i)      = 0.
	   mstn_u_new(i)         = u_c(k)
	   mtrend_u_new(i)       = 0.
	   mstn_v_new(i)         = v_c(k)
	   mtrend_v_new(i)       = 0.
	   mstn_alt_new(i)       = store_4(k,1)
	   mtrend_alt_new(i)     = 0.
 410	   continue
	enddo !k
c
c.....  This section updates with stations previously in the master file,
c.....  but not in the current set of obs.
c
	icount_m = istn_count  !n_obs_b
	do i=1,n_master
	   if(ifound(i) .eq. 0) then
	      n_up = n_updates(i) + 1
	      if(n_up .le. num_keep) then
		 icount_m = icount_m + 1
		 mstn_name_new(icount_m)(1:5) = mstn_name(i)(1:5)
		 mstn_lat_new(icount_m)       = mstn_lat(i)
		 mstn_lon_new(icount_m)       = mstn_lon(i)
		 mstn_elev_new(icount_m)      = mstn_elev(i)
		 mstn_ii_new(icount_m)        = mstn_ii(i)
		 mstn_jj_new(icount_m)        = mstn_jj(i)
		 n_updates_new(icount_m)      = n_up
		 mstn_t_new(icount_m)         = t_est(i)
		 mtrend_t_new(icount_m)       = trend_t_est(i)
		 mstn_td_new(icount_m)        = td_est(i)
		 mtrend_td_new(icount_m)      = trend_td_est(i)
		 mstn_u_new(icount_m)         = u_est(i)
		 mtrend_u_new(icount_m)       = trend_u_est(i)
		 mstn_v_new(icount_m)         = v_est(i)
		 mtrend_v_new(icount_m)       = trend_v_est(i)
		 mstn_alt_new(icount_m)       = alt_est(i)
		 mtrend_alt_new(icount_m)     = trend_alt_est(i)
	      endif
	   endif
	enddo !i
	write(6,901) n_obs_b, istn_count, icount_m
 901	format(' Update master file (all/no moving/total): ',3i8)
c
c.....  Write out the updated master file.  It will only have the current
c.....  obs if there were problems with the read earlier.
c
	call write_master_file(master_file,maxsta,icount_m,
     &         mstn_name_new,mstn_lat_new,mstn_lon_new,
     &         mstn_elev_new,mstn_ii_new,mstn_jj_new,n_updates_new,
     &         mstn_t_new,mtrend_t_new,mstn_td_new,mtrend_td_new,
     &         mstn_u_new,mtrend_u_new,mstn_v_new,mtrend_v_new,
     &         mstn_alt_new,mtrend_alt_new)
c
c.....  Update the store arrays for the LSO.
c
cc	do i=1,n_obs_b
cc	   stations(i)(1:3) = stn_c(i)(2:4)
cc	enddo !i
c
c.....  Add estimated values for stations in the master file but not
c.....  in the current hour to the current hour arrays for LSO.
c
	time = filename9(6:9)
	read(time,*) rtime
c
	kk = n_obs_b 
	do i=1,n_master
	   if(ifound(i) .eq. 0) then
	      kk = kk + 1
	      stations(kk)(1:5) = mstn_name(i)(1:5)
	      write(up_c,299) n_updates(i) + 1
cc	      obstype(kk)  = 'EST     '
cc	      obstype(kk)(5:7) = up_c
	      reptype(kk)(1:3) = 'EST'
	      reptype(kk)(4:6) = up_c
	      atype(kk)(1:6) = '      '
	      wmoid(kk) = 0
	      provider(kk) = 'AFWA       '
	      weather(kk) = '                         '
c
	      store_1(kk,1)  = mstn_lat(i)
	      store_1(kk,2)  = mstn_lon(i)
	      store_1(kk,3)  = mstn_elev(i)
	      store_1(kk,4)  = rtime
c
	      store_2(kk,1)  = t_est(i)
	      store_2(kk,2)  = td_est(i)
	      store_2(kk,3)  = badflag
c
	      call windconv(u_est(i),v_est(i),dd,ff) 
	      store_3(kk,1)  = dd
	      store_3(kk,2)  = ff
	      store_3(kk,3)  = badflag
	      store_3(kk,4)  = badflag
c
	      store_4(kk,1)  = alt_est(i)
	      store_4(kk,2)  = badflag
	      store_4(kk,3)  = badflag
	      store_4(kk,4)  = badflag
	      store_4(kk,5)  = badflag
c
	      store_5(kk,1)  = badflag
	      store_5(kk,2)  = badflag
	      store_5(kk,3)  = badflag
	      store_5(kk,4)  = badflag
c
	      store_6(kk,1)  = badflag
	      store_6(kk,2)  = badflag
	      store_6(kk,3)  = badflag
	      store_6(kk,4)  = badflag
	      store_6(kk,5)  = badflag
c
	      store_7(kk,1)  = 0
	      store_7(kk,2)  = badflag
	      store_7(kk,3)  = badflag
c
	      do j=1,5
		 store_cldht(kk,j)  = badflag
		 store_cldamt(kk,j) = '    '
	      enddo !j
	    endif
	 enddo !i
	 if(kk .gt. n_obs_b) then
	    n_obs_b = kk
	    n_sao_b = kk
	 endif
 299	 format(i3)
c
c.....  Write the LSO file
c
	 call get_directory('lso',outfile,len)
	 outfile = outfile(1:len) // filename9(1:9) // '.lso'
cc        outfile = './' // filename9(1:9) // '.lso'
c
        call write_surface_obs(atime,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,weather,reptype,atype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c.....	That's about it...let's go home.
c
	write(6,*)'Normal completion of OBS_DRIVER_SUB_AF'

	end
c
c
	subroutine barnes_af(trend,imax,jmax,i_stn,j_stn,ii,jj,t_ob,
     &                    numsta,badflag,fnorm)
c
c.....	Routine to do a Barnes analysis that will consider stations in
c.....	the 't_ob' array that are outside the boundaries of the domain.
c
	real fnorm(0:imax-1,0:jmax-1)
	real t_ob(imax*jmax) 
	real val(imax*jmax)
	integer iob(imax*jmax), job(imax*jmax)
	integer ii(imax*jmax), jj(imax*jmax)
	integer dx, dy, i_stn, j_stn
c
c.....	loop over field only once
c
	ncnt = 0
	do n=1,numsta
          if (t_ob(n).ne.badflag) then
	    ncnt = ncnt + 1
	    iob(ncnt) = ii(n)
	    job(ncnt) = jj(n)
	    val(ncnt) = t_ob(n)
          endif
	enddo !n 
c
	if(ncnt .eq. 0) then
	  print *,' +++ NCNT = 0 in BARNES_WIDE. +++'
	  trend = badflag
	  return
	endif
cc	write(6,900) ncnt, numsta
cc900	format('   Selected ',i4,' obs out of ',i4,' total.')
c
	    sum2 = 0.
	    sumwt2 = 0.
	    imaxm1 = imax - 1
	    jmaxm1 = jmax - 1
	    do n=1,ncnt
	      dy = min(abs(j_stn - job(n)), jmaxm1) 
	      dx = min(abs(i_stn - iob(n)), imaxm1) 
	      sum2 = fnorm(dx,dy) * val(n) + sum2
	      sumwt2 = sumwt2 + fnorm(dx,dy)
	    enddo !n
c
c smart/mcginley 7-1-98 
c
	    if(sumwt2 .eq. 0.) then
                trend=badflag
                goto 500
            endif
c               print *,' barneswide wierd loop...........'
c               sum2 = 0.
c       sumwt2 = 0.
c       do n=1,ncnt
c          dx = min(abs(i_stn - iob(n)), imaxm1) 
c          dy = min(abs(j_stn - job(n)), jmaxm1) 
c          sum2 = fnorm(dx,dy) * val(n) + sum2
c          sumwt2 = sumwt2 + fnorm(dx,dy) 
c       enddo !n
c      if(sumwt2 .ne. 0.) then
c   go to 490 
c  else
c  trend = badflag
c  go to 500
c  endif
c   else
c       go to 490
c   endif 
c
490	    continue 
	    trend = sum2 / sumwt2
c
500 	  continue
c
!	print *,'   leaving barnes_wide'
	return
	end
c
c
	subroutine dynamic_wts_af(imax,jmax,n_obs_var,rom2,d,fnorm)
c
c=====================================================================
c
c     Routine to calculate the weights to be used in the Barnes
c     analysis.  The data density is used to set the cutoff for
c     the response function.  Then that cutoff is used to calculate
c     the exp, based on differences so that no additional distance
c     calculations are required in the Barnes routine.  All of this
c     is done in gridpoint space.
c
c     Original:  07-14-95  P. Stamus, NOAA/FSL
c     Changes:   
c                02-28-97  Air Force version
c
c     Notes:
c
c       1.  If variable 'rom2' is passed in as zero, it is calculated
c           from the data density.  Otherwise, the value passed in is
c           used in the weight calculation.
c
c       2.  The response for 2d waves is hard-wired in this routine.
c           This is the 'con' variable, and comes from the eqn:
c                     D = exp -(pi**2 R**2)/lamba**2
c           If we set D (the response) to our desired cutoff, set 
c           lamba to the desired wavelength in gridpt space (2d),
c           then solve for R in terms of d, we get the 'con' value
c           (i.e.,  R = (con)*d).  Here are some values for different
c           cutoffs:
c                     D = 0.01     R = 1.36616d
c                         0.10     R = 0.96602d
c                         0.25     R = 0.74956d
c                         0.50     R = 0.53002d
c
c=====================================================================
c
	real fnorm(0:imax-1,0:jmax-1)
        integer dx,dy
c
c.... First, find the area that each ob covers in gridpt space (this
c.... of course assumes a uniform coverage).
c
cc	con = 0.96602     ! resp of 0.10
cc	con = 0.74956     ! resp of 0.25
	con = 0.53002     ! resp of 0.50
	if(rom2 .eq. 0.) then
	   area = float(imax * jmax) / n_obs_var
	   d = sqrt( area )
	   rom2 = 1. / ((con * con) * (d * d))
	   write(6,900) n_obs_var, area, d, rom2
 900	   format(1x,'Num obs: ',i5,'  Area: ',f8.2,'  d: ',f8.2,
     &       '  rom2: ',f8.5)
	else
	   d = sqrt(1./(con * con * rom2))
	   write(6,902) rom2, d
 902	   format(' Using preset rom2 of: ',f8.5,'  Calc d: ',f8.2)
	endif
c
c.... Now calculate the weights for all the possible distances.
c
	pi = 4. * atan(1.)
	fno = 1. / (sqrt(2. * pi))
c
	do dy=0,jmax-1
	do dx=0,imax-1
	   rr = dx*dx + dy*dy
           fnorm(dx,dy) = fno * (exp( -(rr * rom2)))
	enddo !dx
	enddo !dy
c
c.... That's it.
c
	return
	end
