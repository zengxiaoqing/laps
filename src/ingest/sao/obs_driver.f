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
	program obs_driver
c
c******************************************************************************
c
c	Driver program for the LAPS surface data collection.  This 
c	program gets the correct time and other stuff, then calls
c	the routines that read the different data files.
c
c	History:
c	   P. Stamus  11-09-94  Original version (Interactive from obs_driver)
c                     01-17-95  Turned on mesonet data.
c                     10-30-96  Changes for METARs and CDOTs.
c                     11-13-96  Porting improvments.
c                     12-11-96  Change to FD CDOT files.
c                     01-15-97  Check CDOT filename (pub vs wfo).
c                     03-27-97  Add ability to do interactive runs (removes
c                                 the need for 'obs_driveri').  Remove equivs.
c
c          J. Edwards 07-14-97  Made dynamic and moved data paths 
c                               to nest7grid.parms
c
c          P. Stamus  03-23-98  Changes for stand-alone QC; LS2 format.
c                     05-01-98  Added soil moisture variables.
c                     08-28-98  Added buoy and LDAD mesonet reads.
c                     09-04-98  Install as LSO, using new 'ls2' format.
c                     09-30-98  Housekeeping changes.
c	              12-03-98  Increase dir_s to 256 characters.
c                     06-21-99  Pass lat/lon, grid size and grid spacing
c                                 to get_ routines to calculate box size.
c
c       Notes:
c         1. When run "operationally", 'obs_driver.x' uses the time from
c            the 'systime.dat' file.  Running 'obs_driver.x -i' allows the
c            user to enter the run time.
c
c******************************************************************************
c
        include 'lapsparms.cmn'
        call get_laps_config('nest7grid',istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
	   stop
	endif
        call obs_driver_sub(NX_L_CMN,NY_L_CMN,maxstations_cmn,
     &             maxobs_cmn,path_to_metar_data_cmn,
     &             path_to_local_data_cmn,path_to_buoy_data_cmn,
     &             min_to_wait_for_metars_cmn,laps_cycle_time_cmn)

        end

        subroutine obs_driver_sub(ni,nj,maxsta,maxobs
     &          ,path_to_metar,path_to_local_data,path_to_buoy_data,
     &          minutes_to_wait_for_metars,laps_cycle_time)
c        
        include 'surface_obs.inc'
        integer ni, nj, maxsta, maxobs 
c
	real*4  lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real*4  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
        integer*4  wmoid(maxobs), jstatus, grid_spacing 
        integer    dpchar(maxobs), narg, iargc
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  weather(maxobs)*25 
        character  reptype(maxobs)*6, atype(maxobs)*6
        character  store_cldamt(maxsta,5)*4
	character  atime*24, outfile*200
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13
        character  fname9_to_wfo_fname13*13
	character  data_file_m*150, data_file_l*150, data_file_b*150
c
        character* (*) path_to_metar
        character* (*) path_to_local_data
        character* (*) path_to_buoy_data
c
        integer cnt, minutes_to_wait_for_metars
c        parameter(minutes_to_wait_for_metars=10)
	logical exists
        data exists/.false./
        data cnt/0/
c 
c
c.....	Start here.  Check to see if this is an interactive run.
c
	narg = iargc()
cc      print *,' narg = ', narg
c
c.....  Get the time from the scheduler or from the user if interactive.
c
	if(narg .eq. 0) then
           call get_systime(i4time,filename9,istatus)
	   call cv_i4tim_asc_lp(i4time,atime,istatus)
c
	else
c
 970	   write(6,973)
 973	   format(' Enter input filename (yydddhhmm): ',$)
	   read(5,972) filename9
 972	   format(a9)
	   call i4time_fname_lp(filename9(1:9),i4time,istatus)
	   i4time = i4time / laps_cycle_time * laps_cycle_time
	   call cv_i4tim_asc_lp(i4time, atime, istatus) !find the atime
	endif
c
        call get_directory('lso',outfile,len)
	outfile = outfile(1:len)//filename9(1:9)//'.lso'
cc      outfile = filename9(1:9)//'.lso'
c
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c
        call get_directory('static',dir_s,len)
	ext_s = 'nest7grid'
	var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat,  grid_spacing,istatus)
c
	var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon,  grid_spacing,istatus)
c
	var_s = 'AVG'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      topo, grid_spacing,istatus)
c
c.....	Find east/west and north/south sides of grid (max extension of grid)
c
	grid_east = -999.
	grid_west = 0.
	grid_north = 0.
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
	      store_3ea(i,j) = badflag
	      store_4ea(i,j) = badflag
	      store_6ea(i,j) = badflag
	   enddo !j
c
	   do j=1,3
	      store_2(i,j) = badflag
	      store_7(i,j) = badflag
	      store_2ea(i,j) = badflag
	   enddo !j
c
	   do j=1,4
	      store_1(i,j) = badflag
	      store_3(i,j) = badflag
	      store_5(i,j) = badflag
	      store_5ea(i,j) = badflag
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
c.....  Figure out if the data files are there, paths, etc.
c
        do while(.not. exists .and. 
     &            cnt .lt. minutes_to_wait_for_metars)
c        
	   len_path = index(path_to_METAR,' ') - 1
	   data_file_m = 
     &	      path_to_METAR(1:len_path)//filename9(1:9)// '0100o'
c
	   len_path = index(path_to_local_data,' ') - 1
	   filename13=fname9_to_wfo_fname13(filename9(1:9))
	   data_file_l = 
     &	      path_to_local_data(1:len_path)//filename13
c
 	   len_path = index(path_to_buoy_data,' ') - 1
	   filename13=fname9_to_wfo_fname13(filename9(1:9))
	   data_file_b = 
     &	      path_to_buoy_data(1:len_path)//filename13  
c
	   INQUIRE(FILE=data_file_m,EXIST=exists)
	   if(.not. exists) then
	      filename13=fname9_to_wfo_fname13(filename9(1:9))
	      len_path = index(path_to_METAR,' ') - 1
	      data_file_m = 
     &           path_to_METAR(1:len_path) // filename13
	      len_path = index(path_to_local_data,' ') - 1
	      data_file_l = 
     &           path_to_local_data(1:len_path) // filename13
	      len_path = index(path_to_buoy_data,' ') - 1
	      data_file_b = 
     &           path_to_buoy_data(1:len_path) // filename13
	      INQUIRE(FILE=data_file_m,EXIST=exists)
	      if(.not. exists) then
                 print*,'Waiting for file ', data_file_m
                 call waiting_c(60)
                 cnt = cnt+1               
	      endif
	   endif
	enddo
	if(.not.exists) then
	   print *,' ERROR. File not Found: ', data_file_m
	   stop 'Config error'
        endif
c
c.....  Call the routine that reads the METAR data files, then get
c.....  the data.
c
	print*,'Getting METAR data ', data_file_m
c
        call get_metar_obs(maxobs,maxsta,i4time,data_file_m,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_sao_g,n_sao_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)
c
	if(jstatus .ne. 1) then
	   print *, ' WARNING. Bad status return from GET_METAR_OBS'
	   print *,' '
	endif
c
c.....  Call the routine that reads the LDAD mesonet data files, then get
c.....  the data.
c
	print*,'Getting mesonet data ', data_file_l
c
        call get_local_obs(maxobs,maxsta,i4time,data_file_l,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_local_g,n_local_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, laps_cycle_time, jstatus)
c
	if(jstatus .ne. 1) then
	   print *, ' WARNING. Bad status return from GET_LOCAL_OBS'
	   print *,' '
	endif
c
c.....  Call the routine that reads the Buoy data files, then get
c.....  the data.
c
cc	data_file_b = '/data/fxa/point/maritime/netcdf/' // filename13
	print*,'Getting buoy/ship data ', data_file_b
c
        call get_buoy_obs(maxobs,maxsta,i4time,data_file_b,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_buoy_g,n_buoy_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)
c
	if(jstatus .ne. 1) then
	   print *, ' WARNING. Bad status return from GET_BUOY_OBS'
	   print *,' '
	endif
c
c.....  Count up the obs.
c
	n_obs_g = n_sao_g + n_local_g + n_buoy_g
	n_obs_b = nn
c
c.....  Call the routine to write the LSO file.
c
c
        call write_surface_obs(atime,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,weather,reptype,atype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c.....	That's about it...let's go home.
c
	stop 'Normal completion of OBS_DRIVER'
	end











