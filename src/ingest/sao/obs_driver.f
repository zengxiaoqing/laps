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
c      J. Edwards     07-14-97  Made dynamic and moved data paths 
c                               to nest7grid.parms
c
c       Notes:
c         1. When run "operationally", 'obs_driver.x' uses the time from
c            '../sched/systime.dat'.  Running 'obs_driver.x -i' allows the
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
        call obs_driver_sub(NX_L_CMN,NY_L_CMN,maxstations_cmn
     +             ,maxobs_cmn,path_to_metar_data_cmn
     +             ,path_to_local_data_cmn,laps_cycle_time_cmn)

        end

        subroutine obs_driver_sub(ni,nj,maxsta,maxobs
     +          ,path_to_metar,path_to_local_data,laps_cycle_time)
        
        include 'surface_obs.inc'
        integer ni,nj,maxsta,maxobs
        character* (*) path_to_metar
        character* (*) path_to_local_data

	parameter (num_var = 20)        !Number of variables stored
c
	real*4 store(maxsta,num_var), store_hgt(maxsta,5)
	real*4 lat(ni,nj), lon(ni,nj), topo(ni,nj)
c
	integer*4 jstatus, grid_spacing 
	integer narg, iargc
c
	character store_amt(maxsta,5)*4, store_emv(maxsta,5)*1
	character stations(maxsta)*3, wx(maxsta)*8, obstype(maxsta)*8
c
	character atime*24, outfile*200
	character dir_s*50,ext_s*31,units*10,comment*125,var_s*3
c
	character filename*13
        character fname9_to_wfo_fname13*13
	character data_file_m*80, data_file_c*80
	logical exists
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
           call get_systime(i4time,filename,istatus)
c           call get_directory('etc',outfile,len)
c           open(11,file=outfile(1:len)//'systime.dat',status='unknown')
ccc	   open(11,file='../sched/systime.dat',status='unknown')
c	   read(11,21) i4time
c 21	   format(1x,i11)
c	   read(11,22) filename
c 22	   format(1x,a9)
c	   close(11)
	   call cv_i4tim_asc_lp(i4time,atime,istatus)
c
	else
c
 970	   write(6,973)
 973	   format(' Enter input filename (yydddhhmm): ',$)
	   read(5,972) filename
 972	   format(a9)
	   call i4time_fname_lp(filename(1:9),i4time,istatus)
	   i4time = i4time / laps_cycle_time * laps_cycle_time
	   call cv_i4tim_asc_lp(i4time, atime, istatus) !find the atime
	endif
c
cc	outfile = filename//'.lso'
        call get_directory('lso',outfile,len)
	outfile = outfile(1:len)//filename(1:9)//'.lso'

cc	outfile = '/home/peaks1/stamus/laps/obs/' // filename // '.lso'
c
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c
        call get_directory('static',dir_s,len)
c	dir_s = '../static/' 
cc	dir_s = '/data/laps/nest7grid/static/' 
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
c.....	Set up the counters.
c
	nn = 0
	n_obs_b = 0
	n_obs_g = 0
	n_obs_pos_b = 0
	n_obs_pos_g = 0
	n_meso_old = 0
	n_meso_pos = 0
	n_cdot = 0
c
c.....  Call the routine that reads the METAR data files, then get
c.....  the METAR data.
c

        
c	if(idata_config .eq. 1) then   ! /public at FSL
        len_path = index(path_to_METAR,' ') - 1
	data_file_m = 
     1    path_to_METAR(1:len_path)//filename(1:9)// '0100o'
        len_path = index(path_to_local_data,' ') - 1
	data_file_c = 
     &      path_to_local_data(1:len_path)//filename(1:9)//'0015r'
        INQUIRE(FILE=data_file_m,EXIST=exists)
        if(.not. exists) then
	  filename=fname9_to_wfo_fname13(filename(1:9))
	  len_path = index(path_to_METAR,' ') - 1
	  data_file_m = 
     &        path_to_METAR(1:len_path) // filename
	  len_path = index(path_to_local_data,' ') - 1
	  data_file_c = 
     &        path_to_local_data(1:len_path) // filename
	  INQUIRE(FILE=data_file_m,EXIST=exists)
	  if(.not. exists) then
	     print *,' ERROR. File not Found: ',data_file_m
	     stop 'Config error'
	  endif
	endif
	print*,'Getting surface data ',data_file_m
c	elseif(idata_config .eq. 2) then ! WFO-adv data
c	   data_file_m = 
c     &        path_to_METAR(1:len_path) // filename_wfo
c	else
c
	call get_metar_obs(maxobs,maxsta,i4time,data_file_m,
     &                   grid_east,grid_west,grid_north,
     &                   grid_south,nn,n_sao_g,n_sao_pos_g,n_sao_b,
     &                   n_sao_pos_b,stations,store,wx,obstype,
     &                   store_emv,store_amt,store_hgt,
     &                   num_var,jstatus) 
c
c.....  Set up the local_data filename, then call the routine that reads 
c.....  the local_data data files.
c
        

        
	print*, 'Getting local data',data_file_c

c
	call get_cdot_obs(maxobs,maxsta,i4time,data_file_c,
     &                   grid_east,grid_west,grid_north,
     &                   grid_south,nn,n_cdot,
     &                   stations,store,wx,obstype,
     &                   store_emv,store_amt,store_hgt,
     &                   num_var,jstatus) 
c
c.....  Finish the obs counts, then call the routine to write the LSO file.
c
	n_obs_g = n_sao_g + n_cdot
	n_obs_b = nn
	n_obs_pos_g = n_sao_pos_g + n_meso_pos
	n_obs_pos_b = n_sao_pos_b + n_meso_pos
c
	call write_surface_obs(atime,outfile,n_meso_old,n_meso_pos,
     &    n_sao_g,n_sao_pos_g,n_sao_b,n_sao_pos_b,n_obs_g,n_obs_pos_g,
     &    n_obs_b,n_obs_pos_b,stations,store,wx,obstype,store_emv,
     &    store_amt,store_hgt,maxsta,num_var,badflag,jstatus)

c
c.....	That's about it...let's go home.
c
	stop 'Normal completion of OBS_DRIVER'
	end











