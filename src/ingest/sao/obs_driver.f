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
c                               to parameter file
c
c          P. Stamus  03-23-98  Changes for stand-alone QC; LS2 format.
c                     05-01-98  Added soil moisture variables.
c                     08-28-98  Added buoy and LDAD mesonet reads.
c                     09-04-98  Install as LSO, using new 'ls2' format.
c                     09-30-98  Housekeeping changes.
c	              12-03-98  Increase dir_s to 256 characters.
c                     06-21-99  Pass lat/lon, grid size and grid spacing
c                                 to get_ routines to calculate box size.
c                     09-20-99  Add blacklist test (on KFCS w/known bad P).
c                     11-10-99  Upgrade blacklist code.
c	              01-05-00  Hardwire 903 format statement..LINUX compile
c                                 didn't like variable field length.
c
c       Notes:
c         1. When run "operationally", 'obs_driver.x' uses the time from
c            the 'systime.dat' file.  Running 'obs_driver.x -i' allows the
c            user to enter the run time.
c
c******************************************************************************
c
        integer n_madis_dirs
        parameter (n_madis_dirs = 10)

        character*200 path_to_metar
        character*200 path_to_local_data
        character*200 path_to_buoy_data
        character*200 path_to_gps_data
        character*200 path_to_tower_data
        character*200 path_to_synop_data
        character*8   metar_format, c8_project
	character     atime*24, filename9*9
        character*10  madis_dirs(n_madis_dirs)
        character*10 c_mode

        logical l_allow_empty_lso,l_multiple_reports,l_dupe_names

        ISTAT = INIT_TIMER()

        call get_config(istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error returned from get_config'
	   stop
	endif

        call get_grid_dim_xy(nx,ny,istatus)
	if (istatus .ne. 1) then
           write (6,*) 'Error getting horizontal domain dimensions'
	   stop
	endif

        call get_max_stations(maxsta, istatus)
        if(istatus .ne. 1)stop

        call get_laps_cycle_time(laps_cycle_time,istatus)
        if(istatus .ne. 1)stop

!       Get default value for metar_format
        call get_c8_project(c8_project,istatus)
        if(istatus .ne. 1)stop

        metar_format = 'default'
        madis_dirs = ' '

!       Note that metar_format is updated only when fully specified in namelist
        call get_obs_driver_parms(
     1                            path_to_metar
     1                           ,path_to_local_data
     1                           ,path_to_buoy_data
     1                           ,path_to_gps_data
     1                           ,path_to_tower_data
     1                           ,path_to_synop_data
     1                           ,metar_format
     1                           ,minutes_to_wait_for_metars
     1                           ,ick_metar_time
     1                           ,itime_before,itime_after
     1                           ,madis_dirs
     1                           ,itest_madis_qc
     1                           ,maxobs
     1                           ,l_allow_empty_lso,l_multiple_reports
     1                           ,l_dupe_names
     1                           ,local_obs_thresh, i4wait_local_obs_max
     1                           ,n_cycles,nominal_latency
     1                           ,istatus)
        if(istatus .ne. 1)stop

        call GETENV('MODE_VERIF',c_mode)
        if(trim(c_mode) .eq. '1')then
            write(6,*)
     1          ' MODE_VERIF = 1: setting l_multiple_reports to false'
            l_multiple_reports = .false.
        endif

        call s_len(metar_format,len_metar_format)
        if(metar_format(1:len_metar_format) .eq. 'default')then
            metar_format = c8_project
        endif

        call get_systime(i4time_sys,filename9,istatus)

        do i_cycle = n_cycles,1,-1
            I4_elapsed = ishow_timer()

            if(nominal_latency .ge. 0 .and. i_cycle .eq. 1)then
                i4time_now = i4time_now_gg()
                i4_wait = (i4time_sys + nominal_latency) - i4time_now
                if(i4_wait .gt. 0)then
                    write(6,*)' Waiting for the nominal latency time '
     1                       ,i4_wait     
                    CALL sleep(i4_wait)
                endif
                I4_elapsed = ishow_timer()
            endif

            i4time_proc = i4time_sys - ((i_cycle-1) * laps_cycle_time)       
            call make_fnam_lp(i4time_proc,filename9,istatus)
            call cv_i4tim_asc_lp(i4time_proc,atime,istatus)

            write(6,*)
            write(6,*)' calling obs_driver_sub for cycle at ',filename9       

            call obs_driver_sub(  nx,ny
     1                           ,maxobs,laps_cycle_time
     1                           ,path_to_metar
     1                           ,path_to_local_data
     1                           ,path_to_buoy_data
     1                           ,path_to_gps_data
     1                           ,path_to_tower_data
     1                           ,path_to_synop_data
     1                           ,c8_project,metar_format
     1                           ,minutes_to_wait_for_metars
     1                           ,ick_metar_time
     1                           ,itime_before,itime_after
     1                           ,madis_dirs,n_madis_dirs
     1                           ,itest_madis_qc
     1                           ,maxsta
     1                           ,l_allow_empty_lso,l_multiple_reports
     1                           ,l_dupe_names
     1                           ,local_obs_thresh, i4wait_local_obs_max
     1                           ,i4time_proc,filename9,atime
     1                           ,istatus)
        enddo ! i_cycle

        end

        subroutine obs_driver_sub(ni,nj
     1                           ,maxobs,laps_cycle_time
     1                           ,path_to_metar
     1                           ,path_to_local_data
     1                           ,path_to_buoy_data
     1                           ,path_to_gps_data
     1                           ,path_to_tower_data
     1                           ,path_to_synop_data
     1                           ,c8_project,metar_format
     1                           ,minutes_to_wait_for_metars
     1                           ,ick_metar_time
     1                           ,itime_before,itime_after
     1                           ,madis_dirs,n_madis_dirs
     1                           ,itest_madis_qc
     1                           ,maxsta
     1                           ,l_allow_empty_lso,l_multiple_reports
     1                           ,l_dupe_names
     1                           ,local_obs_thresh, i4wait_local_obs_max
     1                           ,i4time_sys,filename9,atime
     1                           ,istatus)
c        
        integer ni, nj, maxsta, maxobs, n_madis_dirs
c
	real    lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
        integer    wmoid(maxsta), jstatus
        integer    dpchar(maxsta), narg, iargc
        integer    local_obs_thresh_switch(n_madis_dirs)
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  weather(maxsta)*25 
        character  reptype(maxsta)*6, atype(maxsta)*6
        character  store_cldamt(maxsta,5)*4
	character  atime*24, outfile*200
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_m*150, data_file_l*150 
c
        character*200 path_to_metar
        character*200 path_to_local_data
        character*200 path_to_madis_data
        character*200 path_to_buoy_data
        character*200 path_to_gps_data
        character*200 path_to_tower_data
        character*200 path_to_synop_data
        character*8   metar_format, c8_project
        character*8   a9_to_a8, a8_time
        character*10  madis_dirs(n_madis_dirs)

!       data madis_dirs /'mesonet','urbanet','hfmetar'/

        logical l_allow_empty_lso,l_string_contains
        logical l_identical_a(maxsta),l_multiple_reports,l_dupe_names
c
c.....	Start here.  
c
        local_obs_thresh_switch = 0
        local_obs_thresh_switch(1) = 1
c
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

        write(6,*)' systime = ',filename9
c
	call get_directory('lso',outfile,len)
	outfile = outfile(1:len)//filename9(1:9)//'.lso'
cc	outfile = filename9(1:9)//'.lso'
c
c.....	Get the LAPS lat/lon and topo data here so we can pass them to the 
c.....	routines that need them.
c
        call get_directory('static',dir_s,len)
	ext_s = 'nest7grid'
	var_s = 'LAT'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lat,  grid_spacing,istatus)
        if(istatus .ne. 1)return
c
	var_s = 'LON'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      lon,  grid_spacing,istatus)
        if(istatus .ne. 1)return
c
	var_s = 'AVG'
        call rd_laps_static(dir_s,ext_s,ni,nj,1,var_s,units,comment,
     &                      topo, grid_spacing,istatus)
        if(istatus .ne. 1)return
c
c.....	Find east/west and north/south sides of grid (max extension of grid)
c
	grid_east = -999.
	grid_west = 0.
	grid_north = 0.
	grid_south = 90.
	do i=1,ni
	  if(lat(i,nj) .gt. grid_north) grid_north = lat(i,nj)
	  if(lat(i,1)  .lt. grid_south) grid_south = lat(i,1)
	enddo !i
	do j=1,nj	
	  if(lon(ni,j) .gt. grid_east) grid_east = lon(ni,j)
	  if(lon(1,j) .lt. grid_west) grid_west = lon(1,j)
	enddo !j
c
c.....	Set up the counters, and zero/blank arrays.
c
	nn = 0
	n_obs_g = 0 ! # of obs in grid (domain)
	n_obs_b = 0 ! # of obs in box that includes perimeter around domain
	n_sao_g = 0
	n_sao_b = 0
	n_local_g = 0
	n_local_b = 0
	n_buoy_g = 0
	n_buoy_b = 0
	n_gps_g = 0
	n_gps_b = 0
c
	do i=1,maxsta
!          Initialize station i
           call init_station(i
     1                      ,stations,provider,weather,reptype,atype      
     1                      ,store_1,store_2,store_3,store_4,store_5
     1                      ,store_6,store_7
     1                      ,store_2ea,store_3ea,store_4ea,store_5ea
     1                      ,store_6ea,dpchar,wmoid
     1                      ,store_cldht,store_cldamt,maxsta,badflag)
	enddo !i

c
c.....  Figure out if the data files are there, paths, etc.
c
        data_file_m = ' '
        data_file_l = ' '

        call s_len(metar_format,len_metar_format)

        if(    metar_format(1:len_metar_format) .eq. 'NIMBUS'
     1    .or. metar_format(1:len_metar_format) .eq. 'WFO'   
     1    .or. metar_format(1:len_metar_format) .eq. 'MADIS'    )then       

!           Select the hourly METAR file best suited to our obs time window
!           Note that an hourly raw file contains obs from 15 before to 45 after
            i4time_midwindow = i4time_sys + 
     1                         (itime_after - itime_before) / 2      
            i4time_metar_file = ((i4time_midwindow+900) / 3600) * 3600

            call make_fnam_lp(i4time_metar_file,a9time_metar_file
     1                       ,istatus)
            if(istatus .ne. 1)return

            if(metar_format(1:len_metar_format) .eq. 'NIMBUS')then
                len_path = index(path_to_METAR,' ') - 1
	        data_file_m = path_to_METAR(1:len_path)
     1                        //a9time_metar_file// '0100o'       
c        
            elseif(metar_format(1:len_metar_format) .eq. 'WFO'
     1    .or.     metar_format(1:len_metar_format) .eq. 'MADIS')then
                filename13=fname9_to_wfo_fname13(a9time_metar_file)       

                len_path = index(path_to_METAR,' ') - 1
                data_file_m = path_to_METAR(1:len_path)//filename13       

            else
                write(6,*)' ERROR: unknown metar format ',metar_format  
                istatus = 0
                return

            endif

        elseif(metar_format(1:len_metar_format) .eq. 'CWB')then
            continue

        elseif(metar_format(1:len_metar_format) .eq. 'AFWA')then
            continue

        else
            write(6,*)' ERROR: unknown metar format ',metar_format
            istatus=0
            return
       
        endif ! FSL format
c
c.....  Call the routine that reads the METAR data files, then get
c.....  the data.
c
	print*,'Getting METAR data ', data_file_m
c
        if(metar_format(1:len_metar_format) .ne. 'AFWA')then
           call get_metar_obs(maxobs,maxsta,i4time_sys,
     &                        path_to_metar,metar_format,   
     &                        minutes_to_wait_for_metars,
     &                        ick_metar_time,itime_before,itime_after,
     &                        grid_east,grid_west,grid_north,grid_south,       
     &                        lat,lon,ni,nj,grid_spacing,
     &                        nn,n_sao_g,n_sao_b,stations,
     &                        reptype,atype,weather,wmoid,
     &                        store_1,store_2,store_2ea,
     &                        store_3,store_3ea,store_4,store_4ea,
     &                        store_5,store_5ea,store_6,store_6ea,
     &                        store_7,store_cldht,store_cldamt,
     &                        provider, jstatus)
	   if(jstatus .ne. 1) then
	      print *, ' WARNING: Bad status return from GET_METAR_OBS'       
	      print *,' '
	   endif

        else
           call get_sao_obs_af(filename9,path_to_metar,maxsta,
     &                        grid_east,grid_west,grid_north,grid_south,    
     &                        nn,lat,lon,ni,nj,grid_spacing,
     &                        n_sao_g,n_sao_b,stations,
     &                        reptype,atype,weather,wmoid,
     &                        store_1,store_2,store_2ea,
     &                        store_3,store_3ea,store_4,store_4ea,
     &                        store_5,store_5ea,store_6,store_6ea,
     &                        store_7,store_cldht,store_cldamt,
     &                        provider,jstatus)

	   if(jstatus .ne. 1) then
	      print *, ' WARNING. Bad status return from GET_SAO_OBS_AF'       
	      print *,' '
	   endif

        endif

        if(nn .gt. maxsta)then
           write(6,*)' ERROR: nn > maxsta ',nn,maxsta
           return
        endif
c
c.....  Call the routine that reads the mesonet data files, then get the data.
c
        write(6,*)
	write(6,*)'Getting Mesonet Data...'
c
        if(metar_format(1:len_metar_format) .ne. 'CWB')then ! LDAD Netcdf

            if(l_string_contains(path_to_local_data,'madis',istatus)
     1                                                             )then      
	        call s_len(path_to_local_data,len_path)

                n_local_g = 0

                do imadis = 1,n_madis_dirs

                  call s_len(madis_dirs(imadis),len_madis)
                  if(len_madis .gt. 0)then
               
	            path_to_madis_data = path_to_local_data(1:len_path)      
     1                           //'/'//madis_dirs(imadis)(1:len_madis)
     1                           //'/'//'netCDF/'

!                   Wait for 'mesonet' with no waiting for 'urbanet'
                    local_obs_thresh_madis = local_obs_thresh * 
     1              local_obs_thresh_switch(imadis)

                    write(6,*)' MADIS case: ',path_to_madis_data

                    if(madis_dirs(imadis)(1:len_madis) .ne. 'hydro')then       
                      call get_local_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_madis_data,metar_format,
     &                      itime_before,itime_after,
     &                      itest_madis_qc,l_multiple_reports,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_local_gx,n_local_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, laps_cycle_time, 
     &                      local_obs_thresh_madis, 
     &                      i4wait_local_obs_max, 
     &                      jstatus)       
                    else
                      call get_hydro_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_madis_data,metar_format,
     &                      itime_before,itime_after,
     &                      itest_madis_qc,l_multiple_reports,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_local_gx,n_local_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,!store_2,store_2ea,
!    &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
!    &                      store_7,store_cldht,store_cldamt,
     &                      provider, laps_cycle_time, 
     &                      local_obs_thresh_madis, 
     &                      i4wait_local_obs_max, 
     &                      jstatus)       
                    endif

	            if(jstatus .ne. 1) then
	               print *, 
     1                  ' WARNING. Bad status return from GET_LOCAL_...'       
	               print *,' '
	            endif

                    n_local_g = n_local_g + n_local_gx

                    if(nn .gt. maxsta)then
                       write(6,*)' ERROR: nn > maxsta ',nn,maxsta
                       return
                    else
                       write(6,*)' Total obs so far is ',nn
                    endif

                  endif ! directory length > 0

                enddo ! imadis

            else 
                write(6,*)' non-MADIS case'
                call get_local_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_local_data,metar_format,
     &                      itime_before,itime_after,
     &                      itest_madis_qc,l_multiple_reports,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_local_g,n_local_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, laps_cycle_time, 
     &                      local_obs_thresh, i4wait_local_obs_max, 
     &                      jstatus)       

	        if(jstatus .ne. 1) then
	           print *, 
     1                  ' WARNING. Bad status return from GET_LOCAL_...'       
	           print *,' '
	        endif

                if(nn .gt. maxsta)then
                   write(6,*)' ERROR: nn > maxsta ',nn,maxsta
                   return
                endif

            endif

        else ! CWB data
            call get_local_cwb(maxobs,maxsta,i4time_sys,
     &                      path_to_local_data,metar_format,
     &                      itime_before,itime_after,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_local_g,n_local_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, laps_cycle_time, jstatus)

	    if(jstatus .ne. 1) then
	       print *, ' WARNING. Bad status return from GET_LOCAL_...'
	       print *,' '
	    endif

            if(nn .gt. maxsta)then
               write(6,*)' ERROR: nn > maxsta ',nn,maxsta
               return
            endif

        endif
c
c
c.....  Call the routine that reads the Buoy data files, then get
c.....  the data.
c
c
        call get_maritime_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_buoy_data,metar_format,
     &                      itime_before,itime_after,
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
	   print *, ' WARNING. Bad status return from GET_MARITIME_OBS'       
	   print *,' '
	endif

        if(nn .gt. maxsta)then
           write(6,*)' ERROR: nn > maxsta ',nn,maxsta
           return
        endif
c
c.....  Call the routine that reads the GPS data files, then get
c.....  the data.
c
        if(.true.)then

	    print*,'Getting GPS data...'
c
            call get_gps_obs(maxobs,maxsta,i4time_sys,
     &                      path_to_gps_data,metar_format,
     &                      itime_before,itime_after,
     &                      grid_east,grid_west,grid_north,grid_south,       
     &                      lat,lon,ni,nj,grid_spacing,
     &                      nn,n_gps_g,n_gps_b,stations,
     &                      reptype,atype,weather,wmoid,
     &                      store_1,store_2,store_2ea,
     &                      store_3,store_3ea,store_4,store_4ea,
     &                      store_5,store_5ea,store_6,store_6ea,
     &                      store_7,store_cldht,store_cldamt,
     &                      provider, jstatus)
c
	    if(jstatus .ne. 1) then
	       print *, ' WARNING. Bad status return from GET_GPS_OBS'
	       print *,' '
	    endif

            if(nn .gt. maxsta)then
               write(6,*)' ERROR: nn > maxsta ',nn,maxsta
               return
            endif

        endif

        if(c8_project(1:3) .eq. 'RSA')then ! Tower soil moisture
	    print*,'Getting Tower soil moisture data...'
            call tower_sfc_driver(  maxsta,i4time_sys                     ! I
     1                             ,path_to_tower_data                    ! I
     1                             ,lat,lon,ni,nj,grid_spacing            ! I
     1                             ,laps_cycle_time                       ! I
     1                             ,itime_before,itime_after              ! I
     1                             ,nn,n_local_g,n_local_b,stations       ! I/O
     1                             ,store_1,store_2,store_2ea             ! O
     1                             ,store_3,store_3ea,store_4,store_4ea   ! O    
     1                             ,store_5,store_5ea,store_6,store_6ea   ! O
     1                             ,store_7,store_cldht,store_cldamt      ! O
     1                             ,provider,istatus)                     ! O
        endif
c
c.....  Count up the obs.
c
	n_obs_g = n_sao_g + n_local_g + n_buoy_g + n_gps_g
	n_obs_b = nn

c       Call subroutine to blacklist the stations in the "store" arrays
        call apply_blacklist(      maxsta,n_obs_b,stations,provider
     1                            ,store_1,store_2,store_3
     1                            ,store_4,store_5,store_6
     1                            ,store_7,badflag)

        I4_elapsed = ishow_timer()

c       Call subroutine to check for duplicate obs
        if(l_multiple_reports .eqv. .true.)then
            l_dupe_names = .false.
        endif

        call check_for_dupes(      maxsta,n_obs_b,stations
     1                            ,store_1,store_2,store_3
     1                            ,store_4,store_5,store_6
     1                            ,store_7,badflag,l_identical_a
     1                            ,l_dupe_names)

        I4_elapsed = ishow_timer()

!       Remove identical stations by calling 'init_station'
        do i = 1,n_obs_b
            if(l_identical_a(i) .and. .true.)then
                write(6,*)' Removing identical station',i
                call init_station(i
     1                      ,stations,provider,weather,reptype,atype      
     1                      ,store_1,store_2,store_3,store_4,store_5
     1                      ,store_6,store_7
     1                      ,store_2ea,store_3ea,store_4ea,store_5ea
     1                      ,store_6ea,dpchar,wmoid
     1                      ,store_cldht,store_cldamt,maxsta,badflag)
            endif
        enddo
c
!       Final QC check 
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

!       Replace blank/UNK station names with wmoid if possible, else set blank
        iblank = 0
        do i = 1,n_obs_b
            call s_len(stations(i),lensta)
            if(lensta .eq. 0 .or. stations(i)(1:3) .eq. 'UNK')then
                if(wmoid(i) .ne. ibadflag .and. wmoid(i) .ne. 0)then
                    write(stations(i),511,err=512)wmoid(i)
 511		    format(i8)
 512		    continue
                else
                    stations(i) = 'UNK                 '
                    iblank = iblank + 1
                endif
            endif
        enddo

        if(iblank .gt. 0)then
            write(6,*)' Warning: number of UNK stanames = ',iblank       
        endif

!       Count pressure obs
        nalt = 0
        nstp = 0
        nmsl = 0
        nmsl_only = 0
        nalt_and_msl = 0
        nalt_or_msl = 0
        nalt_or_stp = 0
        sum_alt_or_stp = 0.
        do i = 1,n_obs_b
            if(store_4(i,1) .ne. badflag)then
                nalt = nalt+1
            endif
            if(store_4(i,2) .ne. badflag)then
                nstp = nstp+1
            endif
            if(store_4(i,3) .ne. badflag)then
                nmsl = nmsl+1
            endif
            if(store_4(i,3) .ne. badflag .and.
     1         store_4(i,1) .eq. badflag .and.                         
     1         store_4(i,2) .eq. badflag            )then
                nmsl_only = nmsl_only+1
            endif
            if(store_4(i,1) .ne. badflag .and. 
     1         store_4(i,3) .ne. badflag            )then
                nalt_and_msl = nalt_and_msl+1
            endif
            if(store_4(i,1) .ne. badflag .or. 
     1         store_4(i,3) .ne. badflag            )then
                nalt_or_msl = nalt_or_msl+1
            endif
            if(store_4(i,1) .ne. badflag .or. 
     1         store_4(i,2) .ne. badflag            )then
                nalt_or_stp = nalt_or_stp+1
                sum_alt_or_stp = sum_alt_or_stp + store_1(i,3)
            endif
        enddo ! i

        if(nalt_or_stp .gt. 0)then
            ave_elev_pres = sum_alt_or_stp / float(nalt_or_stp)
        else
            ave_elev_pres = badflag
        endif

        write(6,*)
        write(6,*)' Checking pressure reports in box...'
        write(6,*)' # of stations reporting altimeter         ',nalt
        write(6,*)' # of stations reporting station pressure  ',nstp
        write(6,*)' # of stations reporting MSL pressure      ',nmsl
        write(6,*)' # of stations reporting MSL pressure only '
     1                                                 ,nmsl_only     
        write(6,*)' # of stations reporting altimeter and MSLP'
     1                                                 ,nalt_and_msl    
        write(6,*)' # of stations reporting altimeter or MSLP '
     1                                                 ,nalt_or_msl    
        write(6,*)' # of stations reporting altimeter or STNP '
     1                                                 ,nalt_or_stp    
        write(6,*)' Mean elevation of alt/stnp reports = ' 
     1                                                 ,ave_elev_pres       

!       Check for no obs
        if(nn .eq. 0 .and. .not. l_allow_empty_lso)then
            write(6,*)' WARNING: no LSO written due to no obs'
            return
        endif

c
c.....  Call the routine to write the LSO file.
c

        print *
	print *,'  Writing LSO file, # of obs (in box) = ',n_obs_b
c
        call write_surface_obs(atime,outfile,n_obs_g,
     &    n_obs_b,wmoid,stations,provider,weather,reptype,atype,
     &    store_1,store_2,store_3,store_4,store_5,store_6,store_7,
     &    store_2ea,store_3ea,store_4ea,store_5ea,store_6ea,
     &    store_cldamt,store_cldht,maxsta,jstatus)
c
c.....	That's about it...let's go home.
c
	write(6,*)' Normal completion of OBS_DRIVER'

	end


 
       subroutine get_obs_driver_parms(path_to_metar
     1                         ,path_to_local_data
     1                         ,path_to_buoy_data
     1                         ,path_to_gps_data
     1                         ,path_to_tower_data
     1                         ,path_to_synop_data
     1                         ,metar_format
     1                         ,minutes_to_wait_for_metars
     1                         ,ick_metar_time
     1                         ,itime_before,itime_after
     1                         ,madis_dirs
     1                         ,itest_madis_qc
     1                         ,maxobs
     1                         ,l_allow_empty_lso,l_multiple_reports
     1                         ,l_dupe_names
     1                         ,local_obs_thresh, i4wait_local_obs_max
     1                         ,n_cycles,nominal_latency
     1                         ,istatus)

       integer n_madis_dirs
       parameter (n_madis_dirs = 10)

       character*200 path_to_metar
       character*200 path_to_local_data
       character*200 path_to_buoy_data
       character*200 path_to_gps_data
       character*200 path_to_tower_data
       character*200 path_to_synop_data
       character*8   metar_format
       character*10  madis_dirs(n_madis_dirs)
       logical l_allow_empty_lso,l_multiple_reports,l_dupe_names

       namelist /obs_driver_nl/ path_to_metar
     1                         ,path_to_local_data
     1                         ,path_to_buoy_data
     1                         ,path_to_gps_data
     1                         ,path_to_tower_data
     1                         ,path_to_synop_data
     1                         ,metar_format
     1                         ,minutes_to_wait_for_metars
     1                         ,ick_metar_time
     1                         ,itime_before,itime_after
     1                         ,madis_dirs
     1                         ,l_allow_empty_lso
     1                         ,l_multiple_reports
     1                         ,l_dupe_names
     1                         ,maxobs
     1                         ,local_obs_thresh
     1                         ,i4wait_local_obs_max
     1                         ,n_cycles,nominal_latency
     1                         ,itest_madis_qc

       character*150 static_dir,filename
 
       local_obs_thresh = 0 ! default value
       i4wait_local_obs_max = 0   ! default value

       call get_directory('static',static_dir,len_dir)

       filename = static_dir(1:len_dir)//'/obs_driver.nl'
 
       open(1,file=filename,status='old',err=900)
       read(1,obs_driver_nl,err=901)
       close(1)

       istatus = 1
       return

  900  print*,'error opening file ',filename
       istatus = 0
       return

  901  print*,'error reading obs_driver_nl in ',filename
       write(*,obs_driver_nl)
       istatus = 0
       return
 
       end


        subroutine apply_blacklist(maxsta,n_obs_b,stations,provider
     1                            ,store_1,store_2,store_3
     1                            ,store_4,store_5,store_6
     1                            ,store_7,badflag)

	integer    num_varb(maxsta), max_bvar
	parameter  (max_bvar=20)  !max number of variables to blacklist...
                                  !change 903 format statement if you make
                                  !this greater than 20.

	character  dir_b*256, black_path*256, stations_b(maxsta)*20
	character  var_b(maxsta,max_bvar)*3 
        character  stations(maxsta)*20,provider(maxsta)*11
	real    store_1(maxsta,4), 
     &          store_2(maxsta,3), 
     &          store_3(maxsta,4), 
     &          store_4(maxsta,5), 
     &          store_5(maxsta,4), 
     &          store_6(maxsta,5), 
     &          store_7(maxsta,3)

	logical exists, l_provider, l_match
        data exists/.false./
c
c.....  Check for a blacklist file.  If one exists, read it
c.....  and flag the variables listed for each blacklist station
c.....  as bad.  Otherwise, skip over this section and write what
c.....  we have.
c
	call get_directory('static',dir_b,len)
	black_path = dir_b(1:len) // 'Blacklist.dat'
cc	black_path = './Blacklist.dat'
	call s_len(black_path, len)
c
	print *,' '
	print *,'  Looking for blacklist file...'
        INQUIRE(FILE=black_path(1:len), EXIST=exists)
        if( .not. exists ) then
	   print *,'    No blacklist file found.'
	   go to 505		!no blacklist file...skip this stuff
	endif
	print *,'  Found one.  Checking blacklist file.'
c
c.....  Have blacklist file...read it.
c
	open(11, file=black_path, status='unknown')
	read(11,901) num_black
 901	format(1x,i5)
c
	do i=1,num_black
	   read(11,903) stations_b(i), num_varb(i), 
     &                  (var_b(i,j), j=1,max_bvar)
	enddo !i
 903	format(1x,a20,2x,i3,20(2x,a3))
        close(11)
c
c.....  Have blacklist info.  Now flag the var_b's at each station_b as bad.
c
	do ibl=1,num_black ! loop over the blacklist stations

           if(num_varb(ibl) .gt. 100)then
              num_varb(ibl) = num_varb(ibl) - 100
              l_provider = .true.
           else
              l_provider = .false.
           endif           

	   do i=1,n_obs_b ! loop over all stations

              if(l_provider)then
	          if(provider(i) .eq. stations_b(ibl)) then ! found a match
                      l_match = .true.
                  else
                      l_match = .false.
                  endif

              else
	          if(stations(i) .eq. stations_b(ibl)) then ! found a match
                      l_match = .true.
                  else
                      l_match = .false.
                  endif

              endif

              if(l_match)then

                 write(6,911)stations(i),provider(i),l_provider
 911		 format(a,a,l2)

		 do j=1,num_varb(ibl)

		    if( var_b(ibl,j) .eq. 'ALL' ) then
		       store_2(i,1) = badflag ! temperature
		       store_2(i,2) = badflag ! dew point
		       store_2(i,3) = badflag ! relative humidity
		       store_3(i,1) = badflag ! wind dir
		       store_3(i,2) = badflag ! wind speed
		       store_3(i,3) = badflag ! wind gust dir
		       store_3(i,4) = badflag ! wind gust speed
		       store_4(i,1) = badflag ! altimeter
		       store_4(i,2) = badflag ! station pressure
		       store_4(i,3) = badflag ! MSL pressure
		       store_5(i,1) = badflag ! visibility
		       store_5(i,2) = badflag ! solar
		       store_5(i,3) = badflag ! soil/water temp
		       store_5(i,4) = badflag ! soil moisture
		       store_6(i,1) = badflag !  1-h precip
		       store_6(i,2) = badflag !  3-h precip
		       store_6(i,3) = badflag !  6-h precip
		       store_6(i,4) = badflag ! 12-h precip
		       store_6(i,5) = badflag ! snow depth
		       store_7(i,1) = 0       ! clouds (set num cld layers to 0)
c
		    elseif( var_b(ibl,j) .eq. 'TMP' ) then
		       store_2(i,1) = badflag ! temperature
		    elseif( var_b(ibl,j) .eq. 'DEW' ) then
		       store_2(i,2) = badflag ! dew point
		    elseif( var_b(ibl,j) .eq. 'HUM' ) then
		       store_2(i,3) = badflag ! relative humidity
c
		    elseif( var_b(ibl,j) .eq. 'WND' ) then
		       store_3(i,1) = badflag ! wind dir
		       store_3(i,2) = badflag ! wind speed
		       store_3(i,3) = badflag ! wind gust dir
		       store_3(i,4) = badflag ! wind gust speed
c
		    elseif( var_b(ibl,j) .eq. 'ALT' ) then
		       store_4(i,1) = badflag ! altimeter
		    elseif( var_b(ibl,j) .eq. 'STP' ) then
		       store_4(i,2) = badflag ! station pressure
		    elseif( var_b(ibl,j) .eq. 'MSL' ) then
		       store_4(i,3) = badflag ! MSL pressure
c
		    elseif( var_b(ibl,j) .eq. 'VIS' ) then
		       store_5(i,1) = badflag ! visibility
		    elseif( var_b(ibl,j) .eq. 'SOL' ) then
		       store_5(i,2) = badflag ! solar
		    elseif( var_b(ibl,j) .eq. 'SWT' ) then
		       store_5(i,3) = badflag ! soil/water temp
		    elseif( var_b(ibl,j) .eq. 'SWM' ) then
		       store_5(i,4) = badflag ! soil moisture
c
		    elseif( var_b(ibl,j) .eq. 'PCP' ) then
		       store_6(i,1) = badflag !  1-h precip
		       store_6(i,2) = badflag !  3-h precip
		       store_6(i,3) = badflag !  6-h precip
		       store_6(i,4) = badflag ! 12-h precip
		    elseif( var_b(ibl,j) .eq. 'SNW' ) then
		       store_6(i,5) = badflag ! snow depth
c
		    elseif( var_b(ibl,j) .eq. 'CLD' ) then
		       store_7(i,1) = 0       ! clouds (set num cld layers to 0)
c
		    else
		       print *,' '
		       print *,' WARNING. Invalid blacklist variable: ',
     &                   var_b(ibl,j),' for station ',stations_b(ibl)

		    endif
c
		 enddo !j
	      endif ! match found
	   enddo ! i (station loop)
	enddo ! ibl (blacklist loop)

	print *,' '
	print *,'  Done with blacklisting.'
 505	continue

        return
        end


        subroutine check_for_dupes(maxsta,n_obs_b,stations
     1                            ,store_1,store_2,store_3
     1                            ,store_4,store_5,store_6
     1                            ,store_7,badflag,l_identical_a
     1                            ,l_dupe_names)

	character  stations(maxsta)*20
        logical l_identical_a(maxsta),l_identical,l_dupe_names

	real    store_1(maxsta,4), 
     &          store_2(maxsta,3), 
     &          store_3(maxsta,4), 
     &          store_4(maxsta,5), 
     &          store_5(maxsta,4), 
     &          store_6(maxsta,5), 
     &          store_7(maxsta,3)

        write(6,*)
        write(6,*)' Checking for identical stations...'

        l_identical_a = .false.

        do i=1,n_obs_b-1 ! loop over all stations
            if(.not. l_identical_a(i))then
	        do j=i+1,n_obs_b ! loop over all stations
                    l_identical = .true.

                    do k = 1,4
                        if(store_1(i,k) .ne. store_1(j,k)) then ! different loc/time
                            l_identical = .false.
                        endif
	            enddo ! k

                    if(l_identical)then
                        write(6,*)' Location/Time are identical: '
     1                       ,i,stations(i),j,stations(j)

                        if(store_2(i,1) .eq. store_2(j,1))then
                            write(6,*)' Temp is also identical'
     1                               ,store_2(i,1),store_2(j,1)
                        else
                            write(6,*)' Temp is not identical'
     1                               ,store_2(i,1),store_2(j,1)
                            l_identical = .false.
                        endif
                    endif ! l_identical

                    call s_len(stations(i),leni)
!                   call s_len(stations(j),lenj)
                    if(l_dupe_names .and.
     1                 stations(i) .eq. stations(j) .and. 
     1                 leni .gt. 0                  .and. 
     1                 stations(i)(1:3) .ne. 'UNK'       )then
                        write(6,*)' Names are identical: '
     1                           ,i,j,stations(i)     
                        l_identical = .true.
                    endif

                    if(l_identical)then
                        l_identical_a(j) = .true.
                    endif

	        enddo ! j
            endif ! not already flagged as identical

	enddo ! i

        n_identical = 0

        do i=1,n_obs_b ! loop over all stations
            if(l_identical_a(i))then
                n_identical = n_identical + 1
                write(6,*)' Identical station ',i
            endif
        enddo ! i

	print *,' '
	print *,'  Duplicate check: # of identical stations = '
     1         ,n_identical     

        return
        end


        subroutine init_station(i
     1                      ,stations,provider,weather,reptype,atype      
     1                      ,store_1,store_2,store_3,store_4,store_5
     1                      ,store_6,store_7
     1                      ,store_2ea,store_3ea,store_4ea,store_5ea
     1                      ,store_6ea,dpchar,wmoid
     1                      ,store_cldht,store_cldamt,maxsta,badflag)

	real  store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2),
     &          store_4(maxsta,5), store_4ea(maxsta,2),
     &          store_5(maxsta,4), store_5ea(maxsta,4),
     &          store_6(maxsta,5), store_6ea(maxsta,2),
     &          store_7(maxsta,3),
     &          store_cldht(maxsta,5)
c
        integer    wmoid(maxsta)
        integer    dpchar(maxsta)
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  weather(maxsta)*25 
        character  reptype(maxsta)*6, atype(maxsta)*6
        character  store_cldamt(maxsta,5)*4

        if(i .le. 0 .or. i .gt. maxsta)then
            write(6,*)' ERROR: init_station index out of bounds ',i       
            stop
        endif            

        dpchar(i) = 0
        wmoid(i) = 0

	stations(i) = '                    '
	provider(i) = '           '
	weather(i)  = '                         '
	reptype(i)  = '      '
	atype(i)    = '      '
c
	do j=1,2
	    store_3ea(i,j) = 0. ! badflag
	    store_4ea(i,j) = 0. ! badflag
	    store_6ea(i,j) = 0. ! badflag
	enddo !j
c
	do j=1,3
	    store_2(i,j) = badflag
	    store_7(i,j) = badflag
	    store_2ea(i,j) = 0. ! badflag
	enddo !j
c
	do j=1,4
	    store_1(i,j) = badflag
	    store_3(i,j) = badflag
	    store_5(i,j) = badflag
	    store_5ea(i,j) = 0. ! badflag
	enddo !j
c
	do j=1,5
	    store_4(i,j) = badflag
	    store_6(i,j) = badflag
	    store_cldht(i,j) = badflag
	    store_cldamt(i,j) = '    '
	enddo !j

        return
        end
