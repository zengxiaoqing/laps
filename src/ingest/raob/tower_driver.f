
        subroutine tower_driver_sub(ni,nj
     1                           ,maxobs,laps_cycle_time
     1                           ,path_to_local_data
     1                           ,itime_before,itime_after
     1                           ,maxsta
     1                           ,istatus)
c        
        integer ni, nj, maxsta, maxobs 
c
	real    lat(ni,nj), lon(ni,nj), topo(ni,nj)
	real    store_1(maxsta,4), 
     &          store_2(maxsta,3), store_2ea(maxsta,3),
     &          store_3(maxsta,4), store_3ea(maxsta,2) 
!    &          store_4(maxsta,5), store_4ea(maxsta,2),
!    &          store_5(maxsta,4), store_5ea(maxsta,4),
!    &          store_6(maxsta,5), store_6ea(maxsta,2),
!    &          store_7(maxsta,3),
!    &          store_cldht(maxsta,5)
c
        integer    wmoid(maxsta), jstatus, grid_spacing 
        integer    dpchar(maxsta), narg, iargc
c
        character  stations(maxsta)*20, provider(maxsta)*11
        character  weather(maxsta)*25 
        character  reptype(maxsta)*6, atype(maxsta)*6
        character  store_cldamt(maxsta,5)*4
	character  atime*24, outfile*200
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_l*150
c
        character*200 path_to_metar
        character*200 path_to_local_data
        character*200 path_to_buoy_data
        character*200 path_to_gps_data
        character*8   metar_format
        character*8   a9_to_a8, a8_time
c
        integer cnt
        data cnt/0/
c 
c
c.....	Start here.  
c
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

c.....  Get the time from the scheduler or from the user if interactive.
c
        call get_systime(i4time_sys,filename9,istatus)
	call cv_i4tim_asc_lp(i4time_sys,atime,istatus)
c
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
	n_obs_g = 0
	n_obs_b = 0
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
	   stations(i) = '                    '
	   provider(i) = '           '
	   weather(i)  = '                         '
	   reptype(i)  = '      '
	   atype(i)    = '      '
c
	   do j=1,2
	      store_3ea(i,j) = badflag
!	      store_4ea(i,j) = badflag
!	      store_6ea(i,j) = badflag
	   enddo !j
c
	   do j=1,3
	      store_2(i,j) = badflag
!	      store_7(i,j) = badflag
	      store_2ea(i,j) = badflag
	   enddo !j
c
	   do j=1,4
	      store_1(i,j) = badflag
	      store_3(i,j) = badflag
!	      store_5(i,j) = badflag
!	      store_5ea(i,j) = badflag
	   enddo !j
c
!	   do j=1,5
!	      store_4(i,j) = badflag
!	      store_6(i,j) = badflag
!	      store_cldht(i,j) = badflag
!	      store_cldamt(i,j) = '    '
!	   enddo !j
	enddo !i

c
c.....  Figure out if the data files are there, paths, etc.
c
        data_file_l = ' '
        metar_format = 'WFO'

        call s_len(metar_format,len_metar_format)

        if(    metar_format(1:len_metar_format) .eq. 'NIMBUS'
     1    .or. metar_format(1:len_metar_format) .eq. 'WFO'        )then

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
c        
            elseif(metar_format(1:len_metar_format) .eq. 'WFO')then
                filename13=fname9_to_wfo_fname13(a9time_metar_file)       

            else
                write(6,*)' ERROR'
                stop

            endif

        else
            write(6,*)' ERROR: unknown metar format ',metar_format
            stop
       
        endif ! FSL format

c
c.....  Call the routine that reads the mesonet data files, then get the data.
c
        write(6,*)
	write(6,*)'Getting Mesonet Tower Data...'
c
        call get_local_towerobs(maxobs,maxsta,i4time_sys,
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
c
c.....  Count up the obs.
c
	n_obs_g = n_local_g 
	n_obs_b = nn

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

!       Check for no obs
        if(nn .eq. 0)then
            write(6,*)' NOTE: no SND appended due to no tower obs'
            return
        endif

c
c.....  Call the routine to write the SND file.
c

        print *
	print *,'  Appending SND file, # of obs (in box) = ',n_obs_b
c
c
c.....	That's about it...let's go home.
c
	write(6,*)' Normal completion of TOWER_DRIVER'

        return
	end


 


