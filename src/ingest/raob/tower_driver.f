
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
c
        integer    wmoid(maxsta), istatus
        integer    dpchar(maxsta), narg, iargc
c
        character  stations(maxsta)*20
        character  reptype(maxsta)*6, atype(maxsta)*6
	character  atime*24, outfile*200
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_l*150
c
        character*200 path_to_metar
        character*200 path_to_local_data
        character*8   metar_format
        character*8   a9_to_a8, a8_time
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
	call get_directory('snd',outfile,len)
	outfile = outfile(1:len)//filename9(1:9)//'.snd'
cc	outfile = filename9(1:9)//'.snd'
c
c.....  Read in lat/lon/topo
c.....	Find east/west and north/south sides of grid (max extension of grid)
c
        call get_latlon_perimeter(ni,nj,0.0
     1                           ,lat,lon,topo
     1                           ,grid_north,grid_south
     1                           ,grid_east,grid_west,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading LAPS perimeter'
            return
        endif
c
c.....  Call the routine that reads the mesonet data files, then get the data.
c
        write(6,*)
	write(6,*)'Getting Mesonet Tower Data...'
c
        call get_local_towerobs(maxsta,i4time_sys,
     &                      path_to_local_data,metar_format,
     &                      itime_before,itime_after,
     &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,
     &                      nobs,stations,
     &                      reptype,atype,wmoid,
     &                      laps_cycle_time, istatus)

	if(istatus .ne. 1) then
	   print *, ' istatus=0 returned from get_local_towerobs...'
	   return
	endif

        if(nobs .gt. maxsta)then
           write(6,*)' ERROR: nobs > maxsta ',nobs,maxsta
           return
        endif
c
!       Final QC check 
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

!       Replace blank/UNK station names with wmoid if possible, else set blank
        iblank = 0
        do i = 1,nobs
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
        if(nobs .eq. 0)then
            write(6,*)' NOTE: no SND appended due to no tower obs'
            return
        endif

c
c
c.....	That's about it...let's go home.
c
	write(6,*)' Normal completion of TOWER_DRIVER'

        return
	end


 


