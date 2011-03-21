
        subroutine tower_driver_sub(ni,nj,lun_out
     1                           ,maxobs,laps_cycle_time
     1                           ,path_to_local_data
     1                           ,itime_before,itime_after
     1                           ,maxsta                          ! I
     1                           ,istatus)
c        
        integer ni, nj, maxsta, maxobs 
        integer maxlvls ! raw/processed stations for SND file
        parameter (maxlvls=10)
c
	real    lat(ni,nj), lon(ni,nj), topo(ni,nj)
c
        integer    wmoid(maxsta), istatus
        integer    dpchar(maxsta), narg, iargc
c
!       character  stations(maxsta)*20
        character  reptype(maxsta)*6, atype(maxsta)*6
	character  atime*24, outfile*200
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_l*150
c
        character*200 path_to_metar
        character*200 path_to_local_data
        character*8   tower_format
        character*8   a9_to_a8, a8_time

c       Declared then used in 'get_local_towerobs' for SND purposes
        real     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real     staelev_s(maxsta)
        real     soilmoist_p(maxsta)       
	character  stname_s(maxsta)*5
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
        ext_s = 'snd'
        tower_format = 'NIMBUS'
 
        call get_local_towerobs(maxsta,maxlvls,                          ! I
     &                      i4time_sys,lun_out,
     &                      path_to_local_data,tower_format,
     &                      ext_s,
     &                      itime_before,itime_after,
!    &                      grid_east,grid_west,grid_north,grid_south,
     &                      lat,lon,ni,nj,                               ! I
     &                      nobs,                                        ! O
!    &                      stations,
!    &                      reptype,atype,wmoid,
!    &                      laps_cycle_time, 
     &                      stalat_s,stalon_s,staelev_s,                 ! O
     &                      stname_s,                                    ! O
     &                      soilmoist_p,                                 ! O
     &                      istatus)

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

!       Check for no obs
        if(nobs .eq. 0)then
            write(6,*)' NOTE: no SND appended due to no tower obs'
            return
        endif
c
c
c.....	That's about it...let's go home.
c
	write(6,*)' Normal completion of TOWER_DRIVER, nobs = ',nobs

        return
	end


 


