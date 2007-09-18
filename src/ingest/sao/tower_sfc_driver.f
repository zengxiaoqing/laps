
        subroutine tower_sfc_driver(maxsta,i4time_sys                     ! I
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
        character  stations(maxsta)*20, provider(maxsta)*11
        character  reptype(maxsta)*6, atype(maxsta)*6
	character  atime*24
	character  dir_s*256, ext_s*31, units*10, comment*125,var_s*3
	character  filename9*9, filename13*13, a9time_metar_file*9
        character  fname9_to_wfo_fname13*13
	character  data_file_l*150
c
        character*200 path_to_metar
        character*200 path_to_tower_data
        character*8   metar_format
        character*8   a9_to_a8, a8_time

c       Declared then used in 'get_local_towerobs' for SND purposes
        real     stalat_s(maxsta,maxlvls),stalon_s(maxsta,maxlvls)
        real     staelev_s(maxsta)
        real     soilmoist_p(maxsta)       
	character  stname_s(maxsta)*5
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
c.....	Start here.  
c
        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

c.....  Get the time from the scheduler or from the user if interactive.
c
        call make_fnam_lp(i4time_sys,filename9,istatus)
c
        write(6,*)' systime = ',filename9
c
c.....  Call the routine that reads the mesonet data files, then get the data.
c
        write(6,*)
	write(6,*)'Getting Mesonet Tower Data...'
c
        ext_s = 'lso'

        call get_local_towerobs(maxsta,maxlvls,                          ! I
     &                      i4time_sys,lun_out,
     &                      path_to_tower_data,metar_format,
     &                      ext_s,                                       ! I
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


 


