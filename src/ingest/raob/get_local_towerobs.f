c
        subroutine get_local_towerobs(maxobs,maxsta,i4time_sys,
     &                 path_to_local_data,local_format,
     &                 itime_before,itime_after,
     &                 eastg,westg,anorthg,southg,
     &                 lat,lon,ni,nj,
     &                 nobs,stations,
     &                 reptype,atype,wmoid,
     &                 laps_cycle_time, istatus)

c
c.....  Input variables/arrays
c
        integer maxobs ! raw data file
        integer maxsta ! processed stations for SND file
        integer maxlvls ! raw/processed stations for SND file

        parameter (maxlvls=100)

        character*(*) path_to_local_data, local_format
c
c.....  Local variables/arrays
c
        real    lat(ni,nj), lon(ni,nj)

!       Obs arrays (raw files)
	real*8  rh_time(maxobs), p_time(maxobs)
	real*8  t_time(maxobs), dd_time(maxobs), gust_time(maxobs)
	real*8  ff_time(maxobs), timeobs(maxobs)
	real*4  lats(maxobs), lons(maxobs), elev(maxobs)
	real*4  t(maxobs), td(maxobs), rh(maxobs), stnp(maxobs)
	real*4  dd(maxobs), ff(maxobs), ddg(maxobs), ffg(maxobs)
	real*4  mslp(maxobs), alt(maxobs), vis(maxobs)
	character  pro(maxobs)*11
        character*9 a9time_before, a9time_after, a9time_a(maxobs)

        real*4  lvls_m(maxlvls,maxobs)

!       Station arrays (snd file)
        real*4  stalat(maxsta),stalon(maxsta),staelev(maxsta)
        real*4  tsta_c(maxsta),tdsta_c(maxsta)
        real*4  ddsta(maxsta),ffsta(maxsta)

        integer*4  i4time_ob_a(maxobs), before, after
        integer*4  i4_timeobs(maxobs)
	character  provider(maxsta)*11
        logical l_dupe(maxobs)
c
c.....  Output arrays.
c
	integer*4  wmoid(maxsta)
	integer    rtime
	integer    recNum, nf_fid, nf_vid, nf_status
c
	character  stname(maxobs)*6, save_stn(maxobs)*6
	character  timech*9, time*4
	character  stations(maxsta)*20
	character  reptype(maxsta)*6, atype(maxsta)*6
	character  stn_type(maxobs)*11
        character*13 filename13, cvt_i4time_wfo_fname13
        character*150 data_file 
c
c.....  Start.
c
c
c.....	Set istatus flag for the local data to bad until we find otherwise.
c
	istatus = -1

        call get_ibadflag(ibadflag,istatus)
        if(istatus .ne. 1)return

        call get_sfc_badflag(badflag,istatus)
        if(istatus .ne. 1)return

!       call get_box_size(box_size,istatus)
!       if(istatus .ne. 1)return

        box_size = 1.0

c
c.....  Figure out the size of the "box" in gridpoints.  User defines
c.....  the 'box_size' variable in degrees, then we convert that to an
c.....  average number of gridpoints based on the grid spacing.
c

        call get_grid_spacing_cen(grid_spacing_m,istatus)
        box_length = box_size * 111.137 !km/deg lat (close enough for lon)
        ibox_points = box_length / (grid_spacing_m / 1000.) !in km
c
c.....	Zero out the counters.
c
        nobs = 0	        ! # of local obs in the laps grid
c
c.....  Get the data from the NetCDF file.  First, open the file.
c.....  If not there, return to obs_driver.
c
        ix = 1
c
c.....  Set up the time window.
c
	before = i4time_sys - itime_before
	after  = i4time_sys + itime_after

!       Ob times contained in each file
        i4_contains_early = 0 
        i4_contains_late = 3599

        call get_filetime_range(before,after                
     1                         ,i4_contains_early,i4_contains_late       
     1                         ,3600                                     
     1                         ,i4time_file_b,i4time_file_a)              

        do i4time_file = i4time_file_a, i4time_file_b, -3600

            call s_len(path_to_local_data,len_path)
            filename13= cvt_i4time_wfo_fname13(i4time_file)

            if(len_path .lt. 1)goto590
 	    data_file = path_to_local_data(1:len_path)//filename13

            write(6,*)' LDAD tower file = ',data_file(1:len_path+13)
c
c.....  Call the read routine.
c
	    call read_local_tower(data_file, maxobs, maxlvls,     ! I
     &         r_missing_data,                                    ! I
     &         nobs, nlvls, lvls_m(1,ix),                         ! O
     &         elev(ix), lats(ix), lons(ix),                      ! O
     &         t(ix), td(ix), dd(ix), ff(ix),                     ! O
     &         i4_timeobs(ix), stname(ix), wmoid(ix),             ! O
     &         istatus)                                           ! O
 
	    if(istatus .ne. 1)then
                write(6,*)
     1          '     Warning: bad status return from read_local_tower'       
                n_local_file = 0

            else
                n_local_file = nobs
                write(6,*)'     n_local_file = ',n_local_file

            endif

            ix = ix + n_local_file

590     enddo ! i4time_file

        n_local_all = ix - 1
        write(6,*)' n_local_all = ',n_local_all
c
        max_write = 100
      
c
c
c..................................
c.....	Second QC loop over all the obs.
c..................................
c
	jfirst = 1
        box_low = 1. - float(ibox_points)    !buffer on west/south side
        box_idir = float( ni + ibox_points)  !buffer on east
        box_jdir = float( nj + ibox_points)  !buffer on north
c
	do i=1,n_local_all
c
c.....  Bounds check: is station in the box?  Find the ob i,j location
c.....  on the LAPS grid, then check if outside past box boundary.
c
!          Test for badflag 
           call latlon_to_rlapsgrid(lats(i),lons(i),lat,lon,ni,nj,       
     &                              ri_loc,rj_loc,istatus)
           if(ri_loc.lt.box_low .or. ri_loc.gt.box_idir
     1   .or. rj_loc.lt.box_low .or. rj_loc.gt.box_jdir) then
               if(i .le. max_write)then
                   write(6,81,err=125)i,wmoid(i),stname(i)
     1                               ,nint(ri_loc),nint(rj_loc)
 81                format(i6,i7,1x,a8,' out of box ',2i12)
               endif
               go to 125
           endif
c
c.....  Elevation ok?
c
	   if(elev(i).gt.5200. .or. elev(i).lt.-400.) go to 125       
c
c.....  Check to see if its in the desired time window.
c
	   if(i4time_ob_a(i) .lt. before 
     1   .or. i4time_ob_a(i) .gt. after) then
               if(i .le. max_write)then
                   write(6,91,err=125)i,wmoid(i),stname(i)
     1                               ,a9time_a(i),before
     1                               ,after
 91		   format(i6,i7,1x,a8,' out of time ',a11,2i12)
               endif
               go to 125
           endif
c
c.....  Right time, right location...

           timech = a9time_a(i)
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
             if(stname(i) .eq. save_stn(k)) then
                 write(6,*)' Rejecting duplicate ',i,k,stname(i)
     1                    ,' ',a9time_a(i),' ',a9time_a(k)
                 go to 125
             endif
	   enddo !k
c
	   icount = icount + 1
	   save_stn(icount) = stname(i)  ! only one...save for checking
c
 150	   nobs = nobs + 1

           if(nobs .gt. maxsta)then
              write(6,*)' ERROR in get_local_obs: increase maxsta '
     1                 ,nobs,maxsta
              stop
           endif
 
c
c.....  Check if its in the LAPS grid.
c
           if(ri_loc.lt.1. .or. ri_loc.gt.float(ni)) go to 151 !off grid
           if(rj_loc.lt.1. .or. rj_loc.gt.float(nj)) go to 151 !off grid
           nobs = nobs + 1  !on grid...count it
 151	   continue
c
c.....	Figure out the cloud data.
c.....     NOTE: Not currently reading cloud data from mesonets.
c
c
c
 125     continue
       enddo !i
c
c
c.....  That's it...lets go home.
c
c
c.....  Call the routine to write the SND file.
c

        print *
	print *,'  Appending SND file, # of obs (in grid) = ',nobs

        call write_snd(    lun_out                         ! I
     1                    ,maxsnd,maxlvl,nsnd              ! I
     1                    ,iwmostanum                      ! I
     1                    ,stalat,stalon,staelev           ! I
     1                    ,c5_staid,a9time_ob,c8_obstype   ! I
     1                    ,nlvl                            ! I
     1                    ,height_m                        ! I
     1                    ,pressure_pa                     ! I
     1                    ,temp_c                          ! I
     1                    ,dewpoint_c                      ! I
     1                    ,dir_deg                         ! I
     1                    ,spd_mps                         ! I
     1                    ,istatus)                        ! O


         istatus = 1            ! everything's ok...
         return
c
 990     continue               ! no data available
         istatus = 0
         print *,' No data available from GET_LOCAL_TOWEROBS'
         return
c
         end

         subroutine read_local_tower(data_file, maxobs, maxlvls,  ! I
     &         r_missing_data,                                    ! I
     &         nobs, nlvls, lvls_m,                               ! O
     &         elev, lats, lons,                                  ! O
     &         t, td, dd, ff,                                     ! O
     &         i4_timeobs, stname, wmoid,                         ! O
     &         istatus)                                           ! O

         return
         end
