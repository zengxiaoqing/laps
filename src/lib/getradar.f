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


!       File: getradar.f
!
cdoc    Module Summary:
cdoc
cdoc    get_multiradar_vel
cdoc        now called from wind/lplot
cdoc
cdoc    get_radar_ref
cdoc        now called from lplot
cdoc
cdoc    read_radar_3dref
cdoc        now called from deriv
cdoc        now called from accum
cdoc        now called from lplot
cdoc        now called from get_radar_ref
cdoc        now called from mosaic_radar
cdoc
cdoc    read_multiradar_3dref
cdoc        now called from cloud
cdoc        now called from read_radar_3dref
cdoc
cdoc    read_radar_2dref
cdoc        now called directly from wind-derived
cdoc        now called from (wind-derived via) get_radar_max_pd
cdoc
cdoc    read_radar_vel
cdoc        now called from (wind/lplot via) get_multiradar_vel
cdoc
cdoc    read_nowrad_3dref
cdoc        now called from read_multiradar_3dref
cdoc
cdoc
!       1996 Aug    S. Albers FSL


        subroutine get_multiradar_vel(
     1   i4time_ref,i4time_tol,i4time_radar_a
     1  ,max_radars,n_radars,ext_a,r_missing_data
     1  ,l_apply_map,imax,jmax,kmax
     1  ,grid_ra_vel,grid_ra_nyq,v_nyquist_in_a,n_vel_a
     1  ,rlat_radar_a,rlon_radar_a,rheight_radar_a,radar_name_a
     1  ,istatus_multi_vel,istatus_multi_nyq)


cdoc    Returns Velocity from multiple radars.
cdoc    Called from wind/lapsplot


!       i4time_ref          Input   Desired i4time
!       i4time_tol          Input   Half Width of allowable time window
!       i4time_radar_a      Output  Actual times of data you are getting
!       max_radars          Input   Dimensioning for maximum # of radars
!       n_radars            Output  Actual number of radars returned with at
!                                   least one valid velocity measurement
!       ext_a               Local   Array: Possible extensions
!       r_missing_data      Input
!       mode                Input   (1) normal radar data, (2) return clutter map
!       l_apply_map         Input   .true. - remove 3D ground clutter
!       imax,jmax,kmax      Input   LAPS 3D Grid Dimensions
!       lat,lon             Input   2D latitude and longitude arrays (degrees)
!       topo                Input   2D terrain array (meters)
!       grid_ra_vel         Output  4D Velocity Grid
!       grid_ra_nyq         Output  4D Nyquist Velocity Grid
!       v_nyquist_in_a      Output  Array: volume nyquist velocity of the radars
!       n_vel_a             Output  Array: # of grid points with measurable velocity
!       rlat_radar_a        Output  Array: Radar Latitude (Degrees)
!       rlon_radar_a        Output  Array: Radar Longitude (Degrees)
!       rheight_radar_a     Output  Array: Radar Height (m MSL)
!       radar_name_a        Output  Array: Radar Name (Character*4)
!       istatus_multi_vel   Output  Data is useable for 3D vel applications
!       istatus_multi_nyq   Output  Data is useable for 3D nyq applications

!
!       The domain must be specified using grid_fnam_common

        character*31 ext_a(max_radars)

        character*255 c_filespec
        character*150  directory

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        character asc9_tim_radar*9

        logical l_apply_map

        real*4 grid_ra_vel(imax,jmax,kmax,max_radars)
        real*4 grid_ra_nyq(imax,jmax,kmax,max_radars)

        real*4 rlat_radar_a(max_radars),rlon_radar_a(max_radars)
     1    ,rheight_radar_a(max_radars),v_nyquist_in_a(max_radars)

        character*4 radar_name_a(max_radars)

        integer*4 n_vel_a(max_radars)
        integer*4 i4time_radar_a(max_radars)

        logical l_conus_vel

        l_conus_vel = .false.

!       Initialize the flags
        n_radars = 0
        n_ref = 0
        istatus_multi_vel = 1
        istatus_multi_nyq = 1

!       Initialize this stuff
        do i = 1,max_radars
            if(i .lt. 10)then
                write(ext_a(i),91)i
 91             format('v0',i1)
            else
                write(ext_a(i),92)i
 92             format('v',i2)
            endif
            n_vel_a(i) = 0
        enddo

        if(l_conus_vel)then
!           Read CONUS velocity data
!           call read_raw_conus_vel()
        endif

!       Loop through the potential radars
        do i_radar_pot = 1,max_radars

            write(6,*)
            write(6,*)' Looking for potential radar # ',i_radar_pot

            call get_directory(ext_a(i_radar_pot),directory,len_dir)
            c_filespec = directory(1:len_dir)//'*.'//ext_a(i_radar_pot)

            call get_file_time(c_filespec,i4time_ref,i4time_radar)

            call make_fnam_lp(i4time_radar,asc9_tim_radar,istatus)

            if(abs(i4time_radar - i4time_ref) .gt. i4time_tol)then
                write(6,101)i_radar_pot,i4time_tol,asc9_tim_radar
101             format(' Radar',i3,' not available within',i6,
     1          ' sec of analysis time, nearest = ',a9)

            else
                n_radars = n_radars + 1

                do k = 1,kmax
                do j = 1,jmax
                do i = 1,imax
                    grid_ra_vel(i,j,k,n_radars) = r_missing_data
                    grid_ra_nyq(i,j,k,n_radars) = r_missing_data
                enddo
                enddo
                enddo

                i4time_radar_a(n_radars) = i4time_radar

                write(6,*)' Reading vel/nyq for actual radar # '
     1                    ,n_radars

                call read_radar_vel(i4time_radar,l_apply_map,
     1           imax,jmax,kmax,ext_a(i_radar_pot),
     1           grid_ra_vel(1,1,1,n_radars),
     1           grid_ra_nyq(1,1,1,n_radars),v_nyquist_in_a(n_radars),       
     1           rlat_radar_a(n_radars),rlon_radar_a(n_radars)
     1                      ,rheight_radar_a(n_radars)
     1                      ,radar_name_a(n_radars)
     1                      ,n_vel_a(n_radars),
     1                       istatus_vel,istatus_nyq)

                if(n_vel_a(n_radars) .eq. 0 .or. istatus_vel .ne. 1)then       
                    write(6,*)' No valid velocities for radar ',n_radars
                  ! Don't count in a valid radar
                    n_radars = n_radars - 1
                endif

                if(istatus_vel .ne. 1)then
                    istatus_multi_vel = 0
                endif

                if(istatus_nyq .ne. 1)then
                    istatus_multi_nyq = 0
                endif

            endif ! We scored a radar

        enddo ! i_radar_pot

        if(n_radars .gt. max_radars)then
            write(6,*)
     1           ' ERROR in get_multiradar_vel: n_radars > max_radars'
     1           ,n_radars,max_radars
            n_radars = 0
            istatus_multi_vel = 0
            istatus_multi_nyq = 0
        endif

        if(n_radars .eq. 0)then
            write(6,*)' WARNING in get_multiradar_vel: n_radars = 0'
            istatus_multi_vel = 0
            istatus_multi_nyq = 0
        endif

        if(l_conus_vel)then ! Try to supersede 3-D ref data with CONUS
!           Read CONUS reflectivity data
!           call read_raw_conus_vel()
!           istatus_multi_vel = 0
        endif

        return
        end


        subroutine get_radar_ref(i4time_ref,i4time_tol,i4time_radar
     1        ,mode,l_apply_map,imax,jmax,kmax
     1        ,lat,lon,topo,l_low_fill,l_high_fill
     1        ,heights_3d
     1        ,grid_ra_ref,n_ref
     1        ,rlat_radar,rlon_radar,rheight_radar
     1        ,istatus_2dref,istatus_3dref)

cdoc    This routine returns 3D radar ref data from the LAPS radar files
cdoc    Called from lapsplot

!       i4time_ref          Input   Desired i4time
!       i4time_tol          Input   Half Width of allowable time window
!       i4time_radar        Output  Actual time of data you are getting
!       mode                Input   (1) normal radar data, (2) return clutter map
!       l_apply_map         Input   .true. - remove 3D ground clutter
!       imax,jmax,kmax      Input   LAPS 3D Grid Dimensions
!       lat,lon             Input   2D latitude and longitude arrays (degrees)
!       topo                Input   2D terrain array (meters)
!       l_low_fill          Input   Flag to fill vertical gaps in
!                                   low level reflectivities (normally .true.)
!       l_high_fill         Input   Flag to fill vertical gaps in
!                                   high level reflectivities (normally .true.)
!       lat,lon,topo are needed if l_low_fill or l_high_fill are true
!       grid_ra_ref         Output  3D Reflectivity Grid
!       n_ref               Output  # of grid points with measurable reflectivity
!       rlat_radar          Output  Radar Latitude (Degrees)
!       rlon_radar          Output  Radar Longitude (Degrees)
!       rheight_radar       Output  Radar Height (m MSL)
!       istatus_2dref       Output  Data is useable for 2D ref applications
!       istatus_3dref       Output  Data is useable for 3D ref applications

!       Steve Albers           Dec 7 1995

!       USER NOTE: Only one source of radar data per domain is currently
!       permitted with get_radar. The constraint on this is setting up
!       the inputted time window to work with multiple sources of radar data.
!       This may be changed in the future.
!
!       The domain must be specified using grid_fnam_common

        character*255 c_filespec
        character*150  directory

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        character asc9_tim_radar*9

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)
        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

        logical l_low_fill,l_high_fill,l_apply_map,l_clutter_found

        logical l_conus_ref, l_nowrad_3dref

        character*31 radarext
        character*3  ext_gc
        character*4 radar_name

        l_conus_ref = .false.
        l_nowrad_3dref = .false.

        radarext = 'vrc'
        call get_directory(radarext,directory,len_dir)

!       Build the filespec names
        c_filespec = directory(1:len_dir)//'*.'//radarext(1:3)

        IF(mode .eq. 2)then ! Get clutter Map only

        l_clutter_found = .false.

        open(16,
     1    file='clutter.'
     1                                  //grid_fnam_common
     1                  ,status='old',err=1992)


!       Check Dimensions
        read(16,555,end=1990,err=1992)i_dim,j_dim,k_dim ! Read data
        if(i_dim .ne. imax .or. j_dim .ne. jmax
     1                   .or. k_dim .ne. kmax)then
            write(6,*)' Clutter map has wrong dimensions'
            istatus_2dref = 0
            istatus_3dref = 0
            return
        endif

1400    read(16,1401,end=1990,err=1992)ext_gc,i4time_1,i4time_2 ! Read the time
1401    format(1x,a3,2(1x,i10))

        if(.not. l_clutter_found)n_gc_pot = 0

1500    read(16,555,end=1990,err=1992)i_grid,j_grid,k_grid,l_data ! Read data
555     format(1x,3i3,i4)

        if(l_data .eq. 9999)goto1400 ! Update the time if possible (EOData)

        if(i4time_1 .le. i4time_radar .and. i4time_radar .le. i4time_2
     1                    .and. ext_gc .eq. radarext(1:3)
     1                                                          )then

            n_gc_pot = n_gc_pot + 1
            grid_ra_ref(i_grid,j_grid,k_grid) = l_data - 900
            l_clutter_found = .true.

        endif

        goto1500

1990    close(16)

1992    write(6,*)' # Grid boxes = ',n_gc_pot

        if(l_clutter_found)then
            istatus_2dref = 1
            istatus_3dref = 1
        else
            write(6,*)' No Ground Clutter Map Found'
            istatus_2dref = 0
            istatus_3dref = 0
        endif

        return

        ENDIF ! MODE

        call get_file_time(c_filespec,i4time_ref,i4time_radar)

        call make_fnam_lp(i4time_radar,asc9_tim_radar,istatus)

        if(abs(i4time_radar - i4time_ref) .gt. i4time_tol)then
            write(6,101)i4time_tol,asc9_tim_radar
101         format('  No radar data available within',i6,
     1          ' sec of analysis time, nearest = ',a9)
            n_ref = 0
            n_vel = 0
            istatus_2dref = 0
            istatus_3dref = 0

        else ! Radar data is in time window

            if(l_conus_ref .or. l_nowrad_3dref)then

                if(l_conus_ref)then
!                   Read CONUS reflectivity data
!                   call read_raw_conus_ref()
!                   istatus_3dref = 1
                elseif(l_nowrad_3dref)then
!                   Read NOWRAD 3D reflectivity data
!                   call read_nowrad_3dref()
!                   istatus_3dref = 1
                endif


            else ! Attempt to get reflectivity from V01/VRC

                call get_ref_base(ref_base,istatus)
                if(istatus .ne. 1)then
                    istatus_2dref = 0
                    istatus_3dref = 0
                    return
                endif

                call read_radar_3dref(i4time_radar,l_apply_map,       ! I
     1               ref_base,                                        ! I
     1               imax,jmax,kmax,radarext,                         ! I
     1               lat,lon,topo,l_low_fill,l_high_fill,             ! I
     1               heights_3d,                                      ! I
     1               grid_ra_ref,                                     ! O
     1               rlat_radar,rlon_radar,rheight_radar,radar_name,  ! O
     1               n_ref,istatus_2dref,istatus_3dref)               ! O

            endif

        endif ! Radar data is not in time window

        return
        end


        subroutine read_radar_3dref(i4time_radar,                        ! I
!    1   i4_tol,i4_ret,
     1   l_apply_map,ref_missing,                                        ! I
     1   imax,jmax,kmax,radarext,                                        ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                            ! I
     1   heights_3d,                                                     ! I
     1   grid_ra_ref,                                                    ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                 ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                        ! O

cdoc    Steve Albers Feb 1998   This routine will read in a 3D radar
cdoc                            reflectivity field. It is a jacket that
cdoc                            calls read_multiradar_3dref.


        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

!       If a grid "overall" has the following...
!       2dref=1, 3dref=1 - echo top from radar has better confidence than
!                          the associated cloud top from satellite
!       2dref=1, 3dref=0 - echo top from radar has less confidence
!       2dref=0, 3dref=0 - missing data

        integer*4 istatus_2dref
        integer*4 istatus_3dref
        integer*4 istatus_2dref_a(imax,jmax)
        integer*4 istatus_3dref_a(imax,jmax)

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        character*31 radarext

        character*4 radar_name

        logical l_low_fill,l_high_fill,l_apply_map

        write(6,*)' Subroutine read_radar_3dref'

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            n_2dref = 0
            n_3dref = 0
            goto900
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            n_2dref = 0
            n_3dref = 0
            goto900
        endif

        call read_multiradar_3dref(i4time_radar,
     1   l_apply_map,ref_missing,
     1   imax,jmax,kmax,radarext,
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   heights_3d,
     1   grid_ra_ref,
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,
     1   n_ref_grids,n_2dref,n_3dref,istatus_2dref_a,istatus_3dref_a)       

 900    if(    n_2dref .eq. imax*jmax)then    ! Full coverage
            istatus_2dref = +1                
        elseif(n_2dref .gt. 0        )then    ! Partly missing
            istatus_2dref = -1                
        else                                  ! All missing
            istatus_2dref = 0                 
        endif        

        if(    n_3dref .eq. imax*jmax)then    ! Full coverage
            istatus_3dref = +1                
        elseif(n_3dref .gt. 0        )then    ! Partly missing
            istatus_3dref = -1                
        else                                  ! All missing
            istatus_3dref = 0                 
        endif        

        write(6,*)'n_2dref,n_3dref=',n_2dref,n_3dref

        return
        end


        subroutine read_multiradar_3dref(i4time_radar,                  ! I
     1   l_apply_map,ref_missing,                                       ! I
     1   imax,jmax,kmax,radarext,                                       ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                           ! I
     1   heights_3d,                                                    ! I
     1   grid_ra_ref,                                                   ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                ! O
     1   n_ref_grids,n_2dref,n_3dref,istatus_2dref_a,istatus_3dref_a)   ! O 

!       Steve Albers Nov 1998   This routine will read in a 3D radar
!                               reflectivity field. It will either read
!                               in 3D radar refs and do the vertical filling
!                               procedure, or read in the 2D data and do
!                               a quick and dirty vertical fill. Information
!                               is passed back about which data is 2d and which
!                               is 3d.
!
!       SA           Mar 2001   The 'all' option is now being tested that 
!                               merges 3D (v01) and 2D (vrc) data.

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)
        real*4 radar_2dref(imax,jmax)
        real*4 closest_vxx(imax,jmax)

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

!       If a grid element has the following...
!       2dref=1, 3dref=1 - echo top from radar has better confidence than
!                          the associated cloud top from satellite
!       2dref=1, 3dref=0 - echo top from radar has less confidence
!       2dref=0, 3dref=0 - missing data

        integer*4 istatus_2dref
        integer*4 istatus_3dref
        integer*4 istatus_2dref_a(imax,jmax)
        integer*4 istatus_3dref_a(imax,jmax)

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        character*31 radarext

        character*4 radar_name

        logical l_low_fill,l_high_fill,l_apply_map, l_parse

        write(6,*)' Subroutine read_multiradar_3dref'

        i4_tol = 1200

        istatus_2dref_a = 0
        istatus_3dref_a = 0

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            goto900
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            goto900
        endif

        closest_vxx = r_missing_data

!       Initialize 3d reflectivity array with default value
        grid_ra_ref = r_missing_data ! ref_base 

        if(radarext(1:2) .eq. 'v0' .or.
     1     radarext(1:2) .eq. 'v1' .or.
     1     radarext(1:2) .eq. 'v2' .or.
     1     radarext(1:3) .eq. 'all' )then     ! Read Doppler radar ref data
                                              ! from NetCDF files

            write(6,*)' Reading Reflectivity Data from 3D file '
     1                                                 ,radarext

            if(radarext(1:3) .ne. 'all')then
                ext = radarext
            else
                ext = 'v01'
            endif

!           Read Reflectivity
            var_2d = 'REF'
            call get_laps_3dgrid(i4time_radar,i4_tol,i4_ret
     1                    ,imax,jmax,kmax,ext,var_2d
     1                    ,units_2d,comment_2d,grid_ra_ref,istatus)

            if(istatus .eq. 1)then
                read(comment_2d,558)rlat_radar,rlon_radar,rheight_radar
     1                             ,n_ref_grids
558             format(2f9.3,f8.0,i7)

                if(l_low_fill .or. l_high_fill)then
                    call ref_fill_vert(grid_ra_ref,imax,jmax,kmax
     1                          ,l_low_fill,l_high_fill,lat,lon,topo
     1                          ,heights_3d
     1                          ,rlat_radar,rlon_radar,rheight_radar
     1                          ,istatus_rfill)

                    if(istatus_rfill .eq. 1)then
                        write(6,*)' Reflectivity data filled in'
                    else
                        write(6,*)' Reflectivity data fill error'
                    endif

                else
                    write(6,*)' Reflectivity not filled in'

                endif
 
                if(.false.)then
                    call ref_fill_horz(grid_ra_ref,imax,jmax,kmax
     1                ,lat,lon,rlat_radar,rlon_radar,rheight_radar
     1                ,istatus)
                    if(istatus .ne. 1)then
                        istatus_2dref_a = 0
                        istatus_3dref_a = 0
                    endif
                endif

                if(l_low_fill .or. l_high_fill)then
                    do j = 1,jmax
                    do i = 1,imax
                        istat_g = 0
                        do k = 1,kmax
                            if(grid_ra_ref(i,j,k) .ne. r_missing_data
     1                                                             )then       
                                istat_g = 1
                            endif
                        enddo ! k
                      
                        if(istat_g .eq. 1)then ! valid radar data at this point
                            istatus_2dref_a(i,j) = 1
                            istatus_3dref_a(i,j) = 1
                            call latlon_to_radar(
     1                          lat(i,j),lon(i,j),topo(i,j)
     1                         ,azimuth,closest_vxx(i,j),elev
     1                         ,rlat_radar,rlon_radar,rheight_radar)       

                        endif

                    enddo ! i
                    enddo ! j

                endif

            else
                write(6,*)' Radar reflectivity data cannot be read in: '       
     1                    ,ext

            endif ! Success as reflectivity

        endif ! 'vxx'

        if(radarext(1:3) .eq. 'vrc' .or. radarext(1:3) .eq. 'all')then       
50          write(6,*)' Reading NOWRAD/vrc data' ! lumped together for now?

            var_2d = 'REF'
            ext = 'vrc'
            call get_laps_2dgrid(i4time_radar,i4_tol,i4_ret,ext,var_2d
     1                          ,units_2d,comment_2d,imax,jmax
     1                          ,radar_2dref,0,istatus_vrc)

            write(6,*)' istatus_vrc = ',istatus_vrc

            if(istatus_vrc .eq. 1 .or. istatus_vrc .eq. -1)then       
                if(l_parse(comment_2d,'WSI'))then
                    radar_name = 'WSI '
                else
                    len_comment = 37
                    radar_name = comment_2d(len_comment-3:len_comment)       
                endif 
                    
                write(6,*)' Read radar ',radar_name             

                do i = 1,imax
                do j = 1,jmax
                    if(radar_2dref(i,j) .ne. r_missing_data)then
                        istatus_2dref_a(i,j) = 1
                    endif
                enddo ! j
                enddo ! i                 

!               Conditionally fill up the 3D reflectivity array 
!                          (useful for precip type, get_low_ref)

                if(l_low_fill .or. l_high_fill)then
                    do i = 1,imax
                    do j = 1,jmax

!                       Set distant 3D radar points to msg to select NOWRAD/vrc
                        if(closest_vxx(i,j) .ne. r_missing_data
     1               .and. closest_vxx(i,j) .gt. 260000.     )then       
                            istatus_3dref_a(i,j) = 0
                        endif

                        if(istatus_3dref_a(i,j) .eq. 0)then 
                            do k = 1,kmax ! Fill column with NOWRAD/vrc
!                               if(grid_ra_ref(i,j,k).eq.r_missing_data      
!    1                                                             )then       
                                    grid_ra_ref(i,j,k)=radar_2dref(i,j)      
!                               endif
                            enddo ! k
                        endif ! istatus_3dref_a = 0
                    enddo ! j
                    enddo ! i

                endif ! l_low_fill OR l_high_fill

            else
                istatus_2dref_a = 0
                istatus_3dref_a = 0

            endif

        endif ! vrc

        if(radarext(1:3) .eq. 'vrz')then
!           call read_raw_conus_ref()
            write(6,*)' Radar vrz reflectivity data cannot be read in'
            istatus_2dref_a = 0
            istatus_3dref_a = 0

        endif ! 'vrz'

        if(radarext(1:3) .eq. 'ln3')then
            call read_nowrad_3dref(i4time_radar
     1                            ,imax,jmax,kmax
     1                            ,grid_ra_ref,heights_3d
     1                            ,istatus_2dref_a,istatus_3dref_a)

        endif ! 'ln3'

      ! Summary stats
 900    n_2dref = 0
        n_3dref = 0

        do i=1,imax
        do j=1,jmax
            n_2dref = n_2dref + istatus_2dref_a(i,j)
            n_3dref = n_3dref + istatus_3dref_a(i,j)
          
            if(istatus_2dref_a(i,j) .eq. 0)then
                do k = 1,kmax
                    grid_ra_ref(i,j,k) = ref_missing
                enddo ! k
            endif 

        enddo ! j
        enddo ! i

        n_missing = imax*jmax - n_2dref

        write(6,*)'n_2dref,n_3dref=',n_2dref,n_3dref
        write(6,*)'n_missing/ref_missing=',n_missing,ref_missing

        return
        end



        subroutine read_radar_2dref(i4time_radar,radar_name,
     1                  imax,jmax,
     1                  ref_2d,istatus_2dref)

!       Steve Albers Feb 1996   Reads Radar data into 2d REF arrays
!           now called from (wind-derived via) get_radar_max_pd


        real*4 ref_2d(imax,jmax)

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        character*4 radar_name

        write(6,*)' Subroutine read_radar_2dref'

50      write(6,*)' Reading VRC/NOWRAD data'

        radar_name = 'WSI '

        var_2d = 'REF'
        ext = 'vrc'
        call get_laps_2d(i4time_radar,ext,var_2d
     1          ,units_2d,comment_2d,imax,jmax,ref_2d,istatus_vrc)

        istatus_2dref = istatus_vrc

        return
        end


        subroutine read_radar_vel(i4time_radar,l_apply_map,
     1   imax,jmax,kmax,radarext,
     1   grid_ra_vel,grid_ra_nyq,v_nyquist_in,
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,
     1   n_vel_grids,istatus_vel,istatus_nyq)

!       Steve Albers Dec 7 1995 Reads Radar data into VEL array
!                               This routine will be able to select between
!                               various sources of radar data over various
!                               domains

!       USER NOTE: Only one source of radar data per domain is currently
!       permitted with read_radar. The constraint on this is setting up
!       the inputted time to work with multiple sources of radar data.
!       This may be changed in the future. The likely place for the selection
!       of multiple radars is in the calling routine 'get_radar_data'.
!       The domain must be specified using grid_fnam_common


        real*4 grid_ra_vel(imax,jmax,kmax)
        real*4 grid_ra_nyq(imax,jmax,kmax)

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        character*31 radarext

        character*4 radar_name

        logical l_apply_map,l_unfold

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        write(6,*)' Subroutine read_radar_vel'

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            istatus_vel = 0
            istatus_nyq = 0
            return
        endif

        if(radarext(1:3) .eq. 'vrc')then
            write(6,*)' Error NOWRAD has no velocity data'
            istatus_vel = 0

        else ! Read Doppler radar data from NetCDF files
            write(6,*)' Reading Doppler Data from 3D file ',radarext

            ext = radarext

!           Read Velocity
            var_2d = 'VEL'

            write(6,*)

            call get_laps_3d(i4time_radar,imax,jmax,kmax,ext,var_2d
     1                    ,units_2d,comment_2d,grid_ra_vel,istatus_vel)
            write(6,*)' Status of velocity = ',istatus_vel

            read(comment_2d,559,err=560)rlat_radar,rlon_radar
     1            ,rheight_radar
     1            ,n_vel_grids,radar_name,v_nyquist_in,l_unfold
559         format(2f9.3,f8.0,i7,a4,3x,e12.4,l2)
560         print *, ' comment_string: ',comment_2d(1:60)
            print *, ' v_nyquist_in = ',v_nyquist_in


            write(6,*)' Header info for radar: ',radar_name
            write(6,*)rlat_radar,rlon_radar,rheight_radar
     1                  ,n_vel_grids,v_nyquist_in

!           Read Nyquist Velocity
            var_2d = 'NYQ'

            write(6,*)

            call get_laps_3d(i4time_radar,imax,jmax,kmax,ext,var_2d
     1                    ,units_2d,comment_2d,grid_ra_nyq,istatus_nyq)
            write(6,*)' Status of Nyquist velocity = ',istatus_nyq

        endif ! Test of extension (hence radar type )

        return
        end


        subroutine ground_clutter(
     1              imax,jmax,kmax
     1             ,istatus_2dref,istatus_3dref,istatus_vel
     1             ,radarext
     1             ,grid_ra_ref
     1             ,ref_base
     1             ,l_apply_map
     1             ,i4time_radar
     1                                                  )

        logical l_apply_map,l_clutter_found

        character*80 grid_fnam_common
        character*3  ext_gc
        character*31 radarext,ext
        real*4 grid_ra_ref(imax,jmax,kmax)
        character*150 directory

        common / grid_fnam_cmn / grid_fnam_common
        common /laps_diag/ no_laps_diag

        n_gc_grids = 0
        n_gc_pot = 0
        l_clutter_found = .false.

        if(          (.not. l_apply_map)
     1                                                  )goto1992

!       Pull out static directory
        ext = grid_fnam_common
        call get_directory(ext,directory,len_dir)

        open(16,file=directory(1:len_dir)//'clutter.'//grid_fnam_common
     1  ,status='old',err=1991)


!       Check Dimensions
        read(16,555,end=1990,err=1992)i_dim,j_dim,k_dim ! Read data
        if(i_dim .ne. imax .or. j_dim .ne. jmax
     1                   .or. k_dim .ne. kmax)then
            write(6,*)' Clutter map has wrong dimensions'
            istatus_2dref = 0
            istatus_3dref = 0
            istatus_vel = 0
            return
        endif

1400    read(16,1401,end=1990,err=1992)ext_gc,i4time_1,i4time_2 ! Read the time
1401    format(1x,a3,2(1x,i10))

        if(.not. l_clutter_found)n_gc_pot = 0

1500    read(16,555,end=1990,err=1992)i_grid,j_grid,k_grid,l_data ! Read data
555     format(1x,3i3,i4)

        if(l_data .eq. 9999)goto1400 ! Update the time if possible (EOData)

        if(i4time_1 .le. i4time_radar .and. i4time_radar .le. i4time_2
     1                    .and. ext_gc .eq. radarext(1:3)
     1                                                          )then

            n_gc_pot = n_gc_pot + 1

            if(grid_ra_ref(i_grid,j_grid,k_grid) -
     1       float(l_data ) - 900. .lt. 5    .and.
     1        grid_ra_ref(i_grid,j_grid,k_grid) .gt. ref_base)then
                grid_ra_ref(i_grid,j_grid,k_grid) = ref_base
                n_gc_grids = n_gc_grids + 1

            endif

            l_clutter_found = .true.

        endif

        goto1500

1990    close(16)

1991    if(l_clutter_found)then
!           istatus = 1
        else
            write(6,*)' No Ground Clutter Map Found'
            istatus_2dref = 0
            istatus_3dref = 0
            istatus_vel = 0
            return
        endif

1992    if(no_laps_diag .eq. 0)
     1  write(6,*)' # Grid boxes (gc/pot)  = ',n_gc_grids,n_gc_pot

        return
        end



        subroutine read_nowrad_3dref(i4time_radar
     1                              ,imax,jmax,kmax
     1                              ,ref_3d,heights_3d
     1                              ,istatus_2dref_a,istatus_3dref_a)       

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        real*4 ref_3d(imax,jmax,kmax)
        real*4 heights_3d(imax,jmax,kmax)

        real*4 ref_2d_04(imax,jmax)
        real*4 ref_2d_48(imax,jmax)
        real*4 ref_2d_8c(imax,jmax)
        real*4 ref_2d_et(imax,jmax)

        integer*4 istatus_2dref_a(imax,jmax)
        integer*4 istatus_3dref_a(imax,jmax)

        write(6,*)' Subroutine read_nowrad_3dref'

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading ref_base parameter'
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading r_missing_data parameter'
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

        icount_ref1 = 0
        icount_ref2 = 0

        var_2d = 'R04'
        ext = 'ln3'
        call get_laps_2d(i4time_radar,ext,var_2d
     1       ,units_2d,comment_2d,imax,jmax,ref_2d_04,istatus_nowrad)

        if(istatus_nowrad .ne. 1 .and. istatus_nowrad .ne. -1)then
            write(6,*)' Could not read in ',var_2d
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

        var_2d = 'R48'
        ext = 'ln3'
        call get_laps_2d(i4time_radar,ext,var_2d
     1       ,units_2d,comment_2d,imax,jmax,ref_2d_48,istatus_nowrad)

        if(istatus_nowrad .ne. 1 .and. istatus_nowrad .ne. -1)then
            write(6,*)' Could not read in ',var_2d
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

        var_2d = 'R8C'
        ext = 'ln3'
        call get_laps_2d(i4time_radar,ext,var_2d
     1       ,units_2d,comment_2d,imax,jmax,ref_2d_8c,istatus_nowrad)

        if(istatus_nowrad .ne. 1 .and. istatus_nowrad .ne. -1)then
            write(6,*)' Could not read in ',var_2d
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

        var_2d = 'ET'
        ext = 'ln3'
        call get_laps_2d(i4time_radar,ext,var_2d
     1       ,units_2d,comment_2d,imax,jmax,ref_2d_et,istatus_nowrad)

        if(istatus_nowrad .ne. 1 .and. istatus_nowrad .ne. -1)then
            write(6,*)' Could not read in ',var_2d
            istatus_2dref_a = 0
            istatus_3dref_a = 0
            return
        endif

!       Now we process the layer reflectivities into the 3-D grid
        do k = 1,kmax
        do j = 1,jmax
        do i = 1,imax

            if    (heights_3d(i,j,k) .le. 3000.  )then

                frac1 = 1.0
                frac2 = 0.0
                frac3 = 0.0

            elseif(heights_3d(i,j,k) .ge. 3000. .and.
     1             heights_3d(i,j,k) .le. 5000.  )then

                frac2 = (heights_3d(i,j,k) - 3000.) / 2000.
                frac1 = 1.0 - frac2
                frac3 = 0.0

            elseif(heights_3d(i,j,k) .ge. 5000. .and.
     1             heights_3d(i,j,k) .le. 7000.  )then

                frac1 = 0.0
                frac2 = 1.0
                frac3 = 0.0

            elseif(heights_3d(i,j,k) .ge. 7000. .and.
     1             heights_3d(i,j,k) .le. 9000.  )then

                frac3 = (heights_3d(i,j,k) - 7000.) / 2000.
                frac2 = 1.0 - frac3
                frac1 = 0.0

            else ! heights_3d(i,j,k) .ge. 9000.

                frac1 = 0.0
                frac2 = 0.0
                frac3 = 1.0

            endif

            if(ref_2d_04(i,j) .eq. r_missing_data)then
                ref_2d_04(i,j) = ref_base
            endif

            if(ref_2d_48(i,j) .eq. r_missing_data)then
                ref_2d_48(i,j) = ref_base
            endif

            if(ref_2d_8c(i,j) .eq. r_missing_data)then
                ref_2d_8c(i,j) = ref_base
            endif

            ref_3d(i,j,k) = ref_2d_04(i,j) * frac1
     1                    + ref_2d_48(i,j) * frac2
     1                    + ref_2d_8c(i,j) * frac3

            if(ref_3d(i,j,k) .gt. 0)then
                icount_ref1 = icount_ref1 + 1
            endif

!           Determine effective echo top
            if(ref_2d_et(i,j) .eq. r_missing_data)then
                echo_top = 10000.
            elseif(ref_2d_et(i,j) .ge. 10000.)then
                echo_top = ref_2d_et(i,j)
            else
                echo_top = 10000.
            endif

!           Clear echo out if we are above echo top
            if(heights_3d(i,j,k) .gt. echo_top)then ! set echo to ref_base
                ref_3d(i,j,k) = -10.
            endif

            if(ref_3d(i,j,k) .gt. 0)then
                icount_ref2 = icount_ref2 + 1
            endif

        enddo ! i
        enddo ! j
        enddo ! k

        write(6,*)' # of ref points > 0 before/after echo top check = '
     1            ,icount_ref1,icount_ref2

        istatus_2dref_a = 1
        istatus_3dref_a = 0

        return
        end
