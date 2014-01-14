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
cdoc    read_vrz_3dref
cdoc        now called from read_multiradar_3dref
cdoc
cdoc
!       1996 Aug    S. Albers FSL


        subroutine get_multiradar_vel(
     1   i4time_ref,i4time_tol,i4time_radar_a
     1  ,max_radars,n_radars,ext_a,r_missing_data
     1  ,imax,jmax,kmax,lat,lon                                      ! I
     1  ,nx_r,ny_r,igrid_r                                           ! I
     1  ,grid_ra_vel,grid_ra_nyq,idx_radar,v_nyquist_in_a
     1  ,ioffset,joffset                                             ! O
     1  ,l_offset_radar                                              ! I
     1  ,n_vel_a
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
!       imax,jmax,kmax      Input   LAPS 3D Grid Dimensions
!       lat,lon             Input   2D latitude and longitude arrays (degrees)
!       topo                Input   2D terrain array (meters)
!       grid_ra_vel         Output  4D Velocity Grid
!       grid_ra_nyq         Output  4D Nyquist Velocity Grid
!       idx_radar           Output  1D grid of radar 'vxx' numbers
!       v_nyquist_in_a      Output  Array: volume nyquist velocity of the radars
!       i_offset,joffset    Input   offset arrays (1D)
!       l_offset_radar      Input   use offset arrays?
!       n_vel_a             Output  Array: # of grid points with measurable velocity
!       rlat_radar_a        Output  Array: Radar Latitude (Degrees)
!       rlon_radar_a        Output  Array: Radar Longitude (Degrees)
!       rheight_radar_a     Output  Array: Radar Height (m MSL)
!       radar_name_a        Output  Array: Radar Name (Character*4)
!       istatus_multi_vel   Output  If 1 - data is useable for 3D vel applications
!       istatus_multi_nyq   Output  If 1 - data is useable for 3D nyq applications

!
!       The domain must be specified using grid_fnam_common

        character*31 ext_a(max_radars)

        character*255 c_filespec
        character*150  directory

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        character asc9_tim_radar*9

        Integer       ioffset(max_radars)
        Integer       joffset(max_radars)
        logical l_apply_map, l_offset_radar

        real grid_ra_vel(nx_r,ny_r,kmax,max_radars)
        real grid_ra_nyq(nx_r,ny_r,kmax,max_radars)
        real grid_ra_vel_3d(imax,jmax,kmax)
        real grid_ra_nyq_3d(imax,jmax,kmax)
        integer idx_radar(max_radars)

        real lat(imax,jmax),lon(imax,jmax)

        real rlat_radar_a(max_radars),rlon_radar_a(max_radars)
     1    ,rheight_radar_a(max_radars),v_nyquist_in_a(max_radars)

        character*4 radar_name_a(max_radars)

        integer n_vel_a(max_radars)
        integer i4time_radar_a(max_radars)

        logical l_conus_vel

        l_conus_vel = .false.

!       Initialize the flags
        n_radars = 0
        n_ref = 0
        istatus_multi_vel = 1
        istatus_multi_nyq = 1
        idx_radar = 0

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

                idx_radar(n_radars) = i_radar_pot

                grid_ra_vel(:,:,:,n_radars) = r_missing_data
                grid_ra_nyq(:,:,:,n_radars) = r_missing_data

                i4time_radar_a(n_radars) = i4time_radar

                write(6,*)' Reading vel/nyq for actual radar # '
     1                    ,n_radars

                call read_radar_vel(i4time_radar,l_apply_map,
     1           imax,jmax,kmax,ext_a(i_radar_pot),
     1           grid_ra_vel_3d,
     1           grid_ra_nyq_3d,v_nyquist_in_a(n_radars),       
     1           rlat_radar_a(n_radars),rlon_radar_a(n_radars)
     1                      ,rheight_radar_a(n_radars)
     1                      ,radar_name_a(n_radars)
     1                      ,n_vel_a(n_radars),
     1                       istatus_vel,istatus_nyq)

                if(n_vel_a(n_radars) .eq. 0 .or. istatus_vel .ne. 1)then       
                    write(6,*)' No valid velocities for radar ',n_radars
                  ! Don't count in a valid radar
                    idx_radar(n_radars) = 0
                    n_radars = n_radars - 1

                else ! valid radar

                    if(l_offset_radar)then
                      call get_ij_offset_radars(imax,jmax,1,              ! I
     1                              igrid_r,l_offset_radar,               ! I   
     1                              lat,lon,                              ! I
     1                              rlat_radar_a(n_radars),               ! I
     1                              rlon_radar_a(n_radars),               ! I
     1                              ioffset(n_radars),                    ! I
     1                              joffset(n_radars) )                   ! O  


!                     if(rlat_radar_a(n_radars) .eq. r_missing_data .or.
!    1                   rlon_radar_a(n_radars) .eq. r_missing_data  
!    1                                                             )then
!                       write(6,*)
!    1                        ' No valid or single lat/lon for radar '       
!    1                           ,n_radars       
!    1                           ,' ',radar_name(n_radars)
!                       l_valid_latlon(n_radars) = .false.

!                       ioffset(n_radars) = 0
!                       joffset(n_radars) = 0

!                     else
!                       call latlon_to_rlapsgrid(rlat_radar_a(n_radars),
!    &                                           rlon_radar_a(n_radars),
!    &                                           lat,lon,
!    &                                           imax,jmax,
!    &                                           ri,rj,
!    &                                           jstatus)
!                       if(jstatus.ne.1)then
!                           write(6,*)
!    1                   'computing ri/rj for radar (outside domain)'    
!                       endif
!                       write(6,*)'Name: ',radar_name(n_radars)
!     1                          ,ri(n_radars),rj(n_radars),n_radars
!                       l_valid_latlon(n_radars) = .true.

!                       Offset is location of lower left corner of small array in the large array
!                       ioffset(n_radars) = (nint(ri) - igrid_r) - 1
!                       joffset(n_radars) = (nint(rj) - igrid_r) - 1
!                     endif

                      I4_elapsed = ishow_timer()

                      write(6,*)' get_multiradar_vel - offset info '
     1           ,'ri,rj,ioffset(n_radars),joffset(n_radars),igrid_r : '
     1            ,ri,rj,ioffset(n_radars),joffset(n_radars),igrid_r

                      nfill = 0  

                      do jo = 1,ny_r
                        j = jo + joffset(n_radars)
                        if(j .ge. 1 .and. j .le. jmax)then
                          do io = 1,nx_r
                            i = io + ioffset(n_radars)

                            if(i .ge. 1 .and. i .le. imax)then
                              grid_ra_vel(io,jo,:,n_radars) = 
     1                        grid_ra_vel_3d(i,j,:)

                              grid_ra_nyq(io,jo,:,n_radars) = 
     1                        grid_ra_nyq_3d(i,j,:)
                              nfill = nfill + 1
                            endif ! in i bounds

                          enddo ! io
                        endif ! in j bounds
                      enddo ! j

                      write(6,*)' nfill = ',nfill

                    else ! fill 4D array with 3D array contents
                      grid_ra_vel(:,:,:,n_radars) = grid_ra_vel_3d
                      grid_ra_nyq(:,:,:,n_radars) = grid_ra_nyq_3d

                      ioffset(n_radars) = 0
                      joffset(n_radars) = 0

                    endif ! l_offset_radar

                endif ! valid radar

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

        real grid_ra_ref(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)
        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

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

                i4_tol = 1200

                call read_radar_3dref(i4time_radar,                   ! I
     1               i4_tol,i4_ret,                                   ! I/O
     1               l_apply_map,ref_base,                            ! I
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
     1   i4_tol,i4_ret,                                                  ! I/O
     1   l_apply_map,ref_missing,                                        ! I
     1   imax,jmax,kmax,radarext,                                        ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                            ! I
     1   heights_3d,                                                     ! I
     1   grid_ra_ref,                                                    ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                 ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                        ! O

!       Now a jacket routine without 'closest_radar' in the argument list

        real grid_ra_ref(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

        real closest_radar(imax,jmax)

        character*31 radarext
        character*4 radar_name

        logical l_low_fill,l_high_fill,l_apply_map

        call read_radar_3dref_new(i4time_radar,                          ! I
     1   i4_tol,i4_ret,                                                  ! I/O
     1   l_apply_map,ref_missing,                                        ! I
     1   imax,jmax,kmax,radarext,                                        ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                            ! I
     1   heights_3d,                                                     ! I
     1   grid_ra_ref,                                                    ! O
     1   closest_radar,                                                  ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                 ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                        ! O

        return
        end

        subroutine read_radar_3dref_new(i4time_radar,                    ! I
     1   i4_tol,i4_ret,                                                  ! I/O
     1   l_apply_map,ref_missing,                                        ! I
     1   imax,jmax,kmax,radarext,                                        ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                            ! I
     1   heights_3d,                                                     ! I
     1   grid_ra_ref,                                                    ! O
     1   closest_radar,                                                  ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                 ! O
     1   n_ref_grids,istatus_2dref,istatus_3dref)                        ! O

cdoc    Steve Albers Feb 1998   This routine will read in a 3D radar
cdoc                            reflectivity field. It is a jacket that
cdoc                            calls read_multiradar_3dref.


        real grid_ra_ref(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

        real closest_radar(imax,jmax)

!       If a grid "overall" has the following...
!       2dref=1, 3dref=1 - echo top from radar has better confidence than
!                          the associated cloud top from satellite
!       2dref=1, 3dref=0 - echo top from radar has less confidence
!       2dref=0, 3dref=0 - missing data

        integer istatus_2dref
        integer istatus_3dref
        integer istatus_2dref_a(imax,jmax)
        integer istatus_3dref_a(imax,jmax)

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
     1   i4_tol,i4_ret,
     1   l_apply_map,ref_missing,
     1   imax,jmax,kmax,radarext,
     1   lat,lon,topo,l_low_fill,l_high_fill,
     1   heights_3d,
     1   grid_ra_ref,
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,
     1   iqc_2dref,closest_radar,                                       ! O
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
     1   i4_tol,i4_ret,                                                 ! I
     1   l_apply_map,ref_missing,                                       ! I
     1   imax,jmax,kmax,radarext,                                       ! I
     1   lat,lon,topo,l_low_fill,l_high_fill,                           ! I
     1   heights_3d,                                                    ! I
     1   grid_ra_ref,                                                   ! O
     1   rlat_radar,rlon_radar,rheight_radar,radar_name,                ! O
     1   iqc_2dref,closest_radar,                                       ! O
     1   n_ref_grids,n_2dref,n_3dref,istatus_2dref_a,istatus_3dref_a)   ! O 

!       Steve Albers Nov 1998   This routine will read in a 3D radar
!                               reflectivity field. It will either read
!                               in 3D radar refs and do the vertical filling
!                               procedure, or read in the 2D data and do
!                               a quick and dirty vertical fill. Information
!                               is passed back about which data is 2d and which
!                               is 3d.
!
!       SA           Jan 2003   The 'all' option is now being tested that 
!                               merges 3D (vrz) and 2D (vrc) data.
!                               The 'a01' option should merge (v01) with (vrc)

        real grid_ra_ref(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)
        real radar_2dref(imax,jmax)
        real closest_vxx(imax,jmax)
        real closest_vrc(imax,jmax)
        real closest_radar(imax,jmax)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

!       If a grid element has the following...
!       2dref=1, 3dref=1 - echo top from radar has better confidence than
!                          the associated cloud top from satellite
!       2dref=1, 3dref=0 - echo top from radar has less confidence
!       2dref=0, 3dref=0 - missing data

!       iqc_2dref = 0 ! less confidence in 2dref data QC
!       iqc_2dref = 1 ! more confidence in 2dref data (e.g. from NOWRAD)

        integer istatus_2dref
        integer istatus_3dref
        integer istatus_2dref_a(imax,jmax)
        integer istatus_3dref_a(imax,jmax)

        character*3 var_2d
        character*31  readext
        character*10  units_2d
        character*125 comment_2d

        character*31 radarext

        character*4 radar_name

        logical l_low_fill,l_high_fill,l_apply_map, l_parse

        write(6,*)' Subroutine read_multiradar_3dref'

        istatus_2dref_a = 0
        istatus_3dref_a = 0
        iqc_2dref = 0

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
        closest_radar = r_missing_data

!       Initialize 3d reflectivity array with default value
        grid_ra_ref = r_missing_data ! ref_base 

        if(radarext(1:2) .eq. 'v0' .or.
     1     radarext(1:2) .eq. 'v1' .or.
     1     radarext(1:2) .eq. 'v2' .or.
     1     radarext(1:2) .eq. 'v3' .or.
     1     radarext(1:2) .eq. 'v4' .or.
     1     radarext(1:2) .eq. 'v5' .or.
     1     radarext(1:2) .eq. 'v6' .or.
     1     radarext(1:2) .eq. 'v7' .or.
     1     radarext(1:2) .eq. 'v8' .or.
     1     radarext(1:2) .eq. 'v9' .or.
     1     radarext(1:3) .eq. 'vrz' .or.
     1     radarext(1:1) .eq. 'a'   )then     ! Read Doppler 3-D radar ref data
                                              ! from NetCDF VXX or VRZ file

            write(6,*)' Reading Reflectivity Data from 3D file '
     1                                                 ,radarext

            if(radarext(1:3) .eq. 'all')then
                readext = 'vrz' 
            elseif(radarext(1:3) .eq. 'a01')then
                readext = 'v01' 
            else
                readext = radarext
            endif

!           Read Reflectivity
            var_2d = 'REF'

            if(readext(1:3) .ne. 'vrz')then ! read vxx file
                call get_laps_3dgrid(i4time_radar,i4_tol,i4_ret
     1                              ,imax,jmax,kmax,readext,var_2d
     1                              ,units_2d,comment_2d,grid_ra_ref
     1                              ,istatus)

            else                            ! read vrz file
                call get_max_radars(max_radars,istatus)
                if(istatus .ne. 1)return

                if(i4_tol .lt. 0)then ! Bypass reading of vrz
                    istatus = 0
                else
                    call read_vrz_3dref(i4time_radar                     ! I
     1                                 ,imax,jmax,kmax,max_radars        ! I
     1                                 ,readext,var_2d                   ! I
     1                                 ,lat,lon,topo                     ! I
     1                                 ,grid_ra_ref                      ! O
     1                                 ,closest_vxx                      ! O
     1                                 ,istatus)                         ! O
                    closest_radar = closest_vxx ! initialize

                endif
 
            endif

            if(istatus .eq. 1)then
                if(readext(1:3) .ne. 'vrz')then ! vxx
                    read(comment_2d,558)
     1                             rlat_radar,rlon_radar,rheight_radar       
     1                            ,n_ref_grids,radar_name
558                 format(2f9.3,f8.0,i7,a4)

                    if(l_low_fill .or. l_high_fill)then
                        I4_elapsed = ishow_timer()
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

                    endif ! vert_fill

                    write(6,*)' Read radar ',radar_name,' Volume'

                else ! vrz mosaic
                    n_ref_grids = 0
                    rlat_radar = r_missing_data
                    rlon_radar = r_missing_data
                    rheight_radar = r_missing_data
                    radar_name = 'MOSC' 

                endif

!               if(.false.)then
!                   call ref_fill_horz(grid_ra_ref,imax,jmax,kmax
!    1                ,lat,lon,rlat_radar,rlon_radar,rheight_radar
!    1                ,istatus)
!                   if(istatus .ne. 1)then
!                       istatus_2dref_a = 0
!                       istatus_3dref_a = 0
!                   endif
!               endif

!               Set status flags and closest 'vxx'
!               if(l_low_fill .or. l_high_fill .or. 
!    1             radarext(1:3) .eq. 'vrz'                )then
                if(.true.)then
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

                            if(readext(1:3) .ne. 'vrz')then ! vxx radar case
                                call latlon_to_radar(
     1                              lat(i,j),lon(i,j),topo(i,j)
     1                             ,azimuth,closest_vxx(i,j),elev
     1                             ,rlat_radar,rlon_radar,rheight_radar)       
                                closest_radar(i,j) = closest_vxx(i,j) ! initialize
                            endif

                        endif

                    enddo ! i
                    enddo ! j

                endif

            else
                write(6,*)' Radar reflectivity data cannot be read in: '       
     1                    ,readext

            endif ! Success as reflectivity

        endif ! 'vxx' or 'vrz'

        if(radarext(1:3) .eq. 'vrc' .or. radarext(1:1) .eq. 'a')then       
50          write(6,*)' Reading NOWRAD/vrc data' ! lumped together for now?

            readext = 'vrc'
            var_2d = 'REF'
            call get_laps_2dgrid(i4time_radar,i4_tol,i4_ret,readext
     1                          ,var_2d,units_2d,comment_2d,imax,jmax
     1                          ,radar_2dref,0,istatus_vrc)

            write(6,*)' istatus_vrc = ',istatus_vrc

            if(istatus_vrc .eq. 1 .or. istatus_vrc .eq. -1)then       
                var_2d = 'DIS'
                call get_laps_2d(i4_ret,readext
     1                          ,var_2d,units_2d,comment_2d,imax,jmax      
     1                          ,closest_vrc,istatus_dis)

!               closest_vrc = 180000. ! This can be activated by commenting out

                if(istatus_dis .ne. 1 .and. istatus_dis .ne. -1)then
                    write(6,*)' ERROR in get_radar: istatus_dis = '
     1                       ,istatus_dis
                    goto900
                endif

                if(l_parse(comment_2d,'WSI'))then
                    radar_name = 'WSI '
                    write(6,*)' Read radar ',radar_name
     1                       ,' Low Level Mosaic'    
                    iqc_2dref = 1 ! better quality QC
                elseif(l_parse(comment_2d,'Radar mosaic'))then
                    radar_name = '    '
                    write(6,*)
     1                ' Read radar Narrowband Mosaic: # radars = '
     1                ,comment_2d(26:28)
                else
                    len_comment = 37
                    radar_name = comment_2d(len_comment-3:len_comment)       
                    write(6,*)' Read radar ',radar_name,' Low Level'
                endif 
                    
                do i = 1,imax
                do j = 1,jmax
                    if(radar_2dref(i,j) .ne. r_missing_data)then
                        istatus_2dref_a(i,j) = 1
                    endif
                enddo ! j
                enddo ! i                 

!               Conditionally fill up the 3D reflectivity array 
!                          (useful for precip type, get_low_ref)

!               if(l_low_fill .or. l_high_fill)then
                if(.true.)then
                    
                    do i = 1,imax
                    do j = 1,jmax

                        closest_radar(i,j) = closest_vxx(i,j)

!                       Set distant 3D radar points to msg to select NOWRAD/vrc
                        if(closest_vxx(i,j) .ne. r_missing_data
     1               .and. closest_vxx(i,j) .gt. 
     1                                   closest_vrc(i,j) + 25000.)then       
!    1               .and. closest_vxx(i,j) .gt. closest_vrc(i,j) )then       
                            istatus_3dref_a(i,j) = 0
                        endif

                        if(istatus_3dref_a(i,j) .eq. 0)then 
                            do k = 1,kmax ! Fill column with NOWRAD/vrc
!                               if(grid_ra_ref(i,j,k).eq.r_missing_data      
!    1                                                             )then       
                                    grid_ra_ref(i,j,k)=radar_2dref(i,j)      
!                               endif
                            enddo ! k

                            closest_radar(i,j) = closest_vrc(i,j)

                        endif ! istatus_3dref_a = 0

                    enddo ! j
                    enddo ! i

                endif ! l_low_fill OR l_high_fill

            endif ! good vrc status

        endif ! vrc

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

        pct_radar_coverage =  float(n_2dref)/(imax*jmax) * 100.
        pct_volume_coverage = float(n_3dref)/(imax*jmax) * 100.

        write(6,901)pct_radar_coverage,pct_volume_coverage
 901    format(' Radar coverage is '      ,f6.2,'%    '
     1        ,' Full volume coverage is ',f6.2,'%')

        return
        end



        subroutine read_radar_2dref(i4time_radar,radar_name,
     1                  imax,jmax,
     1                  ref_2d,istatus_2dref)

!       Steve Albers Feb 1996   Reads Radar data into 2d REF arrays
!           now called from (wind-derived via) get_radar_max_pd


        real ref_2d(imax,jmax)

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

        subroutine read_radar_raw2d(i4time_radar,             ! I
     1                  imax,jmax,                            ! I
     1                  iradar                                ! I
     1                  radar_name,                           ! O
     1                  ref_2d,radar_dist_2d,istatus_2dref)   ! O

!       Read in unmosaiced radar data from 'vrc' files (under construction)

        real ref_2d(imax,jmax)
        real radar_dist_2d(imax,jmax) ! Not yet filled

        character*3 var_2d
        character*31  ext
        character*10  units_2d
        character*125 comment_2d

        character*4 radar_name

        write(6,*)' Subroutine read_radar_raw2d'

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


        real grid_ra_vel(imax,jmax,kmax)
        real grid_ra_nyq(imax,jmax,kmax)

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


        subroutine get_radar_horizon(i,j,ni,nj,lat,lon,topo)

        include 'trigd.inc'

        use mem_namelist, ONLY: grid_spacing_m

!       Consider line segment from radar location to present location
!       This version approximates the great circle from the radar to the 
!       grid point using a straight line in map projection space

        real lat(ni,nj)
        real lon(ni,nj)
        real topo(ni,nj)

        ri_radar = 10.
        rj_radar = 50.
        rheight_radar = 1000.

        ri_grid = i
        rj_grid = j

        dist = sqrt ( (ri_grid-ri_radar)**2 + (rj_grid-rj_radar)**2 )
        dir = atan3(ri_grid-ri_radar,rj_grid-rj_radar) 
        sindir = sind(dir)
        cosdir = cosd(dir)

        delta_dist = 1.

        do ds = 0.,dist,delta_dist
            dist_i = ds * cosdir
            dist_j = ds * sindir

            ri = ri_radar + dist_i
            rj = rj_radar + dist_j

            i = nint(ri)
            j = nint(rj)

            rheight = topo(i,j)

            call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                          ,azimuth,slant_range,elev_angle
     1                          ,rlat_radar,rlon_radar,rheight_radar)

            angle_max = max(angle_max,elev_angle)

        enddo ! ds                

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
        real grid_ra_ref(imax,jmax,kmax)
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

        real ref_3d(imax,jmax,kmax)
        real heights_3d(imax,jmax,kmax)

        real ref_2d_04(imax,jmax)
        real ref_2d_48(imax,jmax)
        real ref_2d_8c(imax,jmax)
        real ref_2d_et(imax,jmax)

        integer istatus_2dref_a(imax,jmax)
        integer istatus_3dref_a(imax,jmax)

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


        subroutine read_vrz_3dref(i4time_radar                      ! I
     1                           ,imax,jmax,kmax,max_radars         ! I
     1                           ,ext,var_2d                        ! I
     1                           ,lat,lon,topo                      ! I
     1                           ,ref_3d                            ! O
     1                           ,closest_vxx                       ! O
     1                           ,istatus)                          ! O

        character*(*) EXT, var_2d

        real ref_3d(imax,jmax,kmax)
        real closest_vxx(imax,jmax)

        real rlat_radar(max_radars)
        real rlon_radar(max_radars)
        real rheight_radar(max_radars)

        real lat(imax,jmax)
        real lon(imax,jmax)
        real topo(imax,jmax)

        character c_radar_id(max_radars)*4
        character*9 asc9_tim
        character*150 directory

        character*200 comment_3d(kmax)
        character*40 comment_tmp
        character*10 units_3d(kmax)           
        character*3 var_3d(kmax)
        integer LVL_3d(kmax)
        character*4 LVL_COORD_3d(kmax)

        logical ltest_vertical_grid

        write(6,*)' Subroutine read_vrz_3dref'

        call get_directory(ext,directory,len_dir)

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' read_vrz_3dref: bad istatus, return'
            return
        endif

        call s_len(ext,len_ext)

        call make_fnam_lp(i4time_radar,asc9_tim,istatus)

        len_high = max(45,len_dir)
        len_low = len_high - 44

        write(6,11)directory(len_low:len_high),asc9_tim,ext,var_2d
11      format('  read_vrz_3dref: ',a,1x,a,1x,a,1x,a)

        do k = 1,kmax
            if(ltest_vertical_grid('HEIGHT'))then
                lvl_3d(k) = zcoord_of_level(k)/10
                lvl_coord_3d(k) = 'MSL'
            elseif(ltest_vertical_grid('PRESSURE'))then
                lvl_3d(k) = nint(zcoord_of_level(k))/100
                lvl_coord_3d(k) = 'MB'
            else
                write(6,*)' Error, vertical grid not supported,'
     1                   ,' this routine supports PRESSURE or HEIGHT'
                istatus = 0
                return
            endif

            var_3d(k) = var_2d

        enddo ! k

        CALL READ_LAPS_DATA(i4time_radar,DIRECTORY,EXT,imax,jmax,
     1  kmax,kmax,VAR_3D,LVL_3D,LVL_COORD_3D,UNITS_3D,
     1                     COMMENT_3D,ref_3d,ISTATUS)
        if(istatus .ne. 1)return

!       Read radar info from comments
        read(comment_3d(1),101)n_radars
 101    format(25x,i3)

        if(n_radars .le. max_radars)then
            write(6,*)' n_radars/max_radars = ',n_radars,max_radars       
        else
            write(6,*)' ERROR in read_vrz_3dref: n_radars > max_radars'       
     1           ,n_radars,max_radars
            istatus = 0
            return
        endif

!       Read comments from 3 columns each 40 characters wide
        nch=40
        do i_radar = 1,n_radars
            if(i_radar .le. (kmax-1) )then       ! read in 1st column
                ii = i_radar + 1
                read(comment_3d(ii),1)rlat_radar(i_radar)
     1                               ,rlon_radar(i_radar)
     1                               ,rheight_radar(i_radar)
     1                               ,n_ref
     1                               ,c_radar_id(i_radar)
1               format(2f9.3,f8.0,i7,a4)

                write(6,*)' Read radar ',c_radar_id(i_radar)
     1                   ,' Volume (via 3d-mosaic)'

            elseif(i_radar .le. 2*(kmax-1) )then ! read in 2nd column
                ii = i_radar - (kmax-1) + 1
                comment_tmp = comment_3d(ii)(1*nch+1:nch*2)
                read(comment_tmp,1)rlat_radar(i_radar)
     1                            ,rlon_radar(i_radar)
     1                            ,rheight_radar(i_radar)
     1                            ,n_ref
     1                            ,c_radar_id(i_radar)

                write(6,*)' Read radar ',c_radar_id(i_radar)
     1                   ,' Volume (via 3d-mosaic)'

            elseif(i_radar .le. 3*(kmax-1) )then ! read in 3rd column
                ii = i_radar - (2*(kmax-1)) + 1
                comment_tmp = comment_3d(ii)(2*nch+1:nch*3)
                read(comment_tmp,1)rlat_radar(i_radar)
     1                            ,rlon_radar(i_radar)
     1                            ,rheight_radar(i_radar)
     1                            ,n_ref
     1                            ,c_radar_id(i_radar)

                write(6,*)' Read radar ',c_radar_id(i_radar)
     1                   ,' Volume (via 3d-mosaic)'

            elseif(i_radar .le. 4*(kmax-1) )then ! read in 4th column
                ii = i_radar - (3*(kmax-1)) + 1
                comment_tmp = comment_3d(ii)(3*nch+1:nch*4)
                read(comment_tmp,1)rlat_radar(i_radar)
     1                            ,rlon_radar(i_radar)
     1                            ,rheight_radar(i_radar)
     1                            ,n_ref
     1                            ,c_radar_id(i_radar)

                write(6,*)' Read radar ',c_radar_id(i_radar)
     1                   ,' Volume (via 3d-mosaic)'

            elseif(i_radar .le. 5*(kmax-1) )then ! read in 5th column
                ii = i_radar - (4*(kmax-1)) + 1
                comment_tmp = comment_3d(ii)(4*nch+1:nch*5)
                read(comment_tmp,1)rlat_radar(i_radar)
     1                            ,rlon_radar(i_radar)
     1                            ,rheight_radar(i_radar)
     1                            ,n_ref
     1                            ,c_radar_id(i_radar)

                write(6,*)' Read radar ',c_radar_id(i_radar)
     1                   ,' Volume (via 3d-mosaic)'

            else
                write(6,*)
     1          ' Error: too many radars for comment output'
                write(6,*)' Limit is ',i_radar-1
                istatus = 0
                return

            endif

        enddo ! i_radar

!       Create closest radar array (assumes lat/lon projection) 
        do j=1,jmax
        do i=1,imax
           distmin=r_missing_data
           do k=1,n_radars

              if(.false.)then
                  rlatdif=(lat(i,j)-rlat_radar(k))*111100.      !m
                  rlondif=(lon(i,j)-rlon_radar(k))*111100.
                  dist=sqrt(rlatdif*rlatdif + rlondif*rlondif)

              else
                  call latlon_to_radar(
     1                             lat(i,j),lon(i,j),topo(i,j)
     1                            ,azimuth,dist,elev
     1                            ,rlat_radar(k),rlon_radar(k)
     1                            ,rheight_radar(k))       
              endif

              if(dist.lt.distmin)distmin=dist
           enddo
           closest_vxx(i,j)=distmin
        enddo
        enddo

        return
        end

        subroutine get_l_offset_radar(nx_l,ny_l,grid_spacing_cen_m,     ! I
     1                                nx_r,ny_r,igrid_r,l_offset_radar) ! O    

        use mem_namelist, only: i_offset_radar

        logical l_offset_radar

        radius_r = 500000. ! 500km max radar radius
        igrid_r = int(radius_r / grid_spacing_cen_m) + 1
        nx_r_pot = (2 * igrid_r) + 1
        ny_r_pot = (2 * igrid_r) + 1
        write(6,1)i_offset_radar,nx_r_pot,ny_r_pot,nx_l,ny_l
 1      format(' Potential offset radar arrays:',i3,2i6,3x,2i6)

        if(i_offset_radar .eq. 1)then
            l_offset_radar = .true.
            if(nx_r_pot**2 .gt. nx_l*ny_l)then ! offset arrays are larger                     
                write(6,*)' Warning: offset arrays are larger than grid'
            endif
            nx_r = nx_r_pot
            ny_r = ny_r_pot
        elseif(i_offset_radar .eq. 0)then
            if(nx_r_pot**2 .gt. nx_l*ny_l)then ! offset arrays are larger                     
                write(6,*)' Offset arrays are larger than grid'
                nx_r = nx_l
                ny_r = ny_l
                igrid_r = 0 ! dummy value
                l_offset_radar = .false.
            else
                write(6,*)' Offset arrays are within grid'
                nx_r = nx_r_pot
                ny_r = ny_r_pot
                l_offset_radar = .true.
            endif
        else ! i_offset_radar = -1
            nx_r = nx_l
            ny_r = ny_l
            igrid_r = 0 ! dummy value
            l_offset_radar = .false.
        endif

        write(6,*)' l_offset_radar = ',l_offset_radar,nx_r,ny_r

        return
        end


        subroutine get_ij_offset_radars(nx_l,ny_l,n_radars,               ! I
     1                                  igrid_r,l_offset_radar,           ! I   
     1                                  lat,lon,rlat_radar,rlon_radar,    ! I
     1                                  ioffset,joffset)                  ! O

        real lat(nx_l,ny_l),lon(nx_l,ny_l)
        real rlat_radar(n_radars),rlon_radar(n_radars)
        integer ioffset(n_radars),joffset(n_radars)
        logical l_offset_radar

        do k = 1,n_radars
          if(l_offset_radar)then
            if(rlat_radar(k) .eq. r_missing_data .or.
     1         rlon_radar(k) .eq. r_missing_data      )then

                ioffset(k) = 0
                joffset(k) = 0

            else
                call latlon_to_rlapsgrid(rlat_radar(k),
     &                                   rlon_radar(k),
     &                                   lat,lon,
     &                                   nx_l,ny_l,
     &                                   ri,rj,
     &                                   jstatus)
                if(jstatus.ne.1)then
                    write(6,*)
     1               'computing ri/rj for radar (outside domain)'
                endif

!               Offset is location of lower left corner of small array in the large array
                ioffset(k) = (nint(ri) - igrid_r) - 1
                joffset(k) = (nint(rj) - igrid_r) - 1

            endif

          else
            ioffset(k) = 0
            joffset(k) = 0

          endif

        enddo ! k

        return
        end
