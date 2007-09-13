

        subroutine wind_post_process(i4time_sys,EXT,var_3d
     1                              ,units_3d,comment_3d
     1                              ,uanl,vanl                            ! I
     1                              ,NX_L,NY_L,NZ_L                       ! I
     1                              ,N_3D_FIELDS                          ! I
     1                              ,uanl_sfcitrp,vanl_sfcitrp            ! I
     1                              ,topo,lat,lon,grid_spacing_m          ! I
     1                              ,rk_terrain                           ! I
     1                              ,r_missing_data,l_grid_north_out      ! I
     1                              ,istat_lw3)

        real uanl(NX_L,NY_L,NZ_L),vanl(NX_L,NY_L,NZ_L) ! WRT True North ! I
        real wanl(NX_L,NY_L,NZ_L)                                       ! L
        real uanl_sfcitrp(NX_L,NY_L),vanl_sfcitrp(NX_L,NY_L)            ! I

        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)

        real rk_terrain(NX_L,NY_L)

        character*125 comment_3D(N_3D_FIELDS)
        character*10 units_3D(N_3D_FIELDS)
        character*3 var_3D(N_3D_FIELDS)
        character*3 EXT

        logical l_grid_north_out

csms$ignore begin
        write(6,*)' Subroutine wind_post_process...'

        write(6,*)' Computing Omega'
        call vert_wind(uanl,vanl,uanl_sfcitrp,vanl_sfcitrp                ! I
     1                ,NX_L,NY_L,NZ_L                                     ! I
     1                ,wanl                                               ! O
     1                ,topo,lat,lon,grid_spacing_m                        ! I
     1                ,rk_terrain,r_missing_data,l_grid_north_out         ! I
     1                ,istatus)                                           ! O

        if(istatus .ne. 1)then
            write(6,*)' Error: bad data detected by vert_wind'
            write(6,*)
     1      ' Check for missing data on one or more levels of uanl/vanl'
            return
        endif

        I4_elapsed = ishow_timer()


!       Header information for 3D wind
        EXT = 'lw3'

        units_3D(1) = 'M/S'
        units_3D(2) = 'M/S'
        units_3D(3) = 'PA/S'

        var_3D(1) = 'U3'
        var_3D(2) = 'V3'
        var_3D(3) = 'OM'

        comment_3D(1) = '3DWIND'
        comment_3D(2) = '3DWIND'
        comment_3D(3) = '3DWIND'

        I4_elapsed = ishow_timer()

        write(6,*)' Calling write routine for all grids ',ext(1:3)
     1                                  ,i4time_sys

!       call move_3d(uanl,outarray_4D(1,1,1,1),NX_L,NY_L,NZ_L)
!       call move_3d(vanl,outarray_4D(1,1,1,2),NX_L,NY_L,NZ_L)
!       call move_3d(wanl,outarray_4D(1,1,1,3),NX_L,NY_L,NZ_L)

!       call put_laps_multi_3d(i4time_sys,EXT,var_3d,units_3d,
!    1     comment_3d,outarray_4D,NX_L,NY_L,NZ_L,N_3D_FIELDS,istat_lw3)

        call put_laps_multi_3d_jacket(i4time_sys,EXT,var_3d
     1                               ,units_3d,comment_3d
     1                               ,uanl,vanl,wanl
     1                               ,NX_L,NY_L,NZ_L,N_3D_FIELDS
     1                               ,istat_lw3)

csms$ignore end
        return
        end


        subroutine put_laps_multi_3d_jacket(i4time_sys,EXT,var_3d
     1                                     ,units_3d,comment_3d
     1                                     ,uanl,vanl,wanl
     1                                     ,NX_L,NY_L,NZ_L
     1                                     ,N_3D_FIELDS,istat_lw3)

        real uanl(NX_L,NY_L,NZ_L),vanl(NX_L,NY_L,NZ_L) ! WRT True North
        real wanl(NX_L,NY_L,NZ_L)

        real outarray_4D(NX_L,NY_L,NZ_L,N_3D_FIELDS)
        character*125 comment_3D(N_3D_FIELDS)
        character*10 units_3D(N_3D_FIELDS)
        character*3 var_3D(N_3D_FIELDS)
        character*3 EXT

csms$ignore begin
        write(6,*)' Subroutine put_laps_multi_3d_jacket...'

        call move_3d(uanl,outarray_4D(1,1,1,1),NX_L,NY_L,NZ_L)
        call move_3d(vanl,outarray_4D(1,1,1,2),NX_L,NY_L,NZ_L)
        call move_3d(wanl,outarray_4D(1,1,1,3),NX_L,NY_L,NZ_L)

        call put_laps_multi_3d(i4time_sys,EXT,var_3d,units_3d,
     1     comment_3d,outarray_4D,NX_L,NY_L,NZ_L,N_3D_FIELDS,istat_lw3)

csms$ignore end
        return
        end
