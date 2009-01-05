

        subroutine wind_post_process(i4time_sys
     1                              ,uanl,vanl                            ! I
     1                              ,wanl                                 ! O
     1                              ,NX_L,NY_L,NZ_L                       ! I
     1                              ,N_3D_FIELDS                          ! I
     1                              ,heights_3d                           ! I
     1                              ,uanl_sfcitrp,vanl_sfcitrp            ! I
     1                              ,topo,lat,lon,grid_spacing_m          ! I
     1                              ,r_missing_data,l_grid_north_out      ! I
     1                              ,istat_lw3)

        real uanl(NX_L,NY_L,NZ_L),vanl(NX_L,NY_L,NZ_L) ! WRT True North ! I
        real heights_3d(NX_L,NY_L,NZ_L)                                 ! I
        real wanl(NX_L,NY_L,NZ_L)                                       ! O
        real uanl_sfcitrp(NX_L,NY_L),vanl_sfcitrp(NX_L,NY_L)            ! O

        real lat(NX_L,NY_L)
        real lon(NX_L,NY_L)
        real topo(NX_L,NY_L)

        real rk_terrain(NX_L,NY_L)

        logical l_grid_north_out

csms$ignore begin
        write(6,*)' Subroutine wind_post_process...'

!  **** Generate Interpolated SFC analysis ****

        write(6,*)' Generating interpolated laps surface wind'

        i_sfc_bad = 0

        do j = 1,NY_L
        do i = 1,NX_L

!           Interpolate from three dimensional grid to terrain surface
            zlow = height_to_zcoord2(topo(i,j),heights_3d,NX_L,NY_L,NZ_L
     1                                                  ,i,j,istatus)
            if(istatus .ne. 1)then
                write(6,*)' lapswind_anal: error in height_to_zcoord2'
     1                   ,' in sfc wind interpolation',istatus
                write(6,*)i,j,zlow,topo(i,j),
     1                    (heights_3d(i,j,k),k=1,NZ_L)
                return
            endif

            rk_terrain(i,j) = zlow

            klow = max(zlow,1.)
            khigh = klow + 1
            fraclow = float(khigh) - zlow
            frachigh = 1.0 - fraclow

            if( uanl(i,j,klow)  .eq. r_missing_data
     1     .or. vanl(i,j,klow)  .eq. r_missing_data
     1     .or. uanl(i,j,khigh) .eq. r_missing_data
     1     .or. vanl(i,j,khigh) .eq. r_missing_data        )then

                write(6,3333)i,j
3333            format(' Warning: cannot interpolate to sfc at ',2i3)
                i_sfc_bad = 1
                uanl_sfcitrp(i,j) = r_missing_data
                vanl_sfcitrp(i,j) = r_missing_data

            else
                uanl_sfcitrp(i,j) = uanl(i,j,klow ) * fraclow
     1                            + uanl(i,j,khigh) * frachigh

                vanl_sfcitrp(i,j) = vanl(i,j,klow ) * fraclow
     1                            + vanl(i,j,khigh) * frachigh

            endif

        enddo ! j
        enddo ! i

        I4_elapsed = ishow_timer()

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

        return
        end


        subroutine write_wind_output(i4time_sys,EXT,var_3d
     1                              ,uanl,vanl                            ! I
     1                              ,wanl                                 ! I
     1                              ,uanl_sfcitrp,vanl_sfcitrp            ! I
     1                              ,NX_L,NY_L,NZ_L                       ! I
     1                              ,N_3D_FIELDS                          ! I
     1                              ,r_missing_data                       ! I
     1                              ,istat_lw3)

        use mem_wind, ONLY: num_wind_obs       

!       Stuff for 3D winds
        real uanl(NX_L,NY_L,NZ_L),vanl(NX_L,NY_L,NZ_L) ! WRT True North ! I
        real wanl(NX_L,NY_L,NZ_L)                                       ! I

        character*125 comment_3D(N_3D_FIELDS)
        character*10 units_3D(N_3D_FIELDS)
        character*3 var_3D(N_3D_FIELDS)
        character*3 EXT

!       Stuff for SFC Winds
        real uanl_sfcitrp(NX_L,NY_L),vanl_sfcitrp(NX_L,NY_L)            ! I
        real out_sfc_3D(NX_L,NY_L,2)

        character*125 comment_a(2)
        character*10 units_a(2)
        character*3 var_a(2)

        write(6,*)' Subroutine write_wind_output...'

!       Header information for 3D wind
        EXT = 'lw3'

        units_3D(1) = 'M/S'
        units_3D(2) = 'M/S'
        units_3D(3) = 'PA/S'

        var_3D(1) = 'U3'
        var_3D(2) = 'V3'
        var_3D(3) = 'OM'

        do i = 1,3
            write(comment_3D(i),1)num_wind_obs
 1          format('3DWIND num_wind_obs = ',i9)
        enddo ! i

        I4_elapsed = ishow_timer()

        write(6,*)' Calling write routine for all grids ',ext(1:3)
     1                                  ,i4time_sys

        call check_nan3(uanl,NX_L,NY_L,NZ_L,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' ERROR: Nan Detected in uanl field'
            return
        endif

        call check_nan3(vanl,NX_L,NY_L,NZ_L,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' ERROR: Nan Detected in vanl field'
            return
        endif

        call check_nan3(wanl,NX_L,NY_L,NZ_L,istat_lw3)
        if(istat_lw3 .ne. 1)then
            write(6,*)' ERROR: Nan Detected in wanl field'
            return
        endif

        call put_laps_multi_3d_jacket(i4time_sys,EXT,var_3d
     1                               ,units_3d,comment_3d
     1                               ,uanl,vanl,wanl
     1                               ,NX_L,NY_L,NZ_L,N_3D_FIELDS
     1                               ,istat_lw3)
        if(istat_lw3 .eq. 1)then
            write(6,*)' Success writing out LW3 file'
        else
            write(6,*)' Error writing out LW3 file'
        endif

!       Write out derived winds file (sfc wind)
        call move(uanl_sfcitrp,out_sfc_3D(1,1,1),NX_L,NY_L)
        call move(vanl_sfcitrp,out_sfc_3D(1,1,2),NX_L,NY_L)

        ext = 'lwm'

        var_a(1) = 'SU'
        var_a(2) = 'SV'

        do i = 1,2
            units_a(i) = 'm/s'
            comment_a(i) = 'SFCWIND'
        enddo

        call put_laps_multi_2d(i4time_sys,ext,var_a
     1      ,units_a,comment_a,out_sfc_3d,NX_L,NY_L,2,istat_lwm)

        if(istat_lwm .eq. 1)then
            write(6,*)' Success in writing out lwm file (SU,SV)'
        else
            write(6,*)' Error writing out lwm file (SU,SV)'
        endif

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
