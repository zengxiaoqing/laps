  

        subroutine get_wind_3d_obs(
     1            NX_L,NY_L,NZ_L,                                   ! I
     1            r_missing_data,i2_missing_data,                   ! I
     1            i4time_lapswind,heights_3d,heights_1d,            ! I
     1            lat_pr,lon_pr,                                    ! O
     1            lat,lon,                                          ! I
     1            u_mdl_diff,v_mdl_diff,                            ! I
     1            grid_laps_u,grid_laps_v,grid_laps_wt,             ! O
     1            ob_pr_ht,                                         ! O
     1            ob_pr_di, ob_pr_sp,                               ! O
     1            ob_pr_u , ob_pr_v ,                               ! O
     1            ob_pr_r , ob_pr_t ,                               ! O
     1            nlevels_obs_pr,                                   ! O
     1            rlat_radar,rlon_radar,rheight_radar,              ! I
     1            istat_radar_vel,n_vel_grids,                      ! I
     1            istatus_remap_pro,                                ! O
     1            istatus                )                          ! O

!       1997 Jun     Ken Dritz       Added NX_L, NY_L, NZ_L as dummy arguments,
!                                    making arrays with those dimensions
!                                    automatic.
!       1997 Jun     Ken Dritz       Added r_missing_data and i2_missing_data
!                                    as dummy arguments.
!       1997 Jun     Ken Dritz       Removed include of 'lapsparms.for'.

        include 'windparms.inc'
        include 'get_wind_obs.inc'

        dimension u_mdl_diff(NX_L,NY_L,NZ_L),v_mdl_diff(NX_L,NY_L,NZ_L)
        real*4 heights_3d(NX_L,NY_L,NZ_L)
        real*4 heights_1d(NZ_L)

!  ***  Read in Profiler Data  ********************************************

        call get_laps_cycle_time(ilaps_cycle_time,istatus)
        if(istatus .eq. 1)then
            write(6,*)' ilaps_cycle_time = ',ilaps_cycle_time
        else
            write(6,*)' Error getting laps_cycle_time'
            return
        endif

        write(6,*)' Calling read_profiles'

        call read_profiles(
     1            i4time_lapswind,heights_3d,lat_pr,lon_pr,         ! I
     1            lat,lon,                                          ! I
     1            ob_pr_ht,
     1            ob_pr_di, ob_pr_sp,                               ! O
     1            ob_pr_u , ob_pr_v ,                               ! O
     1            ob_pr_r , ob_pr_t ,                               ! O
     1            nlevels_obs_pr,                                   ! O
     1            ob_pr_ht_obs,                                     ! O
     1            ob_pr_di_obs, ob_pr_sp_obs,                       ! O
     1            ob_pr_u_obs , ob_pr_v_obs ,                       ! O
     1            rlat_radar,rlon_radar,rheight_radar,              ! I
     1            n_vel_grids,                                      ! I
     1            u_mdl_diff,v_mdl_diff,                            ! I
     1            ilaps_cycle_time,r_missing_data,                  ! I
     1            NX_L,NY_L,NZ_L,                                   ! I
     1            istatus                )                          ! O

        if(istatus .ne. 1)then
            write(6,*)' Abort read_profiles'
            return
        endif

        I4_elapsed = ishow_timer()

! ***   Remapping + Barnes Analysis of Profiler Data in u & v ******************************

        call remap_profiler(
     1           ob_pr_u,ob_pr_v                                  ! I
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt             ! O
     1          ,lat,lon                                          ! I
     1          ,NX_L,NY_L,NZ_L,MAX_PR                            ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr                     ! I
     1          ,r_missing_data                                   ! I
     1          ,weight_prof                                      ! O
     1          ,l_profiler                                       ! I
     1          ,istatus_remap_pro)                               ! O



!       It would be interesting to see if this routine affects the LW3 file in
!       the presence of radar velocity data. If it doesn't as expected, then it
!       would simplify things to eliminate this call.
        call prof_anal(istat_radar_vel                             ! I
     1          ,n_vel_grids                                       ! I
     1          ,sum_pr_u,sum_pr_v,sum_wt_pr,max_weights_pr        ! I
     1          ,weights_pr                                        ! L
     1          ,ob_pr_u,ob_pr_v                                   ! I
     1          ,grid_pr_u,grid_pr_v,grid_pr_r                     ! L
     1          ,grid_ra_vel,lat,lon,rlat_radar,rlon_radar,rheight_radar ! I
     1          ,i2_missing_data                                   ! I
     1          ,i_profiler_nearest                                ! L
     1          ,NX_L,NY_L,NZ_L,MAX_PR                             ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr                      ! I
     1          ,r_missing_data,weight_prof                        ! I
     1          ,heights_1d                                        ! I
     1          ,istatus)                                          ! O

        if(istatus .ne. 1)then
            write(6,*)' Error in prof_anal'
            return
        endif


        return
        end


        subroutine remap_profiler(
     1           ob_pr_u,ob_pr_v                                     ! I
     1          ,grid_laps_u,grid_laps_v,grid_laps_wt                ! O
     1          ,lat,lon                                             ! I
     1          ,ni,nj,nk,MAX_PR                                     ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr                        ! I
     1          ,r_missing_data                                      ! I
     1          ,weight_prof                                         ! O
     1          ,l_profiler                                          ! I
     1          ,istatus)                                            ! O

!          Perform horizontal remapping of profile obs onto LAPS grid
!          They have already been vertically remapped

!       Profile Observations

        integer nlevels_obs_pr(MAX_PR)
        real*4 lat_pr(MAX_PR)
        real*4 lon_pr(MAX_PR)
        real*4 ob_pr_u (MAX_PR,nk) ! Vertically interpolated Profile wind
        real*4 ob_pr_v (MAX_PR,nk) ! Vertically interpolated Profile wind

!       Barnes Profile analysis

        real*4 lat(ni,nj),lon(ni,nj)

        real*4 grid_laps_u(ni,nj,nk)
        real*4 grid_laps_v(ni,nj,nk)
        real*4 grid_laps_wt(ni,nj,nk)

        logical l_profiler

        write(6,*)
        write(6,*)' Subroutine remap profiler'

        do i_pr = 1,MAX_PR
            if(nlevels_obs_pr(i_pr) .gt. 0)then
                call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr)
     1                                  ,lat,lon,ni,nj,ri,rj,istatus)       
                if(istatus .ne. 1)then
                    write(6,*)
     1                ' ERROR... Profile apparently outside domain'
                    return
                endif

                i_ob = nint(ri)
                j_ob = nint(rj)

                write(6,*)' Remapping profile',i_pr,i_ob,j_ob

                do k = 1,nk

                     if(ob_pr_u(i_pr,k) .ne. r_missing_data)then

                         ob_u = ob_pr_u (i_pr,k)
                         ob_v = ob_pr_v (i_pr,k)

!                 ***    Map Observation onto LAPS grid   ***
                         if(l_profiler)then
                             grid_laps_u(i_ob,j_ob,k) = ob_u
                             grid_laps_v(i_ob,j_ob,k) = ob_v
                             grid_laps_wt(i_ob,j_ob,k) = weight_prof
                             write(6,11)k,ob_u,ob_v
 11                          format(10x,i4,2f8.1)
                         endif

                     endif ! In bounds vertically (of profile data)
                enddo ! level

            endif ! data present
        enddo ! i_pr
        I4_elapsed = ishow_timer()

        istatus = 1

        return
        end



        subroutine prof_anal(istat_radar_vel                            ! I
     1          ,n_vel_grids                                            ! I
     1          ,sum_pr_u,sum_pr_v,sum_wt_pr,max_weights_pr             ! I
     1          ,weights_pr                                             ! O
     1          ,ob_pr_u,ob_pr_v                                        ! I
     1          ,grid_pr_u,grid_pr_v,grid_pr_r ! ,grid_pr_t             ! L
     1          ,grid_ra_vel,lat,lon,rlat_radar,rlon_radar,rheight_radar! I
     1          ,i2_missing_data,i_profiler_nearest                     ! I/O
     1          ,ni,nj,nk,MAX_PR                                        ! I
     1          ,nlevels_obs_pr,lat_pr,lon_pr                           ! I
     1          ,r_missing_data,weight_prof                             ! I/O
     1          ,heights_1d                                             ! I
     1          ,istatus)                                               ! I

!          Do a rough analysis of profile data for comparison with radial
!          velocities

!          It is questionable whether this routine actually passes in
!          the nyquist velocity - or does anything else useful for that matter

!       Profile Observations

        integer nlevels_obs_pr(MAX_PR)
        real*4 lat_pr(MAX_PR)
        real*4 lon_pr(MAX_PR)
        real*4 ob_pr_u (MAX_PR,nk) ! Vertically interpolated Profile wind
        real*4 ob_pr_v (MAX_PR,nk) ! Vertically interpolated Profile wind

!       Barnes Profile analysis

        real*4 weights_pr(ni,nj)
        real*4 lat(ni,nj),lon(ni,nj)
        real*4 max_weights_pr(ni,nj,nk)
        integer*2 i_profiler_nearest(ni,nj,nk)

        real*4 sum_pr_u(ni,nj,nk)
        real*4 sum_pr_v(ni,nj,nk)
        real*4 sum_wt_pr(ni,nj,nk)

        real*4 grid_pr_u(ni,nj,nk) ! Barnes anal. Profile (u)
        real*4 grid_pr_v(ni,nj,nk) ! Barnes anal. Profile (v)
!       real*4 grid_pr_t(ni,nj,nk) ! Barnes anal. Profile (u)
        real*4 grid_pr_r(ni,nj,nk) ! Barnes anal. Profile (v)

        real*4 grid_ra_vel(ni,nj,nk)

        real*4 heights_1d(nk)

        write(6,*)
        write(6,*)' Subroutine prof_anal'

        if(istat_radar_vel .eq. 1 .and. n_vel_grids .gt. 0)then
            write(6,*)' Initializing profile weight arrays'
            do k = 1,nk
                do j = 1,nj
                do i = 1,ni
                    sum_pr_u(i,j,k) = 0.
                    sum_pr_v(i,j,k) = 0.
                    sum_wt_pr(i,j,k) = 0.
                    max_weights_pr(i,j,k) = 0.
                    i_profiler_nearest(i,j,k) = 0
                enddo ! i
                enddo ! j
            enddo ! k
        endif ! Valid radar data (then initialize profile weight arrays)

        I4_elapsed = ishow_timer()

        do i_pr = 1,MAX_PR
            if(nlevels_obs_pr(i_pr) .gt. 0)then
                call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr)
     1                                  ,lat,lon,ni,nj,ri,rj,istatus)
                if(istatus .ne. 1)return

                i_ob = nint(ri)
                j_ob = nint(rj)

                write(6,*)' Analyzing profile',i_pr,i_ob,j_ob

                if(istat_radar_vel .eq. 1 .and. n_vel_grids .gt. 0)the
     1n
                    write(6,*)' Determining weights for profile',i_pr
                    call weights_hor(i_ob,j_ob,weights_pr,i_pr,ni,nj)
                endif

                do k = 1,nk

                     if(ob_pr_u(i_pr,k) .ne. r_missing_data)then

                         ob_u = ob_pr_u (i_pr,k)
                         ob_v = ob_pr_v (i_pr,k)

                         if(istat_radar_vel .eq. 1 .and. n_vel_grids
     1                                                      .gt. 0)then

!d                         write(6,201)i_pr,k,ob_pr_ht(i_pr,k),
!d      1                    ob_pr_di(i_pr,k),ob_pr_sp(i_pr,k)
!d201                      format(1x,' Summing Wts',2i4,2f8.0,f8.1)

                           do j = 1,nj
                           do i = 1,ni
                             sum_pr_u(i,j,k) = sum_pr_u(i,j,k) +
     1                  weights_pr(i,j) * ob_u
                             sum_pr_v(i,j,k) = sum_pr_v(i,j,k) +
     1                  weights_pr(i,j) * ob_v
                             sum_wt_pr(i,j,k) = sum_wt_pr(i,j,k) +
     1                  weights_pr(i,j)

                             if(weights_pr(i,j) .gt. max_weights_pr(i,j,
     1k))then
                                 max_weights_pr(i,j,k) = weights_pr(i,j)
                                 i_profiler_nearest(i,j,k) = i_pr
                             endif

                           enddo ! i
                           enddo ! j

                         endif ! Valid radar data present
                     endif ! In bounds vertically (of profile data)
                enddo ! level

            endif ! data present
        enddo ! i_pr
        I4_elapsed = ishow_timer()

        if(istat_radar_vel .eq. 1 .and. n_vel_grids .gt. 0)then
          write(6,*)' Dividing weights'
          do k = 1,nk
            do j = 1,nj
                 do i = 1,ni

                      if(sum_wt_pr(i,j,k) .gt. 0.)then
                          grid_pr_u(i,j,k)
     1                  = sum_pr_u(i,j,k) / sum_wt_pr(i,j,k)
                          grid_pr_v(i,j,k)
     1                  = sum_pr_v(i,j,k) / sum_wt_pr(i,j,k)
                      else
                          grid_pr_u(i,j,k) = r_missing_data
                          grid_pr_v(i,j,k) = r_missing_data
                          i_profiler_nearest(i,j,k) = i2_missing_data
                      endif

                enddo ! i
            enddo ! j
          enddo ! k
        endif ! Valid velocity data

        I4_elapsed = ishow_timer()

! ***   Convert profile to radial and tangential radar coordinates   ********************************
        write(6,*)' n_vel_grids/istat = ',n_vel_grids,istat_radar_vel       

        if(istat_radar_vel .eq. 1 .and. n_vel_grids .gt. 0)then
          write(6,*)' Generating radar coordinates'
!         do l = 1,n_radars
          do k = 1,nk
            do j = 1,nj
            do i = 1,ni

              if(grid_ra_vel(i,j,k) .ne. r_missing_data)then ! for efficiency

                call latlon_to_radar(lat(i,j),lon(i,j),height_of_level(k
     1)
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                if(grid_pr_u(i,j,k) .ne. r_missing_data)then

                    call uvtrue_to_radar(grid_pr_u(i,j,k),
     1                           grid_pr_v(i,j,k),
     1                           dum_t, !            grid_pr_t(i,j,k),
     1                           grid_pr_r(i,j,k),
     1                           azimuth,
     1                           lon(i,j))

                else
                    dum_t = r_missing_data
                    grid_pr_r(i,j,k) = r_missing_data


                endif


              endif ! efficiency check

            enddo ! j
            enddo ! i
          enddo ! k
!         enddo ! l

          I4_elapsed = ishow_timer()

!         Compare Remapped Radial Velocities to Profile Obs
!         do l = 1,n_radars
              write(6,*)' Radar # ',1
              call comp_vr_prof(grid_ra_vel
     1                 ,grid_pr_r,r_missing_data,ni,nj,nk
     1                 ,lat,lon,max_weights_pr,i_profiler_nearest
     1                       ,v_nyquist_2,unfolding_thresh
     1                 ,heights_1d,rlat_radar,rlon_radar,rheight_radar)
!         enddo ! l

        endif ! We have radar data
        I4_elapsed = ishow_timer()


        return
        end

