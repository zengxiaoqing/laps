 
        subroutine read_profiles(i4time_sys,heights_3d,             ! I
     1                   lat_pr,lon_pr,                             ! O
     1                   lat,lon,                                   ! I
     1                   MAX_PR,MAX_PR_LEVELS,                      ! I
     1                   l_use_raob,                                ! I
     1                   ob_pr_u , ob_pr_v ,                        ! O
     1                   nlevels_obs_pr,                            ! O
     1                   rlat_radar,rlon_radar,rheight_radar,       ! I
     1                   n_vel_grids,                               ! I
     1                   u_mdl_bkg_4d,v_mdl_bkg_4d,NTMIN,NTMAX,     ! I
     1                   ilaps_cycle_time,r_missing_data,           ! I
     1                   imax,jmax,kmax,                            ! I
     1                   istatus                )                   ! O

!       1992 Steve Albers
!       Note that the profiler data in the .PRO files are in knots...
!       1994 Keith Brewster   Added reading of sounding data
c       1995 Keith Brewster   Re-Added reading of sounding data, improved
c                             error handling
c       1996 Steve Albers     Added read of ob times from PRO files
c       1996 Steve Albers     Read nearest PRO file, even if its time does
c                             not exactly match the LAPS analysis time. 
c                             Accept only those profiler obs whose 
c                             mid-window times are within one LAPS cycle
c                             time of the current LAPS analysis time.


!*****************************************************************************

!       LAPS Grid Dimensions

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)
        real rcycles_pr(MAX_PR)
        integer i4time_ob_pr(MAX_PR)

!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)
        character*8 obstype(MAX_PR)
        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)                             ! L
        real ob_pr_pr_obs(MAX_PR,MAX_PR_LEVELS)                             ! L
!       real ob_pr_di_obs(MAX_PR_LEVELS)                                    ! L
!       real ob_pr_sp_obs(MAX_PR_LEVELS)                                    ! L
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)                              ! L
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)                              ! L
        real ob_pr_ht(MAX_PR,kmax)                                          ! L
        real ob_pr_di(MAX_PR,kmax)                                          ! L
        real ob_pr_sp(MAX_PR,kmax)                                          ! L
        real ob_pr_u (MAX_PR,kmax) ! Vertically interpolated Profiler wind  ! O
        real ob_pr_v (MAX_PR,kmax) ! Vertically interpolated Profiler wind  ! O
        real ob_pr_r (MAX_PR,kmax) ! Vertically interpolated Profiler wind  ! L
        real ob_pr_t (MAX_PR,kmax) ! Vertically interpolated Profiler wind  ! L
        real sfc_t(MAX_PR), sfc_p(MAX_PR), sfc_rh(MAX_PR)                   ! L
        real sfc_u(MAX_PR), sfc_v(MAX_PR)                                   ! L

!*****************************************************************************

        real*4 heights_3d(imax,jmax,kmax)

        dimension u_mdl_bkg_4d(imax,jmax,kmax,NTMIN:NTMAX)
        dimension v_mdl_bkg_4d(imax,jmax,kmax,NTMIN:NTMAX)

        character*255 c_filespec
        character ext*31
        character*5 c5_name, c5_name_a(MAX_PR)
        character*9 a9time_ob

        logical l_use_raob

        r_mspkt = .518

        write(6,*)' Subroutine read_profiles: i4time = ',i4time_sys

!       Initialize

        do i_pr = 1,MAX_PR
            nlevels_obs_pr(i_pr) = 0
        enddo

        do i_pr = 1,MAX_PR
            do level = 1,kmax

                ob_pr_ht(i_pr,level) = r_missing_data
                ob_pr_di(i_pr,level) = r_missing_data
                ob_pr_sp(i_pr,level) = r_missing_data
                ob_pr_u(i_pr,level)  = r_missing_data
                ob_pr_v(i_pr,level)  = r_missing_data
                ob_pr_r(i_pr,level)  = r_missing_data
                ob_pr_t(i_pr,level)  = r_missing_data

            enddo
        enddo

        i4time_prg = i4time_sys

        ext = 'prg'
        call open_lapsprd_file(32,i4time_prg,ext,istatus)
        if(istatus .ne. 1)go to 880


! ***   Read in profiler data    ***************************************

!       Open nearest PRO file to the LAPS analysis time
        ext = 'pro'
        call get_filespec(ext,2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time_sys,i4time_prof)

        lun = 12
        call read_pro_data(lun,i4time_prof,ext                        ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                  ! I
     1                         ,n_profiles                            ! O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr  ! O
     1                         ,c5_name_a,i4time_ob_pr,obstype        ! O
     1                         ,ob_pr_ht_obs                          ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs               ! O
     1                         ,sfc_t,sfc_p,sfc_rh,sfc_u,sfc_v        ! O
     1                         ,istatus)                              ! O

        n_profilers = n_profiles

c ***   Read in sonde data    ***************************************
c

      write(6,*)

      if(.not. l_use_raob)then
          write(6,*)' Not using raobs, l_use_raob = ',l_use_raob
          go to 600
      endif

      i4time_raob_file_window = 0

      ext = 'snd'
      call get_filespec(ext,2,c_filespec,istatus)
      call get_file_time(c_filespec,i4time_sys,i4time_nearest)

      i4time_diff = abs(i4time_sys - i4time_nearest)
      if(i4time_diff .le. i4time_raob_file_window)then
          write(6,*)' Nearest SND file is within time window'
     1                ,i4time_diff,i4time_raob_file_window
      else
          write(6,*)' Warning: Nearest SND file is outside time window'       
     1                ,i4time_diff,i4time_raob_file_window
          go to 600
      endif

      i4time_snd = i4time_nearest

      lun = 12
      call read_snd_data(lun,i4time_snd,ext                             ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                    ! I
     1                         ,n_profiles                              ! O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr    ! O
     1                         ,c5_name_a,i4time_ob_pr,obstype          ! O
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs               ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs                 ! O
     1                         ,istatus)                                ! O

 600  continue

      n_snd=n_profiles-n_profilers

      write(6,*)
      write(6,*) ' Read ',n_profilers,' wind profiler(s).'
      write(6,*) ' Read ',n_snd,' sounding(s).'
      write(6,*)
c
c     Process all wind profiles.  Interpolate heights to LAPS levels.
c
      DO i_pr=1,n_profiles

            rcycles_pr(i_pr) = float(i4time_sys - i4time_ob_pr(i_pr))       
     1                                      / float(ilaps_cycle_time)

            if(i_pr .le. 200 .or. i_pr .eq. (i_pr/10)*10)then
                iwrite = 1
            else
                iwrite = 0
            endif

!           Determine if profile is in the LAPS domain

            call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)
            if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1         j_ob .ge. 1 .and. j_ob .le. jmax      )then
                if(iwrite .eq. 1)
     1             write(6,*)'Profile  # ',i_pr,' In Bounds - Doing '       
     1                   ,'Vertical Interpolation'
            else
                if(iwrite .eq. 1)
     1             write(6,*)'Profile  # ',i_pr,' Out of Domain Bounds'       
                nlevels_obs_pr(i_pr)=0 ! This effectively throws out the profile
            endif

            call get_windob_time_window(obstype(i_pr),i4_window_ob
     1                                               ,istatus)

            rcyc_thresh = float(i4_window_ob)
     1                   /float(ilaps_cycle_time)

            rcyc_thresh = min(1.0,rcyc_thresh)

!           Determine if profile was obtained close enough in time....
            if(abs(rcycles_pr(i_pr)) .gt. rcyc_thresh)then
                if(iwrite .eq. 1)
     1             write(6,*)'Profile  # ',i_pr,' Out of time bounds:'       
     1                                    ,rcycles_pr(i_pr)
                nlevels_obs_pr(i_pr)=0 ! This effectively throws out the profile
            endif

!  ***  Interpolate the profiles to the LAPS grid levels  *******

            if(nlevels_obs_pr(i_pr) .gt. 0)then

                if(n_vel_grids .gt. 0)then ! Calculate azimuth relative to radar
                   ht = 0.
                   call latlon_to_radar(lat_pr(i_pr),
     1                          lon_pr(i_pr),
     1                          ht,
     1                          azimuth,
     1                          slant_range,
     1                          elevation_angle,
     1                          rlat_radar,rlon_radar,rheight_radar)       
                endif


                do level = 1,kmax

                    ht = heights_3d(i_ob,j_ob,level)

                    ob_pr_ht(i_pr,level) = ht

                    call get_time_term(u_mdl_bkg_4d,imax,jmax,kmax
     1                                ,NTMIN,NTMAX
     1                                ,i_ob,j_ob,level
     1                                ,i4time_sys,i4time_ob_pr(i_pr)
     1                                ,u_time_interp,u_diff_term
     1                                ,istatus)
                    u_diff = -u_diff_term

                    call get_time_term(v_mdl_bkg_4d,imax,jmax,kmax
     1                                ,NTMIN,NTMAX
     1                                ,i_ob,j_ob,level
     1                                ,i4time_sys,i4time_ob_pr(i_pr)
     1                                ,v_time_interp,v_diff_term
     1                                ,istatus)
                    v_diff = -v_diff_term

                    call interp_prof(ob_pr_ht_obs,ob_pr_u_obs,
     1                             ob_pr_v_obs,
     1                             u_diff,
     1                             v_diff,
     1                             ob_pr_u(i_pr,level),
     1                             ob_pr_v(i_pr,level),
     1                             ob_pr_r(i_pr,level),
     1                             ob_pr_t(i_pr,level),
     1                             ob_pr_di(i_pr,level),
     1                             ob_pr_sp(i_pr,level),
     1                             i_pr,ht,level,nlevels_obs_pr,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             azimuth,r_missing_data,
     1                             heights_3d,imax,jmax,kmax,
     1                             MAX_PR,MAX_PR_LEVELS,
     1                             n_vel_grids,istatus)

c                   write(6,411,err=412)ista,i_pr,level
c       1                ,ob_pr_ht(i_pr,level)
c       1                ,ob_pr_di(i_pr,level)
c       1                ,ob_pr_sp(i_pr,level)
c       1                ,ob_pr_u(i_pr,level)
c       1                ,ob_pr_v(i_pr,level)
c       1                ,u_diff
c       1                ,v_diff
c       1                ,ob_pr_r(i_pr,level)
c       1                ,ob_pr_t(i_pr,level)
411                 format(1x,i6,2i4,f8.1,8f7.1)

412                 write(32,*)ri-1.,rj-1.,level-1
     1         ,ob_pr_di(i_pr,level),ob_pr_sp(i_pr,level)

                enddo ! level

            endif ! # levels > 0

        enddo  ! i_pr

        close(32)
        istatus=1
        return

  880   CONTINUE
        write(6,*) ' Error opening PRG file'
        istatus=0
        return
        end

