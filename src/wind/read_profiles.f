 
        subroutine read_profiles(i4time,heights_3d,lat_pr,lon_pr,   ! I
     1                   lat,lon,                                   ! I
     1                   MAX_PR,MAX_PR_LEVELS,                      ! I
     1                   ob_pr_ht,                                  ! O
     1                   ob_pr_di, ob_pr_sp,                        ! O
     1                   ob_pr_u , ob_pr_v ,                        ! O
     1                   ob_pr_r , ob_pr_t ,                        ! O
     1                   nlevels_obs_pr,                            ! O
     1                   ob_pr_ht_obs,                              ! O
     1                   ob_pr_di_obs,ob_pr_sp_obs,                 ! O
     1                   ob_pr_u_obs ,ob_pr_v_obs ,                 ! O
     1                   rlat_radar,rlon_radar,rheight_radar,       ! I
     1                   n_vel_grids,                               ! I
     1                   u_maps_inc,v_maps_inc,                     ! I
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

!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)
        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_di_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_sp_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_ht(MAX_PR,kmax)
        real ob_pr_di(MAX_PR,kmax)
        real ob_pr_sp(MAX_PR,kmax)
        real ob_pr_u (MAX_PR,kmax) ! Vertically interpolated Profiler wind
        real ob_pr_v (MAX_PR,kmax) ! Vertically interpolated Profiler wind
        real ob_pr_r (MAX_PR,kmax) ! Vertically interpolated Profiler wind
        real ob_pr_t (MAX_PR,kmax) ! Vertically interpolated Profiler wind

!*****************************************************************************

        real*4 heights_3d(imax,jmax,kmax)

!       These two arrays (not used yet) serve for incrementing the out of
!       date profiler obs according to the model rates of change.
        real*4 u_maps_inc(imax,jmax,kmax)
        real*4 v_maps_inc(imax,jmax,kmax)

        character*13 filename13
        character*255 c_filespec
        character ext*31
        character*5 c5_name
        character*9 a9time_ob

        logical l_use_raob, l_use_cdw, l_use_radial_vel

        r_mspkt = .518

        call get_wind_parms(l_use_raob,l_use_cdw,l_use_radial_vel
     1                     ,weight_bkg_const,istatus)     
        if(istatus .ne. 1)then
            write(6,*)' Error getting wind parms'
            return
        endif

        write(6,*)' Subroutine read_profiles: i4time = ',i4time

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

        i4time_prg = i4time

        ext = 'prg'
        call open_lapsprd_file(32,i4time_prg,ext,istatus)
        if(istatus .ne. 1)go to 880


! ***   Read in profiler data    ***************************************

!       Open nearest PRO file to the LAPS analysis time
        ext = 'pro'
        call get_filespec(ext,2,c_filespec,istatus)
        call get_file_time(c_filespec,i4time,i4time_prof)
        call open_lapsprd_file(12,i4time_prof,ext,istatus)
        if(istatus .ne. 1)go to 420

        do i_pr = 1,max_pr

            read(12,401,err=430,end=450)
     1         ista,nlevels_obs_pr(i_pr),lat_pr(i_pr),lon_pr(i_pr)
     1         ,elev_pr(i_pr),c5_name,a9time_ob
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,3x,a9)

            call i4time_fname_lp(a9time_ob,i4time_ob,istatus)

            rcycles_pr(i_pr) = float(i4time - i4time_ob)       
     1                                      / float(ilaps_cycle_time)


            write(6,407)i_pr,ista,nlevels_obs_pr(i_pr),lat_pr(i_pr)
     1                 ,lon_pr(i_pr)
     1                 ,elev_pr(i_pr),rcycles_pr(i_pr),c5_name,a9time_ob       
407         format(/' Profiler #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a6,3x,a9)

            do level = 1,nlevels_obs_pr(i_pr)

                read(12,*,err=312)ob_pr_ht_obs(i_pr,level)
     1           ,ob_pr_di_obs(i_pr,level)
     1           ,ob_pr_sp_obs(i_pr,level)

                ob_pr_sp_obs(i_pr,level) = ob_pr_sp_obs(i_pr,level)
!    1                                          * r_mspkt

                call disp_to_uv(ob_pr_di_obs(i_pr,level),
     1                  ob_pr_sp_obs(i_pr,level),
     1                  ob_pr_u_obs(i_pr,level),
     1                  ob_pr_v_obs(i_pr,level))

!d              write(6,311,err=312)ista,i_pr
!d      1                ,ob_pr_ht_obs(i_pr,level)
!d      1                ,ob_pr_di_obs(i_pr,level)
!d      1                ,ob_pr_sp_obs(i_pr,level)
!d      1                ,ob_pr_u_obs(i_pr,level)
!d      1                ,ob_pr_v_obs(i_pr,level)
311             format(1x,2i4,5f8.1)
312             continue
            enddo ! level
        enddo ! profiler site

        write(6,*)' WARNING: # of profilers has reached'
     1           ,' dimensions of MAX_PR ',MAX_PR

        n_profiles=MAX_PR
        GO TO 500
c
c       Profiler reading error handling
c
  420   write(6,*) ' Error opening PRO file'
        n_profiles=0
        GO TO 500

  430   write(6,*) ' Error during read of PRO file'
        write(6,*) ' While reading profiler number ',i_pr
        n_profiles=i_pr-1
        GO TO 500
c
  450   CONTINUE ! Used for end of file
        n_profiles=i_pr-1

        close(12)

  500   CONTINUE

        n_profilers = n_profiles
c
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
      call get_file_time(c_filespec,i4time,i4time_nearest)

      i4time_diff = abs(i4time - i4time_nearest)
      if(i4time_diff .le. i4time_raob_file_window)then
          write(6,*)' Nearest SND file is within time window'
     1                ,i4time_diff,i4time_raob_file_window
      else
          write(6,*)' Warning: Nearest SND file is outside time window'       
     1                ,i4time_diff,i4time_raob_file_window
          go to 600
      endif

      i4time_snd = i4time_nearest

      call open_lapsprd_file(12,i4time_snd,ext,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' Error opening SND file'
          GO TO 600
      endif

      DO i_pr = n_profiles+1,max_pr
        read(12,511,err=530,end=550)
     1  ista,nlevels_in,
     1  lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),c5_name,a9time_ob
  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9)

        call i4time_fname_lp(a9time_ob,i4time_ob,istatus)
        rcycsnd = float(i4time - i4time_ob) / float(ilaps_cycle_time)

        rcycles_pr(i_pr) = max(min(rcycsnd,1.0),-1.0)

        write(6,512)i_pr,ista,nlevels_in,lat_pr(i_pr)
     1             ,lon_pr(i_pr)
     1             ,elev_pr(i_pr),rcycles_pr(i_pr),c5_name,a9time_ob     
 512    format(/' RAOB #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9)       


        n_good_levels = 0

        DO level = 1,nlevels_in

          read(12,*,err=515)ht_in,pr_in,t_in,td_in,di_in,sp_in ! (sp = m/s)

!         Test this by deliberately setting ht_in to missing
!         ht_in = r_missing_data

!         Determine whether we need to supply our own height (only pres given)
          if(ht_in .eq. r_missing_data .and. 
     1       pr_in .ne. r_missing_data                      )then

              call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon       
     1                              ,imax,jmax,ri,rj,istatus)

              if(istatus .ne. 1)goto505

              i_ob = nint(ri)
              j_ob = nint(rj)
              if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1           j_ob .ge. 1 .and. j_ob .le. jmax      )then
                  pr_in_pa = pr_in * 100.
                  call pressure_to_height(pr_in_pa,heights_3d
     1                     ,imax,jmax,kmax,i_ob,j_ob,ht_buff,istatus)
                  if(istatus .ne. 1)goto505
                  ht_in = ht_buff
                  write(6,*)' Pressure was given, height was derived:'
     1                     ,pr_in,ht_in
              endif
          endif

 505      if(ht_in .ne. r_missing_data .and.
     1       di_in .ne. r_missing_data .and.
     1       sp_in .ne. r_missing_data          )then           

              n_good_levels = n_good_levels + 1
              nlevels_obs_pr(i_pr) = n_good_levels
              ob_pr_ht_obs(i_pr,n_good_levels) = ht_in
              ob_pr_di_obs(i_pr,n_good_levels) = di_in
              ob_pr_sp_obs(i_pr,n_good_levels) = sp_in

              call disp_to_uv(ob_pr_di_obs(i_pr,n_good_levels),
     1                        ob_pr_sp_obs(i_pr,n_good_levels),
     1                        ob_pr_u_obs(i_pr,n_good_levels),
     1                        ob_pr_v_obs(i_pr,n_good_levels))

!             write(6,311,err=312)ista,i_pr
!    1        ,ob_snd_ht_obs(i_pr,level)
!    1        ,ob_snd_t_obs(i_pr,level)
!311          format(1x,i6,i4,5f8.1)

          endif ! Good data at the level

          go to 516

  515     write(6,*)' Error reading RAOB, raw level # =',level
          write(6,*)' While reading sounding # ',(i_pr-n_profiles)
          n_profiles=i_pr-1
          go to 600

  516   END DO                   ! level
      END DO ! i_pr

      write(6,*)' WARNING: # of soundings+profilers has reached'
     1         ,' dimensions of MAX_PR ',MAX_PR

      n_profiles=MAX_PR
      GO TO 600
c
c     Sounding reading error handling
c
  530 write(6,*) ' Error during read of SND file'
      write(6,*) ' While reading sounding number ',(i_pr-n_profiles)
      n_profiles=i_pr-1
      GO TO 600

  550 CONTINUE ! Used for end of file
      n_profiles=i_pr-1

  600 CONTINUE 

      close(12)

      n_snd=n_profiles-n_profilers

      write(6,*)
      write(6,*) ' Read ',n_profilers,' wind profiler(s).'
      write(6,*) ' Read ',n_snd,' sounding(s).'
      write(6,*)
c
c     Process all wind profiles.  Interpolate heights to LAPS levels.
c
      DO i_pr=1,n_profiles

!           Determine if profile is in the LAPS domain

            call latlon_to_rlapsgrid(lat_pr(i_pr),lon_pr(i_pr),lat,lon
     1                              ,imax,jmax,ri,rj,istatus)

            i_ob = nint(ri)
            j_ob = nint(rj)
            if(i_ob .ge. 1 .and. i_ob .le. imax .and.
     1         j_ob .ge. 1 .and. j_ob .le. jmax      )then
                write(6,*)'Profile  # ',i_pr,' In Bounds - Doing '
     1                   ,'Vertical Interpolation'
            else
                write(6,*)'Profile  # ',i_pr,' Out of Domain Bounds'       
                nlevels_obs_pr(i_pr)=0 ! This effectively throws out the profile
            endif

!           Determine if profile was obtained close enough in time....
            if(abs(rcycles_pr(i_pr)) .gt. 1.0)then
                write(6,*)'Profile  # ',i_pr,' Out of time bounds:'
     1                                      ,rcycles_pr(i_pr)
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

                    u_diff = u_maps_inc(i_ob,j_ob,level) 
     1                                              * rcycles_pr(i_pr)       
                    v_diff = v_maps_inc(i_ob,j_ob,level) 
     1                                              * rcycles_pr(i_pr)

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
