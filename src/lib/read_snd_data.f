
        subroutine read_snd_data(lun,i4time_snd,ext                    ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                   ! I
     1                         ,n_profiles                             ! O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! O
     1                         ,c5_name,i4time_ob_pr,obstype           ! O
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! O
     1                         ,istatus)                               ! O

cdoc    Returns sounding data from the SND file

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_pr_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_di_obs(MAX_PR_LEVELS)                               ! L
        real ob_pr_sp_obs(MAX_PR_LEVELS)                               ! L
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)                         ! O

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob
        character ext*(*)

        istatus = 0

        rcycles_pr = 0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' Error getting r_missing_data'
          GO TO 600
      endif

      call open_lapsprd_file_read(lun,i4time_snd,ext,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' Error opening SND file'
          GO TO 600
      endif

      DO i_pr = n_profiles+1,max_pr

        if(i_pr .le. 200 .or. i_pr .eq. (i_pr/10)*10)then
            iwrite = 1
        else
            iwrite = 0
        endif

        read(lun,511,err=530,end=550)ista,nlevels_in,
     1       lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1       c5_name(i_pr),a9time_ob,obstype(i_pr)
  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

!       obstype(i_pr) = 'RAOB'

        call i4time_fname_lp(a9time_ob,i4time_ob,istatus)
        i4time_ob_pr(i_pr) = i4time_ob

        if(iwrite .eq. 1)write(6,512)i_pr,ista,nlevels_in,lat_pr(i_pr)       
     1             ,lon_pr(i_pr)
     1             ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr),a9time_ob       
     1             ,obstype(i_pr)     
 512    format(/' SND #',i4,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9,1x,a8)       

        n_good_levels = 0

        DO level = 1,nlevels_in

          read(lun,*,err=515)ht_in,pr_in,t_in,td_in,di_in,sp_in ! (sp = m/s)

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
                  if(iwrite .eq. 1)
     1            write(6,*)' Pressure was given, height was derived:'       
     1                     ,pr_in,ht_in
              endif
          endif

 505      if(ht_in .ne. r_missing_data .and.
     1       di_in .ne. r_missing_data .and.
     1       sp_in .ne. r_missing_data          )then           

              n_good_levels = n_good_levels + 1

              if(n_good_levels .gt. MAX_PR_LEVELS)then
                  write(6,*)' ERROR: too many sounding (.snd) levels '       
     1                     ,i_pr,n_good_levels,MAX_PR_LEVELS
                  goto515
              endif

              nlevels_obs_pr(i_pr) = n_good_levels
              ob_pr_ht_obs(i_pr,n_good_levels) = ht_in
              ob_pr_pr_obs(i_pr,n_good_levels) = pr_in
              ob_pr_di_obs(n_good_levels) = di_in
              ob_pr_sp_obs(n_good_levels) = sp_in

              call disp_to_uv(ob_pr_di_obs(n_good_levels),
     1                        ob_pr_sp_obs(n_good_levels),
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

      write(6,*)' ERROR: # of soundings+profilers has reached'
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

      close(lun)

      return
      end



        subroutine read_snd_metadata(lun,i4time_snd,ext                ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                   ! I
     1                         ,n_profiles                             ! O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! O
     1                         ,c5_name,i4time_ob_pr,obstype           ! O
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! O
!    1                         ,ob_pr_u_obs,ob_pr_v_obs                ! O
     1                         ,istatus)                               ! O

cdoc    Returns sounding metadata from the SND file

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_pr_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_di_obs(MAX_PR_LEVELS)                               ! L
        real ob_pr_sp_obs(MAX_PR_LEVELS)                               ! L
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)                         ! O

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob
        character ext*(*)

        istatus = 0

        rcycles_pr = 0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' Error getting r_missing_data'
          GO TO 600
      endif

      call open_lapsprd_file_read(lun,i4time_snd,ext,istatus)
      if(istatus .ne. 1)then
          write(6,*) ' Error opening SND file'
          GO TO 600
      endif

      DO i_pr = n_profiles+1,max_pr

        if(i_pr .le. 200 .or. i_pr .eq. (i_pr/10)*10)then
            iwrite = 1
        else
            iwrite = 0
        endif

        read(lun,511,err=530,end=550)ista,nlevels_in,
     1       lat_pr(i_pr),lon_pr(i_pr),elev_pr(i_pr),
     1       c5_name(i_pr),a9time_ob,obstype(i_pr)
  511   format(i12,i12,f11.4,f15.4,f15.0,1x,a5,3x,a9,1x,a8)

!       obstype(i_pr) = 'RAOB'

        call i4time_fname_lp(a9time_ob,i4time_ob,istatus)
        i4time_ob_pr(i_pr) = i4time_ob

        if(iwrite .eq. 1)write(6,512)i_pr,ista,nlevels_in,lat_pr(i_pr)       
     1             ,lon_pr(i_pr)
     1             ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr),a9time_ob
     1             ,obstype(i_pr)     
 512    format(/' SND #',i4,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9,1x,a8)       

        n_good_levels = 0

        DO level = 1,nlevels_in

          read(lun,*,err=515)ht_in,pr_in,t_in,td_in,di_in,sp_in ! (sp = m/s)

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
                  if(iwrite .eq. 1)
     1            write(6,*)' Pressure was given, height was derived:'       
     1                     ,pr_in,ht_in
              endif
          endif

 505      if(ht_in .ne. r_missing_data .and.
     1       di_in .ne. r_missing_data .and.
     1       sp_in .ne. r_missing_data          )then           

              n_good_levels = n_good_levels + 1

              if(n_good_levels .gt. MAX_PR_LEVELS)then
                  write(6,*)' ERROR: too many sounding (.snd) levels '       
     1                     ,i_pr,n_good_levels,MAX_PR_LEVELS
                  goto515
              endif

              nlevels_obs_pr(i_pr) = n_good_levels
              ob_pr_ht_obs(i_pr,n_good_levels) = ht_in
              ob_pr_pr_obs(i_pr,n_good_levels) = pr_in
              ob_pr_di_obs(n_good_levels) = di_in
              ob_pr_sp_obs(n_good_levels) = sp_in

              call disp_to_uv(ob_pr_di_obs(n_good_levels),
     1                        ob_pr_sp_obs(n_good_levels),
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

      write(6,*)' ERROR: # of soundings+profilers has reached'
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

      close(lun)

      return
      end
