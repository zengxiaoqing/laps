
        subroutine read_snd_data(lun,i4time_snd,ext                    ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                   ! I
     1                         ,lat,lon,imax,jmax,kmax                 ! I
     1                         ,heights_3d                             ! I
     1                         ,mode                                   ! I
     1                         ,n_profiles                             ! I/O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! O
     1                         ,c5_name,i4time_ob_pr,obstype           ! O
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! O
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! O
     1                         ,istatus)                               ! O

cdoc    Returns sounding wind data from the SND file

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
        real ob_pr_t_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_td_obs(MAX_PR,MAX_PR_LEVELS)                        ! O

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob
        character ext*(*)

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 heights_3d(imax,jmax,kmax)

        logical l_good_level

        IPR_TIMES_LVLS = MAX_PR*MAX_PR_LEVELS

        write(6,*)' Subroutine read_snd_data: ',MAX_PR,MAX_PR_LEVELS
     1                                         ,IPR_TIMES_LVLS

        if(IPR_TIMES_LVLS .gt. 1000000)then
            write(6,*)' Warning, large allocation for IPR_TIMES_LVLS' 
     1               ,IPR_TIMES_LVLS
        endif

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
 512    format(1x,' SND #',i4,i6,i5,2f8.2,e10.3,f8.2,1x,a5,3x,a9,1x,a8)       

        n_good_levels = 0

        ht_prev = -9999.

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

 505      if(mode .eq. 1 .or. mode .eq. 2)then   ! key good levels to wind data
              if(ht_in .ne. r_missing_data .and.
     1           di_in .ne. r_missing_data .and.
     1           sp_in .ne. r_missing_data          )then           
                  l_good_level = .true.
              else
                  l_good_level = .false.
              endif
          endif ! mode

          if(l_good_level)then           

              n_good_levels = n_good_levels + 1

              if(n_good_levels .gt. MAX_PR_LEVELS)then
                  write(6,*)' ERROR: too many sounding (.snd) levels '       
     1                     ,i_pr,n_good_levels,MAX_PR_LEVELS
                  goto515
              endif

              if(ht_in .le. ht_prev)then
                  write(6,*)
     1              ' ERROR: sounding (.snd) ht levels out of sequence '      
     1                     ,i_pr,n_good_levels,ht_prev,ht_in
                  goto515
              endif

              ht_prev = ht_in

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

              ob_pr_t_obs(i_pr,n_good_levels)  = t_in
              ob_pr_td_obs(i_pr,n_good_levels) = td_in

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
     1                         ,lat,lon,imax,jmax,kmax                 ! I
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

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)

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

	subroutine read_sfc_snd(i4time,atime_s,n_obs_g,n_obs_b, ! regular SND
     &       obstime,wmoid,stations,provider,wx_s,reptype,autostntype,
     &       lat_s,lon_s,elev_s,t_s,td_s,rh_s,dd_s,ff_s,ddg_s,ffg_s,
     &       alt_s,pstn_s,pmsl_s,delpch,delp,vis_s,solar_s,sfct,sfcm,
     &       pcp1,pcp3,pcp6,pcp24,snow,kloud_s,max24t,min24t,t_ea,
     &       td_ea,rh_ea,dd_ea,ff_ea,alt_ea,p_ea,vis_ea,solar_ea,
     &       sfct_ea,sfcm_ea,pcp_ea,snow_ea,store_amt,store_hgt,maxsta,       
     &       lat,lon,imax,jmax,kmax,                                     ! I
     &       istatus)

        include 'constants.inc'

        include 'read_sfc.inc'

!       Declarations for 'read_snd_data' call
        integer MAX_PR,MAX_PR_LEVELS
        parameter (MAX_PR = 1500)
        parameter (MAX_PR_LEVELS = 300)

        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_pr_obs(MAX_PR,MAX_PR_LEVELS)                        ! O
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_t_obs(MAX_PR,MAX_PR_LEVELS)                         ! O
        real ob_pr_td_obs(MAX_PR,MAX_PR_LEVELS)                        ! O

        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob

!       Other declarations
        logical l_good_snd
        character*3 ext

        write(6,*)' Subroutine read_sfc_snd...'

        call get_r_missing_data(r_missing_data,istatus)
        call get_sfc_badflag(badflag,istatus)

        n_good_snd = 0
        ext = 'snd'
        lun = 44
        mode = 2

        call read_snd_data(lun,i4time,ext                              ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                   ! I
     1                         ,lat,lon,imax,jmax,kmax                 ! I
     1                         ,heights_3d                             ! I
     1                         ,mode                                   ! I
     1                         ,n_profiles                             ! I/O
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr   ! O
     1                         ,c5_name,i4time_ob_pr,obstype           ! O
     1                         ,ob_pr_ht_obs,ob_pr_pr_obs              ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs                ! O
     1                         ,ob_pr_t_obs,ob_pr_td_obs               ! O
     1                         ,istatus)                               ! O
        if(istatus .ne. 1)then
            write(6,*)' Bad istatus from read_snd_data'
            go to 999
        endif

        if(n_profiles .eq. 0)then
            write(6,*)' No profiles available via read_snd_data'
            go to 999
        endif

!       Append sounding obs to surface arrays
        do i_pr = 1,n_profiles
            l_good_snd = .false.

            if(nlevels_obs_pr(i_pr) .gt. 0)then
                rlevel_best_wind_diff = 99999.
                rlevel_best_wind      = r_missing_data
                ob_best_wind          = r_missing_data
                ilevel_best_wind      = 0

                rlevel_best_temp_diff = 99999.
                rlevel_best_temp      = r_missing_data
                ob_best_temp          = r_missing_data
                ob_best_dwpt          = r_missing_data
                ilevel_best_temp      = 0

                do il = 1,nlevels_obs_pr(i_pr)
                    rlevel_temp_diff = abs( 
     1                   (ob_pr_ht_obs(i_pr,il)-elev_pr(i_pr)) - 2.0 )

                    rlevel_wind_diff = abs( 
     1                   (ob_pr_ht_obs(i_pr,il)-elev_pr(i_pr)) - 10.0 )

                    if(rlevel_temp_diff .lt. rlevel_best_temp_diff)then
                      if(ob_pr_t_obs(i_pr,il) .ne. r_missing_data)then
                        rlevel_best_temp_diff = rlevel_temp_diff
                        ilevel_best_temp      = il
                        rlevel_best_temp      = ob_pr_ht_obs(i_pr,il)
                      endif
                    endif

                    if(rlevel_wind_diff .lt. rlevel_best_wind_diff)then       
                      if(ob_pr_u_obs(i_pr,il) .ne. r_missing_data
     1             .and. ob_pr_v_obs(i_pr,il) .ne. r_missing_data)then
                        rlevel_best_wind_diff = rlevel_wind_diff
                        rlevel_best_wind      = ob_pr_ht_obs(i_pr,il)
                        ilevel_best_wind      = il
                      endif
                    endif

                enddo ! il

                if(rlevel_best_wind_diff .le. 6.0 .or. 
     1             rlevel_best_temp_diff .le. 2.0)then
                    write(6,*)' Good Sounding Elev     ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,rlevel_best_temp,rlevel_best_temp_diff       
     1                       ,rlevel_best_wind,rlevel_best_wind_diff
                    l_good_snd = .true.
                    n_good_snd = n_good_snd + 1
                elseif(i_pr .le. 100)then
                    write(6,*)' Unusable Sounding Elev ',i_pr
     1                       ,c5_name(i_pr)
     1                       ,rlevel_best_temp,rlevel_best_temp_diff       
     1                       ,rlevel_best_wind,rlevel_best_wind_diff
                endif
            endif

!           if(.false.)then
            if(l_good_snd)then ! This ob is good enough to append
                n_obs_g = n_obs_g + 1
                n_obs_b = n_obs_b + 1

                lat_s(n_obs_b) = lat_pr(i_pr)
                lon_s(n_obs_b) = lon_pr(i_pr)
                elev_s(n_obs_b) = elev_pr(i_pr)
                stations(n_obs_b) = c5_name(i_pr)
                reptype(n_obs_b) = obstype(i_pr)

!               Initialize values to badflag
                t_s(n_obs_b) = badflag
                td_s(n_obs_b) = badflag
                rh_s(n_obs_b) = badflag
                dd_s(n_obs_b) = badflag
                ff_s(n_obs_b) = badflag
                ddg_s(n_obs_b) = badflag
                ffg_s(n_obs_b) = badflag
                alt_s(n_obs_b) = badflag
                pstn_s(n_obs_b) = badflag
                pmsl_s(n_obs_b) = badflag
                delpch(n_obs_b) = badflag
                delp(n_obs_b) = badflag
                vis_s(n_obs_b) = badflag
                solar_s(n_obs_b) = badflag
                sfct(n_obs_b) = badflag
                sfcm(n_obs_b) = badflag
                pcp1(n_obs_b) = badflag
                pcp3(n_obs_b) = badflag
                pcp6(n_obs_b) = badflag
                pcp24(n_obs_b) = badflag
                snow(n_obs_b) = badflag
                kloud_s(n_obs_b) = badflag
                max24t(n_obs_b) = badflag
                min24t(n_obs_b) = badflag

                t_ea(n_obs_b) = badflag
                td_ea(n_obs_b) = badflag
                rh_ea(n_obs_b) = badflag
                dd_ea(n_obs_b) = badflag
                ff_ea(n_obs_b) = badflag
                alt_ea(n_obs_b) = badflag
                p_ea(n_obs_b) = badflag
                vis_ea(n_obs_b) = badflag
                solar_ea(n_obs_b) = badflag
                sfct_ea(n_obs_b) = badflag
                sfcm_ea(n_obs_b) = badflag
                pcp_ea(n_obs_b)  = badflag
                snow_ea(n_obs_b) = badflag
                store_amt(n_obs_b,:) = '     '
                store_hgt(n_obs_b,:) = badflag

!               Add temperature ob into arrays 
                if(ob_pr_t_obs(i_pr,ilevel_best_temp) 
     1                                         .ne. r_missing_data)then
                    t_ob_c = ob_pr_t_obs(i_pr,ilevel_best_temp)     
                    t_s(n_obs_b) = c_to_f(t_ob_c)
                    write(6,*)' Good temp ob',t_s(n_obs_b)
                endif      

!               Add dewpoint ob into arrays 
                if(ob_pr_td_obs(i_pr,ilevel_best_temp) 
     1                                         .ne. r_missing_data)then       
                    td_ob_c = ob_pr_td_obs(i_pr,ilevel_best_temp)     
                    td_s(n_obs_b) = c_to_f(td_ob_c)
                    write(6,*)' Good dwpt ob',td_s(n_obs_b)
                endif      

!               Add wind ob into arrays
                if(ob_pr_u_obs(i_pr,ilevel_best_wind).ne.r_missing_data      
     1       .and. ob_pr_v_obs(i_pr,ilevel_best_wind).ne.r_missing_data
     1                                                         )then
                    call uv_to_disp(ob_pr_u_obs(i_pr,ilevel_best_wind)
     1                             ,ob_pr_v_obs(i_pr,ilevel_best_wind)
     1                             ,dd_s(n_obs_b)
     1                             ,speed_ms )
                    ff_s(n_obs_b) = speed_ms / mspkt    
                    write(6,*)' Good wind ob',dd_s(n_obs_b)
     1                                       ,ff_s(n_obs_b)      
                endif

            endif ! good sounding

        enddo ! i_pr

 999	write(6,*)' # of good soundings = ',n_good_snd,n_obs_g,n_obs_b       

        return
        end

