
        subroutine read_pro_metadata(lun,i4time_prof,ext              ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                  ! I
     1                         ,n_profiles                            ! O
     1                         ,nlevels_obs_pr                        ! I/O
     1                         ,lat_pr,lon_pr,elev_pr                 ! O
     1                         ,c5_name,i4time_ob_pr,obstype          ! O
     1                         ,ob_pr_ht_obs                          ! O
     1                         ,istatus)                              ! O

cdoc    Returns wind profile metadata from the PRO file. 
cdoc    Jacket for read_pro_data.

        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_rms_obs(MAX_PR,MAX_PR_LEVELS)

        real sfc_t(MAX_PR), sfc_p(MAX_PR), sfc_rh(MAX_PR)
        real sfc_u(MAX_PR), sfc_v(MAX_PR) 

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character ext*(*)

        call read_pro_data(     lun,i4time_prof,ext                   ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                  ! I
     1                         ,n_profiles                            ! O
     1                         ,nlevels_obs_pr                        ! I/O
     1                         ,lat_pr,lon_pr,elev_pr                 ! O
     1                         ,c5_name,i4time_ob_pr,obstype          ! O
     1                         ,ob_pr_ht_obs                          ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs               ! O
     1                         ,ob_pr_rms_obs                         ! O
     1                         ,sfc_t,sfc_p,sfc_rh,sfc_u,sfc_v        ! O
     1                         ,istatus)                              ! O

        return
        end 


        subroutine read_pro_data(lun,i4time_prof,ext                  ! I
     1                         ,MAX_PR,MAX_PR_LEVELS                  ! I
     1                         ,n_profiles                            ! O
     1                         ,nlevels_obs_pr                        ! I/O
     1                         ,lat_pr,lon_pr,elev_pr                 ! O
     1                         ,c5_name,i4time_ob_pr,obstype          ! O
     1                         ,ob_pr_ht_obs                          ! O
     1                         ,ob_pr_u_obs,ob_pr_v_obs               ! O
     1                         ,ob_pr_rms_obs                         ! O
     1                         ,sfc_t,sfc_p,sfc_rh,sfc_u,sfc_v        ! O
     1                         ,istatus)                              ! O

cdoc    Returns wind profile data from the PRO file

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_di_obs(MAX_PR_LEVELS)
        real ob_pr_sp_obs(MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS) ! includes sfc wind when valid
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS) ! includes sfc wind when valid
        real ob_pr_rms_obs(MAX_PR,MAX_PR_LEVELS)

!       Note that surface data access is still under development and may not
!       be reliable yet
        real sfc_t(MAX_PR), sfc_p(MAX_PR), sfc_rh(MAX_PR)
        real sfc_u(MAX_PR), sfc_v(MAX_PR) 

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)  ! Includes sfc wind when valid
                                        ! Optionally can be initialized to 0.
                                        ! by the calling program even though 
                                        ! this isn't necessary for the proper
                                        ! operation of this routine
        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character ext*(*)
        character*9 a9time_ob
        character*132 line

        logical l_sfc, l_goodwind

        istatus = 0

        call get_sfc_badflag(badflag,istat_badflag)

        rcycles_pr = 0

        call open_lapsprd_file_read(lun,i4time_prof,ext,istatus)
        if(istatus .ne. 1)go to 420

        do i_pr = 1,max_pr

            read(lun,401,err=430,end=450)
     1           ista,nlevels_in,lat_pr(i_pr),lon_pr(i_pr)
     1          ,elev_pr(i_pr),c5_name(i_pr),a9time_ob,obstype(i_pr)       
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,3x,a9,1x,a8)

            call i4time_fname_lp(a9time_ob,i4time_ob_pr(i_pr),istatus)

            write(6,407)i_pr,ista,nlevels_in,lat_pr(i_pr)
     1                 ,lon_pr(i_pr)
     1                 ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr)
     1                 ,a9time_ob,obstype(i_pr)
407         format(/' Profile #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a6,3x,a9
     1                          ,1x,a8)

            level_out = 0

            do level = 1,nlevels_in
                read(lun,402,end=430)line
 402            format(a)

                read(line,*,err=410,end=410)                 ! Rd sfc if there
     1                              ob_pr_ht_obs_in 
     1                             ,ob_pr_di_obs_in
     1                             ,ob_pr_sp_obs_in
     1                             ,sfc_p_in,sfc_t_in,sfc_rh_in

                sfc_p(i_pr) = sfc_p_in
                sfc_t(i_pr) = sfc_t_in
                sfc_rh(i_pr) = sfc_rh_in
                
                ob_pr_rms_obs_in = 1.0 ! Hardwired until available in PRO file

                l_sfc = .true.

                if(level .gt. 1)then
                    write(6,*)' ERROR: sfc data exists past level 1'       
                    istatus = 0
                    return
                endif

                if( abs(ob_pr_ht_obs_in - elev_pr(i_pr)) .gt. 1.)then
                    write(6,*)' ERROR: elevation disagrees with sfc'     
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)
                    istatus = 0
                    return
                endif                    

                goto415

 410            read(line,*,err=430)ob_pr_ht_obs_in          ! Rd upr lvl only
     1                             ,ob_pr_di_obs_in
     1                             ,ob_pr_sp_obs_in
     1                             ,ob_pr_rms_obs_in
                l_sfc = .false.

                if( abs(ob_pr_ht_obs_in - elev_pr(i_pr)) .lt. 1.)then
                    write(6,*)' Note: elevation agrees with sfc'
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)

                    if(level .gt. 1)then
                        write(6,*)' ERROR: elevation agrees with sfc'
     1                       ,ob_pr_ht_obs_in,elev_pr(i_pr)
                        istatus = 0
                        return
                    else
                        write(6,*)' Note: no P, T, RH present at sfc'       
                        sfc_p(i_pr) = badflag
                        sfc_t(i_pr) = badflag
                        sfc_rh(i_pr) = badflag
                    endif

                endif ! gate elevation matches the surface

 415            if(    ob_pr_di_obs_in .ne. badflag 
     1           .and. ob_pr_sp_obs_in .ne. badflag )then    ! good wind
                    level_out = level_out + 1

                    if(level_out .gt. MAX_PR_LEVELS)then
                        write(6,*)
     1                      ' ERROR: too many profiler (.pro) levels '
     1                      ,i_pr,level_out,MAX_PR_LEVELS
                        goto430
                    endif

                    ob_pr_ht_obs(i_pr,level_out)  = ob_pr_ht_obs_in 
                    ob_pr_di_obs(level_out)       = ob_pr_di_obs_in 
                    ob_pr_sp_obs(level_out)       = ob_pr_sp_obs_in 
                    ob_pr_rms_obs(i_pr,level_out) = ob_pr_rms_obs_in 

                    nlevels_obs_pr(i_pr) = level_out

                    call disp_to_uv(ob_pr_di_obs(level_out),
     1                              ob_pr_sp_obs(level_out),
     1                              ob_pr_u_obs(i_pr,level_out),
     1                              ob_pr_v_obs(i_pr,level_out))

                    if(l_sfc)then
                        sfc_u(i_pr) = ob_pr_u_obs(i_pr,level_out)
                        sfc_v(i_pr) = ob_pr_v_obs(i_pr,level_out)
                    endif

                else                                         ! bad wind
                    if(l_sfc)then
                        sfc_u(i_pr) = badflag
                        sfc_v(i_pr) = badflag
                    else
                        write(6,*)' ERROR: non-sfc wind should be good'       
                        istatus = 0
                        return
                    endif

                endif                

311             format(1x,2i4,5f8.1)
312             continue
            enddo ! level
        enddo ! profiler site

        write(6,*)' ERROR: # of profilers has reached'
     1           ,' dimensions of MAX_PR ',MAX_PR

        n_profiles=MAX_PR
        GO TO 500
c
c       Profiler reading error handling
c
  420   write(6,*) ' Warning, could not open PRO file'
        n_profiles=0
        GO TO 500

  430   write(6,*) ' Error during read of PRO file'
        write(6,*) ' While reading profiler number ',i_pr
        n_profiles=i_pr-1
        GO TO 500
c
  450   CONTINUE ! Used for end of file
        n_profiles=i_pr-1

        istatus = 1

        close(lun)

  500   CONTINUE

        return
        end








