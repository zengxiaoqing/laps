
        subroutine read_pro_data(lun,i4time_prof,ext
     1                         ,MAX_PR,MAX_PR_LEVELS      
     1                         ,n_profiles
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr
     1                         ,c5_name,i4time_ob_pr,obstype
     1                         ,ob_pr_ht_obs
     1                         ,ob_pr_u_obs,ob_pr_v_obs
     1                         ,istatus)

cdoc    Returns wind profile data from the PRO file

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_di_obs(MAX_PR_LEVELS)
        real ob_pr_sp_obs(MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob
        character ext*(*)

        istatus = 0

        rcycles_pr = 0

        call open_lapsprd_file_read(lun,i4time_prof,ext,istatus)
        if(istatus .ne. 1)go to 420

        do i_pr = 1,max_pr

            read(lun,401,err=430,end=450)
     1         ista,nlevels_obs_pr(i_pr),lat_pr(i_pr),lon_pr(i_pr)
     1         ,elev_pr(i_pr),c5_name(i_pr),a9time_ob,obstype(i_pr)
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,3x,a9,1x,a8)

            if(nlevels_obs_pr(i_pr) .gt. MAX_PR_LEVELS)then
                write(6,*)' ERROR: too many profiler (.pro) levels '
     1                   ,i_pr,nlevels_obs_pr(i_pr),MAX_PR_LEVELS
                goto430
            endif

            call i4time_fname_lp(a9time_ob,i4time_ob_pr(i_pr),istatus)

            write(6,407)i_pr,ista,nlevels_obs_pr(i_pr),lat_pr(i_pr)
     1                 ,lon_pr(i_pr)
     1                 ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr)
     1                 ,a9time_ob,obstype(i_pr)
407         format(/' Profile #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a6,3x,a9
     1                          ,1x,a8)

            do level = 1,nlevels_obs_pr(i_pr)

                read(lun,*,err=430,end=430)ob_pr_ht_obs(i_pr,level)
     1           ,ob_pr_di_obs(level)
     1           ,ob_pr_sp_obs(level)

                call disp_to_uv(ob_pr_di_obs(level),
     1                  ob_pr_sp_obs(level),
     1                  ob_pr_u_obs(i_pr,level),
     1                  ob_pr_v_obs(i_pr,level))

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


        subroutine read_pro_metadata(lun,i4time_prof,ext
     1                         ,MAX_PR,MAX_PR_LEVELS      
     1                         ,n_profiles
     1                         ,nlevels_obs_pr,lat_pr,lon_pr,elev_pr
     1                         ,c5_name,i4time_ob_pr,obstype
     1                         ,ob_pr_ht_obs
!    1                         ,ob_pr_u_obs,ob_pr_v_obs
     1                         ,istatus)

cdoc    Returns wind profile metadata from the PRO file

!       Profile Stuff
        real lat_pr(MAX_PR)
        real lon_pr(MAX_PR)
        real elev_pr(MAX_PR)

        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_di_obs(MAX_PR_LEVELS)
        real ob_pr_sp_obs(MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)

        integer i4time_ob_pr(MAX_PR)
        integer nlevels_obs_pr(MAX_PR)

        character*5 c5_name(MAX_PR)
        character*8 obstype(MAX_PR)
        character*9 a9time_ob
        character ext*(*)

        istatus = 0

        rcycles_pr = 0

        call open_lapsprd_file_read(lun,i4time_prof,ext,istatus)
        if(istatus .ne. 1)go to 420

        do i_pr = 1,max_pr

            read(lun,401,err=430,end=450)
     1         ista,nlevels_obs_pr(i_pr),lat_pr(i_pr),lon_pr(i_pr)
     1         ,elev_pr(i_pr),c5_name(i_pr),a9time_ob,obstype(i_pr)
401         format(i12,i12,f11.0,f15.0,f15.0,5x,a5,1x,3x,a9,1x,a8)

            if(nlevels_obs_pr(i_pr) .gt. MAX_PR_LEVELS)then
                write(6,*)' ERROR: too many profiler (.pro) levels '
     1                   ,i_pr,nlevels_obs_pr(i_pr),MAX_PR_LEVELS
                goto430
            endif

            call i4time_fname_lp(a9time_ob,i4time_ob_pr(i_pr),istatus)

            write(6,407)i_pr,ista,nlevels_obs_pr(i_pr),lat_pr(i_pr)
     1                 ,lon_pr(i_pr)
     1                 ,elev_pr(i_pr),rcycles_pr,c5_name(i_pr)
     1                 ,a9time_ob,obstype(i_pr)
407         format(/' Profile #',i3,i6,i5,2f8.2,e10.3,f8.2,1x,a6,3x,a9
     1                          ,1x,a8)

            do level = 1,nlevels_obs_pr(i_pr)

                read(lun,*,err=430,end=430)ob_pr_ht_obs(i_pr,level)
     1           ,ob_pr_di_obs(level)
     1           ,ob_pr_sp_obs(level)

                call disp_to_uv(ob_pr_di_obs(level),
     1                  ob_pr_sp_obs(level),
     1                  ob_pr_u_obs(i_pr,level),
     1                  ob_pr_v_obs(i_pr,level))

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






