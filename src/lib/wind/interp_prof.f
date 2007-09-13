cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
cdis   
cdis
cdis
cdis   
cdis

        subroutine interp_prof(ob_pr_ht_obs,ob_pr_u_obs,ob_pr_v_obs, ! I
     1                             u_diff       , v_diff,            ! I
     1                             u_interp     , v_interp,          ! O
     1                             di_interp    , sp_interp,         ! O
     1                             i_pr,ht,level,nlevels_obs_pr,     ! I
     1                             lat_pr,lon_pr,i_ob,j_ob,          ! I
     1                             r_missing_data,                   ! I
     1                             heights_3d,ni,nj,nk,              ! I
     1                             MAX_PR,MAX_PR_LEVELS,             ! I
     1                             n_vel_grids,istatus)              ! I/O

!************************ARRAYS.FOR******************************************

!       Profiler Observations

        integer nlevels_obs_pr(MAX_PR)
        real ob_pr_ht_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_u_obs(MAX_PR,MAX_PR_LEVELS)
        real ob_pr_v_obs(MAX_PR,MAX_PR_LEVELS)

!**************************** END ARRAYS.FOR ********************************
        real heights_3d(ni,nj,nk)

        u_interp = r_missing_data
        v_interp = r_missing_data
        di_interp = r_missing_data
        sp_interp = r_missing_data

!  ***  Interpolate the profiler observations to the input height *******
        do i_obs = 1,nlevels_obs_pr(i_pr)

          if(i_obs .gt. 1)then

            if(ob_pr_ht_obs(i_pr,i_obs-1) .le. ht .and.
     1         ob_pr_ht_obs(i_pr,i_obs  ) .ge. ht       )then

                h_lower  = ob_pr_ht_obs(i_pr,i_obs-1)
                h_upper =  ob_pr_ht_obs(i_pr,i_obs  )

                height_diff = h_upper - h_lower

                fracl = (h_upper - ht) / height_diff
                frach = 1.0 - fracl

                u_interp = ob_pr_u_obs(i_pr,i_obs) * frach
     1                    +       ob_pr_u_obs(i_pr,i_obs-1) * fracl

                v_interp = ob_pr_v_obs(i_pr,i_obs) * frach
     1                    +       ob_pr_v_obs(i_pr,i_obs-1) * fracl

!               Correct for the time lag
!               u_diff      = du/dt * [t(anal) - t(ob)]
                u_interp = u_interp + u_diff

!               v_diff      = dv/dt * [t(anal) - t(ob)]
                v_interp = v_interp + v_diff

!               Calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

             endif
          endif

        enddo ! level

        if(.true.)return


!       Lower Tail

        if( float(level)
     1    .lt. height_to_zcoord2(ob_pr_ht_obs(i_pr,1),heights_3d
     1                          ,ni,nj,nk,i_ob,j_ob,istatus)
     1                          .and.
     1  (height_to_zcoord2(ob_pr_ht_obs(i_pr,1),heights_3d
     1                          ,ni,nj,nk,i_ob,j_ob,istatus)
     1       - float(level)) .le. 0.5)then

                u_interp  = ob_pr_u_obs(i_pr,1)
                v_interp  = ob_pr_v_obs(i_pr,1)

!               Correct for the time lag
                u_interp = u_interp + u_diff
                v_interp = v_interp + v_diff

!               Calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

        endif

        if(istatus .ne. 1)then
            write(6,*)
     1        ' Warning, out of bounds in height_to_zcoord2/interp_prof'       
        endif

!       Upper Tail

        if( float(level) .gt.
     1  height_to_zcoord2(ob_pr_ht_obs(i_pr,nlevels_obs_pr(i_pr))
     1                  ,heights_3d,ni,nj,nk,i_ob,j_ob,istatus)
     1                           .and.
     1    (float(level) -
     1     height_to_zcoord2(ob_pr_ht_obs(i_pr,nlevels_obs_pr(i_pr))
     1                  ,heights_3d,ni,nj,nk,i_ob,j_ob,istatus))
     1                                          .le. 0.5)then

                u_interp  = ob_pr_u_obs(i_pr,nlevels_obs_pr(i_pr))
                v_interp  = ob_pr_v_obs(i_pr,nlevels_obs_pr(i_pr))

!               Correct for the time lag
                u_interp = u_interp + u_diff
                v_interp = v_interp + v_diff

!               Calculate direction and speed
                call uv_to_disp( u_interp,
     1                           v_interp,
     1                           di_interp,
     1                           sp_interp)

        endif

        if(istatus .ne. 1)then
            write(6,*)' Warning, out of bounds in height_to_zcoord2/inte
     1rp_prof'
        endif

        return

        end

