cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 

        subroutine interp_rass_to_laps(ob_pr_ht_obs,ob_pr_t_obs,
     1                             t_diff,
     1                             t_interp,
     1                             i_pr,level,
     1                             nlevels_obs,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             ni,nj,nk,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             heights_3d)

!       Profiler Stuff
        real lat_pr(max_rs)
        real lon_pr(max_rs)

!       RASS Observations
        integer nlevels_obs(max_rs)
        real ob_pr_ht_obs(max_rs,max_rs_levels)
        real ob_pr_t_obs(max_rs,max_rs_levels)

        real*4 heights_3d(ni,nj,nk)

        t_interp = r_missing_data

!  ***  Interpolate the rass observations to the input LAPS heights *******
        do i_obs = 1,nlevels_obs(i_pr)

          if(i_obs .gt. 1)then

            if(ob_pr_ht_obs(i_pr,i_obs-1) .le. 
     1                                       heights_3d(i_ob,j_ob,level) 
     1                                .AND.
     1       ob_pr_ht_obs(i_pr,i_obs  )   .ge. 
     1                                       heights_3d(i_ob,j_ob,level)
     1                                                             )then

                h_lower  = ob_pr_ht_obs(i_pr,i_obs-1)
                h_upper  = ob_pr_ht_obs(i_pr,i_obs  )

                height_diff = h_upper - h_lower

                fracl = (h_upper - heights_3d(i_ob,j_ob,level)) 
     1                 / height_diff
                frach = 1.0 - fracl

                t_interp = ob_pr_t_obs(i_pr,i_obs)   * frach
     1                   + ob_pr_t_obs(i_pr,i_obs-1) * fracl

!               Correct for the time lag
                t_interp = t_interp + t_diff

             endif
          endif

        enddo ! level

        return

        end


        subroutine interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                             t_interp,p_interp,
     1                             i_pr,k_rass,
     1                             nlevels_obs,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             ni,nj,nk,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             temp_3d,heights_3d,pres_3d)

!       Profiler Stuff
        real lat_pr(max_rs)
        real lon_pr(max_rs)

!       RASS Observations
        integer nlevels_obs(max_rs)
        real ob_pr_ht_obs(max_rs,max_rs_levels)
        real ob_pr_t_obs(max_rs,max_rs_levels)

        real*4 heights_3d(ni,nj,nk)
        real*4 temp_3d(ni,nj,nk)
        real*4 pres_3d(ni,nj,nk)

        t_interp = r_missing_data

!  ***  Interpolate the LAPS temps to the input RASS heights *******
        do k_laps = 2,nk

            if( heights_3d(i_ob,j_ob,k_laps-1) .le. 
     1                                         ob_pr_ht_obs(i_pr,k_rass)
     1                                    .AND.
     1          heights_3d(i_ob,j_ob,k_laps)   .ge. 
     1                                         ob_pr_ht_obs(i_pr,k_rass)
     1                                                             )then

                h_lower = heights_3d(i_ob,j_ob,k_laps-1)
                h_upper = heights_3d(i_ob,j_ob,k_laps)

                height_diff = h_upper - h_lower

                fracl = (h_upper - ob_pr_ht_obs(i_pr,k_rass)) 
     1                 / height_diff
                frach = 1.0 - fracl

                t_interp = temp_3d(i_ob,j_ob,k_laps)   * frach
     1                   + temp_3d(i_ob,j_ob,k_laps-1) * fracl

                p_interp = pres_3d(i_ob,j_ob,k_laps)   * frach
     1                   + pres_3d(i_ob,j_ob,k_laps-1) * fracl

             endif

        enddo ! level

        return

        end

