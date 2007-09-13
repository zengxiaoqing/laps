
        subroutine interp_laps_to_rass(ob_pr_ht_obs,ob_pr_t_obs,
     1                             t_interp,p_interp,
     1                             i_pr,k_rass,
     1                             nlevels_obs,
     1                             lat_pr,lon_pr,i_ob,j_ob,
     1                             ni,nj,nk,
     1                             max_rs,max_rs_levels,r_missing_data,
     1                             temp_bkg_3d,heights_3d,pres_3d)

!       Profiler Stuff
        real lat_pr(max_rs)
        real lon_pr(max_rs)

!       RASS Observations
        integer nlevels_obs(max_rs)
        real ob_pr_ht_obs(max_rs,max_rs_levels)
        real ob_pr_t_obs(max_rs,max_rs_levels)

        real heights_3d(ni,nj,nk)
        real temp_bkg_3d(ni,nj,nk)
        real pres_3d(ni,nj,nk)

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

                t_interp = temp_bkg_3d(i_ob,j_ob,k_laps)   * frach
     1                   + temp_bkg_3d(i_ob,j_ob,k_laps-1) * fracl

                p_interp = pres_3d(i_ob,j_ob,k_laps)   * frach
     1                   + pres_3d(i_ob,j_ob,k_laps-1) * fracl

             endif

        enddo ! level

        return

        end

