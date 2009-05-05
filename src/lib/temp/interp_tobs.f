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

        subroutine interp_tobs_to_laps(ob_pr_ht_obs,ob_pr_t_obs,        ! I
     1                             ob_pr_pr_obs,                        ! I
     1                             t_diff,temp_bkg_3d,                  ! I
     1                             t_interp,                            ! O
     1                             i_pr,iwrite,level,l_3d,              ! I
     1                             nlevels_obs,                         ! I
     1                             lat_pr,lon_pr,i_ob,j_ob,             ! I
     1                             ni,nj,nk,                            ! I
     1                             max_rs,max_rs_levels,r_missing_data, ! I
     1                             pres_3d,                             ! I
     1                             heights_3d)                          ! I

!       Profiler Stuff
        real lat_pr(max_rs)
        real lon_pr(max_rs)

!       Sounding Observations
        integer nlevels_obs(max_rs)
        real ob_pr_ht_obs(max_rs,max_rs_levels)
        real ob_pr_pr_obs(max_rs,max_rs_levels) ! mb
        real ob_pr_t_obs(max_rs,max_rs_levels)

        logical l_3d ! If true, we require obs to be near a vertical level
                     ! and single (or multiple) level obs are allowed

                     ! If false, obs aren't required to be near a vertical
                     ! level. Only multiple level obs are allowed.

        real heights_3d(ni,nj,nk)
        real pres_3d(ni,nj,nk)
        real temp_bkg_3d(ni,nj,nk)

        character*2 c2_obtype

        t_interp = r_missing_data

        rk_delt_min = 1e10

!  ***  Interpolate/map the tsnd observations to the input LAPS level *******
        if(l_3d)then

            do i_obs = 1,nlevels_obs(i_pr)

              height_ob = ob_pr_ht_obs(i_pr,i_obs)
              if(height_ob .eq. r_missing_data)then
                  write(6,*)
     1                  ' ERROR: ob ht missing in interp_tobs_to_laps'    
                  stop
              endif

!             Experimental use of observation pressure
              if(.false.)then
!             if(ob_pr_pr_obs(i_pr,i_obs) .ne. r_missing_data)then
                c2_obtype = 'pr'
                pres_ob_pa = ob_pr_pr_obs(i_pr,i_obs) * 100.
                rk_ob = rlevel_of_logfield(pres_ob_pa,pres_3d,ni,nj,nk
     1                                    ,i_ob,j_ob,istatus)


              else
                c2_obtype = 'ht'
                rk_ob = height_to_zcoord2(height_ob,heights_3d,ni,nj,nk       
     1                                   ,i_ob,j_ob,istatus)

              endif

              rk_delt = rk_ob - float(level)

              if(abs(rk_delt) .le. 0.5 .and.        ! Ob is near grid point
     1           abs(rk_delt) .lt. rk_delt_min)then ! Closest ob      

                 rk_delt_min = abs(rk_delt)

                 k1 = level

!                Determine 2nd level that will be used to bracket the ob
                 if(level .eq. 1)then
                     k2 = level + 1
                 elseif(level .eq. nk)then
                     k2 = level - 1
                 elseif(rk_delt .lt. 0.)then
                     k2 = level - 1
                 elseif(rk_delt .ge. 0.)then
                     k2 = level + 1
                 endif

                 if(   (float(level) - rk_ob) 
     1               * (float(k2)    - rk_ob) .gt. 0.)then
                     write(6,*)' Error: k1/k2 does not bracket ob in '
     1                        ,'interp_tsnd_to_laps'
                     return
                 endif

                 if(c2_obtype .eq. 'ht')then

                     t_lapse = ( temp_bkg_3d(i_ob,j_ob,k2) - 
     1                           temp_bkg_3d(i_ob,j_ob,level) ) 
     1                       / ( heights_3d(i_ob,j_ob,k2)     - 
     1                           heights_3d(i_ob,j_ob,level)     ) 

                     temp_ob = ob_pr_t_obs(i_pr,i_obs)
                     height_ob_diff = height_ob 
     1                              - heights_3d(i_ob,j_ob,level)   

!                    The temperature ob now has an "interpolated" or corrected
!                    value by applying the model lapse rate between the ob and 
!                    the laps grid level. The new value is the estimated 
!                    temperature at the laps grid point.

                     t_interp = temp_ob - height_ob_diff * t_lapse

                     if(iwrite .eq. 1)
     1                   write(6,1)level,rk_ob,height_ob,temp_ob,t_lapse
     1                            ,t_interp      
 1                   format(' level rk_ob height_ob temp_ob t_lapse '
     1                     ,'t_interp'     
     1                      /i5,f8.3,f10.1,f8.3,f10.6,f8.3)

                 else ! 'pr' obtype 
                     t_lapse = ( temp_bkg_3d(i_ob,j_ob,k2) - 
     1                           temp_bkg_3d(i_ob,j_ob,level) ) 
     1                       / ( pres_3d(i_ob,j_ob,k2)     - 
     1                           pres_3d(i_ob,j_ob,level)     ) 

                     temp_ob = ob_pr_t_obs(i_pr,i_obs)
                     pres_ob_diff = pres_ob_pa 
     1                              - pres_3d(i_ob,j_ob,level)   

!                    The temperature ob now has an "interpolated" or corrected
!                    value by applying the model lapse rate between the ob and 
!                    the laps grid level. The new value is the estimated 
!                    temperature at the laps grid point.

                     t_interp = temp_ob - pres_ob_diff * t_lapse

                     if(iwrite .eq. 1)
     1                   write(6,2)level,rk_ob,pres_ob_pa,temp_ob
     1                            ,t_lapse,t_interp      
 2                   format(' level rk_ob height_ob temp_ob t_lapse '
     1                     ,'t_interp'     
     1                      /i5,f8.3,f10.1,f8.3,f10.6,f8.3)

                 endif

!                Correct for the time lag
                 t_interp = t_interp + t_diff

               endif ! Ob is near to desired level

            enddo ! level (iobs)

        else ! not l_3d

          if(nlevels_obs(i_pr) .eq. 1)then
!           To avoid this warning set l_3d to TRUE or supply multiple levels
            write(6,*)
     1        'WARNING in interp_tobs - ignoring single level profile'
          endif

          do i_obs = 1,nlevels_obs(i_pr)

            if(i_obs .gt. 1)then

              if(ob_pr_ht_obs(i_pr,i_obs-1) .le. 
     1                                       heights_3d(i_ob,j_ob,level) 
     1                                .AND.
     1           ob_pr_ht_obs(i_pr,i_obs  )   .ge. 
     1                                       heights_3d(i_ob,j_ob,level)
     1                                                             )then

                 h_lower  = ob_pr_ht_obs(i_pr,i_obs-1)
                 h_upper  = ob_pr_ht_obs(i_pr,i_obs  )

                 height_diff = h_upper - h_lower

                 fracl = (h_upper - heights_3d(i_ob,j_ob,level)) 
     1                  / height_diff
                 frach = 1.0 - fracl

                 t_interp = ob_pr_t_obs(i_pr,i_obs)   * frach
     1                    + ob_pr_t_obs(i_pr,i_obs-1) * fracl


!                Correct for the time lag
                 t_interp = t_interp + t_diff

               endif ! obs bracket laps level

            endif  ! (iobs > 1)

          enddo ! level (iobs)
       
        endif ! l_3d

        return

        end

