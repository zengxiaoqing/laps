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
        subroutine comp_vr_prof(grid_ra_vel,grid_pr_r,r_missing_data,ni,
     1nj,nk
     1                  ,lat,lon,max_weights_pr,i_profiler_nearest
     1                        ,v_nyquist_2,unfolding_thresh
     1                  ,heights_1d,rlat_radar,rlon_radar,rheight_radar)

!       This routine dealiases the radar data by comparing to the profilers
!       No other QC thresholding is presently done, but could be easily added

!       Steve Albers
!       Jan 1992, A fix now sets the radar velocity to missing if the
!       Profiler reference analysis has missing data. If no de-aliasing
!       of the radar data can be done, it is thrown out.

!       A future fix might include using the first guess analysis as the
!       reference for de-aliasing the radar instead of the preliminary
!       profiler analysis.

!       Dec 1992 - De-aliasing moved to subroutine 'qc_radar_obs'. Turned off
!       here, but still being diagnosed.

!       Apr 1993 - Removed arrays include

        real grid_ra_vel(ni,nj,nk)
        real grid_pr_r(ni,nj,nk)
        real heights_1d(nk)
        real lat(ni,nj)
        real lon(ni,nj)
        real max_weights_pr(ni,nj,nk)
        integer i_profiler_nearest(ni,nj,nk)

        write(6,2)
2       format(/'                  Comparing Radial Velocities to Profil
     1er'
     1  /'    i   j   k   radar  prof   diff    wt   prof  nobs')

        do k = 1,nk
        heights_1d(k) = height_of_level(k)
        do j = 1,nj
        do i = 1,ni

!           Check for presence of radar ob at this grid point
            if(grid_pr_r(i,j,k) .ne. r_missing_data
     1  .and. grid_ra_vel(i,j,k) .ne. r_missing_data)then

                nobs1 = nobs1 + 1

!               Compare radar radial velocity to profiler first guess
                diff = grid_pr_r(i,j,k) - grid_ra_vel(i,j,k)

!               Check for folded radar data (up to +/- 1 Nyquist velocity)
!               Test if residual is more than about 1.3 V Nyquist and
!                                   less than about 2.7 V Nyquist
                if(abs(abs(diff)-v_nyquist_2) .lt. unfolding_thresh)then
                    call latlon_to_radar(lat(i,j),lon(i,j),heights_1d(k)
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                    velold = grid_ra_vel(i,j,k)

!                   Adjust the velocity value
                    if(.false.)then
                        r_nyquist_number = nint(diff/v_nyquist_2)
                        diff = diff - r_nyquist_number * v_nyquist_2
                        grid_ra_vel(i,j,k) =
     1              grid_ra_vel(i,j,k)
     1            + r_nyquist_number * v_nyquist_2
                    endif

                    write(6,102)i,j,k
     1                                  ,grid_pr_r(i,j,k)
     1                                  ,velold
     1                                  ,grid_ra_vel(i,j,k)
     1                                  ,diff
     1                                  ,nint(azimuth)
     1                                  ,nint(slant_range/1000.)
     1                                  ,elev

102                 format(1x,'Folding at'
     1            ,i3,i3,i3,' vp,vr,vrnw',3f6.1,' df',f6.1
     1            ,' azran=',i3,'/',i3,' el=',f4.1)

                endif


                residual1 = residual1 + diff**2
!d              if( float(nobs1)/50. .eq. int(float(nobs1)/50.) )
!d                  write(6,101)i,j,k,grid_ra_vel(i,j,k),grid_pr_r(i,j,k),
!d      1                       diff,max_weights_pr(i,j,k),
!d      1                       i_profiler_nearest(i,j,k),nobs2

                if(max_weights_pr(i,j,k) .ge. .999)then
                    sumsqpr = sumsqpr + grid_pr_r(i,j,k) ** 2
                    sumsqra = sumsqra + grid_ra_vel(i,j,k) ** 2
                    nobs2 = nobs2 + 1
                    residual2 = residual2 + diff ** 2
                    write(6,101)i,j,k,grid_ra_vel(i,j,k),grid_pr_r(i,j,k
     1),
     1                  diff,max_weights_pr(i,j,k),
     1                  i_profiler_nearest(i,j,k),nobs2
101                 format(1x,3i4,3f7.1,f7.3,2i6)
                endif


            else ! Missing data in profiler anal or radar data
!               if(k .eq. 13)then
!                   write(6,*)i,j,k,' Radar Present, profiler?, k = 13'
!               endif

!               if(grid_pr_r(i,j,k) .eq. r_missing_data)then
!                   write(6,*)i,j,k,' Radar Present, profiler absent'
!               endif

!               Make sure radar data is not used if there is no profiler
!               data present to de-alias it.
                if(.false.)then ! QCing mode
                    grid_ra_vel(i,j,k) = r_missing_data
                endif

            endif

        enddo ! i
        enddo ! j
        enddo ! k

        if(nobs1 .gt. 0)then
            rms1 = sqrt(residual1/nobs1)
        else
            rms1 = 0.
        endif

        if(nobs2 .gt. 0)then
            rms2 = sqrt(residual2/nobs2)
            rmspr2 = sqrt(sumsqpr/nobs2)
            rmsra2 = sqrt(sumsqra/nobs2)
        else
            rms2 = 0.
            rmspr2 = 0.
            rmsra2 = 0.
        endif

        write(6,*)' RMS between radial velocities & profiler anal= ',nob
     1s1,rms1
        write(6,*)' RMS Of profiler,radar = ',rmspr2,rmsra2
        write(6,*)' RMS between radial velocities & profiler obs = ',nob
     1s2,rms2

        write(15,*)' RMS between radial velocities & profiler anal = ',n
     1obs1,rms1
        write(15,*)' RMS Of profiler,radar = ',rmspr2,rmsra2
        write(15,*)' RMS between radial velocities & profiler obs = ',no
     1bs2,rms2

        return

        end

