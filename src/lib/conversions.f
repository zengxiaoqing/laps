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

!       Miscellaneous Conversion Routines

!       Steve Albers           1990

        function omega_to_w(omega,pressure_pa)

        real*4 omega_to_w

        real*4 scale_height
        parameter (scale_height = 8000.)

        omega_to_w = - (omega / pressure_pa) * scale_height

        return
        end

        function w_to_omega(w,pressure_pa)

        real*4 w_to_omega

        real*4 scale_height
        parameter (scale_height = 8000.)

        w_to_omega = - (w * pressure_pa) / scale_height

        return
        end

        subroutine latlon_to_radar(lat_grid,lon_grid,height_grid
     1                  ,azimuth,slant_range,elev
     1                  ,rlat_radar,rlon_radar,rheight_radar)
        include 'trigd.inc'

        implicit real*4 (a-z)

        if(rlat_radar .eq. 0.0)then
            write(6,*)' Warning, Radar Coords NOT Initialized'
        endif

        rpd = 3.141592653589/180.
        mpd = 111194.
        radius_earth = 6371.e3
        radius_earth_8_thirds = 6371.e3 * 2.6666666
        diff = 0.

        delta_lat = lat_grid - rlat_radar
        delta_lon = lon_grid - rlon_radar


        cos_factor =  cosd( 0.5 * (rlat_radar + lat_grid ) )


        delta_y = delta_lat * mpd

        delta_x = delta_lon * mpd * cos_factor

        height_factor =
     1   (radius_earth + 0.5 * (rheight_radar + height_grid))
     1  /radius_earth

        hor_dist = sqrt(delta_x**2 + delta_y**2) * height_factor

        delta_z = height_grid - rheight_radar

        slant_range = sqrt(hor_dist**2 + delta_z**2)

        if(hor_dist .gt. 0.0)then
            azimuth = atan3(delta_x,delta_y)/rpd
        else
            azimuth = 0.0
        endif

        curvature = hor_dist **2 / radius_earth_8_thirds
        height_diff = delta_z - curvature


        if(slant_range .gt. 0.0)then
            elev = atan2 ( height_diff , hor_dist ) / rpd
        else
            elev = 0.0
        endif


c       write(6,1)slant_range,delta_x,delta_y,azimuth
c       1                       ,hor_dist,curvature
1       format(1x,7f10.0)

        return

        end



        function height_to_zcoord(height_m,istatus) ! Using standard atmosphere

!       Steve Albers FSL
!       Note that this routine works with the standard atmosphere.
!       When the vertical grid is pressure, the height is converted to
!       pressure, then the interpolation to the vertical grid is performed.
!       Thus if the height is midway between two LAPS levels in height space,
!       the value of height_to_zcoord will not have a fraction of 0.5.
!       If the pressure is midway between two LAPS levels, then the
!       value of height_to_zcoord will have a fraction of 0.5.

        implicit real*4 (a-z)

        integer*4 istatus

        include 'lapsparms.cmn'

        logical ltest_vertical_grid

        call get_laps_config('nest7grid',istatus)
        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            return
        endif

        if(ltest_vertical_grid('HEIGHT'))then
!           height_to_zcoord = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            pressure_pa = ztopsa(height_m) * 100.
            if(pressure_pa .eq. 9999900.)pressure_pa = 1e-2
            height_to_zcoord = (PRESSURE_0_L - pressure_pa)
     1                          / PRESSURE_INTERVAL_L

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function height_to_zcoord2(height_m,heights_3d
     1                                          ,ni,nj,nk,i,j,istatus)

!       WARNING: This routine is designed to be efficient. As a result, it
!       will fail if the input heights differ too much from the standard
!       atmosphere. In such a case, the istatus will be returned as zero.

!       1994 Steve Albers FSL

!       Note that this routine works with the real atmosphere.
!       The type of interpolation is similar to that in 'height_to_zcoord'.
!       When the vertical grid is pressure, the height is converted to
!       pressure, then the interpolation to the vertical grid is performed.
!       Thus if the height is midway between two LAPS levels in height space,
!       the value of 'height_to_zcoord2' will not have a fraction of 0.5.
!       If the pressure is midway between two LAPS levels, then the
!       value of 'height_to_zcoord2' will have a fraction of 0.5.

        implicit real*4 (a-z)

        integer i,j,k,ni,nj,nk,kref,istatus

        real*4 heights_3d(ni,nj,nk)

        logical ltest_vertical_grid

        istatus = 1

        if(ltest_vertical_grid('HEIGHT'))then
            height_to_zcoord2 = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_to_zcoord2 = nk+1 ! Default value is off the grid

          ! Standard Atmosphere Guess + a cushion
!           This must always be >= height_to_zcoord2
            kref = min(int(height_to_zcoord((height_m+600.)*1.2,istatus)
     1),nk)

            heights_above = heights_3d(i,j,kref)

            if(height_m .gt. heights_above)then
                istatus = 0
                goto999
            endif

            do k = kref-1,1,-1
                if(heights_above     .ge. height_m .and.
     1           heights_3d(i,j,k) .le. height_m         )then
                    thickness = heights_above - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low  = zcoord_of_level(k)
                    pressure_high = zcoord_of_level(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    pressure = pressure_low * exp(diff_log_space*fractio
     1n)
                    height_to_zcoord2 = zcoord_of_pressure(pressure)

!                   if(j .eq. 29)then
!                       write(6,*)' height_to_zcoord2: kref,k,kref-k+1'
!       1                                             ,kref,k,kref-k+1
!                   endif

                    goto999

                endif

                heights_above = heights_3d(i,j,k)

            enddo ! k

            istatus = 0
            height_to_zcoord2 = 0
            write(6,101)kref,height_m,heights_3d(i,j,1)
101         format('  Error: below domain in height_to_zcoord2, kref,h,h
     1(1)',
     1             i3,2e11.4)

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

999     return
        end


        function height_to_zcoord2_lin(height_m,heights_3d
     1                                  ,ni,nj,nk,i,j,istatus)

!       WARNING: This routine is designed to be efficient. As a result, it
!       will fail if the input heights differ too much from the standard
!       atmosphere. In such a case, the istatus will be returned as zero.

!       1994 Steve Albers

        implicit real*4 (a-z)

        integer i,j,k,ni,nj,nk,kref,istatus

        real*4 heights_3d(ni,nj,nk)

        logical ltest_vertical_grid

        istatus = 1

        if(ltest_vertical_grid('HEIGHT'))then
            height_to_zcoord2_lin = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_to_zcoord2_lin = nk+1 ! Default value is off the grid

          ! Standard Atmosphere Guess + a cushion
!           This must always be >= height_to_zcoord2_lin
            kref = min(int(height_to_zcoord((height_m+600.)*1.2,istatus)
     1),nk)

            heights_above = heights_3d(i,j,kref)

            if(height_m .gt. heights_above)then
                istatus = 0
                goto999
            endif

            do k = kref-1,1,-1
                if(heights_above     .ge. height_m .and.
     1           heights_3d(i,j,k) .le. height_m         )then
                    thickness = heights_above - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness

                    height_to_zcoord2_lin = k + fraction

!                   if(j .eq. 29)then
!                       write(6,*)' height_to_zcoord2_lin: kref,k,kref-k+1'
!       1                                             ,kref,k,kref-k+1
!                   endif

                    goto999

                endif

                heights_above = heights_3d(i,j,k)

            enddo ! k

            istatus = 0
            height_to_zcoord2_lin = 0
            write(6,101)kref,height_m,heights_3d(i,j,1)
101         format(
     1    ' Error: below domain in height_to_zcoord2_lin, kref,h,h(1)'
     1             ,i3,2e11.4)

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

999     return
        end

        function height_to_zcoord3(height_m,heights_3d,zcoords_1d
     1                                  ,ni,nj,nk,i,j,istatus)


!       Note that this routine works with the real atmosphere.
!       The type of interpolation is similar to that in 'height_to_zcoord'.
!       When the vertical grid is pressure, the height is converted to
!       pressure, then the interpolation to the vertical grid is performed.
!       Thus if the height is midway between two LAPS levels in height space,
!       the value of 'height_to_zcoord3' will not have a fraction of 0.5.
!       If the pressure is midway between two LAPS levels, then the
!       value of 'height_to_zcoord3' will have a fraction of 0.5.

        implicit real*4 (a-z)

        integer i,j,k,ni,nj,nk,kref,istatus

        real*4 heights_3d(ni,nj,nk)
        real*4 zcoords_1d(nk)

        logical ltest_vertical_grid

        istatus = 1

        if(ltest_vertical_grid('HEIGHT'))then
            height_to_zcoord3 = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_to_zcoord3 = nk+1 ! Default value is off the grid

          ! Standard Atmosphere Guess + a cushion
            kref = min(int(height_to_zcoord((height_m+600.)*1.2,istatus)
     1),nk)

            heights_above = heights_3d(i,j,kref)

            if(height_m .gt. heights_above)then
                istatus = 0
                goto999
            endif

            do k = kref-1,1,-1
                if(heights_above     .ge. height_m .and.
     1           heights_3d(i,j,k) .le. height_m         )then
                    thickness = heights_above - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low  = zcoords_1d(k)
                    pressure_high = zcoords_1d(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    pressure = pressure_low * exp(diff_log_space*fractio
     1n)
                    height_to_zcoord3 = zcoord_of_pressure(pressure)

!                   if(j .eq. 29)then
!                       write(6,*)' height_to_zcoord3: kref,k,kref-k+1'
!       1                                             ,kref,k,kref-k+1
!                   endif

                    goto999

                endif

                heights_above = heights_3d(i,j,k)

            enddo ! k

            istatus = 0
            height_to_zcoord3 = 0
            write(6,101)kref,height_m,heights_3d(i,j,1)
101         format('  Error: below domain in height_to_zcoord3, kref,h,h
     1(1)',
     1             i3,2e11.4)

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

999     return
        end

        function height_to_zcoord4(height_m,heights_3d,zcoords_1d,kref
     1                                  ,ni,nj,nk,i,j,istatus)


!       Note that this routine works with the real atmosphere.
!       The type of interpolation is similar to that in 'height_to_zcoord'.
!       When the vertical grid is pressure, the height is converted to
!       pressure, then the interpolation to the vertical grid is performed.
!       Thus if the height is midway between two LAPS levels in height space,
!       the value of 'height_to_zcoord4' will not have a fraction of 0.5.
!       If the pressure is midway between two LAPS levels, then the
!       value of 'height_to_zcoord4' will have a fraction of 0.5.

        implicit real*4 (a-z)

        integer i,j,k,ni,nj,nk,kref,istatus

        real*4 heights_3d(ni,nj,nk)
        real*4 zcoords_1d(nk)

        logical ltest_vertical_grid

        istatus = 1

        if(ltest_vertical_grid('HEIGHT'))then
            height_to_zcoord4 = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_to_zcoord4 = nk+1 ! Default value is off the grid

            heights_above = heights_3d(i,j,kref)

            if(height_m .gt. heights_above)then
                istatus = 0
                goto999
            endif

            do k = kref-1,1,-1
                if(heights_above     .ge. height_m .and.
     1           heights_3d(i,j,k) .le. height_m         )then
                    thickness = heights_above - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low  = zcoords_1d(k)
                    pressure_high = zcoords_1d(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    pressure = pressure_low * exp(diff_log_space*fractio
     1n)
                    height_to_zcoord4 = zcoord_of_pressure(pressure)

!                   if(j .eq. 29)then
!                       write(6,*)' height_to_zcoord4: kref,k,kref-k+1'
!       1                                             ,kref,k,kref-k+1
!                   endif

                    goto999

                endif

                heights_above = heights_3d(i,j,k)

            enddo ! k

            istatus = 0
            height_to_zcoord4 = 0
            write(6,101)kref,height_m,heights_3d(i,j,1)
101         format('  Error: below domain in height_to_zcoord4, kref,h,h
     1(1)',
     1             i3,2e11.4)

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

999     return
        end


        function height_to_pressure(height_m,heights_3d
     1                          ,pressures_1d,ni,nj,nk,i,j)

        implicit real*4 (a-z)
        integer i,j,k,ni,nj,nk

        logical ltest_vertical_grid

        real*4 heights_3d(ni,nj,nk)
        real*4 pressures_1d(nk)

        if(ltest_vertical_grid('HEIGHT'))then
            height_to_zcoord2 = height_m / HEIGHT_INTERVAL

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_to_pressure = -999. ! Default value is off the grid
            do k = 1,nk-1
                if(heights_3d(i,j,k+1) .ge. height_m .and.
     1           heights_3d(i,j,k)   .le. height_m         )then
                    thickness = heights_3d(i,j,k+1) - heights_3d(i,j,k)
                    fraction = (height_m - heights_3d(i,j,k))/thickness
                    pressure_low = pressures_1d(k)
                    pressure_high = pressures_1d(k+1)
                    diff_log_space = log(pressure_high/pressure_low)
                    height_to_pressure
     1          = pressure_low * exp(diff_log_space*fraction)
                endif

            enddo ! k

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        return
        end


        subroutine pressure_to_height(pres_pa,heights_3d
     1                             ,ni,nj,nk,i,j,height_out,istatus)

        real*4 heights_3d(ni,nj,nk)

        logical ltest_vertical_grid

        if(ltest_vertical_grid('HEIGHT'))then
            istatus = 0
            return

        elseif(ltest_vertical_grid('PRESSURE'))then
            if(pres_pa .lt. pressure_of_level(nk) .or.
     1         pres_pa .gt. pressure_of_level(1)       )then
                istatus = 0
                return
            endif


            rk = zcoord_of_logpressure(pres_pa)
            if(rk .lt. 1. .or. rk .gt. float(nk))then
                istatus = 0
                return

            else
                if(rk .eq. float(nk))then
                    klow = nk - 1
                    khigh = nk
                else
                    klow = int(rk)
                    khigh = klow + 1
                endif

                frac = rk - klow

                height_out = heights_3d(i,j,klow)  * (1. - frac)
     1                     + heights_3d(i,j,khigh) * (frac)

                istatus = 1
                return

            endif ! rk is within domain

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        istatus = 0
        return
        end


        function height_of_level(level)

        implicit real*4 (a-z)

        integer*4 level,istatus

        logical ltest_vertical_grid

        if(ltest_vertical_grid('HEIGHT'))then
            height_of_level = HEIGHT_INTERVAL * level

        elseif(ltest_vertical_grid('PRESSURE'))then
            height_of_level = psatoz(pressure_of_level(level) * .01)

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function zcoord_of_level(level)

        implicit real*4 (a-z)

        integer*4 level, istatus

        include 'lapsparms.cmn'

        logical ltest_vertical_grid

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        if(ltest_vertical_grid('HEIGHT'))then
!           zcoord_of_level = height_interval * level

        elseif(ltest_vertical_grid('PRESSURE'))then
            zcoord_of_level = PRESSURE_0_L
     1                  - PRESSURE_INTERVAL_L * level

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        function pressure_of_level(level)

        implicit real*4 (a-z)

        integer*4 level, istatus

        include 'lapsparms.cmn'

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        pressure_of_level = PRESSURE_0_L
     1                  - PRESSURE_INTERVAL_L * level

        istatus = 1
        return
        end


        function pressure_of_rlevel(rlevel)

        implicit real*4 (a-z)

        integer istatus

        include 'lapsparms.cmn'

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        pressure_of_rlevel = PRESSURE_0_L
     1                  - PRESSURE_INTERVAL_L * rlevel

        istatus = 1
        return
        end


        function zcoord_of_pressure(pres_pa)

        implicit real*4 (a-z)

        include 'lapsparms.cmn'

!       logical ltest_vertical_grid

        integer istatus

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

!       if(ltest_vertical_grid('HEIGHT'))then

!       elseif(ltest_vertical_grid('PRESSURE'))then
            zcoord_of_pressure =
     1  (PRESSURE_0_L - pres_pa) / PRESSURE_INTERVAL_L

!       else
!           write(6,*)' Error, vertical grid not supported,'
!    1               ,' this routine supports PRESSURE or HEIGHT'
!           istatus = 0
!           return

!       endif

        istatus = 1
        return
        end

        function zcoord_of_logpressure(pres_pa)

!       This routine interpolates between LAPS levels in logp space

        implicit real*4 (a-z)

        include 'lapsparms.cmn'

        logical ltest_vertical_grid

        integer istatus

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        if(ltest_vertical_grid('HEIGHT'))then

        elseif(ltest_vertical_grid('PRESSURE'))then
            rz = (PRESSURE_0_L - pres_pa) / PRESSURE_INTERVAL_L

            pres_low  = pressure_of_level(int(rz))
            pres_high = pressure_of_level(int(rz)+1)

            frac = log(pres_pa/pres_low) / log(pres_high/pres_low)

            zcoord_of_logpressure = int(rz) + frac

        else
            write(6,*)' Error, vertical grid not supported,'
     1               ,' this routine supports PRESSURE or HEIGHT'
            istatus = 0
            return

        endif

        istatus = 1
        return
        end


        subroutine   uvgrid_to_radar(u_grid,
     1                       v_grid,
     1                       t_radar,
     1                       r_radar,
     1                       azimuth,
     1                       longitude)

        real longitude

        call   uvgrid_to_disptrue      (u_grid,
     1                          v_grid,
     1                          di_true,
     1                          speed,
     1                          longitude)

        call disptrue_to_radar( di_true,
     1                  speed,
     1                  t_radar,
     1                  r_radar,
     1                  azimuth)

        return
        end


        subroutine   uvtrue_to_radar(u_true,
     1                       v_true,
     1                       t_radar,
     1                       r_radar,
     1                       azimuth)

        call   uv_to_disp(u_true,
     1            v_true,
     1            di_true,
     1            speed)

        call disptrue_to_radar( di_true,
     1                  speed,
     1                  t_radar,
     1                  r_radar,
     1                  azimuth)

        return
        end


        subroutine   uvgrid_to_disptrue(u_grid,
     1                          v_grid,
     1                          di_true,
     1                          speed,
     1                          longitude)

        real longitude

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        speed = sqrt( u_grid **2 + v_grid **2 )

        if(speed .gt. 0)then
            di_grid = atan3d(-u_grid,-v_grid)
            di_true = di_grid + projrot_laps(longitude)
            di_true = mod(di_true+360.,360.)
        else
            di_true = 0.
        endif

        return
        end



        subroutine   disptrue_to_radar(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)

        include 'trigd.inc'

!       real longitude

        angle_diff = di_true - azimuth
        t_radar = speed * (-sind(angle_diff))
        r_radar = speed * (-cosd(angle_diff))

        return
        end


        subroutine   radar_to_uvgrid(t_radar,
     1                       r_radar,
     1                       u_grid,
     1                       v_grid,
     1                       azimuth,
     1                       longitude)

        real longitude

        call radar_to_disptrue( di_true,
     1                  speed,
     1                  t_radar,
     1                  r_radar,
     1                  azimuth)

        call   disptrue_to_uvgrid      (di_true,
     1                          speed,
     1                          u_grid,
     1                          v_grid,
     1                          longitude)


        return
        end

        subroutine   radar_to_uvtrue(t_radar,
     1                       r_radar,
     1                       u_true,
     1                       v_true,
     1                       azimuth)

        call radar_to_disptrue( di_true,
     1                  speed,
     1                  t_radar,
     1                  r_radar,
     1                  azimuth)

        call   disp_to_uv      (di_true,
     1                  speed,
     1                  u_true,
     1                  v_true)


        return
        end

        subroutine   radar_to_disptrue(di_true,
     1                         speed,
     1                         t_radar,
     1                         r_radar,
     1                         azimuth)


!       real longitude

        speed = sqrt(t_radar**2 + r_radar**2)

        if(speed .gt. 0)then
            di_radar = atan3d(-t_radar,-r_radar)
            di_true = mod(di_radar + azimuth + 360.,360.)
        else
            di_true = 0.
        endif

        return
        end

        subroutine   disptrue_to_uvgrid(di_true,
     1                          speed,
     1                          u_grid,
     1                          v_grid,
     1                          longitude)

        real longitude

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        di_grid = di_true - projrot_laps(longitude)

        call         disp_to_uv(di_grid,
     1                  speed,
     1                  u_grid,
     1                  v_grid)

        return
        end

        subroutine   uvtrue_to_uvgrid(  u_true,
     1                          v_true,
     1                          u_grid,
     1                          v_grid,
     1                          longitude)

        real longitude

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        angle = projrot_laps(longitude)

        call         rotate_vec(u_true,
     1                  v_true,
     1                  u_grid,
     1                  v_grid,
     1                  angle)

        return
        end

        subroutine   uvgrid_to_uvtrue(  u_grid,
     1                          v_grid,
     1                          u_true,
     1                          v_true,
     1                          longitude)

        real longitude

        call get_laps_config('nest7grid',istatus)

        if(istatus .ne. 1)then
            write(6,*)' ERROR, get_laps_config not successfully called'       
            stop
        endif

        angle = -projrot_laps(longitude)


        call         rotate_vec(u_grid,
     1                  v_grid,
     1                  u_true,
     1                  v_true,
     1                  angle)

        return
        end



      subroutine rotate_vec(u1,v1,u2,v2,angle)
      include 'trigd.inc'
      u2 =  u1 * cosd(angle) + v1 * sind(angle)
      v2 = -u1 * sind(angle) + v1 * cosd(angle)

      return
      end


        subroutine   disp_to_uv(dir,
     1                  speed,
     1                  u,
     1                  v)
        include 'trigd.inc'
        u  = - sind(dir) * speed
        v  = - cosd(dir) * speed

        return
        end


        subroutine   uv_to_disp(u,
     1                  v,
     1                  dir,
     1                  speed)

        speed = sqrt( u**2 + v**2 )

        if(speed .gt. 0)then
            dir = atan3d(-u,-v)
!           dir = mod(dir+360.,360.)
        else
            dir = 0.
        endif

        return
        end

        function k_to_f(x)
        real*4 k_to_f
        k_to_f = (x - 273.15) * 1.8 + 32.
        return
        end

        function f_to_k(x)
        real*4 f_to_k
        f_to_k = (x - 32.) / 1.8 + 273.15
        return
        end

c
c

c
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

        subroutine conv_meters_to_inches(data_in, numarr,prodno,
     1  imax, jmax)
c   subroutine to convert accumulations meters to inches
        REAL*4 data_in(imax, jmax, numarr)
        Integer Prodno
c   begin
        do j=1,jmax
        do i=1,imax
          data_in(i,j,prodno) = data_in(i,j,prodno) / 0.0254
        enddo !i
        enddo !j
c
        return
        end
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

