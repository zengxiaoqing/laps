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
        subroutine zs(precip_rate,temp_col_max,ni,nj,s_2d_out)

!           1994    Steve Albers
!       Jan 1998    Steve Albers     Remove lapsparms.inc and other cleanup

        real precip_rate(ni,nj)
        real temp_col_max(ni,nj) ! Deg K
        real s_2d_out(ni,nj)

        call get_ref_base(ref_base,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting ref_base in zs'
            stop
        endif

        n_snow_pts = 0
!       n_warm_pts = 0

        do j = 1,nj
        do i = 1,ni

!           if(dbz .eq. ref_base)then
!               s_2d_out(i,j) = +1e-30

!           else
                ratio = snow_to_rain_ratio(temp_col_max(i,j))
                n_snow_pts = n_snow_pts + 1
                s_2d_out(i,j) = precip_rate(i,j) * ratio

                if(n_snow_pts .eq. (n_snow_pts/200) * 200   )then
                    write(6,*)i,j,temp_col_max(i,j)-273.15, ratio
                endif

!           endif

        enddo
        enddo

        write(6,*)' n_snow_pts = ',n_snow_pts

        return
        end

        function snow_to_rain_ratio(temp_col_max)

!       Note that temp_col_max is in Deg K

        temp_col_c = temp_col_max - 273.15   ! Convert from K to C

        if(temp_col_c .ge. 0.0)then          !  T > 0C, use 10

            snow_to_rain_ratio = 10.

        elseif(temp_col_c .ge. -3.0)then     !  0C > T >  -3C, ramp 10 - 15

            frac = temp_col_c/ (-3.0)
            snow_to_rain_ratio = 10. * (1. - frac) + 15. * frac

        elseif(temp_col_c .ge. -10.0)then    !  -3C > T > -10C, ramp 15 - 25

            frac = (temp_col_c - (-3.0)) / (-7.0)
            snow_to_rain_ratio = 15. * (1. - frac) + 25. * frac

        elseif(temp_col_c .ge. -18.0)then    ! -10C > T > -18C, use 25

            snow_to_rain_ratio = 25.

        elseif(temp_col_c .ge. -22.0)then    ! -18C > T > -22C, ramp 25 - 15

            frac = (temp_col_c - (-18.0)) / (-4.0)
            snow_to_rain_ratio = 25. * (1. - frac) + 15. * frac

        else                                 !        T < -22C, use 15

            snow_to_rain_ratio = 15.

        endif

        return
        end
