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

!       1994    Steve Albers

        include 'lapsparms.inc'   ! ref_base

        real*4 precip_rate(ni,nj)
        real*4 temp_col_max(ni,nj) ! Deg K
        real*4 s_2d_out(ni,nj)

        n_snow_pts = 0
!       n_warm_pts = 0

        do j = 1,nj
        do i = 1,ni

            if(dbz .eq. ref_base)then
                s_2d_out(i,j) = +1e-30

            else
                ratio = snow_to_rain_ratio(temp_col_max(i,j))
                n_snow_pts = n_snow_pts + 1
                s_2d_out(i,j) = precip_rate(i,j) * ratio

                if(n_snow_pts .eq. (n_snow_pts/20) * 20   )then
                    write(6,*)i,j,temp_col_max(i,j)-273.15, ratio
                endif

            endif

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

        subroutine zs_old(z_2d_in,temp_col_max
!       1       ,temp_sfc_k,td_sfc_k,pres_sta_pa,tw_sfc_k
     1                                          ,ni,nj,s_2d_out)

!       1991    Steve Albers
!       1992    J. Smart        Add temperature bias adjustment

        include 'lapsparms.inc'

        integer*4 iTemp(2501)
        integer*4 index_from_temp(24815:27315)  ! index to acquire adjustment
                                                ! factor from temp.
        integer*4 iT
        real*4 ZS_bias_factor(2501)

        real*4 z_2d_in(ni,nj)
        real*4 temp_col_max(ni,nj)
!       real*4 temp_sfc_k(ni,nj)
!       real*4 td_sfc_k(ni,nj)
!       real*4 pres_sta_pa(ni,nj)
!       real*4 tw_sfc_k(ni,nj)
        real*4 s_2d_out(ni,nj)

        logical compute_bias_factor

        real*4 a,b,rate_max
        parameter (a = 50.)         ! Z S relationship
        parameter (b = 1.8)         ! Z S relationship
        parameter (rate_max = 10.0) ! cm/hr

        common/zs_comm/ compute_bias_factor

C
C Snow depth temperature adjustment factor computation.  Resolution
C is to 0.01 degrees.
C Adjustment is based upon linear ramp in which precip at colder temps is
C increased while at warmer temps it is decreased.  Presently the zero
C adjustment (ie., ZS_bias_factor = 1.0) occurs at -8C max col temperature.
C Temperatures colder than -25C are adjusted by the same amount as -25C.
C Precip type does not allow snow accumulation if max col temperature > 0C.
C The temperature is in K.

        if(compute_bias_factor)then
          icnt=0
          do iT=0,2500
            icnt=icnt+1
            iTemp(icnt)=27315-iT
            ZS_bias_factor(icnt)=(27848.33-FLOAT(iTemp(icnt)))/1333.33
            index_from_temp(iTemp(icnt))=icnt
          end do
          compute_bias_factor = .false.
        end if

        write(6,*)
     1   ' Converting from 2D Z to Snowfall Rate field'


        n_snow_pts = 0
!       n_warm_pts = 0

        aterm = alog10(1./a)
        bterm = 1./b
        cterm = .01 / 3600  ! (M/S) / (CM/HR)

        do j = 1,nj
        do i = 1,ni
            dbz = z_2d_in(i,j)

            if(dbz .eq. ref_base)then
                s_2d_out(i,j) = +1e-30

            else
                n_snow_pts = n_snow_pts + 1

              ! Generate S (cm/hr) in Z=a*S**b
                s_cm_hr = 10.**(bterm*(aterm + dBZ/10.))
                s_2d_out(i,j) = min(s_cm_hr,rate_max) * cterm

C adjust snow rate by max col temperature
C
                if((temp_col_max(i,j).ge.248.15) .and.
     1     (temp_col_max(i,j).le.273.15))then

                  iT = (temp_col_max(i,j)*100.)
                  factor = ZS_bias_factor(index_from_temp(iT))
                  s_2d_out(i,j) = s_2d_out(i,j)*factor

                else if(temp_col_max(i,j).gt.273.15)then

                  factor = ZS_bias_factor(1)
                  s_2d_out(i,j) = s_2d_out(i,j)*factor

                else    !temp < -25C

                  s_2d_out(i,j) = s_2d_out(i,j)*ZS_bias_factor(2501)

                end if

!               if(tw_sfc_k(i,j) .gt. 273.65)then ! Wet Bulb > 0.5 Deg C
!                   s_2d_out(i,j) = 0.
!                   n_warm_pts = n_warm_pts + 1
!               endif

            endif

        enddo
        enddo

        write(6,*)' n_snow_pts = ',n_snow_pts

        return
        end

c        block data
c        logical compute_bias_factor
c        common/zs_comm/ compute_bias_factor
c        data compute_bias_factor /.true./
c        end

