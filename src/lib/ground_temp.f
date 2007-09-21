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

        subroutine get_ground_temperature(i4time,laps_cycle_time
     1                    ,imax,jmax,lat,lon,r_missing_data
     1                                      ,cvr_snow,t_sfc_k,t_gnd_k)

!       This routine returns ground temp under the assumption of clear skies
!       After this routine is called from clouds, it may be worth filtering
!       the band 8 temps to expand the cloudy areas slightly. The kernel
!       of the filter would be the coldest pixel of the nine neighbors

        character var*3,comment*125,ext*31,units*10


        integer imax,jmax         ! input
        integer i4time            ! input
        integer laps_cycle_time   ! input
        real lat(imax,jmax)       ! input
        real lon(imax,jmax)       ! input
        real r_missing_data       ! input
        real t_sfc_k(imax,jmax)   ! input
        real t_gnd_k(imax,jmax)   ! output
        real cvr_snow(imax,jmax)  ! local

        var = 'SC'
        ext = 'lm2'
        call get_laps_2d(i4time-laps_cycle_time,ext,var,units,comment
     1                            ,imax,jmax,cvr_snow,istat_cvr_snow)

        if(istat_cvr_snow .ne. 1)then    ! Try using snow cover field
            write(6,*)' Error reading snow cover'
            do i = 1,imax
            do j = 1,jmax
                cvr_snow(i,j) = r_missing_data
            enddo
            enddo
        endif



!       Calculate solar altitude and ground temperature
        do i = 1,imax
        do j = 1,jmax
            call solar_position(lat(i,j),lon(i,j),i4time,solar_alt
     1                                     ,solar_dec,solar_ha)

!           Calculate ground temperature (for now equate to sfc air temp)
!           Note that function t_ground_k needs to be moved to the lapslib
            t_gnd_k(i,j) = t_ground_k(t_sfc_k(i,j)
     1                          ,solar_alt,solar_ha
     1              ,solar_dec,lat(i,j),cvr_snow(i,j),r_missing_data
     1                          ,i,j,imax,jmax)


        enddo
        enddo


        return
        end

        function t_ground_k(t_sfc_k,solar_alt,solar_ha
     1    ,solar_dec,rlat,cvr_snow,r_missing_data,i,j,ni,nj)

        real high_alt,low_alt,corr_low,corr_high,corr_ramp

        solar_transit_alt = 90. - abs(rlat - solar_dec)

        if(solar_ha .lt. 0.)then ! morning
            high_alt = solar_transit_alt * .66
            low_alt  = solar_transit_alt * .36
!           low_alt  = high_alt - 11.
        else                     ! afternoon
            high_alt = solar_transit_alt * .79
            low_alt  = solar_transit_alt * .49
!           low_alt  = high_alt - 11.
        endif

        corr_high =+4.
        corr_low = -4.
        corr_ramp = (corr_high - corr_low) / (high_alt - low_alt)

!       Warmer at day, colder at night
        if(solar_alt .lt. low_alt)then
            corr = corr_low

        elseif(solar_alt .gt. high_alt)then
            corr = corr_high

        else ! ramp the function
            corr = corr_low + (solar_alt - low_alt) * corr_ramp

        endif

        if(i .eq. ni/2 .and. j .eq. nj/2)then
            write(6,*)' solar_ha,solar_alt,solar_transit_alt,ratio = '
     1             ,solar_ha,solar_alt,solar_transit_alt
     1             ,solar_alt/solar_transit_alt
            write(6,*)' high/low alt, corr = ',high_alt,low_alt,corr
        endif

        corr_negonly = min(corr,0.0)            ! Only colder at night

        if(cvr_snow .ne. r_missing_data)then
            t_ground_fullcorr = t_sfc_k + corr  ! Warmer at day, colder at night
            t_snow_k = t_sfc_k + corr_negonly   ! Only colder at night
            t_snow_k = min(t_snow_k,273.15)     ! Snow can't be above 0C
            t_ground_k = t_ground_fullcorr * (1.-cvr_snow) + t_snow_k *
     1cvr_snow
        else
            t_ground_k = t_sfc_k + corr_negonly ! Only colder at night
        endif

        return
        end

























