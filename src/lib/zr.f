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
        subroutine zr(z_2d_in
!       1       ,temp_sfc_k,td_sfc_k,pres_sta_pa,tw_sfc_k
     1                                          ,ni,nj,r_2d_out)

!           1991                                        Steve Albers, J. Smart
!       Apr 1994 rate_max parameter disabled            SA
!       Jan 1998 Remove lapsparms.inc                   SA

        real z_2d_in(ni,nj)
!       real temp_sfc_k(ni,nj)
!       real td_sfc_k(ni,nj)
!       real pres_sta_pa(ni,nj)
!       real tw_sfc_k(ni,nj)
        real r_2d_out(ni,nj)

        real a,b,rate_max
        parameter (a = 200.)        ! Z-R relationship
        parameter (b = 1.6)         ! Z-R relationship
!       parameter (rate_max = 10.0) ! mm/hr; equiv to 39 dBZ with Z=200R**1.6
        parameter (rate_max = 1000.0) ! disabled

        write(6,*)
     1   ' Converting from 2D Z to Rainfall Rate field'

        call get_ref_base_useable(ref_base_useable,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting ref_base_useable in zr'
            stop
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error getting r_missing_data in zr'
            stop
        endif

        n_leq_pts = 0
!       n_warm_pts = 0

        aterm = alog10(1./a)
        bterm = 1./b
        cterm = .001 / 3600  ! (M/S) / (MM/HR)

        do j = 1,nj
        do i = 1,ni
            dbz = z_2d_in(i,j)

            if(dbz .eq. r_missing_data)then
                r_2d_out(i,j) = r_missing_data

            elseif(dbz .lt. ref_base_useable)then
!               r_2d_out(i,j) = +1e-30
                r_2d_out(i,j) = 0.

            else
                n_leq_pts = n_leq_pts + 1

              ! Generate R (mm/hr) in Z=a*R**b
                r_mm_hr = 10.**(bterm*(aterm + dBZ/10.))
                r_2d_out(i,j) = min(r_mm_hr,rate_max) * cterm

!               if(tw_sfc_k(i,j) .gt. 273.65)then ! Wet Bulb > 0.5 Deg C
!                   r_2d_out(i,j) = 0.
!                   n_warm_pts = n_warm_pts + 1
!               endif

            endif

        enddo
        enddo

        write(6,*)' n_leq_pts = ',n_leq_pts

        return
        end
