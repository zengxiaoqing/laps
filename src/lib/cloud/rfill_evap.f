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
        subroutine rfill_evap(ref_3d,ni,nj,nk,l_low_fill,l_high_fill
     1  ,lat,lon,topo,rlat_radar,rlon_radar,rheight_radar
     1  ,temp_3d,rh_3d_pct,cldpcp_type_3d,istatus)

!       This Routine Fills in the gaps in a 3D radar volume that has been
!       stored as a sparse array. A linear interpolation in the vertical
!       direction is used. An additional option 'l_low_fill' can be set to
!       .true.. This will fill in low level echoes that might have been
!       missed because the grid points are below the radar horizon or too
!       close to the terrain. This is done if radar echo is detected not
!       too far above such grid points.

!       Steve Albers            1992          Original rfill subroutine
!       Steve Albers            1993          Add low level evaporation

        include 'lapsparms.for'

!       ni,nj,nk are input LAPS grid dimensions
!       rlat_radar,rlon_radar,rheight_radar are input radar coordinates

        real*4 ref_3d(ni,nj,nk)                  ! Input/Output 3D reflctvy grid
        real*4 temp_3d(ni,nj,nk)                 ! Input 3D temp grid
        real*4 rh_3d_pct(ni,nj,nk)               ! Input 3D rh grid
!       real*4 heights_1d(nk)                    ! Input
        integer*2 cldpcp_type_3d(ni,nj,nk)       ! Input 3D pcp type grid
        real*4 lat(ni,nj),lon(ni,nj),topo(ni,nj) ! Input 2D grids

        integer*4 isum_ref_2d(NX_L_MAX,NY_L_MAX) ! Local array
        real*4 heights_1d(NZ_L_MAX)              ! Local array

        logical l_low_fill,l_high_fill,l_test

        write(6,*)' Interpolating vertically through gaps with evaporati
     1on'
        write(6,*)'i,j,k_upper,dbz_upper,dbz_lower'
     1           ,',rh_3d_pct(i,j,k_upper),r_upper'
     1           ,',delta_t,r_lower,ptype'

        n_low_fill = 0

        isum_test = nint(ref_base) * nk

        do j = 1,nj
        do i = 1,ni
            isum_ref_2d(i,j) = 0
        enddo
        enddo

        do k = 1,nk
        do j = 1,nj
        do i = 1,ni
            isum_ref_2d(i,j) = isum_ref_2d(i,j) + nint(ref_3d(i,j,k))
        enddo
        enddo
        enddo

        do k = 1,nk
            heights_1d(k) = psatoz(pressure_of_level(k)/100.)
        enddo ! k

        if(l_low_fill)then
          if(lat(1,1) .eq. 0. .or. lon(1,1) .eq. 0. .or. topo(1,1) .eq. 
     10.)then
            write(6,*)' Error in RFILL, lat/lon/topo has zero value'
            return
          endif
        endif

        do j = 1,nj
c       write(6,*)' Doing Column ',j

        do i = 1,ni

          if(isum_ref_2d(i,j) .ne. isum_test)then ! Test for presence of echo

            if(l_high_fill)then ! Fill in between high level echoes
                k = 1

!               Test for presence of top
15              l_test = .false.
                ref_below = ref_3d(i,j,k)

                do while (k .lt. nk .and. (.not. l_test))
                    ref_above = ref_3d(i,j,k+1)

                    if(ref_below .gt. ref_base
     1         .and. ref_above .eq. ref_base)then
                        k_top = k
                        l_test = .true.
                    endif

                    k = k + 1
                    ref_below = ref_above

                enddo

                if(.not. l_test)goto100 ! No Top exists

!               Top exists, Search for next bottom
                l_test = .false.
                ref_below = ref_3d(i,j,k)

                do while(k .lt. nk .and. (.not. l_test))
                    k = k + 1
                    ref_above = ref_3d(i,j,k)

                    if(ref_above .gt. ref_base
     1         .and. ref_below .eq. ref_base)then
                        k_bottom = k
                        l_test = .true.
                    endif

                    ref_below = ref_above

                enddo

                if(.not. l_test)goto100 ! No Bottom exists

!               Fill in gap if it exists and is small enough
!               if(.true.)then
                if(k_bottom .le. k_top+5)then ! Fill gaps smaller than const-1
c                   write(6,*)' Filling gap between',i,j,k_bottom,k_top
c                   write(6,101)(nint(max(ref_3d(i,j,kwrt),ref_base)),kwrt=1,nk)
101                 format(1x,17i4)

                    do k_bet = k_top+1,k_bottom-1
                        frac = float(k_bet-k_bottom)/float(k_top-k_botto
     1m)
                        ref_3d(i,j,k_bet) = ref_3d(i,j,k_bottom)*(1.-fra
     1c)
     1                            + ref_3d(i,j,k_top)*(frac)
                    enddo ! k

c                   write(6,101)(nint(max(ref_3d(i,j,kwrt),ref_base)),kwrt=1,nk)

                endif

                if(k .le. nk-2)then ! Look for another gap
                    goto15
                endif

            endif ! l_high_fill

100         if(l_low_fill)then ! Fill below low level echoes

!               Search for bottom of detected echo
                k_bottom = 0
                l_test = .false.

                k = 2

                ref_below = ref_3d(i,j,1)

                do while((.not. l_test) .and. k .le. nk)
                    ref_above = ref_3d(i,j,k)

                    if(ref_above .gt. ref_base
     1         .and. ref_below .eq. ref_base)then
                        k_bottom = k
                        dbz_bottom = ref_above
                        l_test = .true.
                    endif

                    k = k + 1
                    ref_below = ref_above

                enddo

                istatus = 1

                k_topo = max(int(height_to_zcoord(topo(i,j),istatus)),1)

                if(istatus .ne. 1)then
                    write(6,*)' ERROR return in rfill'
                    return
                endif

                if(k_bottom .gt. k_topo)then ! Echo above ground

                    height_bottom = height_of_level(k_bottom)

                    call latlon_to_radar(lat(i,j),lon(i,j),height_bottom
     1                  ,azimuth,slant_range,elev_bottom
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                    call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                  ,azimuth,slant_range,elev_topo
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                    if(k_bottom .le. k_topo+2 ! Echo base near ground
!       1               .or. elev_bottom .lt. (elev_topo + 1.5) ! Echo below radar horizon
     1          .or. elev_bottom .lt. 1.5 ! Echo base near radar horizon
     1                                                  )then ! Fill in

                        do k = k_bottom-1,k_topo,-1
!                           ref_3d(i,j,k) = ref_3d(i,j,k_bottom)
                            call modify_dbz_level(i,j,k+1,k,ref_3d
     1                           ,temp_3d,rh_3d_pct,cldpcp_type_3d
     1                           ,heights_1d,ni,nj,nk,ref_base,istatus)        
                            if(istatus .ne. 1)then
                                return
                            endif
                        enddo ! k

                        write(6,*)

!                       write(6,211)i,j,k_topo,k_bottom,elev_topo,elev_bottom
!       1                       ,slant_range,ref_3d(i,j,k_bottom)
211                     format(' low fill ',4i5,2f6.1,2f8.0)

                        n_low_fill = n_low_fill + 1

                    endif ! Fill in from ground to bottom of echo

                endif ! Bottom of radar echo above ground

            endif ! l_low_fill

          endif ! echo present

        enddo ! i
        enddo ! j

        write(6,*)' n_low_fill = ',n_low_fill

        return
        end

        subroutine modify_dbz_level(i,j,k_upper,k_lower,ref_3d,temp_3d
     1  ,rh_3d_pct,cldpcp_type_3d,heights_1d,ni,nj,nk,ref_base,istatus)       

        real*4 ref_3d(ni,nj,nk)                  ! Input/Output 3D reflctvy grid
        real*4 temp_3d(ni,nj,nk)                 ! Input 3D temp grid
        real*4 rh_3d_pct(ni,nj,nk)               ! Input 3D rh grid
        real*4 heights_1d(nk)                    ! Input
        integer*2 cldpcp_type_3d(ni,nj,nk)       ! Input 3D pcp type grid

        real*4 maxrate_R2V,maxrate_S2V,maxrate_I2V
!       parameter (maxrate_R2V = 0.83e-5) ! s**-1     (Original Schultz value)
!       parameter (maxrate_S2V = 1.67e-5) ! s**-1     (Original Schultz value)
!       parameter (maxrate_I2V = 3.33e-6) ! s**-1     (Original Schultz value)

        parameter (maxrate_R2V = 0.83e-7) ! s**-1     (New Empirical value)
        parameter (maxrate_S2V = 1.67e-7) ! s**-1     (New Empirical value)
        parameter (maxrate_I2V = 3.33e-8) ! s**-1     (New Empirical value)

        integer*4 i_debug
        data i_debug /0/
        save i_debug

        i_debug = i_debug + 1

!       This routine computes the change in radar reflectivity from one
!       LAPS level to the next lower level due to evaporation.

!       Find precip rate at the higher level
        dbz_upper = ref_3d(i,j,k_upper)

        if(dbz_upper .gt. ref_base)then
            pres_upper = pressure_of_level(k_upper)
            ipcp_type_upper = cldpcp_type_3d(i,j,k_upper) / 16   ! Pull out precip type
            pcp_rate_upper = dbz_to_rate(dbz_upper,temp_3d(i,j,k_upper)
     1  ,pres_upper,ipcp_type_upper,istatus)

            thickness = heights_1d(k_upper) - heights_1d(k_lower)

            if(ipcp_type_upper .ne. 0)then
                call cpt_fall_velocity(ipcp_type_upper,temp_3d(i,j,k_upp
     1er)
     1                                  ,pres_upper,fall_velocity)
!               Get concentration in g/m**3
                call cpt_concentration(pcp_rate_upper,fall_velocity
     1                                                  ,pcpcnc_upper)
            else
                write(6,*)' Error in modify_dbz_level, no precip type'
     1                          ,i,j,k_upper,dbz_upper
                istatus = 0
                return
            endif

            delta_t = thickness / fall_velocity

            if(ipcp_type_upper .eq. 1 .or. ipcp_type_upper .eq. 3)then ! R,ZR

!               Calculate saturation mixing ratio
                rvsatliq = W_fast(temp_3d(i,j,k_upper)-273.15
     1                         ,pressure_of_level(k)/100.     )       * 
     1.001

                rv = rvsatliq * rh_3d_pct(i,j,k_upper)/100.   ! Ambient mixing ratio

                call ConvR2V (maxrate_R2V, rv, rvsatliq, rate_evap)

            elseif(ipcp_type_upper .eq. 2)then                         ! S

                rvsatliq = W_fast(temp_3d(i,j,k_upper)-273.15
     1                         ,pressure_of_level(k)/100.     )       * 
     1.001

                if(temp_3d(i,j,k_upper) .le. 273.15)then
                    rvsatice = Wice_fast(temp_3d(i,j,k_upper)-273.15
     1                         ,pressure_of_level(k)/100.     )       * 
     1.001

                else
                    rvsatice = rvsatliq

                endif

                rv = rvsatliq * rh_3d_pct(i,j,k_upper)/100.   ! Ambient mixing ratio

                call ConvS2V (maxrate_S2V, rv, rvsatice, rate_evap)

            elseif(ipcp_type_upper .eq. 4 .or. ipcp_type_upper .eq. 5)th
     1en ! IP,A

                rvsatliq = W_fast(temp_3d(i,j,k_upper)-273.15
     1                         ,pressure_of_level(k)/100.     )       * 
     1.001

                if(temp_3d(i,j,k_upper) .le. 273.15)then
                    rvsatice = Wice_fast(temp_3d(i,j,k_upper)-273.15
     1                         ,pressure_of_level(k)/100.     )       * 
     1.001

                else
                    rvsatice = rvsatliq

                endif

                rv = rvsatliq * rh_3d_pct(i,j,k_upper)/100.   ! Ambient mixing ratio

                call ConvI2V (maxrate_I2V, rv, rvsatice, rate_evap)

            endif


            evaporation = rate_evap * delta_t

!           Calculate Mixing Ratio at the lower level
            r_upper = pcpcnc_upper * .001 ! (g/m**3 to kg/kg)
            r_lower = max(r_upper - evaporation,0.)

            pcpcnc_lower = pcpcnc_upper * (r_lower / r_upper)

!           In absence of evaporation, pcp_rate is conserved through the level
            pcp_rate_lower = (pcpcnc_lower / 1e6) * fall_velocity

            pres_lower = pressure_of_level(k_lower)
            ipcp_type_lower = cldpcp_type_3d(i,j,k_lower) / 16   ! Pull out precip type
            dbz_lower = rate_to_dbz(pcp_rate_lower,temp_3d(i,j,k_lower)
     1                  ,pres_lower,ref_base,ipcp_type_lower,istatus)

            if(i_debug .lt. 200)write(6,1)i,j,k_upper,dbz_upper
     1          ,dbz_lower,rh_3d_pct(i,j,k_upper),r_upper
     1          ,nint(delta_t),r_lower,ipcp_type_upper
1           format(1x,3i4,2f7.1,f6.0,f10.7,i6,f10.7,i2)

        else
            dbz_lower = dbz_upper
            if(i_debug .lt. 200)write(6,1)i,j,k_upper,dbz_upper,dbz_lowe
     1r

        endif

        ref_3d(i,j,k_lower) = dbz_lower

        istatus = 1
        return
        end

        function dbz_to_rate(dbz,t,p,itype,istatus)

        real*4 a,b,rate_max
        parameter (a = 200.)        ! Z-R relationship
        parameter (b = 1.6)         ! Z-R relationship
        parameter (rate_max = 1000.0) ! Currently disabled

        aterm = alog10(1./a)
        bterm = 1./b
        cterm = .001 / 3600.  ! (M/S) / (MM/HR)

        r_mm_hr = 10.**(bterm*(aterm + dBZ/10.))

        dbz_to_rate = (min(r_mm_hr,rate_max) * cterm)

        if(p .ne. 101300. .or. t .ne. 273.15)then
            density_norm = (p / 101300.) * (273.15 / t)
            sqrt_density_norm = sqrt(density_norm)
        else
            sqrt_density_norm = 1.
        endif

        dbz_to_rate = dbz_to_rate / sqrt_density_norm

        return
        end

        function rate_to_dbz(rate,t,p,ref_base,itype,istatus)

        real*4 a,b
        parameter (a = 200.)        ! Z-R relationship
        parameter (b = 1.6)         ! Z-R relationship

        if(rate .le. 0)then
            rate_to_dbz = ref_base
            return
        endif

        cterm = .001 / 3600.  ! (M/S) / (MM/HR)
        cterm_inv = 1. / cterm

        if(p .ne. 101300. .or. t .ne. 273.15)then
            density_norm = (p / 101300.) * (273.15 / t)
            sqrt_density_norm = sqrt(density_norm)
        else
            sqrt_density_norm = 1.
        endif

        rate = rate * sqrt_density_norm

        z = a * (rate*cterm_inv)**b
        rate_to_dbz = 10. * alog10(z)

        rate_to_dbz = max(rate_to_dbz,ref_base)

        return
        end
