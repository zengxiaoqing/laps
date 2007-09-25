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
        subroutine rfill_evap(ref_3d,ni,nj,nk,cloud_base
     1  ,lat,lon,topo,mode_evap
     1  ,temp_3d,rh_3d_pct,cldpcp_type_3d,heights_3d,istatus,ref_base)       

!       This Routine works with reflectivity read in from the LPS file. This
!       is originally either 2D or 3D and may have been further processed
!       within the cloud analysis.

!       ni,nj,nk are input LAPS grid dimensions

        real ref_3d(ni,nj,nk)                  ! I/O 3D reflctvy grid
        real temp_3d(ni,nj,nk)                 ! I
        real rh_3d_pct(ni,nj,nk)               ! I
        real heights_3d(ni,nj,nk)              ! I
        integer cldpcp_type_3d(ni,nj,nk)         ! I
        real lat(ni,nj),lon(ni,nj),topo(ni,nj) ! I
        real cloud_base(ni,nj)                 ! I

        real heights_1d(nk)                    ! L

        logical l_low_fill,l_high_fill,l_test

        write(6,*)' subroutine rfill_evap'

        n_low_attenuate = 0

        isum_test = nint(ref_base) * nk

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        do j = 1,nj
c       write(6,*)' Doing Column ',j

        do i = 1,ni

!           Test for vertical existance/variance to characterize reflectivity
            reference = r_missing_data
            iref_dim = 0

            do k = 1,nk
                if(ref_3d(i,j,k) .ne. ref_base 
     1       .and. ref_3d(i,j,k) .ne. r_missing_data)then
                    if(reference .eq. r_missing_data)then ! initial ref
                        iref_dim = 2
                        reference = ref_3d(i,j,k)
                    else                                  ! subsequent ref
                        if(ref_3d(i,j,k) .ne. reference)then
                            iref_dim = 3
                        endif
                    endif
                endif
            enddo ! k

            if(iref_dim .eq. 2 .and. mode_evap .eq. 3)then 
                goto 900                              ! don't evaporate 2d data
            endif

!           Determine "Radar Base" to be Cloud Base for 2dref points
!           and Radar Horizon for 3dref points.

            if(iref_dim .eq. 2)then 
                k_cloud_base = height_to_zcoord2(cloud_base(i,j)
     1                               ,heights_3d,ni,nj,nk,i,j,istatus)

                if(cloud_base(i,j) .ne. r_missing_data .and.
     1             k_cloud_base    .ne. r_missing_data .and. 
     1             istatus         .eq. 1                    )then     

                    k_bottom = k_cloud_base

                else
                    write(6,*)' Warning, radar echo without cloud base'
                    k_bottom = 0
 
                endif

            elseif(iref_dim .eq. 3)then

!               Search for bottom of detected echo
                k_bottom = 0
                l_test = .false.

                k = 2

                ref_below = ref_3d(i,j,1)

                do while((.not. l_test) .and. k .le. nk)
                    ref_above = ref_3d(i,j,k)

                    if(ref_above .gt. ref_base
     1           .and. ref_below .eq. ref_base)then
                        k_bottom = k
                        dbz_bottom = ref_above
                        l_test = .true.
                    endif

                    k = k + 1
                    ref_below = ref_above

                enddo

            else
                k_bottom = 0

            endif

            if(iref_dim .gt. 0)then    ! We have an echo
                istatus = 1

                k_topo = max(int(height_to_zcoord2(topo(i,j),heights_3d       
     1                           ,ni,nj,nk,i,j,istatus)),1)

                if(istatus .ne. 1)then
                    write(6,*)' ERROR return in rfill'
                    return
                endif

                if(k_bottom .gt. k_topo)then ! Echo above ground
                    height_bottom = heights_3d(i,j,k_bottom)

                    call latlon_to_radar(lat(i,j),lon(i,j),height_bottom
     1                  ,azimuth,slant_range,elev_bottom
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                    call latlon_to_radar(lat(i,j),lon(i,j),topo(i,j)
     1                  ,azimuth,slant_range,elev_topo
     1                  ,rlat_radar,rlon_radar,rheight_radar)

                    if(k_bottom .le. k_topo+2 ! Echo base near ground
!       1         .or. elev_bottom .lt. (elev_topo + 1.5) ! Echo below radar horizon
     1            .or. elev_bottom .lt. 1.5 ! Echo base near radar horizon
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

                        n_low_attenuate = n_low_attenuate + 1

                    endif ! Fill in from ground to bottom of echo

                endif ! Bottom of radar echo above ground

            endif ! l_low_fill

 900        continue ! Skipped this grid point

        enddo ! i
        enddo ! j

        write(6,*)' n_low_attenuate = ',n_low_attenuate

        return
        end

        subroutine modify_dbz_level(i,j,k_upper,k_lower,ref_3d,temp_3d
     1  ,rh_3d_pct,cldpcp_type_3d,heights_1d,ni,nj,nk,ref_base,istatus)       

        real ref_3d(ni,nj,nk)                  ! Input/Output 3D reflctvy grid
        real temp_3d(ni,nj,nk)                 ! Input 3D temp grid
        real rh_3d_pct(ni,nj,nk)               ! Input 3D rh grid
        real heights_1d(nk)                    ! Input
        integer cldpcp_type_3d(ni,nj,nk)       ! Input 3D pcp type grid

        real maxrate_R2V,maxrate_S2V,maxrate_I2V
!       parameter (maxrate_R2V = 0.83e-5) ! s**-1     (Original Schultz value)
!       parameter (maxrate_S2V = 1.67e-5) ! s**-1     (Original Schultz value)
!       parameter (maxrate_I2V = 3.33e-6) ! s**-1     (Original Schultz value)

        parameter (maxrate_R2V = 0.83e-7) ! s**-1     (New Empirical value)
        parameter (maxrate_S2V = 1.67e-7) ! s**-1     (New Empirical value)
        parameter (maxrate_I2V = 3.33e-8) ! s**-1     (New Empirical value)

        integer i_debug
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

        real a,b,rate_max
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

        real a,b
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
