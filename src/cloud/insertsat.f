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
c
c
        subroutine insert_sat(i4time,cldcv,
     1  cldcv_sao,cld_hts,rlat,rlon,                                     ! I
     1  pct_req_lvd_s8a,default_clear_cover,                             ! I
     1  tb8_k,istat_tb8,                                                 ! I
     1  sst_k,istat_sst,                                                 ! I
     1  istat_39_a, l_use_39,                                            ! I
     1  istat_39_add_a,                                                  ! O
     1  tb8_cold_k,                                                      ! O
     1  grid_spacing_m,surface_sao_buffer,                               ! I
!    1  cloud_frac_vis_a,istat_vis,
     1  solar_alt,solar_ha,solar_dec,                                    ! I
     1  cloud_frac_co2_a,                                                ! O
     1  rlaps_land_frac,topo,heights_3d,temp_3d,t_sfc_k,pres_sfc_pa,     ! I
     1  cvr_snow,imax,jmax,kcld,klaps,r_missing_data,                    ! I
     1  t_gnd_k,                                                         ! O
     1  cldtop_m_co2,cldtop_m_tb8,cldtop_m,                              ! O
     1  istatus)                                                         ! O
c
c*************************************************************************
c
c       Routine to process satellite data for clouds and to modify the
c       3d cloud cover array.
c       Currently Band 8 (11.2 micron) brightness temps are used with a
c       dummy call for the CO2 slicing method.
c
c       1993        Steve Albers
c       1995 Dec 12 Steve Albers   QC check added prior to call of subroutine
c                                  rad_to_temp
c       1996 Sep    Steve Albers   Fix QC check comparing cloud heights
c                                  to heights_3d(i,j,klaps)
!       1997 Aug 01 Ken Dritz      Changed NX_L to imax, NY_L to jmax, and
!                                  NZ_L to klaps
!       1997 Aug 01 Ken Dritz      Added r_missing_data as dummy argument
!       1997 Aug 01 Ken Dritz      Removed include of lapsparms.for
c
c*************************************************************************
c
!       Prevents clearing out using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real*4 surface_sao_buffer

!       Prevents adding cloud using satellite (hence letting SAOs dominate)
!       below this altitude (M AGL)
        real*4 surface_ir_buffer
        parameter (surface_ir_buffer = 5000.) ! 3000.

!       Default Thickness of Clouds Inserted by Satellite
        real*4 thk_def
        parameter (thk_def            = 1500.)

!       Cloud cover threshold in evaluating presence of SAO/PIREP layers
        real*4 thr_sao_cvr
        parameter (thr_sao_cvr = 0.1)

!       Threshold for IR cloud detection (SFC temp - IR temp)
        real*4 thresh_ir_diff1
        parameter (thresh_ir_diff1 = 8.)

!       Second threshold for IR cloud detection (SFC temp - IR temp)
        real*4 thresh_ir_diff2
        parameter (thresh_ir_diff2 = 21.)

        character*3 lvd_ext
        data lvd_ext /'lvd'/

!       Input/Output
        real*4 cldcv(imax,jmax,kcld)       ! 3D Cloud cover array

        include 'laps_cloud.inc'

!       Input
!       real*4 cld_hts(kcloud)
        real*4 cldcv_sao(imax,jmax,kcld)       ! 3D Cloud cover array
        real*4 rlat(imax,jmax),rlon(imax,jmax)
        real*4 tb8_k(imax,jmax)
        real*4 sst_k(imax,jmax)
        real*4 tb8_cold_k(imax,jmax)
        real*4 topo(imax,jmax)
        real*4 rlaps_land_frac(imax,jmax)
!       real*4 cloud_frac_vis_a(imax,jmax)
        real*4 solar_alt(imax,jmax)
        real*4 solar_ha(imax,jmax)
        real*4 temp_3d(imax,jmax,klaps)
        real*4 t_sfc_k(imax,jmax)
        real*4 cvr_snow(imax,jmax)
        real*4 pres_sfc_pa(imax,jmax)
        real*4 heights_3d(imax,jmax,klaps)

        integer*4 istat_39_a(imax,jmax)
        integer*4 istat_39_add_a(imax,jmax)

!       Output
        real*4 t_gnd_k(imax,jmax)
        real*4 cldtop_m(imax,jmax)
        real*4 cldtop_m_co2(imax,jmax)
        real*4 cldtop_m_tb8(imax,jmax)
        real*4 cloud_frac_co2_a(imax,jmax)

!       Local

        real*4 k_terrain(imax,jmax)
        real*4 zcoords_1d(klaps)
        real*4 cldcv_1d(kcloud)
        real*4 laps_p(klaps)

        character*31 ext
        character var*3,comment*125,units*10

        logical  l_tb8,l_cloud_present,l_use_39
        logical l_co2
        data l_co2 /.false./ ! Attempt to use co2 slicing method?

        real*4 k_to_f

!       Control search box for SAO analyzed data
        integer*4 idelt(3)
        integer*4 jdelt(3)
!       data idelt/-5,0,+5/
!       data jdelt/0,+5,-5/

        integer*4 nidelt,njdelt
        data nidelt/3/,njdelt/3/

        idelt_max = nint(50000. / grid_spacing_m)

        idelt(1) = -idelt_max
        idelt(2) = 0
        idelt(3) = +idelt_max

        jdelt(1) = 0
        jdelt(2) = +idelt_max
        jdelt(3) = -idelt_max

        do k = 1,klaps
            zcoords_1d(k) = zcoord_of_level(k)
            laps_p(k)     = pressure_of_level(k)
        enddo ! k

        if(istat_tb8 .ne. 1)then
            if(pct_req_lvd_s8a .gt. 0.)then
                write(6,*)
     1              ' Error reading tb8_k - return from insert_sat'            
                istatus = 0
                return

            else ! pct_req_lvd_s8a .eq. 0
                write(6,*)' No lvd files, however pct_req_lvd_s8a = 0.'       
                istatus = 1
                return

            endif

        endif

!       Calculate valid percentage
        icount = 0
        do j = 1,jmax
        do i = 1,imax
            if(tb8_k(i,j) .eq. r_missing_data)then
                icount = icount + 1
            endif
        enddo
        enddo

        valid_pct = (1.- float(icount)/float(imax*jmax) ) * 100.

        if(valid_pct .lt. pct_req_lvd_s8a)then
            write(6,51)pct_req_lvd_s8a,valid_pct
 51         format(' WARNING: insufficient LVD S8A data ',2f8.2)
            istatus = 0
            return

        else
            write(6,52)pct_req_lvd_s8a,valid_pct
 52         format(' Pct LVD S8A data required/actual...',2f8.2)
            
        endif

!       Calculate cold filtered temperatures (over ~20km box)
        i_delt = max(1,nint(10000./grid_spacing_m))
        do j = 1,jmax
            jl = max(1   ,j-i_delt)
            jh = min(jmax,j+i_delt)

            do i = 1,imax
                il = max(1   ,i-i_delt)
                ih = min(imax,i+i_delt)
                tb8_cold_k(i,j) = tb8_k(i,j)

                do jj = jl,jh,1
                do ii = il,ih,1
                  if(tb8_k(ii,jj) .ne. r_missing_data)then
                    tb8_cold_k(i,j) = min(tb8_cold_k(i,j),tb8_k(ii,jj))
                  endif
                enddo ! ii
                enddo ! jj
            enddo ! i
        enddo ! j


!       Calculate ground temperature
        if(istat_sst .ne. 1)then
            write(6,*)' Note: sst_k not available'
        endif

        do j = 1,jmax
        do i = 1,imax
          if(rlaps_land_frac(i,j) .ge. 0.5)then
            t_gnd_k(i,j) = t_ground_k(t_sfc_k(i,j),solar_alt(i,j)
     1                               ,solar_ha(i,j),solar_dec,rlat(i,j)       
     1                               ,cvr_snow(i,j),r_missing_data,i,j
     1                               ,imax,jmax)
          else ! water environment
            if(istat_sst .eq. 1 
     1                   .and. sst_k(i,j) .ne. r_missing_data)then     
                t_gnd_k(i,j) = sst_k(i,j)
            else
                t_gnd_k(i,j) = t_sfc_k(i,j)
            endif
          endif
        enddo
        enddo

!       Check number of points thrown out as a function of offset to check
!       for IR navigation errors
        thresh1 = 5.
!       call correlation(t_gnd_k,tb8_k,thresh1,imax,jmax)
        call correlation(t_gnd_k,tb8_k,thresh_ir_diff1,imax,jmax)
        thresh3 = 15.
!       call correlation(t_gnd_k,tb8_k,thresh3,imax,jmax)

!       Find terrain locations
        do j=1,jmax
        do i=1,imax
             k_terrain(i,j) = 1 ! Minimum value allowed
        enddo ! i
        enddo ! j

        do kl=1,klaps
        do j=1,jmax
        do i=1,imax
             if(heights_3d(i,j,kl) .lt. topo(i,j))k_terrain(i,j) = kl
        enddo ! i
        enddo ! j
        enddo ! k


        I4_elapsed = ishow_timer()


        write(6,*)
        write(6,*)'  i  k   frac_k   cldtop_temp_f  cldtop_m',
     1                  '     sfc T (f)   Topo'

        init_co2 = 0 ! Initializes CO2 Slicing routine
        n_valid_co2 = 0
        n_missing_co2 = 0
        n_no_sao1 = 0
        n_no_sao2 = 0
        n_no_sao3 = 0

        do j=1,jmax
        do i=1,imax

         if(tb8_k(i,j) .ne. r_missing_data)then

!         Compare brightness temp to surface temperature
          if(j .eq. 29)then
              write(6,111,err=112)i,k_to_f(tb8_k(i,j))
     1                             ,k_to_f(t_gnd_k(i,j))
111           format(1x,i3,11x,f14.1,8x,f14.1)
112       endif

!         Calculate cloud top height from Band 8 and/or CO2 slicing method
          call cloud_top( init_co2,i4time,tb8_k(i,j)
!    1     ,cloud_frac_vis_a(i,j),istat_vis,cloud_frac_vis_a(i,j)
     1     ,t_gnd_k,pres_sfc_pa
     1     ,thresh_ir_diff1,topo(i,j),r_missing_data
     1     ,i,j,imax,jmax,klaps,heights_3d,temp_3d,k_terrain(i,j),laps_p       
     1     ,istat_39_a(i,j), l_use_39                                     ! I
     1     ,istat_39_add_a(i,j)                                           ! O
     1     ,l_co2                                                         ! I
     1     ,n_valid_co2,n_missing_co2,cldtop_m_co2(i,j),istat_co2         ! O
     1     ,cldtop_m_tb8(i,j),l_tb8                                       ! O
     1     ,cldtop_m(i,j),l_cloud_present                                 ! O
     1     ,sat_cover)                                                    ! O

          if(l_co2 .and. istat_co2 .eq. 1)then ! Using CO2 slicing method

!           Clear out those levels higher than what the satellite is showing
!           with the co2 slicing method

            do k=kcld,1,-1
              if(cldcv(i,j,k) .gt. .04)then ! Efficiency test
                if(cldtop_m_co2(i,j) .ne. r_missing_data)then
                  if(cld_hts(k) .gt. cldtop_m_co2(i,j))then
                    cldcv(i,j,k) = default_clear_cover ! Include surface buffer?
                  endif
                endif
              endif
            enddo ! k

          endif

          if(.true.)then ! Use Band 8 (11.2 mm)

!           Modify those levels where the satellite shows warmer than the
!           calculated brightness temp from the analysis of SAO/Pireps

            do k=kcld,1,-1

              if(cldcv(i,j,k) .gt. .04)then ! Efficiency test

                z_temp = height_to_zcoord2(cld_hts(k),heights_3d,
     1                                     imax,jmax,klaps,i,j,istatus1)       

                if(istatus1 .ne. 1)then
                    if(cld_hts(k) .gt. heights_3d(i,j,klaps))then
!                       We cannot perform the cloud clearing when the cloud
!                       height is above the top of the pressure grid
                        go to 500
                    else
                        write(6,*)
     1                   ' Error status returned from height_to_zcoord2'       
     1                  ,k,klaps,cld_hts(k),heights_3d(i,j,klaps)
                        istatus = 0
                        return
                    endif
                endif
                z_temp = max(1.,min(z_temp,float(klaps)-.001))
                iz_temp = int(z_temp)
                frac = z_temp - iz_temp
                temp_grid_point = temp_3d(i,j,iz_temp)    * (1. - frac)  
     1                          + temp_3d(i,j,iz_temp+1)  * frac

                if(cldcv(i,j,k) .ne. r_missing_data)then
                  if(cldcv(i,j,k) .gt. 1.050)then   ! excessively over 1.0
                      write(6,*)' Error in insert_sat, cldcv > 1.050'
     1                         ,i,j,k,cldcv(i,j,k)
                      istatus = 0
                      return
                  elseif(cldcv(i,j,k) .gt. 1.0)then ! slightly over 1.0
                      write(6,*)' Warning in insert_sat, cldcv > 1.0'       
     1                         ,i,j,k,cldcv(i,j,k)
                      write(6,*)' Resetting cloud cover to 1.0'
                      cldcv(i,j,k) = 1.0
                  endif

                  tb8_calculated_rad = 
     1                temp_to_rad(temp_grid_point) * cldcv(i,j,k) +
     1                temp_to_rad(t_gnd_k(i,j))    * (1.-cldcv(i,j,k))       

                  if(tb8_calculated_rad .lt. 0.)then
                      write(6,*)' Error, tb8_calculated_rad < 0 '
     1                                  ,tb8_calculated_rad
                      write(6,*)'i,j,k,cldcv',i,j,k,cldcv(i,j,k)
                      write(6,*)'temp_grid_point,t_gnd_k'
     1                          ,temp_grid_point,t_gnd_k(i,j)
                      istatus = 0
                      return
                  endif

                  tb8_calculated = rad_to_temp(tb8_calculated_rad)

                  tb8_calculated_test = temp_grid_point *  cldcv(i,j,k) 
     1                                + t_gnd_k(i,j) * (1.-cldcv(i,j,k))
                else
                  tb8_calculated = t_gnd_k(i,j)
                endif


!               Test if clouds detected by SAO/PIREP should have been
!               detected by the satellite (if satellite is warmer than analysis)
                if(tb8_calculated - tb8_k(i,j) .lt. -0.)then

!                 Don't touch points within buffer of surface
                  if(cld_hts(k) - topo(i,j) .gt. surface_sao_buffer)then
                    if(.false.)then
                      cldcv(i,j,k)=default_clear_cover
                    else ! .true.
!                     Does satellite still imply at least some cloud?
                      if(t_gnd_k(i,j) - tb8_k(i,j)  .gt. 8.0)then ! Some cloud
                        if(cldcv(i,j,k) .gt. 0.9)then ! Lower top of solid cld
                            cldcv(i,j,k)=default_clear_cover
                        else                          ! Cover < 0.9, correct it
                            cldcv(i,j,k) = band8_cover( tb8_k(i,j)
     1                          ,t_gnd_k(i,j),temp_grid_point)
                            if(cldcv(i,j,k) .gt. 1.0 .or.
     1                         cldcv(i,j,k) .lt. 0.0       )then
                                write(6,*)' ERROR--cover out of bounds'
                                istatus = 0
                                return
                            endif
                        endif
                      else ! Band 8 nearly matches ground, clear it
!                       Insure that "black (or grey) stratus" is not present
                        temp_thresh = 
     1                            min(t_gnd_k(i,j),t_sfc_k(i,j)-10.0)       
                        if(temp_grid_point .lt. temp_thresh)then
!                       if(temp_grid_point .lt. t_gnd_k(i,j))then
                            cldcv(i,j,k)=default_clear_cover ! not in inversion, 
                                                             ! clear it out
                        endif
                      endif ! IR signature present
                    endif ! .true.
                  endif ! Analysis has cloud above surface buffer

                endif ! SAO grid point is colder than SAT

              endif ! Current Cloud Cover is significant (> .04)

 500        enddo ! k (for clearing clouds)

          endif ! Using Band 8

!         Test if we're confident that a cloud is present and we know where
!         the cloud top is.
          IF(l_cloud_present) then ! Insert satellite clouds

!           Set initial satellite cloud base
            htbase_init=cldtop_m(i,j) - thk_def
            htbase = htbase_init

!           Locate lowest SAO cloud base
            ht_sao_base = 1e30

            cldcv_above = cldcv_sao(i,j,kcld)

            do k=kcld-1,1,-1
              if(cldcv_sao(i,j,k) .le. thr_sao_cvr .and.
     1         cldcv_above      .gt. thr_sao_cvr)then
                    ht_sao_base = cld_hts(k+1)
              endif

              cldcv_above = cldcv_sao(i,j,k)

            enddo ! k

            i_sao = i
            j_sao = j

!           if(.false.)then ! Search for nearby SAO cloud layers
            if(ht_sao_base .eq. 1e30)then ! Search for nearby SAO cloud layers
                                          ! because the satellite says cloud
                                          ! and the SAO doesn't
              n_no_sao1 = n_no_sao1 + 1

              do jdelt_index = 1,njdelt
              jj = j + jdelt(jdelt_index)

              do idelt_index = 1,nidelt
              ii = i + idelt(idelt_index)
                if(ii .ge. 1 .and. ii .le. imax .and.
     1             jj .ge. 1 .and. jj .le. jmax          )then
                    cldcv_above = cldcv_sao(ii,jj,kcld)
                    do k=kcld-1,1,-1
                      if(cldcv_sao(ii,jj,k) .le. thr_sao_cvr .and.
     1                 cldcv_above        .gt. thr_sao_cvr)then
                          ht_sao_base = cld_hts(k+1)
                      endif
                      cldcv_above = cldcv_sao(ii,jj,k)
                    enddo ! k
                    if(ht_sao_base .ne. 1e30)then
                        i_sao = ii
                        j_sao = jj
                        goto201
                    endif
                endif  ! In bounds
              enddo ! ii
              enddo ! jj
            endif

201         if(ht_sao_base .eq. 1e30)then ! Satellite cloud but no SAO cloud
              n_no_sao2 = n_no_sao2 + 1
              cover=sat_cover
              htbase_init = ht_sao_base

              if(tb8_k(i,j) - t_gnd_k(i,j) .lt. -thresh_ir_diff2)then 
                  buffer = 2100.             ! We more likely have a cloud
              else                            
                  buffer = surface_ir_buffer ! Weed out IR tops w/higher buffer
              endif

!             Calculate new cloud top and cover (based on filtered tb8)
!             This gives a cloud edge with more uniform height for an isolated
!             cloud when the edge has a "soft" appearance in the imagery.

              cldtop_m_avg = cldtop_m(i,j)
              call cloud_top(init_co2,i4time,tb8_cold_k(i,j)
!    1            ,cloud_frac_vis_a(i,j),istat_vis,cloud_frac_co2_dum
     1            ,t_gnd_k,pres_sfc_pa
     1            ,thresh_ir_diff1,topo(i,j),r_missing_data
     1            ,i,j,imax,jmax,klaps,heights_3d,temp_3d
     1            ,k_terrain(i,j),laps_p
     1            ,istat_39_a(i,j), l_use_39                             ! I
     1            ,istat_39_add_dum                                      ! O
     1            ,l_co2                                                 ! I
     1            ,n_valid_co2,n_missing_co2,cldtop_m_co2(i,j),istat_co2 ! O
     1            ,cldtop_m_tb8(i,j),l_tb8                               ! O
     1            ,cldtop_m(i,j),l_cloud_present                         ! O
     1            ,sat_cover)                                            ! O



!             Calculate the cover (opacity) given the brightness temperature,
!             ground temperature, and assumed ambient cloud-top temperature.
              cover = 
     1              band8_cover(tb8_k(i,j),t_gnd_k(i,j),tb8_cold_k(i,j))       

              htbase = max( topo(i,j) + buffer , cldtop_m(i,j)-thk_def )

              if(htbase .gt. cldtop_m(i,j))then
                  n_no_sao3 = n_no_sao3 + 1
              else
                  write(6,211,err=212)i,j,tb8_k(i,j),tb8_cold_k(i,j)
     1                   ,cldtop_m_avg,cldtop_m(i,j)
!       1                        ,cldtop_m_avg,cldtop_m_cold
     1                   ,cover
211               format(' AVG/COLD',2i4,2f7.1,2x,2f7.0,f9.3)
212           endif

            elseif(ht_sao_base .gt. cldtop_m(i,j))then ! Satellite top below ceiling
              cover=sat_cover
              htbase_init = ht_sao_base
              htbase = htbase_init
              cldtop_old = cldtop_m(i,j)
              cldtop_m(i,j) = htbase_init + thk_def

!             Find a thinner value for cloud cover consistent with the new
!             higher cloud top and the known brightness temperature.
              if(.true.)then ! Should this depend on co2?

!                 Note that cover is not really used here as an input
                  call correct_cover(cover,cover_new,cldtop_old
     1                              ,cldtop_m(i,j)
     1                              ,temp_3d,tb8_k(i,j),t_gnd_k(i,j)
     1                              ,heights_3d
     1                              ,imax,jmax,klaps,i,j,istatus)
                  if(istatus .ne. 1)then
                      write(6,*)' Error in correct_cover'
                      write(6,*)cldtop_old,cldtop_m(i,j)
     1                         ,htbase_init,thk_def,cover
                      write(6,*)(heights_3d(i,j,k),k=1,klaps)
                      return
                  endif
                  cover = cover_new
              endif ! l_co2

            else ! Normal use of satellite data
              cover=sat_cover

!             Locate SAO cloud base below satellite cloud top, modify
!             satellite cloud base. Highest SAO ceiling within default thickness
!             range of satellite layer is used.
              do k=kcld,1,-1

                if(       cld_hts(k) .ge. htbase_init
     1             .and.  cld_hts(k) .le. cldtop_m(i,j)      )then

                  if(cldcv_sao(i_sao,j_sao,k)   .le. thr_sao_cvr .and.
     1               cldcv_sao(i_sao,j_sao,k+1) .gt. thr_sao_cvr
     1                                       )then ! We have an SAO base
                    htbase = cld_hts(k+1)

!                   If SAO (hence satellite) base is above the satellite
!                   cloud top, lower the satellite base by one grid level
                    if(htbase .gt. cldtop_m(i,j))then
                        htbase = cld_hts(k)
                    endif

                    goto301

                  endif ! We have an SAO base

                endif ! in satellite layer

              enddo ! k

            endif ! Satellite top below SAO ceiling


301         if(htbase .ne. htbase_init)then
!                write(6,*)' Satellite ceiling reset by SAO',i,j,htbase
!       1       ,htbase_init,cldtop_m(i,j)
            endif

!           Add satellite cloud to array
            do k=kcld,1,-1
              if(cld_hts(k) .ge. htbase  .and.
     1           cld_hts(k) .le. cldtop_m(i,j) )then ! in satellite layer
                 cldcv(i,j,k)=cover
              endif
            enddo

          ENDIF ! l_cloud_present (Cloudy)

         endif ! tb8_k(i,j) .ne. r_missing_data

        enddo ! imax
        enddo ! jmax

!       Write stats on CO2 and Band 8 (11.2mm) methods
        write(6,*)' n_valid_co2 = '  ,n_valid_co2
     1           ,' n_missing_co2 = ',n_missing_co2
        write(6,*)' n_no_sao (1/2/3) = ',n_no_sao1,n_no_sao2,n_no_sao3

        I4_elapsed = ishow_timer()

        call compare_radiation(kcld,temp_3d,klaps,imax,jmax
     1      ,cldcv,cldcv_1d,cld_hts,t_sfc_k,t_gnd_k,tb8_k
     1               ,r_missing_data,cvr_snow,heights_3d,nlyr,istatus)

999     return
        end

        subroutine cloud_top( init_co2,i4time,tb8_k
!    1  ,cloud_frac_vis,istat_vis,cloud_frac_co2
     1  ,t_gnd_k,pres_sfc_pa,thresh_ir_diff1,topo,r_missing_data
     1  ,i,j,imax,jmax,klaps,heights_3d,temp_3d,k_terrain,laps_p
     1  ,istat_39, l_use_39                                            ! I
     1  ,istat_39_add                                                  ! O
     1  ,l_co2                                                         ! I
     1  ,n_valid_co2,n_missing_co2,cldtop_m_co2,istat_co2              ! O
     1  ,cldtop_m_tb8,l_tb8                                            ! O
     1  ,cldtop_m,l_cloud_present                                      ! O
     1  ,sat_cover)                                                    ! O


!       This routine computes the cloud top height given a band 8 brightness
!       temperature and 3D fields of temp and height. The CO2 method is also
!       employed to yield a cloud top pressure and height. This is not yet
!       fully integrated into the final cloud analysis.

        real*4 zeros
        parameter (zeros = 1.e-30)

!       Argument list
        integer*4 init_co2                      ! Input/Output
        integer*4 i4time                        ! Input
        real*4 tb8_k                            ! Input
!       real*4 cloud_frac_vis                   ! Input
        integer*4 i,j,imax,jmax,klaps           ! Input
        real*4 t_gnd_k(imax,jmax)               ! Input
        real*4 pres_sfc_pa(imax,jmax)           ! Input
        real*4 thresh_ir_diff1                  ! Input
        real*4 topo                             ! Input
        real*4 r_missing_data                   ! Input
!       integer*4 istat_vis                     ! Input
        real*4 heights_3d(imax,jmax,klaps)      ! Input
        real*4 temp_3d(imax,jmax,klaps)         ! Input
        real*4 k_terrain                        ! Input
        real*4 laps_p(klaps)                    ! Input
        integer*4 n_valid_co2,n_missing_co2     ! Input/Output
        real*4 cldtop_m_co2                     ! Output
        logical l_co2,l_use_39                  ! Input
        integer istat_co2                       ! Output
        real*4 cldtop_m_tb8                     ! Output
        logical l_tb8                           ! Output
        real*4 cldtop_m                         ! Output
        logical l_cloud_present                 ! Output
        real*4 sat_cover                        ! Output

!       Local
!       real*4 dum_3d(imax,jmax,klaps)          ! Local (Dummy array for Q)
        real*4 arg,frac_k,temp_above,cldtop_temp_k
        integer*4 kl
!       real*4 ppcc(8)

!       Function call
        real*4 k_to_f

!       Call the CO2 slicing method to get cloud tops

        if(init_co2 .eq. 0)istatus_co2 = 0
        istatus_co2 = -3

!       init_co2 is a flag telling whether it is the first time the routine
!       is being called so data/array initialization or satellite data
!       ingest can be performed.

!       if(istatus_co2 .ne. -3)
!       1       call co2slice (i,j,imax,jmax,klaps,i4time,init_co2,
!       1       laps_p,temp_3d,t_gnd_k,pres_sfc_pa,dum_3d,
!       1       ppcc, iipcld, iiacld, iitcld, iiflg, istatus_co2)

c       input parameters

!       integer i,j,imax,jmax,klaps,i4time,init ! coords and dimensions, time
c       and init parameter.  init has the following meaning.
c
c       init = 0 this is the first call of the routine, it must get the grids
c       init = 1 this is a subsequent call to the routine and therefore, it
c                will skip getting the grids  -- furthermore, this module will
c               affect changes in the value of init insofar as changing it
c               from 0 to 1.  This module will not change the value from 1 to
c               0.. that is up to the calling program.

c       output parameters

!       real ppcc(8) ! defined by co2 slicing code
!       integer iipcld, iiacld, iitcld, istatus, iiflg

c       iiflg is the co2 code's return status flag.  it is as follows..

c       -2 = not enough channels to run the code for this point
c       -1 = bad data encountered in the data. for this point
c       0  = all has run ok

c       in addition to the above nominal returns from the routine, we have
c       added the following parameter

c       -3 = data not available for this call ... do not call again for this
c       time period...( this is to allow you to save resources and continue on
c       in processing.)

C                  IIPCLD = PRESSURE CLOUD TOP (MB)
C                  IIACLD = CLOUD AMOUNT*EMISSIVITY (%)
C                  IITCLD = CLOUD TOP TEMP (K) FROM PROFILE
C                  PPCC(1) = CTP FOR RATIO 1 (3/5)
C                  PPCC(2) = CTP FOR RATIO 2 (3/4)
C                  PPCC(3) = CTP FOR RATIO 3 (4/5)
C                  PPCC(4) = CTP FOR RATIO 4 (5/8)
C                  PPCC(5) = CTP FOR WINDOW CHANNEL
C                  PPCC(6) = NUMBER OF BEST CTP
C                  PPCC(7) = EFFECTIVE CLOUD AMOUNT
C                  PPCC(8) = EFFECTIVE CLOUD AMOUNT FROM 5/8 RATIO

!       Convert CO2 cloud top from pressure to height
        if(istatus_co2 .eq. 0)then
            istat_co2 = 1
            cldtop_p_co2 = iipcld * 100.
            cldtop_z_co2 = zcoord_of_logpressure(cldtop_p_co2)

            iarg = int(cldtop_z_co2)
            frac = cldtop_z_co2 - float(iarg)
            if(iarg .lt. klaps)then
                cldtop_m_co2 = heights_3d(i,j,iarg) * (1.0 - frac) +
     1                 heights_3d(i,j,iarg+1) * frac
            else
                cldtop_m_co2 = r_missing_data
                write(6,*)' WARNING, CO2 pressure is out of bounds '
     1                   ,iipcld
            endif

        else
            istat_co2 = 0
            cldtop_m_co2 = r_missing_data

        endif

        if(cldtop_m_co2 .ne. r_missing_data)then
            n_valid_co2 = n_valid_co2 + 1
        else
            n_missing_co2 = n_missing_co2 + 1
            cldtop_m_co2 = r_missing_data 
        endif

!       This section finds the cloud top using Band 8 data and temperatures
!       Estimate whether tb8_k - t < threshold
        cldtop_m_tb8 = r_missing_data ! zeros
        if(tb8_k - t_gnd_k(i,j) .lt. -thresh_ir_diff1) then ! probably clds
            l_tb8 = .true.

        else ! No clouds according to SATELLITE (Band 8 - 11.2mm)
            l_tb8 = .false.

        endif

        if(l_tb8 .OR. (l_use_39 .and. istat_39 .eq. 1))then ! get 11u cloud top
            cldtop_temp_k = tb8_k

!           Correct the cloud top temperature for thin clouds using VIS data
!!!!!       call correct(tb8_k,tb8_cold_k,cloud_frac_vis,istat_vis,cldtop_temp_k)

!           Locate cloud top in 3-D Temperature Grid (Using lowest crossing point)
            temp_above = temp_3d(i,j,klaps)

            do kl = klaps-1,k_terrain,-1

                if( (temp_3d(i,j,kl) - cldtop_temp_k) *
     1              (temp_above      - cldtop_temp_k) .lt. 0.)then ! Crossing Pt

                    frac_k = (cldtop_temp_k - temp_3d(i,j,kl))
     1                    /  (temp_above    - temp_3d(i,j,kl))

                    arg = heights_3d(i,j,kl) + frac_k *
     1                   (heights_3d(i,j,kl+1) - heights_3d(i,j,kl))

                    if(arg .ge. topo)then
                        cldtop_m_tb8 = arg
                    endif

                    if(j .eq. 29)then
                        if(arg .lt. topo)then
                            write(6,*)
     1                      ' Cloud Top Below Ground - not used'
                        endif

                        write(6,111,err=121)cldtop_m_tb8,cldtop_m_co2
111                     format(1x,f10.0,1x,f10.0)

121                     write(6,122,err=123)i,kl,frac_k
     1                           ,k_to_f(cldtop_temp_k)
     1                           ,arg,topo,k_to_f(temp_3d(i,j,kl))
     1                           ,k_to_f(temp_above)
122                     format(1x,2i3,f8.3,f14.1,f11.1,f21.1,2f6.1)
123                 endif

                endif

                temp_above = temp_3d(i,j,kl)

            enddo ! kl

        endif ! We will want to use a 11u determined cloud top

        istat_39_add = 0

!       Set variables depending on whether in Band 8 or CO2 mode
        if(l_co2 .and. istat_co2 .eq. 1)then ! Using CO2 method
            if(cldtop_m_co2 .ne. r_missing_data)then
                l_cloud_present = .true.
            else
                l_cloud_present = .false.
            endif

            cldtop_m = cldtop_m_co2
!           sat_cover = PPCC(7)

        elseif( (.not. l_tb8) .AND. (l_use_39 .and. istat_39 .eq. 1) 
     1                        .AND. cldtop_m_tb8 .ne. r_missing_data 
     1                                                            ) then      

!           Band 8 (11mm) threshold says no but 3.9 micron says yes
!           We did still get a valid Band 8 derived cloud top

            l_cloud_present = .true.
            cldtop_m = cldtop_m_tb8
            sat_cover = 1.0 
            istat_39_add = 1

        else                          ! Using Band 8 (11.2mm) data only
            l_cloud_present = l_tb8
            cldtop_m = cldtop_m_tb8
            sat_cover = 1.0 

        endif

        return
        end


        subroutine correlation(t,tb8_k,thresh,ni,nj)

        integer*4 IBOX
        parameter (IBOX = 3)

        real*4 t(ni,nj),tb8_k(ni,nj)
        integer*4 i_corr_array(-IBOX:IBOX,-IBOX:IBOX)

        write(6,*)
        write(6,*)' Checking Navigation of IR satellite, thresh = '
     1                                          ,thresh

        thresh_minus = -thresh
        min_count = 999999

        do joff = -IBOX,IBOX
        do ioff = -IBOX,IBOX
            icount = 0

            do j = 1,nj
              jj = max(min(j+joff,nj),1)

              do i = 1,ni

                ii = max(min(i+ioff,ni),1)

                if(tb8_k(ii,jj) - t(i,j) .lt. thresh_minus)then
                    icount = icount+1 ! Clouds Indicated
                endif

              enddo ! i

            enddo ! j

            if(icount .lt. min_count)then
                min_count = icount
                ioff_min = ioff
                joff_min = joff
            endif

            i_corr_array(ioff,joff) = icount

        enddo ! ioff
        enddo ! joff

        do j = -ibox,ibox
            write(6,101)(i_corr_array(i,j),i=-ibox,ibox)
101         format(10i8)
        enddo

        write(6,201)min_count,ioff_min,joff_min
201     format('  Minimum of',i5,' points flagged with an offset of',2i4
     1)

        return
        end

        subroutine correct_cover(cover_in,cover_new_f,cldtop_old
     1             ,cldtop_new,temp_3d,tb8_k,t_gnd_k,heights_3d
     1             ,imax,jmax,klaps,i,j,istatus)

!       Find a thinner value for cloud cover consistent with the new
!       higher cloud top and the known brightness temperature.

        real*4 temp_3d(imax,jmax,klaps)
        real*4 heights_3d(imax,jmax,klaps)

        integer iwrite
        data iwrite /0/
        save iwrite

        istatus = 1

!       Find Temperature of old cloud top
        z_temp = height_to_zcoord2(cldtop_old,heights_3d,imax,jmax,klaps
     1                                          ,i,j,istatus)
        if(istatus .ne. 1)then
            return
        endif

        z_temp = max(1.,min(z_temp,float(klaps)-.001))
        iz_temp = int(z_temp)
        frac = z_temp - iz_temp
        temp_old = temp_3d(i,j,iz_temp)    * (1. - frac)
     1          +  temp_3d(i,j,iz_temp+1)  * frac

!       Find Temperature of new cloud top
        z_temp = height_to_zcoord2(cldtop_new,heights_3d,imax,jmax,klaps
     1                                          ,i,j,istatus)
        if(istatus .ne. 1)then
            return
        endif

        z_temp = max(1.,min(z_temp,float(klaps)-.001))
        iz_temp = int(z_temp)
        frac = z_temp - iz_temp
        temp_new = temp_3d(i,j,iz_temp)    * (1. - frac)
     1          +  temp_3d(i,j,iz_temp+1)  * frac

!       This one utilizes a linear approximation to the sigma T**4 relationship
        cover_old = min(cover_in,1.0)
        cover_new = cover_old *
     1  (tb8_k - t_gnd_k) / (temp_new - t_gnd_k)

!       This one utilizes the sigma T**4 relationship
        cover_new_f = band8_cover(tb8_k,t_gnd_k,temp_new)

        if(j .eq. int(j-9)/10*10+9)then
            iwrite = iwrite + 1
            if(iwrite .lt. 60)then
                write(6,1,err=2)i,j,t_gnd_k,temp_old,temp_new,cldtop_old
     1                         ,cldtop_new,cover_new,cover_new_f
1               format(1x,'Corr-cvr ',2i3,3f7.0,2f8.0,2f8.2)
2           endif
        endif

        return
        end


        function band8_cover(tb8_k,t_gnd_k,t_cld)

        r_sfc = temp_to_rad(t_gnd_k)
        r_sat = temp_to_rad(tb8_k)
        r_cld = temp_to_rad(t_cld)

        band8_cover = (r_sat - r_sfc) / (r_cld - r_sfc)

        if(band8_cover .gt. 1.0)then
            write(6,*)' WARNING: resetting band8_cover down to 1.0'
            write(6,*)' tb8_k,t_gnd_k,t_cld,band8_cover'
     1                 ,tb8_k,t_gnd_k,t_cld,band8_cover
            band8_cover = 1.0 
        endif

        return
        end

        function temp_to_rad(temp)

        data init /0/
        save init
        if(init .eq. 0)then
            init = 1
            NSAT = 3
            call PLNKIV(NSAT)
        endif

        temp_to_rad = temp
        temp_to_rad = VPLANC(temp,8)

        return
        end

        function rad_to_temp(rad)

        data init /0/
        save init
        if(init .eq. 0)then
            init = 1
            NSAT = 3
            call PLNKIV(NSAT)
        endif

        rad_to_temp = rad
        rad_to_temp = VBRITE(rad,8)

        return
        end

        subroutine compare_radiation(kcld,temp_3d,klaps,imax,jmax
     1      ,cldcv,cldcv_1d,cld_hts,t_sfc_k,t_gnd_k,tb8_k
     1           ,r_missing_data,cvr_snow,heights_3d,nlyr,istatus)

        include 'laps_cloud.inc'

        real*4 temp_3d(imax,jmax,klaps)
        real*4 heights_3d(imax,jmax,klaps)
        real*4 cvr_snow(imax,jmax)
        real*4 tb8_k(imax,jmax),t_gnd_k(imax,jmax),t_sfc_k(imax,jmax)
        real*4 a(100),f(100)
        integer*4 ilyr(KCLOUD) ! Dimension needs to be changed to KCLOUD
        real*4 a_new(100),f_new(100)
        integer*4 ilyr_new(KCLOUD) ! Dimension needs to be changed to KCLOUD
        real*4 cldcv(imax,jmax,kcld)
        real*4 cldcv_1d(kcld) ! ,cld_hts(kcld)
        logical l_correct,l_output

!       This routine compares the analyzed clouds to the 11.2mm radiation
!       and determines adjusted cloud fractions of the cloud layers to yield
!       a better fit.

        write(6,*)' Comparing radiation'

        corr_thr = 4.
        cover_step = 0.05
        iter_max = 10

        iwrite = 0
        n_clear = 0
        n_clear_ns = 0
        n_clear_sn = 0
        n_cloudy = 0
        n_total = 0
        n_correct = 0
        n_iter = 0
        tdiff_sumsq = 0.
        tdiff_cld_sumsq = 0.
        tdiff_corr_sumsq = 0.
        tdiff_corr_sumsq2 = 0.
        tdiff_corr_sumsq3 = 0.
        tdiff_corr_cld_sumsq = 0.
        tb8_g_clr_sum   = 0.
        tb8_a_clr_sum   = 0.
        tb8_g_clr_sumsq = 0.
        tb8_a_clr_sumsq = 0.
        tb8_g_clr_ns_sum   = 0.
        tb8_a_clr_ns_sum   = 0.
        tb8_g_clr_ns_sumsq = 0.
        tb8_a_clr_ns_sumsq = 0.
        tb8_g_clr_sn_sum   = 0.
        tb8_a_clr_sn_sum   = 0.
        tb8_g_clr_sn_sumsq = 0.
        tb8_a_clr_sn_sumsq = 0.

        do j = 1,jmax
          if(j .eq. int(j-9)/10*10+9)then
            l_output = .true.
          else
            l_output = .false.
          endif

          do i = 1,imax
           if(tb8_k(i,j) .ne. r_missing_data)then
            do k = 1,kcld
                cldcv_1d(k) = cldcv(i,j,k)
            enddo

            call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j,imax,jmax
     1      ,a,f,ilyr,cldcv_1d,cld_hts,t_gnd_k(i,j),heights_3d
     1                                  ,t_effective,nlyr,istatus)
            tdiff = tb8_k(i,j)-t_effective ! Band 8 - calculated
            frac_clouds = 1.-f(nlyr)

            if(nlyr .eq. 1)then
                tb8_g_clr_sum   = tb8_g_clr_sum   + tdiff
                tb8_g_clr_sumsq = tb8_g_clr_sumsq + tdiff**2
                tb8_a_clr_sum   = tb8_a_clr_sum   + tb8_k(i,j) 
     1                                            - t_sfc_k(i,j)
                tb8_a_clr_sumsq = tb8_a_clr_sumsq +
     1                                  (tb8_k(i,j) - t_sfc_k(i,j))**2
                n_clear = n_clear + 1

                if(cvr_snow(i,j) .ne. r_missing_data)then
                    if(cvr_snow(i,j) .lt. 0.5)then
                        tb8_g_clr_ns_sum   = tb8_g_clr_ns_sum   + tdiff
                        tb8_g_clr_ns_sumsq = tb8_g_clr_ns_sumsq + tdiff*
     1*2
                        tb8_a_clr_ns_sum   = tb8_a_clr_ns_sum
     1                             + tb8_k(i,j) - t_sfc_k(i,j)
                        tb8_a_clr_ns_sumsq = tb8_a_clr_ns_sumsq +
     1                              (tb8_k(i,j) - t_sfc_k(i,j))**2
                        n_clear_ns = n_clear_ns + 1
                    else
                        tb8_g_clr_sn_sum   = tb8_g_clr_sn_sum   + tdiff
                        tb8_g_clr_sn_sumsq = tb8_g_clr_sn_sumsq + tdiff*
     1*2
                        tb8_a_clr_sn_sum   = tb8_a_clr_sn_sum
     1                             + tb8_k(i,j) - t_sfc_k(i,j)
                        tb8_a_clr_sn_sumsq = tb8_a_clr_sn_sumsq +
     1                              (tb8_k(i,j) - t_sfc_k(i,j))**2
                        n_clear_sn = n_clear_sn + 1
                    endif
                endif
            else
                n_cloudy = n_cloudy + 1
                tdiff_cld_sumsq = tdiff_cld_sumsq + tdiff**2

            endif

            n_total = n_total + 1
            tdiff_sumsq = tdiff_sumsq + tdiff**2

!           Apply corrections?
            if(nlyr .ge. 2 .and. abs(tdiff) .gt. corr_thr
     1                   .and. t_gnd_k(i,j) - tb8_k(i,j) .gt. 8.
     1             .and. frac_clouds .gt. 0.4             )then
                l_correct = .true.
                n_correct = n_correct + 1
                tdiff_corr_sumsq3 = tdiff_corr_sumsq3 + tdiff**2
            else
                l_correct = .false.
            endif

            if(abs(tdiff) .gt. corr_thr .and. l_output)then
                write(6,901,err=902)i,j,nlyr-1,tb8_k(i,j),t_effective
     1               ,tdiff,frac_clouds,l_correct,(a(l),l=nlyr-1,1,-1)       
901             format(1x,2i4,i3,2f7.0,f7.1,f7.3,l2,' *',10f6.2)
902         endif

!           This corrective section is now turned on
            iter = 0
            delta_cover = 0.

905         if(l_correct)then
                iter = iter + 1
                tdiff_ref = tdiff
                delta_cover_ref = delta_cover

                if(tdiff .lt. 0.)then
                    delta_cover = delta_cover + cover_step
                else
                    delta_cover = delta_cover - cover_step
                endif

                call apply_correction(cldcv,i,j,imax,jmax
     1                       ,cldcv_1d,1,1,1,1,kcld,f,ilyr,delta_cover)

                call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j
     1              ,imax,jmax,a_new,f_new,ilyr_new,cldcv_1d,cld_hts       
     1              ,t_gnd_k(i,j),heights_3d
     1              ,t_effective,nlyr_new,istatus)

                tdiff = tb8_k(i,j)-t_effective
                frac_clouds = 1.-f_new(nlyr_new)

!               Continue to apply corrections?
                if(nlyr_new .ge. 2 .and. frac_clouds .gt. 0.4)then
                    i_correct = 1
                else
                    i_correct = 0
                endif


                if(iwrite .le. 100 .and. l_output)then
                    iwrite = iwrite + 1
                    write(6,911,err=912)i,j,nlyr_new-1,tb8_k(i,j)
     1               ,t_effective
     1               ,tdiff,frac_clouds,i_correct
     1               ,(a_new(l),l=nlyr_new-1,1,-1)
911                 format(1x,2i3,i3,2f7.0,f7.1,f7.3,i2,'  ',10f6.2)
912             endif

                if(i_correct .eq. 1 .and. iter .lt. iter_max
     1          .and. abs(tdiff) .lt. abs(tdiff_ref)
     1          .and.     tdiff * tdiff_ref .gt. 0.
     1                          )goto905 ! Loop back & increment cover

                n_iter = n_iter + iter

!               Final iteration
                if(.false. .or. tdiff*tdiff_ref .ge. 0.)then
                  ! Select best of last two increments
                    if(abs(tdiff) .ge. abs(tdiff_ref))then
                        delta_cover = delta_cover_ref
                        tdiff = tdiff_ref
                    endif

                else           ! Do one Newton iteration
                    frac = tdiff / (tdiff - tdiff_ref)
                    delta_cover = delta_cover
     1                  + frac * (delta_cover_ref - delta_cover)

                    call apply_correction(cldcv,i,j,imax,jmax
     1                        ,cldcv_1d,1,1,1,1,kcld,f,ilyr,delta_cover)

                    call cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j
     1                                       ,imax,jmax
     1                                       ,a_new,f_new,ilyr_new
     1                                       ,cldcv_1d,cld_hts
     1                                       ,t_gnd_k(i,j),heights_3d       
     1                                       ,t_effective,nlyr_new
     1                                       ,istatus)
                    tdiff = tb8_k(i,j)-t_effective
                    n_iter = n_iter + 1
                endif

!               if(iwrite .le. 100 .and. l_output)write(6,913)delta_cover,tdiff
                if(iter .gt. 0 .and. l_output)write(6,913)delta_cover,td
     1iff
913             format(70x,'Best delta_cover/tdiff = ',f6.2,f6.1)

!               Apply correction to 3D cloud cover field
                call apply_correction(cldcv,i,j,imax,jmax
     1                   ,cldcv,i,j,imax,jmax,kcld,f,ilyr,delta_cover)

            endif ! Corrections were made

            tdiff_corr_sumsq     = tdiff_corr_sumsq + tdiff**2
            if(nlyr .gt. 1)then ! Cloudy (at least before adjustment)
                tdiff_corr_cld_sumsq = tdiff_corr_cld_sumsq + tdiff**2
            endif

            if(l_correct)then
                tdiff_corr_sumsq2 = tdiff_corr_sumsq2 + tdiff**2
            endif

           endif ! tb8_k(i,j) is valid
          enddo ! i
        enddo ! j

!       Write out statistics on consistency of Band 8 data and cloud cover
        if(n_clear .gt. 0)then
            write(6,951,err=960)n_clear, tb8_g_clr_sum  /float(n_clear)
     1                 ,sqrt    (tb8_g_clr_sumsq/float(n_clear))
951         format(/
     1      ' Mean/RMS band 8/gnd  temp residual in clear skies =  '
     1                                                ,i5,2f9.3)
960         write(6,961,err=970)n_clear, tb8_a_clr_sum  /float(n_clear)
     1                 ,sqrt    (tb8_a_clr_sumsq/float(n_clear))
961         format(' Mean/RMS band 8/air  temp residual in clear skies =
     1  '
     1                                                ,i5,2f9.3)
970     endif

        if(n_clear_ns .gt. 0)then
            write(6,956,err=965)n_clear_ns
     1                         ,tb8_g_clr_ns_sum  /float(n_clear_ns)
     1                       ,sqrt(tb8_g_clr_ns_sumsq/float(n_clear_ns))
956         format(/
     1      ' Mean/RMS band 8/gnd temp resid in clear/nsnow skies ='
     1                                                ,i5,2f9.3)
965         write(6,966,err=975)n_clear_ns
     1                         ,tb8_a_clr_ns_sum  /float(n_clear_ns)
     1                       ,sqrt(tb8_a_clr_ns_sumsq/float(n_clear_ns))
966         format(
     1      ' Mean/RMS band 8/air temp resid in clear/nsnow skies ='
     1                                                ,i5,2f9.3)
975     endif

        if(n_clear_sn .gt. 0)then
            write(6,976,err=985)n_clear_sn
     1                         ,tb8_g_clr_sn_sum  /float(n_clear_sn)
     1                       ,sqrt(tb8_g_clr_sn_sumsq/float(n_clear_sn))
976         format(/
     1      ' Mean/RMS band 8/gnd temp resid in clear/snow skies = '
     1                                                ,i5,2f9.3)
985         write(6,986,err=995)n_clear_sn
     1                         ,tb8_a_clr_sn_sum  /float(n_clear_sn)
     1                       ,sqrt(tb8_a_clr_sn_sumsq/float(n_clear_sn))
986         format(
     1      ' Mean/RMS band 8/air temp resid in clear/snow skies = '
     1                                                ,i5,2f9.3)
995     endif

        if(n_total .gt. 0)then
            write(6,971,err=980)n_total,sqrt(tdiff_sumsq/float(n_total))
971         format(/
     1      ' RMS  band 8/teff residual (bfr corr) in all  skies = '
     1                                                ,i5, f9.3)
980         write(6,981,err=990)n_total
     1                         ,sqrt(tdiff_corr_sumsq/float(n_total))
981         format(
     1      ' RMS  band 8/teff residual (aft corr) in all  skies = '
     1                                                ,i5, f9.3)
990     endif

        if(n_cloudy .gt. 0)then
            write(6,991,err=1000)n_cloudy
     1                          ,sqrt(tdiff_cld_sumsq/float(n_cloudy))
991         format(
     1      ' RMS  band 8/teff residual (bfr corr) in cldy skies = '
     1                                                ,i5, f9.3)
1000        write(6,1001,err=1010)
     1              n_cloudy,sqrt(tdiff_corr_cld_sumsq/float(n_cloudy))
1001        format(
     1      ' RMS  band 8/teff residual (aft corr) in cldy skies = '
     1                                                ,i5, f9.3)
1010    endif

        if(n_correct .gt. 0)then
            write(6,1011,err=1020)
     1          n_correct,sqrt(tdiff_corr_sumsq3/float(n_correct))
1011        format(/
     1      ' RMS  band 8/teff resid (bfr corr - corrected pts)  = '
     1                                                ,i5, f9.3)
1020        write(6,1021,err=1030)
     1          n_correct,sqrt(tdiff_corr_sumsq2/float(n_correct))
1021        format(
     1      ' RMS  band 8/teff resid (aft corr - corrected pts)  = '
     1                                                ,i5, f9.3)
            write(6,*)' Total/Average # of iterations = ',n_iter,
     1                          float(n_iter)/float(n_correct)
            write(6,*)
1030    endif

        return
        end

        subroutine cvr_to_tb8_effective(kcld,temp_3d,klaps,i,j,ni,nj,a
     1                                 ,f,ilyr,cvr,cld_hts,t_gnd_k
     1                                 ,heights_3d,t_effective,nlyr
     1                                 ,istatus)

        real*4 thresh_cvr
        parameter (thresh_cvr = 0.1)

        real*4 a(100)          ! Cloud fractions of layers
        real*4 f(100)          ! Apparent "cross-section" of cloud layers seen from above
        integer*4 ik(100)      ! Height level representative of cloud layers
        integer*4 ilyr(kcld)   ! Layer index for each cloud lvl (needs KCLOUD)
        real*4 cvr(kcld)       ! Cloud cover from analysis
        real*4 temp_3d(ni,nj,klaps),temp_lyr(100)
        real*4 heights_3d(ni,nj,klaps)
        real*4 cld_hts(kcld)

!       Convert from cloud cover to discreet cloud layer indices (cvr to a)
        nlyr = 0
        do k = kcld-1,1,-1
            if(cvr(k)   .ge. thresh_cvr .and.
     1       cvr(k+1) .lt. thresh_cvr)then      ! Top of layer
                nlyr = nlyr + 1
                a(nlyr) = cvr(k)
                ik(nlyr) = k

             else
                if(nlyr .ge. 1)then
                    if(cvr(k) .gt. a(nlyr))then   ! Max within layer
                        a(nlyr) = cvr(k)
                        ik(nlyr) = k
                    endif
                endif
             endif

             if(cvr(k) .ge. thresh_cvr)then       ! Still within layer
                 ilyr(k) = nlyr
             else                                 ! Below layer
                 ilyr(k) = 0
             endif

        enddo ! k

!       Get temperatures of the layers
        do n = 1,nlyr
            k = ik(n)
            z_temp = height_to_zcoord2(cld_hts(k),heights_3d,ni,nj,klaps
     1                                          ,i,j,istatus)
            if(istatus .ne. 1)then
                return
            endif

            z_temp = max(1.,min(z_temp,float(klaps)-.001))
            iz_temp = int(z_temp)
            frac = z_temp - iz_temp
            temp_lyr(n) = temp_3d(i,j,iz_temp)    * (1. - frac)
     1               +  temp_3d(i,j,iz_temp+1)  * frac
        enddo ! n

!       Add a layer for the ground
        nlyr = nlyr + 1
        a(nlyr) = 1.0
        temp_lyr(nlyr) = t_gnd_k

!       Convert cloud layer fractions to "cross-section" seen from satellite
!       This solves for the f array given the a array
        a(1) = min(a(1),1.0)
        f(1) = a(1)
        sumf = f(1)
        if(nlyr .ge. 2)then
            do n = 2,nlyr
                a(n) = min(a(n),1.0)
                f(n) = a(n) * (1.0 - sumf)
                sumf = sumf + f(n)
            enddo ! n
        endif ! nlyr

!       Calculate total radiance from all cloud layers + ground
        rsum = 0
        do n = 1,nlyr
            rsum = rsum + temp_to_rad(temp_lyr(n)) * f(n)
        enddo ! n

!       Convert to effective temperature and compare to observed brightness temp
        t_effective = rad_to_temp(rsum)

        return
        end


        subroutine apply_correction(cldcv_in,i_in,j_in,imax_in,jmax_in
     1                             ,cldcv_out,i_out,j_out,imax_out
     1                             ,jmax_out,kcld,f,ilyr,delta_cover)

        real*4 cldcv_in(imax_in,jmax_in,kcld)
        real*4 cldcv_out(imax_out,jmax_out,kcld)
        real*4 f(100)
        integer*4 ilyr(kcld)    ! Dimension needs to be changed to KCLOUD

!       Apply correction to 3D cloud cover field
        do k = 1,kcld
            if(ilyr(k) .gt. 0)then
                if(f(ilyr(k)) .gt. 0.)then
                    cldcv_out(i_out,j_out,k) =
     1          min(cldcv_in(i_in,j_in,k),1.0) + delta_cover
                    cldcv_out(i_out,j_out,k) =
     1          max(min(cldcv_out(i_out,j_out,k),1.0),0.0)
                endif
            endif
        enddo

        return
        end
