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

        subroutine get_cloud_deriv(ni,nj,nk,clouds_3d,cld_hts,            ! I
     1                          temp_3d,rh_3d_pct,heights_3d,pres_3d,     ! I
     1                          istat_radar,                              ! I/O
     1                          radar_3d,grid_spacing_cen_m,              ! I  
     1                          l_mask_pcptype,                           ! O
     1                          ibase_array_lwc,itop_array_lwc,           ! O
     1                          iflag_slwc,slwc_3d,cice_3d,
     1                          thresh_cvr_cty_vv,thresh_cvr_lwc,
     1                          l_flag_cloud_type,cldpcp_type_3d,         ! I/O
     1                          l_flag_mvd,mvd_3d,
     1                          l_flag_icing_index,icing_index_3d,
     1                          vv_to_height_ratio_Cu,                    ! I
     1                          vv_to_height_ratio_Sc,                    ! I
     1                          vv_for_St,                                ! I
     1                          l_flag_bogus_w,omega_3d,l_bogus_radar_w,
     1                          l_deep_vv,                                ! I
     1                          twet_snow,                                ! I
     1                          l_flag_pcp_type,                          ! I
     1                          istatus)                                  ! O

!       Steve Albers
cdoc    This routine calculates SLWC, Cloud Type, MVD, and Icing Index
cdoc    This routine also does the Cloud Bogussed Omega and the Snow Potential.
!       These have been combined to make more efficient use of looping through
!       large arrays.
!       The input flags can be adjusted to allow only some of these to be returned
!
!       Steve Albers Oct 1992 Add Haines cloud-ice modification
!                    Apr 1994 Remove byte usage
!                    Dec 1996 Clean up error handling
!                    Jun 2002 Added dx to calling arguments for
!                             cloud_bogus_w

!       LWC Flags
!       iflag_slwc =  1 : Adiabatic LWC
!       iflag_slwc =  2 : Adjusted  LWC
!       iflag_slwc =  3 : Adjusted  SLWC
!       iflag_slwc = 11 : Smith-Feddes LWC
!       iflag_slwc = 12 : Smith-Feddes SLWC
!       iflag_slwc = 13 : New Smith-Feddes LWC
!       iflag_slwc = 14 : New Smith-Feddes SLWC

        real temp_3d(ni,nj,nk)         ! Input
        real rh_3d_pct(ni,nj,nk)       ! Input
        real heights_3d(ni,nj,nk)      ! Input
        real pres_3d(ni,nj,nk)         ! Input
        real radar_3d(ni,nj,nk)        ! Input

        real omega_3d(ni,nj,nk)        ! Input / Output (if l_flag_bogus_w = .true.)
                                         !                (Omega is in Pa/S)
        real slwc_3d(ni,nj,nk)         ! Output
        real cice_3d(ni,nj,nk)         ! Output
        integer cldpcp_type_3d(ni,nj,nk) ! Output (1st 4bits pcp, 2nd 4 are cld)

        real mvd_3d(ni,nj,nk)          ! Output
        integer icing_index_3d(ni,nj,nk) ! Output
!       real lwc_res_3d(ni,nj,nk)      ! Output

!       real snow_2d(ni,nj)            ! Output

        real temp_1d(nk)               ! Local
        real slwc_1d(nk)
        real cice_1d(nk)
        real heights_1d(nk)
        real pressures_mb(nk)
        real pressures_pa(nk)
        real rlwc_laps1d(nk)
        real w_1d(nk)                  ! Units are M/S
        real lwc_res_1d(nk)
        real prob_laps(nk)
        real d_thetae_dz_1d(nk)

        integer cloud_type_1d(nk)

        integer iarg

        integer iflag_slwc
        logical   l_flag_cloud_type,l_flag_mvd,l_flag_icing_index
     1           ,l_flag_bogus_w,l_cloud,l_bogus_radar_w,l_deep_vv
     1           ,l_flag_pcp_type

      ! Used for "Potential" Precip Type
        logical l_mask_pcptype(ni,nj)
        integer ibase_array_lwc(ni,nj)
        integer itop_array_lwc(ni,nj)
        integer ibase_array_cty_vv(ni,nj)
        integer itop_array_cty_vv(ni,nj)


!       EXTERNAL        LIB$INIT_TIMER,
!    1                  LIB$SHOW_TIMER,
!    1                  my_show_timer

        include 'laps_cloud.inc'

        real clouds_3d(ni,nj,KCLOUD)

        Integer  KCLOUD_M1,  KCLOUD_P1
        parameter (KCLOUD_M1 = KCLOUD - 1)
        parameter (KCLOUD_P1 = KCLOUD + 1)

        real thresh_cvr_cty_vv,thresh_cvr_lwc
!       parameter (THRESH_CVR = 0.65)
!       parameter (THRESH_CVR = 0.75)

        character*2 c2_type

        Real vv_to_height_ratio_Cu
        Real vv_to_height_ratio_Sc
        Real vv_for_St

        write(6,*)' Start LWC/Omega/Snow Potential Routine'

!       thresh_cvr_cty_vv = thresh_cvr
!       thresh_cvr_lwc = thresh_cvr

        zero = 0.

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)then
            write(6,*)' Error reading r_missing_data'
            stop
        endif

!       Check consistency of input flags
        if(l_flag_icing_index)then
            iflag_slwc = 13
            l_flag_mvd = .true.
        endif

        if(l_flag_mvd .or. l_flag_bogus_w
!       1             .or. iflag_slwc .ge. 10
     1                .or. iflag_slwc .ge. 1
     1                                          )then
            l_flag_cloud_type = .true.
        endif

        if(l_flag_bogus_w)then
            write(6,*)' Computing Cloud Bogused Omega'
        endif

        write(6,*)' Generating Lowest Base and Highest Top Arrays'
        ibase_array_lwc = KCLOUD_P1
        itop_array_lwc = 0
        ibase_array_cty_vv = KCLOUD_P1
        itop_array_cty_vv = 0

        do k = KCLOUD,1,-1
            do j = 1,nj
            do i = 1,ni
                if(clouds_3d(i,j,k) .ge. THRESH_CVR_CTY_VV)then
                    ibase_array_cty_vv(i,j) = k
                endif
                if(clouds_3d(i,j,k) .ge. THRESH_CVR_LWC)then
                    ibase_array_lwc(i,j) = k
                endif
            enddo
            enddo
        enddo

        do k = 1,KCLOUD
            do j = 1,nj
            do i = 1,ni
                if(clouds_3d(i,j,k) .ge. THRESH_CVR_CTY_VV)then
                    itop_array_cty_vv(i,j) = k
                endif
                if(clouds_3d(i,j,k) .ge. THRESH_CVR_LWC)then
                    itop_array_lwc(i,j) = k
                endif
            enddo
            enddo
        enddo

        I4_elapsed = ishow_timer()
        write(6,*)

        if(iflag_slwc .ne. 0)then
            write(6,*)' Initializing SLWC/CICE arrays, iflag_slwc = '
     1               ,iflag_slwc
            slwc_3d = zero
            cice_3d = zero
        endif

        if(l_flag_bogus_w)then
            write(6,*)' Initializing Omega array'
            omega_3d = r_missing_data
        endif

        if(l_flag_mvd)then
            write(6,*)' Initializing MVD array'
            mvd_3d = zero
        endif

        if(l_flag_cloud_type)then
            write(6,*)' Initializing Cloud Type Array'
            cldpcp_type_3d = 0
        endif

        if(l_flag_icing_index)then
            write(6,*)' Initializing Icing Index Array'
            icing_index_3d = 0
        endif

        I4_elapsed = ishow_timer()

        n_cloud_columns_cty_vv = 0
        n_cloud_columns_slwc = 0
        n_slwc_call = 0

        write(6,*)' Finding Cloud Layers and Computing Output Field(s)'
        do i = 1,ni
        do j = 1,nj

            l_cloud = .false.

!           Generate vertical sounding at this grid point and call SLWC routine

            if(ibase_array_cty_vv(i,j) .ne. KCLOUD_P1)then ! At least one layer exists
                l_cloud = .true.
                n_cloud_columns_cty_vv = n_cloud_columns_cty_vv + 1

                do k = 1,nk ! Initialize
                    temp_1d(k) = temp_3d(i,j,k)
                    heights_1d(k) = heights_3d(i,j,k)
                    pressures_mb(k) = pres_3d(i,j,k) / 100.
                    pressures_pa(k) = pres_3d(i,j,k)
                    cloud_type_1d(k) = 0
                enddo

!               cloud_ceil(i,j) = 3000. ! Dummied in for testing
!               cloud_top(i,j)  = 5000. ! Dummied in for testing

!               Get Base and Top
                k = max(ibase_array_cty_vv(i,j) - 1,1)
                k_highest = itop_array_cty_vv(i,j) - 1

!               First time around has lower threshold
!                                     keep stability & cloud type
!                                     exclude SLWC 
!                                     keep MVD
                                      
                do while (k .le. k_highest)
                  if(clouds_3d(i,j,k+1) .ge. THRESH_CVR_CTY_VV .and.
     1               clouds_3d(i,j,k  ) .lt. THRESH_CVR_CTY_VV
     1                              .OR.
     1             k .eq. 1 .and. clouds_3d(i,j,k).ge.THRESH_CVR_CTY_VV      
     1                                                  )then
                    cld_base_m = 0.5 * (cld_hts(k) + cld_hts(k+1))
                    k_base = k + 1

c                   if(i .eq. 1)write(6,*)i,j,k,' Cloud Base'

                    k = k + 1

                    do while (k .le. kcloud-1)
                      if(clouds_3d(i,j,k  ) .gt. THRESH_CVR_CTY_VV .and.
     1                   clouds_3d(i,j,k+1) .le. THRESH_CVR_CTY_VV)then       
                        cld_top_m = 0.5 * (cld_hts(k) + cld_hts(k+1))

!                       Constrain cloud top to top of domain
                        cld_top_m = min(cld_top_m,heights_3d(i,j,nk))

                        k_top = k

!                       We have now defined a cloud base and top
                        k_1d_base = int(height_to_zcoord3(
     1                              cld_base_m,heights_3d
     1                         ,pressures_pa,ni,nj,nk,i,j,istatus)) + 1       
                        k_1d_top  = int(height_to_zcoord3(
     1                              cld_top_m ,heights_3d
     1                         ,pressures_pa,ni,nj,nk,i,j,istatus))

                        if(istatus .ne. 1)then
                            write(6,*)' Returning from '
     1                       ,'get_cloud_deriv (height_to_zcoord3 call)'  
     1                       ,cld_base_m,cld_top_m
     1                       ,(heights_3d(i,j,kk),kk=1,nk)
                            return
                        endif

!                       Make sure cloud base and top stay in LAPS domain
                        k_1d_base = min(k_1d_base,nk)
                        k_1d_top  = min(k_1d_top ,nk)

c                       if(i .eq. 1)write(6,*)i,j,k,' Cloud Top',k_base,k_top

                        if(l_flag_cloud_type)then ! Get 1d Stability & Cloud type
                            call get_stability(nk,temp_1d,heights_1d
     1                                         ,pressures_mb,k_1d_base
     1                                         ,k_1d_top,d_thetae_dz_1d)

                            do k_1d = k_1d_base,k_1d_top
                                call get_cloudtype(temp_1d(k_1d)
     1                           ,d_thetae_dz_1d(k_1d),cld_base_m
     1                           ,cld_top_m,itype,c2_type)

                                if(radar_3d(i,j,k_1d) .gt. 45.
     1                       .and. radar_3d(i,j,k_1d) .ne. 
     1                                              r_missing_data)then       
                                    itype = 10 ! CB
                                endif

                                cloud_type_1d(k_1d)     = itype

                                cldpcp_type_3d(i,j,k_1d) = itype
                            enddo
                        endif

!                       We have now defined a cloud base and top
                        do k_1d = k_1d_base,k_1d_top ! Loop through the cloud layer
                            if(l_flag_mvd)then
                                call get_mvd(cloud_type_1d(k_1d)
     1                                      ,rmvd_microns)
                                mvd_3d(i,j,k_1d) = rmvd_microns * 1e-6
                            endif

                        enddo ! k_1d

                        goto1000

                      endif ! Found Cloud Top

                      k = k + 1

                    enddo ! k

                  endif ! Found Cloud Base

1000              k = k + 1

                enddo ! k (cloud layer loop)

!               Get Base and Top
                k = max(ibase_array_lwc(i,j) - 1,1)
                k_highest = itop_array_lwc(i,j) - 1

!               Second time around has higher threshold
!                   remove redundant stability & cloud type
!                   keep SLWC
!                   remove MVD
                do while (k .le. k_highest)
                  if(clouds_3d(i,j,k+1) .ge. THRESH_CVR_LWC .and.
     1               clouds_3d(i,j,k  ) .lt. THRESH_CVR_LWC
     1                              .OR.
     1             k .eq. 1 .and. clouds_3d(i,j,k) .ge. THRESH_CVR_LWC
     1                                                  )then
                    cld_base_m = 0.5 * (cld_hts(k) + cld_hts(k+1))
                    k_base = k + 1

c                   if(i .eq. 1)write(6,*)i,j,k,' Cloud Base'

                    k = k + 1

                    do while (k .le. kcloud-1)
                      if(clouds_3d(i,j,k  ) .gt. THRESH_CVR_LWC .and.
     1                   clouds_3d(i,j,k+1) .le. THRESH_CVR_LWC)then
                        cld_top_m = 0.5 * (cld_hts(k) + cld_hts(k+1))

!                       Constrain cloud top to top of domain
                        cld_top_m = min(cld_top_m,heights_3d(i,j,nk))

                        k_top = k

!                       We have now defined a cloud base and top
                        k_1d_base = int(height_to_zcoord3(
     1                              cld_base_m,heights_3d
     1                         ,pressures_pa,ni,nj,nk,i,j,istatus)) + 1       
                        k_1d_top  = int(height_to_zcoord3(
     1                              cld_top_m ,heights_3d
     1                         ,pressures_pa,ni,nj,nk,i,j,istatus))

                        if(istatus .ne. 1)then
                            write(6,*)' Returning from '
     1                       ,'get_cloud_deriv (height_to_zcoord3 call)'  
     1                       ,cld_base_m,cld_top_m
     1                       ,(heights_3d(i,j,kk),kk=1,nk)
                            return
                        endif

!                       Make sure cloud base and top stay in LAPS domain
                        k_1d_base = min(k_1d_base,nk)
                        k_1d_top  = min(k_1d_top ,nk)

c                       if(i .eq. 1)write(6,*)i,j,k,' Cloud Top',k_base,k_top

!                       We have now defined a cloud base and top
                        if(iflag_slwc .ne. 0)then

                          n_slwc_call = n_slwc_call + 1

                          if(iflag_slwc .lt. 10)then
                            call get_slwc1d(nk,cld_base_m,cld_top_m
     1                                     ,k_1d_base,k_1d_top
     1                                     ,heights_1d,temp_1d
     1                                     ,pressures_pa
     1                                     ,iflag_slwc,zero,slwc_1d)       
                          else ! Get Smith-Feddes Output
                            mode = 1

                            do k_1d = 1,nk ! Initialize
                              slwc_1d(k_1d) = zero
                              cice_1d(k_1d) = zero
                              prob_laps(k_1d) = zero
                            enddo

!                           QC the data going into SMF
                            if(cld_top_m .gt. heights_1d(nk) - 110.)then
                                cld_top_qc_m = heights_1d(nk) - 110.
                                cld_base_qc_m =
     1                          min(cld_base_m,cld_top_qc_m - 110.)
                                write(6,501)i,j,nint(heights_1d(nk))
     1                          ,nint(cld_base_qc_m),nint(cld_top_qc_m)
501                             format(1x,'QC data for SMF, ht/bs/tp'
     1                                ,2i4,3i6)

                            else ! normal case
                                cld_top_qc_m = cld_top_m
                                cld_base_qc_m = cld_base_m

                            endif

                            if(iflag_slwc .le. 12)then
!                               Performance is 200/350 calls/sec on PROFS3/FSL9
                                call get_smf_1d(nk,cld_base_qc_m
     1                              ,cld_top_qc_m
!    1                              ,cloud_type_1d((k_1d_base + k_1d_top)/2)
     1                              ,1
     1                              ,heights_1d,pressures_mb,temp_1d
     1                              ,slwc_1d,prob_laps,mode)
                            else ! iflag_slwc = 13 or 14

!                               Performance is 450/700 calls/sec on PROFS3/FSL9
                                call get_sfm_1d(nk,cld_base_qc_m
     1                          ,cld_top_qc_m
     1                          ,cloud_type_1d((k_1d_base + k_1d_top)/2)
!    1                          ,1
     1                          ,heights_1d,pressures_mb,temp_1d
     1                          ,slwc_1d,cice_1d,prob_laps,mode)
!    1                          ,slwc_1d,prob_laps,mode)
                            endif


                            if(iflag_slwc .eq. 12 .or.
     1                         iflag_slwc .eq. 14           )then
                              do k_1d = k_1d_base,k_1d_top
                                if(temp_1d(k_1d) .gt. 273.15)then
                                  slwc_1d(k_1d) = 0.
                                endif
                              enddo ! k_1d
                            endif

                          endif ! iflag <> 10

                        endif ! iflag .ne. 0


                        do k_1d = k_1d_base,k_1d_top ! Loop through the cloud layer
                            if(iflag_slwc .ne. 0)then
                              if(slwc_1d(k_1d) .gt. 0.)then
                                if(istat_radar .eq. 1
     1                                     .AND. 
     1                             temp_3d(i,j,k_1d) .lt. 273.15
     1                                            )then ! Apply Radar Depletion

!                               if(.false.)then ! Don't Apply Radar Depletion

                                 if(cloud_type_1d(k_1d) .ne. 10)then ! Not CB

                                  depl1 = 25. 
                                  depl2 = 30.

                                  if(radar_3d(i,j,k_1d) .le. depl1)then ! No Depletion
                                    continue

                                  elseif(radar_3d(i,j,k_1d) .gt. depl2 
     1                                                  )then ! Total Depletion
                                    cice_1d(k_1d) = cice_1d(k_1d) 
     1                                            + slwc_1d(k_1d)
                                    slwc_1d(k_1d) = zero

                                  else ! Ramped Depletion
                                    ramp = 1.0 - 
     1                                  ((radar_3d(i,j,k_1d) - depl1) 
     1                                                / (depl2-depl1))       
                                    cice_1d(k_1d) = cice_1d(k_1d) 
     1                                       + slwc_1d(k_1d) * (1.-ramp)       
                                    slwc_1d(k_1d) = slwc_1d(k_1d) * ramp

                                  endif ! Dbz Sufficient for Depletion

                                 endif ! Not a CB

                                endif ! Valid Radar Data and below freezing

                              endif ! LWC Returned by Algorithm

!                             Pull lwc/ice from 1d array to 3d array
                              if(slwc_1d(k_1d) .gt. 0.)
     1                           slwc_3d(i,j,k_1d) = slwc_1d(k_1d)
                              if(cice_1d(k_1d) .gt. 0.)
     1                           cice_3d(i,j,k_1d) = cice_1d(k_1d)

                            endif ! iflag_slwc

                        enddo ! k_1d

                        goto2000

                      endif ! Found Cloud Top

                      k = k + 1

                    enddo ! k

                  endif ! Found Cloud Base

2000              k = k + 1

                enddo ! k (cloud layer loop)

            endif ! ibase_array_cty_vv (At least one layer exists)

            if(l_flag_bogus_w)then
              if(l_cloud)then
                  call cloud_bogus_w
     1            (grid_spacing_cen_m,cloud_type_1d,heights_1d,nk  ! I
     1                           ,vv_to_height_ratio_Cu            ! I
     1                           ,vv_to_height_ratio_Sc            ! I
     1                           ,vv_for_St                        ! I
     1                           ,l_deep_vv                        ! I
     1                           ,w_1d)                            ! O

                  do k = 1,nk ! Transfer the 1D w (M/S) into the output 
                              ! omega (Pa/S) array
                      if(w_1d(k) .ne. r_missing_data)then
                          omega_3d(i,j,k) = w_to_omega(w_1d(k)
     1                                       ,pressures_pa(k))
!                     else
!                         w_1d(k) = 0. ! Condition the array for the snow_diag
                      endif
                  enddo

              endif ! l_cloud
            endif ! l_flag_bogus_w

        enddo ! j
        enddo ! i

        write(6,*)
     1    ' N_CLOUD_COLUMNS_CTY_VV,N_CLOUD_COLUMNS_SLWC,N_SLWC_CALL = '
     1                                       ,n_cloud_columns_cty_vv     
     1                                       ,n_cloud_columns_slwc
     1                                       ,n_slwc_call

        I4_elapsed = ishow_timer()

!       Simulate Radar Data
        istat_radar = 1
        l_mask_pcptype = .false.

        if(istat_radar .eq. 1)then ! Get 3D precip type

            write(6,*)' Computing 3D Precip type'

            call cpt_pcp_type_3d(temp_3d,rh_3d_pct,pres_3d
     1                  ,radar_3d,l_mask_pcptype,grid_spacing_cen_m       
     1                  ,ni,nj,nk,twet_snow,cldpcp_type_3d,istatus)
            if(istatus .ne. 1)then
                write(6,*)'Bad status returned from cpt_pcp_type_3d'
                return
            endif

            I4_elapsed = ishow_timer()

        endif ! Valid Radar Data

        if(l_flag_icing_index)then

            write(6,*)' Computing Icing Severity Index field'

            do i = 1,ni
            do j = 1,nj
            do k = 1,nk
                iarg = cldpcp_type_3d(i,j,k)
                i_precip_type = iarg/16
                i_cloud_type  = iarg - i_precip_type*16

                if(iarg .gt. 0)then ! clouds or precip are present
                    call isi3(slwc_3d(i,j,k),temp_3d(i,j,k)-273.15
     1                          ,i_cloud_type,i_precip_type,xindex)
                    iarg = xindex
                    icing_index_3d(i,j,k) = iarg

                else
                    icing_index_3d(i,j,k) = 0

                endif

            enddo ! k
            enddo ! j
            enddo ! i

            I4_elapsed = ishow_timer()

        endif ! iflag icing index

        istatus = 1
        return

        end


        subroutine cpt_pcp_type_3d(temp_3d,rh_3d_pct,pres_3d
     1  ,radar_3d,l_mask,grid_spacing_cen_m
     1  ,ni,nj,nk,twet_snow,cldpcp_type_3d,istatus)

!       1991    Steve Albers
!       1997    Steve Albers - Allow for supercooled precip generation
!                            - Misc streamlining of logic

cdoc    Compute 3D Precip Type given profiles of T, RH, Reflectivity
!       This program modifies the most significant 4 bits of the integer
!       array by inserting multiples of 16.

        real temp_3d(ni,nj,nk)                          ! Input
        real rh_3d_pct(ni,nj,nk)                        ! Input
        real pres_3d(ni,nj,nk)                          ! Input
        integer cldpcp_type_3d(ni,nj,nk)                  ! Output
        real radar_3d(ni,nj,nk)
        integer itype
        logical l_mask(ni,nj) ! Used for "Potential" Precip Type

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return
       
        call get_ref_base_useable(ref_base_useable,istatus)
        if(istatus .ne. 1)return

!       Stuff precip type into cloud type array
!       0 - No Precip
!       1 - Rain
!       2 - Snow
!       3 - Freezing Rain
!       4 - Sleet
!       5 - Hail / Graupel

        zero_c = 273.15
        rlayer_refreez_max = 0.0

!       Ramp the hail_ref_thresh between 0-10km grid spacing
        arg = min(grid_spacing_cen_m,10000.)
        hail_ref_thresh1 = 55. - arg/1000.

        n_zr = 0
        n_sl = 0
        n_last = 0

        do j = 1,nj
        do i = 1,ni

            iflag_melt = 0       ! Equivalent to iprecip_type = liquid while
                                 ! t_wb_c > 0 (present or past), depending on 
                                 ! where we are in the loop.
            iflag_refreez = 0
            rlayer_refreez = 0.0

            iprecip_type_last = 0

            do k = nk,1,-1

                if( (radar_3d(i,j,k) .ge. ref_base_useable .and.  ! Above Noise
     1               radar_3d(i,j,k) .ge. 0.               .and.  ! Sig Echo
     1               radar_3d(i,j,k) .ne. r_missing_data         )
     1                               .OR. l_mask(i,j)            )then        

!                   Set refreezing flag
                    t_c         = temp_3d(i,j,k) - zero_c
                    td_c        = DWPT(t_c,rh_3d_pct(i,j,k))
                    pressure_mb = pres_3d(i,j,k) / 100.

                    thresh_melt_c =
     1               wb_melting_threshold(t_c,radar_3d(i,j,k),twet_snow)

!                   This function call here is fast but returns a t_wb_c
!                   equal to t_c if pres < 500mb. This approximation should
!                   not hurt the algorithm. Same if abs(t_c) > 60.

                    t_wb_c = twet_fast(t_c,td_c,pressure_mb)

                    if(t_wb_c .lt. 0.    .and. 
     1                 iflag_melt .eq. 1)then ! integrate refreezing eq.

!                       Integrate below freezing temperature times column
!                       thickness - ONLY for portion of layer below freezing

                        temp_lower_c = t_wb_c
                        k_upper = min(k+1,nk)

!                       For simplicity and efficiency, the assumption is
!                       here made that the wet bulb depression is constant
!                       throughout the level.

                        temp_upper_c = t_wb_c +
     1                  (temp_3d(i,j,k_upper) - temp_3d(i,j,k))

                        if(temp_upper_c .le. 0.)then
                            frac_below_zero = 1.0
                            tbar_c = 0.5 * (temp_lower_c + temp_upper_c)

                        else ! Layer straddles the freezing level
                            frac_below_zero = temp_lower_c
     1                                   / (temp_lower_c - temp_upper_c)
                            tbar_c = 0.5 * temp_lower_c

                        endif

                        pressure_interval_pa = pres_3d(i,j,k) 
     1                                       - pres_3d(i,j,k+1)

                        rlayer_refreez = rlayer_refreez
     1            + abs(tbar_c * pressure_interval_pa * frac_below_zero)     

!                       if(rlayer_refreez .ge. 75000.)then
                        if(rlayer_refreez .ge. 25000.)then
                            iflag_refreez = 1
                        endif

                        rlayer_refreez_max =
     1                          max(rlayer_refreez_max,rlayer_refreez)

                    else ! Temp > 0C or iflag_melt = 0
                        iflag_refreez = 0
                        rlayer_refreez = 0.0

                    endif ! Temp is freezing

!                   Calculate Hail/Graupel threshold with temperature factored in
                    h_thr_depr = max(10. - abs(-5. - temp_3d(i,j,k)),0.)         
                    hail_ref_thresh = hail_ref_thresh1 - 3. * h_thr_depr
                    hail_ref_thresh = max(hail_ref_thresh,30.)                    

                    if(radar_3d(i,j,k) .lt. hail_ref_thresh)then ! Not Hail
                        if(t_wb_c .ge. thresh_melt_c)then     ! Warm, got rain
                            iprecip_type = 1
                            iflag_melt = 1

                        elseif(t_wb_c .ge. 0.0)then ! Between 0C and Melt threshold
                            if(iprecip_type_last .eq. 0       ! Generating lyr
     1                    .or. iprecip_type_last .eq. 3       ! Freezing Rain
     1                                                 )then  ! Set to Rain
                                iprecip_type = 1              
                                iflag_melt = 1

                            else                              ! Unchanged pcp
                                iprecip_type = iprecip_type_last
                                n_last = n_last + 1
                                if(n_last .lt. 5)then
                                    write(6,*)'Unchanged Precip'
     1                                        ,i,j,k,t_wb_c
                                endif
                            endif

                        else ! below 0C (Freezing Precip or Snow)
                            if(iprecip_type_last .eq. 0)then  ! Generating lyr
!                               if(t_wb_c .ge. -6.)then       ! Supercooled pcp
                                if(.false.)then               ! Supercooled pcp
                                    iflag_melt = 0
                                    iprecip_type = 3          ! Freezing Rain
 
                                else
                                    iprecip_type = 2          ! Snow

                                endif

                            elseif(iprecip_type_last .eq. 1 .or.
     1                             iprecip_type_last .eq. 3
     1                                                 )then  ! Liquid
                                                              ! Now Supercooled

                              ! Is it freezing rain or sleet?
                                if(iflag_refreez .eq. 0)then  ! Freezing Rain
                                    n_zr = n_zr + 1
                                    if(n_zr .lt. 30)then
                                        write(6,5)i,j,k,t_wb_c
     1                                  ,temp_3d(i,j,k),rh_3d_pct(i,j,k) 
5                                       format('zr',3i3,2f8.2,f8.1)
                                    endif
                                    iprecip_type = 3

                                else  ! (iflag_refreez = 1)      ! Sleet
                                    n_sl = n_sl + 1
                                    iprecip_type = 4
                                    iflag_melt = 0               ! Reset flags
                                    iflag_refreez = 0
                                    rlayer_refreez = 0.0

                                endif

                            elseif(iprecip_type_last .eq. 5)then ! Hail
                                iprecip_type = 2                 ! Snow

                            else                                 ! Not Hail
                                iprecip_type = iprecip_type_last ! Unchanged

                            endif ! iflag_melt = 1

                        endif ! t_wb_c

                    else ! >= hail_ref_thresh
                        iprecip_type = 5                         ! Hail

                    endif ! Intense enough for hail

                else ! No Radar Echo
                    iprecip_type = 0
                    iflag_melt = 0                               ! Reset
                    iflag_refreez = 0   
                    rlayer_refreez = 0.

                endif ! Radar Echo?

!               Insert most sig 4 bits into array
                itype = cldpcp_type_3d(i,j,k)
                itype = itype - (itype/16)*16     ! Initialize the 4 bits
                itype = itype + iprecip_type * 16 ! Add in the new value
                cldpcp_type_3d(i,j,k) = itype

                iprecip_type_last = iprecip_type

            enddo ! k
        enddo ! j
        enddo ! i

        write(6,*)' rlayer_refreez_max = ',rlayer_refreez_max
        write(6,*)' n_zr/n_sl = ',n_zr,n_sl

        istatus = 1
        return
        end


        subroutine get_sfc_preciptype(pres_2d,t_sfc_k,td_sfc_k
     1            ,cldpcp_type_3d,twet_snow
     1            ,dbz_2d,pcp_type_2d,ni,nj,nk)

!       Steve Albers 1991

cdoc    Compute Sfc Precip Type, given both sfc and 3D fields

        real pres_2d(ni,nj)             ! Input
        real t_sfc_k(ni,nj)             ! Input
        real td_sfc_k(ni,nj)            ! Input
        integer cldpcp_type_3d(ni,nj,nk)  ! Input
        real dbz_2d(ni,nj)              ! Input (Low Level reflectivity)
        integer pcp_type_2d(ni,nj)        ! Output
                                       ! Leftmost 4 bits contain the precip type

        integer iarg

        n_precip = 0
        n_chg_frz = 0
        n_chg_mlt = 0
        n_hail = 0

        do j = 1,nj
        do i = 1,ni

!           Interpolate from three dimensional grid to terrain surface
            rksfc = zcoord_of_pressure(pres_2d(i,j))

!           If no precip at the NEAREST grid point, check the point just above
!           the surface
            ksfc_close = max(nint(rksfc),1)
            if(float(ksfc_close) .lt. rksfc)then
                ksfc_above = ksfc_close + 1
            else
                ksfc_above = ksfc_close
            endif

!           Pull out the precip type
            iarg = cldpcp_type_3d(i,j,ksfc_close)
            iprecip_type = iarg/16

            if(iprecip_type .ne. 0)then
                pcp_type_2d(i,j) = cldpcp_type_3d(i,j,ksfc_close)
            else
                pcp_type_2d(i,j) = cldpcp_type_3d(i,j,ksfc_above)
                iarg = cldpcp_type_3d(i,j,ksfc_above)
                iprecip_type = iarg/16
            endif

!           0 - No Precip
!           1 - Rain
!           2 - Snow
!           3 - Freezing Rain
!           4 - Sleet
!           5 - Hail

            zero_c = 273.15

!           Given no hail, if sfc wetbulb is below freezing, change the
!           precip type to the lowest frozen type in the 3D column.
!           If sfc wetbulb is above freezing change the precip type to rain

            if(iprecip_type .ne. 0 .and. iprecip_type .ne. 5)then ! Not Hail

                t_sfc_c  = t_sfc_k(i,j)  - 273.15
                td_sfc_c = td_sfc_k(i,j) - 273.15
                pres_mb  = pres_2d(i,j) / 100.

                tw_sfc_c = twet_fast(t_sfc_c,td_sfc_c,pres_mb)

!               Note that the dbz value has not been passed in yet
                thresh_melt_c = 
     1              wb_melting_threshold(t_sfc_c,dbz_2d(i,j),twet_snow)       

                n_precip = n_precip + 1
                if(iprecip_type .ne. 1)then       ! Not Rain
                    if(tw_sfc_c .gt. thresh_melt_c)then ! Above freezing
                        iprecip_type = 1
                        n_chg_mlt = n_chg_mlt + 1
                    endif

                elseif(iprecip_type .eq. 1)then   ! Rain
                    if(tw_sfc_c .lt. thresh_melt_c)then ! Below freezing

!                       Get lowest frozen precip type
                        iprecip_frozen = 0
                        do k = nk,ksfc_close,-1
                            iarg = cldpcp_type_3d(i,j,k)
                            iprecip = iarg/16

                            if(iprecip .eq. 2 .or. iprecip .eq. 3 .or.
     1                                             iprecip .eq. 4)then
                                iprecip_frozen = iprecip
                            endif

                        enddo ! k
                        iprecip_type = iprecip_frozen
                        n_chg_frz = n_chg_frz + 1

                    endif ! Below melting thresh at sfc

                endif ! 3D Precip is rain

!               Ensure that ZR is diagnosed only if sfc dry bulb < 0C
                if(iprecip_type .eq. 3 .and. t_sfc_c .gt. 0.)then
                    iprecip_type = 1
                    n_chg_mlt = n_chg_mlt + 1
                endif

                iarg = iprecip_type * 16
                pcp_type_2d(i,j) = iarg

            endif ! Precip

            if(iprecip_type .eq. 5)n_hail = n_hail + 1

        enddo ! i
        enddo ! j

        write(6,*)' # Sfc pts with precip =                   ',n_precip       
        write(6,*)' # Sfc pts changed to frozen by sfc data = '
     1                                                        ,n_chg_frz
        write(6,*)' # Sfc pts changed to rain by sfc data  =  '
     1                                                        ,n_chg_mlt      
        write(6,*)' # Sfc pts with hail =                     ',n_hail

        return
        end


        subroutine get_slwc1d(nk,cbase_m,ctop_m,kbase,ktop
     1                       ,heights_1d,temp_1d,pressures_pa
     1                       ,iflag_slwc,zero,slwc_1d)

cdoc    Determine cloud liquid profile within a given cloud layer

        real temp_1d(nk)
        real heights_1d(nk)
        real pressures_pa(nk)
        real slwc_1d(nk)  ! Output

        do k = 1,nk ! Initialize
            slwc_1d(k) = zero
        enddo

        if(ctop_m .gt. cbase_m)then

!           Determine Lowest and Highest Grid Points within the cloud
            if(ktop .ge. kbase .and. kbase .ge. 2)then

!               Get cloud base pressure and temperature
                cbase_pa = height_to_pressure(cbase_m,heights_1d
     1                          ,pressures_pa,1,1,nk,1,1)

!               cbase_pa = ztopsa(cbase_m) * 100.

                frac_k = (cbase_m - heights_1d(kbase-1))  /
     1         (heights_1d(kbase) - heights_1d(kbase-1))

                cbase_k = temp_1d(kbase-1) * (1.0 - frac_k)
     1                  + temp_1d(kbase)   * frac_k


!               Get cloud top temperature
                frac_k = (ctop_m - heights_1d(ktop-1))  /
     1         (heights_1d(ktop) - heights_1d(ktop-1))

                ctop_k = temp_1d(ktop-1) * (1.0 - frac_k)
     1                 + temp_1d(ktop)   * frac_k


!               Calculate SLWC at each vertical grid point. For each level
!               we use an assumed cloud extending from the actual cloud base
!               to the height of the grid point in question.

                do k=kbase,ktop
                    grid_top_pa = zcoord_of_level(k)

                    grid_top_k = temp_1d(k)

                    call laps_slwc_revb(cbase_pa,cbase_k
     1                ,grid_top_pa,grid_top_k,ctop_k
     1                ,ADIABATIC_LWC,ADJUSTED_LWC,ADJUSTED_SLWC
     1                ,I_STATUS1,I_STATUS2)

!                   These three lines are dummied in
!                   i_status = .true.
!                   slwc_adiabatic = -(ctop_k-cbase_k) / 10.
!                   slwc_adjusted =  -(ctop_k-cbase_k) / 10.

                    if(i_status2 .eq. 1)then
                        if(iflag_slwc .eq. 1)then
                            slwc_1d(k) = ADIABATIC_LWC
                        elseif(iflag_slwc .eq. 2)then
                            slwc_1d(k) = ADJUSTED_LWC
                        elseif(iflag_slwc .eq. 3)then
                            slwc_1d(k) = ADJUSTED_SLWC
                        endif
                    else
                        write(6,*)' Error Detected in SLWC'
                    endif
                enddo ! k
            endif ! thick enough cloud exists
        endif ! cloud exists

        return
        end

        subroutine integrate_slwc(slwc,heights_3d,imax,jmax,kmax
     1                           ,slwc_int)            

cdoc    Integrates cloud liquid through the column

        real slwc(imax,jmax,kmax)       ! Input in g/m**3
        real heights_3d(imax,jmax,kmax) ! Input 
        real slwc_int(imax,jmax)  ! Output in m (metric tons of water / m**2)
        real depth(kmax)          ! Local

        do j = 1,jmax
        do i = 1,imax
            slwc_int(i,j) = 0.
            do k = 1,kmax-1
                depth(k) = heights_3d(i,j,k+1) - heights_3d(i,j,k) ! meters
                                                         ! convert units
                slwc_ave = (slwc(i,j,k) + slwc(i,j,k+1)) * 0.5 * 1e-6 
                slwc_int(i,j) = slwc_int(i,j) + slwc_ave * depth(k)
            enddo ! k
        enddo ! i
        enddo ! j

        return
        end

        subroutine get_mvd(itype,rmvd)

cdoc    rmvd        Output Mean Volume Diameter in microns, given cloud type

!                                                 ! no cloud
        IF(itype.EQ.0) THEN
            rmvd = 0.0
!                                                 ! St
        ELSEIF (itype.EQ. 1) THEN
            rmvd = 12.0
!                                                 ! Sc
        ELSEIF (itype.EQ. 2) THEN
            rmvd = 10.0
!                                                 ! Cu
        ELSEIF (itype.EQ. 3) THEN
            rmvd = 18.0
!                                                 ! Ns
        ELSEIF (itype.EQ. 4) THEN
            rmvd = 12.0
!                                                 ! Ac
        ELSEIF (itype.EQ. 5) THEN
            rmvd = 18.0
!                                                 ! As
        ELSEIF (itype.EQ. 6) THEN
            rmvd = 10.0
!                                                 ! Cb
        ELSEIF (itype.EQ.10) THEN
            rmvd = 25.0
!                                                 ! cirriform
        ELSE ! itype = 7,8,9
            rmvd = 10.0
!
        ENDIF

        return
        end


        subroutine get_cloudtype(temp_k,d_thetae_dz
     1                          ,cbase_m,ctop_m,itype,c2_type)

cdoc    Determine Cloud Type, given temperature and stability d(theta(e))/dz

        character*2 c2_type,c2_cloudtypes(10)

        data c2_cloudtypes
     1  /'St','Sc','Cu','Ns','Ac','As','Cs','Ci','Cc','Cb'/

        temp_c = temp_k - 273.15
        depth_m = ctop_m - cbase_m

!       Go from Stability to Cloud Type
        If ( temp_c.ge.-10.) then

            if (d_thetae_dz.ge. +.001) then
                itype = 1 ! St
            elseif (d_thetae_dz .lt. +.001 .and. d_thetae_dz .ge. -.001)
     1 then
                itype = 2 ! Sc
            elseif (d_thetae_dz .lt. -.001 .and. d_thetae_dz .ge. -.005)
     1 then
                itype = 3 ! Cu
            else ! d_thetae_dz .lt. -.005
                if(depth_m .gt. 5000)then
                    itype = 10 ! Cb
                else
                    itype = 3  ! Cu
                endif
            endif

        Elseif (temp_c.lt.-10. .and. temp_c.ge.-20.) then

            if (d_thetae_dz .lt. 0.) then
                if(depth_m .gt. 5000)then
                    itype = 10 ! Cb
                else
                    itype = 5 ! Ac
                endif
            else
                itype = 6 ! As
            endif

        Else !if ( temp_c.lt.-20.) then

            if (d_thetae_dz .ge. +.0005) then
                itype = 7 ! Cs
            elseif (d_thetae_dz .lt. +.0005 .and. 
     1              d_thetae_dz .ge. -.0005) then
                itype = 8 ! Ci
            else ! d_thetae_dz .lt. -.0005
                itype = 9 ! Cc
            endif

            if(depth_m .gt. 5000 .and. d_thetae_dz .lt. -.0000)then
                itype = 10 ! Cb
            endif

        Endif

        c2_type = c2_cloudtypes(itype)

!       write(6,*)' ',c2_type

        return
        end

        subroutine get_stability(nk,temp_1d,heights_1d,pressures_mb
     1                           ,kbottom,ktop,d_thetae_dz_1d)

cdoc    This routine returns stability at a given level given 1D array inputs
cdoc    Saturation is assumed and d(theta(e))/dz is returned

        real temp_1d(nk)                  ! Input
        real heights_1d(nk)               ! Input
        real pressures_mb(nk)             ! Input
        real d_thetae_dz_1d(nk)           ! Output
        real thetae_1d(nk)                ! Local

!       Calculate Stability
        klow  = max(kbottom-1,1)
        khigh = min(ktop+1,nk)

        do k = klow,khigh
            thetae_1d(k)  = OS_fast(temp_1d(k), pressures_mb(k))
        enddo ! k

        do k = kbottom,ktop
            km1  = max(k-1,1)
            kp1  = min(k+1,nk)

            d_thetae_dz_1d(k) = (thetae_1d(kp1) - thetae_1d(km1))
     1                / (heights_1d(kp1) - heights_1d(km1))

        enddo ! k

        return
        end

c
      subroutine isi3(slw,temp,ictype,iptype,xindex)
c     this version 2.0 dated 30 Dec 1991
c       m. politovich, rap, ncar
c       "dammit jim, I'm a scientist, not a programmer!"
cdoc  routine to calculate icing severity index
c
c     input:  slw     liquid water content in g/m3
c             temp    temperature in deg C
c             ictype  cloud type
c             iptype  precip type
c     output: index   severity index - real number from 0-6
c
c     note: there are really 3 categories, with
c            continuous or intermittent designation
c           1-3 = cats 1, 2 and 3 continuous
c           4-6 = cats 1, 2 and 3 intermittent
c           continuous:   cloud types st,as,cs
c           intermittent: cloud types cu,cs,cb,ac,sc
c
c     precip type: if zr, category 3 is designated
c
c     version 2.0  17 December 1991
c
c
c     scat and tcat are arrays of thresholds for slw and temp
      parameter (islw=3, itemp=3, ityp=2)
      dimension scat(islw),tcat(itemp)
c     iaray is the 3D severity index matrix
      dimension iaray(islw+1, itemp+1, ityp)
c
c     These are the lw and temp categories
      data scat /0.01,0.1,0.5/
      data tcat /-10.,-5.,0./
c     iaray is the array of severity index values
c       iaray (slw, temp, type)
c     this version 2.0 dated 30 Dec 1991
c
      data iaray/  0, 1, 2, 3, 0, 2, 2, 3,
     j             0, 2, 3, 3, 0, 0, 0, 0,
     b             0, 4, 5, 6, 0, 5, 5, 6,
     k             0, 5, 6, 6, 0, 0, 0, 0/
c
c
c     sort input values
c    **note that - test is for "ge" rather than "gt"
c
      i = 1
      do 460 ii = 1, islw
        if (slw.gt.scat(ii)) i = ii + 1
460     continue
761   j = 1
      do 461 ii = 1, itemp
        if (temp.gt.tcat(ii)) j = ii+1
461     continue
c     cloud type section (continuous/intermittent)
c      types: 0/none 1/st 2/sc 3/cu 4/ns 5/ac
c             6/as   7/cs 8/ci 9/cc 10/cb
c      continuous, k=1
c      intermittent k=2
c     ** note that - if there is no cloud type (=0) but there
c      IS slw calculated, the slw overrides**
      k = 1
      if (ictype.eq.2) k = 2
      if (ictype.eq.3) k = 2
      if (ictype.eq.5) k = 2
      if (ictype.ge.9) k = 2
c
c     precip type section (freezing rain/drizzle=cat 3)
c      types: 0/none 1/rn 2/sn 3/zr 4/sl 5/ha 6-10/no assignment
c     ** note that - freezing rain sets slw to category 4,
c      and, sets temperature to category 3 (even if it is
c      somehow diagnosed to "above freezing"
c      also - if no cloud type given (eg, below cloud), zr
c       is considered "continuous"***
      if (iptype.eq.3) j = 3
      if (iptype.eq.3) i = 4
c
!     write (6,5995) slw,temp,type
5995  format (' l,t,ty:',2f7.2, a2)
!     write (6,5996) i,j,k
5996  format (' i,j,k:',3i3)
c
c     assign severity index from sorted values
c
      index = iaray(i,j,k)
      xindex = index
c
c     ahhh, there's an easy solution to everything...almost
c
      return
      end

      function wb_melting_threshold(t_c,dbz,twet_snow)

cdoc  This function calculates the wet-bulb threshold for melting snow into
cdoc  rain as a function of dbz and t_c. In the absence of radar a default
cdoc  value is used.

      wb_melting_threshold = twet_snow  ! Units are C

      return
      end


        subroutine nowrad_virga_correction(r_pcp_type_2d,
     1                                     r_pcp_type_thresh_2d,
     1                                     t_sfc_k,
     1                                     td_sfc_k,
     1                                     istat_radar_3dref)

cdoc    Correct precip type for snow drying in sub-cloud layer based on sfc rh

        r_pcp_type_thresh_2d = r_pcp_type_2d

        if(r_pcp_type_2d .eq. 2.0)then

!           Threshold by surface dewpoint depression
            if(t_sfc_k - td_sfc_k .gt. 10.0)then
                if(istat_radar_3dref .eq. 0)then ! orig 2d data
                    r_pcp_type_thresh_2d = 0.  ! No Snow
                endif
            endif

        endif


        return
        end

        subroutine integrate_tpw(sh,sh_sfc,p,psfc,ni,nj,nk,tpw)

!       Calculate tpw using pressures of the levels

        use constants_laps, ONLY: grav ! m/s**2

        implicit none

        real sh(ni,nj,nk)     ! (dimensionless)          Input
        real sh_sfc           ! (dimensionless)          Input
        real p(ni,nj,nk)      ! (Pa)                     Input
        real psfc             ! (Pa)                     Input
        real tpw(ni,nj)       ! (meters)                 Output

        integer ni,nj,nk,i,j,k
        real pbot,ptop,shbot,shtop,shave,dp,rho_water

        rho_water = 1e3       ! mass in kilograms of a cubic meter of water

        do i = 1,ni
        do j = 1,nj
            tpw(i,j) = 0.

            do k = 2,nk
                if(p(i,j,k-1) .lt. psfc)then    ! layer is all above the sfc 
                    pbot = p(i,j,k-1)
                    ptop = p(i,j,k)
                    shbot = sh(i,j,k-1)
                    shtop = sh(i,j,k)
                elseif(p(i,j,k) .lt. psfc) then ! layer straddles the sfc
                    pbot = psfc          
                    ptop = p(i,j,k)
                    shbot = sh_sfc          
                    shtop = sh(i,j,k)
                else                            ! layer is all below the sfc
                    pbot = 0.
                    ptop = 0.
                endif

                if(pbot .gt. 0.)then
                    shave = 0.5 * (shbot + shtop)
                    dp = pbot - ptop
                    tpw(i,j) = tpw(i,j) + shave * dp
                endif

            enddo ! k

            tpw(i,j) = (tpw(i,j) / grav) / rho_water 

!           if(i .eq. 1 .and. j .eq. 1)then
!               write(6,*)' integrate_tpw: ',tpw(i,j)
!           endif

        enddo ! j
        enddo ! i
 
        return
        end
