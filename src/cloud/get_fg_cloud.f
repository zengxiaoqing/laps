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

        subroutine get_modelfg(cf_modelfg,t_modelfg                      ! O
     1                        ,model_q_3d                                ! O
     1                        ,model_ref_3d                              ! O
     1                        ,default_clear_cover,r_missing_data        ! I
     1                        ,temp_3d,heights_3d,cld_hts                ! I
     1                        ,i4time_needed,ilaps_cycle_time            ! I
     1                        ,ni,nj,klaps,KCLOUD                        ! I
     1                        ,istatus)                                  ! O

!       Obtain model first guess cloud cover and temperature fields
!       This should probably be free of very small scale horizontal structures
!       for best results when combining with the satellite data

!                Steve Albers   Use model first guess info to guess cloud cover
!       1994     Steve Albers   Read in SH from model background (instead of
!                                                                         RH)
!       1995 Dec Steve Albers   QC check to prevent cf_modelfg > 1.0
!       1999     Steve Albers   Added in LWC/ICE model first guess
!                               Simple usage for the moment.

        real heights_3d(ni,nj,klaps)       ! Input
        real temp_3d(ni,nj,klaps)          ! Input
        real cf_modelfg(ni,nj,KCLOUD)      ! Output
        real t_modelfg(ni,nj,KCLOUD)       ! Output
        real cld_hts(KCLOUD)               ! Input
        real model_q_3d(ni,nj,klaps)       ! Output
        real model_ref_3d(ni,nj,klaps)     ! Output

        real model_lwc_3d(ni,nj,klaps)     ! Local (units are mixing ratio)
        real model_ice_3d(ni,nj,klaps)     ! Local (units are mixing ratio)

        real make_rh,lwc_modelfg,ice_modelfg

        real icethresh, liqthresh

        character*31 ext
        character*3 var_2d

        write(6,*)
        write(6,*)' Getting first guess cloud cover'

!       mode = 1 ! Prefer cloud derived from model RH 
        mode = 2 ! use maximum of first guess LWC and RH implied cloud

!       Initialize model first guess cover field with default value
        do k = 1,KCLOUD
        do j = 1,nj
        do i = 1,ni
            cf_modelfg(i,j,k) = default_clear_cover
            t_modelfg(i,j,k) = 0.
        enddo
        enddo
        enddo

        i_hum_high = 0
        i_hum_low = 0
        i_hum_ok = 0

        i_condensate = 0

        istat_sh  = 0
        istat_lwc = 0
        istat_ice = 0

        write(6,*)' Getting MODEL LWC background'

!       Get Model First Guess LWC
        var_2d = 'LWC'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_lwc_3d,istat_lwc)

        if(istat_lwc .ne. 1)then
            write(6,*)' No first guess available for ',var_2d
        endif

        write(6,*)' Getting MODEL ICE background'

!       Get Model First Guess ICE
        var_2d = 'ICE'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_ice_3d,istat_ice)

        if(istat_ice .ne. 1)then
            write(6,*)' No first guess available for ',var_2d
        endif

        write(6,*)' Getting MODEL SH background'

!       Get Model First Guess SH
        var_2d = 'SH'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_q_3d,istat_sh)

        if(istat_sh .ne. 1)then
            write(6,*)' No first guess available for ',var_2d
            istatus = 0
            return
        endif

        write(6,*)' Getting MODEL REF background'

!       Get Model First Guess REF
        var_2d = 'REF'
        call get_modelfg_3d(i4time_needed,var_2d,ni,nj,klaps
     1                     ,model_ref_3d,istat_ref)

        if(istat_ref .ne. 1)then
            write(6,*)' No first guess available for ',var_2d
            model_ref_3d = r_missing_data
        endif

!       Check status of LWC/ICE from model first guess
        if(istat_lwc .ne. 1 .or. istat_ice .ne. 1)then
            continue

        else ! We are using model LWC/ICE fields
            call get_grid_spacing_cen(grid_spacing_cen_m,istatus)
            if(istatus .ne. 1)return

            IF(grid_spacing_cen_m .LE. 10000.)THEN  ! Grid spacing in meters
                icethresh  = 0.000005 ! These are in kg kg{-1} (mixing ratio)
                snowthresh = 0.000003
                liqthresh  = 0.000003
            ELSE
                icethresh  = 0.000005
                snowthresh = 0.000025
                liqthresh  = 0.000025
            ENDIF

        endif ! Status for MODEL data

!       Remap to cloud height grid and convert to cloud cover
        t_ref = -10. ! colder than this ice saturation is assumed

        do k = 1,KCLOUD
        do j = 1,nj
          istat_z = 0
          do i = 1,ni

!           Find the model pressure at this location in the cloud height grid
            if(i-1 .eq. (i-1)/10*10)then ! Update every 10th grid point
                z_laps = height_to_zcoord2(cld_hts(k),heights_3d
     1                  ,ni,nj,klaps,i,j,istat_z)

                if(istat_z .ne. 1)then
!                   We're OK if cloud height grid is above pressure grid
                    if(cld_hts(k) .le. heights_3d(i,j,klaps))then
                        write(6,*)' Error: Bad status from '
     1                           ,'height_to_zcoord2'
                        istatus = 0
                        return
                    endif

                else
                    z_laps = max(1.,min(z_laps,float(klaps)-.001))
                    iz_laps = int(z_laps)
                    frac = z_laps - iz_laps

                    p_modelfg = pressure_of_level(iz_laps) * (1. - frac)       
     1                        + pressure_of_level(iz_laps+1)  * frac

                    p_modelfg_mb = p_modelfg * .01

                endif ! istat_z .ne. 1

            endif

            if(istat_z .ne. 1)then
                i_grid_high = i_grid_high + 1
                cf_modelfg(i,j,k) = default_clear_cover
                go to 1000
            endif

            if(istat_lwc .eq. 1 .and. istat_lwc .eq. 1)then
                lwc_modelfg =  model_lwc_3d(i,j,iz_laps)   * (1. - frac)
     1                      +  model_lwc_3d(i,j,iz_laps+1)       * frac

                ice_modelfg =  model_ice_3d(i,j,iz_laps)   * (1. - frac)
     1                      +  model_ice_3d(i,j,iz_laps+1)       * frac

                if(lwc_modelfg .ge. liqthresh .or. 
     1             ice_modelfg .ge. icethresh        )then
                    cf_modelfg(i,j,k) = 1.
                    i_condensate = i_condensate + 1
                else
                    cf_modelfg(i,j,k) = default_clear_cover
                endif

            endif

            if(istat_sh .eq. 1)then
              if(.NOT. (istat_lwc .eq. 1 .and. istat_lwc .eq. 1) 
     1                          .OR.   mode .eq. 2
     1                                                           )then
!               Find the model temp at this location in the cloud height grid
                t_modelfg(i,j,k) =  temp_3d(i,j,iz_laps)   * (1. - frac)      
     1                           +  temp_3d(i,j,iz_laps+1) * frac

                t_modelfg_c = t_modelfg(i,j,k) - 273.15

!               Find the model sh at this location in the cloud height grid
                q_modelfg =  model_q_3d(i,j,iz_laps)    * (1. - frac)       
     1                    +  model_q_3d(i,j,iz_laps+1)  * frac

                q_modelfg_gkg = q_modelfg * 1000.

                rh_modelfg = make_rh(p_modelfg_mb              ! fractional rh
     1                     ,t_modelfg_c,q_modelfg_gkg,t_ref)

!               QC the rh
                rh_qc = rh_modelfg                             ! fractional rh

                if(rh_qc .gt. 1.0)then
                    rh_qc = 1.0
                    i_hum_high = i_hum_high + 1
                elseif(rh_qc .lt. 0.0)then
                    rh_qc = 0.0
                    i_hum_low = i_hum_low + 1
                else
                    i_hum_ok = i_hum_ok + 1
                endif

                if(cld_hts(k) .gt. 11000.)rh_qc = .01   ! set upper lvls to dry
                                                        ! counters model (ruc) 
                                                        ! moist bias

                if(mode .eq. 1)then
                    cf_modelfg(i,j,k) = rh_to_cldcv(rh_qc)       ! fractional_rh
                else
                    cf_modelfg(i,j,k) = 
     1              max(cf_modelfg(i,j,k),rh_to_cldcv(rh_qc))    ! fractional_rh
                endif

              endif

            endif

 1000     enddo ! i
        enddo ! j
        enddo ! k (cloud height array level)

        write(6,*)' # RH values QCed  high/low/ok'
     1            ,i_hum_high,i_hum_low,i_hum_ok
        write(6,*)' # cloud height grids above pressure grid'
     1            ,i_grid_high
        write(6,*)' # points set to cloud based on condensate = '
     1           ,i_condensate

        istatus = 1

        return
        end
