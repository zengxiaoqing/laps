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

        subroutine get_precip_accum(i4time_beg,i4time_end   ! Input
     1          ,imax,jmax,kmax                   ! Input
     1          ,lat,lon,topo      ! Input
     1          ,ilaps_cycle_time  ! Input
     1          ,radarext_3d_accum ! Input
     1          ,snow_accum,precip_accum,frac_sum      ! Outputs
     1          ,istatus)                 ! Output

!       Steve Albers 1991
!       Steve Albers 1995 Dec  Modify radar call to read_radar_ref
!       Returns Accumulated Snow and Liquid Precip
!       This routine uses 3D precip type (using temp data from the
!       closest cycle) as a cutoff for snow/no snow. This is calculated and
!       applied every time we are leaving the valid time window for a
!       particular set of environmental data. Typically, the first half of the
!       cycle time comes from the previous cycle's temp data, the second half of
!       the period uses current data. A mask is also used so that the 3D
!       precip type calculations are performed only where there are accumulating
!       echoes. In other words, this code is more complex that it would
!       otherwise be so that real-time speed is optimized.

        integer*4 MAX_FILES
        parameter (MAX_FILES = 3000)

!       Input
        real*4 lat(imax,jmax)
        real*4 lon(imax,jmax)
        real*4 topo(imax,jmax)

!       Output
        real*4 snow_accum(imax,jmax) ! M
        real*4 precip_accum(imax,jmax) ! M

        real*4 snow_accum_pd(imax,jmax)
        real*4 snow_rate(imax,jmax) ! M/S
        real*4 precip_rate(imax,jmax) ! M/S
        real*4 dbz_2d(imax,jmax)
        real*4 t_sfc_k(imax,jmax)
        real*4 td_sfc_k(imax,jmax)
        real*4 pres_sfc_pa(imax,jmax)
        real*4 tw_sfc_k(imax,jmax)
        real*4 temp_3d(imax,jmax,kmax)
        real*4 height_3d(imax,jmax,kmax)
        real*4 temp_col_max(imax,jmax)
        real*4 rh_3d(imax,jmax,kmax)
        real*4 pressures_mb(kmax)
        logical l_mask(imax,jmax)
        integer*2 i2_pcp_type_2d(imax,jmax)
        integer*2 i2_cldpcp_type_3d(imax,jmax,kmax)
        integer ipcp_1d(kmax)

        real*4 grid_ra_ref(imax,jmax,kmax)
        real*4 grid_ra_vel(imax,jmax,kmax)
        real*4 grid_ra_nyq(imax,jmax,kmax)

        character*9 asc_tim_9,asc_tim_9_beg,asc_tim_9_end
        integer*4 i4time_file(MAX_FILES)
        real*4 frac(MAX_FILES)
        character c_fnames(MAX_FILES)*80

        character*255 c_filespec

        character*4  radar_name ! Local

        character*3 var_2d
        character*31  ext, radarext_3d_accum
        character*10  units_2d
        character*125 comment_2d

        integer*4 iarg

        logical l_first_sfc_update_completed

        character*80 grid_fnam_common
        common / grid_fnam_cmn / grid_fnam_common

        data mode_radar/1/

        radarext_3d_accum = 'vrc'

!       istatus = .true.

        max_radar_gap = float(ilaps_cycle_time) + 1200. ! * 1.33334

        im = imax/2 ! 15
        jm = jmax/2 ! 56

        do k = 1,kmax
            pressures_mb(k) = pressure_of_level(k) / 100.
        enddo ! k

        call make_fnam_lp(i4time_beg,asc_tim_9_beg,istatus)
        call make_fnam_lp(i4time_end,asc_tim_9_end,istatus)

        write(6,*)' Radar accumulation from ',asc_tim_9_beg,
     1                                 ' to ',asc_tim_9_end

!       Get File Times
        call get_filespec(radarext_3d_accum(1:3),2,c_filespec,istatus)

        call    Get_file_names(  c_filespec,
     1                   i_nbr_files_ret,
     1                   c_fnames,
     1                   max_files,
     1                   i_status )

        min_diff = 1999999999

        if(i_nbr_files_ret .gt. 0)then
            call get_directory_length(c_fnames(1),lenf)
        else ! Error Condition
            write(6,*)' No Radar Data Available for Snow/Precip Accumula
     1tion'
            istatus = 0
            return
        endif

10      do i=1,i_nbr_files_ret
            asc_tim_9 = c_fnames(i)(lenf+1:lenf+9)
            call i4time_fname_lp(asc_tim_9,I4time_file(i),istatus)
            if(istatus .ne. 1)then
                write(6,*)' Bad return from i4time_fname_lp,',
     1            ' called from get_precip_accum: ',asc_tim_9
                return
            endif
        enddo

        call get_fracs(i4time_beg,i4time_end,max_radar_gap,i_nbr_files_r
     1et
     1            ,i4time_file,frac,frac_sum,istatus)

        if(istatus .ne. 1)then
            return
        endif

!       Initialize both period and total precip/snow

!       Snow is accumulated in "periods" corresponding to when the ancillary
!       LAPS data is read in. Each time new ancillary data is read in, the
!       "period" snow accumulation is added to the overall snow total for the
!       window.

!       This is done so that precip type (which determines snow/no snow) only
!       has to be calculated only once for every time ancillary data is read in.
!       The precip type is calculated over a masked area corresponding to where
!       precip (potentially snow) has occurred over the "period".

!       This approximation may lead to slightly inaccurate accumulations
!       over the period if the precip type is changing over the period.
!       Note that precip type is a function of temperature, rh, and even
!       reflectivity - all of which can change in various ways over the
!       laps cycle. The approximation is there to increase the speed of
!       execution. When computers get better, it may be desirable to
!       redesign this so the precip type (and snow/rain ratio) gets calculated
!       for each radar scan. Of course this can also be taken care of by
!       decreasing the overall LAPS analysis cycle time.

        do j = 1,jmax
        do i = 1,imax
            snow_accum(i,j) = 0.
            snow_accum_pd(i,j) = 0.
            precip_accum(i,j) = 0.
            l_mask(i,j) = .false.
        enddo ! i
        enddo ! j

        i4_interval = (i4time_end - i4time_beg)
        i4time_mid = i4time_beg + i4_interval / 2
        i4_tol = max(float(ilaps_cycle_time) * 0.6,20.)
        i4time_temp = 0

        l_first_sfc_update_completed = .false.

!       Loop through all the radar scan times
        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then ! This scan is needed for the window at hand

            i4time_radar = i4time_file(ifile)

            call make_fnam_lp(i4time_radar,asc_tim_9,istatus)

!           write(6,101)asc_tim_9,frac(ifile)
!101        format(' Time, Frac = ',a9,2x,f6.3)

            write(6,*)

!           Determine whether we need to update the sfc data to match the radar
!           and if we should add the period snow accum onto the total
            if(abs(i4time_radar - i4time_temp) .gt. ilaps_cycle_time/2       
     1                                                             )then

!               Determine whether to add to snow accum total based on precip 
!               type using ancillary LAPS data valid for the interval just 
!               ended. Note that the first time around, we are just 
!               initializing and don't need to add to the snow accum, we just 
!               read in the initial ancillary LAPS data.

                if(l_first_sfc_update_completed .eqv. .true.)then


                    do j = 1,jmax
                    do i = 1,imax
                        if(snow_accum_pd(i,j) .gt. 1e-10)l_mask(i,j) = .
     1true.
                    enddo
                    enddo


                    write(6,*)' Compute 3D precip type over masked area'

!                   Note that the reflectivities here are valid at the end
!                   of the accumulating period. This means that the wb threshold!                   is calculated using reflectivities that may not be
!                   representative for the entire accumulation subperiod
!                   (typically one half of the laps cycle time).

                    call cpt_pcp_type_3d(temp_3d,rh_3d,pressures_mb
     1                  ,grid_ra_ref,l_mask
     1                  ,imax,jmax,kmax,i2_cldpcp_type_3d,istatus)
                    if(istatus .ne. 1)then
                        return
                    endif

                    I4_elapsed = ishow_timer()

                    write(6,*)' Compute sfc precip type'
                    call get_sfc_preciptype(pres_sfc_pa,t_sfc_k,td_sfc_k
     1                                  ,i2_cldpcp_type_3d,dbz_2d
     1                                  ,i2_pcp_type_2d,imax,jmax,kmax)

                    I4_elapsed = ishow_timer()

                    write(6,*)' Adding in accumulation for intervening p
     1eriod'

                    n_pcp_pts = 0
                    n_snw_pts = 0
                    n_nopcp_pts = 0
                    n_zr_pts = 0

                    do j = 1,jmax
                    do i = 1,imax

                        if(l_mask(i,j))then
                            n_pcp_pts = n_pcp_pts + 1
                            iarg = i2_pcp_type_2d(i,j)/16

                            if(iarg .eq. 2)then     ! precip type is snow

                                r_pcp_type        = iarg
                                r_pcp_type_thresh = iarg

                                if(.true.)then  ! testcode
                                    call nowrad_virga_correction(
     1                              r_pcp_type,
     1                              r_pcp_type_thresh,
     1                              t_sfc_k(i,j),
     1                              td_sfc_k(i,j),
     1                              istatus_3dref)
                                endif

                                if(r_pcp_type_thresh .eq. 2.0)then
                                    snow_accum(i,j)
     1                          = snow_accum(i,j) + snow_accum_pd(i,j)
                                    n_snw_pts = n_snw_pts + 1
                                else
                                    write(6,*)' Nowrad_virga thresh out 
     1#1',i,j
                                endif

                            elseif(iarg .eq. 0)then ! precip type is no precip
                                n_nopcp_pts = n_nopcp_pts + 1
                                if(n_nopcp_pts .le. 20)
     1                          write(6,*)' No2dPcpType pt at',i,j
                            elseif(iarg .eq. 3)then ! precip type is freezing rain
                                n_zr_pts = n_zr_pts + 1
                            endif

                        endif
                    enddo ! i
                    enddo ! j

!                   This is useful for debugging
                    do k = 1,kmax
                        iarg = i2_cldpcp_type_3d(im,jm,k) / 16
                        ipcp_1d(k) = iarg
                    enddo ! k
                    write(6,101)l_mask(im,jm),snow_accum_pd(im,jm)
     1                          ,snow_accum(im,jm)
     1                          ,t_sfc_k(im,jm),td_sfc_k(im,jm)
     1                          ,i2_pcp_type_2d(im,jm)
     1                          ,(ipcp_1d(k),k=1,min(kmax,10))
101                 format(1x,'accum/tw/type2,3d',l2,2e11.3,2f6.1,i3,1x,
     110i2)

                    write(6,*)' # of Points Snow/Precip/ZR = '
     1                          ,n_snw_pts,n_pcp_pts,n_zr_pts

                    if(n_nopcp_pts .gt. 0)
     1          write(6,*)' WARNING: n_nopcp_pts = ',n_nopcp_pts

                endif ! l_first_sfc_update_completed

!               Initialize
                do j = 1,jmax
                do i = 1,imax
                    snow_accum_pd(i,j) = 0.
                    l_mask(i,j) = .false.
                enddo ! i
                enddo ! j

!               Read in surface data that is time matched to the radar data
                write(6,*)
                write(6,*)' Updating Surface and 3D Temp Information'

!               Read in surface temp data
                var_2d = 'T'
                ext = 'lsx'
                call get_laps_2dgrid(i4time_radar,i4_tol,i4time_temp
     1   ,ext,var_2d,units_2d,comment_2d,imax,jmax,t_sfc_k,0,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS Sfc Temp not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

!               Read in surface dewpoint data
                var_2d = 'TD'
                ext = 'lsx'
                call get_laps_2d(i4time_temp,ext,var_2d
     1                ,units_2d,comment_2d,imax,jmax,td_sfc_k,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS Sfc Dewpoint not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

!               Read in surface pressure data
                var_2d = 'PS'
                ext = 'lsx'
                call get_laps_2d(i4time_temp,ext,var_2d
     1             ,units_2d,comment_2d,imax,jmax,pres_sfc_pa,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS Sfc Pressure not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

                I4_elapsed = ishow_timer()

                write(6,*)' Getting LT1 Height/Temp' 
                var_2d = 'HT' 
                ext = 'lt1'

                call get_laps_3d(i4time_temp
     1                  ,imax,jmax,kmax,ext,var_2d
     1                  ,units_2d,comment_2d,height_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS 3D Height not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

                var_2d = 'T3' 
                ext = 'lt1'

                call get_laps_3d(i4time_temp
     1                  ,imax,jmax,kmax,ext,var_2d
     1                  ,units_2d,comment_2d,temp_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS 3D Temp not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

!               Calculate column max temperatures
                do j = 1,jmax
                do i = 1,imax
                    temp_col_max(i,j) = t_sfc_k(i,j)
                    k_sfc = int(zcoord_of_pressure(pres_sfc_pa(i,j)))
                    do k = k_sfc+1,kmax
                        temp_col_max(i,j) =
     1                  max(temp_col_max(i,j),temp_3d(i,j,k))
                    enddo ! k
                enddo ! i
                enddo ! j


                write(6,*)' Getting LH3 file (or equivalent)'
                var_2d = 'RHL'
                ext = 'lh3'
                call get_laps_3dgrid(i4time_temp,ilaps_cycle_time ! *2
     1                                  ,i4time_rh
     1          ,imax,jmax,kmax,ext,var_2d
     1                  ,units_2d,comment_2d,rh_3d,istatus)
                if(istatus .ne. 1)then
                    write(6,*)' LAPS 3D RH not available'
                    frac_sum = -1.0 ! Turns off the wait loop for more radar
                    return
                endif

                l_first_sfc_update_completed = .true.

                write(6,*)' Ancillary LAPS data update completed'
                write(6,*)

                I4_elapsed = ishow_timer()

            else ! l_first_sfc_update_completed = .true.
                write(6,*)' No Ancillary LAPS data update needed'

            endif ! We need to update the LAPS information

!           Get LAPS reflectivities at the surface (or immediately above it)
            write(6,*)

!           Repeat read in radar data with low level reflectivities filled in
            call read_radar_3dref(i4time_radar,
!    1       0,i4_dum,
     1       .true.,imax,jmax,kmax,
     1       radarext_3d_accum,lat,lon,topo,
     1       .true.,.false.,
     1       height_3d,
     1       grid_ra_ref,
     1       rlat_radar,rlon_radar,rheight_radar,radar_name,
     1       n_ref,istatus_2dref,istatus_3dref)

            if(istatus_2dref .eq. 0)then
                write(6,*)' Error in reading radar data in PRECIP ACCUM'
                frac_sum = -1.0 ! Turns off the wait loop for more radar
                istatus = 0
                return
            endif

            write(6,*)' Call get_low_ref'

            call get_low_ref(grid_ra_ref,pres_sfc_pa,imax,jmax,kmax,dbz_
     12d)

            write(6,*)' Incrementing Precip Accumulation '
     1               ,'rate for this scan (call zr)'

            call zr(dbz_2d,imax,jmax,precip_rate)

            do j = 1,jmax
            do i = 1,imax
                precip_accum(i,j) = precip_accum(i,j)
     1                            + precip_rate(i,j) * frac(ifile)
            enddo ! i
            enddo ! j


            write(6,*)' Incrementing Snow Accumulation_pd for this scan'

            call zs(precip_rate,temp_col_max,imax,jmax,snow_rate)

            do j = 1,jmax
            do i = 1,imax
                snow_accum_pd(i,j) = snow_accum_pd(i,j)
     1                             + snow_rate(i,j) * frac(ifile)
            enddo ! i
            enddo ! j

            write(6,202)dbz_2d(im,jm),snow_rate(im,jm),snow_accum_pd(im,
     1jm)
202         format(1x,'dbz/rate/accum_pd',3e12.4)

          endif ! Frac > 0 (This radar file is within the accumulating window

!         write(6,*)' Cycle to next file'

        enddo ! ifile (Loop through all radar file times)


!       Add in snow accumulation from final period; first define the mask
        do j = 1,jmax
        do i = 1,imax
            if(snow_accum_pd(i,j) .gt. 1e-10)l_mask(i,j) = .true.
        enddo
        enddo


        write(6,*)' Compute 3D precip type over masked area'

        call cpt_pcp_type_3d(temp_3d,rh_3d,pressures_mb,grid_ra_ref,l_ma
     1sk
     1          ,imax,jmax,kmax,i2_cldpcp_type_3d,istatus)
        if(istatus .ne. 1)then
            return
        endif

        write(6,*)' Compute sfc precip type'
        call get_sfc_preciptype(pres_sfc_pa,t_sfc_k,td_sfc_k
     1        ,i2_cldpcp_type_3d,dbz_2d,i2_pcp_type_2d,imax,jmax,kmax)

        write(6,*)' Adding in accumulation for the last period'

        n_pcp_pts = 0
        n_snw_pts = 0
        n_nopcp_pts = 0
        n_zr_pts = 0

        do j = 1,jmax
        do i = 1,imax

            if(l_mask(i,j))then
                n_pcp_pts = n_pcp_pts + 1
                iarg = i2_pcp_type_2d(i,j) / 16

                if(iarg .eq. 2)then     ! precip type is snow

                    r_pcp_type        = iarg
                    r_pcp_type_thresh = iarg

                    if(.true.)then ! testcode
                        call nowrad_virga_correction(
     1                              r_pcp_type,
     1                              r_pcp_type_thresh,
     1                              t_sfc_k(i,j),
     1                              td_sfc_k(i,j),
     1                              istatus_3dref)
                    else
                        write(6,*)' Nowrad_virga thresh out #2',i,j
                    endif

                    if(r_pcp_type_thresh .eq. 2.0)then
                        snow_accum(i,j) = snow_accum(i,j) + snow_accum_p
     1d(i,j)
                        n_snw_pts = n_snw_pts + 1
                    endif

                elseif(iarg .eq. 0)then ! precip type is no precip
                    n_nopcp_pts = n_nopcp_pts + 1
                    if(n_nopcp_pts .le. 20)write(6,*)' No2dPcpType pt at
     1',i,j
                elseif(iarg .eq. 3)then ! precip type is freezing rain
                    n_zr_pts = n_zr_pts + 1
                endif

            endif

        enddo ! i
        enddo ! j


!       This is useful for debugging
        do k = 1,kmax
            iarg = i2_cldpcp_type_3d(im,jm,k) / 16
            ipcp_1d(k) = iarg
        enddo ! k
        write(6,101)l_mask(im,jm),snow_accum_pd(im,jm)
     1                          ,snow_accum(im,jm)
     1                          ,t_sfc_k(im,jm),td_sfc_k(im,jm)
     1                          ,i2_pcp_type_2d(im,jm)
     1                          ,(ipcp_1d(k),k=1,min(kmax,10))

        write(6,*)' # of Points Snow/Precip/ZR = ',n_snw_pts,n_pcp_pts,n
     1_zr_pts

        if(n_nopcp_pts .gt. 0)write(6,*)' WARNING: n_nopcp_pts = ',n_nop
     1cp_pts

        write(6,*)' Converting from average rate to actual accumulation'

!       Convert from time averaged rate to accumulation
        do j = 1,jmax
        do i = 1,imax
            snow_accum(i,j) = snow_accum(i,j) * i4_interval
            precip_accum(i,j) = precip_accum(i,j) * i4_interval
        enddo ! i
        enddo ! j

        write(6,*)' Final snow_accum(im,jm) = ',snow_accum(im,jm)
        write(6,*)' Final precip_accum(im,jm) = ',precip_accum(im,jm)

        return
        end

        subroutine get_fracs(i4time_beg,i4time_end,max_radar_gap
     1          ,i_nbr_files_ret
     1          ,i4time_file,frac,frac_sum,istatus)

!       Steve Albers    1991    This routine calculates linear combination
!                               coefficients for radar scans which can be used
!                               to arrive at an integrated precipitation rate
!                               over a specified time window

        integer*4 MAX_FILES
        parameter (MAX_FILES = 3000)

        character*9 asc_tim_9
        real*4 frac(MAX_FILES)
        integer*4 i4time_file(MAX_FILES)

        i4_interval = i4time_end - i4time_beg

        do i = 1,i_nbr_files_ret
!           write(6,301)i4time_beg,i4time_file(i),i4time_end
!301        format(1x,3i12)
            frac(i) = 0.
        enddo

        do i = 1,i_nbr_files_ret-1
            ibeg = i
            iend = i+1

            interval_between_files = i4time_file(iend) - i4time_file(ibe
     1g)

            frac_between_files = float(interval_between_files)
     1                 /     float(i4_interval)

            if(i4time_file(ibeg) .ge. i4time_beg .and.
     1       i4time_file(iend) .le. i4time_end        )then ! Fully Within Pd
                frac(ibeg) = frac(ibeg) + 0.5 * frac_between_files
                frac(iend) = frac(iend) + 0.5 * frac_between_files

                if(interval_between_files .gt. max_radar_gap)then
                    write(6,*)' ERROR: Gap in radar files (min) >'
     1                                          ,max_radar_gap/60
                    istatus = 0
                    return
                endif

            endif

            if(i4time_file(ibeg) .lt. i4time_beg .and.
     1       i4time_file(iend) .gt. i4time_beg .and.
     1       i4time_file(iend) .le. i4time_end        )then ! Straddle Beginning
                a = i4time_beg - i4time_file(ibeg)
                b = i4time_file(iend) - i4time_beg
                partial_frac = b/(a+b)
                frac(ibeg) = frac(ibeg) + (      0.5 * partial_frac)
     1                          * frac_between_files * partial_frac
                frac(iend) = frac(iend) + (1.0 - 0.5 * partial_frac)
     1                          * frac_between_files * partial_frac

                if(interval_between_files .gt. max_radar_gap)then
                    write(6,*)' ERROR: Gap in radar files (min) >'
     1                                          ,max_radar_gap/60
                    istatus = 0
                    return
                endif

            endif

            if(i4time_file(ibeg) .lt. i4time_end   .and.
     1       i4time_file(ibeg) .ge. i4time_beg   .and.
     1       i4time_file(iend) .gt. i4time_end       )then ! Straddle End
                a = i4time_end - i4time_file(ibeg)
                b = i4time_file(iend) - i4time_end
                partial_frac = a/(a+b)
                frac(ibeg) = frac(ibeg) + (1.0 - 0.5 * partial_frac)
     1                  * frac_between_files * partial_frac
                frac(iend) = frac(iend) + (      0.5 * partial_frac)
     1                  * frac_between_files * partial_frac

                if(interval_between_files .gt. max_radar_gap)then
                    write(6,*)' ERROR: Gap in radar files (min) >'
     1                                          ,max_radar_gap/60
                    istatus = 0
                    return
                endif

            endif

            if(i4time_file(ibeg) .lt. i4time_beg   .and.
     1       i4time_file(iend) .gt. i4time_end    )then ! Brackets the Pd
                i4time_mid = i4time_beg + (i4time_end - i4time_beg) / 2
                frac_mid = float(i4time_mid        - i4time_file(ibeg))
     1           /       float(i4time_file(iend) - i4time_file(ibeg))
                frac(ibeg) = 1.0 - frac_mid
                frac(iend) = frac_mid

                if(interval_between_files .gt. max_radar_gap)then
                    write(6,*)' ERROR: Gap in radar files (min) >'
     1                                          ,max_radar_gap/60
                    istatus = 0
                    return
                endif

            endif

        enddo ! i

        frac_sum = 0.
        do ifile = 1,i_nbr_files_ret
          if(frac(ifile) .gt. 0.)then
            call make_fnam_lp(i4time_file(ifile),asc_tim_9,istat_fnam)
            write(6,101)asc_tim_9,frac(ifile)
 101        format('        FileTime, Frac = ',a9,2x,f8.5)
            frac_sum = frac_sum + frac(ifile)
          endif
        enddo

        if(abs(frac_sum - 1.0) .gt. 1e-5)then
!           Note: we can here potentially wait for more radar data in driver
            write(6,*)' ERROR: Fractions do not add up to 1.0',frac_sum
            istatus = 0
        endif

        if(i_nbr_files_ret .gt. 0)then
            if(i4time_file(1) .gt. i4time_beg)then
                write(6,*)' Radar files begin after start of accumulatio
     1n window'
                frac_sum = -1.0 ! Turns off the wait loop for more radar
            endif
        endif

999     if(istatus .ne. 1)then
            write(6,*)' Insufficient files within time window'
            return
        endif

        istatus = 1
        return
        end
