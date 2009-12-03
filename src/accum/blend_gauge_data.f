
        subroutine blend_gauge_data(    i4time,ni,nj,maxsta        ! I
     1                                 ,r_missing_data             ! I
     1                                 ,lat,lon                    ! I
     1                                 ,pcp_bkg_m                  ! I
     1                                 ,ilaps_cycle_time           ! I
     1                                 ,closest_radar              ! I
     1                                 ,precip_accum_m)            ! I/O


!       Blend the radar derived and background first guess to provide a 
!       spatially continuous precip estimate. This estimate is then corrected
!       using rain gauges. If the radar/background field has any missing values
!       a gauge only analysis is done instead.

        use mem_namelist, only: l_accum_fg, l_accum_radar, l_accum_gauge

        include 'read_sfc.inc'
        include 'constants.inc'

        real precip_accum_m(ni,nj)! input radar / output precip analysis (M)
        real pcp_2d_in(ni,nj)     ! gauge corrected analysis (IN)
        real pcp_bkg_m(ni,nj)     ! background field for gauge analysis (M)
        real pcp_cmb_m(ni,nj)     ! background + radar field combined (M)
        real pcp_cmb_in(ni,nj)    ! background + radar field combined (IN)
        real wt_bkg_a(ni,nj)      ! background weight for gauge analysis
        real one(ni,nj)           ! unity array for bias analysis
        real zero(ni,nj)          ! zero array for gauge only analysis
        real bias_anal(ni,nj)     ! bias analysis
        real closest_radar(ni,nj) ! M
        real lat(ni,nj)
        real lon(ni,nj)

        real pcp_gauge(maxsta)    ! relevant time of accumulation is selected
        character*5 c_field

        integer ilaps(maxsta),jlaps(maxsta)

        logical l_accum_bias_ratio ! Apply bias correction as a constant ratio
                                   ! (given we are analyzing radar/fg & gauges)

        logical l_gauge_only ! Gauge only analysis

        l_gauge_only = (.not. l_accum_fg) .AND. (.not. l_accum_radar) 
     1                                    .AND.        l_accum_gauge
        l_accum_bias_ratio = .true. 

!       Combine background and radar field given radar gap areas
        n_radar = 0
        n_bkg = 0
        n_msg_rdr_bkg = 0

        do i = 1,ni
        do j = 1,nj
            if(precip_accum_m(i,j) .ne. r_missing_data)then 
                n_radar = n_radar + 1 ! ! we have radar input
                pcp_cmb_m(i,j) = precip_accum_m(i,j) 
   
            elseif(pcp_bkg_m(i,j) .ne. r_missing_data)then
                n_bkg = n_bkg + 1
                pcp_cmb_m(i,j) = pcp_bkg_m(i,j)
            
            else
                pcp_cmb_m(i,j) = r_missing_data
                n_msg_rdr_bkg = n_msg_rdr_bkg + 1

            endif
        enddo ! j
        enddo ! i

        write(6,*)
        write(6,*)' Subroutine blend_gauge_data (1hr pcp inches)...'
 
        write(6,*)' Number of radar points = ',n_radar
        write(6,*)' Number of background points = ',n_bkg
        write(6,*)' Number of missing points = ',n_msg_rdr_bkg

        write(6,*)
        write(6,*)'   #  Name        Gauge Analyzed  Range'

        call read_sfc_precip(i4time,btime,n_obs_g,n_obs_b,
     &           stations,provider,lat_s,lon_s,elev_s,
     &           pcp1,pcp3,pcp6,pcp24,
     &           snow,maxsta,jstatus)

        call get_sfc_badflag(badflag,istatus)

        n_gauge_noradar = 0

!       Loop through obs and write out precip values (when gauge reports precip)
        do iob = 1,n_obs_b

!           Fill gauge array according to cycle time
            if(ilaps_cycle_time .eq. 3600)then
                pcp_gauge(iob) = pcp1(iob)
                c_field = 'pcp1'
            elseif(ilaps_cycle_time .eq. 10800)then
                pcp_gauge(iob) = pcp3(iob)
                c_field = 'pcp3'
            elseif(ilaps_cycle_time .eq. 21600)then
                pcp_gauge(iob) = pcp6(iob)
                c_field = 'pcp6'
            elseif(ilaps_cycle_time .eq. 86400)then
                pcp_gauge(iob) = pcp24(iob)
                c_field = 'pcp24'
            else
                pcp_gauge(iob) = badflag
                c_field = 'none'
            endif

!           Obtain LAPS i,j at ob location
            call latlon_to_rlapsgrid(lat_s(iob),lon_s(iob),lat,lon
     1                              ,ni,nj,ri,rj
     1                              ,istatus)
            if(istatus.ne.1)goto20

            ilaps(iob) = nint(ri)
            jlaps(iob) = nint(rj)

            if(ilaps(iob) .ge. 1 .and. ilaps(iob) .le. ni .and. 
     1         jlaps(iob) .ge. 1 .and. jlaps(iob) .le. nj      )then

!               Convert from meters to inches
                if(pcp_cmb_m(ilaps(iob),jlaps(iob)) .ne. r_missing_data
     1                                                            )then
                   pcp_laps_in = pcp_cmb_m(ilaps(iob),jlaps(iob)) 
     1                          * in_per_m        
                else
                   pcp_laps_in = r_missing_data
                endif

                if(closest_radar(ilaps(iob),jlaps(iob)) 
     1                          .ne. r_missing_data ) then
                    closest_radar_km = 
     1              closest_radar(ilaps(iob),jlaps(iob)) / 1000.
                else
                    closest_radar_km = -999.
                endif

                if(pcp_gauge(iob) .ge. 0.)then
                   if(pcp_laps_in .ne. r_missing_data)then
                      write(6,11)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),pcp_laps_in     
     1                          ,closest_radar_km,provider(iob)
11                    format(i6,1x,a,2f7.3,f10.1,1x,a11,' RADAR/FG')           
                   else
                      write(6,12)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),lat_s(iob),lon_s(iob)    
     1                          ,closest_radar_km,provider(iob)
12                    format(i6,1x,a,f7.3,2f8.2,f10.1,1x,a11
     1                                               ,' NORADAR/FG')          
                      n_gauge_noradar = n_gauge_noradar + 1
                   endif
                endif

!               Other QC can be done here if needed by setting to badflag
                if(l_accum_bias_ratio .and. n_msg_rdr_bkg .eq. 0)then
                   if(pcp_gauge(iob) .gt. 0. .AND. 
     1                pcp_laps_in .ne. r_missing_data .AND.
     1                pcp_laps_in .gt. 0.)then
                      bias_ratio = pcp_gauge(iob) / pcp_laps_in
                      if(bias_ratio .ge. 0.5 .AND. 
     1                   bias_ratio .le. 2.0      )then
                         pcp_gauge(iob) = bias_ratio
                      else
                         pcp_gauge(iob) = badflag
                      endif
                   else
                      pcp_gauge(iob) = badflag
                   endif
                endif   

            endif

20          continue

        enddo ! iob

        wt_bkg_a = 5e28

        if(n_msg_rdr_bkg .gt. 0 .or. l_gauge_only)then ! do gauge only analysis
            if(n_msg_rdr_bkg .gt. 0)then            
                write(6,*)' Background/radar field has missing points'
            endif

            write(6,*)' Performing gauge only analysis'  
            call precip_barnes_jacket(           c_field              ! I
     1                                           ,ilaps,jlaps         ! I
     1                                           ,pcp_gauge           ! I
     1                                           ,maxsta              ! I
     1                                           ,zero                ! I
     1                                           ,badflag,ni,nj       ! I
     1                                           ,topo,ldf            ! I
     1                                           ,wt_bkg_a            ! I
     1                                           ,pcp_2d_in,istatus)  ! O

            precip_accum_m = pcp_2d_in * meters_per_inch        

            precip_accum_m = max(precip_accum_m,0.)
 
        elseif(.not. l_accum_bias_ratio)then ! analyze gauge values via increments

            write(6,*)' Performing gauge increment analysis'

            pcp_cmb_in = pcp_cmb_m * in_per_m

            call precip_barnes_jacket(           c_field              ! I
     1                                           ,ilaps,jlaps         ! I
     1                                           ,pcp_gauge           ! I
     1                                           ,maxsta              ! I
     1                                           ,pcp_cmb_in          ! I
     1                                           ,badflag,ni,nj       ! I
     1                                           ,topo,ldf            ! I
     1                                           ,wt_bkg_a            ! I
     1                                           ,pcp_2d_in,istatus)  ! O

!           
            precip_accum_m = pcp_2d_in * meters_per_inch        

            precip_accum_m = max(precip_accum_m,0.)
                                                                        
        elseif(l_accum_bias_ratio)then ! perform gauge bias analysis
            write(6,*)' Performing gauge bias analysis'
            one = 1.0
            call precip_barnes_jacket(           c_field              ! I
     1                                           ,ilaps,jlaps         ! I
     1                                           ,pcp_gauge           ! I
     1                                           ,maxsta              ! I
     1                                           ,one                 ! I
     1                                           ,badflag,ni,nj       ! I
     1                                           ,topo,ldf            ! I
     1                                           ,wt_bkg_a            ! I
     1                                           ,bias_anal,istatus)  ! O
            write(6,*)' Max Bias Anal: ', MAXVAL(bias_anal)
            write(6,*)' Min Bias Anal: ', MINVAL(bias_anal)
            precip_accum_m = pcp_cmb_m * bias_anal


        else                ! return blended radar/background (no gauges)
            write(6,*)' Returning radar/background blend (no gauges)'
            precip_accum_m = pcp_cmb_m

        endif

        return
        end

