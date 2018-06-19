
        subroutine blend_gauge_data(    i4time,ni,nj,maxsta        ! I
     1                                 ,r_missing_data             ! I
     1                                 ,lat,lon,topo,ldf           ! I
     1                                 ,pcp_bkg_m                  ! I
     1                                 ,ilaps_cycle_time           ! I
     1                                 ,closest_radar,istat_radar  ! I
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
        real topo(ni,nj)
        real ldf(ni,nj)

        real pcp_gauge(maxsta)    ! relevant time of accumulation is selected
        real pcp_laps_in_a(maxsta)! radar/first guess precip of paired obs        
        real gauge_pairs_a(maxsta)! rain gauge reading of paired observations
        character*5 c_field
        character*24 btime

        integer ilaps(maxsta),jlaps(maxsta)

        logical l_regression       ! Apply regression of radar/fg & gauges

        logical l_accum_bias_ratio ! Apply bias correction as a constant ratio
                                   ! (given we are analyzing radar/fg & gauges)

        logical l_gauge_only       ! Gauge only analysis

        logical l_qc               ! QC of radar gauge pairs

        l_gauge_only = (.not. l_accum_fg) .AND. (.not. l_accum_radar) 
     1                                    .AND.        l_accum_gauge

!       If both 'l_regression' and 'l_bias_ratio' are F then the gauge increment
!       analysis can be performed.
        l_regression = .false.
        l_accum_bias_ratio = .false. 
!       l_accum_bias_ratio = (.not. l_regression) 
        
!       Combine background and radar field given radar gap areas
        n_radar = 0
        n_bkg = 0
        n_msg_rdr_bkg = 0

        idebug = 0 ! extra output for debugging (0,1)

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

        if(n_radar + n_bkg .eq. 0)l_gauge_only = .true.

        write(6,*)
        write(6,*)' Subroutine blend_gauge_data (1hr pcp inches)...'

        write(6,*)' l_accum_fg    = ',l_accum_fg
        write(6,*)' l_accum_radar = ',l_accum_radar
        write(6,*)' l_accum_gauge = ',l_accum_gauge
        write(6,*)' l_gauge_only  = ',l_gauge_only
        write(6,*)
        write(6,*)' l_accum_bias_ratio  = ',l_accum_bias_ratio
        write(6,*)' l_regression        = ',l_regression
        write(6,*) 
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
        n_pairs = 0

!       Loop through obs and write out precip values (when gauge reports precip)
        do iob = 1,n_obs_b

!           Fill gauge array according to cycle time
            if(ilaps_cycle_time .le. 3600)then
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

                if(idebug .gt. 0 .and. iob .le. 100)then
                   write(6,*)'gauge in domain',iob,pcp_gauge(iob)
                endif

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
                      if(pcp_laps_in    .gt. 0. .and. 
     1                   pcp_gauge(iob) .gt. 0.      )then

                        l_qc = .true.
                        if(pcp_gauge(iob) .gt. .20 .and.
     1                    pcp_laps_in / pcp_gauge(iob) .lt. 0.1)then
                          l_qc = .false.
                        endif

                        if(l_qc)then
                          n_pairs = n_pairs + 1
                          pcp_laps_in_a(n_pairs) = pcp_laps_in
                          gauge_pairs_a(n_pairs) = pcp_gauge(iob)

                          write(6,12)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),pcp_laps_in     
     1                          ,closest_radar_km,provider(iob)
12                        format(i6,1x,a,2f7.3,f10.1,1x,a11
     1                          ,' RADAR/FG VALID PAIR')   

                        else
                          write(6,13)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),pcp_laps_in     
     1                          ,closest_radar_km,provider(iob)
13                        format(i6,1x,a,2f7.3,f10.1,1x,a11
     1                          ,' RADAR/FG QC FAILED')           

                        endif

                      else
                        write(6,14)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),pcp_laps_in     
     1                          ,closest_radar_km,provider(iob)
14                      format(i6,1x,a,2f7.3,f10.1,1x,a11,' RADAR/FG')           

                      endif
                   else
                      write(6,15)iob,stations(iob)(1:10)
     1                          ,pcp_gauge(iob),lat_s(iob),lon_s(iob)    
     1                          ,closest_radar_km,provider(iob)
15                    format(i6,1x,a,f7.3,2f8.2,f10.1,1x,a11
     1                                               ,' NORADAR/FG')          
                      n_gauge_noradar = n_gauge_noradar + 1
                   endif
                endif

!               Other QC can be done here if needed by setting to badflag
                if(l_accum_bias_ratio .and. n_msg_rdr_bkg .eq. 0
     1             .and. (.not. l_regression)                    )then
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

!       60% threshold allowed to be missing
        frac_msg_rdr_bkg = float(n_msg_rdr_bkg) / float(ni*nj)

        write(6,*)' n_msg_rdr_bkg = ',n_msg_rdr_bkg
        write(6,*)' frac_msg_rdr_bkg = ',frac_msg_rdr_bkg
        write(6,*)' l_gauge_only  = ',l_gauge_only

!       if(n_msg_rdr_bkg .gt. 0 .or. l_gauge_only)then ! do gauge only analysis
        if(frac_msg_rdr_bkg .gt. .60 .or. l_gauge_only)then ! do gauge only analysis
            if(istat_radar .eq. 0)then
                write(6,*)' Radar data unavailable'
            endif              
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
 
        elseif( (.not. l_accum_bias_ratio) .and. l_accum_gauge
     1           .and. (.not. l_regression) )then 

!           Analyze gauge values via increments

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
                                                                        
        elseif(l_accum_gauge)then 

!           Perform gauge bias analysis

            write(6,*)' Performing gauge bias analysis'

            if(l_accum_bias_ratio)then
                one = 1.0
                call precip_barnes_jacket(       c_field              ! I
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
                write(6,*)' multiplying precip by bias analysis'
                precip_accum_m = pcp_cmb_m * bias_anal

            else ! do regression of radar/gauge pairs (pcp_gauge,pcp_laps_in_a)
                if(n_pairs .gt. 0)then
                    write(6,*)' Do regression of radar/gauge pairs '
     1                       ,n_pairs
                    call regress_precip(n_pairs,pcp_laps_in_a
     1                                 ,gauge_pairs_a
     1                                 ,a_t,b_t,rbar,gbar,istatus)
                    if(istatus .eq. 1)then
                        write(6,*)' Apply regression to radar/fg'
                        precip_accum_m = (pcp_cmb_m * a_t) 
     1                                 + (b_t * meters_per_inch)
                        precip_accum_m = max(precip_accum_m,0.)
                        where(pcp_cmb_m .eq. 0.)precip_accum_m = 0.
                    else
                        write(6,*)' regression not applied'
                        write(6,*)
     1                   ' Returning radar/background blend (no gauges)'
                        precip_accum_m = pcp_cmb_m
                    endif
                else
                    write(6,*)' no gauge/radar pairs for regression'
                endif

            endif

        else                ! return blended radar/background (no gauges)
            write(6,*)' Returning radar/background blend (no gauges)'
            precip_accum_m = pcp_cmb_m

        endif

        return
        end

