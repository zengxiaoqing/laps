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
        subroutine insert_radar(i4time,cldcv,cld_hts
     1         ,temp_3d,temp_sfc_k,grid_spacing_m,ni,nj,nk,kcloud
     1         ,cloud_base,cloud_base_buf,ref_base
     1         ,grid_ra_ref,dbz_max_2d,vis_radar_thresh_dbz,l_unresolved
     1         ,heights_3d,istatus)

!       Insert radar data into cloud grid

!       Feb 1992        Steve Albers

        real*4 temp_3d(ni,nj,nk)
        real*4 heights_3d(ni,nj,nk)
        real*4 temp_sfc_k(ni,nj)

        real*4 cldcv(ni,nj,kcloud)
        real*4 cld_hts(kcloud)

        real*4 grid_ra_ref(ni,nj,nk)
        real*4 dbz_max_2d(ni,nj)
        real*4 cloud_base(ni,nj)
        real*4 cloud_base_buf(ni,nj) ! Lowest SAO/IR base within search radius

!       Cloud not filled in unless radar echo is higher than base calculated
!       with THIS threshold.
        real*4     thresh_cvr
        parameter (thresh_cvr = 0.20)

        logical l_below_base
        logical l_unresolved(ni,nj)

!       Calculate Cloud Bases
        unlimited_ceiling = 200000.

        do j = 1,nj
        do i = 1,ni
            cloud_base(i,j) = unlimited_ceiling

            do k = kcloud-1,1,-1
                if(cldcv(i,j,k  ) .lt. thresh_cvr .and.
     1           cldcv(i,j,k+1) .ge. thresh_cvr)then
                    cloud_base(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
            enddo ! k

            cloud_base_buf(i,j) = cloud_base(i,j)

            l_unresolved(i,j) = .false.

        enddo ! i
        enddo ! j

        icount_below = 0
        isearch_base = 0
        insert_count_tot = 0

        isearch_radius = nint(120000. / grid_spacing_m)
        intvl = max((isearch_radius / 4),1)

        write(6,*)' isearch_radius/intvl = ',isearch_radius,intvl

        do k = nk,1,-1 ! Essential that this go downward to detect radar tops
                       ! in time to search for a new cloud base

            kp1 = min(k+1,nk)

            icount_radar_lvl = 0
            insert_count_lvl = 0

            klow = 10000
            khigh = -1

           !Find range of cloud grid levels for this LAPS grid level
           !Find lowest level
            do kk = 1,kcloud
c               write(6,*)k,kk,height_to_zcoord(cld_hts(kk),istatus)
                if(nint(height_to_zcoord(cld_hts(kk),istatus)) .eq. k
     1                                                        )then
c                   write(6,*)' klow = ',kk
                    klow = kk
                    goto400
                endif
            enddo
400         do kk = kcloud,1,-1
                if(nint(height_to_zcoord(cld_hts(kk),istatus)) .eq. k
     1                                                        )then
c                   write(6,*)' khigh = ',kk
                    khigh = kk
                    goto500
                endif
            enddo

            if(istatus .ne. 1)then
                write(6,*)' Error in insert_radar'
                return
            endif

            if(klow .eq. 10000 .or. khigh .eq. -1)goto600

500         do i = 1,ni
            do j = 1,nj

                if(temp_3d(i,j,k) .lt. 278.)then
                    ref_thresh = 0.0
                else
                    ref_thresh = 0.0 ! Lowered from 20. and 13. 9/94
                endif

                if(grid_ra_ref(i,j,k) .gt. ref_thresh)then
                    icount_radar_lvl = icount_radar_lvl + 1
                    l_below_base = .false.

!                   Test if we are at echo top
                    if(k .eq. nk   .or. 
     1                 grid_ra_ref(i,j,kp1) .lt. ref_thresh)then

!                       Test if we are below the cloud base
                        if(heights_3d(i,j,k) .lt. cloud_base_buf(i,j)
     1                                                         )then

!                           Radar Echo Top is below analyzed cloud base
!                           Search for Modified Cloud Base, i.e. other neighboring
!                           bases lower than the current buffer value

                            ilow =  max(i-isearch_radius,1)
                            ihigh = min(i+isearch_radius,ni)
                            jlow =  max(j-isearch_radius,1)
                            jhigh = min(j+isearch_radius,nj)

                            do jj = jlow,jhigh,intvl
                            do ii = ilow,ihigh,intvl
                                cloud_base_buf(i,j)
     1                    = min(cloud_base(ii,jj),cloud_base_buf(i,j))
                            enddo ! i
                            enddo ! j

                            if(cloud_base_buf(i,j) .lt. 
     1                                  heights_3d(i,j,k))then

                              isearch_base = isearch_base + 1
                              if(isearch_base .lt. 50)then ! limit log output
                                write(6,71)i,j,k
     1                                    ,nint(heights_3d(i,j,k))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
71                              format(' Rdr Top > Bse ',3i3,3i7
     1                                ,' Resolved')
                              endif

                            else ! Probably Unresolved base
                                write(6,72)i,j,k
     1                                    ,nint(heights_3d(i,j,k))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
72                              format(' Rdr Top < Bse ',3i3,3i7
     1                                ,' Prob Unresolved - CLD_RDR')

                                if(cloud_base(i,j)
     1                                      .eq. unlimited_ceiling)then
!                                   We will want to reconcile cloud/radar
                                    l_unresolved(i,j) = .true.
                                endif

                            endif

                        endif ! Below Cloud Base

                    endif ! At Echo Top

!                   Loop through range of cloud grid levels for this LAPS level
                    do kk = klow,khigh
!                       Insert radar if we are above cloud base
                        if(cld_hts(kk) .gt. cloud_base_buf(i,j)
!    1                                .and. .not. l_unresolved(i,j)
     1                                                           )then
                            cldcv(i,j,kk) = 1.0
                            insert_count_lvl = insert_count_lvl + 1
                            insert_count_tot = insert_count_tot + 1
                        else ! Radar Echo below cloud base
                            l_below_base = .true.
                        endif
                    enddo ! kk

                    if(l_below_base)then
                        icount_below = icount_below + 1

                        if(icount_below .le. 100)then
                            write(6,81)i,j,k,nint(cld_hts(klow))
     1                                  ,nint(cloud_base_buf(i,j))
81                          format(' Rdr     < Bse ',3i3,2i7)
                        endif

                    endif

                endif ! Reflectivity > thresh

            enddo ! j
            enddo ! i

            write(6,591)k,klow,khigh
     1          ,icount_radar_lvl,insert_count_lvl,insert_count_tot
591         format(' Inserted radar',3i4,3i8)

600         continue
        enddo ! k

        do i = 1,ni
        do j = 1,nj
            if(l_unresolved(i,j))then

!               Reconcile radar and satellite
                if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then
!                   Blank out radar
                    write(6,601)i,j,dbz_max_2d(i,j)
601                 format(' CLD_RDR - insert_radar: '
     1                             ,'Blank out radar < '
     1                             ,2x,2i4,f6.1,' dbz')

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                else

!                   This situation is undesirable because we want to
!                   avoid cases where the radar reflectivity is high
!                   but there are no analyzed (SAO/IR) clouds at this point in
!                   the analysis. We should check as to why there were no
!                   clouds analyzed. For now we blank out the radar, but we
!                   could also create a cloud at the radar echo location.

                    write(6,602)i,j,dbz_max_2d(i,j)
602                 format(' CLD_RDR - insert_radar: '
     1                             ,'Blank out radar >  *'
     1                             ,1x,2i4,f6.1,' dbz')

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                endif

            endif
        enddo ! j
        enddo ! i

        write(6,*)' Total cloud grid points modified by radar = '
     1                                          ,insert_count_tot

999     return

        end
