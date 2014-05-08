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
        subroutine insert_radar(i4time,cldcv,cld_hts
     1         ,temp_3d,temp_sfc_k,td_sfc_k                          ! I
     1         ,grid_spacing_m,ni,nj,nk,kcloud                       ! I
     1         ,cloud_base,ref_base                                  ! I
!    1         ,tb8_k,sh_3d                                          ! I
     1         ,topo,solar_alt,r_missing_data                        ! I
     1         ,grid_ra_ref,dbz_max_2d                               ! I/O
     1         ,vis_radar_thresh_dbz                                 ! I
     1         ,l_unresolved                                         ! O
     1         ,heights_3d                                           ! I
     1         ,istatus)                                             ! O

!       Insert radar data into cloud grid and reconcile cloud/radar fields.

!       For wideband radar, "resolved" echoes can be added in as cloud if the 
!       echo top is >2000m AGL and the echo is located above an existing 
!       "nominal" cloud base. The nominal cloud base is taken (if present) from 
!       the pre-existing cloud field, either at the horizontal grid location 
!       in question or from a neighboring grid point that has existing cloud.

!       High echo tops with consistent cloud bases are treated as "resolved".
!       Low echo tops are treated as "unresolved" as are echoes that lack
!       pre-existing cloud.

!       For all "unresolved" echoes no cloud is added and the echo is blanked 
!       out.

!       THIS PARAGRAPH SHOULD BE REVIEWED TO ENSURE IT IS CURRENT
!       Narrowband radar is treated similarly, except that the echo top
!       threshold does not apply. Instead, this subroutine will be actively
!       adding clouds (between the pre-existing cloud layers) only for 
!       "trusted" narrowband radars such as NOWRAD. For other narrowband 
!       radars, the radar is blanked out if there is no pre-existing cloud 
!       in the calling routine 'laps_cloud_sub'.

!       Further reconciling of radar and satellite/cloud information is done 
!       later on in subroutine 'insert_vis'.

!       1992 - 2003        Steve Albers

        use mem_namelist, ONLY: echotop_thr_a

        real temp_3d(ni,nj,nk)
        real heights_3d(ni,nj,nk)
        real temp_sfc_k(ni,nj)
        real td_sfc_k(ni,nj)

        real cldcv(ni,nj,kcloud)
        real cld_hts(kcloud)

        real grid_ra_ref(ni,nj,nk)
        real dbz_max_2d(ni,nj)
        real dbz_max_thr(ni,nj)
        real cloud_base(ni,nj)
        real cloud_base_buf(ni,nj) ! Lowest SAO/IR base within search radius
        real topo(ni,nj)
        real solar_alt(ni,nj)

        real echo_top(ni,nj)           ! L
        real echo_top_agl(ni,nj)       ! L
        real echo_top_temp(ni,nj)      ! L
        real echo_agl_thr(ni,nj)       ! L
        real lcl_agl(ni,nj)            ! L

        integer idebug_a(ni,nj)        ! L

!       Cloud not filled in unless radar echo is higher than base calculated
!       with THIS threshold.
        real     thresh_cvr
        parameter (thresh_cvr = 0.10)

        logical l_below_base, l_inserted
        logical l_resolved(ni,nj)
        logical l_unresolved(ni,nj)

!       Calculate Cloud Bases
        unlimited_ceiling = 200000.

!       Echo top threshold
!       echo_agl_thr = 2000.

!       Search radius for cloud bases
        search_radius_m = 60000.

        call get_ref_base_useable(ref_base_useable,istatus)
        if(istatus .ne. 1)return

        itest = min(13,ni)

        do j = 1,nj
        do i = 1,ni

!           Approximation based on a 4 deg F dewpoint depression for each
!           1000 ft of cloud base above the ground
            lcl_agl(i,j) = (temp_sfc_k(i,j) - td_sfc_k(i,j)) * 137.16 ! meters

            if(i .eq. (i/5)*5 .and. j .eq. (j/5)*5)then
                idebug_a(i,j) = 1
            else
                idebug_a(i,j) = 0
            endif

            if(solar_alt(i,j) .lt. 15. .and. 
     1          temp_sfc_k(i,j) .gt. echotop_thr_a(3))then ! warmer evenings
                echo_agl_thr(i,j) = echotop_thr_a(2)
                dbz_max_thr(i,j) = 20.
            else                                           ! colder evenings
                echo_agl_thr(i,j) = echotop_thr_a(1)
                dbz_max_thr(i,j) = ref_base_useable
            endif

            cloud_base(i,j) = unlimited_ceiling
            echo_top(i,j) = r_missing_data
            echo_top_agl(i,j) = r_missing_data

            do k = kcloud-1,1,-1
                if(cldcv(i,j,k  ) .lt. thresh_cvr .and.
     1             cldcv(i,j,k+1) .ge. thresh_cvr       )then
                    cloud_base(i,j) = 0.5 * (cld_hts(k) + cld_hts(k+1))
                endif
            enddo ! k

            cloud_base_buf(i,j) = cloud_base(i,j)

            l_resolved(i,j) = .true.

        enddo ! i
        enddo ! j

        icount_below = 0
        isearch_base = 0
        insert_count_tot = 0
        iwrite_inserted = 0

        isearch_radius = nint(search_radius_m / grid_spacing_m)
        intvl = max((isearch_radius / 4),1)

        write(6,*)' isearch_radius/intvl = ',isearch_radius,intvl

        kp1 = nk

        do k = nk,1,-1 ! Essential that this go downward to detect radar tops
                       ! in time to search for a new cloud base

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

            if(klow .eq. 10000 .or. khigh .eq. -1)then   
                write(6,*)' Skip looping, no cloud levels'
                goto600
            endif

500         write(6,*)' Looping through level ',k,klow,khigh

            do i = 1,ni
            do j = 1,nj

                l_inserted = .false.

                if(temp_3d(i,j,k) .lt. 278.)then
                    ref_thresh = ref_base_useable
                else
                    ref_thresh = ref_base_useable ! Lowered from 20. and 13. 9/94
!                   ref_thresh = 20.              ! Raised back to 13. in 2012          
                endif

!               if(idebug_a(i,j) .eq. 1 .and. k .eq. 29)then
!               if(idebug_a(i,j) .eq. 1)then
!                   write(6,*)'etop 1:',i,j,k
!    1                    ,grid_ra_ref(i,j,k),grid_ra_ref(i,j,kp1)
!               endif

                if(grid_ra_ref(i,j,k) .ge. ref_thresh)then
                    icount_radar_lvl = icount_radar_lvl + 1
                    l_below_base = .false.

!                   if(idebug_a(i,j) .eq. 1)then
!                       write(6,*)'etop 2:',i,j,k
!    1                        ,grid_ra_ref(i,j,k),grid_ra_ref(i,j,kp1)
!                   endif

!                   Test if we are at echo top
                    if(k .eq. nk   .or. 
     1                 grid_ra_ref(i,j,kp1) .lt. ref_thresh)then

                        echo_top(i,j) = heights_3d(i,j,k)
                        echo_top_agl(i,j) = echo_top(i,j) - topo(i,j)
                        echo_top_temp(i,j) = temp_3d(i,j,k)                     

!                       if(idebug_a(i,j) .eq. 1)then
!                           write(6,*)'Echo top found:',i,j,k
!    1                                                  ,echo_top(i,j)
!                       endif

!                       Test if we are below the cloud base
                        if(echo_top(i,j) .lt. cloud_base(i,j))then

!                           Radar Echo Top is below analyzed cloud base
!                           Various tests attempting to resolve...
!                           Search for Modified Cloud Base, i.e. other neighboring
!                           bases lower than the current buffer value

                            ilow =  max(i-isearch_radius,1)
                            ihigh = min(i+isearch_radius,ni)
                            jlow =  max(j-isearch_radius,1)
                            jhigh = min(j+isearch_radius,nj)

                            do jj = jlow,jhigh,intvl
                            do ii = ilow,ihigh,intvl
                                base_buf_lcl = max(cloud_base(ii,jj)
     1                                            ,lcl_agl(i,j))
                                cloud_base_buf(i,j)
     1                    = min(base_buf_lcl,cloud_base_buf(i,j))
                            enddo ! ii
                            enddo ! jj

                            if(cloud_base_buf(i,j) .lt. echo_top(i,j)
     1                                       .AND. 
     1                         echo_top_agl(i,j) .gt. echo_agl_thr(i,j)
     1                                       .AND. 
     1                         echo_top_temp(i,j) .le. 278.15   ! Cold Precip          
     1                                       .AND. 
     1                         dbz_max_2d(i,j) .ge. dbz_max_thr(i,j)
     1                                                             )then

                              isearch_base = isearch_base + 1
                              if(idebug_a(i,j) .eq. 1)then 
                                write(6,71)i,j,k
     1                                    ,nint(echo_top(i,j))
     1                                    ,nint(cloud_base(i,j))
     1                                    ,nint(cloud_base_buf(i,j))
71                              format(' Rdr Top > Bse ',2i4,i3,3i7
     1                                ,' Resolved')
                              endif

!                             l_resolved is still TRUE

                            else ! Unresolved base
                              l_resolved(i,j) = .false.

                            endif ! Resolved Base tests (nearby base, high echo
                                  !                      top, and cold)

                        endif ! Below Cloud Base (try to resolve)

                    endif ! At Echo Top

!                   Loop through range of cloud grid levels for this LAPS level
                    do kk = klow,khigh
!                       Insert radar if we are above cloud base
                        if(cld_hts(kk) .gt. cloud_base_buf(i,j)
     1                                .and. (l_resolved(i,j).eqv..true.)
     1                                                           )then
                            cldcv(i,j,kk) = 1.0
                            insert_count_lvl = insert_count_lvl + 1
                            insert_count_tot = insert_count_tot + 1
                            l_inserted = .true.
                        else ! Radar Echo below cloud base
                            l_below_base = .true.
                        endif
                    enddo ! kk

                    if(l_below_base)then
                        icount_below = icount_below + 1
                        if(idebug_a(i,j) .eq. 1)then 
                            write(6,81)i,j,k,nint(cld_hts(klow))
     1                                  ,nint(cloud_base_buf(i,j))
81                          format(' Rdr     < Bse ',2i4,i3,2i7)
                        endif
                    endif ! below base

                endif ! Reflectivity > thresh

                if( (l_inserted .eqv. .true.) .AND. 
     1              (idebug_a(i,j) .eq. 1)          )then
                    write(6,591)i,j,k,l_resolved(i,j)
     1                         ,dbz_max_2d(i,j),echo_top(i,j)
     1                         ,echo_top_temp(i,j)
591                 format(' Inserted radar',2i4,i3,l2
     1                                      ,1x,f8.1,f8.0,f8.2)
                    if(echo_top(i,j) .eq. r_missing_data)then
                        write(6,*)grid_ra_ref(i,j,:)
                    endif
                    iwrite_inserted = iwrite_inserted + 1
                endif

            enddo ! j
            enddo ! i

            write(6,592)k,klow,khigh
     1          ,icount_radar_lvl,insert_count_lvl,insert_count_tot
592         format(' Inserted radar',3i3,3i8)
 
            kp1 = k

600         continue

        enddo ! k (LAPS grid level)

        write(6,*)' Total cloud grid points modified by radar = '
     1                                          ,insert_count_tot

        do i = 1,ni
        do j = 1,nj

            if(echo_top(i,j) .ne. r_missing_data)then
                if(echo_top_agl(i,j) .lt. echo_agl_thr(i,j) .and. 
     1             (l_resolved(i,j).eqv..true.)       )then

!                   Should we set these to unresolved here or better yet above?
                    if(echo_top_temp(i,j) .ge. 278.15)then
                      iwarm_resolved = 1
                      l_resolved(i,j) = .false.
                    else
                      iwarm_resolved = 0
                    endif

                    if(idebug_a(i,j) .eq. 1)then 
                      write(6,610)i,j       
     1                        ,nint(echo_top_agl(i,j))
     1                        ,nint(echo_agl_thr(i,j))
     1                        ,nint(echo_top(i,j))
     1                        ,nint(cloud_base(i,j))
     1                        ,nint(cloud_base_buf(i,j))
     1                        ,echo_top_temp(i,j),iwarm_resolved
 610                  format('CLD_RDR - Low echo top yet resolved '
     1                    ,3i5,' ET',2i5,' Base',2i7,f7.1,i3)
                    endif
                endif
            endif

            if(l_resolved(i,j) .eqv. .false.)then

!               Reconcile radar and satellite 
!               (Block out "unresolved" radar echoes)

                if(dbz_max_2d(i,j) .lt. vis_radar_thresh_dbz)then
!                   Blank out radar
                    if(idebug_a(i,j) .eq. 1)then 
                      write(6,601)i,j,dbz_max_2d(i,j)
     1                         ,nint(echo_top_agl(i,j))
     1                         ,nint(echo_agl_thr(i,j))
     1                         ,echo_top_temp(i,j)
601                   format('CLD_RDR - insert_radar: '
     1                         ,'Blank out radar < '
     1                         ,2x,2i4,f6.1,' dbz',i6,1x,i6,' agl',f8.1)       
                    endif

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                else

!                   This situation is undesirable because we want to
!                   avoid cases where the radar reflectivity is high
!                   but is "unresolved" (e.g. analyzed (SAO/IR) clouds at 
!                   this point in) the analysis. We should check as to why 
!                   it is unresolved (e.g. no clouds analyzed). For now we 
!                   blank out the radar, but we could also create a cloud at 
!                   the radar echo location.

                    if(idebug_a(i,j) .eq. 1)then 
                      write(6,602)i,j,dbz_max_2d(i,j)
     1                         ,nint(echo_top_agl(i,j))       
     1                         ,nint(echo_agl_thr(i,j))
     1                         ,echo_top_temp(i,j)
602                   format('CLD_RDR - insert_radar: '
     1                         ,'Blank out radar > *'
     1                         ,1x,2i4,f6.1,' dbz',i6,1x,i6,' agl',f8.1)
                    endif

                    if(.true.)then ! Block out the radar
                        do k = 1,nk
                            grid_ra_ref(i,j,k) = ref_base
                        enddo ! k
                        dbz_max_2d(i,j) = ref_base
                    endif

                endif

            else
                if(idebug_a(i,j) .eq. 1)then 
                    if(echo_top(i,j) .ne. r_missing_data)then
                        write(6,612)i,j,nint(cloud_base(i,j))
     1                             ,nint(echo_top(i,j))
612                     format('Resolved ',2i4,1x,2i6)                    
                    endif
                endif

            endif ! l_resolved

            if(echo_top(i,j) .ne. r_missing_data)then
                if(echo_top_temp(i,j) .ge. 278.15  .AND.  
     1             dbz_max_2d(i,j)    .ge. ref_base_useable)then   
                    if(idebug_a(i,j) .eq. 1)then
                        write(6,*)' Warm echo remaining at',i,j
     1                           ,dbz_max_2d(i,j),echo_top_agl(i,j)
     1                           ,echo_top_temp(i,j)
                    endif
                endif
            endif

        enddo ! j
        enddo ! i

!       Blank out echoes less than "useable"
        where (grid_ra_ref .lt. ref_base_useable)grid_ra_ref = ref_base       
        where (dbz_max_2d  .lt. ref_base_useable)dbz_max_2d  = ref_base       

999     l_unresolved = .not. (l_resolved)

        return

        end
