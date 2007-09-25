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

        subroutine get_istat_39(t39_k,tb8_k,solar_alt,r_missing_data    ! I
     1                         ,rlaps_land_frac,ni,nj                   ! I
     1                         ,static_albedo                           ! I
     1                         ,istat_39_a)                             ! O

!       This routine returns the status of stratus cloud detection for each 
!       grid point (containing liquid water).
!       -1 = cloud is not present (of any kind)
!        0 = indeterminate cloud existence
!       +1 = cloud is present (and is stratus)

!       We think this works even with thin cloud scenarios that could cause
!       an underestimation of cloud-top temperature when using 'tb8_k'.

        real    t39_k(ni,nj)
        real    tb8_k(ni,nj)
        real    solar_alt(ni,nj)
        real    rlaps_land_frac(ni,nj)
        real    static_albedo(ni,nj)          ! Static albedo database
        integer istat_39_a(ni,nj)
        integer istat_39_buff(ni,nj)
        integer icount(-1:+1)

        real k_to_c

        write(6,*)' Subroutine get_istat_39...'

        do ic = -1,+1
            icount(ic) = 0
        enddo ! ic

        do i = 1,ni
        do j = 1,nj
            if(tb8_k(i,j) .ne. r_missing_data .and.
     1         t39_k(i,j) .ne. r_missing_data      )then
                t39_c = k_to_c(t39_k(i,j))
                tb8_c = k_to_c(tb8_k(i,j))

                if(tb8_c          .gt. -10. ! Cloud-top "Definitely" warm
     1                                      ! enough for possibility of liquid
     1                                      ! water
     1                     .AND.
     1             solar_alt(i,j) .lt. 0.        )then

                    t39_diff = t39_k(i,j) - tb8_k(i,j)

                    istat_39_a(i,j) = 0                       ! Indeterminate

                    if(
     1                  (t39_diff .le. -2.5 .and. 
     1                                    rlaps_land_frac(i,j) .gt. 0.5)      
     1                                 .OR.
     1                  (t39_diff .le. -2.0 .and. 
     1                                    rlaps_land_frac(i,j) .le. 0.5)      
     1                                             )then ! Sufficient diff     
                        istat_39_a(i,j) = +1

                        if(icount(1) .le. 20)then             ! Log output
                            write(6,11)i,j,t39_c,tb8_c,t39_c-tb8_c
 11                         format(' Stratus present t39/tb8/df: '
     1                            ,2i5,3f8.1)
                        endif 
                    endif

                    if(t39_diff .gt. +0.0 .and. t39_diff .lt. +2.5
     1                          .AND. tb8_c .gt. 0.    
     1                          .AND. rlaps_land_frac(i,j) .gt. 0.5)then ! land
                        istat_39_a(i,j) = 0 ! -1              ! Warm - no cloud?
                    endif

                    if(t39_diff .gt. +0.0 .and. t39_diff .lt. +2.5
     1                          .AND. tb8_c .gt. 0.    
     1                          .AND. rlaps_land_frac(i,j) .lt. 0.5)then ! sea
                        istat_39_a(i,j) = 0 ! -1              ! Warm - no cloud?
                    endif

                else ! Sun above horizon or possibly cold cloud top: ambiguous
                    istat_39_a(i,j) = 0 ! Indeterminate

                endif

            else ! we have missing data
                istat_39_a(i,j) = 0     ! Indeterminate

            endif

!           Screen points according to static albedo
            if(static_albedo(i,j) .eq. r_missing_data .or.
     1         static_albedo(i,j) .gt. 0.18                  )then
                istat_39_a(i,j) = 0     ! Indeterminate
            endif

            icount(istat_39_a(i,j)) = icount(istat_39_a(i,j)) + 1

        enddo 
        enddo

        write(6,12)icount(-1),icount(0),icount(+1)
 12     format('  3.9u cloud mask status vals [-1,0,+1] =',3i8)       

!       Filter out isolated points
        istat_39_buff = istat_39_a

        do i = 2,ni-1
        do j = 2,nj-1
            if(istat_39_buff(i,j) .eq. 1)then
!               Count neighbors
                n_neighbors = 0
                do ii = i-1,i+1
                do jj = j-1,j+1
                    if(istat_39_buff(ii,jj) .eq. 1)then
                        n_neighbors = n_neighbors + 1
                    endif
                enddo ! jj
                enddo ! ii

                if(n_neighbors .gt. 2)then
!                   Keep if at least 2 neighbors (plus central pt) are there
                    istat_39_a(i,j) = istat_39_buff(i,j) ! 1
                else
!                   Throw out if isolated point or just one neighbor
                    istat_39_a(i,j) = 0
                endif
            else
                istat_39_a(i,j) = istat_39_buff(i,j)
            endif
        enddo ! j
        enddo ! i

        return
        end 


        subroutine get_istat_39_lwc(t39_k,tb8_k,solar_alt
     1                             ,r_missing_data,ni,nj,istat_39_lwc_a)

!       This routine returns the status of lwc detection GIVEN that a cloud 
!       has been shown to be present by OTHER means.
!       -1 = lwc is present at cloud-top
!       0  = indeterminate
!       +1 = lwc is absent at cloud-top (cloud-ice at cloud-top)

!       We think this works even with thin cloud scenarios that could cause
!       an underestimation of cloud-top temperature when using 'tb8_k'.

        real    t39_k(ni,nj)
        real    tb8_k(ni,nj)
        real    solar_alt(ni,nj)
        integer istat_39_lwc_a(ni,nj)
        integer icount(-1:+1)

        real k_to_c

        do ic = -1,+1
            icount(ic) = 0
        enddo ! ic

        do i = 1,ni
        do j = 1,nj
            if(tb8_k(i,j)  .ne. r_missing_data .and.
     1         t39_k(i,j)  .ne. r_missing_data      )then
                t39_c = k_to_c(t39_k(i,j))
                tb8_c = k_to_c(tb8_k(i,j))

                if(tb8_c .gt. 0.)then ! "Definitely" warm cloud top
                    istat_39_lwc_a(i,j) = +1

                elseif(solar_alt(i,j) .lt. 0.        )then
                    if(t39_k(i,j) - tb8_k(i,j) .lt. -3.)then ! Sufficient diff
                        istat_39_lwc_a(i,j) = +1
                    else                                     ! Insuff diff
                        istat_39_lwc_a(i,j) = -1
                    endif

                else ! Sun above horizon AND possibly cold cloud top: ambiguous
                    istat_39_lwc_a(i,j) = 0

                endif

            else ! we have missing data
                istat_39_lwc_a(i,j) = 0

            endif

            icount(istat_39_lwc_a(i,j)) = 
     1      icount(istat_39_lwc_a(i,j)) + 1       

        enddo 
        enddo

        write(6,12)icount(-1),icount(0),icount(+1)
 12     format('  3.9u cloud lwc status vals [-1,0,+1] =',3i8)       

        return
        end 

