
        subroutine get_istat_39(t39_k,tb8_k,solar_alt
     1                         ,r_missing_data,ni,nj,istat_39_a)

!       This routine returns the status of stratus cloud detection for each 
!       grid point (containing liquid water).
!       -1 = cloud is not present (of any kind)
!        0 = indeterminate cloud existence
!       +1 = cloud is present (and is stratus)

!       We think this works even with thin cloud scenarios that could cause
!       an underestimation of cloud-top temperature when using 'tb8_k'.

        real*4    t39_k(ni,nj)
        real*4    tb8_k(ni,nj)
        real*4    solar_alt(ni,nj)
        integer*4 istat_39_a(ni,nj)
        integer*4 icount(-1:+1)

        real*4 k_to_c

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

                    if(t39_k(i,j) - tb8_k(i,j) .lt. -4.)then ! Sufficient diff
                        istat_39_a(i,j) = +1

                        if(icount(1) .le. 20)then            ! Log output
                            write(6,11)i,j,t39_c,tb8_c,t39_c-tb8_c
 11                         format(' Stratus present t39/tb8/df: '
     1                            ,2i5,3f8.1)
                        endif 

                    else                                     ! Insuff diff
                        istat_39_a(i,j) = -1

                    endif

                else ! Sun above horizon or possibly cold cloud top: ambiguous
                    istat_39_a(i,j) = 0

                endif

            else ! we have missing data
                istat_39_a(i,j) = 0

            endif

            icount(istat_39_a(i,j)) = icount(istat_39_a(i,j)) + 1

        enddo 
        enddo

        write(6,12)icount(-1),icount(0),icount(+1)
 12     format('  3.9u cloud mask status vals [-1,0,+1] =',3i8)       

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

        real*4    t39_k(ni,nj)
        real*4    tb8_k(ni,nj)
        real*4    solar_alt(ni,nj)
        integer*4 istat_39_lwc_a(ni,nj)
        integer*4 icount(-1:+1)

        real*4 k_to_c

        do ic = -1,+1
            icount(ic) = 0
        enddo ! ic

        do i = 1,ni
        do j = 1,nj
            if(tb8_k(i,j)  .ne. r_missing_data .and.
     1         t39_k(i,j) .ne. r_missing_data      )then
                t39_c = k_to_c(t39_k(i,j))
                tb8_c = k_to_c(tb8_k(i,j))

                if(tb8_c .gt. 0.)then ! "Definitely" warm cloud top
                    istat_39_lwc_a(i,j) = +1

                elseif(solar_alt(i,j) .lt. 0.        )then
                    if(t39_k(i,j) - tb8_k(i,j) .lt. -4.)then ! Sufficient diff
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

        function k_to_c(x)
        real*4 k_to_c
        k_to_c = (x - 273.15)
        return
        end

