
        subroutine get_istat_39(t39_k,tb8_k,solar_alt
     1                         ,r_missing_data,ni,nj,istat_39_a)

!       This routine returns the status of cloud detection for each grid point.
!       -1 = cloud is not present
!        0 = indeterminate cloud existence
!       +1 = cloud is present

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

                    if(tb8_k(i,j) - t39_k(i,j) .lt. -4.)then ! Sufficient diff
                        istat_39_a(i,j) = +1
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

        do ic = -1,+1
            write(6,*)' status value / number of grid points ',ic
     1               ,icount(ic)
        enddo ! ic

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
                    if(tb8_k(i,j) - t39_k(i,j) .lt. -4.)then ! Sufficient diff
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

            icount(istat_39_lwc_a(i,j)) = icount(istat_39(i,j)) + 1       

        enddo 
        enddo

        do ic = -1,+1
            write(6,*)' status value / number of grid points ',ic
     1               ,icount(ic)
        enddo ! ic

        return
        end 

        function k_to_c(x)
        real*4 k_to_c
        k_to_c = (x - 273.15)
        return
        end

