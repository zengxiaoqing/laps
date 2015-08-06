
        subroutine mean_wind_bunkers(uanl,vanl,topo,imax,jmax,kmax  ! I
     1                              ,heights_3d                     ! I
     1                              ,umean,vmean                    ! O
     1                              ,ushear,vshear                  ! O
     1                              ,ustorm,vstorm,istatus)         ! O

        logical ltest_vertical_grid

        real umean(imax,jmax),vmean(imax,jmax)                    ! O
        real ustorm(imax,jmax),vstorm(imax,jmax)                  ! O
        real ushear(imax,jmax),vshear(imax,jmax)                  ! O
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)            ! I
        real heights_3d(imax,jmax,kmax)                           ! I

        real topo(imax,jmax)                                      ! I

        real sum(imax,jmax)                                       ! L
        real usum(imax,jmax)                                      ! L
        real vsum(imax,jmax)                                      ! L
        integer klow(imax,jmax)                                   ! L
        integer khigh(imax,jmax)                                  ! L

        write(6,*)
        write(6,*)' Calculating Mean Wind (BSM)'

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        do j = 1,jmax
          do i = 1,imax
!            Layer is 0-6 km AGL, denoted from "sfc" to "top"

             klow(i,j) = nint(height_to_zcoord2(topo(i,j)      
     1                       ,heights_3d,imax,jmax,kmax,i,j,istatus))
             if(istatus .ne. 1)return

             khigh(i,j) = nint(height_to_zcoord2(topo(i,j)+6000.
     1                        ,heights_3d,imax,jmax,kmax,i,j,istatus))
             if(istatus .ne. 1)return

             sum(i,j) = 0.
             usum(i,j) = 0.
             vsum(i,j) = 0.

          enddo ! j
        enddo ! i

        do j = 1,jmax
          do i = 1,imax
            do k = 1,khigh(i,j)
              if(uanl(i,j,k) .ne. r_missing_data .and.
     1           vanl(i,j,k) .ne. r_missing_data
     1                        .and. k .ge. klow(i,j))then
                sum(i,j) = sum(i,j) + 1.
                usum(i,j) = usum(i,j) + uanl(i,j,k)
                vsum(i,j) = vsum(i,j) + vanl(i,j,k)
              endif
            enddo ! k
          enddo ! i
        enddo ! j

        do j = 1,jmax
          do i = 1,imax

C            COMPUTE STORM MOTION VECTOR (a la the RUC code)
C            IT IS DEFINED AS 7.5 M/S TO THE RIGHT OF THE 0-6 KM MEAN
C            WIND CONSTRAINED ALONG A LINE WHICH IS BOTH PERPENDICULAR
C            TO THE 0-6 KM MEAN VERTICAL WIND SHEAR VECTOR AND PASSES
C            THROUGH THE 0-6 KM MEAN WIND.  

!            Mean wind through the layer
             umean(i,j) = usum(i,j) / sum(i,j)
             vmean(i,j) = vsum(i,j) / sum(i,j)

!            Shear Vector through the layer
             if(uanl(i,j,klow(i,j)) .ne. r_missing_data .and.
     1          vanl(i,j,klow(i,j)) .ne. r_missing_data       )then
                 ushear(i,j) = uanl(i,j,khigh(i,j))-uanl(i,j,klow(i,j))
                 vshear(i,j) = vanl(i,j,khigh(i,j))-vanl(i,j,klow(i,j))

                 shearspeed = sqrt(ushear(i,j)*ushear(i,j)
     1                            +vshear(i,j)*vshear(i,j))

                 ustorm(i,j) = umean(i,j) + (7.5*vshear(i,j)/shearspeed)
                 vstorm(i,j) = vmean(i,j) - (7.5*ushear(i,j)/shearspeed)

             else
                 write(6,*)' Error in meanwind, missing low level wind'       

             endif


          enddo ! i
        enddo ! j

        return
        end

        subroutine mean_wind(uanl,vanl,topo,imax,jmax,kmax
     1                                  ,umean,vmean
     1                                  ,ustorm,vstorm,istatus)

        logical ltest_vertical_grid

        real umean(imax,jmax),vmean(imax,jmax)                  ! Output
        real ustorm(imax,jmax),vstorm(imax,jmax)                ! Output
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! Input

        real topo(imax,jmax)                                    ! Input

        real sum(imax,jmax)                                     ! Local
        real usum(imax,jmax)                                    ! Local
        real vsum(imax,jmax)                                    ! Local
        integer klow(imax,jmax)                                 ! Local

        write(6,*)
        write(6,*)' Calculating Mean Wind (LSM)'

        if(ltest_vertical_grid('HEIGHT'))then
            khigh = nint(height_to_zcoord(5000.,istatus))
        elseif(ltest_vertical_grid('PRESSURE'))then
            pres_mb = 300.
            pres_pa = pres_mb * 100.
            khigh = nint(zcoord_of_pressure(pres_pa))
        else
            write(6,*)' mean_wind: unknown vertical grid'
            istatus = 0
            return
        endif

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        write(6,*)' Top level of mean wind computation = ',khigh

!       Mean wind (mass weighted) is calculated for whole levels within the 
!       range.
        do j = 1,jmax
          do i = 1,imax
             klow(i,j) =
     1            max(nint(height_to_zcoord(topo(i,j),istatus)),1)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: ERROR in height_to_zcoord'
                 return
             endif
             sum(i,j) = 0.
             usum(i,j) = 0.
             vsum(i,j) = 0.
          enddo ! j
        enddo ! i

        do k = 1,khigh
          do j = 1,jmax
            do i = 1,imax
              if(uanl(i,j,k) .ne. r_missing_data .and.
     1           vanl(i,j,k) .ne. r_missing_data
     1                        .and. k .ge. klow(i,j))then
                sum(i,j) = sum(i,j) + 1.
                usum(i,j) = usum(i,j) + uanl(i,j,k)
                vsum(i,j) = vsum(i,j) + vanl(i,j,k)
              endif
             enddo
          enddo
        enddo

        do j = 1,jmax
          do i = 1,imax

!            Mean wind through the layer
             umean(i,j) = usum(i,j) / sum(i,j)
             vmean(i,j) = vsum(i,j) / sum(i,j)

!            Shear Vector through the layer
             ushear = uanl(i,j,khigh) - uanl(i,j,klow(i,j))
             vshear = vanl(i,j,khigh) - vanl(i,j,klow(i,j))

!            Estimate storm motion of a right moving storm
!            Rotate shear vector by 90 deg, multiply by .15, add to mean wind
             ustorm(i,j) = umean(i,j) ! + 0.15 * vshear
             vstorm(i,j) = vmean(i,j) ! - 0.15 * ushear

          enddo ! i
        enddo ! j

        return
        end
