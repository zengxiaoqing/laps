c
c**************new routine as adapted at FSL**************************

      subroutine cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d  ! I
     1                          ,istatus)                           ! O

      real*4 pres_3d_pa(ni,nj,nk)                                   ! I
!     real*4 heights_3d(ni,nj,nk)                                   ! I
      real*4 temp_3d(ni,nj,nk)                                      ! I
      real*4 td_3d(ni,nj,nk)                                        ! I
      real*4 u_3d(ni,nj,nk)
      real*4 v_3d(ni,nj,nk)

      real*4 haines_mid_2d(ni,nj)                                   ! L
      real*4 vent_2d(ni,nj)                                         ! L
      real*4 pbl_top_pa(ni,nj)                                      ! L
      real*4 pbl_top_m(ni,nj)                                       ! L
      real*4 pbl_depth_m(ni,nj)                                     ! L

      real*4 pres_3d_mb(ni,nj,nk)                                   ! L

      write(6,*)' Subroutine cpt_fire_fields'

!     Calculate Haines Index

!     pres_3d_mb = pres_3d_pa / 100.
!     call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_mid_2d,ni,nj,nk
!    1                ,850.,700.) 

!     Calculate Fosberg Fireweather Index
!     call fireweatherindex()      

!     Read PBL info from file

!     Calculate Ventilation Index
      call ventilation_index(u_3d,v_3d,pbl_top_pa,ni,nj,nk,heights_3d    ! I
     1                      ,vent_2d,istatus)                            ! O

      return
      end

      subroutine ventilation_index(u_3d,v_3d,pbl_top_pa,ni,nj,nk         ! I
     1                            ,heights_3d                            ! I
     1                            ,vent_2d,istatus)                      ! O

      real*4 heights_3d(ni,nj,nk)                                        ! I
      real*4 u_3d(ni,nj,nk)                                              ! I
      real*4 v_3d(ni,nj,nk)                                              ! I

      real*4 pbl_top_pa(ni,nj)                                           ! I
      real*4 vent_2d(ni,nj)                                              ! O

      write(6,*)' Subroutine ventilation_index'

!     Find boundary-layer depth in meters

!     Call mean_wind subroutine

!     Multiply to obtain ventilation index

      return
      end



        subroutine pbl_mean_wind(uanl,vanl,topo,imax,jmax,kmax    ! I
     1                          ,umean,vmean,istatus)             ! O

        logical ltest_vertical_grid

        real*4 umean(imax,jmax),vmean(imax,jmax)                  ! Output
        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! Input

        real*4 topo(imax,jmax)                                    ! Input

        real*4 sum(imax,jmax)                                     ! Local
        real*4 usum(imax,jmax)                                    ! Local
        real*4 vsum(imax,jmax)                                    ! Local
        integer*4 klow(imax,jmax)                                 ! Local

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

          enddo ! i
        enddo ! j

        return
        end
