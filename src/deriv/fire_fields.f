c
c**************new routine as adapted at FSL**************************

      subroutine fire_fields(ni,nj,nk,temp_3d,td_3d                 ! I
!    1                      ,heights_3d                             ! I
     1                      ,u_3d,v_3d                              ! I
     1                      ,r_missing_data,i4time                  ! I
     1                      ,istatus)                               ! O

      integer MAX_FIELDS
      parameter (MAX_FIELDS = 3)

      character*3 var_a(MAX_FIELDS)
      character*125 comment_a(MAX_FIELDS)
      character*10  units_a(MAX_FIELDS)
      character*31 ext

!     real*4 heights_3d(ni,nj,nk)                                   ! I
      real*4 temp_3d(ni,nj,nk)                                      ! I
      real*4 td_3d(ni,nj,nk)                                        ! I
      real*4 u_3d(ni,nj,nk)
      real*4 v_3d(ni,nj,nk)

      real*4 t_sfc_k(ni,nj)                                         ! I?
      real*4 rh_sfc(ni,nj)                                          ! I?
      real*4 p_sfc_pa(ni,nj)                                        ! I?

      real*4 pbl_top_pa(ni,nj)                                      ! L
      real*4 pres_3d_pa(ni,nj,nk)                                   ! L
      real*4 haines_mid_2d(ni,nj)                                   ! L
      real*4 vent_2d(ni,nj)                                         ! L
      real*4 fosberg_2d(ni,nj)                                      ! L

      write(6,*)' Subroutine fire_fields (under construction)'

      write(6,*)' Read in pressure of PBL top'
      pbl_top_pa = r_missing_data      ! We can read PBL in here
      istat_pbl = 0

      call get_pres_3d(i4time,ni,nj,nk,pres_3d_pa,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Could not obtain pres_3d'
          return
      endif

      call cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d        ! I
!    1                          ,heights_3d                         ! I
     1                          ,u_3d,v_3d                          ! I
     1                          ,pbl_top_pa,istat_pbl               ! I
     1                          ,t_sfc_k                            ! I
     1                          ,rh_sfc                             ! I
     1                          ,p_sfc_pa                           ! I
     1                          ,fosberg_2d                         ! O
     1                          ,haines_mid_2d                      ! O
     1                          ,vent_2d                            ! O
     1                          ,istatus)                           ! O

      if(istatus .eq. 1)then ! write out LFR fire fields file
          write(6,*)' Write out LFR file'
!         Fields are: fosberg_2d, haines_mid_2d, vent_2d
          ext = 'lfr'
!         call put_laps_multi_2d(i4time,ext,var_a,units_a,
!    1                           comment_a,out_array_3d,NX_L,NY_L,3,
!    1                           istatus)

      else
          write(6,*)' Skipping write of LFR file'

      endif

      return
      end

       subroutine cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d  ! I
!    1                           ,heights_3d                         ! I
     1                           ,u_3d,v_3d                          ! I
     1                           ,pbl_top_pa,istat_pbl               ! I
     1                           ,t_sfc_k                            ! I
     1                           ,rh_sfc                             ! I
     1                           ,p_sfc_pa                           ! I
     1                           ,fosberg_2d                         ! O
     1                           ,haines_mid_2d                      ! O
     1                           ,vent_2d                            ! O
     1                           ,istatus)                           ! O

       real*4 pres_3d_pa(ni,nj,nk)                                   ! I
!      real*4 heights_3d(ni,nj,nk)                                   ! I
       real*4 temp_3d(ni,nj,nk)                                      ! I
       real*4 td_3d(ni,nj,nk)                                        ! I
       real*4 u_3d(ni,nj,nk)
       real*4 v_3d(ni,nj,nk)

       real*4 pbl_top_pa(ni,nj)                                      ! I
       real*4 pbl_top_m(ni,nj)                                       ! L
       real*4 pbl_depth_m(ni,nj)                                     ! L

       real*4 t_sfc_k(ni,nj)
       real*4 rh_sfc(ni,nj)
       real*4 p_sfc_pa(ni,nj)
       real*4 p_sfc_mb(ni,nj)

       real*4 haines_mid_2d(ni,nj)                                   ! O
       real*4 vent_2d(ni,nj)                                         ! O
       real*4 fosberg_2d(ni,nj)                                      ! O

       real*4 pres_3d_mb(ni,nj,nk)                                   ! L

       write(6,*)' Subroutine cpt_fire_fields...'

!      Calculate Haines Index
       write(6,*)' Calculate Mid-Level Haines Index'
       pres_3d_mb = pres_3d_pa / 100.
       call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_mid_2d,ni,nj,nk
     1                 ,850.,700.) 

!      Calculate Fosberg Fireweather Index
       write(6,*)' Calculate Fosberg Fireweather Index'
!      Note that 3d inputs are presently anticipated as sigma coordinates
!      call fireweatherindex(t_sfc_k,rh_sfc,p_sfc_mb,u10,v10
!    1                      ,fosberg_2d)      

!      Calculate Ventilation Index
       if(istat_pbl .eq. 1)then
           write(6,*)' Calculate Ventilation Index'
           call ventilation_index(u_3d,v_3d,pbl_top_pa,ni,nj,nk           ! I
!    1                           ,heights_3d                              ! I
     1                           ,vent_2d,istatus)                        ! O
       else
           write(6,*)' Skip ventilation index due to no PBL'

       endif

       return
       end

       subroutine ventilation_index(u_3d,v_3d,pbl_top_pa,ni,nj,nk         ! I
!    1                             ,heights_3d                            ! I
     1                             ,vent_2d,istatus)                      ! O

!      real*4 heights_3d(ni,nj,nk)                                        ! I
       real*4 u_3d(ni,nj,nk)                                              ! I
       real*4 v_3d(ni,nj,nk)                                              ! I

       real*4 pbl_top_pa(ni,nj)                                           ! I
       real*4 vent_2d(ni,nj)                                              ! O

       write(6,*)' Subroutine ventilation_index'

!      Find boundary-layer depth in meters

!      Call mean_wind subroutine

!      Multiply to obtain ventilation index

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
