c
c**************new routine as adapted at FSL**************************

      subroutine fire_fields(ni,nj,nk,temp_3d,td_3d                 ! I
!    1                      ,heights_3d                             ! I
     1                      ,u_3d,v_3d                              ! I
     1                      ,t_sfc_k,p_sfc_pa                       ! I
!    1                      ,rh_sfc,u_sfc,v_sfc                     ! I
     1                      ,r_missing_data,i4time                  ! I
     1                      ,istatus)                               ! O

      character*10  units_2d
      character*125 comment_2d
      character*3 var_2d

      integer MAX_FIELDS
      parameter (MAX_FIELDS = 3)
      real*4 field_array(ni,nj,MAX_FIELDS)
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
      real*4 u_sfc(ni,nj)
      real*4 v_sfc(ni,nj)

      real*4 pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! L
      real*4 pres_3d_pa(ni,nj,nk)                                   ! L
      real*4 haines_mid_2d(ni,nj)                                   ! L
      real*4 vent_2d(ni,nj)                                         ! L
      real*4 fosberg_2d(ni,nj)                                      ! L

      write(6,*)' Subroutine fire_fields (under construction)'

      ext = 'pbl'
      istat_pbl = 1

      write(6,*)' Read in pressure of PBL top'
      var_2d = 'PTP'
      call get_laps_2dgrid(i4time,0,i4time_nearest
     1                    ,ext,var_2d,units_2d,comment_2d,ni,nj
     1                    ,pbl_top_pa,0,istatus)

      if(istatus .ne. 1)then
          write(6,*)' LAPS PBL top not available'
          istat_pbl = 0
      endif

      write(6,*)' Read in PBL depth'
      var_2d = 'PDM'
      call get_laps_2dgrid(i4time,0,i4time_nearest
     1                    ,ext,var_2d,units_2d,comment_2d,ni,nj
     1                    ,pbl_depth_m,0,istatus)

      if(istatus .ne. 1)then
          write(6,*)' LAPS PBL depth not available'
          istat_pbl = 0
      endif

      call get_pres_3d(i4time,ni,nj,nk,pres_3d_pa,istatus)
      if(istatus .ne. 1)then
          write(6,*)' Could not obtain pres_3d'
          return
      endif

      call cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d        ! I
!    1                          ,heights_3d                         ! I
     1                          ,u_3d,v_3d                          ! I
     1                          ,pbl_top_pa,pbl_depth_m,istat_pbl   ! I
     1                          ,t_sfc_k                            ! I
     1                          ,rh_sfc                             ! I
     1                          ,u_sfc                              ! I
     1                          ,v_sfc                              ! I
     1                          ,p_sfc_pa                           ! I
     1                          ,fosberg_2d                         ! O
     1                          ,haines_mid_2d                      ! O
     1                          ,vent_2d                            ! O
     1                          ,istatus)                           ! O

      if(istatus .eq. 1)then ! write out LFR fire fields file
          write(6,*)' Write out LFR file'

          call move(vent_2d      ,field_array(1,1,1),ni,nj)
          call move(haines_mid_2d,field_array(1,1,2),ni,nj)
          call move(fosberg_2d   ,field_array(1,1,3),ni,nj)

          ext = 'lfr'
          var_a(1) = 'VNT'
          var_a(2) = 'HAM'
          var_a(3) = 'FWI'
          units_a(1) = 'M**2/S'
          units_a(2) = '      '
          units_a(3) = '      '
          comment_a(1) = 'Ventilation Index'
          comment_a(2) = 'Mid-Level Haines Index'
          comment_a(3) = 'Fosberg Fire Wx Index'
          call put_laps_multi_2d(i4time,ext,var_a,units_a
     1                          ,comment_a,field_array,ni,nj
     1                          ,3,istatus)

      else
          write(6,*)' Skipping write of LFR file'

      endif

      I4_elapsed = ishow_timer()

      return
      end

       subroutine cpt_fire_fields(ni,nj,nk,pres_3d_pa,temp_3d,td_3d  ! I
!    1                           ,heights_3d                         ! I
     1                           ,u_3d,v_3d                          ! I
     1                           ,pbl_top_pa,pbl_depth_m,istat_pbl   ! I
     1                           ,t_sfc_k                            ! I
     1                           ,rh_sfc                             ! I
     1                           ,u_sfc                              ! I
     1                           ,v_sfc                              ! I
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

       real*4 pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! I
       real*4 pbl_top_m(ni,nj)                                       ! L

       real*4 t_sfc_k(ni,nj)
       real*4 rh_sfc(ni,nj)
       real*4 u_sfc(ni,nj)
       real*4 v_sfc(ni,nj)
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
!      call fireweatherindex(t_sfc_k,rh_sfc,p_sfc_mb,u_sfc,v_sfc
!    1                      ,fosberg_2d)      

!      Calculate Ventilation Index
       if(istat_pbl .eq. 1)then
           write(6,*)' Calculate Ventilation Index'
           call ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m        ! I
     1                           ,ni,nj,nk                                ! I
!    1                           ,heights_3d                              ! I
     1                           ,vent_2d,istatus)                        ! O
       else
           write(6,*)' Skip ventilation index due to no PBL'

       endif

       return
       end

       subroutine ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m      ! I
     1                             ,ni,nj,nk                              ! I
!    1                             ,heights_3d                            ! I
     1                             ,vent_2d,istatus)                      ! O

!      real*4 heights_3d(ni,nj,nk)                                        ! I
       real*4 u_3d(ni,nj,nk)                                              ! I
       real*4 v_3d(ni,nj,nk)                                              ! I

       real*4 pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                        ! I
       real*4 umean(ni,nj),vmean(ni,nj)                                   ! L
       real*4 vent_2d(ni,nj)                                              ! O

       real*4 topo(ni,nj)           ! Switch to sfc_pres_pa?

       write(6,*)' Subroutine ventilation_index'

!      Calculate mean wind within the PBL
       call pbl_mean_wind(u_3d,v_3d,topo,pbl_top_pa,ni,nj,nk              ! I
     1                   ,umean,vmean,istatus)                            ! O

!      Multiply PBL depth by mean wind to obtain ventilation index
       do i = 1,ni
       do j = 1,nj
           spmean = sqrt(umean(i,j)**2 + vmean(i,j)**2)
           vent_2d(i,j) = pbl_depth_m(i,j) * spmean
       enddo ! j
       enddo ! i

       return
       end



        subroutine pbl_mean_wind(uanl,vanl,topo,pbl_top_pa        ! I
     1                          ,imax,jmax,kmax                   ! I
     1                          ,umean,vmean,istatus)             ! O

        logical ltest_vertical_grid

        real*4 umean(imax,jmax),vmean(imax,jmax)                  ! Output
        real*4 uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! Input

        real*4 topo(imax,jmax)                                    ! Input
        real*4 pbl_top_pa(imax,jmax)                              ! Input

        real*4 sum(imax,jmax)                                     ! Local
        real*4 usum(imax,jmax)                                    ! Local
        real*4 vsum(imax,jmax)                                    ! Local
        integer*4 klow(imax,jmax)                                 ! Local

        topo = 0.             ! Just for testing (also switch to pres_sfc_pa)?

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        do j = 1,jmax
          do i = 1,imax

             khigh = nint(zcoord_of_pressure(pbl_top_pa(i,j)))

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
