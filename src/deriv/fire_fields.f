c
c**************new routine as adapted at FSL**************************

      subroutine fire_fields(ni,nj,nk,temp_3d,td_3d                 ! I
!    1                      ,heights_3d                             ! I
     1                      ,u_3d,v_3d                              ! I
     1                      ,t_sfc_k,p_sfc_pa                       ! I
     1                      ,rh_sfc,u_sfc,v_sfc                     ! I
     1                      ,r_missing_data,i4time                  ! I
     1                      ,istatus)                               ! O

      character*10  units_2d
      character*125 comment_2d
      character*3 var_2d

      integer MAX_FIELDS
      parameter (MAX_FIELDS = 7)
      real field_array(ni,nj,MAX_FIELDS)
      character*3 var_a(MAX_FIELDS)
      character*125 comment_a(MAX_FIELDS)
      character*10  units_a(MAX_FIELDS)
      character*31 ext

!     real heights_3d(ni,nj,nk)                                   ! I
      real temp_3d(ni,nj,nk)                                      ! I
      real td_3d(ni,nj,nk)                                        ! I
      real u_3d(ni,nj,nk)
      real v_3d(ni,nj,nk)

      real t_sfc_k(ni,nj)                                         ! I 
      real rh_sfc(ni,nj)                                          ! I 
      real p_sfc_pa(ni,nj)                                        ! I 
      real u_sfc(ni,nj)                                           ! I
      real v_sfc(ni,nj)                                           ! I

      real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! L
      real pres_3d_pa(ni,nj,nk)                                   ! L
      real haines_mid_2d(ni,nj)                                   ! L
      real haines_hi_2d(ni,nj)                                    ! L
      real vent_2d(ni,nj)                                         ! L
      real fosberg_2d(ni,nj)                                      ! L
      real umean_2d(ni,nj),vmean_2d(ni,nj)                        ! L
      real cfwi(ni,nj)                                            ! L

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
     1                          ,r_missing_data                     ! I
     1                          ,i4time                             ! I
     1                          ,fosberg_2d                         ! O
     1                          ,haines_mid_2d                      ! O
     1                          ,haines_hi_2d                       ! O
     1                          ,vent_2d                            ! O
     1                          ,umean_2d,vmean_2d                  ! O
     1                          ,cfwi                               ! O
     1                          ,istatus)                           ! O

      if(istatus .eq. 1)then ! write out LFR fire fields file
          write(6,*)' Write out LFR file'

          call move(vent_2d      ,field_array(1,1,1),ni,nj)
          call move(haines_mid_2d,field_array(1,1,2),ni,nj)
          call move(haines_hi_2d, field_array(1,1,3),ni,nj)
          call move(fosberg_2d   ,field_array(1,1,4),ni,nj)
          call move(umean_2d     ,field_array(1,1,5),ni,nj)
          call move(vmean_2d     ,field_array(1,1,6),ni,nj)
          call move(cfwi         ,field_array(1,1,7),ni,nj)

          ext = 'lfr'
          var_a(1) = 'VNT'
          var_a(2) = 'HAM'
          var_a(3) = 'HAH'
          var_a(4) = 'FWI'
          var_a(5) = 'UPB'
          var_a(6) = 'VPB'
          var_a(7) = 'CWI'
          units_a(1) = 'M**2/S'
          units_a(2) = '      '
          units_a(3) = '      '
          units_a(4) = '      '
          units_a(5) = 'M/S'
          units_a(6) = 'M/S'
          units_a(7) = '      '
          comment_a(1) = 'Ventilation Index'
          comment_a(2) = 'Haines Index (850-700hPa)'
          comment_a(3) = 'Haines Index (700-500hPa)'
          comment_a(4) = 'Fosberg Fire Wx Index'
          comment_a(5) = 'Boundary Layer Mean U Component'
          comment_a(6) = 'Boundary Layer Mean V Component'
          comment_a(7) = 'Critical Fire Weather Index'
          call put_laps_multi_2d(i4time,ext,var_a,units_a
     1                          ,comment_a,field_array,ni,nj
     1                          ,7,istatus)

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
     1                           ,r_missing_data                     ! I
     1                           ,i4time                             ! I
     1                           ,fosberg_2d                         ! O
     1                           ,haines_mid_2d                      ! O
     1                           ,haines_hi_2d                       ! O
     1                           ,vent_2d                            ! O
     1                           ,umean_2d,vmean_2d                  ! O
     1                           ,cfwi                               ! O
     1                           ,istatus)                           ! O

       real pres_3d_pa(ni,nj,nk)                                   ! I
!      real heights_3d(ni,nj,nk)                                   ! I
       real temp_3d(ni,nj,nk)                                      ! I
       real td_3d(ni,nj,nk)                                        ! I
       real u_3d(ni,nj,nk)
       real v_3d(ni,nj,nk)

       real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                   ! I
       real pbl_top_m(ni,nj)                                       ! L

       real t_sfc_k(ni,nj)
       real rh_sfc(ni,nj)
       real u_sfc(ni,nj)
       real v_sfc(ni,nj)
       real p_sfc_pa(ni,nj)                                        ! I
       real p_sfc_mb(ni,nj)

       real haines_mid_2d(ni,nj)                                   ! O
       real haines_hi_2d(ni,nj)                                    ! O
       real vent_2d(ni,nj)                                         ! O
       real umean_2d(ni,nj),vmean_2d(ni,nj)                        ! O
       real fosberg_2d(ni,nj)                                      ! O
       real cfwi(ni,nj)                                            ! O

       real pres_3d_mb(ni,nj,nk)                                   ! L

       write(6,*)' Subroutine cpt_fire_fields...'

!      Calculate Haines Index
       write(6,*)' Calculate Mid-Level Haines Index'
       pres_3d_mb = pres_3d_pa / 100.
       call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_mid_2d,ni,nj,nk       
     1                 ,850.,700.) 
       do i = 1,ni
       do j = 1,nj
           if(p_sfc_pa(i,j) .lt. 85000.)then
               haines_mid_2d(i,j) = r_missing_data
           endif
       enddo ! j
       enddo ! i

       call hainesindex(pres_3d_mb,temp_3d,td_3d,haines_hi_2d,ni,nj,nk
     1                 ,700.,500.) 
       do i = 1,ni
       do j = 1,nj
           if(p_sfc_pa(i,j) .lt. 70000.)then
               haines_hi_2d(i,j) = r_missing_data
           endif
       enddo ! j
       enddo ! i

!      Calculate Fosberg Fireweather Index
       write(6,*)' Calculate Fosberg Fireweather Index'
       p_sfc_mb = p_sfc_pa / 100.
!      Argument list has been changed to use sfc inputs instead of 3D
       call fireweatherindex(t_sfc_k,rh_sfc,p_sfc_mb,u_sfc,v_sfc          ! I
     1                      ,ni,nj                                        ! I
     1                      ,fosberg_2d)                                  ! O

!      Calculate Ventilation Index
       if(istat_pbl .eq. 1)then
           write(6,*)' Calculate Ventilation Index'
           call ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m        ! I
     1                           ,pres_3d_pa,p_sfc_pa                     ! I
     1                           ,ni,nj,nk                                ! I
     1                           ,r_missing_data                          ! I
!    1                           ,heights_3d                              ! I
     1                           ,umean_2d,vmean_2d                       ! O
     1                           ,vent_2d,istatus)                        ! O
       else
           write(6,*)' Skip ventilation index due to no PBL'

       endif

!      Calculate Critical Fire Weather Index.
!         (RH<15% and speed>20mph for any 3 consecutive hours during the
!          past 24 hours)
       write(6,*)' Calculate Critical Fireweather Index',i4time
       call critical_fwi(rh_sfc,u_sfc,v_sfc                               ! I
     1                  ,ni,nj,i4time                                     ! I
     1                  ,cfwi)                                            ! O

       return
       end

       subroutine ventilation_index(u_3d,v_3d,pbl_top_pa,pbl_depth_m      ! I
     1                             ,pres_3d_pa,p_sfc_pa                   ! I
     1                             ,ni,nj,nk                              ! I
     1                             ,r_missing_data                        ! I
!    1                             ,heights_3d                            ! I
     1                             ,umean_2d,vmean_2d                     ! O
     1                             ,vent_2d,istatus)                      ! O

!      real heights_3d(ni,nj,nk)                                        ! I
       real pres_3d_pa(ni,nj,nk)                                        ! I
       real u_3d(ni,nj,nk)                                              ! I
       real v_3d(ni,nj,nk)                                              ! I

       real pbl_top_pa(ni,nj),pbl_depth_m(ni,nj)                        ! I
       real umean_2d(ni,nj),vmean_2d(ni,nj)                             ! O
       real vent_2d(ni,nj)                                              ! O

       real topo(ni,nj)           ! Switch to sfc_pres_pa?
       real p_sfc_pa(ni,nj)                                             ! I

       write(6,*)' Subroutine ventilation_index'

       vent_2d = r_missing_data

!      Calculate mean wind within the PBL
       call pbl_mean_wind(u_3d,v_3d,topo,pbl_top_pa,ni,nj,nk              ! I
     1                   ,pres_3d_pa,p_sfc_pa                             ! I
     1                   ,umean_2d,vmean_2d,istatus)                      ! O
       if(istatus .ne. 1)then
           write(6,*)' WARNING: Bad status returned from pbl_mean_wind'       
           write(6,*)' Returning vent_2d field as missing data'
           return
       endif

       write(6,*)' Compute VI from mean wind speed and PBL Depth'

!      Multiply PBL depth by mean wind to obtain ventilation index
       do i = 1,ni
       do j = 1,nj
           if(abs(umean_2d(i,j)) .gt. 1000.)then
               write(6,*)' ERROR, umean out of bounds',i,j,umean_2d(i,j)
               istatus = 0
               return
           endif

           if(abs(vmean_2d(i,j)) .gt. 1000.)then
               write(6,*)' ERROR, vmean out of bounds',i,j,vmean_2d(i,j)       
               istatus = 0
               return
           endif

           spmean = sqrt(umean_2d(i,j)**2 + vmean_2d(i,j)**2)
           vent_2d(i,j) = pbl_depth_m(i,j) * spmean

       enddo ! j
       enddo ! i

       return
       end



        subroutine pbl_mean_wind(uanl,vanl,topo,pbl_top_pa        ! I
     1                          ,imax,jmax,kmax                   ! I
     1                          ,pres_3d_pa,p_sfc_pa              ! I
     1                          ,umean_2d,vmean_2d,istatus)       ! O

        logical ltest_vertical_grid

        real umean_2d(imax,jmax),vmean_2d(imax,jmax)            ! Output
        real uanl(imax,jmax,kmax),vanl(imax,jmax,kmax)          ! Input
        real pres_3d_pa(imax,jmax,kmax)                         ! Input

        real topo(imax,jmax)                                    ! Input
        real p_sfc_pa(imax,jmax)                                ! Input
        real pbl_top_pa(imax,jmax)                              ! Input

        real sum(imax,jmax)                                     ! Local
        real usum(imax,jmax)                                    ! Local
        real vsum(imax,jmax)                                    ! Local
        integer klow(imax,jmax)                                 ! Local
        integer khigh(imax,jmax)                                ! Local

!       topo = 0.             ! Just for testing (also switch to pres_sfc_pa)?

!       Calculate the mass weighted mean wind for the PBL. Inputs are the
!       lower and upper bounds in terms of pressure. Fractional levels are 
!       not accounted for - only whole levels between the pressure bounds
!       are integrated.

        call get_r_missing_data(r_missing_data,istatus)
        if(istatus .ne. 1)return

        umean_2d = r_missing_data
        vmean_2d = r_missing_data

        do j = 1,jmax
          do i = 1,imax
             if(pbl_top_pa(i,j) .gt. p_sfc_pa(i,j))then
                 write(6,*)' ERROR in pbl_mean_wind: Pbl Top > Sfc P'
     1                    ,i,j,pbl_top_pa(i,j),p_sfc_pa(i,j)
                 istatus = 0
                 return
             endif

             khigh(i,j) = rlevel_of_field(pbl_top_pa(i,j),pres_3d_pa
     1                              ,imax,jmax,kmax,i,j,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: ERROR in rlevel_of_field'
                 return
             endif

             klow(i,j) = rlevel_of_field(p_sfc_pa(i,j),pres_3d_pa
     1                                  ,imax,jmax,kmax,i,j,istatus)
             if(istatus .ne. 1)then
                 write(6,*)' mean_wind: ERROR in rlevel_of_field'
                 return
             endif

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
             enddo
          enddo
        enddo

        do j = 1,jmax
          do i = 1,imax
             if(sum(i,j) .gt. 0.)then ! Mean wind through the layer
                 umean_2d(i,j) = usum(i,j) / sum(i,j)
                 vmean_2d(i,j) = vsum(i,j) / sum(i,j)

             else 
                 write(6,*)' WARNING: sum <= 0',i,j,klow(i,j),khigh(i,j)       

             endif
          enddo ! i
        enddo ! j

        return
        end


       subroutine critical_fwi(rh_sfc,u_sfc,v_sfc                         ! I
     1                        ,ni,nj,i4time                               ! I
     1                        ,cfwi)                                      ! O

       implicit none

       real spd_thresh
       parameter (spd_thresh = 79.90372)   ! 20mph squared in m/s
 
       integer ni,nj,i4time,pi4time,i,j,n,istatus

       real rh_sfc(ni,nj),u_sfc(ni,nj),v_sfc(ni,nj)
     .       ,rh(ni,nj),u(ni,nj),v(ni,nj)
     .       ,cfwi1(ni,nj,24),cfwi(ni,nj)
     .       ,speed

       character*31  ext
       character*10  units_2d
       character*125 comment_2d
       character*3   var_2d

       do n=1,24
          if (n .eq. 1) then

             do j=1,nj
             do i=1,ni
                u(i,j)=u_sfc(i,j)
                v(i,j)=v_sfc(i,j)
                rh(i,j)=rh_sfc(i,j)
             enddo
             enddo

          else

             pi4time = i4time - (n-1)*3600
             print *,'critical_fwi: i4time, pi4time',i4time,pi4time

             ext = 'lsx'

             var_2d = 'U'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,u,istatus)
             var_2d = 'V'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,v,istatus)
             var_2d = 'RH'
             call get_laps_2d(pi4time,ext,var_2d,units_2d,comment_2d
     1                       ,ni,nj,rh,istatus)

          endif

          do j=1,nj
          do i=1,ni
             if (rh(i,j) .gt. 0 .and. rh(i,j) .lt. 15. .and.
     1           u(i,j) .lt. 100. and. v(i,j) .lt. 100.) then
                speed = u(i,j)**2 + v(i,j)**2
                if (speed .gt. spd_thresh) then
                   cfwi1(i,j,n) = 1.
                else
                   cfwi1(i,j,n) = 0.
                endif
             else
                cfwi1(i,j,n) = 0.
             endif
          enddo
          enddo

       enddo

       do j=1,nj
       do i=1,ni
          cfwi(i,j) = 0.
       enddo
       enddo

       do n=3,24
       do j=1,nj
       do i=1,ni
          if (cfwi1(i,j,n)+cfwi1(i,j,n-1)+cfwi1(i,j,n-2) .gt. 0.) 
     1       cfwi(i,j) = 1.
       enddo
       enddo
       enddo

       return
       end
