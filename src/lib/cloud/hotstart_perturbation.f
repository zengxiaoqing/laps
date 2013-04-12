

        subroutine hotstart_perturbation(temp_3d                       ! I/O
     1                                  ,pres_3d                       ! I
     1                                  ,pres_sfc_2d                   ! I
     1                                  ,ht_3d                         ! I
     1                                  ,sh_3d                         ! I
     1                                  ,cld_base_2d                   ! I
     1                                  ,ref_3d                        ! I
     1                                  ,cld_liq_3d                    ! I
     1                                  ,cld_ice_3d                    ! I
     1                                  ,ni,nj,nk                      ! I
     1                                  ,rh_liq_3d                     ! O
     1                                  ,wt_sat_3d                     ! O
     1                                  ,l_lb_3d                       ! O
     1                                  ,l_ub_3d                       ! O
     1                                          )

        use mem_namelist, ONLY: r_missing_data

!       Arrays from argument list
        real temp_3d(ni,nj,nk)            ! Temperature (K)
        real pres_3d(ni,nj,nk)            ! Pressure (Pa)
        real pres_sfc_2d(ni,nj)           ! Surface Pressure (Pa)
        real ht_3d(ni,nj,nk)              ! Height in Meters MSL
        real sh_3d(ni,nj,nk)              ! Specific Humidity (dimensionless)
        real cld_base_2d(ni,nj)           ! Cloud base in Meters MSL
        real ref_3d(ni,nj,nk)             ! Radar Reflectivity (dBZ)
        real cld_liq_3d(ni,nj,nk)         ! Analyzed Cloud Liquid
        real cld_ice_3d(ni,nj,nk)         ! Analyzed Cloud Ice
        real rh_liq_3d(ni,nj,nk)          ! Saturation humidity value (can be less than 1.0 for ice)
        real wt_sat_3d(ni,nj,nk)          ! Weight for use in variational RH constraint
        logical l_lb_3d(ni,nj,nk)         ! Flag to use as lower bound for T
        logical l_ub_3d(ni,nj,nk)         ! Flag to use as upper bound for T

!       Local arrays
        real ref_max_400(ni,nj)

        real k_to_c

        logical l_valid_col

        write(6,*)' Subroutine hotstart_perturbation...'

!       Initialize
        l_lb_3d = .false.                   
        l_ub_3d = .false.                   

!       Determine level above 400mb
        do k = nk,1,-1
            if(pres_3d(ni/2,nj/2,k) .le. 40000.)then
                k400 = k
            endif
        enddo 

        ref_max_400(:,:) = maxval(ref_3d(:,:,k400:nk))

        n_valid_cols = 0

!       Amplitude of temperature perturbation
        do i = 1,ni
        do j = 1,nj

            l_valid_col = .false.

!           Significant reflectivity in the column?
            if(ref_max_400(i,j) .ge. 10. .AND. 
     1         ref_max_400(i,j) .ne. r_missing_data)then

!               Obtain k of cloud base
                k_cld_base = 0
                do k = 2,nk
                    if(ht_3d(i,j,k)   .gt. cld_base_2d(i,j) .AND.
     1                 ht_3d(i,j,k-1) .lt. cld_base_2d(i,j)      )then
                        k_cld_base = k
                    endif
                enddo ! k

!               Obtain k of echo top
                k_echo_top = 0
                do k = 2,nk    
                    if(ref_3d(i,j,k  ) .gt. 10. .AND.                   
     1                 ref_3d(i,j,k-1) .le. 10.       )then
                        k_echo_top = k
                    endif
                enddo ! k

!               Test for valid cloud base and echo top
                if(k_cld_base .gt. 0          .AND. 
     1             k_cld_base .le. k400       .AND.
     1             k_echo_top .ge. k_cld_base .AND.
     1             k_echo_top .ge. k400            ) then
                    l_valid_col = .true.
                    n_valid_cols = n_valid_cols + 1
                    amp_t = (ref_max_400(i,j)-10.) / 10.

!                   Perturb at k400 and below
                    pres_cld_base = pres_3d(i,j,k_cld_base)
                    do k = k_cld_base,k400
                        fracp = (pres_3d(i,j,k) - pres_cld_base)
     1                        / (40000.         - pres_cld_base)
                        fracp2 = sqrt(fracp)
                        pert_t = amp_t * fracp2
                        temp_3d(i,j,k) = temp_3d(i,j,k) + pert_t
                        l_lb_3d(i,j,k) = .true.
                    enddo ! k

!                   Perturb above k400            
                    pres_echo_top = pres_3d(i,j,k_echo_top)
                    if(k_echo_top .gt. k400)then
                      do k = k400+1,k_echo_top
                        fracp = (pres_3d(i,j,k) - pres_echo_top)
     1                        / (40000.         - pres_echo_top)
                        fracp2 = sqrt(fracp)
                        pert_t = amp_t * fracp2
                        temp_3d(i,j,k) = temp_3d(i,j,k) + pert_t
                        l_lb_3d(i,j,k) = .true.
                      enddo ! k
                    endif
   
                endif
            endif

!           Consider wet-bulb perturbation below cloud base
            if(l_valid_col .eqv. .true.)then
                do k = 1,k_cld_base
                    if(pres_3d(i,j,k) .lt. pres_sfc_2d(i,j) .AND.
     1                 ref_3d(i,j,k)  .gt. 10.                   )then

!                       Calculate wet bulb temp
                        t_c = k_to_c(temp_3d(i,j,k))
                        p_mb = pres_3d(i,j,k) / 100.
                        td_c = tdew(p_mb,sh(i,j,k))
                        tw_c = tw(t_c,td_c,p_mb)
                        temp_3d(i,j,k) = c_to_k(tw_c)
                        l_ub_3d(i,j,k) = .true.

                    endif
                enddo ! k
            endif

        enddo ! j
        enddo ! i

        write(6,*)' n_valid_cols = ',n_valid_cols

        return
        end

