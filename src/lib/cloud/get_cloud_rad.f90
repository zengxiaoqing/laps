

      subroutine get_cloud_rad(sol_alt,sol_azi,clwc_3d,cice_3d,rain_3d, &
            snow_3d,topo_a,heights_3d,transm_3d,transm_4d,idb,jdb,ni,nj,nk)

      use mem_namelist, ONLY: r_missing_data, earth_radius
      use cloud_rad ! Cloud Radiation and Microphysics Parameters
      include 'trigd.inc'

!     airmass(z) = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z))) ! if z <= 91
      airmassn(cosz) =           1.002432 * cosz**2 + 0.148386  * cosz + 0.0096467
      airmassd(cosz) = cosz**3 + 0.149864 * cosz**2 + 0.0102963 * cosz + .000303978
      airmassf(cosz) = airmassn(cosz) / airmassd(cosz)         ! if z <= 93
      trans(od) = exp(-min(od,80.))

      parameter (nc = 3)

      real ext_g(nc),trans_c(nc)  ! od per airmass, transmissivity
      data ext_g /.07,.14,.28/    ! refine via Schaeffer

      real clwc_3d(ni,nj,nk) ! kg/m**3
      real cice_3d(ni,nj,nk) ! kg/m**3
      real rain_3d(ni,nj,nk) ! kg/m**3
      real snow_3d(ni,nj,nk) ! kg/m**3
      real cond_3d(ni,nj,2)
      real heights_3d(ni,nj,nk)
      real transm_3d(ni,nj,nk) ! direct transmission plus forward scattered
      real transm_4d(ni,nj,nk,nc) ! adding 3 color information, account for
                                  ! solar intensity at top of cloud
      real sol_alt(ni,nj)
      real sol_azi(ni,nj)

      real clwc_int(ni,nj)
      real cice_int(ni,nj)
      real rain_int(ni,nj)
      real snow_int(ni,nj)
      real backscatter_int(ni,nj)
      real aod_2d(ni,nj)          ! aerosol optical depth (tau per airmass) 
      real topo_a(ni,nj)

!     n                                    (number concentration:   m**-3)
!     sigma                                (cross-section:          m**2)
!     kappa = n * sigma / rho              (opacity:                m**2 per kg)
!     K or alpha = sigma * n = kappa * rho (extinction coefficient: m**-1)
!     tau = K * s = kappa * rho * s        (optical depth:          dimensionless)

!     Statement functions
!     od_to_albedo_lwc(tau) = 1. - exp(-tau*.07) ! lwc backscattering efficiency term 
!     od_to_albedo_ice(tau) = 1. - exp(-tau*.14) ! ice backscattering efficiency term 
!     od_to_albedo(tau_bks) = 1. - exp(-tau_bks)

      clwc_int = 0.
      cice_int = 0.
      rain_int = 0.
      snow_int = 0.
      backscatter_int = 0.

      call get_grid_spacing_cen(grid_spacing_m,istatus)

!     Note that idb,jdb is a "nominal" grid point location from which we derive a constant solar alt/az for some purposes
      sol_alt_eff = max(sol_alt(idb,jdb),1.5)
      if(sol_alt(idb,jdb) .le. -2.0)sol_alt_eff = +5.0 ! twilight arch light source
      ds_dh = 1. / sind(sol_alt_eff)
      dxy_dh = 1. / tand(sol_alt_eff)

      transm_3d(:,:,:) = 1.

      sinazi = sind(sol_azi(idb,jdb))
      cosazi = cosd(sol_azi(idb,jdb))

      terr_max = maxval(topo_a)

      do k = nk-1,1,-1

        patm_k = ztopsa(heights_3d(idb,jdb,k)) / 1013.

        ku = k+1 ; kl = k

        cond_3d(:,:,1) = clwc_3d(:,:,kl) + cice_3d(:,:,kl) &
                       + rain_3d(:,:,kl) + snow_3d(:,:,kl)

        cond_3d(:,:,2) = clwc_3d(:,:,ku) + cice_3d(:,:,ku) &
                       + rain_3d(:,:,ku) + snow_3d(:,:,ku)

        dh = heights_3d(1,1,ku) - heights_3d(1,1,kl)
        ds = dh * ds_dh
        dij = (dh * dxy_dh) / grid_spacing_m           
        di =  sinazi * dij
        dj =  cosazi * dij
        dil =  sinazi * (heights_3d(1,1,kl) * dxy_dh) / grid_spacing_m
        diu =  sinazi * (heights_3d(1,1,ku) * dxy_dh) / grid_spacing_m
        djl =  cosazi * (heights_3d(1,1,kl) * dxy_dh) / grid_spacing_m
        dju =  cosazi * (heights_3d(1,1,ku) * dxy_dh) / grid_spacing_m

        if(heights_3d(1,1,kl) .lt. terr_max)then
            nsub = max(nint(dij),1)
        else
            nsub = 1
        endif
        write(6,*)
        write(6,*)'k,nsub = ',k,nsub

!       Loop through array of slant columns passing through i,j at sea level
        do i = 1,ni
!        il,ih,jl,jh are offset in reference to MSL height
         il = i + nint(dil) ; il = max(il,1) ; il = min(il,ni)
         iu = i + nint(diu) ; iu = max(iu,1) ; iu = min(iu,ni)

         do j = 1,nj

          jl = j + nint(djl) ; jl = max(jl,1) ; jl = min(jl,nj)
          ju = j + nint(dju) ; ju = max(ju,1) ; ju = min(ju,nj)

          if(il .eq. idb .AND. jl .eq. jdb)then
              idebug = 1
          else
              idebug = 0
          endif

          cond_m = 0.5 * (cond_3d(iu,ju,2 ) + cond_3d(il,jl, 1))   

!         Check slant ray & terrain within grid box coming toward MSL observer
          ihit_terrain = 0
          do isub = nsub,1,-1 
            fracs = float(isub) / float(nsub)
            is = float(iu) * fracs + float(il) * (1. - fracs)
            js = float(ju) * fracs + float(jl) * (1. - fracs)
            height_int = heights_3d(is,js,ku) * fracs + heights_3d(is,js,kl) * (1. - fracs)
            if(height_int .lt. topo_a(is,js))then
              ihit_terrain = ihit_terrain + 1
            endif
            if(idebug .eq. 1)then
              rks = float(kl) + fracs
              write(6,91)is,js,fracs,rks,height_int,topo_a(is,js),ihit_terrain
91            format(' Compare height with terrain: ',2i6,2f8.1,2f10.2,i4)
            endif
          enddo ! isub

          if(ihit_terrain .ge. 1)then
            transm_3d(il,jl,kl) = 0. ! terrain shadow
            if(idebug .eq. 1)then
                ht_agl = heights_3d(il,jl,kl)-topo_a(il,jl)
                write(6,92)il,jl,kl,ht_agl
92              format(' Set terrain shadow for ',2i5,i4,f8.1,'M AGL')
            endif

!         elseif(clwc_m+cice_m+rain_m+snow_m .gt. 0.)then ! for efficiency
          elseif(cond_m .gt. 0.)then ! for efficiency

!           Calculate cloud integrated species along slant path
            clwc_m = 0.5 * (clwc_3d(iu,ju,ku) + clwc_3d(il,jl,kl))
            cice_m = 0.5 * (cice_3d(iu,ju,ku) + cice_3d(il,jl,kl))
            rain_m = 0.5 * (rain_3d(iu,ju,ku) + rain_3d(il,jl,kl))
            snow_m = 0.5 * (snow_3d(iu,ju,ku) + snow_3d(il,jl,kl))

            clwc_lyr_int = clwc_m * ds ! LWP                             
            cice_lyr_int = cice_m * ds                                   
            rain_lyr_int = rain_m * ds                                   
            snow_lyr_int = snow_m * ds                                   
 
!           Convert to optical depth via particle size (Stephens' or LAPS thin cloud equation)
!           od = 3./2. * LWP / reff
            od_lyr_clwc = (1.5 * (clwc_lyr_int / rholiq )) / reff_clwc
            od_lyr_cice = (1.5 * (cice_lyr_int / rholiq )) / reff_cice
            od_lyr_rain = (1.5 * (rain_lyr_int / rholiq )) / reff_rain
            od_lyr_snow = (1.5 * (snow_lyr_int / rhosnow)) / reff_snow

!           Convert to reflectance (using od_to_albedo)
!           od_int = od_int + (od_lyr_clwc + od_lyr_cice + od_lyr_rain + od_lyr_snow)
            backscatter_int(i,j) = backscatter_int(i,j) &
                                           + od_lyr_clwc * bksct_eff_clwc &
                                           + od_lyr_cice * bksct_eff_cice &
                                           + od_lyr_rain * bksct_eff_rain &
                                           + od_lyr_snow * bksct_eff_snow  

!           albedo_int = od_to_albedo(backscatter_int(i,j))                                                            
            albedo_int = 1.0 - exp(-backscatter_int(i,j))                                                            

!           Convert to transmittance
            transm_3d(il,jl,kl) = 1. - albedo_int

            if( idebug .eq. 1 .OR. &
               (transm_3d(il,jl,kl) .eq. 0. .and. jl .eq. jdb) )then
              write(6,101)k,il,clwc_m,cice_m,rain_m,snow_m &
                     ,clwc_lyr_int,cice_lyr_int,rain_lyr_int,snow_lyr_int &
                     ,od_lyr_clwc,od_lyr_cice,od_lyr_rain,od_lyr_snow &
                     ,backscatter_int(i,j),transm_3d(il,jl,kl)
101           format('k/il/clwc/lwp/od/bks/transm: ',i3,i4,2x,4f9.6,2x,4f7.4,2x,4f7.4,1x,2f9.5)
            endif

          else ! for efficiency
            transm_3d(il,jl,kl) = transm_3d(iu,ju,ku)
            if( idebug .eq. 1 .OR. &
               (transm_3d(il,jl,kl) .eq. 0. .and. jl .eq. jdb) )then
              write(6,102)k,heights_3d(il,jl,kl),transm_3d(il,jl,kl),dh,ds&
                         ,sol_alt(il,jl),sol_azi(il,jl),di,dj,iu,ju,il,jl,i,j
102           format('k/ht/trans/dhds/solaltaz/dij/u/l/ij ' &
                    ,i5,f8.0,f8.5,1x,2f9.2,1x,2f7.2,1x,2f7.3,3(2x,2i4))
            endif
          endif

          topo = 1500.
          ht_agl = heights_3d(il,jl,k) - topo

!         See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
          if(ht_agl .gt. 0.)then                               
              horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
          else
              horz_dep_d = 0.
          endif

          refraction = 0.5 ! typical value near horizon

          sol_alt_cld = sol_alt(il,jl) + horz_dep_d + refraction

!         Estimate solar extinction/reddening by Rayleigh scattering at this cloud altitude
          if(sol_alt_cld .lt. 0.)then    ! (early) twilight cloud lighting
            twi_int = .1 * 10.**(+sol_alt_cld * 0.4) ! magnitudes per deg
            rint = twi_int
            grn_rat = 1.0 ; blu_rat = 1.0            
          elseif(sol_alt_cld .ge. 0.)then            ! low daylight sun
            if(.true.)then
!             Direct illumination of the cloud is calculated here
!             Indirect illumination is factored in via 'scat_frac'
              am = airmassf(cosd(90. - max(sol_alt(il,jl),-3.0)))
              scat_frac = 0.75
              do ic = 1,nc
                trans_c(ic) = trans(am*ext_g(ic)*patm_k*scat_frac)
              enddo
              rint = trans_c(1)
              grn_rat = trans_c(2) / trans_c(1)
              blu_rat = trans_c(3) / trans_c(1)
            endif
          endif  

          if(idebug .eq. 1)then
              write(6,103)k,sol_alt(il,jl),horz_dep_d,sol_alt_cld,am*patm_k,rint,rint*blu_rat**0.3,rint*blu_rat
103           format('k/salt/hdep/salt_cld/amk/R/G/B',43x,i4,3f6.2,f8.2,2x,3f6.2)                                   
          endif

!         Modify transm array for each of 3 colors depending on solar intensity and color
          transm_4d(il,jl,kl,1) = transm_3d(il,jl,kl) * rint
          transm_4d(il,jl,kl,2) = transm_3d(il,jl,kl) * rint * grn_rat      
          transm_4d(il,jl,kl,3) = transm_3d(il,jl,kl) * rint * blu_rat   

         enddo ! j
        enddo ! i
        I4_elapsed = ishow_timer()
        write(6,*)' Finished loop for k = ',k
      enddo ! k

      I4_elapsed = ishow_timer()

!     Check the presence of terrain shadow grid points
      n_terr_shadow = 0.
      do k = 1,nk
      do i = 1,ni
      do j = 1,nj
          if(heights_3d(i,j,k) .gt. topo_a(i,j) .AND. transm_3d(i,j,k) .eq. 0.)then
              n_terr_shadow = n_terr_shadow + 1
              if(n_terr_shadow .le. 10)then
                  write(6,*)' terrain shadow is at ',i,j,k
              endif
          endif 
      enddo ! j
      enddo ! i
      enddo ! k      

      write(6,*)' Number of points above ground and in shadow is',n_terr_shadow

      write(6,*)' Returning from get_cloud_rad...'

      return
      end

