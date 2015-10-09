

     subroutine get_cloud_rad(obj_alt,obj_azi,solalt,solazi,clwc_3d,cice_3d,rain_3d, &
           snow_3d,topo_a,lat,lon,heights_3d,transm_3d,transm_4d,idb,jdb,ni,nj,nk,twi_alt,sfc_glow)

     use mem_namelist, ONLY: r_missing_data, earth_radius
     use cloud_rad ! Cloud Radiation and Microphysics Parameters
     include 'trigd.inc'

!    airmass(z) = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z))) ! if z <= 91
!    airmassn(cosz) =           1.002432 * cosz**2 + 0.148386  * cosz + 0.0096467
!    airmassd(cosz) = cosz**3 + 0.149864 * cosz**2 + 0.0102963 * cosz + .000303978
!    airmassf(cosz) = airmassn(cosz) / airmassd(cosz)         ! if z <= 93
     include 'rad.inc'

     trans(od) = exp(-min(od,80.))

!    parameter (nc = 3)

!    real ext_g(nc)           ! od per airmass
     real trans_c(nc)         ! transmissivity
!    data ext_g /.07,.14,.28/ ! refine via Schaeffer

     real clwc_3d(ni,nj,nk) ! kg/m**3
     real cice_3d(ni,nj,nk) ! kg/m**3
     real rain_3d(ni,nj,nk) ! kg/m**3
     real snow_3d(ni,nj,nk) ! kg/m**3
     real cond_3d(ni,nj,2)
     real heights_3d(ni,nj,nk)
     real transm_3d(ni,nj,nk) ! direct transmission plus forward scattered
     real transm_4d(ni,nj,nk,nc) ! adding 3 color information, account for
                                  ! solar intensity at top of cloud
     real obj_alt(ni,nj)
     real obj_azi(ni,nj)
     real sfc_glow(ni,nj)        ! surface lighting intensity (nl)                 

     real clwc_int(ni,nj)
     real cice_int(ni,nj)
     real rain_int(ni,nj)
     real snow_int(ni,nj)
     real backscatter_int(ni,nj)
     real aod_2d(ni,nj)          ! aerosol optical depth (tau per airmass) 
     real topo_a(ni,nj)
     real lat(ni,nj)
     real lon(ni,nj)
     real bi_coeff(2,2)
     real eclipse(ni,nj) ; logical l_solar_eclipse /.false./

!    n                                    (number concentration:   m**-3)
!    sigma                                (cross-section:          m**2)
!    kappa = n * sigma / rho              (opacity:                m**2 per kg)
!    K or alpha = sigma * n = kappa * rho (extinction coefficient: m**-1)
!    tau = K * s = kappa * rho * s        (optical depth:          dimensionless)

     eclipse = 1.0 ! Default is no solar eclipse

     clwc_int = 0.
     cice_int = 0.
     rain_int = 0.
     snow_int = 0.
     backscatter_int = 0.

     call get_grid_spacing_cen(grid_spacing_m,istatus)

!    Note that idb,jdb is a "nominal" grid point location from 
!    which we derive a constant object alt/az for some purposes
     obj_alt_eff = max(obj_alt(idb,jdb),1.5)
     if(obj_alt(idb,jdb) .le. -2.0)obj_alt_eff = +5.0 ! twilight arch light source
     ds_dh = 1. / sind(obj_alt_eff)
     dxy_dh = 1. / tand(obj_alt_eff)

!    Initialize
     transm_3d(:,:,:) = 1.
     transm_4d(:,:,:,:) = r_missing_data ! Zero gives slightly better results than
                              ! one. We might try a more explicit
                              ! calculation for those points that aren't
                              ! covered by the slant rays. In that case
                              ! 'r_missing_data' should be initialized
                              ! here.                

     sinazi = sind(obj_azi(idb,jdb))
     cosazi = cosd(obj_azi(idb,jdb))

     terr_max = maxval(topo_a); terr_min = minval(topo_a)

     I4_elapsed = ishow_timer()

     ht_ref = 5000. ! 0.

     do k = nk-1,1,-1

       patm_k = ztopsa(heights_3d(idb,jdb,k)) / 1013.

       ku = k+1 ; kl = k

       dh = heights_3d(1,1,ku) - heights_3d(1,1,kl)
       ds = dh * ds_dh
       dij = (dh * dxy_dh) / grid_spacing_m           
       di =  sinazi * dij
       dj =  cosazi * dij
       dil =  sinazi * ((heights_3d(1,1,kl)-ht_ref) * dxy_dh) &
                                                       / grid_spacing_m
       diu =  sinazi * ((heights_3d(1,1,ku)-ht_ref) * dxy_dh) &
                                                       / grid_spacing_m
       djl =  cosazi * ((heights_3d(1,1,kl)-ht_ref) * dxy_dh) &
                                                       / grid_spacing_m
       dju =  cosazi * ((heights_3d(1,1,ku)-ht_ref) * dxy_dh) &
                                                       / grid_spacing_m

!      Convert hydrometeor concentration to backscatter optical depth
       const_clwc = ((1.5 / rholiq ) / reff_clwc) * bksct_eff_clwc * ds
       const_cice = ((1.5 / rholiq ) / reff_cice) * bksct_eff_cice * ds
       const_rain = ((1.5 / rholiq ) / reff_rain) * bksct_eff_rain * ds
       const_snow = ((1.5 / rhosnow) / reff_snow) * bksct_eff_snow * ds

       if(ku .eq. nk)then
           cond_3d(:,:,2) &
               = clwc_3d(:,:,ku)*const_clwc + cice_3d(:,:,ku)*const_cice &
               + rain_3d(:,:,ku)*const_rain + snow_3d(:,:,ku)*const_snow
       else
           cond_3d(:,:,2) = cond_3d(:,:,1)  
       endif

       cond_3d(:,:,1) &
               = clwc_3d(:,:,kl)*const_clwc + cice_3d(:,:,kl)*const_cice &
               + rain_3d(:,:,kl)*const_rain + snow_3d(:,:,kl)*const_snow

       write(6,*)

!      Compare heights to terrain zone (move inside ij loops?)
!      if(heights_3d(1,1,kl) .le. terr_max .AND. &
!         heights_3d(1,1,kl) .ge. terr_min       )then ! between
!          nsub = max(nint(dij),1)
!      elseif(heights_3d(1,1,kl) .gt. terr_max)then    ! above
!          nsub = 0
!          ihit_terrain_ref = 0
!      else                                            ! below
!          nsub = 0
!          ihit_terrain_ref = 1
!      endif

       write(6,*)'k = ',k        
       write(6,51)dij,heights_3d(1,1,kl),terr_max,terr_min
51     format(' dij,height,terr range',f8.3,3f8.1)

!      Loop through array of slant columns passing through i,j at ht_ref
       do j = 1,nj

!       il,ih,jl,jh are offset in reference to MSL height
        rjl = j + djl ; rjl = max(rjl,1.) ; rjl = min(rjl,float(nj))
                                             jl = nint(rjl)
        rju = j + dju ; rju = max(rju,1.) ; rju = min(rju,float(nj))
                                             ju = nint(rju)

        do i = 1,ni
         ril = i + dil ; ril = max(ril,1.) ; ril = min(ril,float(ni))
                                              il = nint(ril)
         riu = i + diu ; riu = max(riu,1.) ; riu = min(riu,float(ni))
                                              iu = nint(riu)

          if(il .eq. idb .AND. jl .eq. jdb)then
              idebug = 1
          else
              idebug = 0
          endif

!         Compare heights to terrain zone                          
          heights_3d_il_jl_kl = heights_3d(il,jl,kl)
          if(heights_3d_il_jl_kl .gt. terr_max)then      ! above
            nsub = 0
            ihit_terrain_ref = 0
            cond_m = 0.5 * (cond_3d(iu,ju,2 ) + cond_3d(il,jl, 1))   
          elseif(heights_3d(iu,ju,ku) .lt. terr_min)then ! below
            nsub = 0
            ihit_terrain_ref = 1
            cond_m = 0.0                                             
          else                                           ! between                                                  
            nsub = max(nint(dij),1)
            cond_m = 0.5 * (cond_3d(iu,ju,2 ) + cond_3d(il,jl, 1))   
          endif

!         Check slant ray & terrain within grid box coming toward MSL observer
          if(nsub .gt. 0)then
           ihit_terrain = 0
           isub = nsub
           do while (isub .ge. 1 .AND. ihit_terrain .eq. 0)
            fracs = float(isub) / float(nsub)
            ris = riu * fracs + ril * (1. - fracs); is = nint(ris)
            rjs = rju * fracs + rjl * (1. - fracs); js = nint(rjs)
            height_int = heights_3d(is,js,ku) * fracs + heights_3d(is,js,kl) * (1. - fracs)

!           Interpolate to get topography at fractional grid point
            i1 = min(int(ris),ni-1); fi = ris - i1; i2=i1+1
            j1 = min(int(rjs),nj-1); fj = rjs - j1; j2=j1+1

            bi_coeff(1,1) = (1.-fi) * (1.-fj)
            bi_coeff(2,1) = fi      * (1.-fj)
            bi_coeff(1,2) = (1.-fi) *     fj 
            bi_coeff(2,2) = fi      *     fj 

            topo_bilin = sum(bi_coeff(:,:) * topo_a(i1:i2,j1:j2))

!           if(height_int .lt. topo_a(is,js))then
            if(height_int .lt. topo_bilin)then ! use topo bilin 
              ihit_terrain = ihit_terrain + 1
            endif
            if(idebug .eq. 1)then
              rks = float(kl) + fracs
              write(6,91)is,js,fracs,rks,height_int,topo_a(is,js),ihit_terrain
91            format(' Compare height with terrain: ',2i6,2f8.1,2f10.2,i4)
            endif
            isub = isub - 1
           enddo ! isub
          else ! simpler case (above or below terrain zone)
           ihit_terrain = ihit_terrain_ref
          endif

          if(ihit_terrain .ge. 1)then
            transm_3d(il,jl,kl) = 0. ! terrain shadow
            if(idebug .eq. 1)then
                ht_agl = heights_3d_il_jl_kl-topo_a(il,jl)
                write(6,92)il,jl,kl,ht_agl
92              format(' Set terrain shadow for ',2i5,i4,f8.1,'M AGL')
            endif

!         elseif(clwc_m+cice_m+rain_m+snow_m .gt. 0.)then ! for efficiency
          elseif(cond_m .gt. 0.)then ! for efficiency

!           Convert to reflectance (using od_to_albedo)
            backscatter_int(i,j) = backscatter_int(i,j) &
!                                          + od_lyr_clwc * bksct_eff_clwc &
!                                          + od_lyr_cice * bksct_eff_cice &
!                                          + od_lyr_rain * bksct_eff_rain &
!                                          + od_lyr_snow * bksct_eff_snow  
                                           + cond_m

!           albedo_int = od_to_albedo(backscatter_int(i,j))                                                            
!           albedo_int = 1.0 - exp(-backscatter_int(i,j))                                                            

!           New albedo relationship
            albedo_int = backscatter_int(i,j) / (backscatter_int(i,j) + 1.)

!           Convert to transmittance
            transm_3d(il,jl,kl) = 1. - albedo_int

            if( idebug .eq. 1 .OR. &
               (transm_3d(il,jl,kl) .eq. 0. .and. jl .eq. jdb) )then

!             Calculate cloud integrated species along slant path
              clwc_m = 0.5 * (clwc_3d(iu,ju,ku) + clwc_3d(il,jl,kl))
              cice_m = 0.5 * (cice_3d(iu,ju,ku) + cice_3d(il,jl,kl))
              rain_m = 0.5 * (rain_3d(iu,ju,ku) + rain_3d(il,jl,kl))
              snow_m = 0.5 * (snow_3d(iu,ju,ku) + snow_3d(il,jl,kl))

              clwc_lyr_int = clwc_m * ds ! LWP                             
              cice_lyr_int = cice_m * ds                                   
              rain_lyr_int = rain_m * ds                                   
              snow_lyr_int = snow_m * ds                                   
 
!             Convert to optical depth via particle size (Stephens' or LAPS thin cloud equation)
!             od = 3./2. * LWP / reff
              od_lyr_clwc = (1.5 * (clwc_lyr_int / rholiq )) / reff_clwc
              od_lyr_cice = (1.5 * (cice_lyr_int / rholiq )) / reff_cice
              od_lyr_rain = (1.5 * (rain_lyr_int / rholiq )) / reff_rain
              od_lyr_snow = (1.5 * (snow_lyr_int / rhosnow)) / reff_snow

              write(6,101)k,il,clwc_m,cice_m,rain_m,snow_m &
                     ,clwc_lyr_int,cice_lyr_int,rain_lyr_int,snow_lyr_int &
                     ,od_lyr_clwc,od_lyr_cice,od_lyr_rain,od_lyr_snow &
                     ,backscatter_int(i,j),transm_3d(il,jl,kl)
101           format('k/il/clwc/lwp/od/bks/transm: ',i3,i4,2x,4f9.6,2x,4f7.4,2x,4f7.4,1x,2f9.5)
            endif

          else ! for efficiency
            transm_3d(il,jl,kl) = transm_3d(iu,ju,ku)
!           if( idebug .eq. 1 .OR. &
!              (transm_3d(il,jl,kl) .eq. 0. .and. jl .eq. jdb) )then
            if( idebug .eq. 1 )then 
              write(6,102)k,heights_3d_il_jl_kl,transm_3d(il,jl,kl),dh,ds&
                         ,obj_alt(il,jl),obj_azi(il,jl),di,dj,iu,ju,il,jl,i,j
102           format('k/ht/trans/dhds/solaltaz/dij/u/l/ij ' &
                    ,i4,f8.0,f8.5,1x,2f9.2,1x,2f7.2,1x,2f7.3,3(2x,2i4))
            endif
          endif

          if(ihit_terrain .eq. 0)then
            topo = 1500. ! generic topo value (possibly redp_lvl?)
            ht_agl = heights_3d_il_jl_kl - topo

!           See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
!           Use new statement function?
            if(ht_agl .gt. 0.)then                               
              horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
            else
              horz_dep_d = 0.
            endif

            refraction = 0.5 ! typical value near horizon

            obj_alt_cld = obj_alt(il,jl) + horz_dep_d + refraction

!           Estimate solar extinction/reddening by Rayleigh scattering at this cloud altitude
            if(obj_alt_cld .lt. 0.)then    ! (early) twilight cloud lighting
              twi_int = .1 * 10.**(+obj_alt_cld * 0.4) ! magnitudes per deg
              rint = twi_int
              grn_rat = 1.0 ; blu_rat = 1.0            
            elseif(obj_alt_cld .ge. 0.)then            ! low daylight sun
!             Direct illumination of the cloud is calculated here
!             Indirect illumination is factored in via 'scat_frac'
!             am = airmassf(cosd(90. - max(obj_alt(il,jl),-3.0)))
              am = airmassf(90.-obj_alt(il,jl),patm_k)
              scat_frac = 1.00
              do ic = 1,nc
!               trans_c(ic) = trans(am*ext_g(ic)*patm_k*scat_frac)
                trans_c(ic) = trans(am*ext_g(ic)       *scat_frac)
              enddo
              rint = trans_c(1)
              grn_rat = trans_c(2) / trans_c(1)
              blu_rat = trans_c(3) / trans_c(1)
            endif  

            if(idebug .eq. 1)then
              write(6,103)k,obj_alt(il,jl),horz_dep_d,obj_alt_cld,am*patm_k,rint,rint*blu_rat**0.3,rint*blu_rat
103           format('k/salt/hdep/salt_cld/amk/R/G/B',43x,i4,3f6.2,f8.2,2x,3f6.2)                                   
            endif

!           Modify transm array for each of 3 colors depending on solar intensity and color
            transm_4d(il,jl,kl,1) = transm_3d(il,jl,kl) * rint           * eclipse(i,j)
            transm_4d(il,jl,kl,2) = transm_3d(il,jl,kl) * rint * grn_rat * eclipse(i,j)
            transm_4d(il,jl,kl,3) = transm_3d(il,jl,kl) * rint * blu_rat * eclipse(i,j)

          else
            transm_4d(il,jl,kl,:) = 0.
          endif ! ihit_terrain = 0
         enddo ! i
        enddo ! j
        I4_elapsed = ishow_timer()
        write(6,*)' Finished loop for k = ',k
      enddo ! k

      I4_elapsed = ishow_timer()

!     Check the presence of terrain shadow grid points + add sfc glow
      write(6,*)' heights_3d column = ',heights_3d(idb,jdb,:)
      write(6,*)' transm_3d column = ',transm_3d(idb,jdb,:)
      write(6,*)' transm_4d column = ',transm_4d(idb,jdb,:,1)

!     if(solalt .lt. -4.)then ! use red channel for sfc lighting
!         write(6,*)' Range of transm_4d(red channel nl) = ',minval(transm_4d(:,:,:,1)),maxval(transm_4d(:,:,:,1))
!         write(6,*)' Range of transm_4d(grn channel nl) = ',minval(transm_4d(:,:,:,2)),maxval(transm_4d(:,:,:,2))
!     endif

      write(6,*)' Additional filling of missing areas and at night...'

      n_terr_shadow = 0.
      day_int = 3e9 ! nl
      do k = 1,nk
        patm_k = ztopsa(heights_3d(idb,jdb,k)) / 1013.
        do j = 1,nj
        do i = 1,ni
          if(heights_3d(i,j,k) .gt. topo_a(i,j) .AND. transm_3d(i,j,k) .eq. 0.)then
              n_terr_shadow = n_terr_shadow + 1
              if(n_terr_shadow .le. 10)then
                  write(6,*)' terrain shadow is at ',i,j,k
              endif
          endif 
          if(solalt .lt. twi_alt)then ! use red channel for sfc lighting
              transm_4d(i,j,k,1) = sfc_glow(i,j) * 0.3 ! nominal backsct
              if(transm_4d(i,j,k,2) .eq. r_missing_data)then
!                 write(6,*)' WARNING transm_4d green channel missing',i,j,k
                  transm_4d(i,j,k,2) = 0.
              endif
          elseif(transm_4d(i,j,k,1) .eq. r_missing_data)then 
!             Consider filling in missing transm_4d with low sun
!             values assuming transm_3d = 1
              if(.true.)then
                  topo = 1500. ! generic topo value (possibly redp_lvl?)
                  ht_agl = heights_3d(i,j,k) - topo

!                 See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
                  if(ht_agl .gt. 0.)then                               
                    horz_dep_d = sqrt(2.0 * ht_agl / earth_radius) * 180./3.14
                  else
                    horz_dep_d = 0.
                  endif

                  refraction = 0.5 ! typical value near horizon

                  obj_alt_cld = obj_alt(i,j) + horz_dep_d + refraction

!                 Estimate solar extinction/reddening by Rayleigh scattering at this cloud altitude
                  if(obj_alt_cld .lt. 0.)then    ! (early) twilight cloud lighting
                    twi_int = .1 * 10.**(+obj_alt_cld * 0.4) ! magnitudes per deg
                    rint = twi_int
                    grn_rat = 1.0 ; blu_rat = 1.0            
                    transm_3d_generic = 0.5
                    transm_4d(i,j,k,1) = transm_3d_generic * rint * eclipse(i,j)
                    transm_4d(i,j,k,2) = transm_3d_generic * rint * eclipse(i,j)
                    transm_4d(i,j,k,3) = transm_3d_generic * rint * eclipse(i,j)
                  elseif(obj_alt_cld .ge. 0.)then            ! low daylight sun
!                   Direct illumination of the cloud is calculated here
!                   Indirect illumination is factored in via 'scat_frac'
!                   am = airmassf(cosd(90. - max(obj_alt(i,j),-3.0)))
                    am = airmassf(90.-obj_alt(i,j),patm_k)
                    scat_frac = 0.75
                    do ic = 1,nc
!                     trans_c(ic) = trans(am*ext_g(ic)*patm_k*scat_frac)
                      trans_c(ic) = trans(am*ext_g(ic)       *scat_frac)
                    enddo
                    rint = trans_c(1)
                    grn_rat = trans_c(2) / trans_c(1)
                    blu_rat = trans_c(3) / trans_c(1)
                    transm_3d_generic = 0.5
                    transm_4d(i,j,k,1) = transm_3d_generic * rint    * eclipse(i,j)
                    transm_4d(i,j,k,2) = transm_3d_generic * grn_rat * eclipse(i,j)
                    transm_4d(i,j,k,3) = transm_3d_generic * blu_rat * eclipse(i,j)
                  endif  
                   
                  if(i .eq. idb .AND. j .eq. jdb)then
                      write(6,103)k,obj_alt(i,j),horz_dep_d,obj_alt_cld,am*patm_k,rint,rint*blu_rat**0.3,rint*blu_rat
                  endif

              else ! generic fill in value
                  transm_4d(i,j,k,1) = 0.7
                  transm_4d(i,j,k,2) = 0.5
                  transm_4d(i,j,k,3) = 0.15
              endif
          endif
        enddo ! i
        enddo ! j
      enddo ! k      
      write(6,*)' transm_4d column = ',transm_4d(idb,jdb,:,1)

      write(6,*)' Number of points above ground and in shadow is',n_terr_shadow
      if(solalt .lt. -4.)then ! use red channel for sfc lighting
          write(6,*)' Range of transm_4d(red channel nl) = ',minval(transm_4d(:,:,:,1)),maxval(transm_4d(:,:,:,1))
          write(6,*)' Range of transm_4d(grn channel nl) = ',minval(transm_4d(:,:,:,2)),maxval(transm_4d(:,:,:,2))
      endif

      I4_elapsed = ishow_timer()

      write(6,*)' Returning from get_cloud_rad...'

      return
      end

      subroutine get_sfc_glow(ni,nj,grid_spacing_m,lat,lon,sfc_glow,gnd_glow)

!     Simple surface lighting model for Boulder area based on city populations. 
!     We can later use DMSP or VIIRS imagery for more detailed information
!     on a more global basis.

      real lat(ni,nj)
      real lon(ni,nj)
      real sfc_glow(ni,nj)   ! surface lighting intensity of clouds (nl)
      real gnd_glow(ni,nj)   ! ground lighting intensity (nl)      

      parameter (ncities = 25)
      real ctylat(ncities)
      real ctylon(ncities)
      real ctypop(ncities)   ! population

!     http://www.geonames.org/US/CO/largest-cities-in-colorado.html
!     http://download.geonames.org/export/dump/US.zip
!     http://voices.yahoo.com/largest-cities-colorado-2011-5984666.html?cat=16
!     http://en.wikipedia.org/wiki/List_of_cities_and_towns_in_Colorado

!     Latitude, Longitude, and Population of Cities
      ctylat(1)=40.03;  ctylon(1)=-105.25;  ctypop(1) = 100000 ! Boulder
      ctylat(2)=39.76;  ctylon(2)=-104.88;  ctypop(2) =1000000 ! Denver
      ctylat(3)=39.93;  ctylon(3)=-105.16;  ctypop(3) =  12000 ! Superior
      ctylat(4)=39.97;  ctylon(4)=-105.14;  ctypop(4) =  19000 ! Louisville
      ctylat(5)=40.17;  ctylon(5)=-105.10;  ctypop(5) =  88000 ! Longmont
      ctylat(6)=40.00;  ctylon(6)=-105.10;  ctypop(6) =  25000 ! Lafayette 
      ctylat(7)=38.83;  ctylon(7)=-104.82;  ctypop(7) = 416000 ! Co Springs
      ctylat(8)=39.73;  ctylon(8)=-104.83;  ctypop(8) = 325000 ! Aurora
      ctylat(9) =40.58; ctylon(9) =-105.08; ctypop(9) = 144000 ! Ft Collins
      ctylat(10)=39.71; ctylon(10)=-105.08; ctypop(10)= 143000 ! Lakewood
      ctylat(11)=39.87; ctylon(11)=-104.97; ctypop(11)= 119000 ! Thornton
      ctylat(12)=38.25; ctylon(12)=-104.61; ctypop(12)= 106000 ! Pueblo
      ctylat(13)=39.80; ctylon(13)=-105.09; ctypop(13)= 106000 ! Arvada
      ctylat(14)=39.84; ctylon(14)=-105.03; ctypop(14)= 106000 ! Wstminster
      ctylat(15)=39.58; ctylon(15)=-104.88; ctypop(15)= 100000 ! Centennial
      ctylat(16)=39.95; ctylon(16)=-105.05; ctypop(16)=  58000 ! Broomfield
      ctylat(17)=40.048; ctylon(17)=-105.067; ctypop(17)=  19000 ! Erie
      ctylat(18)=40.084; ctylon(18)=-104.937; ctypop(18)=   4000 ! Dacono
      ctylat(19)=40.083; ctylon(19)=-104.811; ctypop(19)=   7600 ! Ft.Lptn
      ctylat(20)=40.112; ctylon(20)=-104.936; ctypop(20)=  11000 ! Firestn
      ctylat(21)=40.102; ctylon(21)=-104.937; ctypop(21)=   9400 ! Frdrck
      ctylat(22)=39.985; ctylon(22)=-104.815; ctypop(22)=  34000 ! Brghtn
      ctylat(23)=39.823; ctylon(23)=-104.921; ctypop(23)=  48000 ! ComrcCty
      ctylat(24)=40.099; ctylon(24)=-105.161; ctypop(24)=   4000 ! Niwot
      ctylat(25)=40.065; ctylon(25)=-105.191; ctypop(25)=   9300 ! Gunbrl

      sfc_glow = 0.
      gnd_glow = 0.

      do icity = 1,ncities
!       if(icity .eq. 1)then ! Boulder
!         glow_city = 5000.  ! at city edge
!         radius_city = 6000.
!       else                 ! Denver
!         glow_city = 6000.  ! at city edge
!         radius_city = 14000.
!       endif
        glow_city = 5500.
        radius_city = (ctypop(icity)/100000.)**0.38 * 6000.
        call latlon_to_rlapsgrid(ctylat(icity),ctylon(icity),lat,lon,ni,nj &
                                ,ricity,rjcity,istatus)

        write(6,11)ctypop(icity),radius_city,glow_city,ricity,rjcity
11      format(' pop/radius/glow/i/j = ',f10.0,f10.3,f8.1,f8.1,f8.1)

        do j = 1,nj
        do i = 1,ni
          dist_city = sqrt( (float(i)-ricity)**2 + (float(j)-rjcity)**2 ) * grid_spacing_m
          radii_city = dist_city / radius_city
          if(radii_city .gt. 1.0)then
              sfc_glow(i,j) = sfc_glow(i,j) + (glow_city / (radii_city**2.5))
!             gnd_glow(i,j) = 0.
          else
              sfc_glow(i,j) = sfc_glow(i,j) +  glow_city                      ! nl
              gnd_glow(i,j) = gnd_glow(i,j) +  glow_city
          endif
!         if(j .eq. nj/2)then
!           write(6,*)' dist/radii/glow = ',dist_city,radii_city,sfc_glow(i,j)
!         endif
        enddo ! i
        enddo ! j
      enddo ! icity
      write(6,*)' range of sfc_glow (nl) is ',minval(sfc_glow),maxval(sfc_glow)
      write(6,*)' range of gnd_glow (nl) is ',minval(gnd_glow),maxval(gnd_glow)
      return
      end
