

      subroutine get_cloud_rad(sol_alt,sol_azi,clwc_3d,cice_3d,rain_3d,snow_3d, &
                               heights_3d,transm_3d,transm_4d,idb,jdb,ni,nj,nk)

      use mem_namelist, ONLY: r_missing_data, earth_radius
      use cloud_rad ! Cloud Radiation and Microphysics Parameters
      include 'trigd.inc'

      parameter (nc = 3)

      real clwc_3d(ni,nj,nk)
      real cice_3d(ni,nj,nk)
      real rain_3d(ni,nj,nk)
      real snow_3d(ni,nj,nk)
      real heights_3d(ni,nj,nk)
      real transm_3d(ni,nj,nk) ! direct transmission plus forward scattered, account for solar intensity at top of cloud
      real transm_4d(nc,ni,nj,nk) ! adding 3 color information

      real sol_alt(ni,nj)
      real sol_azi(ni,nj)

      real clwc_int(ni,nj)
      real cice_int(ni,nj)
      real rain_int(ni,nj)
      real snow_int(ni,nj)
      real backscatter_int(ni,nj)

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
      ds_dh = 1. / sind(sol_alt_eff)
      dxy_dh = 1. / tand(sol_alt_eff)

      transm_3d(:,:,nk) = 1.

      do k = nk-1,1,-1

        ku = k+1 ; kl = k

        dh = heights_3d(1,1,ku) - heights_3d(1,1,kl)
        ds = dh * ds_dh
        dij = (dh * dxy_dh) / grid_spacing_m           
        di =  sind(sol_azi(idb,jdb)) * dij
        dj = -cosd(sol_azi(idb,jdb)) * dij
        dil =  sind(sol_azi(idb,jdb)) * (heights_3d(1,1,kl) * dxy_dh) / grid_spacing_m
        dih =  sind(sol_azi(idb,jdb)) * (heights_3d(1,1,ku) * dxy_dh) / grid_spacing_m
        djl = -cosd(sol_azi(idb,jdb)) * (heights_3d(1,1,kl) * dxy_dh) / grid_spacing_m
        djh = -cosd(sol_azi(idb,jdb)) * (heights_3d(1,1,ku) * dxy_dh) / grid_spacing_m

        do i = 1,ni
        do j = 1,nj

          if(i .eq. idb .AND. j .eq. jdb)then
              idebug = 1
          else
              idebug = 0
          endif

          il = i + nint(dil) ; il = max(il,1) ; il = min(il,ni)
          iu = i + nint(diu) ; iu = max(iu,1) ; iu = min(iu,ni)
          jl = j + nint(djl) ; jl = max(jl,1) ; jl = min(jl,nj)
          ju = j + nint(dju) ; ju = max(ju,1) ; ju = min(ju,nj)

!         Calculate cloud integrated species along slant path
          clwc_m = 0.5 * (clwc_3d(iu,ju,ku) + clwc_3d(il,jl,kl))
          cice_m = 0.5 * (cice_3d(iu,ju,ku) + cice_3d(il,jl,kl))
          rain_m = 0.5 * (rain_3d(iu,ju,ku) + rain_3d(il,jl,kl))
          snow_m = 0.5 * (snow_3d(iu,ju,ku) + snow_3d(il,jl,kl))
 
          if(clwc_m+cice_m+rain_m+snow_m .gt. 0.)then ! for efficiency

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
            backscatter_int(i,j) = backscatter_int(i,j) + od_lyr_clwc * bksct_eff_clwc &
                                                        + od_lyr_cice * bksct_eff_cice &
                                                        + od_lyr_rain * bksct_eff_rain &
                                                        + od_lyr_snow * bksct_eff_snow  

!           albedo_int = od_to_albedo(backscatter_int(i,j))                                                            
            albedo_int = 1.0 - exp(-backscatter_int(i,j))                                                            

!           Convert to transmittance
            transm_3d(i,j,k) = 1. - albedo_int

            if(idebug .eq. 1)then
              write(6,101)k,il,clwc_m,cice_m,rain_m,snow_m &
                       ,clwc_lyr_int,cice_lyr_int,rain_lyr_int,snow_lyr_int &
                       ,od_lyr_clwc,od_lyr_cice,od_lyr_rain,od_lyr_snow &
                       ,backscatter_int(i,j),transm_3d(i,j,k)
101           format('k/il/clwc/lwp/od/bks/transm: ',i3,i4,2x,4f9.6,2x,4f7.4,2x,4f7.4,1x,2f9.5)
            endif

          else ! for efficiency
            transm_3d(i,j,k) = transm_3d(i,j,ku)
            if(idebug .eq. 1)then
              write(6,102)k,transm_3d(i,j,k),dh,ds,sol_alt(idb,jdb),sol_azi(idb,jdb),di,dj
102           format('k/trans/dhds/solaltaz/dij ',i5,f8.5,1x,2f9.3,1x,2f7.2,1x,2f7.3)
            endif
          endif

!         See http://mintaka.sdsu.edu/GF/explain/atmos_refr/dip.html
          if(heights_3d(i,j,k) .gt. 0.)then ! consider correction for topo
              horz_dep_d = sqrt(2.0 * heights_3d(i,j,k) / earth_radius) * 180./3.14
          else
              horz_dep_d = 0.
          endif

          refraction = 0.5 ! typical value near horizon

          sol_alt_cld = sol_alt(i,j) + horz_dep_d + refraction

!         Estimate solar extinction/reddening by Rayleigh scattering at this cloud altitude
          z = 90. - sol_alt_cld
!         airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))
!         extinction = 0.28 * airmass * patm

          ramp_ang = 5.0
          if(sol_alt_cld .le. ramp_ang)then
            ramp_slope = 1. / ramp_ang
            rint    = max( (1.0 - (ramp_ang - sol_alt_cld) * ramp_slope * 1.0),0.0)
            blu_rat = max( (1.0 - (ramp_ang - sol_alt_cld) * ramp_slope * 2.0),0.0)
          else
            rint = 1.
            blu_rat = 1.
          endif  

          if(idebug .eq. 1)then
              write(6,103)k,sol_alt(i,j),horz_dep_d,sol_alt_cld, rint,rint*blu_rat**0.3,rint*blu_rat
103           format('k/salt/hdep/salt_cld/R/G/B',70x,i4,3f9.3,2x,3f6.2)                                   
          endif

!         Modify transm array for each of 3 colors depending on solar intensity and color at cloud (top)
          transm_4d(1,i,j,k) = transm_3d(i,j,k) * rint
          transm_4d(2,i,j,k) = transm_3d(i,j,k) * rint * blu_rat**0.3 
          transm_4d(3,i,j,k) = transm_3d(i,j,k) * rint * blu_rat   

        enddo ! j
        enddo ! i
      enddo ! k

      return
      end

