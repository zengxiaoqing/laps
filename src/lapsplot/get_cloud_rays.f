         
        subroutine get_cloud_rays(i4time,clwc_3d,cice_3d,heights_3d     ! I
     1                           ,rain_3d,snow_3d                       ! I
     1                           ,pres_3d,aod_3d,topo_sfc,topo_a,swi_2d ! I
     1                           ,topo_albedo_2d                        ! I
     1                           ,topo_swi,topo_albedo                  ! O
     1                           ,aod_ray,aod_2_cloud,aod_2_topo        ! O
     1                           ,aod_ill,aod_ill_dir                   ! O
     1                           ,aod_tot,transm_obs                    ! O
     1                           ,r_cloud_3d,cloud_od,cloud_od_sp       ! O
     1                           ,r_cloud_rad,cloud_rad_c               ! O
     1                           ,clear_rad_c,clear_radf_c,patm         ! O
     1                           ,airmass_2_cloud_3d,airmass_2_topo_3d  ! O
     1                           ,htstart                               ! O
     1                           ,ni,nj,nk,i,j,kstart                   ! I
     1                           ,view_alt,view_az,sol_alt,sol_azi      ! I
     1                           ,alt_norm                              ! I
     1                           ,moon_alt,moon_azi                     ! I
     1                           ,moon_mag,moon_mag_thr                 ! I
     1                           ,l_solar_eclipse,rlat,rlon,lat,lon     ! I
     1                           ,minalt,maxalt,minazi,maxazi           ! I
     1                           ,alt_scale,azi_scale                   ! I
     1                           ,grid_spacing_m,r_missing_data)        ! I

        use mem_namelist, ONLY: earth_radius,aod,aero_scaleht,redp_lvl

        use cloud_rad

        include 'trigd.inc'

!       Statement Functions
        trans(od) = exp(-min(od,80.))
!       brt(a) = (1.0 - exp(-.14 * a))/.14 ! relative sky brightness (max~7)
        brt(a) = (1.0 - exp(-.14 * a))     ! relative sky brightness (max=1)

        angdif(X,Y)=MOD(X-Y+540.,360.)-180.
        angleunitvectors(a1,a2,a3,b1,b2,b3) = acosd(a1*b1+a2*b2+a3*b3)

        parameter (rpd = 3.14159 / 180.)

!       Determine shadow regions of the input 3-D cloud array

        parameter (nc = 3)

        real ext_g(nc),twi_trans_c(nc) ! od per airmass, tramsmissivity
        data ext_g /.07,.14,.28/       ! refine via Schaeffer

        real clwc_3d(ni,nj,nk) ! kg/m**3
        real cice_3d(ni,nj,nk) ! kg/m**3
        real rain_3d(ni,nj,nk) ! kg/m**3
        real snow_3d(ni,nj,nk) ! kg/m**3
        real cond_3d(ni,nj,nk) ! kg/m**3 (effective LWC)
        real aod_3d(ni,nj,nk)  ! aerosol optical depth (per vertical atmosphere)
        real bi_coeff(2,2),tri_coeff(2,2,2)
        real heights_3d(ni,nj,nk) ! MSL
        real transm_3d(ni,nj,nk)
        real transm_4d(ni,nj,nk,nc) ! L
        real transm_4d_m(nc)
        real heights_1d(nk)
        real topo_a(ni,nj)
        real lat(ni,nj)
        real lon(ni,nj)
        real swi_2d(ni,nj)
        real topo_albedo_2d(nc,ni,nj)
        real pres_3d(ni,nj,nk)
        real pres_1d(nk)
        real view_alt(ni,nj)
        real view_az(ni,nj)

        real sol_alt(ni,nj)
        real sol_azi(ni,nj)
        real alt_norm(ni,nj)

        real moon_alt(ni,nj)
        real moon_azi(ni,nj)
        real moon_mag,moon_mag_thr

        real obj_alt(ni,nj)
        real obj_azi(ni,nj)

        real gnd_glow(ni,nj)        ! ground lighting intensity (nl)                 

!       logical l_process(minalt:maxalt,minazi:maxazi)
        logical l_solar_eclipse
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        parameter (nsp = 4)
        real cvr_path_sum_sp(nsp)
        real wt_sp(nsp) ; data wt_sp /1.0,0.5,0.02,0.0714/

        real elong(minalt:maxalt,minazi:maxazi)
        real aod_ray_eff(minalt:maxalt,minazi:maxazi) ! zenithal
        real aod_ray_dir(minalt:maxalt,minazi:maxazi) ! zenithal direct (dummy)
        real r_cloud_3d(minalt:maxalt,minazi:maxazi)  ! cloud opacity
        real cloud_od(minalt:maxalt,minazi:maxazi)    ! cloud optical depth
        real cloud_od_sp(minalt:maxalt,minazi:maxazi,nsp)! cloud species tau
        real r_cloud_rad(minalt:maxalt,minazi:maxazi)    ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,minalt:maxalt,minazi:maxazi) ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi) ! clear sky illumination
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)! integrated 
               ! fraction of air illuminated by the sun along line of sight
               ! (consider Earth's shadow + clouds, used when sun is below
               !  the horizon)
        real ag_2d(minalt:maxalt,minazi:maxazi)       ! dummy
        real airmass_2_cloud_3d(minalt:maxalt,minazi:maxazi)
        real airmass_2_topo_3d(minalt:maxalt,minazi:maxazi)
        real topo_swi(minalt:maxalt,minazi:maxazi)    ! short wave down W/m**2
        real topo_albedo(nc,minalt:maxalt,minazi:maxazi)
        real aod_2_cloud(minalt:maxalt,minazi:maxazi) ! slant path
        real aod_2_topo(minalt:maxalt,minazi:maxazi)  ! slant path
        real aod_ill(minalt:maxalt,minazi:maxazi)     ! slant path
        real aod_ill_dir(minalt:maxalt,minazi:maxazi) ! slant path
        real aod_tot(minalt:maxalt,minazi:maxazi)     ! slant path
        real sum_odrad_c(nc)

        character*1 cslant

        I4_elapsed = ishow_timer()

        write(6,*)' Subroutine get_cloud_rays... ',i,j

!       moon_alt = -10.0
!       moon_azi = 0.

        radius_earth_8_thirds = 6371.e3 * 2.6666666

        pstd = 101325.
        airmass_2_cloud_3d = 0.
        airmass_2_topo_3d = 0.
        topo_swi = 0.
        r_cloud_rad = 1.
        clear_rad_c = 0.
        clear_radf_c = 0.

        idebug_a = 0

        icloud_tot = 0                          

        if(((sol_alt(i,j) .lt. -6.0  .AND. moon_alt(i,j) .gt.  0.0) .OR.
     1      (sol_alt(i,j) .lt. -16.0 .AND. moon_alt(i,j) .gt. -6.0))
     1                         .AND. 
     1                moon_mag .le. moon_mag_thr                  )then
            obj_alt = moon_alt
            obj_azi = moon_azi
            obj_bri = 10. ** ((-26.7 - moon_mag)*0.4)
            write(6,*)' Object is moon, brightness is:',obj_bri
!       else ! add case for twilight light source (depending on handling
!                                                  in get_cloud_rad)
        else
            obj_alt = sol_alt
            obj_azi = sol_azi
            obj_bri = 1.0
            write(6,*)' Object is sun, brightness is:',obj_bri
        endif

        cloud_rad_c = obj_bri

        write(6,*)' range of clwc,cice is ',maxval(clwc_3d)
     1                                     ,maxval(cice_3d)

        write(6,*)' range of heights is ',minval(heights_3d)
     1                                   ,maxval(heights_3d)

        I4_elapsed = ishow_timer()

        write(6,*)' call get_cloud_rad...'
      
        if(.true.)then
            call get_cloud_rad(obj_alt,obj_azi,sol_alt(i,j),sol_azi(i,j)
     1                    ,clwc_3d,cice_3d
     1                    ,rain_3d,snow_3d,topo_a,lat,lon
     1                    ,heights_3d,transm_3d,transm_4d,i,j,ni,nj,nk
     1                    ,gnd_glow)
!           do jj = 1,nj
!           do ii = 1,ni
!               if(gnd_glow(ii,jj) .gt. 0.)then
!                   write(6,*)' ii,jj,gndglow',ii,jj,gnd_glow(ii,jj)
!               endif
!           enddo ! i
!           enddo ! j
        endif

        if(sol_alt(i,j) .lt. -4.)then ! Modify moon glow section in green channel
                                      ! Preserve sfc sky glow section in red channel
            transm_4d(:,:,:,2) = transm_4d(:,:,:,2) * obj_bri ! correct for sun/moon brightness
        else
            transm_4d(:,:,:,:) = transm_4d(:,:,:,:) * obj_bri ! correct for sun/moon brightness
        endif
        write(6,*)' Range of transm_4d(red channel) = '
     1           ,minval(transm_4d(:,:,:,1)),maxval(transm_4d(:,:,:,1))

        I4_elapsed = ishow_timer()

        topo_max_ang = 45.0
        topo_max_ht = maxval(topo_a)

!       Relative values using constants from 'module_cloud_rad.f90'
!       Values are 1.5 / (rho * reff)
!       clwc2alpha = 75. 
        clwc2alpha = 1.5 / (rholiq  * reff_clwc)
        cice2alpha = 1.5 / (rholiq  * reff_cice)
        rain2alpha = 1.5 / (rholiq  * reff_rain)
        snow2alpha = 1.5 / (rhosnow * reff_snow)

        cond_3d = clwc_3d 
     1          + cice_3d * (cice2alpha/clwc2alpha)
     1          + rain_3d * (rain2alpha/clwc2alpha) 
     1          + snow_3d * (snow2alpha/clwc2alpha)

        ri = i
        rj = j

!       kstart = 0 ! 0 means sfc, otherwise level of start

        if(kstart .eq. 0)then ! start from sfc
          htstart = topo_sfc  ! MSL
          write(6,*)' i/j/topo_sfc = ',i,j,topo_sfc
          write(6,*)' height column = ',heights_3d(i,j,:)
          do k = 1,nk-1
            if(heights_3d(i,j,k)   .le. topo_sfc .AND.
     1         heights_3d(i,j,k+1) .ge. topo_sfc      )then
                frach = (topo_sfc - heights_3d(i,j,k)) / 
     1                  (heights_3d(i,j,k+1) - heights_3d(i,j,k))
                ksfc = k
                write(6,*)' Start at sfc, ksfc = ',ksfc
                rksfc = float(ksfc) + frach
            endif
          enddo
          rkstart = rksfc
          patm = pres_3d(i,j,ksfc) / 101325.
        else ! start aloft
          rkstart = kstart
          htstart = heights_3d(i,j,kstart)
          patm = pres_3d(i,j,kstart) / 101325.
          write(6,*)' Start aloft at k = ',kstart
        endif

        aod_ray = aod * exp(-(htstart-redp_lvl)/aero_scaleht)

        write(6,*)' rkstart/htstart/patm = ',rkstart,htstart,patm
        write(6,*)' aero_scaleht = ',aero_scaleht
        write(6,*)' aod/redp_lvl/aod_ray = ',aod,redp_lvl,aod_ray

        aod_ray_eff = aod_ray
        aod_ray_dir = aod_ray

        heights_1d(:) = heights_3d(i,j,:)
        pres_1d(:)    = pres_3d(i,j,:)

        idelt = nint(2. / alt_scale)

!       azid1 = 46. ; azid2 = 226.
        azid1 = 90. ; azid2 = 270.
        if(sol_alt(i,j) .gt. 0.)then
            azid1 = nint(sol_azi(i,j))
!           azid2 = mod(azid1+180.,360.)
            azid2 = azid1
        elseif(moon_alt(i,j) .gt. 0.)then
            azid1 = nint(moon_azi(i,j))
            azid2 = mod(azid1+180.,360.)
        endif

        do ialt = minalt,maxalt

!        l_process = .false.

         call get_val(ialt,minalt,alt_scale,altray)

         if(altray .ge. 20.)then     ! 2 deg alt, 2 deg azi
             if(ialt .eq. (ialt/idelt)*idelt)then
                 jazi_delt = nint(2. / azi_scale)
             else
                 jazi_delt = maxazi
             endif
         elseif(altray .ge. 10.)then ! alt_scale, 1 deg azi
             jazi_delt = nint(1. / azi_scale)
         else                        ! alt_scale, azi_scale
             jazi_delt = 1
         endif

         do jazi = minazi,maxazi,jazi_delt
          view_azi_deg = float(jazi) * azi_scale

          if((view_azi_deg .eq. azid1 .or. 
     1        view_azi_deg .eq. azid2)            .AND.
     1        (abs(altray) .eq. 12  .or. abs(altray) .eq. 9 .or.
     1         abs(altray) .le.  9. .or. ialt .eq. minalt .or.
     1         abs(altray) .eq. 20. .or. altray .eq. 30. .or.
     1             altray  .eq. 40. .or. altray .eq. 50.)      )then
              idebug = 1
              idebug_a(ialt,jazi) = 1
          else
              idebug = 0
          endif

!         Trace towards sky from each grid point
!         view_altitude_deg = max(altray,0.) ! handle Earth curvature later
          view_altitude_deg = altray

!         Get direction cosines based on azimuth
          xcos = sind(view_azi_deg)
          ycos = cosd(view_azi_deg)

          icloud = 0
          cvr_path_sum = 0.
          cvr_path_sum_sp = 0.
          sum_odrad = 0.
          sum_odrad_c = 0.
          sum_clrrad = 0.
          sum_aod = 0.
          sum_aod_ill = 0.
          sum_aod_ill_dir = 0.
          sum_am2cld_num = 0.
          sum_am2cld_den = 0.

          if(.true.)then
!           do k = ksfc,ksfc
              r_cloud_3d(ialt,jazi) = 0.
              cloud_od(ialt,jazi) = 0.
              cloud_od_sp(ialt,jazi,:) = 0.
              
              ri1 = ri
              rj1 = rj

!             ht_sfc = heights_1d(ksfc)
!             ht_sfc = topo_sfc       

              grid_factor = grid_spacing_m / 3000.

              if(view_altitude_deg .lt. -15.)then
                  rkdelt1 = -1.00 * grid_factor
                  rkdelt2 = -1.00 * grid_factor
!             elseif(view_altitude_deg .lt. -4.5)then ! produces red image at -6
!                 rkdelt1 = -0.50 * grid_factor
!                 rkdelt2 = -0.50 * grid_factor
              elseif(view_altitude_deg .lt. 0.0)then
                  rkdelt1 =  0.0
                  rkdelt2 =  0.0
              elseif(view_altitude_deg .le. 0.5)then
                  rkdelt1 = 0.0 ! 0.01 * grid_factor 
                  rkdelt2 = 0.0 ! 0.10 * grid_factor 
              elseif(view_altitude_deg .le. 4.)then
                  rkdelt1 = 0.0 ! 0.10 * grid_factor
                  rkdelt2 = 0.0 ! 0.25 * grid_factor
              elseif(view_altitude_deg .le. 15.)then
                  rkdelt1 = 0.0 ! 0.50 * grid_factor
                  rkdelt2 = 0.0 ! 0.50 * grid_factor
              elseif(view_altitude_deg .le. 45.)then
                  rkdelt1 = 1.00 * grid_factor
                  rkdelt2 = 1.00 * grid_factor
              else
                  rkdelt1 = 1.0
                  rkdelt2 = 1.0
              endif

!             arg = max(view_altitude_deg,1.0)
!             rkdelt = tand(90. - arg)                     
!             rkdelt = max(min(rkdelt,2.0),0.5)
              
              if(idebug .eq. 1)write(6,*)
              if(idebug .eq. 1)then
                write(6,11)altray,view_azi_deg,ialt,jazi,jazi_delt
     1                    ,rkdelt,i,j
11              format(
     1      'Trace the slant path (alt/azi/ialt/jazi/jdelt/rkdelt/i/j):'
     1                ,2f6.1,3i5,f6.2,2i5)                                     
                write(6,12)                                   
12              format('   dz1_l   dz1_h   dxy1_l   dxy1_h  rin',
     1           'ew  rjnew   rk    ht_m   topo_m  ',
     1           ' path     lwc    ice    rain   snow      slant',
     1           '  cvrpathsum  cloudfrac  airmass cloud_rad',
     1           ' aeroext  transm aod_sum aod_sum_ill')
              endif

!             Initialize ray
              dxy1_h = 0.
              dz1_h = 0.
              slant1_h = 0.
              airmass1_h = 0.
              ihit_topo = 0
              ihit_bounds = 0
              rk_h = rkstart

              rkdelt = rkdelt1

              rk = rkstart          

!             Up front check for going below the domain (and topo)
!             if((rk + rkdelt) .lt. 1.0 .AND. rkdelt .lt. 0.)then
!                 ihit_topo = 1
!             endif

              do while(rk .le. float(nk) - rkdelt 
     1           .AND. ihit_topo .eq. 0
     1           .AND. ihit_bounds .eq. 0
     1           .AND. cvr_path_sum .le. 1.0) ! Tau < ~75

                if(rkdelt .ne. 0.)then ! trace by pressure levels

                  if(.false.)then ! optimize step size
!                 if(view_altitude_deg .gt. 0.)then ! optimize step size
                    rkdelt_vert = 1.0 - (rk - int(rk))
                    kk_ref = min(int(rk),nk-1)
                    delta_grid_height = heights_1d(kk_ref+1)
     1                                - heights_1d(kk_ref)
                    aspect_ratio = delta_grid_height / grid_spacing_m
                    rkdelt_horz = aspect_ratio / tand(view_altitude_deg) 
                    rkdelt = max(min(rkdelt_horz,rkdelt_vert),0.2)
                    if(idebug .eq. 1)then
                      write(6,*)' rk,rkdelt_vert,rkdelt_horz,rkdelt',
     1                            rk,rkdelt_vert,rkdelt_horz,rkdelt
                      write(6,*)' view_altitude_deg,delta_grid_height',      
     1                            view_altitude_deg,delta_grid_height
                      write(6,*)' aspect_ratio',
     1                            aspect_ratio
!                     stop
                    endif
                  endif

                  rk = rk + rkdelt

                  rk_l = rk - rkdelt
                  kk_l = int(rk_l)
                  frac_l = rk_l - float(kk_l)

                  if(kk_l .ge. nk .or. kk_l .le. 0)then
                    write(6,*)' ERROR: rk/kk_l',rk,kk_l
                  endif

                  rk_h = rk
                  kk_h = int(rk_h)
                  frac_h = rk_h - float(kk_h)

                  ht_l = heights_1d(kk_l) * (1. - frac_l) 
     1                 + heights_1d(kk_l+1) * frac_l

                  pr_l = pres_1d(kk_l) * (1. - frac_l) 
     1                 + pres_1d(kk_l+1) * frac_l

                  if(rk_h .lt. float(nk) .AND. rk_h .ge. 1.0)then
                    ht_h = heights_1d(kk_h) * (1. - frac_h) 
     1                   + heights_1d(kk_h+1) * frac_h

                    pr_h = pres_1d(kk_h) * (1. - frac_h) 
     1                   + pres_1d(kk_h+1) * frac_h

                  elseif(rk_h .eq. float(nk))then
                    ht_h = heights_1d(kk_h)

                    pr_h = pres_1d(kk_h)
                  else
                    write(6,*)' ERROR: rk_h is out of bounds ',rk_h,nk
     1                                                        ,rkdelt
                    istatus = 0
                    return
                  endif

                  dz1_l = ht_l - htstart        
                  dz1_h = ht_h - htstart          
                  dz2   = dz1_h - dz1_l ! layer

                  if(.false.)then
                    dxy1_l = dz1_l * tand(90. - view_altitude_deg)
                    dxy1_h = dz1_h * tand(90. - view_altitude_deg)
                    dxy2   = dz2   * tand(90. - view_altitude_deg)
                  else ! horzdist call (determine distance from height & elev)
                    aterm = 1. / radius_earth_8_thirds
                    bterm = tand(min(view_altitude_deg,89.9))

                    cterm = dz1_l                                            
                    discriminant = 4.*aterm*cterm + bterm**2.

                    if(discriminant .gt. 0.)then     
                      if(view_altitude_deg .gt. 0.)then ! above horizon
                          dxy1_l = ( sqrt(discriminant) - bterm)   
     1                                              / (2.*aterm)
                      else                           ! below horizon
                          dxy1_l = (-sqrt(discriminant) - bterm)   
     1                                              / (2.*aterm)
                      endif
                    else                             ! catch start point
                      dxy1_l = 0.
                    endif

                    cterm = dz1_h                                          
                    discriminant = 4.*aterm*cterm + bterm**2.

                    if(view_altitude_deg .gt. 0.)then ! above horizon
                        dxy1_h = ( sqrt(discriminant) - bterm) 
     1                                           / (2.*aterm)
                    else                           ! below horizon
                        dxy1_h = (-sqrt(discriminant) - bterm) 
     1                                           / (2.*aterm)
                    endif
                 
                    dxy2 = dxy1_h - dxy1_l
                  endif

                  slant2 = sqrt(dxy2**2 + dz2**2) ! layer
                  slant1_l = slant1_h
                  slant1_h = slant1_h + slant2

                  airmassv = (pr_l - pr_h) / pstd
                  airmass2 = airmassv * (slant2 / dz2)

                  if(dxy2 .gt. grid_spacing_m)then
                    cslant = '*'
                  else
                    cslant = ' '
                  endif

                else  ! Trace by slant range (rkdelt = 0. )
                  slant2 = grid_spacing_m
                  slant1_l = slant1_h
                  slant1_h = slant1_h + slant2

                  dxy1_l = dxy1_h
                  dxy1_h = slant1_h * cosd(altray)

!                 Determine height from distance and elev
                  dz1_l = dz1_h
                  dz1_h = slant1_h * sind(altray) 
     1                  + dxy1_h**2 / (2. * radius_earth_8_thirds)

                  rk_l = rk_h
                  ht_l = htstart + dz1_l
                  ht_h = htstart + dz1_h

                  rk_h = r_missing_data

                  do k = 1,nk-1
                    if(heights_1d(k)   .le. ht_h .AND.
     1                 heights_1d(k+1) .ge. ht_h      )then
                        frach = (ht_h - heights_1d(k)) / 
     1                          (heights_1d(k+1) - heights_1d(k))      
                        rk_h = float(k) + frach
                        pr_h = (pres_1d(k)*(1.-frach)) 
     1                       + (pres_1d(k+1)*frach)
                    endif
                  enddo

                  rk = rk_h

                  airmass2 = (slant2 / 8000.) * (pr_h / pstd)

                endif ! rkdelt .ne. 0.

                cslant = ' '

                dx1_l = dxy1_l * xcos
                dy1_l = dxy1_l * ycos

                dx1_h = dxy1_h * xcos
                dy1_h = dxy1_h * ycos

                airmass1_l = airmass1_h
                airmass1_h = airmass1_h + airmass2

                ridelt_l = dx1_l / grid_spacing_m
                rjdelt_l = dy1_l / grid_spacing_m

                ridelt_h = dx1_h / grid_spacing_m
                rjdelt_h = dy1_h / grid_spacing_m

                rinew_l = ri + ridelt_l
                rinew_h = ri + ridelt_h

                rjnew_l = rj + rjdelt_l
                rjnew_h = rj + rjdelt_h

                rinew_m = 0.5 * (rinew_l + rinew_h)
                rjnew_m = 0.5 * (rjnew_l + rjnew_h)

                rk_m = 0.5 * (rk_l + rk_h)
                ht_m = 0.5 * (ht_l + ht_h)

                if((rk_h-rk_l) .gt. 0.)then
                    grid_dh = (ht_h-ht_l) / (rk_h-rk_l)
                    grid_tan = grid_dh / grid_spacing_m
                else
                    grid_tan = 0.25 ! nominal value
                endif

                rni = ni
                rnj = nj

                inew_m = nint(rinew_m); jnew_m = nint(rjnew_m)
                k_m = nint(rk_m)

                if(rinew_h .ge. 1. .and. rinew_h .le. rni .AND. 
     1             rjnew_h .ge. 1. .and. rjnew_h .le. rnj .AND.
     1             rk_h    .ne. r_missing_data )then 

!                 Cloud present on any of 8 interpolation vertices
                  if(maxval(
     1              cond_3d(int(rinew_m):int(rinew_m)+1
     1                     ,int(rjnew_m):int(rjnew_m)+1
     1                     ,int(rk_m)   :int(rk_m)+1    ) ) .gt. 0.)then        

                      i1 = min(int(rinew_m),ni-1)
                      j1 = min(int(rjnew_m),nj-1) 
                      k1 = min(int(rk_m)   ,nk-1) 

                      fi = rinew_m - i1; i2=i1+1
                      fj = rjnew_m - j1; j2=j1+1
                      fk = rk_m    - k1; k2=k1+1

                      tri_coeff(1,1,1) = (1.-fi) * (1.-fj) * (1.- fk)
                      tri_coeff(2,1,1) = fi      * (1.-fj) * (1.- fk)
                      tri_coeff(1,2,1) = (1.-fi) *     fj  * (1.- fk)
                      tri_coeff(1,1,2) = (1.-fi) * (1.-fj) *      fk
                      tri_coeff(1,2,2) = (1.-fi) *     fj  *      fk
                      tri_coeff(2,1,2) = fi      * (1.-fj) *      fk
                      tri_coeff(2,2,1) = fi      *     fj  * (1.- fk)
                      tri_coeff(2,2,2) = fi      *     fj  *      fk

                      cond_m = sum(tri_coeff(:,:,:) * 
     1                             cond_3d(i1:i2,j1:j2,k1:k2))

                      clwc_m = sum(tri_coeff(:,:,:) * 
     1                             clwc_3d(i1:i2,j1:j2,k1:k2))

                      cice_m = sum(tri_coeff(:,:,:) * 
     1                             cice_3d(i1:i2,j1:j2,k1:k2))

                      rain_m = sum(tri_coeff(:,:,:) * 
     1                             rain_3d(i1:i2,j1:j2,k1:k2))

                      snow_m = sum(tri_coeff(:,:,:) * 
     1                             snow_3d(i1:i2,j1:j2,k1:k2))

!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                   ,cond_3d,cond_m)
!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                   ,clwc_3d,clwc_m)
!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                   ,cice_3d,cice_m)
!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                   ,rain_3d,rain_m)
!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                   ,snow_3d,snow_m)
                  else
                      cond_m = 0.
                      clwc_m = 0.
                      cice_m = 0.
                      rain_m = 0.
                      snow_m = 0.
                  endif

                  cvr_path = cond_m                                 
                  cvr_path_sum_last = cvr_path_sum
                  cvr_path_sum      = cvr_path_sum + cvr_path * slant2
                  cvr_path_sum_sp(1) = 
     1            cvr_path_sum_sp(1) + clwc_m * wt_sp(1) * slant2
                  cvr_path_sum_sp(2) = 
     1            cvr_path_sum_sp(2) + cice_m * wt_sp(2) * slant2
                  cvr_path_sum_sp(3) = 
     1            cvr_path_sum_sp(3) + rain_m * wt_sp(3) * slant2
                  cvr_path_sum_sp(4) = 
     1            cvr_path_sum_sp(4) + snow_m * wt_sp(4) * slant2

!                 if(idebug .eq. 1 .OR. cond_m .gt. 0.)then
                  if(.true.)then

                    if(cond_m .eq. 0.)then ! wasn't computed above
                      i1 = min(int(rinew_m),ni-1)
                      j1 = min(int(rjnew_m),nj-1) 
                      k1 = min(int(rk_m)   ,nk-1) 

                      fi = rinew_m - i1; i2=i1+1
                      fj = rjnew_m - j1; j2=j1+1
                      fk = rk_m    - k1; k2=k1+1
                    endif

                    transm_3d_m = sum(tri_coeff(:,:,:) * 
     1                            transm_3d(i1:i2,j1:j2,k1:k2))

!                   call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                                 ,transm_3d,transm_3d_m)

                    sum_odrad = sum_odrad + 
     1               (cvr_path * slant2 * transm_3d_m)       

                    do ic = 1,nc
                      transm_4d_m(ic) = sum(tri_coeff(:,:,:) * 
     1                              transm_4d(i1:i2,j1:j2,k1:k2,ic))
!                     call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
!    1                            ,transm_4d(:,:,:,ic),transm_4d_m(ic))
                    enddo ! ic

                    sum_odrad_c(:) = sum_odrad_c(:) + 
     1              (cvr_path * slant2 * transm_4d_m(:))
                  endif ! idebug .eq. 1 .OR. cond_m .gt. 0.

                  sum_clrrad = sum_clrrad + transm_3d(inew_m,jnew_m,k_m) 
     1                                    * (airmass1_h - airmass1_l)

                  aero_ext_coeff = aod_3d(inew_m,jnew_m,k_m)             
                  sum_aod     = sum_aod     + aero_ext_coeff * slant2

                  if(rk_m .lt. (rkstart + 1.0))then ! near topo
!                 if(.false.)then ! near topo
                    sum_aod_ill = sum_aod_ill + aero_ext_coeff * slant2
     1                          * transm_3d(inew_m,jnew_m,int(rk_m)+1)  
                    if(transm_3d(inew_m,jnew_m,int(rk_m)+1) .gt. .1)then
                        transm_3d_dir = 
     1                exp(log(transm_3d(inew_m,jnew_m,int(rk_m)+1))*10.)
                    else
                        transm_3d_dir = 0.
                    endif
                    sum_aod_ill_dir = sum_aod_ill_dir 
     1                         + aero_ext_coeff * slant2 * transm_3d_dir
                  else ! typical case far from topo
                    sum_aod_ill = sum_aod_ill + aero_ext_coeff * slant2
     1                          * transm_3d_m
                    if(transm_3d_m .gt. .1)then
                        transm_3d_dir = exp(log(transm_3d_m)*10.)
                    else
                        transm_3d_dir = 0.
                    endif
                    sum_aod_ill_dir = sum_aod_ill_dir 
     1                         + aero_ext_coeff * slant2 * transm_3d_dir
                  endif

                  if((cvr_path          .gt. 0.00) .AND.
     1               (cvr_path_sum      .le.  .013     .OR. ! tau ~1
     1                cvr_path_sum_last .eq. 0.  )            )then 
                      aod_2_cloud(ialt,jazi) = sum_aod
                  endif

!                 Determine backscatter threshold for averaging cloud rad
!                 This can be refined by determining the phase of the
!                 cloud so far in the ray traversal
                  if(cvr_path .gt. 0.00)then  
                    frac_liq = (cvr_path_sum_sp(1) + cvr_path_sum_sp(3))
     1                       / cvr_path_sum
                    tau_thr = 15. * frac_liq + 7. * (1.-frac_liq)
                    bks_thr = tau_thr / clwc2alpha
                    if(cvr_path_sum      .le.  bks_thr .OR. ! tau ~7-15
     1                 cvr_path_sum_last .eq. 0.            )then 
!                     Average value over cloud path (tau ~10)
                      r_cloud_rad(ialt,jazi) = sum_odrad / cvr_path_sum 
                      cloud_rad_c(:,ialt,jazi) 
     1                                = sum_odrad_c(:) / cvr_path_sum 
                    endif
                  endif

!                 Calculated weighted value of airmass to cloud
!                 Summation(airmass*cvrpath*trans) / Summation(cvrpath*trans)
                  if(cvr_path .gt. 0.00)then
                      am2cld_den = cvr_path 
     1                           * trans(clwc2alpha*cvr_path_sum)
                      am2cld_num = am2cld_den 
     1                           * 0.5 * (airmass1_l + airmass1_h)
                      sum_am2cld_num = sum_am2cld_num + am2cld_num
                      sum_am2cld_den = sum_am2cld_den + am2cld_den
                      airmass_2_cloud_3d(ialt,jazi) = sum_am2cld_num
     1                                              / sum_am2cld_den
                  endif

!                 if(ht_m - htstart .le. 1000.)then ! more accurate topo when low
                  if(ht_m   .le. topo_max_ht .AND. 
     1               altray .le. topo_max_ang     )then 
!                 if(ht_m   .le. topo_max_ht)then 
                    i1 = min(int(rinew_m),ni-1); fi = rinew_m - i1; i2=i1+1
                    j1 = min(int(rjnew_m),nj-1); fj = rjnew_m - j1; j2=j1+1

                    bi_coeff(1,1) = (1.-fi) * (1.-fj)
                    bi_coeff(2,1) = fi      * (1.-fj)
                    bi_coeff(1,2) = (1.-fi) *     fj 
                    bi_coeff(2,2) = fi      *     fj 
                    topo_m = sum(bi_coeff(:,:) * topo_a(i1:i2,j1:j2))
!                   call bilinear_laps(rinew_m,rjnew_m,ni,nj
!    1                                ,topo_a,topo_m)
                  else
                    topo_m = topo_a(inew_m,jnew_m)
                  endif

                  if(topo_m .gt. ht_m .AND. ihit_topo .eq. 0)then
                      ihit_topo = 1

!                     Land illumination related to terrain slope
                      if(sol_alt(inew_m,jnew_m)  .gt. 0. )then  
                        alt_norm_int = sum(bi_coeff(:,:) 
     1                               * alt_norm(i1:i2,j1:j2))
                        if(alt_norm_int .gt. 0. )then
                          solar_corr = sind(alt_norm_int) 
     1                               / sind(sol_alt(inew_m,jnew_m))
                          solar_corr = min(max(solar_corr,0.2),2.0) 
                        else ! land is shadowed
                          solar_corr = 0.2
                        endif
                      else ! sun is down     
                          solar_corr = 1.0
                      endif

                      topo_swi(ialt,jazi) = 
     1                  sum(bi_coeff(:,:) * swi_2d(i1:i2,j1:j2))
     1                                    * solar_corr
!                     topo_swi(ialt,jazi) = swi_2d(inew_m,jnew_m) 
!    1                                    * solar_corr

!                     City lights on the ground
                      if(sol_alt(inew_m,jnew_m) .lt. -4.)then  
                        topo_swi(ialt,jazi) = topo_swi(ialt,jazi) + 
     1                    sum(bi_coeff(:,:) * gnd_glow(i1:i2,j1:j2))
                        if(gnd_glow(inew_m,jnew_m) .gt. 0. .OR. 
     1                                                idebug .eq. 1)then
                          write(6,81)ialt,jazi,inew_m,jnew_m
     1                              ,gnd_glow(inew_m,jnew_m)
 81                       format(' gnd glow adding to topo_swi',4i5
     1                          ,f9.1)
                        endif
                      endif

                      do ic = 1,nc
                        topo_albedo(ic,ialt,jazi) = sum(bi_coeff(:,:) 
     1                                 * topo_albedo_2d(ic,i1:i2,j1:j2))
                      enddo ! ic
!                     topo_albedo(:,ialt,jazi) = 
!    1                    topo_albedo_2d(:,inew_m,jnew_m)

!                     Test for first segment
                      if(dxy1_l .eq. 0. .AND. htstart .eq. topo_sfc)then
                          airmass_2_topo_3d(ialt,jazi) = 1e-5
                          aod_2_topo(ialt,jazi)        = 1e-4
                          sum_aod_ill                  = 1e-4 
     1                            * transm_3d(inew_m,jnew_m,int(rk_m)+1)
                          sum_aod_ill_dir              = 1e-4 
     1                            * transm_3d(inew_m,jnew_m,int(rk_m)+1)
                      else
                          airmass_2_topo_3d(ialt,jazi) 
     1                                 = 0.5 * (airmass1_l + airmass1_h)
                          aod_2_topo(ialt,jazi) = sum_aod
                      endif

                      if(idebug .eq. 1)write(6,91)solar_corr  
     1                                    ,topo_swi(ialt,jazi)
     1                                    ,topo_albedo(:,ialt,jazi)
     1                                    ,aod_2_topo(ialt,jazi)
     1                                    ,airmass_2_topo_3d(ialt,jazi)
     1                                    ,gnd_glow(inew_m,jnew_m)
91                    format(' Hit topo ',f8.3,f8.2,3f8.3,3f8.3)
                  endif

                  cloud_od(ialt,jazi) = clwc2alpha*cvr_path_sum   
                  cloud_od_sp(ialt,jazi,:) 
     1                                = clwc2alpha*cvr_path_sum_sp(:)   
                  r_cloud_3d(ialt,jazi) = 1.-(exp(-cloud_od(ialt,jazi)))
                  if(idebug .eq. 1)write(6,101)dz1_l,dz1_h,dxy1_l,dxy1_h
     1                     ,cslant
     1                     ,rinew_h,rjnew_h,rk,ht_m,topo_m 
     1                     ,cvr_path*1e3
     1                     ,clwc_m*1e3,cice_m*1e3,rain_m*1e3,snow_m*1e3
     1                     ,slant2,cvr_path_sum
     1                     ,r_cloud_3d(ialt,jazi)
     1                     ,airmass_2_cloud_3d(ialt,jazi)
     1                     ,r_cloud_rad(ialt,jazi)
     1                     ,aero_ext_coeff
!    1                     ,transm_3d(inew_m,jnew_m,int(rk_m)+1)
     1                     ,transm_3d_m
     1                     ,sum_aod,sum_aod_ill,sum_aod_ill_dir
 101              format(2f8.1,2f9.1,a1,f6.1,f7.1,f6.2,2f8.1
     1                  ,1x,f7.4,2x,4f7.4,f10.1,2f11.4,f9.2,f9.4
     1                  ,f10.5,4f8.2)
                else
                  if(idebug .eq. 1)then
                      write(6,*)' out of bounds ',ialt,jazi
     1                    ,rinew_h,rjnew_h,rk,slant2,dxy1_l,dxy1_h,dz1_h   
                      if(rkdelt .ne. 0.)then
                          write(6,*)' aterm/bterm/cterm/discriminant='
     1                               ,aterm,bterm,cterm,discriminant   
!                         write(6,*)' argd1/argd2=',argd1,argd2
                      endif
                  endif
                  ihit_bounds = 1
                endif ! in horizontal/vertical domain

                if(rk_h - rkstart .gt. 0.0)then
!                 rkdelt = rkdelt2
!                 rkdelt = 0.5 * (rk_h - rkstart)
                  rkdelt = rkdelt * 1.1             
                  rkdelt = min(rkdelt,rkdelt2)
                endif
              enddo ! kk/rk

!           enddo ! k

          endif ! true

          aod_ill(ialt,jazi) = sum_aod_ill
          aod_ill_dir(ialt,jazi) = sum_aod_ill_dir
          aod_tot(ialt,jazi) = sum_aod

          if(r_cloud_3d(ialt,jazi) .gt. .5)then
              icloud = 1
          else
              icloud = 0
          endif

          icloud_tot = icloud_tot + icloud

!         l_process(ialt,jazi) = .true.

!         include '../lib/cloud/skyglow_phys.inc'

          if(sol_alt(i,j) .gt. 0.)then

!           Get clear sky daylight brightness ratio at point
!           (i.e. fraction of atmosphere illuminated by the sun)

            clear_radf_c(:,ialt,jazi) = 0.25 ! secondary scattering in cloud shadow
     1                                + 0.75 * (sum_clrrad/airmass1_h)
          endif

!         end include 'skyglow.inc'

         enddo ! jazi

!        Get clear sky twilight brightness in ring
         if(sol_alt(i,j) .le. 0.)then
!            write(6,*)' args: ',l_solar_eclipse,i4time,rlat,rlon
             call skyglow_phys(ialt,ialt,1,minazi,maxazi,jazi_delt
     1             ,minalt,maxalt,minazi,maxazi,idebug_a
     1             ,sol_alt(i,j),sol_azi(i,j),view_alt,view_az
     1             ,earth_radius,patm,aod_ray_eff,aod_ray_dir
     1             ,aero_scaleht,htstart,redp_lvl             ! I
     1             ,aod_ill                                   ! I (dummy)
     1             ,l_solar_eclipse,i4time,rlat,rlon          ! I
     1             ,clear_radf_c,ag_2d                        ! I (ag2d is dummy)
     1             ,clear_rad_c,elong                       ) ! O
         endif

         if(jazi_delt .eq. 2 .OR. jazi_delt .eq. 4 .OR. 
     1      jazi_delt .eq. 8               )then ! fill missing azimuths
          do jazi = minazi,maxazi
            call get_interp_parms(minazi,maxazi,jazi_delt,jazi       ! I
     1                           ,fm,fp,jazim,jazip,ir)              ! O
            if(ir .ne. 0)then
              r_cloud_3d(ialt,jazi) = 
     1         fm * r_cloud_3d(ialt,jazim) + fp * r_cloud_3d(ialt,jazip)
              cloud_od(ialt,jazi) = 
     1         fm * cloud_od(ialt,jazim)   + fp * cloud_od(ialt,jazip)
              cloud_od_sp(ialt,jazi,:) = 
     1         fm * cloud_od_sp(ialt,jazim,:) 
     1                                  + fp * cloud_od_sp(ialt,jazip,:)
              r_cloud_rad(ialt,jazi) = 
     1         fm * r_cloud_rad(ialt,jazim)+fp * r_cloud_rad(ialt,jazip)      
              cloud_rad_c(:,ialt,jazi) = 
     1         fm * cloud_rad_c(:,ialt,jazim)  
     1                                 + fp *cloud_rad_c(:,ialt,jazip)
              clear_rad_c(:,ialt,jazi) = 
     1         fm * clear_rad_c(:,ialt,jazim)   
     1                                 + fp * clear_rad_c(:,ialt,jazip)
              clear_radf_c(:,ialt,jazi) = 
     1         fm * clear_radf_c(:,ialt,jazim)   
     1                                 + fp * clear_radf_c(:,ialt,jazip)  
              airmass_2_cloud_3d(ialt,jazi) = 
     1                fm * airmass_2_cloud_3d(ialt,jazim) 
     1              + fp * airmass_2_cloud_3d(ialt,jazip)
              airmass_2_topo_3d(ialt,jazi) = 
     1                fm * airmass_2_topo_3d(ialt,jazim) 
     1              + fp * airmass_2_topo_3d(ialt,jazip)
              aod_ill(ialt,jazi) = 
     1                fm * aod_ill(ialt,jazim) 
     1              + fp * aod_ill(ialt,jazip)
              aod_ill_dir(ialt,jazi) = 
     1                fm * aod_ill_dir(ialt,jazim) 
     1              + fp * aod_ill_dir(ialt,jazip)
              aod_tot(ialt,jazi) = 
     1                fm * aod_tot(ialt,jazim) 
     1              + fp * aod_tot(ialt,jazip)
              topo_swi(ialt,jazi) = 
     1                fm * topo_swi(ialt,jazim) 
     1              + fp * topo_swi(ialt,jazip)
              topo_albedo(:,ialt,jazi) = 
     1                fm * topo_albedo(:,ialt,jazim) 
     1              + fp * topo_albedo(:,ialt,jazip)
            endif ! ir .ne. 0
          enddo ! jazi interp
         endif ! fill azimuth

        enddo ! ialt

        call get_idx(20.,minalt,alt_scale,ialt_min)
        call get_idx(90.,minalt,alt_scale,ialt_max)

        do ialt = ialt_min,ialt_max ! fill in missing alt rings
          call get_interp_parms(minalt,maxalt,idelt,ialt           ! I
     1                         ,fm,fp,ialtm,ialtp,ir)              ! O
          if(ir .ne. 0)then
            r_cloud_3d(ialt,:) =
     1       fm * r_cloud_3d(ialtm,:) + fp * r_cloud_3d(ialtp,:)
            cloud_od(ialt,:) =
     1       fm * cloud_od(ialtm,:)   + fp * cloud_od(ialtp,:)
            cloud_od_sp(ialt,:,:) =
     1       fm * cloud_od_sp(ialtm,:,:) + fp * cloud_od_sp(ialtp,:,:)
            r_cloud_rad(ialt,:) =
     1       fm * r_cloud_rad(ialtm,:)  + fp * r_cloud_rad(ialtp,:)
            cloud_rad_c(:,ialt,:) =
     1       fm * cloud_rad_c(:,ialtm,:) + fp * cloud_rad_c(:,ialtp,:)
            clear_rad_c(:,ialt,:) =
     1       fm * clear_rad_c(:,ialtm,:) + fp * clear_rad_c(:,ialtp,:)
            clear_radf_c(:,ialt,:) =
     1       fm * clear_radf_c(:,ialtm,:) + fp * clear_radf_c(:,ialtp,:)
            airmass_2_cloud_3d(ialt,:) =
     1           fm * airmass_2_cloud_3d(ialtm,:) 
     1         + fp * airmass_2_cloud_3d(ialtp,:)
            airmass_2_topo_3d(ialt,:) =
     1           fm * airmass_2_topo_3d(ialtm,:) 
     1         + fp * airmass_2_topo_3d(ialtp,:)
            aod_ill(ialt,:) =
     1           fm * aod_ill(ialtm,:) 
     1         + fp * aod_ill(ialtp,:)
            aod_ill_dir(ialt,:) =
     1           fm * aod_ill_dir(ialtm,:) 
     1         + fp * aod_ill_dir(ialtp,:)
            aod_tot(ialt,:) =
     1           fm * aod_tot(ialtm,:) 
     1         + fp * aod_tot(ialtp,:)
            topo_swi(ialt,:) =
     1         fm * topo_swi(ialtm,:)      + fp * topo_swi(ialtp,:)
            topo_albedo(:,ialt,:) =
     1         fm * topo_albedo(:,ialtm,:) + fp * topo_albedo(:,ialtp,:)
          endif ! ir .ne. 0
        enddo ! ialt

        write(6,*)' Number of rays with cloud = ',icloud_tot
     1           ,' out of ',91*361

        write(6,*)' Range of r_cloud_3d = ',minval(r_cloud_3d)
     1                                     ,maxval(r_cloud_3d)

        write(6,*)' Range of cloud_od = ',minval(cloud_od)
     1                                   ,maxval(cloud_od)

        do isp = 1,nsp
          write(6,*)' Range of cloud_od_sp = ',isp
     1                                    ,minval(cloud_od_sp(:,:,isp))
     1                                    ,maxval(cloud_od_sp(:,:,isp))
        enddo ! isp

        write(6,*)' Range of r_cloud_rad = ',minval(r_cloud_rad)
     1                                      ,maxval(r_cloud_rad)

        write(6,*)' Range of cloud_rad_c R =',minval(cloud_rad_c(1,:,:))
     1                                       ,maxval(cloud_rad_c(1,:,:))

        write(6,*)' Range of cloud_rad_c G =',minval(cloud_rad_c(2,:,:))
     1                                       ,maxval(cloud_rad_c(2,:,:))

        write(6,*)' Range of cloud_rad_c B =',minval(cloud_rad_c(3,:,:))
     1                                       ,maxval(cloud_rad_c(3,:,:))

        write(6,*)' Range of clear_rad_c 3 =',minval(clear_rad_c(3,:,:))
     1                                       ,maxval(clear_rad_c(3,:,:))

        write(6,*)' Range of clear_radf_c 3 ='
     1                                      ,minval(clear_radf_c(3,:,:))
     1                                      ,maxval(clear_radf_c(3,:,:))

        write(6,*)' Range of topo_swi = ',minval(topo_swi)
     1                                   ,maxval(topo_swi)

        write(6,*)' Range of topo_albedo 2 = '
     1                                       ,minval(topo_albedo(2,:,:))
     1                                       ,maxval(topo_albedo(2,:,:))

        write(6,*)' Range of aod_ill = ',minval(aod_ill)
     1                                  ,maxval(aod_ill)       

        write(6,*)' Range of aod_ill_dir = ',minval(aod_ill_dir)
     1                                      ,maxval(aod_ill_dir)       

        transm_obs = transm_3d(i,j,int(rkstart)+1)
        write(6,*)' transm of observer is ',i,j,int(rkstart)+1
     1                                     ,transm_obs

        I4_elapsed = ishow_timer()
 
        return
        end
