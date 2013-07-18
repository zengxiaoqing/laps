        
        subroutine get_cloud_rays(i4time,clwc_3d,cice_3d,heights_3d
     1                           ,rain_3d,snow_3d
     1                           ,pres_3d,topo_sfc,topo_a,swi_2d
     1                           ,topo_swi,topo_albedo_2d,topo_albedo
     1                           ,aod_2_cloud,aod_2_topo
     1                           ,r_cloud_3d,cloud_od
     1                           ,cloud_rad,cloud_rad_c
     1                           ,clear_rad_c
     1                           ,airmass_2_cloud_3d,airmass_2_topo_3d
     1                           ,ni,nj,nk,i,j
     1                           ,view_alt,view_az,sol_alt,sol_azi 
     1                           ,minalt,maxalt
     1                           ,grid_spacing_m,r_missing_data)

        use mem_namelist, ONLY: earth_radius

        include 'trigd.inc'

        angdif(X,Y)=MOD(X-Y+540.,360.)-180.

        parameter (rpd = 3.14159 / 180.)

!       Determine shadow regions of the input 3-D cloud array

        parameter (nc = 3)

        real clwc_3d(ni,nj,nk) ! kg/m**3
        real cice_3d(ni,nj,nk) ! kg/m**3
        real rain_3d(ni,nj,nk) ! kg/m**3
        real snow_3d(ni,nj,nk) ! kg/m**3
        real cond_3d(ni,nj,nk) ! kg/m**3 (effective LWC)
        real aod_3d(ni,nj,nk)  ! aerosol optical depth (per vertical atmosphere)
        real heights_3d(ni,nj,nk) ! MSL
        real transm_3d(ni,nj,nk)
        real transm_4d(nc,ni,nj,nk)
        real heights_1d(nk)
        real topo_a(ni,nj)
        real swi_2d(ni,nj)
        real topo_albedo_2d(nc,ni,nj)
        real pres_3d(ni,nj,nk)
        real pres_1d(nk)
        real view_alt(ni,nj)
        real view_az(ni,nj)

        real sol_alt(ni,nj)
        real sol_azi(ni,nj)

        integer isky(minalt:maxalt,0:360)
        real r_cloud_3d(minalt:maxalt,0:360)        ! cloud opacity
        real cloud_od(minalt:maxalt,0:360)          ! cloud optical depth
        real cloud_rad(minalt:maxalt,0:360)         ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,minalt:maxalt,0:360)    ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,minalt:maxalt,0:360)    ! integrated fraction of air illuminated by the sun along line of sight
        real airmass_2_cloud_3d(minalt:maxalt,0:360)
        real airmass_2_topo_3d(minalt:maxalt,0:360)
        real topo_swi(minalt:maxalt,0:360)
        real topo_albedo(nc,minalt:maxalt,0:360)
        real aod_2_cloud(minalt:maxalt,0:360)
        real aod_2_topo(minalt:maxalt,0:360)
        real sum_odrad_c(nc)

        real*8 xplane,yplane,zplane
        real*8 xray,yray,zray,angle_r 

        character*1 cslant

        I4_elapsed = ishow_timer()

        write(6,*)' Subroutine get_cloud_rays ',i,j

        pstd = 101325.
        airmass_2_cloud_3d = 0.
        airmass_2_topo_3d = 0.
        topo_swi = 0.
        cloud_rad = 1.
        cloud_rad_c = 1.

        ishadow_tot = 0                          

        write(6,*)' range of clwc,cice is ',maxval(clwc_3d)
     1                                     ,maxval(cice_3d)

        write(6,*)' range of heights is ',minval(heights_3d)
     1                                   ,maxval(heights_3d)

        write(6,*)' call get_cloud_rad'

        call get_cloud_rad(sol_alt,sol_azi,clwc_3d,cice_3d
     1                    ,rain_3d,snow_3d
     1                    ,heights_3d,transm_3d,transm_4d,i,j,ni,nj,nk)

        I4_elapsed = ishow_timer()

        cond_3d = clwc_3d +        cice_3d * 0.5 
     1          + rain_3d * 0.02 + snow_3d * 0.05

        ri = i
        rj = j

        kstart = 0 ! 0 means sfc, otherwise level of start

        if(kstart .eq. 0)then ! start from sfc
          htstart = topo_sfc  ! MSL
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
        else ! start aloft
          rkstart = kstart
          htstart = heights_3d(i,j,kstart)
          write(6,*)' Start aloft at k = ',kstart
        endif

        heights_1d(:) = heights_3d(i,j,:)
        pres_1d(:)    = pres_3d(i,j,:)

        do ialt = minalt,maxalt

         if(ialt .ge. 20)then
             if(ialt .eq. (ialt/2)*2)then
                 jazi_delt = 2
             else
                 jazi_delt = 360
             endif
         else
             jazi_delt = 1
         endif

         do jazi = 0,360,jazi_delt

          if((jazi .eq. 225 .or. jazi .eq. 45) .AND.
!    1        (ialt .eq. (ialt/10)*10 .or. ialt .le. 5) )then
     1        (ialt .eq. 12 .or. ialt .le. 6) )then
              idebug = 1
          else
              idebug = 0
          endif

!         Trace towards sky from each grid point
          view_altitude_deg = max(float(ialt),0.) ! handle Earth curvature later
          view_azi_deg = jazi

!         Get direction cosines based on azimuth
          xcos = sind(view_azi_deg)
          ycos = cosd(view_azi_deg)

          ishadow = 0
          cvr_path_sum = 0.
          sum_odrad = 0.
          sum_odrad_c = 0.

          if(.true.)then
!           do k = ksfc,ksfc
              r_cloud_3d(ialt,jazi) = 0.
              cloud_od(ialt,jazi) = 0.
              
              ri1 = ri
              rj1 = rj

!             ht_sfc = heights_1d(ksfc)
!             ht_sfc = topo_sfc       

              grid_factor = grid_spacing_m / 3000.

              if(view_altitude_deg .lt. -4.5)then
                  rkdelt1 = -1.00 * grid_factor
                  rkdelt2 = -1.00 * grid_factor
              elseif(view_altitude_deg .lt. 0.0)then
                  rkdelt1 =  0.0  * grid_factor
                  rkdelt2 =  0.0  * grid_factor
              elseif(view_altitude_deg .le. 0.5)then
                  rkdelt1 = 0.01 * grid_factor
                  rkdelt2 = 0.10 * grid_factor
              elseif(view_altitude_deg .le. 4.)then
                  rkdelt1 = 0.10 * grid_factor
                  rkdelt2 = 0.25 * grid_factor
              elseif(view_altitude_deg .le. 10.)then
                  rkdelt1 = 0.50 * grid_factor
                  rkdelt2 = 0.50 * grid_factor
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
                write(6,11)ialt,jazi,jazi_delt,rkdelt,i,j
11              format(
     1          'Trace the slant path (ialt/jazi/jdelt/rkdelt/i/j):'
     1                ,3i5,f6.2,2i5)                                     
                write(6,12)                                                            
12              format('   dz1_l   dz1_h   dxy1_l   dxy1_h  rin',
     1           'ew  rjnew   rk    ht_m   topo_m  ',
     1           ' path     lwc    ice    rain   snow      slant',
     1           '      sum     cloudfrac  airmass cloud_rad')
              endif

!             Initialize ray
              dz1_h = 0.
              slant1_h = 0.
              airmass1_h = 0.
              ihit_topo = 0
              rk_h = rkstart

              rkdelt = rkdelt1

              rk = rkstart          
              do while(rk .le. float(nk) - rkdelt 
     1                                 .AND. ihit_topo .eq. 0)

                if(rkdelt .ne. 0.)then ! trace by pressure levels
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

                  if(rk_h .lt. float(nk))then
                    ht_h = heights_1d(kk_h) * (1. - frac_h) 
     1                   + heights_1d(kk_h+1) * frac_h

                    pr_h = pres_1d(kk_h) * (1. - frac_h) 
     1                   + pres_1d(kk_h+1) * frac_h

                  elseif(rk_h .eq. float(nk))then
                    ht_h = heights_1d(kk_h)

                    pr_h = pres_1d(kk_h)
                  else
                    write(6,*)' ERROR: rk_h is too high ',rk_h,nk
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
                    radius_earth_8_thirds = 6371.e3 * 2.6666666
                    aterm = 1. / radius_earth_8_thirds
                    bterm = tand(min(view_altitude_deg,89.9))
                    cterm = dz1_l                                            
                    if(cterm .gt. 0.)then ! catch start point
                      dxy1_l = 
     1                       (sqrt(4.*aterm*cterm + bterm**2.) - bterm)   
     1                                  / (2.*aterm)
                    else
                      dxy1_l = 0.
                    endif

                    cterm = dz1_h                                          
                    dxy1_h = (sqrt(4.*aterm*cterm + bterm**2.) - bterm) 
     1                                  / (2.*aterm)
                    dxy2 = dxy1_h - dxy1_l
                  endif

                  slant2 = sqrt(dxy2**2 + dz2**2) ! layer
                  slant1_l = slant1_h
                  slant1_h = slant1_h + slant2

                  airmassv = (pr_l - pr_h) / pstd
                  airmass2 = airmassv * (slant2 / dz2)

                else  ! Trace by slant range (rkdelt = 0. )
                  slant2 = grid_spacing_m
                  slant1_l = slant1_h
                  slant1_h = slant1_h + slant2

!                 Determine height from distance and elev
                  dz1_l = dz1_h
                  dz1_h = slant1_h * sind(float(ialt)) 
     1                  + slant1_h**2 / (2. * radius_earth_8_thirds)

                  rk_l = rk_h
                  ht_l = htstart + dz1_l
                  ht_h = htstart + dz1_h

                  do k = 1,nk-1
                    if(heights_1d(k)   .le. ht_h .AND.
     1                 heights_1d(k+1) .ge. ht_h      )then
                        frach = (ht_h - heights_1d(k)) / 
     1                          (heights_1d(k+1) - heights_1d(k))      
                        rk_h = float(k) + frach
                    endif
                  enddo

                  pr_h = (pres_1d(k)*(1.-frach)) + (pres_1d(k+1)*frach)

                  airmass2 = (slant2 / 8000.) * (pr_h / pstd)

                endif ! rkdelt .ne. 0.

                if(dxy2 .gt. grid_spacing_m)then
!                 if(idebug .eq. 1)then
!                   write(6,*)' dxy2 is more than the grid spacing'
!                 endif
                  cslant = '*'
                else
                  cslant = ' '
                endif

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

                rni = ni
                rnj = nj
              
                if(rinew_h .ge. 1. .and. rinew_h .le. rni .AND. 
     1             rjnew_h .ge. 1. .and. rjnew_h .le. rnj      )then 

                  if(.true.)then
                    call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                 ,cond_3d,cond_m)
                    if(idebug .eq. 1 .OR. cond_m .gt. 0.)then
                      call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                   ,clwc_3d,clwc_m)
                      call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                   ,cice_3d,cice_m)
                      call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                   ,rain_3d,rain_m)
                      call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                   ,snow_3d,snow_m)
                    endif

!                   nf = 1
!                   call trilinear_laps_1d(rinew_m,rjnew_m,rk_m
!    1                                    ,ni,nj,nk,nf
!    1                                    ,cond_3d,cond_m)

                  else
                    clwc_m = 
     1              clwc_3d(nint(rinew_m),nint(rjnew_m),nint(rk_m))
                    cice_m = 
     1              cice_3d(nint(rinew_m),nint(rjnew_m),nint(rk_m))
                  endif

                  if(ht_m - htstart .le. 1000.)then ! more accurate topo when low
                    call bilinear_laps(rinew_m,rjnew_m,ni,nj
     1                                ,topo_a,topo_m)
                  else
                    topo_m = topo_a(nint(rinew_m),nint(rjnew_m))
                  endif

                  cvr_path = cond_m                                 
                  cvr_path_sum_last = cvr_path_sum
                  cvr_path_sum      = cvr_path_sum + cvr_path * slant2       
                  sum_odrad = sum_odrad + 
     1               (cvr_path * slant2 * 
     1                transm_3d(nint(rinew_m),nint(rjnew_m),nint(rk_m)))       

                  sum_odrad_c(:) = sum_odrad_c(:) + 
     1               (cvr_path * slant2 * 
     1              transm_4d(:,nint(rinew_m),nint(rjnew_m),nint(rk_m)))       

                  if((cvr_path          .gt. 0.00) .AND.
     1               (cvr_path_sum      .le.  .013     .OR. ! tau ~1
     1                cvr_path_sum_last .eq. 0.  )            )then 
                      airmass_2_cloud_3d(ialt,jazi) 
     1                                 = 0.5 * (airmass1_l + airmass1_h)

!                     Average value over cloud path (tau ~1)
                      cloud_rad(ialt,jazi) = sum_odrad / cvr_path_sum 
                      cloud_rad_c(:,ialt,jazi) 
     1                                = sum_odrad_c(:) / cvr_path_sum 
                  endif

                  if(topo_m .gt. ht_m .AND. ihit_topo .eq. 0)then
                      ihit_topo = 1
                      airmass_2_topo_3d(ialt,jazi) 
     1                                 = 0.5 * (airmass1_l + airmass1_h)
                      topo_swi(ialt,jazi) = 
     1                    swi_2d(nint(rinew_m),nint(rjnew_m))
                      topo_albedo(:,ialt,jazi) = 
     1                    topo_albedo_2d(:,nint(rinew_m),nint(rjnew_m))
                      if(idebug .eq. 1)write(6,*)' Hit topo '
     1                                ,topo_swi(ialt,jazi)
     1                                ,topo_albedo(:,ialt,jazi)
                  endif

!                 value of 75. from clwc constants in 'get_cloud_rad'
                  cloud_od(ialt,jazi) = 75.*cvr_path_sum   
                  r_cloud_3d(ialt,jazi) = 1.-(exp(-cloud_od(ialt,jazi)))
                  if(idebug .eq. 1)write(6,101)dz1_l,dz1_h,dxy1_l,dxy1_h
     1                     ,cslant
     1                     ,rinew_h,rjnew_h,rk,ht_m,topo_m 
     1                     ,cvr_path*1e3
     1                     ,clwc_m*1e3,cice_m*1e3,rain_m*1e3,snow_m*1e3
     1                     ,slant2,cvr_path_sum
     1                     ,r_cloud_3d(ialt,jazi)
     1                     ,airmass_2_cloud_3d(ialt,jazi)
     1                     ,cloud_rad(ialt,jazi)
 101              format(2f8.1,2f9.1,a1,f6.1,f7.1,f6.2,2f8.1
     1                  ,1x,f7.4,2x,4f7.4,f10.1,2f11.4,f9.2,f9.4)
                else
                  if(idebug .eq. 1)write(6,*)' out of bounds ',ialt,jazi
     1                            ,rinew_h,rjnew_h,rk
                endif

                if(rk_h - rkstart .gt. 0.0)then
!                 rkdelt = rkdelt2
!                 rkdelt = 0.5 * (rk_h - rkstart)
                  rkdelt = rkdelt * 1.1             
                  rkdelt = min(rkdelt,rkdelt2)
                endif
              enddo ! kk/rk

!           enddo ! k

          endif ! true

          if(r_cloud_3d(ialt,jazi) .gt. .5)then
              ishadow = 1
          else
              ishadow = 0
          endif

          ishadow_tot = ishadow_tot + ishadow

!         Estimate portion of ray in illuminated region of atmosphere
          if(sol_alt(i,j) .gt. 0.)then
            clear_rad_c(:,ialt,jazi) = 1. ! assume all of atmosphere is illuminated by the sun
            angle_plane = 0.; dist_ray_plane = 0.; ht_ray_plane = 0.
          else ! sun below horizon
            alt_plane = 90. - abs(sol_alt(i,j))
            azi_plane = sol_azi(i,j)
            xplane = cosd(azi_plane) * cosd(alt_plane)
            yplane = sind(azi_plane) * cosd(alt_plane)
            zplane =                   sind(alt_plane)
            xray   = cosd(float(jazi)) * cosd(float(ialt))
            yray   = sind(float(jazi)) * cosd(float(ialt))
            zray   =                     sind(float(ialt))
            call anglevectors(xplane,yplane,zplane
     1                       ,xray,yray,zray,angle_r)
            angle_plane = 90. - (angle_r / rpd) ! angle between light ray and plane                              
            if(angle_plane .gt. 0.)then ! assume part of atmosphere is illuminated by the sun
              horz_dep_r = -sol_alt(i,j) * rpd                                                  
              dist_pp_plane = horz_dep_r**2 * earth_radius / 2.0 ! approx perpendicular dist
              dist_ray_plane = dist_pp_plane / sind(angle_plane) ! distance along ray
              ht_ray_plane = dist_ray_plane * sind(float(ialt)) 
     1                     + dist_ray_plane**2 / (2. * earth_radius)
              if(ht_ray_plane .le. 99000.)then
                frac_airmass_lit = min(ZtoPsa(ht_ray_plane) / 1013.,1.0)
              else
                frac_airmass_lit = 0.
              endif
              z = 90. - float(ialt)
              airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))
              airmass_lit = airmass * frac_airmass_lit      
              clear_int = max(min(airmass_lit,1.0),.001)                                                
            else ! in Earth's shadow
              dist_ray_plane = -999.9
              ht_ray_plane = -999.9
              airmass_lit = 0.
              clear_int = .001                                                                           
            endif

            sat_ramp = min(-sol_alt(i,j) / 3.0,1.0)
            bluness = (1.0-saturation) / 2.0 
     1              + exp(-airmass_lit*0.5) * saturation ! range from 0.2 to 0.8 
            redness = 1.0 - bluness
!           HSI
            hue = exp(-airmass_lit*0.2) ! 0:R 1:B
            clear_rad_c(1,ialt,jazi) = hue                           ! Hue
            clear_rad_c(2,ialt,jazi) = abs(hue-0.4) * 1.0 * sat_ramp ! Saturation
            clear_rad_c(3,ialt,jazi) = clear_int                     ! Intensity

          endif ! sun above horizon

          if(idebug .eq. 1)then 
            write(6,111)ialt,jazi,airmass_2_cloud_3d(ialt,jazi)
     1                ,cloud_rad(ialt,jazi),sol_alt(i,j)
     1                ,angle_plane,dist_ray_plane,ht_ray_plane
     1                ,clear_rad_c(2,ialt,jazi)
111         format(
     1     'ialt/jazi/airm2cld/cldrad/salt/ang_pln/ds_ray/ht_ray/clrrad'
     1            ,2i5,2f10.6,2f8.2,2f10.1,f8.4)
          endif

         enddo ! jazi

         if(jazi_delt .eq. 2)then ! fill in missing azimuths
          do jazi = 1,359,jazi_delt
              r_cloud_3d(ialt,jazi) = 
     1         0.5 * (r_cloud_3d(ialt,jazi-1) + r_cloud_3d(ialt,jazi+1))
              cloud_od(ialt,jazi) = 
     1         0.5 * (cloud_od(ialt,jazi-1)   + cloud_od(ialt,jazi+1))
              cloud_rad(ialt,jazi) = 
     1         0.5 * (cloud_rad(ialt,jazi-1)  + cloud_rad(ialt,jazi+1))
              cloud_rad_c(:,ialt,jazi) = 
     1         0.5 * (cloud_rad_c(:,ialt,jazi-1)  
     1                                     + cloud_rad_c(:,ialt,jazi+1))
              clear_rad_c(:,ialt,jazi) = 
     1         0.5 * (clear_rad_c(:,ialt,jazi-1)   
     1                                     + clear_rad_c(:,ialt,jazi+1))
              airmass_2_cloud_3d(ialt,jazi) = 
     1         0.5 * (airmass_2_cloud_3d(ialt,jazi-1) 
     1              + airmass_2_cloud_3d(ialt,jazi+1))
              airmass_2_topo_3d(ialt,jazi) = 
     1         0.5 * (airmass_2_topo_3d(ialt,jazi-1) 
     1              + airmass_2_topo_3d(ialt,jazi+1))
              topo_swi(ialt,jazi) = 
     1         0.5 * (topo_swi(ialt,jazi-1) 
     1              + topo_swi(ialt,jazi+1))
          enddo
         endif

        enddo ! ialt

        do ialt = 21,89,2 ! fill in missing alt rings
            r_cloud_3d(ialt,:) =
     1         0.5 * (r_cloud_3d(ialt-1,:) + r_cloud_3d(ialt+1,:))
            cloud_od(ialt,:) =
     1         0.5 * (cloud_od(ialt-1,:)   + cloud_od(ialt+1,:))
            cloud_rad(ialt,:) =
     1         0.5 * (cloud_rad(ialt-1,:)     + cloud_rad(ialt+1,:))
            cloud_rad_c(:,ialt,:) =
     1         0.5 * (cloud_rad_c(:,ialt-1,:) + cloud_rad_c(:,ialt+1,:))
            clear_rad_c(:,ialt,:) =
     1         0.5 * (clear_rad_c(:,ialt-1,:) + clear_rad_c(:,ialt+1,:))
            airmass_2_cloud_3d(ialt,:) =
     1         0.5 * (airmass_2_cloud_3d(ialt-1,:) 
     1              + airmass_2_cloud_3d(ialt+1,:))
            airmass_2_topo_3d(ialt,:) =
     1         0.5 * (airmass_2_topo_3d(ialt-1,:) 
     1              + airmass_2_topo_3d(ialt+1,:))
            topo_swi(ialt,:) =
     1         0.5 * (topo_swi(ialt-1,:) 
     1              + topo_swi(ialt+1,:))
        enddo ! ialt

        write(6,*)' Number of rays with cloud = ',ishadow_tot
     1           ,' out of ',91*361

        write(6,*)' Range of r_cloud_3d = ',minval(r_cloud_3d)
     1                                     ,maxval(r_cloud_3d)

        write(6,*)' Range of cloud_od = ',minval(cloud_od)
     1                                   ,maxval(cloud_od)

        write(6,*)' Range of cloud_rad = ',minval(cloud_rad)
     1                                    ,maxval(cloud_rad)

        write(6,*)' Range of cloud_rad_c R =',minval(cloud_rad_c(1,:,:))
     1                                       ,maxval(cloud_rad_c(1,:,:))

        write(6,*)' Range of cloud_rad_c G =',minval(cloud_rad_c(2,:,:))
     1                                       ,maxval(cloud_rad_c(2,:,:))

        write(6,*)' Range of cloud_rad_c B =',minval(cloud_rad_c(3,:,:))
     1                                       ,maxval(cloud_rad_c(3,:,:))

        I4_elapsed = ishow_timer()
 
        return
        end
