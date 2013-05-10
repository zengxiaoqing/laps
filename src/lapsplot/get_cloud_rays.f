        
        subroutine get_cloud_rays(i4time,clwc_3d,cice_3d,heights_3d
     1                             ,pres_3d,topo_sfc
     1                             ,r_cloud_3d,airmass_2_cloud_3d
     1                             ,ni,nj,nk,i,j
     1                             ,view_alt,view_az 
     1                             ,grid_spacing_m,r_missing_data)

        include 'trigd.inc'

!       Determine shadow regions of the input 3-D cloud array

        real clwc_3d(ni,nj,nk)
        real cice_3d(ni,nj,nk)
        real cond_3d(ni,nj,nk)
        real heights_3d(ni,nj,nk)
        real heights_1d(nk)
        real pres_3d(ni,nj,nk)
        real pres_1d(nk)
        real view_alt(ni,nj)
        real view_az(ni,nj)

        integer isky(0:90,0:360)
        real r_shadow_3d(0:90,0:360)
        real r_cloud_3d(0:90,0:360)
        real airmass_2_cloud_3d(0:90,0:360)

        I4_elapsed = ishow_timer()

        write(6,*)' Subroutine get_cloud_rays ',i,j

        pstd = 101325.
        airmass_2_cloud_3d = 0

        ishadow_tot = 0                          

        write(6,*)' range of clwc,cice is ',maxval(clwc_3d)
     1                                     ,maxval(cice_3d)

        write(6,*)' range of heights is ',minval(heights_3d)
     1                                   ,maxval(heights_3d)

        cond_3d = clwc_3d + cice_3d

        ri = i
        rj = j

        do k = 1,nk-1
            if(heights_3d(i,j,k)   .le. topo_sfc .AND.
     1         heights_3d(i,j,k+1) .ge. topo_sfc      )then
                ksfc = k
                write(6,*)' ksfc = ',ksfc
            endif
        enddo

        heights_1d(:) = heights_3d(i,j,:)
        pres_1d(:)    = pres_3d(i,j,:)

        do ialt = 0,90

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

          if((jazi .eq. 45 .or. jazi .eq. 225) .AND.
     1        ialt .eq. (ialt/10)*10                 )then
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

          if(.true.)then
            do k = ksfc,ksfc
              r_shadow_3d(ialt,jazi) = 0.
              
              ri1 = ri
              rj1 = rj

              ht_sfc = heights_1d(k)

              if(view_altitude_deg .le. 10.)then
                  rkdelt = 0.5
              else
                  rkdelt = 1.0
              endif

!             arg = max(view_altitude_deg,1.0)
!             rkdelt = tand(90. - arg)                     
!             rkdelt = max(min(rkdelt,2.0),0.5)
              
              if(idebug .eq. 1)write(6,*)
              if(idebug .eq. 1)write(6,*)
     1        ' Trace the slant path (ialt/jazi/jdelt,rkdelt) '
     1                               ,ialt,jazi,jazi_delt,rkdelt

              slant1_h = 0.
              airmass1_h = 0.

!             do kk = ksfc+1,nk
              do rk = ksfc+rkdelt,float(nk),rkdelt

!               rk = kk

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
     1               + heights_1d(kk_l+1) * frac_l

                pr_l = pres_1d(kk_l) * (1. - frac_l) 
     1               + pres_1d(kk_l+1) * frac_l

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

                dz1_l = ht_l - ht_sfc        
                dz1_h = ht_h - ht_sfc          
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
                  dxy1_l = (sqrt(4.*aterm*cterm + bterm**2.) - bterm)
     1                                  / (2.*aterm)
                  cterm = dz1_h                                                
                  dxy1_h = (sqrt(4.*aterm*cterm + bterm**2.) - bterm)
     1                                  / (2.*aterm)
                  dxy2 = dxy1_h - dxy1_l
                endif

                dx1_l = dxy1_l * xcos
                dy1_l = dxy1_l * ycos

                dx1_h = dxy1_h * xcos
                dy1_h = dxy1_h * ycos

                slant2 = sqrt(dxy2**2 + dz2**2) ! layer
                slant1_l = slant1_h
                slant1_h = slant1_h + slant2

                airmassv = (pr_l - pr_h) / pstd
                airmass2 = airmassv * (slant2 / dz2)
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

                rni = ni
                rnj = nj
              
                if(rinew_h .ge. 1. .and. rinew_h .le. rni .AND. 
     1             rjnew_h .ge. 1. .and. rjnew_h .le. rnj      )then 

                  if(.true.)then
                    call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                 ,cond_3d,cond_m)
                  elseif(.false.)then
                    call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                 ,clwc_3d,clwc_m)

                    call trilinear_laps(rinew_m,rjnew_m,rk_m,ni,nj,nk
     1                                 ,cice_3d,cice_m)
                  else
                    clwc_m = 
     1              clwc_3d(nint(rinew_m),nint(rjnew_m),nint(rk_m))
                    cice_m = 
     1              cice_3d(nint(rinew_m),nint(rjnew_m),nint(rk_m))
                  endif

                  cvr_path = cond_m                                 
                  cvr_path_sum_last = cvr_path_sum
                  cvr_path_sum      = cvr_path_sum + cvr_path * slant2       

                  if(cvr_path_sum      .gt. .05 .and.
     1               cvr_path_sum_last .lt. .05       )then ! optical depth ~1
                      airmass_2_cloud_3d(ialt,jazi) 
     1                                 = 0.5 * (airmass1_l + airmass1_h)
                  endif

                  r_shadow_3d(ialt,jazi) = 1. - (exp(-20.*cvr_path_sum))
                  if(idebug .eq. 1)write(6,101)dz1_l,dz1_h,dxy1_l,dxy1_h
     1                     ,rinew_h,rjnew_h,rk 
     1                     ,cvr_path,slant2,cvr_path_sum
     1                     ,r_shadow_3d(ialt,jazi)
 101              format(2f8.1,2f9.1,4x,2f7.1,f6.1,f12.7,f10.1,2f11.4)
                else
                  if(idebug .eq. 1)write(6,*)' out of bounds ',ialt,jazi
     1                            ,rinew_h,rjnew_h,rk
                endif
              enddo ! kk/rk

            enddo ! k

          endif ! high enough sun to use vis   

          if(r_shadow_3d(ialt,jazi) .gt. .5)then
!         if(r_shadow_3d(ialt,jazi) .ge. 0.000000)then
              ishadow = 1
              r_cloud_3d(ialt,jazi) = r_shadow_3d(ialt,jazi) ! 1              
          else
              ishadow = 0
              r_cloud_3d(ialt,jazi) = r_shadow_3d(ialt,jazi) ! 0
          endif

          ishadow_tot = ishadow_tot + ishadow

          if(idebug .eq. 1 .OR. 
     1       (ishadow .eq. 1 .AND. ishadow_tot .le. 10) )then
              write(6,*)'ialt,jazi,ishadow',ialt,jazi,ishadow
          endif

         enddo ! jazi

         if(jazi_delt .eq. 2)then ! fill in missing azimuths
          do jazi = 1,359,jazi_delt
              r_cloud_3d(ialt,jazi) = 
     1         0.5 * (r_cloud_3d(ialt,jazi-1) + r_cloud_3d(ialt,jazi+1))
          enddo
         endif

        enddo ! ialt

        do ialt = 21,89,2 ! fill in missing alt rings
            r_cloud_3d(ialt,:) =
     1         0.5 * (r_cloud_3d(ialt-1,:) + r_cloud_3d(ialt+1,:))
        enddo ! ialt

        write(6,*)' Number of rays with cloud = ',ishadow_tot
     1           ,' out of ',91*361

        write(6,*)' Range of r_shadow_3d = ',minval(r_shadow_3d)
     1                                      ,maxval(r_shadow_3d)

        write(6,*)' Range of r_cloud_3d = ',minval(r_cloud_3d)
     1                                     ,maxval(r_cloud_3d)

        I4_elapsed = ishow_timer()
 
        return
        end
