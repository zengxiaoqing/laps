

        subroutine skyglow_phys(ialt_start,ialt_end,ialt_delt &! I
                   ,jazi_start,jazi_end,jazi_delt             &! I
                   ,minalt,maxalt,minazi,maxazi,idebug_a      &! I
                   ,sol_alt,sol_azi,view_alt,view_az          &! I
                   ,earth_radius,patm,aod_ray,aero_scaleht    &! I
                   ,htmsl,redp_lvl                            &! I
                   ,l_solar_eclipse,i4time,rlat,rlon          &! I
                   ,clear_rad_c,elong                       )  ! O

        include 'trigd.inc'

        use mem_namelist, ONLY: aod_bin, aod_asy 

!       Statement Functions
        trans(od) = exp(-od)
!       brt(a) = (1.0 - exp(-.14 * a))/.14 ! relative sky brightness (max~7)
        brt(a) = (1.0 - exp(-.14 * a))     ! relative sky brightness (max=1)
        brtf(am,od_per_am) = 1.0 - exp(-am*od_per_am) ! rel. sky brightness (max=1)
        brto(od) = 1.0 - exp(-od)                     ! rel. sky brightness (max=1)

        angdif(X,Y)=MOD(X-Y+540.,360.)-180.
        angleunitvectors(a1,a2,a3,b1,b2,b3) = acosd(a1*b1+a2*b2+a3*b3)
        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        real linecyl,linecylp
        linecylp(xcos,ycos,x1,r) = (-(2.*xcos*x1) + sqrt((2.*xcos*x1)**2 - 4. * (xcos**2 + ycos**2) * (x1**2 - r**2))) &
                                 /  (2. * (xcos**2 + ycos**2))

        include 'rad.inc'

        real twi_trans_c(nc)           ! transmissivity

        real mie, alpha_g, alpha_a

        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi) ! integrated fraction of air illuminated by the sun along line of sight
                                           ! (consider Earth's shadow + clouds)
        real aod_ray(minalt:maxalt,minazi:maxazi)
        real elong(minalt:maxalt,minazi:maxazi)
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

        logical l_solar_eclipse

        if(sol_alt .gt. 0.)then
            write(6,*)' skyglow_phys: i4time is ',i4time,l_solar_eclipse
            write(6,*)' aod_bin = ',aod_bin
            write(6,*)' aod_asy = ',aod_asy
        endif

        do ialt = ialt_start,ialt_end,ialt_delt
  
         altray = view_alt(ialt,jazi_start)
         aero_refht = redp_lvl
         call get_airmass(altray,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius &            ! I
                         ,ag,ao,aa)                 ! O

         do jazi = jazi_start,jazi_end,jazi_delt

          altray = view_alt(ialt,jazi)
          view_altitude_deg = altray

          view_azi_deg = view_az(ialt,jazi)

!         include 'skyglow_phys.inc'
!         include 'skyglow_phys.inc'

          xs = cosd(sol_alt) * cosd(sol_azi)
          ys = cosd(sol_alt) * sind(sol_azi)
          zs = sind(sol_alt)

          xo = cosd(altray) * cosd(view_azi_deg)
          yo = cosd(altray) * sind(view_azi_deg)
          zo = sind(altray)

          elong(ialt,jazi) = angleunitvectors(xs,ys,zs,xo,yo,zo)

          idebug = idebug_a(ialt,jazi)

          if(sol_alt .gt. 0.)then

!           Fraction of atmosphere illuminated by the sun
!           clear_rad_c(3,ialt,jazi) = 1. 
!           clear_rad_c(3,ialt,jazi) = 0.25&! secondary scattering in cloud shadow
!                                    + 0.75 * (sum_clrrad/airmass1_h)
            angle_plane = 0.; dist_ray_plane = 0.; ht_ray_plane = 0.

!           Illumination using Rayleigh and HG phase functions
            z = (90. - altray) + refractd_app(altray      ,patm)
            airmass_g = ag
            sb_corr = 2.0 * (1.0 - (sind(sol_alt)**0.5))
            day_int = 3e9 / 10.**(0.4 * sb_corr)

!           HG illumination
            hg2 = aod_bin(1) * hg(aod_asy(1),elong(ialt,jazi)) &
                + aod_bin(2) * hg(aod_asy(2),elong(ialt,jazi)) &
                + aod_bin(3) * hg(aod_asy(3),elong(ialt,jazi)) 
            mie = brt(aod_ray(ialt,jazi)*airmass_g) * hg2                       

            if(.false.)then ! sum brightness from gas and aerosols

              do ic = 1,nc
!               Rayleigh illumination
                od_per_am = ext_g(ic)
                rayleigh = brtf(airmass_g,od_per_am) * rayleigh_pf(elong(ialt,jazi)) 

!               Total illumination
                clear_rad_c(ic,ialt,jazi) = day_int * (rayleigh + mie)

                if(idebug .ge. 1 .and. ic .eq. 2)then
                  write(6,71)day_int,elong(ialt,jazi),airmass_g,brtf(airmass_g,od_per_am) &
                            ,rayleigh,mie,hg2,clear_rad_c(2,ialt,jazi)
71                format('day_int/elong/am/brt/rayleigh/mie/hg2/clrrd2',f12.0,6f8.3,f12.0)      
                endif

               enddo ! ic

               if(idebug .ge. 1)then
                 write(6,72)day_int,elong(ialt,jazi),mie,hg2,clear_rad_c(:,ialt,jazi)
72               format('day_int/elong/mie/hg2/clear_rad :',f12.0,3f8.3,3f11.0)      
               endif

            else ! assume two scattering layers (g+a, g)
              do ic = 1,nc
                od_g = ext_g(ic)*airmass_g
                od_a = aod_ray(ialt,jazi)*aa
                alpha_g = od_g / 8000.
                alpha_a = od_a / aero_scaleht
                if(od_a .gt. 0)then
                    od_g1 = od_a * (alpha_g / alpha_a)
                else
                    od_g1 = 0.
                endif
                od_g2 = od_g - od_g1

                if(.true.)then ! available light for Mie scattering, relative to Rayleigh scattering
                    am_sun = airmassf(90. - sol_alt,patm) * (od_g2/od_g)
                    sun_trans_g2 = trans(am_sun*ext_g(ic))
                    solar_int_g2 = sun_trans_g2 * 10.**(0.4*sbcorr)
                    solar_int_g2 = min(solar_int_g2,1.0)
                else
                    solar_int_g2 = 1.0
                endif

                pf_eff1 = (rayleigh_pf(elong(ialt,jazi)) * alpha_g + hg2 * alpha_a * solar_int_g2) &
                        / (alpha_g + alpha_a * solar_int_g2)
                od_1 = od_g1 + od_a
                brt1 = brto(od_1) * pf_eff1
                brt2 = brto(od_g2) * rayleigh_pf(elong(ialt,jazi))
                trans1 = trans(od_1)
                clear_rad_c(ic,ialt,jazi) = day_int * ((1.-trans1) * brt1 + trans1 * brt2)

                if(idebug .ge. 1 .and. ic .eq. 2)then
                  write(6,73)day_int,airmass_g,od_g,aod_ray(ialt,jazi),od_a,alpha_g*1e3,alpha_a*1e3,od_g1,od_g2,clear_rad_c(2,ialt,jazi)
73                format('day_int/airmass_g/od_g/aod_ray/od_a/alpha_g/alpha_a/od_g1/od_g2/clear_rad :' &
                        ,f12.0,4f7.3,2x,2f7.3,2x,2f7.3,f11.0)      
                  write(6,*)'am_sun,solar_int_g2 = ',am_sun,solar_int_g2
                endif

             enddo ! ic

            endif

            if(idebug .ge. 1)then
              rmaglim = b_to_maglim(clear_rad_c(2,ialt,jazi))
            endif

          else ! sun below horizon (twilight / night brightness)

!           See http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.67.2567&rep=rep1&type=pdf 
!           patm = 0.85
            z = min(90. - altray,91.)
            airmass_g = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))

            day_int = 3e9
            twi_int = 3e9 ! / 6.

            alt_plane = 90. - abs(sol_alt)
            azi_plane = sol_azi
            xplane = cosd(azi_plane) * cosd(alt_plane)
            yplane = sind(azi_plane) * cosd(alt_plane)
            zplane =                   sind(alt_plane)
            xray   = cosd(view_azi_deg) * cosd(altray)
            yray   = sind(view_azi_deg) * cosd(altray)
            zray   =                      sind(altray)
            angle_r = angleunitvectors(xplane,yplane,zplane&! normal to plane
                                      ,xray,yray,zray)      ! light ray
!           write(6,*)'xplane/yplane/zplane',xplane,yplane,zplane ! test
!           write(6,*)'xray/yray/zray',xray,yray,zray             ! test
            angle_plane = 90. - angle_r ! angle between light ray and plane                              
!           write(6,*)'angle_r/angle_plane = ',angle_r,angle_plane ! test

            horz_dep_r = -sol_alt * rpd                                                  
            dist_pp_plane = horz_dep_r**2 * earth_radius / 2.0 ! approx perpendicular dist
            dist_pp_plane = dist_pp_plane + 13000.             ! shadow enlargement

            if(.true.)then ! get light ray distance to shadow cylinder
                xcos = cosd(angle_r) ! points along the sun's azimuth at 90 degees elong
                zcos = cosd(elong(ialt,jazi))   ! points towards the sun
                ycos_2 = max(1. - xcos**2 - zcos**2,0.)
                ycos = sqrt(ycos_2)             ! perpendicular to xcos and zcos
                x1 = earth_radius - dist_pp_plane
!               y1 = 0.
!               x3 = 0.           
!               y3 = 0.
                dist_ray_plane = linecylp(xcos,ycos,x1,earth_radius) 
            endif

            skyref = .0000001 ! related to airglow / surface lighting?

            if(.true.)then ! assume part of atmosphere is illuminated by the sun
!             dist_ray_plane = dist_pp_plane / sind(angle_plane) ! distance along ray
              bterm = dist_ray_plane**2 / (2. * earth_radius)
              ht_ray_plane = dist_ray_plane * sind(altray) + bterm
              if(ht_ray_plane .le. 99000.)then
                patm_ray_plane = min(ZtoPsa(ht_ray_plane) / 1013.,1.0) * patm
                aero_ray_plane = aod_ray(ialt,jazi)*exp(-ht_ray_plane/aero_scaleht)
                po3_ray_plane = .02 *ay
              else
                patm_ray_plane = 0.
                aero_ray_plane = 0.
              endif

              if(l_solar_eclipse .eqv. .true.)then ! get lat/lon of light ray hitting shadow cylinder
                  call sun_eclipse_parms(i4time,rlat,rlon,htmsl,idebug &
                                        ,altray,view_azi_deg,dist_ray_plane &
                                        ,earth_radius,elgms,emag,eobsc)
                  ecl_int = 1.0 - eobsc
                  if(idebug .ge. 1)then
                      write(6,*)' eclipse elg,emag,eobs ',elgms,emag,eobsc
                  endif
              else
                  ecl_int = 1.0          
              endif

              if(.false.)then ! This appears to become inaccurate near/below the horizon
                  exp_term = 1.0 - 0.5 * (bterm / ht_ray_plane)
                  frac_airmass_unlit = (1.0 - patm_ray_plane)**exp_term
                  frac_airmass_lit = 1.0 - frac_airmass_unlit
                  airmass_lit = airmass_g * frac_airmass_lit      
                  aero_red = 0.
              else
                  altray_plane = altray + cosd(altray)*dist_ray_plane &
                               / (earth_radius*rpd)
                  z = (90. - altray_plane) &
                    + refractd_app(altray_plane,patm)
                  airmass_lit = airmassf(z,patm_ray_plane)
                  z = (90. - altray) + refractd_app(altray      ,patm)
                  airmass_tot = airmassf(z,patm)
                  frac_airmass_lit = airmass_lit / airmass_tot
                  frac_airmass_unlit = 1.0 - frac_airmass_lit
                  airmass_unlit = airmass_tot - airmass_lit

                  aod_path = aod_ray(ialt,jazi) * airmassf(z,1.0)
                  frac_aero_lit = aero_ray_plane / aod_ray(ialt,jazi)
                  aod_lit = aod_path * frac_aero_lit
                  aod_unlit = aod_path * (1.-frac_aero_lit)
!                 arg = aod_lit * 0.5 * (1. - cosd(elong(ialt,jazi)))
!                 Note: is the 40. factor too large?
                  aero_red = (1.0 - exp(-40. * aod_lit)) &       ! aod factor   (0-1)
                           * 0.5 * (1. - cosd(elong(ialt,jazi)))&! elong factor (0-1) 
                           * min(-sol_alt*0.7,1.0)               ! solalt factor(0-1)
              endif

              airmass_to_twi = airmass_unlit
              twi_trans = exp(-.14 * 0.75 * airmass_to_twi)
              do ic = 1,nc
                twi_trans_c(ic) = trans(ext_g(ic)*airmass_to_twi*0.75) 
              enddo ! ic
!             twi_trans_c(:) = trans(ext_g(:) * airmass_to_twi * 0.75) 
!             twi_trans = twi_trans_c(1)

              if(.false.)then ! fractional illumination of sunrise/set value
!               clear_int = max(min(airmass_lit,32.0),.00001)               
                clear_int = &
                   max(min(frac_airmass_lit *twi_trans     ,1.0),.00001)
!               clear_int = &
!                  max(min(     airmass_lit *twi_trans_c(1),1.0),.00001)
!               clear_int = &
!                  max(     brt(airmass_lit)*twi_trans_c(1)     ,.00001)
                clear_intf = clear_int
              else ! experimental absolute illumination (nl)
                rayleigh = rayleigh_pf(elong(ialt,jazi))  ! phase function
                clear_int = &
                  max(twi_int*ecl_int*rayleigh*brt(airmass_lit)*twi_trans_c(1) &
                                                              ,1e3)       
                if(idebug .ge. 2)then
                  write(6,91)elong(ialt,jazi) &
                            ,rayleigh,brt(airmass_lit),twi_trans_c(1),twi_int,clear_int
91                format('elong/rayleigh/brt/twi_trans/day_int/clear_int ',4f9.5,2f11.0)
                endif
                clear_intf = clear_int / twi_int
              endif

            else ! in Earth's shadow (may need secondary scattering)  
              dist_ray_plane = -999.9
              ht_ray_plane = -999.9
              airmass_lit = 0.
              airmass_unlit = 0.
              airmass_tot = 0.
              clear_intf = skyref
              clear_int  = clear_intf * twi_int                           
              aod_lit = 0.
              aod_unlit = 0.
              aero_red = 0.
            endif ! part of atmosphere is illuminated

!           Apply saturation ramp (=1) at high altitudes or low sun
!                                 (=0) at low altitudes and high sun
            sat_sol_ramp = min(-sol_alt / 3.0,1.0) ! 0-1 (lower sun)
            sat_alt_ramp = sqrt(sind(max(altray,0.)))   ! 0-1 (higher alt)
            sat_ramp = 1.0-((1.0 - sat_alt_ramp) * (1.0 - sat_sol_ramp))

!           Ramp2 is zero with either high sun or high alt
            hue_alt_ramp = (sind(altray))**2.0     ! 0-1 (higher alt)
            hue_ramp2 = 1.0 ! sat_sol_ramp * (1.0 - hue_alt_ramp)
            hue_coeff = max(1.0 + min((sol_alt+3.),0.) /  6.,0.2)

!           The value of 'sat_twi_ramp_nt' should also appear in
!           subroutine 'get_clr_rad_nt', variable 'sat_twi_ramp'
            sat_twi_ramp_nt = 0.4

!           clear_int here represents intensity at the zenith
            if(clear_intf .gt. skyref*10.)then      ! mid twilight
                sat_twi_ramp = 1.0
                rint_alt_ramp = 1.0
            elseif(clear_intf .le. skyref)then ! night or late twilight
                sat_twi_ramp = sat_twi_ramp_nt
                am_term = min(sqrt(airmass_g),5.)
                rint_alt_ramp = am_term                  
                clear_intf = skyref
                clear_int  = clear_intf * twi_int                           
            else                              ! late twilight
                frac_twi_ramp = (clear_intf - skyref) / (9.*skyref)
!               sat_twi_ramp = 0.6 + 0.4 * scurve(frac_twi_ramp)
                sat_twi_ramp = sat_twi_ramp_nt &
                       + (1. - sat_twi_ramp_nt) * scurve(frac_twi_ramp)
                am_term = min(sqrt(airmass_g),5.)
                rint_alt_ramp = (1.0       *        frac_twi_ramp) &
                              +  am_term   * (1.0 - frac_twi_ramp)
            endif

!           HSI

!           Hue
!           Hue tends to red with high airmass and blue with low airmass
!           Higher exp coefficient (or low hue) makes it more red
!           Also set to blue with high alt or high sun (near horizon)
!           Aero_red reddens due to aerosols (if value is 1)
!           hue = exp(-airmass_lit*0.75) ! 0:R 1:B 2:G 3:R
            huea = exp(-airmass_unlit*0.12*hue_coeff) ! 0:R 1:B 2:G 3:R
            hue = min(huea,(1.0 - aero_red))
            hue2 = (2.7 - 1.7 * hue ** 1.6) *      hue_ramp2 &
                                      + 1.0 * (1.0-hue_ramp2)
            hue2 = max(hue2,1.5) ! keep aqua color high up

!           Redden further based on aerosols
!           hue2 = min(max(hue2,aero_red/clear_intf),2.7)

            if(hue2 .lt. 2.0)then
                hue2 = 2.0 - sqrt(2.0 - hue2)
            else
                hue2 = 2.0 + .836 * sqrt(hue2 - 2.0)
            endif

            clear_rad_c(1,ialt,jazi) = hue2                          ! Hue

!           Saturation
            if(hue2 .gt. 2.0)then ! Red End                          
                sat_arg = 0.7*abs((hue2-2.0))**1.5
            else                  ! Blue End
                sat_arg = 0.4*abs((hue2-2.0))**1.5
            endif
            clear_rad_c(2,ialt,jazi) = 0.00 + (sat_arg) &            ! Sat
                                           * sat_ramp &                 
                                           * sat_twi_ramp 

!           Intensity
            clear_rad_c(3,ialt,jazi) = clear_int * rint_alt_ramp     ! Int

            if(clear_rad_c(3,ialt,jazi) .lt. .0000099)then
                write(6,*)' WARNING: low value of clear_rad_c' &
                         ,clear_rad_c(3,ialt,jazi),z,airmass_g &
                         ,rint_alt_ramp,clear_int,angle_plane &
                         ,ht_ray_plane,ZtoPsa(ht_ray_plane)
            endif

            if(idebug .ge. 2)then
              write(6,101) airmass_lit,airmass_tot,airmass_unlit&
                          ,hue_coeff,huea,hue  
101           format('airmass_lit/tot/unlit/huec/huea/hue',6f9.5)
              write(6,102)aod_lit,aod_unlit,aero_red &         
                         ,aero_red/clear_intf &                
                         ,clear_intf,skyref & 
                         ,sat_arg,sat_ramp,sat_twi_ramp
102           format('aod_lit/unlit/red/rat/ref/sat',4f9.5,2x,2f9.5,2x,3f9.5)
              rmaglim = b_to_maglim(clear_rad_c(3,ialt,jazi))
            endif

          endif ! sun above horizon

          if(idebug .ge. 1)then 
              if(sol_alt .gt. 0.)then
                write(6,111)ialt,jazi,sol_alt &
                      ,angle_plane,dist_ray_plane,ht_ray_plane &
                      ,clear_rad_c(:,ialt,jazi)/1e9
              else
                write(6,112)ialt,jazi,sol_alt &
                      ,angle_plane,dist_ray_plane,ht_ray_plane &
                      ,patm_ray_plane,airmass_lit,airmass_unlit &
                      ,twi_trans,clear_intf,rint_alt_ramp &
                      ,clear_rad_c(1,ialt,jazi) &
                      ,clear_rad_c(2,ialt,jazi) &
                      ,clear_rad_c(3,ialt,jazi) ! /1e9
              endif
111           format( &
       'ialt/jazi/salt/ang_pln/ds_ray/ht_ray/clrrd' &
                  ,2i4,2x,2f6.2,2f10.1,2x,3f8.3)
112           format( &
       'ialt/jazi/salt/ang_pln/ds_ray/ht_ray/a/t/clrrad' &
!                 ,2i5,2f7.3,2f8.2,2f10.1,2x,4f8.5,f8.0,f8.5,2x,3f8.5) &
                  ,2i4,2x,2f6.2,2f10.1,2x,3f8.4,2x,f8.4,2f8.5 &
                                      ,2x,f7.4,f7.4,f11.0)

              if(idebug .ge. 2)then
                write(6,121)rmaglim
121             format('rmaglim = ',f8.3)
              endif
              write(6,*)

          endif ! idebug .eq. 1

         enddo ! jazi
        enddo ! ialt

        return
        end

        subroutine get_airmass(alt,htmsl,patm &          ! I
                              ,aero_refht,aero_scaleht & ! I
                              ,earth_radius &            ! I
                              ,ag,ao,aa)                 ! O

!       Airmasses relative to zenith at sea level pressure

        include 'trigd.inc'

        include 'rad.inc'

        zapp = 90. - alt

        ztrue = zapp + refractd_app(alt ,patm)

!       If alt < 0 we might calculate htmin 
        if(alt .lt. 0)then
          patm2 = patm
          niter = 2; iter = 1
          do while (iter .le. niter)
            rk = 0.5*(patm+patm2) * 0.13 ! curvature at midpoint of ray
            erad_eff = earth_radius / (1. - rk)
            htmin = htmsl - erad_eff * (alt*rpd)**2 / 2.
            patm2 = ztopsa(htmin) / 1013.25
!           write(6,1)iter,htmin,patm,patm2
1           format('     iter/htmin/patm/patm2 = ',i4,3f10.2)
            iter = iter + 1
          enddo ! while
        endif

!       Gas component for Rayleigh Scattering
        ag = airmassf(ztrue,patm)

!       Ozone component
        ao = airmasso(zapp,htmsl) * patm_o3(htmsl)

!       Aerosol component

!       Note that near the horizon the aerosol airmass should
!       exceed the gas component by the inverse square root of the
!       respective scale heights.

        patm_aero = exp(-((htmsl-aero_refht) / aero_scaleht))
        ZZ = (min(zapp,90.)) * rpd
        aa=1./(COS(ZZ)+.0123*EXP(-24.5*COS(ZZ))) * patm_aero

        return
        end
        