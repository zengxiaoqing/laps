

        subroutine skyglow_phys_twi(ialt_start,ialt_end,ialt_delt &! I
                   ,jazi_start,jazi_end,jazi_delt                 &! I
                   ,minalt,maxalt,minazi,maxazi,idebug_a          &! I
                   ,sol_alt,sol_azi,view_alt,view_az              &! I
                   ,earth_radius,patm,aod_vrt,aod_ray,aod_ray_dir &! I
                   ,aero_scaleht                                  &! I
                   ,htmsl,redp_lvl                                &! I
                   ,aod_ill                                       &! I
                   ,l_solar_eclipse,i4time,rlat,rlon              &! I
                   ,clear_radf_c,horz_dep_d,ag_2d                 &! I
                   ,clear_rad_c,elong                       )      ! O

!       Sky glow with solar altitude <= 0

        include 'trigd.inc'

        use mem_namelist, ONLY: aod_bin, aod_asy

!       Statement Functions
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)
!       brt(a) = (1.0 - exp(-.14 * a))/.14 ! rel. sky brightness (max~7)
        brt(a) = (1.0 - exp(-.14 * a))     ! rel. sky brightness (max=1)
        brtf(am,od_per_am) = 1.0 - exp(-am*od_per_am) ! rel. sky brightness (max=1)
        brto(od) = 1.0 - exp(-min(od,80.)) ! rel. sky brightness (max=1)

        angdif(X,Y)=MOD(X-Y+540.,360.)-180.
        angleunitvectors(a1,a2,a3,b1,b2,b3) = acosd(a1*b1+a2*b2+a3*b3)
        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        real linecyl,linecylp
        linecylp(xcos,ycos,x1,r) = (-(2.*xcos*x1) + sqrt((2.*xcos*x1)**2 - 4. * (xcos**2 + ycos**2) * (x1**2 - r**2))) &
                                 /  (2. * (xcos**2 + ycos**2))

        include 'rad.inc'

        cosp(b,x) = (1. + b/2.) * cosd(x/2)**b

        real twi_trans_c(nc)           ! transmissivity

        real hg2d(nc), alphav_g, alphav_a ! hg2(nc)
        real srcdir_90(nc),srcdir(nc),clear_int_c(nc)
        real ecl_intd(nc),ecl_dir_rat(nc),ecl_scat(nc)

        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)  ! clear sky illumination
        real clear_rad_c0(nc,minalt:maxalt,minazi:maxazi) ! clear sky (no eclipse)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)! integrated
               ! fraction of air molecules illuminated by the sun along
               ! line of sight (consider Earth's shadow + clouds)
               ! topo and multiple scattering is accounted for
        real ag_2d(minalt:maxalt,minazi:maxazi) ! gas airmass (topo/notopo)
        real aod_ill(minalt:maxalt,minazi:maxazi) ! aerosol illuminated
                                      ! optical depth (slant - topo/notopo)
        real aod_ray(minalt:maxalt,minazi:maxazi) ! aerosol optical 
                                      ! depth (zenithal) may be adjusted
                                      ! for illumination
        real aod_ray_dir(minalt:maxalt,minazi:maxazi) ! aerosol optical 
                                      ! depth (zenithal) may be adjusted
                                      ! for direct illumination / topo
        real patm_ray(minalt:maxalt,minazi:maxazi)! effective patm adjusted
                                      ! for illumination
        real elong(minalt:maxalt,minazi:maxazi)
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

        parameter (nsc = 3)
!       real aod_asy_eff(nsc) ! ,nc (add nc dimension throughout)

        logical l_solar_eclipse

        parameter (nsteps = 100000)
        parameter (nopac = 10)

        real distecl(nopac,nc)
        real tausum_a(nsteps)
        real taumid(nopac),opacmid(nopac),eobsc_a(nopac,nc)
        real eobsc(nc),eobsc_sum(nc)
        real sky_rad_scat(nc,minalt:maxalt,minazi:maxazi)            
        real sky_rad_scata(minazi:maxazi)
        real ags_a(-5:+5), aas_a(-5:+5)

        sfc_alb = 0.15 ! pass this in and account for snow cover?

        aod_bin(1) = .000
        aod_bin(2) = .987
        aod_bin(3) = .013

        aod_asy(1) =  .75
        aod_asy(2) =  .70
        aod_asy(3) = -.65

        patm_ray = patm

        aero_refht = redp_lvl

        if(sol_alt .gt. 0.)then
            write(6,*)' skyglow_phys: i4time is ',i4time,l_solar_eclipse
            write(6,*)' aod_bin = ',aod_bin
            write(6,*)' aod_asy = ',aod_asy
            write(6,*)' htmsl / aero_refht = ',htmsl,aero_refht
            write(6,*)' aod_vrt = ',aod_vrt
            write(6,*)' aod_ray max = ',maxval(aod_ray)
!           write(6,*)' patm_ray max = ',maxval(patm_ray)
        endif

!       Obtain reference values of source term
!       call get_airmass(90.,htmsl,patm &           ! I
!                        ,aero_refht,aero_scaleht & ! I
!                        ,earth_radius,1 &          ! I
!                        ,ag_90,ao_90,aa_90)        ! O

!       Obtain reference values for sun             
!       call get_airmass(sol_alt,htmsl,patm &       ! I
!                        ,aero_refht,aero_scaleht & ! I
!                        ,earth_radius,1 &          ! I
!                        ,ag_s,ao_s,aa_s)           ! O

        do ialt = ialt_start,ialt_end,ialt_delt
  
         altray = view_alt(ialt,jazi_start)
         if(altray .eq. nint(altray))then
           iverbose = 1
         else
           iverbose = 0
         endif
         call get_airmass(altray,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,iverbose &   ! I
                         ,ag,ao,aa)                 ! O

!        Determine aerosol multiple scattering order
!        altscat = 1.00 * altray + 0.00 * sol_alt
!        altscat = max(altray,sol_alt)
!        altscat = sqrt(0.5 * (altray**2 + solalt**2))
!        call get_airmass(altscat,htmsl,patm &      ! I
!                        ,aero_refht,aero_scaleht & ! I
!                        ,earth_radius,1 &          ! I
!                        ,agdum,aodum,aascat)       ! O
!        scatter_order = 1.0 ! f(altray,sol_alt)
!        scatter_order = max(aod_vrt*aascat,1.0)**1.0 ! 0.5-1.5
!        write(6,*)' aod_ray/aascat/sco',aod_ray(ialt,minazi),aascat,scatter_order

!        Determine effective asymmetry parameter from multiple scattering
!        http://www.arm.gov/publications/proceedings/conf15/extended_abs/sakerin_sm.pdf
!        do isc = 1,nsc
!          if(aod_asy(isc) .gt. 0.)then
!             aod_asy_eff(isc) = aod_asy(isc) ** scatter_order
!          elseif(aod_asy(isc) .lt. 0.)then
!             aod_asy_eff(isc) = -(abs(aod_asy(isc)) ** scatter_order)
!          else
!             aod_asy_eff(isc) = 0. 
!          endif
!        enddo ! isc

!        cosp_frac = max(1.-scatter_order,0.)

!        Set this for values of -5 to +5 degrees of solar altitude
         if(.false.)then
           ags_a = ag_s/ag_90
           aas_a = aa_s/aa_90
         elseif(.false.)then
           do isolalt = -5,+5
             sol_alt_a = sol_alt + float(isolalt)
             call get_airmass(sol_alt_a,htmsl,patm &    ! I
                             ,aero_refht,aero_scaleht & ! I
                             ,earth_radius,1 &          ! I
                             ,ag_s,ao_s,aa_s)           ! O
             ags_a(isolalt) = ag_s / ag_90
             aas_a(isolalt) = aa_s / aa_90
           enddo 
         endif

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

          idebug = idebug_a(ialt,jazi) * 2

          if(sol_alt .gt. 0. .or. altray .lt. -horz_dep_d)then
            continue

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
            dist_pp_plane = (horz_dep_r**2 * earth_radius / 2.0) - htmsl ! approx perpendicular dist
            dist_pp_plane_s = dist_pp_plane

!           Enlarge the shadow except when near the sun and near sunrise
!           shadow_enlarge_high = 13000. * (1.0 - cosd(elong(ialt,jazi)))/2
            if(htmsl .le. 14600.)then
                elong_arg = min(elong(ialt,jazi),90.)
                shadow_enlarge_high = (14600.-htmsl) * sind(elong_arg)
                shadow_enlarge_low  =  14600.-htmsl                     
                ramp_enl = min(max(1.0 - (-2. - sol_alt) / 2.0,0.),1.)
                shadow_enlarge = shadow_enlarge_high * ramp_enl + shadow_enlarge_low * (1. - ramp_enl)
            else
                shadow_enlarge = 0.
            endif
            dist_pp_plane = dist_pp_plane + shadow_enlarge ! shadow enlargement
!           dist_pp_plane = dist_pp_plane + 13000.

            if(.true.)then ! distance along light ray to shadow cylinder
                xcos = cosd(angle_r) ! points along sun's azimuth at 90 deg elong
                zcos = cosd(elong(ialt,jazi)) ! points towards the sun
                ycos_2 = max(1. - xcos**2 - zcos**2,0.)
                ycos = sqrt(ycos_2)        ! perpendicular to xcos and zcos
                x1 = earth_radius - dist_pp_plane
                x1_s = earth_radius - dist_pp_plane_s
                dist_ray_plane   = linecylp(xcos,ycos,x1,  earth_radius)
                dist_ray_plane_s = linecylp(xcos,ycos,x1_s,earth_radius) 
            endif

            skyref = .000000003 ! smaller values reduce twilight artifacts
                                ! setting also related to airglow?

            if(.true.)then ! assume part of atmosphere is illuminated by the sun
              bterm = dist_ray_plane**2 / (2. * earth_radius)
              ht_ray_plane = (dist_ray_plane * sind(altray) + bterm) + htmsl

              if(ht_ray_plane .le. 99000.)then
!               Pressure at ray plane intersection in standard atmospheres?
                patm_ray_plane = min(ZtoPsa(ht_ray_plane)/1013.,1.0) !*patm
              else
!               Pressure at ray plane intersection in standard atmospheres?
                patm_ray_plane = min(ZtoPsa(99000.) / 1013.,1.0) !*patm
                scalehts = min((ht_ray_plane - 99000.) / 8000.,10.)
                patm_ray_plane = patm_ray_plane * exp(-scalehts)
              endif

              bterm = dist_ray_plane_s**2 / (2. * earth_radius)
              ht_ray_plane_s = (dist_ray_plane_s * sind(altray) + bterm) &
                             + htmsl

              if(ht_ray_plane_s .le. 99000.)then
                aero_ray_plane = aod_vrt * exp(-ht_ray_plane_s/aero_scaleht)
              else
                aero_ray_plane = 0.
              endif

              if(l_solar_eclipse .eqv. .true.)then ! get lat/lon of light ray hitting shadow cylinder
                  call sun_eclipse_parms(i4time,rlat,rlon,htmsl,idebug &
                                    ,altray,view_azi_deg,dist_ray_plane &
                                    ,earth_radius,elgms,emag,eobscf,eobsc(:))
                  ecl_int = 1.0 - eobsc(2)
                  if(idebug .ge. 1)then
                      write(6,*)' eclipse elg,emag,eobs ',elgms,emag,eobsc(2)
                  endif
              else
                  ecl_int = 1.0          
              endif

              if(.true.)then 
                  altray_plane = altray + cosd(altray)*dist_ray_plane &
                               / (earth_radius*rpd)
                  z = (90. - altray_plane) &
                    + refractd_app(altray_plane,patm)
                  airmass_lit = airmassf(z,patm_ray_plane)

                  za = (90. - altray) + refractd_app(altray      ,patm)
!                 airmass_tot = airmassf(za,patm)
                  airmass_tot = ag

                  frac_airmass_lit = airmass_lit / airmass_tot
                  frac_airmass_unlit = 1.0 - frac_airmass_lit
                  airmass_unlit = airmass_tot - airmass_lit

!                 aod_path = aod_ray(ialt,jazi) * airmassf(z,1.0)
                  aod_path = aod_ray(ialt,jazi)*aa ! *ext_a(2)
                  if(aod_ray(ialt,jazi) .gt. 0.)then
                      frac_aero_lit = aero_ray_plane / aod_ray(ialt,jazi)
                  else
                      frac_aero_lit = 0.
                  endif
                  aod_lit = aod_path * frac_aero_lit
                  aod_unlit = aod_path * (1.-frac_aero_lit)
!                 arg = aod_lit * 0.5 * (1. - cosd(elong(ialt,jazi)))
!                 Note: is the 40. factor too large?
!                 aero_red = (1.0 - exp(-40. * aod_lit)) &  ! aod fct (0-1)
                  aero_red = (1.0 - exp(-240. * aod_lit)) &  ! aod fct (0-1)
                           * 0.5 * (1. - cosd(elong(ialt,jazi)))&! elong factor (0-1) 
                           * min(-sol_alt*0.7,1.0) &     ! solalt fct (0-1)
                           * cosd(altray)**60.           ! alt fact   (0-1)
              endif

              airmass_to_twi = airmass_unlit
              do ic = 1,nc
                twi_od_eff = ext_g(ic)*airmass_to_twi + aod_unlit*0.75
                twi_trans_c(ic) = trans(twi_od_eff)
              enddo ! ic

              rayleigh_gnd = rayleigh_pf(elong(ialt,jazi))
              rayleigh = rayleigh_gnd  ! phase function

              mode_twi = 1
              if(mode_twi .eq. 2)then ! rad twilight
                continue

              else ! absolute illumination (nl), HSI
                if(airmass_lit .gt. 1.e-5)then
                    brta = brt(airmass_lit)

                else
                    brta = airmass_lit * .14
                endif
                clear_int = &
                  max(twi_int*ecl_int*rayleigh*brta*twi_trans_c(1) &
                                                              ,1e2) ! floor
                if(idebug .ge. 2)then
                  write(6,91)elong(ialt,jazi) &
                            ,rayleigh,brta,twi_trans_c(1),twi_int,clear_int
91                format('elong/rayleigh/brt/twi_trans/day_int/clear_int ',2f9.5,f12.9,f9.5,2f12.0)
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
!           if(airmass_lit .lt. .02)then ! Earth shadow
!               sat_ramp = sat_ramp * (airmass_lit / .02)
!           endif

!           Ramp2 is zero with either high sun or high alt
            hue_alt_ramp = (sind(altray))**2.0     ! 0-1 (higher alt)
!           hue_ramp2 = 1.0 ! sat_sol_ramp * (1.0 - hue_alt_ramp)
            hue_coeff = 0.2 ! max(1.0 + min((sol_alt+3.),0.) /  6.,0.2)

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

            if(mode_twi .eq. 1)then ! HSI twilight

!             HSI

!             Hue
!             Hue tends to red with high airmass and blue with low airmass
!             Higher exp coefficient (or low hue) makes it more red
!             Also set to blue with high alt or high sun (near horizon)
!             Aero_red reddens due to aerosols (if value is 1)
!             Increase huea coefficient to redden the horizon?
!             hue = exp(-airmass_lit*0.75) ! 0:R 1:B 2:G 3:R
              huea = exp(-(airmass_unlit-4.0)*0.40*hue_coeff) 
              hue = min(huea,(1.0 - aero_red)) ! set hue to 0 via aero_red
              hue2 = (2.8 - 1.8 * hue ** 1.6)  ! 0:R 1:B 2:G 3:R
              hue2 = max(hue2,1.5) ! keep aqua color high up

!             Redden further based on aerosols
!             hue2 = min(max(hue2,aero_red/clear_intf),2.7)

              if(hue2 .lt. 2.0)then
                  hue2 = 2.0 - sqrt(2.0 - hue2)
              else
                  hue2 = 2.0 + .836 * sqrt(hue2 - 2.0)
              endif

              clear_rad_c(1,ialt,jazi) = hue2                        ! Hue

!             Saturation
              if(hue2 .gt. 2.0)then ! Red End                          
                  sat_arg = 0.7*abs((hue2-2.0))**2.0
              else                  ! Blue End
                  sat_arg = 0.4*abs((hue2-2.0))**2.0
              endif
              sat_aod_ramp = 1.2 - (aod_vrt * 2.)
              clear_rad_c(2,ialt,jazi) = 0.10 + (sat_arg) &          ! Sat
                                           * sat_ramp &                 
                                           * sat_twi_ramp &
                                           * sat_aod_ramp

!             Intensity
              clear_rad_c(3,ialt,jazi) = clear_int * rint_alt_ramp   ! Int

              if(clear_rad_c(3,ialt,jazi) .lt. .0000099)then
                  write(6,*)' WARNING: low value of clear_rad_c' &
                           ,clear_rad_c(3,ialt,jazi),za,airmass_g &
                           ,rint_alt_ramp,clear_int,angle_plane &
                           ,ht_ray_plane,ZtoPsa(ht_ray_plane)
              endif

            else ! rad twilight (experimental)
              do ic = 1,nc
                  clear_rad_c(ic,ialt,jazi) = clear_int_c(ic)
              enddo ! ic                
            endif          

            huecall = clear_rad_c(1,ialt,jazi)
            satcall = min(clear_rad_c(2,ialt,jazi)*2.0, 1.0)
            rincall = clear_rad_c(3,ialt,jazi) * 0.255
            call hsl_to_rgb(huecall,satcall,rincall,  &
                            clear_rad_c(1,ialt,jazi), &
                            clear_rad_c(2,ialt,jazi), &
                            clear_rad_c(3,ialt,jazi) ) 

            if(idebug .ge. 2)then
              write(6,101) airmass_lit,airmass_tot,airmass_unlit&
                          ,hue_coeff,huea,hue,hue2
101           format('airmass_lit/tot/unlit/huec/huea/hue/hue2',f12.8,6f9.5)
              write(6,102)aa,ht_ray_plane_s,aero_ray_plane,aod_lit &         
                         ,aero_red,aero_red/clear_intf &                
                         ,clear_intf,skyref & 
                         ,sat_arg,sat_ramp,sat_twi_ramp
102           format('aa/htryplns/plane/aod_lit/red/rat/cif/skr/sat',f9.5,f10.0,3f9.5,f10.3,f10.3,f14.10,2x,3f9.5)
              rmaglim = b_to_maglim(clear_rad_c(3,ialt,jazi))
            endif

          endif ! sun above horizon

          if(idebug .ge. 1)then 
              if(sol_alt .gt. 0.)then
                write(6,111)ialt,jazi,sol_alt &
                      ,od_a,clear_rad_c(:,ialt,jazi)/1e9 
              else
                write(6,112)altray,view_azi_deg,sol_alt &
                      ,angle_plane,dist_pp_plane,dist_ray_plane,ht_ray_plane,shadow_enlarge &
                      ,patm_ray_plane,airmass_lit,airmass_unlit &
                      ,twi_trans_c(1),clear_intf,rint_alt_ramp &
                      ,clear_rad_c(1,ialt,jazi) &
                      ,clear_rad_c(2,ialt,jazi) &
                      ,clear_rad_c(3,ialt,jazi) ! /1e9
              endif
111           format( &
       'ialt/jazi/salt/od_a/clrrd' &
                  ,i4,i5,2x,2f6.2,2f10.1,2x,f8.3)
112           format( &
       'alt/azi/salt/ang_pln/ds_pp/ds_ray/ht_ray/enl/a/t/clrrad' &
                  ,2f6.1,1x,2f6.2,3f11.1,f9.0,2x,3f8.4,2x,f8.4,2f8.5 &
                                      ,2x,3f12.0)

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
