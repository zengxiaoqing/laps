
        subroutine get_sky_rgb(r_cloud_3d,cloud_od,cloud_od_sp,nsp, &
                               r_cloud_rad,cloud_rad_c, &
                               clear_rad_c, &
                               clear_radf_c,patm,htmsl, &               ! I
                               glow,glow_sun,glow_moon,glow_stars, &
                               od_atm_a,transm_obs,isun,jsun, &         ! I
                               airmass_2_cloud,airmass_2_topo, &
                               topo_swi,topo_albedo,albedo_sfc, &       ! I
                               aod_2_cloud,aod_2_topo,aod_ill, &
                               alt_a,azi_a,elong_a,ni,nj,sol_alt,sol_az, &
                               moon_alt,moon_az,moon_mag, &
                               sky_rgb)                                 ! O

        use mem_namelist, ONLY: r_missing_data,earth_radius,aero_scaleht,redp_lvl
        include 'trigd.inc'

!       Statement functions
        addlogs(x,y) = log10(10.**x + 10.**y)
        trans(od) = exp(-od)
        opac(od) = 1.0 - exp(-od)
        brt(a,ext) = 1.0 - exp(-a*ext)
        rad_to_counts(rad) = (log10(rad)-7.3)*100.
        counts_to_rad(counts) = 10.**((counts/100.)+7.3)

        include 'rad.inc'

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real cloud_od(ni,nj)        ! cloud optical depth
        real cloud_od_sp(ni,nj,nsp) ! cloud species optical depth
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,ni,nj)  ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,ni,nj)  ! clear sky illumination
                                    ! local/input when sun is above/below horizon
        real moon_rad_c(nc,ni,nj)   ! clear sky illumination from moon
        real clear_radf_c(nc,ni,nj) ! integrated fraction of air illuminated 
                                    ! by the sun along line of sight  
                                    ! (accounting for Earth shadow + clouds)
        real clear_rad_c_nt(3)      ! HSV night sky brightness
        real glow(ni,nj)            ! skyglow (log b in nl, sun or moon from 'vi')           
        real glow_sun(ni,nj)        ! sunglow (log b in nl, extended obj)
        real glow_moon(ni,nj)       ! moonglow (log b in nl, extended obj)
        real glow_moon_sc(ni,nj)    ! moonglow (log b in nl, scattered)
        real glow_stars(nc,ni,nj)   ! starglow (log b in nl)           
        real airmass_2_cloud(ni,nj) ! airmass to cloud 
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        real topo_swi(ni,nj)        ! terrain illumination (relative to clr sky)
        real topo_albedo(nc,ni,nj)  ! terrain albedo (* sin solar alt)
        real aod_2_cloud(ni,nj)     ! future use
        real aod_2_topo(ni,nj)      ! aerosol optical depth to topo
        real aod_ill(ni,nj)         ! aerosol illuminated optical depth 
        real od_atm_a_eff(ni,nj)    ! aerosol illuminated tau per airmass
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real elong_a(ni,nj)
        real idebug_a(ni,nj)
        real rintensity(nc), cld_rgb_rat(nc), glow_cld_c(nc)
        real pf_scat(ni,nj), pf_scat1(ni,nj), pf_scat2(ni,nj)
        real bkscat_alb(ni,nj)
       
        real trans_c(nc)            ! transmissivity

        real sky_rgb(0:2,ni,nj)
        real moon_alt,moon_az,moon_mag

        new_skyglow = 1

        write(6,*)' get_sky_rgb: sol_alt = ',sol_alt
        write(6,*)' moon alt/az/mag = ',moon_alt,moon_az,moon_mag

        idebug_a = 0

!       htmsl = psatoz(patm*1013.25)

        if(sol_alt .le. 0.)then

!         Use max sky brightness to calculate a secondary glow (twilight)
          call get_twi_glow_ave(            log10(clear_rad_c(3,:,:)) &
                     ,alt_a,azi_a,ni,nj,sol_alt,sol_az,twi_glow_ave)
          arg = maxval(log10(clear_rad_c(3,:,:))) ! avoid moon

          write(6,*)' twi_glow_ave = ',twi_glow_ave
          write(6,*)' twi_glow_max = ',arg

!         glow_secondary = arg - 1.5  ! one thirtieth of max sky brightness
          glow_diff_cld      = 0.8
          glow_diff_clr      = 1.2 ! 1.5
          glow_secondary_cld = twi_glow_ave - glow_diff_cld
          glow_secondary_clr = twi_glow_ave - glow_diff_clr

          write(6,*)' glow_secondary_cld = ',glow_secondary_cld
          write(6,*)' glow_secondary_clr = ',glow_secondary_clr

          sol_alt_eff = max(sol_alt,-16.)
!         Test at -3.2 solar alt
          fracalt = sol_alt_eff / -16.
          corr1 = 8.6; corr2 = 3.612
          corr1 = 8.6; corr2 = 3.412
!         argref = corr1 + (fracalt**0.83 * (corr2-corr1))
          argref = corr1 + (fracalt**0.68 * (corr2-corr1))
          contrast = 70. + 30. * (abs(sol_alt_eff + 8.)/8.)**2
          write(6,*)' argref = ',argref
          write(6,*)' contrast = ',contrast
        endif

!       Sun/Cloud glow Calculation
!       http://www.astro.umd.edu/~ssm/ASTR620/mags.html#solarabsmag
!       Sun = 96000 lux = 1367 w/m**2 = -26.81 bolometric mag (above the atmosphere)
!       Vega = 2.54 * 1e-6 lux
!       1 lambert = 1 lumen / cm**2 = 1e4/pi candela/m**2
!       Sun is -26.74 mag per sphere (UBVRI: -25.96, -26.09, -26.74, -27.26, -27.55)
!       (180/pi)^2 * 4 * pi = 41250 square degrees in a sphere.
!       5.346e11 square arcsec in a sphere
        glow_cld_day = v_to_b(-26.74 + log10(5.346e11)*2.5)
        write(6,*)' glow_cld_day (nl) = ',glow_cld_day

!       Redness section (sun at low solar altitude)
!       Compare with 'get_cloud_rad' redness calculation?
        sol_alt_red_thr = 9.0 + (od_atm_a * 20.)
        if(sol_alt .le. sol_alt_red_thr .and. sol_alt .gt. -16.0)then
!           Add aerosols to trans calculation (call 'get_airmass')
            am = airmassf(90.-sol_alt,patm)
            do ic = 1,nc
              trans_c(ic) = trans(am*ext_g(ic)*patm)
            enddo
            rob_sun = trans_c(1) / trans_c(3)
            rog_sun = trans_c(1) / trans_c(2)
        else
            rob_sun = 1.
            rog_sun = 1.
            am = -99.9
        endif
        write(6,5)sol_alt_red_thr,od_atm_a,am,rob_sun,rog_sun
5       format('  sol_alt_red_thr/od_atm_a/am/rob_sun/rog_sun = ',5f9.3)

!       Grnness section (clear sky at low solar altitudes)      
        sol_alt_grn_thr = 10.0                       
        if(sol_alt .le. sol_alt_grn_thr .and. sol_alt .gt. -16.0)then
            grnness = min((sol_alt_grn_thr - sol_alt)/sol_alt_grn_thr,1.0)
        else
            grnness = 0.
        endif
        write(6,*)' sol_alt_grn_thr/grnness = ',sol_alt_grn_thr,grnness

!       Brighten resulting image at night
        if(sol_alt .lt. -16.)then
            ramp_night = 1.0
        else
            ramp_night = 1.0
        endif

!       azid1 = 46. ; azid2 = 226.
        azid1 = 90. ; azid2 = 270.
        if(sol_alt .gt. 0.)then
            azid1 = nint(sol_az)
            azid2 = mod(azid1+180.,360.)
        elseif(moon_alt .gt. 0.)then
            azid1 = nint(moon_az)
            azid2 = mod(azid1+180.,360.)
        endif

        do j = 1,nj
        do i = 1,ni

         if(.true.)then                           
          if(.false.)then ! polar
!             if(i .eq. ni/2 .AND. j .eq. (j/5) * 5)then
              if(i .eq.    j .AND. j .eq. (j/5) * 5)then ! SW/NE
                  idebug_a(i,j) = 1
              endif
          else ! cyl
              if(azi_a(i,j) .eq. azid1 .OR. azi_a(i,j) .eq. azid2)then ! constant azimuth
                  idebug_a(i,j) = 1
              endif
          endif
         endif

        enddo ! i
        enddo ! j

        if(sol_alt .gt. 0.)then

            I4_elapsed = ishow_timer()

            write(6,*)' isun,jsun=',isun,jsun

!           Determine if observer is in terrain shadow
            if(airmass_2_topo(isun,jsun) .gt. 0.)then 
                write(6,*)' Sun is behind terrain'
                od_atm_a_eff = 0.
                sun_vis = 0.
            else ! sun is outside of terrain
                if(transm_obs .gt. 0.9)then
                    od_atm_a_eff = od_atm_a
                    sun_vis = 1.
                else ! in cloud/terrain shadow
                    od_atm_a_eff = 0.
                    sun_vis = 0.
                endif
            endif

!           Calculate normalized illuminated aerosol optical depth (zenith)
            iprint = 0
            do i = ni,1,-1
              call get_airmass(alt_a(i,1),htmsl,patm &   ! I
                              ,redp_lvl,aero_scaleht &   ! I
                              ,earth_radius &            ! I
                              ,ag,ao,aa)                 ! O 
              do j = 1,nj
!               Crepuscular ray experiments
!               if(airmass_2_topo(i,j) .gt. 0. .OR. transm_obs .lt. 0.9)then
                if(airmass_2_topo(i,j) .gt. 0.)then ! inside terrain 
                  od_atm_a_eff(i,j) = od_atm_a * aod_ill(i,j) / aod_2_topo(i,j)
                  if(idebug_a(i,j) .eq. 1)then
                    od_atm_test = od_atm_a * aod_ill(i,j) / aod_2_topo(i,j)
                    write(6,6)alt_a(i,j),od_atm_a,aod_2_topo(i,j) &
                             ,aod_ill(i,j),od_atm_test     
6                   format(' alt/od_atm_a/aod2topo/aod_ill/od_atm_test',2f9.3,2f10.5,f9.3)
                  endif

                else ! outside of terrain
                  od_atm_a_eff(i,j) = (od_atm_a * aod_ill(i,j)) / (aa * od_atm_a)
                  if((i .eq. isun .AND. j .eq. jsun) .OR. &
                     (clear_radf_c(2,i,j) .lt. 0.5 .AND. iprint .le. 30) .OR. &
                     (cloud_od(i,j) .eq. 0. .AND. &
                      od_atm_a_eff(i,j)/od_atm_a .le. 0.5 .AND. &
                      alt_a(i,j) .gt. 3.0                 .AND. &
                      iprint .le. 30)                              )then
                      if(i .eq. isun .AND. j .eq. jsun)then   
                          write(6,*)'Calcs at solar location outside of terrain:'
                      elseif(clear_radf_c(2,i,j) .lt. 0.5)then
                          write(6,*)'Clear with cloud shadowing outside of terrain:'
                      else
                          write(6,*)'Clear with aerosol shadowing outside of terrain:'
                      endif
                      write(6,7)alt_a(i,j),azi_a(i,j),od_atm_a,aod_ill(i,j),aa &
                               ,od_atm_a_eff(i,j),clear_radf_c(2,i,j)
7                     format(' altaz/od_atm_a/aod_ill/aa/od_atm_a_eff/clr_radf =   ',2f9.2,3f9.3,2f9.3)
                      iprint = iprint + 1
                  endif
                endif

              enddo ! j
            enddo ! i

            write(6,8)transm_obs,od_atm_a,od_atm_a_eff(isun,jsun),sun_vis
8           format(' transm_obs,od_atm_a,od_atm_a_eff,sun_vis=',4f8.3)

            write(6,*)' od_atm_a_eff range = ',minval(od_atm_a_eff),maxval(od_atm_a_eff)

            write(6,*)' call skyglow_phys:'

            call skyglow_phys(1,ni,1 &                               ! I
                     ,1,nj,1 &                                       ! I
                     ,1,ni,1,nj,idebug_a &                           ! I
                     ,sol_alt,sol_az,alt_a,azi_a &                   ! I
                     ,earth_radius,patm,od_atm_a_eff,aero_scaleht &  ! I
                     ,htmsl,redp_lvl &                               ! I
                     ,.false.,i4time,rlat,rlon &                     ! I
                     ,clear_radf_c &                                 ! I
                     ,clear_rad_c,elong_a                     )      ! O

            write(6,*)' range of clear_rad_c(2) is ',minval(clear_rad_c(2,:,:)),maxval(clear_rad_c(2,:,:))

        elseif(sol_alt .lt. -16. .and. moon_alt .gt. 0.)then ! sun below -16. and moon is up
            od_atm_a_eff = od_atm_a

            write(6,*)' range of clear_radf_c is ',minval(clear_radf_c(2,:,:)),maxval(clear_radf_c(2,:,:))

            write(6,*)' call skyglow_phys for moon testing:'
            call skyglow_phys(1,ni,1 &                               ! I
                     ,1,nj,1 &                                       ! I
                     ,1,ni,1,nj,idebug_a &                           ! I
                     ,moon_alt,moon_az,alt_a,azi_a &                 ! I
                     ,earth_radius,patm,od_atm_a_eff,aero_scaleht &  ! I
                     ,htmsl,redp_lvl &                               ! I
                     ,.false.,i4time,rlat,rlon &                     ! I
                     ,clear_radf_c &                                 ! I
                     ,moon_rad_c,elong_a                     )       ! O

            glow_moon_sc(:,:) = log10(moon_rad_c(2,:,:)) + (-26.7 - moon_mag) * 0.4
            glow(:,:) = glow_moon_sc(:,:) ! experimental

            write(6,*)' range of moon_rad_c(2) is ',minval(moon_rad_c(2,:,:)),maxval(moon_rad_c(2,:,:))
            write(6,*)' range of glow_moon_sc is ',minval(glow_moon_sc(:,:)),maxval(glow_moon_sc(:,:))
            write(6,*)' range of glow is ',minval(glow(:,:)),maxval(glow(:,:))

        endif

        write(6,*)' range of glow_sun is ',minval(glow_sun),maxval(glow_sun)
        write(6,*)' range of glow_moon is ',minval(glow_moon),maxval(glow_moon)

        I4_elapsed = ishow_timer()

        if(ni .eq. nj)then ! polar
            write(6,*)' slice from SW to NE through midpoint'
        else
            write(6,*)' slices at degrees azimuth: ',azid1,azid2
        endif

        if(.true.)then
          if(sol_alt .ge. 0.)then       ! daylight
            write(6,11)
11          format('    i   j   alt   azi  elong   pf_scat   opac    od /        species            alb    cloud airmass   rad   ', &
                   'rinten   airtopo  switopo topoalb aodill  topood topovis cld_visb   glow      skyrgb')
          elseif(sol_alt .ge. -16.)then ! twilight
            write(6,12)
12          format('    i   j      alt      azi     elong   pf_scat   opac       od      alb     cloud  airmass   rad    ', &
                   'rinten glw_cld_nt glw_cld  glw_twi glw2ndclr rmaglim  cld_visb  glow      skyrgb')
          else                          ! night
            write(6,13)
13          format('    i   j      alt      azi     elong   pf_scat   opac       od      alb     cloud  airmass   rade-3 ', &
                   'rinten glw_cld_nt glwcldmn glw_cld glw2ndclr rmaglim  cld_visb  glow      skyrgb')
          endif
        endif

        call get_cld_pf(elong_a,r_cloud_rad,cloud_od,cloud_od_sp,airmass_2_topo,ni,nj & ! I
                       ,pf_scat1,pf_scat2,pf_scat,bkscat_alb) ! O

        do j = 1,nj
        do i = 1,ni

          sky_rgb(:,i,j) = 0.          
          if(azi_a(i,j) .eq. azid1 .OR. azi_a(i,j) .eq. azid2)then ! constant azimuth
              if(i .eq. 1)write(6,*)   
          endif

          idebug = idebug_a(i,j)

!         Obtain cloud brightness
          if(sol_alt .ge. -4.)then ! Day/twilight from cloud_rad_c array
!             Potential intensity of cloud if it is opaque 
!               (240. is nominal intensity of a white cloud far from the sun)
!               (0.25 is dark cloud base value)                                  
              rint_top  = 240.                                * pf_scat(i,j) 
              rint_base = 240. * (  0.3 * (1. + albedo_sfc) ) * pf_scat(i,j) 

!             Gamma color correction applied when sun is near horizon 
              if(cloud_rad_c(1,i,j) .gt. 0.)then
                  cld_rgb_rat(:) = (cloud_rad_c(:,i,j)/cloud_rad_c(1,i,j))**0.35
              else
                  cld_rgb_rat(:) = 1.0
              endif

!             rintensity = min(rintensity,255.)
              rintensity(:) = rint_top  * cld_rgb_rat(:) * r_cloud_rad(i,j) &
                            + rint_base * (1.0 - r_cloud_rad(i,j))
              rintensity = max(rintensity(:),0.)

              if(idebug_a(i,j) .eq. 1 .AND. abs(alt_a(i,j)) .le. 2.0)then
!                 write(6,41)pf_scat1(i,j),elong_a(i,j),pf_snow,snow_factor,pf_scat(i,j),rint_base,r_cloud_rad(i,j),(rint_top  * cld_rgb_rat(1) * r_cloud_rad(i,j)),rint_base * (1.0 - r_cloud_rad(i,j)),rintensity(1)
                  write(6,41)pf_scat1(i,j),elong_a(i,j),cloud_od_snow,cloud_od_tot,snow_bin1,pf_snow,snow_factor,pf_scat2(i,j),pf_scat(i,j),rintensity(1)
 41               format(' pf/rintensity(1) = ',10f9.3)
              endif

!             Apply cloud reddening
              do ic = 1,nc
                  trans_c(ic) = trans(ext_g(ic) * airmass_2_cloud(i,j))              
              enddo
              rintensity(:) = rintensity(:) * (trans_c(:)/trans_c(1))**0.25 ! 0.45

          else ! later twilight (clear_rad_c) and nighttime (surface lighting)
              glow_cld = glow_secondary_cld                                     

              if(sol_alt .lt. -16.)then
!                 Add phys term for scattering by moonlight on the clouds
                  pf_scat_moon = pf_scat(i,j)**1. ! empirically stronger due to linear scale?             
                  glow_cld_moon = log10(glow_cld_day * cloud_rad_c(2,i,j) * pf_scat_moon)
                  glow_cld = addlogs(glow_cld,glow_cld_moon)
              else
                  glow_cld_moon = -999.
              endif

!             Note that cloud_rad_c(1) is surface night lighting of clouds (nl)
              glow_cld_nt = log10(cloud_rad_c(1,i,j))
              glow_cld_c(1) = addlogs(glow_cld,glow_cld_nt+0.20)
              glow_cld_c(2) = addlogs(glow_cld,glow_cld_nt+0.00)
              glow_cld_c(3) = addlogs(glow_cld,glow_cld_nt-0.20)

!             During twilight, compare with clear sky background
!             Note that secondary scattering might also be considered in
!             early twilight away from the bright twilight arch.
!             The result creates a contrast effect for clouds vs twilight            
              rintensity(:) = max(min(((glow_cld_c(:) -argref) * contrast + 128.),455.),0.)

              if(idebug .eq. 1 .AND. r_cloud_3d(i,j) .gt. 0.)then ! cloud present
                  write(6,51)cloud_rad_c(1:2,i,j),r_cloud_rad(i,j),pf_scat_moon,glow_secondary_cld &
                            ,glow_cld_moon,glow_cld_nt,glow_cld_c(2),nint(rintensity(:))
51                format('   cld_rad_c/crad/pf/glow: sec|moon|nt|cld rint',2e11.4,6f8.2,3i5)
              endif

          endif

          cld_red = nint(rintensity(1))                
          cld_grn = nint(rintensity(2))                        
          cld_blu = nint(rintensity(3))                         

!         Obtain brightness (glow) of clear sky
          if(sol_alt .gt. 0.)then ! Daylight from skyglow routine
              if(new_skyglow .eq. 1)then
                glow_tot = log10(clear_rad_c(2,i,j)) ! + log10(clear_radf_c(2,i,j))
                if(sun_vis .eq. 1.0)then
                    glow_tot = addlogs(glow_tot,glow_sun(i,j))
                endif
                glow_tot = addlogs(glow_tot,glow_moon(i,j))
!               glow_tot = log10(clear_rad_c(2,i,j))
                rintensity_glow = min(((glow_tot -7.3) * 100.),255.)

!               Convert RGB brightness into image counts preserving color balance
                gamclr = 0.70
                purple = 0.45
                bog_rad = clear_rad_c(3,i,j) / clear_rad_c(2,i,j) 
                if(bog_rad .gt. 2.)then ! Rayleigh Value
                    rog_exp = purple    ! boost red component
                elseif(bog_rad .lt. 1.)then ! Mie Yellow/Red
                    rog_exp = gamclr
                else
                    rog_exp = gamclr - (gamclr-purple) * (bog_rad - 1.0)
                endif
                rog = (clear_rad_c(1,i,j) / clear_rad_c(2,i,j))**rog_exp
                bog = (clear_rad_c(3,i,j) / clear_rad_c(2,i,j))**gamclr
                clr_red = rintensity_glow * rog
                clr_grn = rintensity_glow
                clr_blu = rintensity_glow * bog
              endif

          elseif(sol_alt .ge. -16.)then ! Twilight from clear_rad_c array
              call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),clear_rad_c_nt)
              glow_nt = log10(clear_rad_c_nt(3)) ! log nL           

              hue = clear_rad_c(1,i,j)
              sat = clear_rad_c(2,i,j)
              if(.false.)then
                glow_twi = glow(i,j) + log10(clear_rad_c(3,i,j)) ! phys val
              else ! new phys values (log nL)
                glow_twi =             log10(clear_rad_c(3,i,j)) ! phys val
              endif
              glow_twi = addlogs(glow_twi,glow_secondary_clr)
              glow_twi = addlogs(glow_twi,glow_nt)

              if(glow_stars(2,i,j) .gt. 1.0)then
!                 write(6,*)'i,j,glow_stars',i,j,glow_stars(2,i,j)
              endif
!             glow_stars(i,j) = 1.0                           ! test
              glow_tot = addlogs(glow_twi,glow_stars(2,i,j))  ! with stars 
              star_ratio = 10. ** ((glow_tot - glow_twi) * 0.45)

!             arg = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15             
              arg = glow_tot                                  ! experiment?            
              rintensity_floor = 0. ! 75. + sol_alt

!             rintensity_glow = max(min(((arg     -7.) * 100.),255.),rintensity_floor)
!             rintensity_glow = max(min(((arg -argref) * 100.),255.),rintensity_floor)
              rintensity_glow = max(min(((arg -argref) * contrast + 128.),255.),rintensity_floor)
!             rintensity_glow = min(rintensity_glow*star_ratio,255.)              
              call hsl_to_rgb(hue,sat,rintensity_glow,clr_red,clr_grn,clr_blu)
!             if(idebug .eq. 1)then
!                 write(6,61)glow_tot, argref, rintensity_glow,nint(clr_red),nint(clr_grn),nint(clr_blu)
61                format('Twi Clr RGB = ',3f9.3,3i5)
!             endif

          else ! Night from clear_rad_c array (default flat value of glow)
              call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),clear_rad_c_nt)
              hue = clear_rad_c_nt(1)
              sat = clear_rad_c_nt(2)
              glow_nt = log10(clear_rad_c_nt(3)) ! log nL           

              if(moon_alt .gt. 0.)then ! add moon mag condition
!                 Glow from Rayleigh, no clear_rad crepuscular rays yet
!                 argm = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15
                  glow_moon_s = glow(i,j)          ! log nL                 
                  glow_tot = addlogs(glow_nt,glow_moon_s)
              else
                  glow_tot = glow_nt
                  glow_moon_s = 0.
              endif

!             Add in stars. Stars have a background glow of 1.0
              glow_tot = addlogs(glow_tot,glow_stars(2,i,j))

!!            if((idebug .eq. 1 .and. moon_alt .gt. 0.) .OR. glow_stars(2,i,j) .gt. 1.0)then
!             if((idebug .eq. 1) .OR. glow_stars(2,i,j) .gt. 3.0)then
!                 write(6,91)i,j,idebug,glow_nt,glow_moon_s,glow_stars(2,i,j),glow_tot
91                format('   glow: nt/moon/stars/tot = ',3i5,4f9.3)
!                 idebug = 1 ! for subsequent writing at this grid point
!             endif

!             rintensity_glow = max(min(((glow_tot - 2.3) * 100.),255.),20.)
              rintensity_glow = max(min(((glow_tot - argref) * contrast + 128.),255.),20.)
              call hsl_to_rgb(hue,sat,rintensity_glow,clr_red,clr_grn,clr_blu)

          endif

!         Apply redness to clear sky / sun
!         Define area around the sun that is reddenned since scattering by aerosols
!         is dominant. Try and preserve both original luminance and desired color         
!         if(sol_alt .gt. 0)then
          elong_red = 0.25 ! only redden solar disk
          if(sol_alt .ge. -2.0)then
            if(elong_a(i,j) .le. elong_red)then
!           if(.false.)then                        
              clr_luma1 = .30 * clr_red + .59 * clr_grn + .11 * clr_blu
!             red_elong = ((elong_red - elong_a(i,j)) / elong_red)**2
              red_elong = 1.0
              red_sun = 255. ! Unless slant trans < 1./160000.
              clr_red =  red_sun 
              clr_grn = (red_sun / (rog_sun**(0.25))) 
              clr_blu = (red_sun / (rob_sun**(0.25))) 
              clr_luma2 = .30 * clr_red + .59 * clr_grn + .11 * clr_blu
              ratio_luma = min(clr_luma1 / clr_luma2,255. / clr_red)
              clr_red = clr_red * ratio_luma               
              clr_grn = clr_grn * ratio_luma              
              clr_blu = clr_blu * ratio_luma               
              write(6,*)' alt/azi/elong_red/elong/redelong: ',alt_a(i,j),azi_a(i,j),elong_red,elong_a(i,j),red_elong,nint(clr_red),nint(clr_grn),nint(clr_blu)
            endif
          endif

!                     Rayleigh  Ozone   Mag per optical depth            
          od_atm_g = (0.1451  + .016) / 1.086
          od_2_cloud = (od_atm_g + od_atm_a) * airmass_2_cloud(i,j)

!         Empirical correction to account for bright clouds being visible
          cloud_visibility = exp(-0.71*od_2_cloud) ! empirical coefficient

!         Use clear sky values if cloud cover is less than 0.5
          frac_cloud = r_cloud_3d(i,j)
!         frac_cloud = 0.0 ; Testing
          frac_cloud = frac_cloud * cloud_visibility  
          frac_clr = 1.0 - frac_cloud                         
          sky_rgb(0,I,J) = clr_red * frac_clr + cld_red * frac_cloud
          sky_rgb(1,I,J) = clr_grn * frac_clr + cld_grn * frac_cloud
          sky_rgb(2,I,J) = clr_blu * frac_clr + cld_blu * frac_cloud
          sky_rgb(:,I,J) = min(sky_rgb(:,I,J),255.)

!         if(idebug_a(i,j) .eq. 1 .AND. abs(alt_a(i,j)) .le. 2.0)then
!             write(6,*)'cld_red/sky_rgb',cld_red,sky_rgb(0,I,J)
!         endif

!         Use topo value if airmass to topo > 0
          if(airmass_2_topo(i,j) .gt. 0.)then

!             The sky RGB values should already include reductions from
!             cloud/terrain shadowing and limited distance to the topography
              if(sol_alt .le. 0.)then 
                od_2_topo = 0. !(od_atm_g + od_atm_a) * airmass_2_topo(i,j)
              else ! eventually use clear_rad influenced by topo?
                od_2_topo = (od_atm_g * airmass_2_topo(i,j)) + aod_2_topo(i,j)
!               od_2_topo = (od_atm_g * airmass_2_topo(i,j)) + aod_ill(i,j)
              endif

              topo_visibility = exp(-1.00*od_2_topo)                    

              if(airmass_2_cloud(i,j) .gt. 0. .AND. airmass_2_cloud(i,j) .lt. airmass_2_topo(i,j)) then
                  topo_visibility = topo_visibility * (1.0 - r_cloud_3d(i,j))
              endif

              if(sol_alt .gt. -4.0)then 
!                 Daytime assume topo is lit by sunlight (W/m**2)
                  topo_swi_frac = (max(topo_swi(i,j),001.) / 1000.) ** 0.45
                  rtopo_red = 120. * topo_swi_frac * (topo_albedo(1,i,j)/.15)**0.45
                  rtopo_grn = 120. * topo_swi_frac * (topo_albedo(2,i,j)/.15)**0.45
                  rtopo_blu = 120. * topo_swi_frac * (topo_albedo(3,i,j)/.15)**0.45
              else
!                 Nighttime assume topo is lit by city lights (nL)
                  topo_swi_frac = (max(topo_swi(i,j),001.) / 5000.) ** 0.45
                  rtopo_red = 120. * topo_swi_frac * 1.1
                  rtopo_grn = 120. * topo_swi_frac * 1.0
                  rtopo_blu = 120. * topo_swi_frac * 0.9
              endif

              sky_rgb(0,I,J) = nint(rad_to_counts(counts_to_rad(rtopo_red)*topo_visibility &
                             + counts_to_rad(sky_rgb(0,I,J))*(1.0-topo_visibility) ) )
              sky_rgb(1,I,J) = nint(rad_to_counts(counts_to_rad(rtopo_grn)*topo_visibility &
                             + counts_to_rad(sky_rgb(1,I,J))*(1.0-topo_visibility) ) )
              sky_rgb(2,I,J) = nint(rad_to_counts(counts_to_rad(rtopo_blu)*topo_visibility &
                             + counts_to_rad(sky_rgb(2,I,J))*(1.0-topo_visibility) ) )
          else
              od_2_topo = 0.
          endif

          sky_rgb(:,i,j) = min(sky_rgb(:,i,j) * ramp_night,255.)

          if(idebug .eq. 1)then
              rmaglim = b_to_maglim(10.**glow_tot)
              call apply_rel_extinction(rmaglim,alt_a(i,j),od_atm_g+od_atm_a)
              if(sol_alt .ge. 0.)then        ! daylight
                  write(6,102)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat(i,j),r_cloud_3d(i,j),cloud_od(i,j),cloud_od_sp(i,j,:),bkscat_alb(i,j) &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                      ,topo_swi(i,j),topo_albedo(1,i,j),aod_ill(i,j),od_2_topo,topo_visibility,cloud_visibility,rintensity_glow &
                      ,nint(sky_rgb(:,i,j)),nint(cld_red),nint(cld_grn),nint(cld_blu)
              elseif(sol_alt .ge. -16.)then ! twilight
                  write(6,103)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat(i,j),r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb(i,j) &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j) &
                      ,rintensity(1),glow_cld_nt,glow_cld,glow_twi &
                      ,glow_secondary_clr,rmaglim,cloud_visibility &
                      ,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c(:,i,j),nint(clr_red),nint(clr_grn),nint(clr_blu)                              
              else ! night
                  write(6,104)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat(i,j),r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb(i,j) &
!                     ,frac_cloud,airmass_2_cloud(i,j),cloud_rad_c(2,i,j)*1e3,rintensity(1) &
                      ,frac_cloud,airmass_2_cloud(i,j),cloud_rad_c(1,i,j)*1e3,rintensity(1) &
                      ,glow_cld_nt,glow_cld_moon,glow_cld,glow_secondary_clr,rmaglim &
                      ,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c_nt(:)
              endif
102           format(2i5,3f6.1,f9.3,f9.4,5f6.2,3f8.3,f8.4,f7.1,f9.3,f9.1,5f8.3,f9.2,2x,3i4,' cldrgb',1x,3i4)
103           format(2i5,3f9.2,f9.3,f9.4,4f9.3,f9.4,f7.1,f9.3,f9.1,4f9.3,f9.2,2x,3i4,' clrrad',2f9.5,f9.0,3i4)
104           format(2i5,3f9.2,f9.3,f9.4,4f9.3,f9.6,f7.1,f9.3,f9.3,4f9.3,f9.2,2x,3i4,' clrrad',3f8.2)
          endif

        enddo ! i
        enddo ! j

        return
        end
