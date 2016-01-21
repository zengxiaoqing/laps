
        subroutine get_sky_rgb(r_cloud_3d,cloud_od,cloud_od_sp,nsp, &
                   r_cloud_rad,cloud_rad_c,cloud_rad_w, &
                   clear_rad_c,l_solar_eclipse,i4time,rlat,rlon,eobs, & ! I
                   clear_radf_c,patm,patm_sfc,htmsl, &                  ! I
                   glow_sun,glow_moon,glow_stars, &                     ! I
                   od_atm_a,aod_ref,transm_obs,obs_glow_zen,isun,jsun, &! I
                   airmass_2_cloud,airmass_2_topo, &                    ! I
                   topo_gti,topo_albedo,gtic,dtic,btic,emic,albedo_sfc,&! I
                   aod_2_cloud,aod_2_topo,aod_ill,aod_ill_dir,aod_tot, &! I
                   dist_2_topo,topo_solalt,trace_solalt,          &     ! I
                   alt_a,azi_a,elong_a,ni,nj,azi_scale,sol_alt,sol_az, &! I
                   minalt,maxalt,minazi,maxazi, &                       ! I
                   twi_0,horz_dep,solalt_limb_true, &                   ! I
                   moon_alt,moon_az,moon_mag,corr1_in, &                ! I
                   sky_rgb)                                             ! O

        use mem_namelist, ONLY: r_missing_data,earth_radius,aero_scaleht,redp_lvl,fcterm
        include 'trigd.inc'

!       Statement functions
        addlogs(x,y) = log10(10.**x + 10.**y)
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)
        brt(a,ext) = 1.0 - exp(-a*ext)
        rad_to_counts(rad) = (log10(rad)-7.3)*contrast ! 7.1 matches 240
        counts_to_rad(counts) = 10.**((counts/contrast)+7.3)
        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        include 'rad.inc'
        include 'wac.inc'

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real cloud_od(ni,nj)        ! cloud optical depth
        real cloud_od_sp(ni,nj,nsp) ! cloud species optical depth
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,ni,nj)  ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real cloud_rad_w(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat) * rad
        real clear_rad_c(nc,ni,nj)  ! clear sky illumination
                                    ! local/input when sun is above/below horizon
        real clear_rad_2nd_c(nc,ni,nj) ! secondary scattering clear sky illumination
        real moon_rad_c(nc,ni,nj)   ! clear sky illumination from moon
        real clear_radf_c(nc,ni,nj) ! integrated fraction of gas illumin- 
                                    ! ated by the sun along line of sight  
                                    ! (accounting for Earth shadow+clouds)
                                    ! though possibly not topo?
                                    ! attenuated behind clouds
!       real clear_radf_c_eff(nc,ni,nj) ! accounts for airmass_2_topo         
        real clear_rad_c_nt(3,ni,nj)! night sky brightness
        real ag_2d(ni,nj)           ! gas airmass (topo/notopo)
        real glow_moon1(ni,nj)      ! glow (experimental)
        real glow_sun(ni,nj)        ! sunglow (log nl, extendd obj, extnct)
        real glow_moon(ni,nj)       ! moonglow (log b in nl, extended obj)
        real glow_moon_sc(ni,nj)    ! moonglow (log b in nl, scattered)
        real rad_moon_sc(nc,ni,nj)  ! moonglow (nl, scattered multispctrl)
        real glow_stars(nc,ni,nj)   ! starglow (log b in nl)           
        real airmass_2_cloud(ni,nj) ! airmass to cloud 
        real airmass_2_topo(ni,nj)  ! airmass to topo (rel to zenith @ std atmos)
        real topo_gti(ni,nj)        ! terrain normal global irradiance
        real topo_albedo(nc,ni,nj)  ! terrain albedo
        real gtic(nc,ni,nj)         ! spectral terrain normal irradiance
        real dtic(nc,ni,nj)         !    "        "    diffuse     "
        real btic(nc,ni,nj)         !    "        "    beam/direct "
        real emic(nc,ni,nj)         ! spectral exitance (from sfc lights)
        real aod_2_cloud(ni,nj)     ! future use
        real aod_2_topo(ni,nj)      ! aerosol optical depth to topo (slant)
        real aod_ill(ni,nj)         ! aerosol illuminated slant optical depth (topo/notopo)
        real aod_ill_dir(ni,nj)     ! aerosol directly slant illuminated optical depth, atten bhd clds 
        real aod_tot(ni,nj)         ! aerosol slant (in domain) optical depth 
        real dist_2_topo(ni,nj)     ! distance to topo (m)
        real topo_solalt(ni,nj)     ! solar altitude from ground
        real trace_solalt(ni,nj)    ! solar altitude from ray trace
        real od_atm_a               ! zenithal AOD seen from observer (aod_vrt)
        real od_atm_a_eff(ni,nj)    ! aerosol illuminated tau per airmass
        real od_atm_a_dir(ni,nj)    ! aerosol directly ill tau per airmass
        real aod_slant(ni,nj)       ! aerosol optical depth (slant path)
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real elong_a(ni,nj)         ! elong of sun or moon
        real elong_m(ni,nj)         ! elong of sun or moon
        real elong_s(ni,nj)         ! elong of sun or moon
        integer idebug_a(ni,nj), idebug_pf(ni,nj)
        real cld_radt(nc), cld_radb(nc), cld_rad(nc), sky_rad(nc)
        real rintensity(nc), cld_rgb_rat(nc), glow_cld_c(nc), ref_nl(nc)
        real rad_sec_cld(nc), albedo_sfc(nc)
        real pf_scat(nc,ni,nj), pf_scat1(nc,ni,nj), pf_scat2(nc,ni,nj)
        real pf_land(nc,ni,nj)
        real bkscat_alb(ni,nj)
        real rint_top(nc),rint_base(nc)
        real trans_c(nc)            ! transmissivity
        real od_g_slant_a(nc,ni)    ! use for sun/moon/star attenuation        
        real od_o_slant_a(nc,ni)    ! use for sun/moon/star attenuation
        real od_a_slant_a(nc,ni)    ! use for sun/moon/star attenuation
        real clr_od(nc), sky_rad_ave(nc), transterm(nc), sph_rad_ave(nc)

        real sky_rgb(0:2,ni,nj)
        real moon_alt,moon_az,moon_mag

        integer new_color /2/ ! sky_rad can be more fully used still

        logical l_solar_eclipse, l_sun_behind_terrain

        write(6,1)sol_alt,horz_dep,solalt_limb_true
1       format(' get_sky_rgb: sol_alt/horz_dep/solalt_limb_true = ',3f9.2)
        write(6,*)' moon alt/az/mag = ',moon_alt,moon_az,moon_mag
        write(6,*)' range of r_cloud_rad is ',minval(r_cloud_rad),maxval(r_cloud_rad)

        idebug_a = 0
        clear_rad_2nd_c(:,:,:) = 0. ! set (initialize) for testing
        day_int = day_int0 ! via includes

        angstrom_exp_a = 2.4 - (fcterm * 15.)
        do ic = 1,nc
            ext_a(ic) = (wa(ic)/.55)**(-angstrom_exp_a)
            write(6,*)' ic/wa/ext_a ',ic,wa(ic),ext_a(ic)
        enddo ! ic

!       htmsl = psatoz(patm*1013.25)

!       Sun/Cloud glow Calculation
!       http://www.astro.umd.edu/~ssm/ASTR620/mags.html#solarabsmag
!       Sun = 96000 lux = 1367 w/m**2 = -26.81 bolometric mag (TOA)
!       Vega = 2.54 * 1e-6 lux
!       1 lambert = 1 lumen / cm**2 = 1e4/pi candela/m**2
!       Sun is -26.74 mag per sphere (UBVRI: -25.96, -26.09, -26.74, -27.26, -27.55)
!       (180/pi)^2 * 4 * pi = 41252.961 square degrees in a sphere.
!       5.3463837e11 square arcsec in a sphere
!       0.76 mag per square arcsec
        rad_cld_day = v_to_b(-26.74 + log10(5.346e11)*2.5)
        write(6,*)' rad_cld_day (nl) = ',rad_cld_day
        rad_cld_day2 =  s10_to_nl(magsecsq_to_s10(0.76))
        write(6,*)' rad_cld_day2 (nl) = ',rad_cld_day2

        twi_alt = -4.5

        if(sol_alt .ge. twi_0)then
          corr1 = corr1_in; corr2 = 3.55  ! darkness of start/end of twilight (gamma)
          corr1 = corr1_in; corr2 = 3.75  ! darkness of start/end of twilight (gamma)
        else
          corr1 = corr1_in; corr2 = 3.55  ! darkness of start/end of twilight (gamma)
          corr1 = corr1_in; corr2 = 3.75  ! darkness of start/end of twilight (gamma)
        endif

        if(sol_alt .le. 0.)then

!         Use max sky brightness to calculate a secondary glow (twilight)
          call get_twi_glow_ave(            log10(clear_rad_c(3,:,:)) &
                     ,alt_a,azi_a,ni,nj,sol_alt,sol_az,twi_glow_ave)
          arg = maxval(log10(clear_rad_c(3,:,:))) ! avoid moon

          write(6,*)' twi_glow_ave = ',twi_glow_ave
          write(6,*)' twi_glow_max = ',arg

!         glow_secondary = arg - 1.5  ! one thirtieth of max sky brightness
          if(sol_alt .gt. twi_alt)then
              glow_diff_cld      = 1.1
          else
              glow_diff_cld      = 1.4 ! 0.8
          endif
          glow_secondary_cld = twi_glow_ave - glow_diff_cld

!         This is used only when sol_alt < twi_0
          glow_diff_clr      = 1.5 ! 1.5 ! 1.5
          glow_secondary_clr = twi_glow_ave - glow_diff_clr 

!         Last coefficient represents altitude of incident twilight glow
!         Try get_sky_glow_ave for twilight and/or daylight?
!         rindirect = 1300. * (10.**glow_secondary_cld / rad_cld_day) * 0.1
          rad_sec_cld(:) = 10.**glow_secondary_cld

          write(6,*)' glow_secondary_cld = ',glow_secondary_cld
          if(sol_alt .lt. twi_0)then
              write(6,*)' glow_secondary_clr = ',glow_secondary_clr
          endif
!         write(6,*)' rindirect = ',rindirect                 
          write(6,2)rad_sec_cld(:)
2         format('  rad_sec_cld (based on glow_secondary_cld) = ',3f12.0)

!         Note that 'argref' is applied to the clear sky from solar alt
!         0 down to -16 and applied to clouds just from -4 down to -16.
          sol_alt_eff = max(sol_alt,-16.)
!         Test at -3.2 solar alt
          fracalt = sol_alt_eff / (-16.)

!         A higher coefficient darkens mid-twilight
!         argref = corr1 + (fracalt**0.83 * (corr2-corr1))
!         argref = corr1 + (fracalt**0.68 * (corr2-corr1))
!         argref = corr1 + (fracalt**0.40 * (corr2-corr1))
!         fraccorr = 0.1 * fracalt + 0.9 * scurve(fracalt)
!         fraccorr = scurve(fracalt**0.40)
          fraccorr = scurve(fracalt**(0.45-fracalt*0.1))
          argref = corr1 + fraccorr * (corr2-corr1)

          argbri = corr1-argref
!         contrast = 70. + 30. * (abs(sol_alt_eff + 8.)/8.)**2 ! each image
          contrast = 54. + 30. * (abs(sol_alt_eff + 8.)/8.)**2 ! each image
          write(6,*)' argref/argbri/fraccorr = ',argref,argbri,fraccorr
          glwlow = 7.3 - argbri

        else ! sun above horizon
          glow_secondary_cld = 1. ! dark default in daytime
!         2.0 is mean airmasses relative to zenith, 0.5 is 90deg phase func
          sb_corr = 2.0 * (1.0 - (sind(sol_alt)**0.5)) 
          rad_sec_cld(:) = (day_int * ext_g(:) * patm_sfc * 2.0 * 0.5) / 10.**(0.4*sb_corr)
          write(6,23)rad_sec_cld(:)
23        format('  rad_sec_cld (based on day_int) = ',3f12.0)
          argbri = sb_corr * 0.30
          glwlow = 7.3 - argbri
          write(6,*)' sb_corr/argbri = ',sb_corr,argbri
          contrast = 100.

        endif

        write(6,*)' rad_sec_cld (top of code) = ',rad_sec_cld(:)

!       First arg (altmidcorr) increase will brighten all altitudes
!       Second arg increase will darken shallow twilight and brighten
!       deep twilight.
!       Higher 'fracerf0' will darken sunrise/set.
!       Brighten the mid-twilight image with thicker tropospheric aerosols.
!       Darken the mid-twilight image with thicker stratospheric aerosols?
        if(.false.)then ! avoid saturation on the bright end
          altmidcorr = -5.95 + (od_atm_a * 2.5) 
          fracerf0 = 0.69 ! value with sun on horizon
          write(6,*)' altmidcorr: ',altmidcorr
        else              ! show more like the camera
          alt_top = alt_a(ni,1)
          if(alt_top .eq. 90.)then ! fisheye lens
            altmidcorr = -3.60
            fracerf0 = 0.59  ! value with sun on horizon
          else                     ! panoramic camera
            altmidcorr = -4.89
            fracerf0 = 0.65  ! value with sun on horizon
          endif
          write(6,*)' alt_top,altmidcorr: ',alt_top,altmidcorr
        endif
        if(sol_alt .gt. altmidcorr)then ! shallow twilight
          fracerf = (sol_alt - altmidcorr) * (-fracerf0/altmidcorr)
        else                            ! deep twilight
          fracerf = (sol_alt - altmidcorr) * 0.220 
        endif
        erfterm = (erf(fracerf) + 1.) / 2.
        glwmid = corr2*(1.-erfterm) + corr1*erfterm

        offset = 0.
        glwlow = glwmid - 128. / 100. ! contrast
        write(6,3)sol_alt,corr1,glwlow,glwmid,contrast,fracerf,erfterm
3       format('  sol_alt/corr1/glwlow/glwmid/contrast/fracerf/erfterm',f9.2,3f9.3,f9.1,2f9.3)

        ref_nl = day_int0
        if(new_color .gt. 0)then
          call nl_to_RGB(ref_nl(:),glwmid,contrast & 
                        ,128.,1,ref_red,ref_grn,ref_blu)
        else
          ref_red = 240.; ref_grn = 240.; ref_blu = 240.
        endif
        refcolmax = max(ref_red,ref_grn,ref_blu)
        refred1 = (ref_red/refcolmax) * 255.
        refgrn1 = (ref_grn/refcolmax) * 255.
        refblu1 = (ref_blu/refcolmax) * 255.
        write(6,4)new_color,ref_red,ref_grn,ref_blu,refred1,refgrn1,refblu1
4       format(' newc/reference RGB = ',i2,3f7.1,3x,3f7.1)

!       Redness section (sun at low solar altitude)
!       Compare with 'get_cloud_rad' redness calculation?
        sol_alt_red_thr = 9.0 + (od_atm_a * 20.)
        if(sol_alt .le. sol_alt_red_thr .and. sol_alt .gt. -3.0)then
!           Add aerosols to trans calculation (call 'get_airmass')
            am = airmassf(90.-sol_alt,patm)
            do ic = 1,nc
              trans_c(ic) = trans(am*ext_g(ic)*patm)
            enddo
!           rob_sun = trans_c(1) / trans_c(3)
!           rog_sun = trans_c(1) / trans_c(2)
        else
!           rob_sun = 1.
!           rog_sun = 1.
            am = -99.9
        endif
!       write(6,5)sol_alt_red_thr,od_atm_a,am,rob_sun,rog_sun
!5      format('  sol_alt_red_thr/od_atm_a/am/rob_sun/rog_sun = ',5f9.3)

!       Determine moon_cond_clr and azid1/d2 
        azid1 = 90. ; azid2 = 270. ! default
!       if(sol_alt .gt. 0. .and. isun .ge. 1  .and. jsun .ge. 1 &
!                          .and. isun .le. ni .and. jsun .le. nj)then
        if(sol_alt .gt. 0. .and. jsun .ge. 1 &
                           .and. jsun .le. nj)then ! solar azimuth
            azid1 = azi_a(1,jsun) ! nint(sol_az)
            azid2 = mod(azid1+180.,360.)
!           azid2 = azid1 ! block antisolar az               
            moon_cond_clr = 0
        elseif(moon_alt .gt. -6. .and. sol_alt .lt. twi_0)then
            azid1 = nint(moon_az)
            azid2 = mod(azid1+180.,360.)
            if(sol_alt .gt. -8.)then ! block moon alt/az
                azid1 = 90. ; azid2 = 270.
!               azid1 = nint(sol_az)-40; azid2 = nint(sol_az)+40
            endif
            moon_cond_clr = 1
        else ! generic night
            azid1 = 90. ; azid2 = 270.
            moon_cond_clr = 0
        endif
!       if(htmsl .gt. 1000e3)then
!           azid1 = 215. ; azid2 = 215. ! custom
!       endif
!       azid1 = 179. ; azid2 = 181. ! custom volcanic

        write(6,*)' azid1/2 are at ',azid1,azid2

        if(moon_cond_clr .eq. 0)then ! sun
            write(6,*)' fill elong_a array based on solar elongation'
            call get_elong_a(1,ni,1 &                                ! I
                     ,1,nj,1 &                                       ! I
                     ,1,ni,1,nj &                                    ! I
                     ,sol_alt,sol_az,alt_a,azi_a &                   ! I
                     ,elong_a                                 )      ! O
            elong_s = elong_a
        else                    ! moon
            write(6,*)' fill elong_a array based on lunar elongation'
            call get_elong_a(1,ni,1 &                                ! I
                     ,1,nj,1 &                                       ! I
                     ,1,ni,1,nj &                                    ! I
                     ,moon_alt,moon_az,alt_a,azi_a &                 ! I
                     ,elong_a                                 )      ! O

            elong_m = elong_a

            call get_elong_a(1,ni,1 &                                ! I
                     ,1,nj,1 &                                       ! I
                     ,1,ni,1,nj &                                    ! I
                     ,sol_alt,sol_az,alt_a,azi_a &                   ! I
                     ,elong_s                                 )      ! O
        endif

        do j = 1,nj
        do i = 1,ni
            if(abs(azi_a(i,j)-azid1) .le. .001 .OR. & 
               abs(azi_a(i,j)-azid2) .le. .001     )then ! constant azimuth
                if(abs(alt_a(i,j)+horz_dep) .le. 3.4)then ! near horizon/limb
                    idebug_a(i,j) = 1
                elseif(abs(alt_a(i,j)) .le. 20. .AND. &
                       alt_a(i,j) .eq. float(nint(alt_a(i,j))))then
                    idebug_a(i,j) = 1
!               elseif(abs(alt_a(i,j)) .le. 21.)then   
!                   idebug_a(i,j) = 1
                elseif(alt_a(i,j) .ge. -75. .AND. &
                       alt_a(i,j) .le. -65. .AND. &
                       alt_a(i,j) .eq. float(nint(alt_a(i,j))))then
                    idebug_a(i,j) = 1
                elseif(alt_a(i,j) .ge. -78. .AND. &
                       alt_a(i,j) .le. -75.           )then
                    idebug_a(i,j) = 1
                elseif(alt_a(i,j) .eq. float((nint(alt_a(i,j))/5)*5))then
                    idebug_a(i,j) = 1
                endif
            endif
!           if(alt_a(i,j) .eq. +45.0 .and. azi_a(i,j) .eq. nint(azi_a(i,j)))then
!               idebug_a(i,j) = 1
!           endif
!           if(alt_a(i,j) .eq. +15.0)then
!               idebug_a(i,j) = 1
!           endif
        enddo ! i
        enddo ! j

        if(isun .gt. 0 .and. isun .le. ni .and. jsun .gt. 0 .and. jsun .le. nj)then
          idebug_a(isun,jsun) = 1
        endif

        if(sol_alt .gt. twi_0)then

            I4_elapsed = ishow_timer()

            write(6,*)' isun,jsun=',isun,jsun

!           Determine if observer is in terrain shadow
!           Note sun_vis is used only for sun glowing below horizon
            if(.not.(isun .ge. 1  .and. jsun .ge. 1 .and. &
                     isun .le. ni .and. jsun .le. nj))then
               l_sun_behind_terrain = .false. ! outside window
            elseif(airmass_2_topo(isun,jsun) .gt. 0.)then 
               l_sun_behind_terrain = .true. 
            else
               l_sun_behind_terrain = .false. 
            endif

            if(solalt_limb_true .lt. 0.)then ! below horizon
               write(6,*)' Sun is below horizon/limb',sol_alt,solalt_limb_true
               l_sun_behind_terrain = .true. 
            endif
                
            if(l_sun_behind_terrain .eqv. .true.)then 
                write(6,*)' Sun is behind terrain'
                od_atm_a_eff = 0.
                od_atm_a_dir = 0.
                sun_vis = 0.
            else ! sun is outside of terrain
                write(6,*)' Sun is outside of terrain'
                if(transm_obs .gt. 0.9)then
                    od_atm_a_eff = od_atm_a
                    od_atm_a_dir = od_atm_a
                    sun_vis = 1.
                else ! in cloud/terrain shadow
                    od_atm_a_eff = 0.
                    od_atm_a_dir = 0.
                    sun_vis = 1. ! (render attenuated sun better)
!                   sun_vis = 0. 
                endif
            endif

!           Calculate normalized illuminated aerosol optical depth (zenith)
            iprint1 = 0
            iprint2 = 0
            write(6,*)' alt  od_atm_a  aod2topo    aod_ill/dir  od_atm_test/dir  clrd'
            do i = ni,1,-1
              call get_airmass(alt_a(i,1),htmsl,patm &   ! I
                              ,redp_lvl,aero_scaleht &   ! I
                              ,earth_radius,1 &          ! I
                              ,ag,ao,aa,refr_deg)        ! O 
              ag_2d(i,:) = ag
              aod_slant(i,:) = aa * od_atm_a ! od_atm_a (should) = aod_vert 
              do j = 1,nj
!               Crepuscular ray experiments
!               if(airmass_2_topo(i,j) .gt. 0. .OR. transm_obs .lt. 0.9)then
                if(airmass_2_topo(i,j) .gt. 0.)then ! inside terrain 
!                 Note that aod_ill_dir discriminates between direct/diffuse
                  od_atm_a_eff(i,j) = od_atm_a * aod_ill(i,j)     / aod_2_topo(i,j)
                  od_atm_a_dir(i,j) = od_atm_a * aod_ill_dir(i,j) / aod_2_topo(i,j)
!                 clear_radf_c_eff(:,i,j) = clear_radf_c(:,i,j) * airmass_2_topo(i,j) / ag
                  ag_2d(i,j) = airmass_2_topo(i,j)
                  if(idebug_a(i,j) .eq. 1)then
                    od_atm_test = od_atm_a * aod_ill(i,j) / aod_2_topo(i,j)
                    iprint1 = iprint1 + 1
                    if(iprint1 .eq. 1)then
                        write(6,*)'    alt    azi     od_atm_a  aod2topo      aod_ill/dir     od_atm_test/dir    clrdf'
                    endif
                    write(6,6)alt_a(i,j),azi_a(i,j),od_atm_a,aod_2_topo(i,j) &
                      ,aod_ill(i,j),aod_ill_dir(i,j),od_atm_test,od_atm_a_dir(i,j) &   
                      ,clear_radf_c(2,i,j) ! ,clear_radf_c_eff(2,i,j)
6                   format(3f9.3,3f10.5,4f9.3)
                  endif

                else ! outside of terrain
                  if(alt_a(i,j) .ge. 0. .and. aa .gt. 0.)then
!                   od_atm_a_eff(i,j) = (od_atm_a * aod_ill(i,j))     / (aa * od_atm_a)
!                   od_atm_a_dir(i,j) = (od_atm_a * aod_ill_dir(i,j)) / (aa * od_atm_a)
                    od_atm_a_eff(i,j) = (           aod_ill(i,j))     / (aa           )
                    od_atm_a_dir(i,j) = (           aod_ill_dir(i,j)) / (aa           )
                  else ! might make use of more accurate aa below horizon
!                   od_atm_a_eff(i,j) = 1.0
!                   od_atm_a_dir(i,j) = 1.0
                    od_atm_a_eff(i,j) = (           aod_ill(i,j))     / (aa           )
                    od_atm_a_dir(i,j) = (           aod_ill_dir(i,j)) / (aa           )
                  endif
                  if((i .eq. isun .AND. j .eq. jsun) .OR. &
                     (alt_a(i,j) .lt. 0. .and. idebug_a(i,j) .eq. 1) .OR. &
                     (clear_radf_c(2,i,j) .lt. 0.5 .AND. iprint2 .le. 30) .OR. &
                     (cloud_od(i,j) .eq. 0. .AND. &
                      od_atm_a_eff(i,j)/od_atm_a .le. 0.5 .AND. &
                      alt_a(i,j) .gt. 10.0                 .AND. &
                      iprint2 .le. 30)                              )then
                      if(i .eq. isun .AND. j .eq. jsun)then   
                          write(6,*)'Calcs at solar location outside of terrain:'
                          write(6,*)'cloud_od (slant/vert) = ',cloud_od(i,j),cloud_od(i,j)*sind(sol_alt)
                      elseif(clear_radf_c(2,i,j) .lt. 0.5)then
                          write(6,*)'Clear with cloud shadowing outside of terrain:'
                      else
                          write(6,*)'Clear with aerosol shadowing outside of terrain:'
                      endif
                      write(6,7)alt_a(i,j),azi_a(i,j),od_atm_a,aod_tot(i,j),aod_ill(i,j),aod_ill_dir(i,j),aa &
                               ,od_atm_a_eff(i,j),od_atm_a_dir(i,j),clear_radf_c(2,i,j)
7                     format(' altaz/od_atm_a/aod_tot/ill/dir/aa/od_atm_a_eff/dir/clr_radf =   ',2f9.2,5f9.3,3f9.3)
                      iprint2 = iprint2 + 1
                  endif
!                 clear_radf_c_eff(:,i,j) = clear_radf_c(:,i,j)
                endif

!               if(idebug_a(i,j) .eq. 1)then
!                   if(aod_ill_dir(i,j) .gt. aod_ill(i,j))then
!                       write(6,71)alt_a(i,j),azi_a(i,j),aod_ill(i,j),aod_ill_dir(i,j) 
!71                     format('WARNING: large aod_ill_dir at',2f9.2,2f9.4)
!                   endif
!               endif

              enddo ! j
            enddo ! i

            if(isun .gt. 0 .and. isun .le. ni .and. jsun .gt. 0 .and. jsun .le. nj)then
              write(6,8)transm_obs,od_atm_a,od_atm_a_eff(isun,jsun),sun_vis
8             format(' transm_obs,od_atm_a,od_atm_a_eff,sun_vis=',4f8.3)
            endif

            write(6,*)' od_atm_a_eff range = ',minval(od_atm_a_eff),maxval(od_atm_a_eff)
            write(6,*)' od_atm_a_dir range = ',minval(od_atm_a_dir),maxval(od_atm_a_dir)

!           if(maxval(od_atm_a_dir) .gt. maxval(od_atm_a_eff))then
!               write(6,*)' WARNING: od_atm_a_dir larger than od_atm_a_eff'
!           endif

            I4_elapsed = ishow_timer()

            write(6,*)' range of azi_a = ',minval(azi_a),maxval(azi_a)
                 
            isolalt_lo=-21; isolalt_hi=+21
            write(6,*)' call skyglow_phys for daytime or solalt > twi_0:'

!           Introduced airmass_2_topo effect in ag_2d
            call skyglow_phys(minalt,maxalt,1 &                        ! I
                     ,minazi,maxazi,1,azi_scale &                      ! I
                     ,minalt,maxalt,minazi,maxazi,idebug_a &           ! I
                     ,sol_alt,sol_az,alt_a,azi_a,twi_0,twi_alt &       ! I
                     ,isolalt_lo,isolalt_hi,topo_solalt,trace_solalt & ! I
                     ,earth_radius,patm &                              ! I
                     ,od_atm_a,od_atm_a_eff,od_atm_a_dir &             ! I
                     ,aod_ref,aero_scaleht,dist_2_topo &               ! I
                     ,htmsl,redp_lvl,horz_dep &                        ! I
                     ,aod_ill,aod_2_topo,aod_tot &                     ! I
                     ,l_solar_eclipse,i4time,rlat,rlon &               ! I
                     ,clear_radf_c,ag_2d &                             ! I
                     ,od_g_slant_a,od_o_slant_a,od_a_slant_a,ext_a &   ! O
                     ,clear_rad_2nd_c &                                ! O
                     ,clear_rad_c,sky_rad_ave,elong_s         )      ! O/I

            write(6,*)' range of clear_radf_c(2) is ',minval(clear_radf_c(2,:,:)),maxval(clear_radf_c(2,:,:))
            write(6,*)' range of clear_rad_c(2) is ',minval(clear_rad_c(2,:,:)),maxval(clear_rad_c(2,:,:))
            if(minval(clear_rad_c(2,:,:)) .le. 0.)then
                write(6,*)' ERROR: clear_rad_c(2,:,:) has min <= 0.'
            endif
            if(jsun .gt. 0 .and. jsun .le. nj)then
              write(6,*)' clear_rad_c(2) in solar column:',clear_rad_c(2,:,jsun)
            endif
!           Returned when 'sol_alt' is between 'twi_0' and 10.
            write(6,*)' sky_rad_ave = ',sky_rad_ave(:)

!           Add sfc_albedo?
            where(sky_rad_ave(:) .ne. r_missing_data)                   
              sph_rad_ave(:) = sky_rad_ave(:) * 0.5 * (1. + 0.15)
            elsewhere
              sph_rad_ave(:) = r_missing_data
            endwhere
            write(6,*)' sph_rad_ave = ',sph_rad_ave(:)

        elseif(moon_cond_clr .eq. 1)then 
            od_atm_a_eff = od_atm_a
            od_atm_a_dir = od_atm_a

            write(6,*)' range of clear_radf_c is ',minval(clear_radf_c(2,:,:)),maxval(clear_radf_c(2,:,:))
            write(6,*)' Calling get_airmass for all altitudes'
            do i = ni,1,-1
              call get_airmass(alt_a(i,1),htmsl,patm &   ! I
                              ,redp_lvl,aero_scaleht &   ! I
                              ,earth_radius,0 &          ! I
                              ,ag,ao,aa,refr_deg)        ! O 
              ag_2d(i,:) = ag
            enddo ! i

            I4_elapsed = ishow_timer()

            isolalt_lo=-11; isolalt_hi=+11
            write(6,*)' call skyglow_phys for moon testing:'
            call skyglow_phys(minalt,maxalt,1 &                        ! I
                     ,minazi,maxazi,1,azi_scale &                      ! I
                     ,minalt,maxalt,minazi,maxazi,idebug_a &           ! I
                     ,moon_alt,moon_az,alt_a,azi_a,twi_0,twi_alt &     ! I
                     ,isolalt_lo,isolalt_hi,topo_solalt,trace_solalt & ! I
                     ,earth_radius,patm &                              ! I
                     ,od_atm_a,od_atm_a_eff,od_atm_a_dir &             ! I
                     ,aod_ref,aero_scaleht,dist_2_topo &               ! I
                     ,htmsl,redp_lvl,horz_dep &                        ! I
                     ,aod_ill,aod_2_topo,aod_tot &                     ! I
                     ,.false.,i4time,rlat,rlon &                       ! I
                     ,clear_radf_c,ag_2d &                             ! I
                     ,od_g_slant_a,od_o_slant_a,od_a_slant_a,ext_a &   ! O
                     ,clear_rad_2nd_c &                                ! O
                     ,moon_rad_c,sky_rad_ave,elong_a         )       ! O/I

            glow_moon_sc(:,:) = 0.; rad_moon_sc(:,:,:) = 0.
            where(moon_rad_c(2,:,:) .gt. 0.)
              glow_moon_sc(:,:) = log10(moon_rad_c(2,:,:)) + (-26.7 - moon_mag) * 0.4
              do ic = 1,nc
                rad_moon_sc(ic,:,:) = log10(moon_rad_c(ic,:,:)) + (-26.7 - moon_mag) * 0.4
              enddo ! ic
            endwhere
            glow_moon1(:,:) = glow_moon_sc(:,:) ! experimental

            write(6,*)' range of moon_rad_c(2) is ',minval(moon_rad_c(2,:,:)),maxval(moon_rad_c(2,:,:))
            write(6,*)' log correction is ',(-26.7 - moon_mag) * 0.4
            write(6,*)' range of glow_moon_sc is ',minval(glow_moon_sc(:,:)),maxval(glow_moon_sc(:,:))
            write(6,*)' range of glow_moon1 is ',minval(glow_moon1(:,:)),maxval(glow_moon1(:,:))
            write(6,*)' moonglow1:',glow_moon1(ni/2,1:nj:10)

        else ! get just solar elongation
            if(isun .gt. 0 .AND. jsun .gt. 0 .AND. isun .le. ni .AND. jsun .le. nj)then
                write(6,*)' isun,jsun,am_2_topo=',isun,jsun,airmass_2_topo(isun,jsun)

!               Determine if observer is in terrain shadow
                if(airmass_2_topo(isun,jsun) .gt. 0.)then 
                    write(6,*)' Sun is behind terrain/horizon'
                    sun_vis = 0.
                else
                    write(6,*)' Sun is above terrain/horizon'
                    sun_vis = 1.
                endif
            else
                write(6,*)' Sun is outside window'
                sun_vis = 0.
            endif

            write(6,*)' Calling get_airmass for all altitudes'
            do i = ni,1,-1
              call get_airmass(alt_a(i,1),htmsl,patm &   ! I
                              ,redp_lvl,aero_scaleht &   ! I
                              ,earth_radius,0 &          ! I
                              ,ag,ao,aa,refr_deg)        ! O 
              od_g_slant_a(:,i) = ext_g(:) * ag 
              od_o_slant_a(:,i) = (o3_du/300.) * ext_o(:) * ao
              od_a_slant_a(:,i) = aod_ref * aa * ext_a(:)
            enddo ! i
        endif

        if(sol_alt .ge. 0.)then
            htmin_sun = htmsl
        else
            htmin_sun = htminf(htmsl,sol_alt,earth_radius)
        endif

        write(6,*)' range of glow_sun is ',minval(glow_sun),maxval(glow_sun)
        write(6,*)' range of glow_moon is ',minval(glow_moon),maxval(glow_moon)
        write(6,*)' sun_vis/isun/jsun = ',sun_vis,isun,jsun
        write(6,*)' htmin_sun = ',htmin_sun

        I4_elapsed = ishow_timer()

        bkscat_alb = -99.9 ! dummy value to initialize for logging
        write(6,*)' max of idebug_a (1) = ',maxval(idebug_a)
        if(sol_alt .lt. twi_alt)then
            idebug_pf = 0
        else
            idebug_pf = idebug_a
        endif
        call get_cld_pf(elong_a,alt_a,r_cloud_rad,cloud_rad_w,cloud_od,cloud_od_sp & ! I
                       ,nsp,airmass_2_topo,idebug_pf,ni,nj &  ! I
                       ,pf_scat1,pf_scat2,pf_scat,bkscat_alb) ! O

        I4_elapsed = ishow_timer()

        call get_lnd_pf(elong_a,alt_a,azi_a,topo_gti,topo_albedo    & ! I
                       ,transm_obs                                  & ! I
                       ,gtic,dtic,btic                              & ! I
                       ,dist_2_topo,topo_solalt,azi_scale           & ! I
                       ,sol_alt,sol_az,nsp,airmass_2_topo,idebug_pf & ! I
                       ,ni,nj,pf_land)                                ! O

        I4_elapsed = ishow_timer()

        call get_clr_rad_nt_2d(alt_a,ni,nj,obs_glow_zen &             ! I
                              ,patm,htmsl,horz_dep &                  ! I
                              ,airmass_2_topo &                       ! I
                              ,clear_rad_c_nt)                        ! O

        I4_elapsed = ishow_timer()

        write(6,*)' albedo_sfc = ',albedo_sfc(:)
        write(6,*)
        if(ni .eq. nj)then ! polar
            write(6,*)' slice from SW to NE through midpoint'
        else
            write(6,*)' slices at degrees azimuth: ',azid1,azid2
        endif
        write(6,*)

        if(.true.)then
!         if(sol_alt .ge. -2.0)then     ! daylight
          if(sol_alt .ge. twi_alt)then  ! daylight
            write(6,11)
11          format('    i   j   alt   azi  elong   pf_scat     opac     od /         species            alb    cloud airmass   rad', &
                   '   rinten  airtopo   gtitopo topoalb aodill  topood topovis cld_visb   glow      skyrgb')
          elseif(sol_alt .ge. -16.)then ! twilight
            write(6,12)
12          format('    i   j      alt      azi     elong   pf_scat     opac       od      alb     cloud  airmass   rad    ', &
                   'rinten  glwcldnt glwcld  glwnt glwtwi gw2clr glwtot maglim cld_visb  glow      skyrgb')
          else                          ! night
            write(6,13)
13          format('    i   j      alt      azi     elong   pf_scat     opac       od      alb     cloud  airmass   rade-3 ', &
                   'rinten glw_cld_nt glwcldmn glw_cld glw2ndclr rmaglim  cld_visb  glow      skyrgb')
          endif
        endif

        do j = 1,nj
        do i = 1,ni

          sky_rgb(:,i,j) = 0.          
          clr_od(:) = od_g_slant_a(:,i) + od_o_slant_a(:,i) &
                    + od_a_slant_a(:,i)
          if(azi_a(i,j) .eq. azid1 .OR. azi_a(i,j) .eq. azid2)then ! constant azimuth
              if(i .eq. 1)write(6,*)   
          endif

          idebug = idebug_a(i,j)

!         Determine relevant solar altitude along ray
          if(htmsl .gt. 150e3)then ! highalt strategy
            if(alt_a(i,j) .gt. 0.)then
              solalt_ref = sol_alt
            else
              if(trace_solalt(i,j) .ne. sol_alt)then ! valid trace
                solalt_ref = trace_solalt(i,j)
              else
                solalt_ref = sol_alt - (altray*cosd(view_azi_deg-sol_azi))
                solalt_ref = min(solalt_ref,+180.-solalt_ref)
                solalt_ref = max(solalt_ref,-180.-solalt_ref)
              endif
            endif
          else
            solalt_ref = topo_solalt(i,j)
          endif

!         Obtain cloud brightness

!         if(sol_alt .ge. twi_alt)then ! Day/twilight from cloud_rad_c array
!         if(topo_solalt(i,j) .ge. twi_alt)then ! Day/twilight from cloud_rad_c array
          if(solalt_ref .ge. twi_alt)then ! Day/twilight from cloud_rad_c array
              if(solalt .le. 0.)then ! between shallow twilight and 0.
                  where(sph_rad_ave .ne. r_missing_data)rad_sec_cld = sph_rad_ave
              endif
!             Potential intensity of cloud if it is opaque 
!               (240. is nominal intensity of a white cloud far from the sun)
!                217 would better match -7.3 value above
!               (0.25 is dark cloud base value)                                  
!             Use surface albedo equation?
              do ic = 1,nc
                  rad = day_int             * pf_scat(ic,i,j)

!                 Improve sunset clouds
!                 rad = rad * cloud_rad_c(1,i,j) / r_cloud_rad(i,j) 
!                 rint_top(ic) = rad_to_counts(rad)

                  cld_radt(ic) = rad * cloud_rad_c(ic,i,j) + rad_sec_cld(ic)
                  rint_top(ic) = rad_to_counts(cld_radt(ic))

                  if(sol_alt .gt. 0.)then
!                     Calculate 'sky_rad_ave' upstream for each color
!                     topo_arg = sky_rad_ave(ic) * albedo_sfc(ic)
                      rad = day_int * (0.02 * (1. + albedo_sfc(ic)))   * pf_scat(ic,i,j)
                  else ! secondary scattering term allowed to dominate
                      rad = 0.
                  endif
                  cld_radb(ic) = rad + rad_sec_cld(ic)

!                 Improve sunset clouds
!                 rad = rad * cloud_rad_c(1,i,j) / r_cloud_rad(i,j)
!                 rint_base(ic) = rad_to_counts(rad)

                  rint_base(ic) = rad_to_counts(cld_radb(ic))
              enddo ! ic

!             Gamma color correction applied when sun is near horizon 
!             Add reduction based on sun's clear air attenutation as well?
              if(cloud_rad_c(1,i,j) .gt. 0.)then
                  cld_rgb_rat(:) = (cloud_rad_c(:,i,j)/cloud_rad_c(1,i,j))**0.45
              else
                  cld_rgb_rat(:) = 1.0
              endif

!             rintensity = min(rintensity,255.)
!             rintensity(:) = rint_top(:)  * cld_rgb_rat(:) * r_cloud_rad(i,j) &
!                           + rint_base(:) * (1.0 - r_cloud_rad(i,j))
!             rintensity = max(rintensity(:),0.)
              cld_rad(:) = cld_radt(:) * r_cloud_rad(i,j) & ! * cld_rgb_rat(:) 
                         + cld_radb(:) * (1.0 - r_cloud_rad(i,j))

!             Apply cloud reddening from viewer intervening airmass
!             This may be redundant with cloud_visibility
              do ic = 1,nc
!                 trans_c(ic) = trans(ext_g(ic) * airmass_2_cloud(i,j)) 
                  trans_c(ic) = 1.0 ! if redundant with cloud_visibility
              enddo
!             rintensity(:) = rintensity(:) * (trans_c(:)/trans_c(1))**0.25 ! 0.45
              cld_rad(:) = cld_rad(:) * trans_c(:)

!             Informational only
              call nl_to_RGB(cld_rad(:),glwmid,contrast & 
                        ,128.,0,rintensity(1),rintensity(2),rintensity(3))

              if( ( idebug_a(i,j) .eq. 1 .AND. alt_a(i,j) .eq. nint(alt_a(i,j)) .AND. & 
                (abs(alt_a(i,j)).le.40.0 .or. abs(alt_a(i,j)).eq.90.0 .or. abs(alt_a(i,j)+70.).lt.7.) )  &
                  .OR. (i .eq. isun .and. j .eq. jsun) )then
               if(r_cloud_3d(i,j) .gt. 0.)then
                  write(6,*)' solalt_ref/twi_alt ',solalt_ref,twi_alt
                  write(6,42)i,j,elong_a(i,j),pf_scat(2,i,j),rintensity(2)&
                   ,trans_c(2),r_cloud_rad(i,j) &
                   ,cloud_rad_c(:,i,j),cld_radt(:)/1e6 &
                   ,cld_radb(:)/1e6,(cld_rad(:)/1e6)/trans_c(:) &
                   ,cld_rad(:)/1e6,rad_sec_cld(:)/1e6
 42               format(&
                  ' elg/pf/rint2/trnsc2/rcldrd/cldrdtba/cldrad/rdsc = ',2i5,5f9.3,' c',3f8.5,2x,3f6.0,2x,3f6.0,2x,3f6.0,2x,3f6.0,2x,3f5.0)
               endif ! cloud present
              endif

          else ! later twilight (clear_rad_c) and nighttime (surface lighting)

!             This might consider the optical thickness and albedo of the cloud as well?
              if(sph_rad_ave(2) .eq. r_missing_data)then
                  glow_cld1 = glow_secondary_cld ! + log10(r_cloud_rad(i,j))                                     
              elseif(sph_rad_ave(2) .gt. 0.)then
                  if(moon_cond_clr .eq. 0)then
                      glow_cld1 = log10(sph_rad_ave(2))
                  else ! moon present - reduce sph_rad_ave
                      ratio_moonsun = 10. ** ((-26.7 - moon_mag) * 0.4)
                      glow_cld1 = log10(sph_rad_ave(2) * ratio_moonsun)
                  endif
              else
                  glow_cld1 = 0.
              endif

              if(moon_cond_clr .eq. 1)then
!                 Add phys term for scattering by moonlight on the clouds
                  pf_scat_moon = pf_scat(2,i,j)
                  rad_cld_moon = rad_cld_day * cloud_rad_c(2,i,j) * pf_scat_moon
                  if(rad_cld_moon .gt. 0.)then
                      glow_cld_moon = log10(rad_cld_day * cloud_rad_c(2,i,j) * pf_scat_moon)
                  else
                      glow_cld_moon = 0.
                  endif
                  glow_cld = addlogs(glow_cld1,glow_cld_moon)
              else
                  glow_cld_moon = -999.
                  glow_cld = glow_cld1
              endif

!             Note that cloud_rad_c(1) is surface night lighting of clouds (nl)
              if(cloud_rad_c(1,i,j) .gt. 0.)then
                  glow_cld_nt = log10(cloud_rad_c(1,i,j))
              else
                  glow_cld_nt = 0.
              endif
              glow_cld_c(1) = addlogs(glow_cld,glow_cld_nt+0.10)
              glow_cld_c(2) = addlogs(glow_cld,glow_cld_nt+0.00)
              glow_cld_c(3) = addlogs(glow_cld,glow_cld_nt-0.20)

!             During twilight, compare with clear sky background
!             Note that secondary scattering might also be considered in
!             early twilight away from the bright twilight arch.
!             The result creates a contrast effect for clouds vs twilight            
!             rintensity(:) = max(min(((glow_cld_c(:) -argref) * contrast + 128.),455.),0.)
!             rintensity(:) = max(min(((glow_cld_c(:) -glwmid) * contrast + 128.),455.),0.)

              if(idebug .eq. 1 .AND. r_cloud_3d(i,j) .gt. 0.)then ! cloud present
                  if(cloud_rad_c(2,i,j) .ge. 1e20)then
                      write(6,*)' WARNING: large cloud_rad_c',glow_cld_c(:),glow_cld,glow_cld_moon,glow_cld_nt
                  endif
                  write(6,51)solalt_ref,twi_alt,cloud_rad_c(1:2,i,j),r_cloud_rad(i,j),pf_scat_moon,glow_cld1 &
                            ,glow_cld_moon,glow_cld_nt,glow_cld,glow_cld_c(:),cld_rad(:)
51                format('   salt-rf-tw/cloud_rad_c/crad/pf/glow: cld1|moon|nt|cld|cldc cldrad',2f8.2,2e11.4,9f8.2,3f10.0)
              endif

              cld_rad(:) = 10. ** glow_cld_c(:) ! newly defined

          endif ! solalt_ref > twi_alt

!         Note cld_rad only newly defined for sol_alt < twi_alt
          call nl_to_RGB(cld_rad(:),glwmid,contrast & 
                        ,128.,0,cld_red,cld_grn,cld_blu)

!         Obtain brightness (glow) of clear sky
          if(sol_alt .gt. twi_0)then ! Daylight to intermediate twilight
                                     ! from skyglow routine
              if(.true.)then
                if(sol_alt .lt. -4.)then ! city lights on
                  clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + clear_rad_c_nt(:,i,j)
                  clear_rad_2nd_c(:,i,j) = clear_rad_2nd_c(:,i,j) + clear_rad_c_nt(:,i,j)
                endif

                if(clear_rad_c(2,i,j) .le. 0.)then
                  if(idebug_a(i,j) .ge. 1)then
                    write(6,53)clear_rad_c(:,i,j),i,j,alt_a(i,j),azi_a(i,j)
53                  format(' ERROR: clear_rad_c(2,i,j) <= 0.',3e12.4,2i5,2f9.2)
                  endif
                  glow_tot = 0.0
                else
                  glow_tot = log10(clear_rad_c(2,i,j)) ! + log10(clear_radf_c(2,i,j))
                endif

                if(idebug_a(i,j) .ge. 1)then
                  write(6,54)clear_rad_c(:,i,j),clear_rad_c_nt(:,i,j)
54                format(' clrradc/nt',3f12.0,3x,3f10.0)
                endif

                if((l_solar_eclipse .eqv. .true.) .and. eobs .gt. 0.9)then
                  ramp_eobs = -log10(1.-min(eobs,.9995)) ! ranges 1-3
                  ramp_eobs = -log10(1.-min(eobs,.9975)) ! test      
!                 arge = 5.0 ! reasonable avoidance of saturation     
!                 arge = 7.3 - ramp_eobs*0.55
!                 arge = 7.3 - ramp_eobs*0.7 ! for 11070[2-8]a version
                  arge = 7.3 - ramp_eobs*0.9
                  fstops = (7.3-arge) / log10(2.)
                  if(i .eq. isun .and. j .eq. jsun)then
                      write(6,55)eobs,ramp_eobs,arge,fstops
55                    format(' eobs/ramp_eobs/arge/fstops',f9.5,3f9.3)
                  endif
                  glwlow = arge
                endif

                call nl_to_RGB(clear_rad_c(:,i,j),glwmid,contrast & 
                              ,128.,0,clr_red,clr_grn,clr_blu)
              endif ! .true.

              if(idebug .eq. 1 .AND. (elong_a(i,j) .lt. 0. .or. &
                 abs(alt_a(i,j)) .le. 2.0 .or. l_solar_eclipse .eqv. .true.) )then
                  write(6,60)elong_a(i,j),alt_a(i,j),azi_a(i,j),clear_rad_c(:,i,j)
60                format(' elong / altaz / Clrrad rgb = ',3f9.2,3f14.1)
              endif

          elseif(sol_alt .ge. -16.)then ! Deep Twi from clear_rad_c array
!             call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),obs_glow_zen,patm,htmsl,clear_rad_c_nt)
              glow_nt = log10(clear_rad_c_nt(3,i,j)) ! log nL           
 
!             Gray out twilight clear_rad_c and add nt component
!             clear_rad_c(:,i,j) = clear_rad_c(3,i,j) + clear_rad_c_nt(3)
              clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + clear_rad_c_nt(:,i,j)

              if(moon_cond_clr .eq. 1)then
                  glow_moon_s = glow_moon1(i,j)
!                 clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + 10.**glow_moon_sc(i,j)
                  clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + rad_moon_sc(:,i,j)
              else
                  glow_moon_s = 1.
              endif
              glow_nt = addlogs(glow_nt,glow_moon_s)

              hue = clear_rad_c(1,i,j)
              sat = clear_rad_c(2,i,j)
              glow_twi = log10(clear_rad_c(3,i,j)) ! phys values (log nL)

              glow_tot = addlogs(glow_twi,glow_secondary_clr)
              glow_tot = addlogs(glow_tot,glow_nt)

              frac_twi  = (10.**glow_twi)           / (10.**glow_tot)
              frac_sec  = (10.**glow_secondary_clr) / (10.**glow_tot)
              frac_nt   = (10.**glow_nt)            / (10.**glow_tot)
              alt_ramp = cosd(alt_a(i,j))**40. ! hits 0.5 at ~10 deg
!             0 if either low secondary or low night lighting dominates
              sat_ramp = (1. - frac_sec*alt_ramp) * (1. - frac_nt*alt_ramp)
              sat = sat * sat_ramp         

              if(glow_stars(2,i,j) .gt. 1.0)then
!                 write(6,*)'i,j,glow_stars',i,j,glow_stars(2,i,j)
              endif
!             glow_stars(i,j) = 1.0                           ! test

              if(airmass_2_topo(i,j) .eq. 0.)then ! free of terrain
                  frac_nt = ((10.**glow_nt) / (10.**glow_tot))**1
                  frac_sat = frac_sec**3. ! (0-1 if frac_sec is high)     
!                 Has been a bit green on the anti-twi horizon for -9 deg
!                 Now trying to reduce above exponent from 6 to 3
                  hueo = hue ; sato = sat
                  hue = 1.3 * frac_sat + hue * (1. - frac_sat)
                  sat = 0.4 * frac_sat + sat * (1. - frac_sat)
                  satmax = 0.10 + abs(hue-2.0)**2 ! chromaticity curve
                  sat = min(sat,satmax)  ! restrict by elongation if needed
!                 glow_tot = addlogs(glow_tot,glow_stars(2,i,j))! add stars 
              else ! bluish twilight 2ndary scattering in front of terrain
                  glow_tot = glow_secondary_clr
                  frac_sec = ((10.**glow_secondary_clr)/(10.**glow_tot))**6
                  hueo = hue ; sato = sat
                  hue = hue * (1.-frac_sec) + 1.3 * frac_sec
                  sat = sat * (1.-frac_sec) + 0.5 * frac_sec
                  satmax = 0.10 + abs(hue-2.0)**2 ! chromaticity curve
                  sat = min(sat,satmax)
              endif

              star_ratio = 10. ** ((glow_tot - glow_twi) * 0.45)

!             if(sun_vis .eq. 1.0)then ! sun glowing below horizon
!                 glow_tot = addlogs(glow_tot,glow_sun(i,j))
!             endif

!             arg = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15             
              arg = glow_tot                                  ! experiment?            

              rintensity_floor = 0. ! 75. + sol_alt
!             rintensity_glow = max(min(((arg     -7.) * 100.),255.),rintensity_floor)
!             rintensity_glow = max(min(((arg -argref) * 100.),255.),rintensity_floor)
              rintensity_glow = max(min(((arg -argref) * contrast + 128.),255.),rintensity_floor)
!             rintensity_glow = min(rintensity_glow*star_ratio,255.)              

          else ! Night from clear_rad_c array (default flat value of glow)
            if(airmass_2_topo(i,j) .eq. 0.)then ! free of terrain
!             call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),obs_glow_zen,patm,htmsl,clear_rad_c_nt)
              hue = 1. 
              sat = 0. 
              glow_nt = log10(clear_rad_c_nt(2,i,j)) ! log nL           
              clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + clear_rad_c_nt(:,i,j)

              if(moon_cond_clr .eq. 1)then ! add moon mag condition
!               Glow from Rayleigh, no clear_rad crepuscular rays yet
!               argm = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15
                glow_moon_s = glow_moon1(i,j)          ! log nL                 
                glow_tot = addlogs(glow_nt,glow_moon_s)
!               clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + 10.**glow_moon_sc(i,j)
                clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + rad_moon_sc(:,i,j)
              else
                glow_tot = glow_nt
                glow_moon_s = 0.
              endif

!             Add in stars. Stars have a background glow of 1.0
!             glow_tot = addlogs(glow_tot,glow_stars(2,i,j))

!!            if((idebug .eq. 1 .and. moon_alt .gt. 0.) .OR. glow_stars(2,i,j) .gt. 1.0)then
!             if((idebug .eq. 1) .OR. glow_stars(2,i,j) .gt. 3.0)then
              if(idebug .eq. 1)then
                  write(6,91)i,j,idebug,elong_a(i,j),glow_nt,glow_moon_s,glow_stars(2,i,j),glow_tot,clear_rad_c(:,i,j)
91                format('   glow: elg/nt/moon/stars/tot/rad = ',3i5,f7.1,4f9.3,3f10.0)
!                 idebug = 1 ! for subsequent writing at this grid point
              endif

!             rintensity_glow = max(min(((glow_tot - 2.3) * 100.),255.),20.)
              rintensity_glow = max(min(((glow_tot - argref) * contrast + 128.),255.),20.)

            else ! terrain present (catch airglow from above)
!             call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),obs_glow_zen,patm,htmsl,clear_rad_c_nt)
              clear_rad_c(:,i,j) = clear_rad_c(:,i,j) + clear_rad_c_nt(:,i,j)

!             Add in moonglow depending of opacity of airmass to terrain
              od_2_topo = (od_atm_g * airmass_2_topo(i,j)) + aod_2_topo(i,j)
              clear_rad_c(:,i,j) = clear_rad_c(:,i,j) &
                                 + rad_moon_sc(:,i,j) * opac(od_2_topo)

            endif
          endif

          do ic = 1,nc 

!                     Rayleigh  Ozone   Mag per optical depth            
            od_atm_g = (0.1451  + .016) / 1.086
!           od_2_cloud = (od_atm_g + od_atm_a*ext_a(ic)) * airmass_2_cloud(i,j)
            od_2_cloud = ext_g(ic) * airmass_2_cloud(i,j) &
                       + ext_a(ic) * aod_2_cloud(i,j)

!           Empirical correction to account for bright clouds being visible
            cloud_visibility = exp(-od_2_cloud) ! empirical coeff

!           Use clear sky values if cloud cover is less than 0.5
            opac_cloud = r_cloud_3d(i,j)
!           opac_cloud = 0.0 ; Testing

!           0 by day and 1 by night
            ramp_cld_nt = min(max((-sol_alt/3.0),0.),1.)               
            ramp_cld_nt = ramp_cld_nt**0.7

!           Always 1 by day and goes to 0 at night near sun
            scurve_term = scurve(scurve(elong_s(i,j)/180.)) * ramp_cld_nt + 1.0 * (1.0-ramp_cld_nt)

!           Scurve can turn off cloud_visibility at night near the sun
!           When scurve=0 and cloud_visibility remains at 1
!           When scurve=1 cloud visibility is active

!           Best results for 'frac_front' may be when using gas+aero at
!           low elevation looking up and just gas at high elevation
!           looking down. It depends if we are automatically considering
!           the clear sky illumination of gas+aero for all types of
!           cloud attenuation. 'frac_front' should be relatively high
!           looking at sunlit clouds from above and should be low
!           looking at shaded clouds from below.
            if(sol_alt .gt. -4.5)then ! distant clouds fade more during day

              if(airmass_2_topo(i,j) .gt. 0.)then ! terrain
                od_2_topo = (ext_g(ic) * airmass_2_topo(i,j)) + aod_2_topo(i,j)
                frac_front = opac(od_2_cloud) / opac(od_2_topo)
              else ! fraction of air/scatterers in front of cloud                
                frac_front = opac(od_2_cloud) / opac(clr_od(ic)) ! gas+aero
              endif
              frac_front = min(frac_front,1.0)

!             Fraction of clr/lit air behind cloud
              frac_behind = 1.0 - frac_front * (1.0-ramp_cld_nt)

            else ! fade only when opposite the twilight arch
              frac_front = 0.
              frac_behind = 1.0              ! frac clr/lit air behind cloud

            endif

!           frac_cloud_2nd = opac_cloud * cloud_visibility

!           Fraction of white cloudy radiance making it to the observer
!           Fade in daytime and if opposite the twilight arch
            frac_cloud = opac_cloud * cloud_visibility**scurve_term  

            btau = 0.10 * cloud_od(i,j)
            cloud_albedo = btau / (1. + btau)

!           Consider frac_clr also in relation to fraction of airmass 
!           between the observer and the cloud contributing to Rayleigh/Mie
!           scattered light. This helps for the case of thin cirrus adding
!           to the clear sky light. Testing a new option to blend clr/cld 
!           more radiatively a la topo.
!           Note cloud may be located in front of or behind clear sky glow
!           It would be in front during twilight for single scattering.

!           Fraction of clear radiance making it to the observer
            if(sol_alt .gt. 0.0)then ! daytime
!             frac_clr = 1.0 - frac_cloud
              frac_clr = frac_front + (1.-frac_front) * (1.-opac_cloud)
            else
!             frac_clr = 1.0 - frac_cloud*frac_behind ! clear sky behind cloud                         
!             frac_clr_2nd = clear_rad_2nd_c(ic,i,j) / (clear_rad_c(ic,i,j) + clear_rad_2nd_c(ic,i,j))
              if(clear_rad_c(ic,i,j) .gt. 0.)then
                frac_clr_2nd = clear_rad_2nd_c(ic,i,j) / clear_rad_c(ic,i,j)
              else
                frac_clr_2nd = 0.
              endif
              arg = ((1. - ramp_cld_nt) * (1. - frac_clr_2nd)) + (1. * frac_clr_2nd)
              frac_front = frac_front * arg
              frac_clr = frac_front + (1.-frac_front) * (1.-opac_cloud)
            endif

!           Why is frac_clr slightly high on a cloudy day at the surface?
!           if(sol_alt .gt. 0.0)then ! daytime
!             frac_clr = frac_clr * clear_radf_c(ic,i,j)
!           endif

!           Clear sky radiance should be modulated by the cloud if it has to
!           go through the cloud. It can be added to the cloud if the cloud
!           is transparent, or if the illuminated air is in front of the cloud.
            radf_through_cloud = 1.0 - (r_cloud_3d(i,j)*scurve(patm/patm_sfc))

!           This should be turned on when we are looking at a dark cloud
!           lying in front of the sun and the sub-cloud air (particularly
!           the aerosol component) is unilluminated. Thus the aureole from 
!           clear air scattering should not show up projected in front of 
!           a cloud, if the line of sight from the observer to the cloud
!           is in shadow.
!           TRUE option works on a cloudy day at any altitude.         
!           It fails at the surface on a sunny day.
!           FALSE option works on a sunny day at any altitude. 
!           It used to fail at the surface on a cloudy day.
            if(sol_alt .gt. 0.0 .and. .false.)then ! daytime
              sky_rad(ic) = clear_rad_c(ic,i,j) * frac_clr * radf_through_cloud &
                          + cld_rad(ic) * frac_cloud
            else ! always
              sky_rad(ic) = clear_rad_c(ic,i,j) * frac_clr &
                          + cld_rad(ic) * frac_cloud
            endif

          enddo ! ic

          call nl_to_RGB(sky_rad(:),glwmid,contrast & 
                        ,128.,0,sky_rgb(0,I,J),sky_rgb(1,I,J),sky_rgb(2,I,J))
!         sky_rgb(:,I,J) = min(sky_rgb(:,I,J),255.)

!         if(idebug_a(i,j) .eq. 1 .AND. alt_a(i,j) .le. 2.0)then
          if(idebug_a(i,j) .eq. 1 .AND. alt_a(i,j) .le. 90.0)then
            write(6,96)od_2_topo,od_2_cloud,clr_od(2),aod_2_cloud(i,j),cloud_albedo
96          format(' od2tpo/od2cld/odclr/aod2cld/cldalb',5f10.4)
            write(6,97)frac_front,frac_behind,r_cloud_3d(i,j),clear_radf_c(2,i,j),frac_clr,frac_clr_2nd,frac_cloud,scurve_term,clear_rad_c(:,i,j),cld_rad,sky_rad
97          format(' ffnt/fbhd/rcld/radf/fclr/fclr2/fcld/scrv/clrrd/cldrd/skyrd',8f7.3,3f13.0,2x,3f13.0,2x,3f13.0,2x,3i4)
          endif

!         Use topo value if airmass to topo > 0
          if(airmass_2_topo(i,j) .gt. 0.)then ! looking at terrain

!             The sky RGB values should already include reductions from
!             cloud/terrain shadowing and limited distance to the topography
!             This may however consider aerosol brightness and not gas
!             component?
              if(sol_alt .le. twi_alt)then 
                od_2_topo = 0. !(od_atm_g + od_atm_a) * airmass_2_topo(i,j)
              else ! eventually use clear_rad influenced by topo?
                od_2_topo = (od_atm_g * airmass_2_topo(i,j)) + aod_2_topo(i,j)
!               od_2_topo = (od_atm_g * airmass_2_topo(i,j)) + aod_ill(i,j)
              endif

              topo_visibility = trans(+1.00*od_2_topo)                    

              if(airmass_2_cloud(i,j) .gt. 0. .AND. airmass_2_cloud(i,j) .lt. airmass_2_topo(i,j)) then
                  topo_visibility = topo_visibility * (1.0 - r_cloud_3d(i,j))
              endif

!             If the view altitude is near -90, we can consider the topo as
!             adding to the cloud via forward scattering
              if(alt_a(i,j) .lt. 0.)then
                  rlnd_fwsc = -sind(alt_a(i,j))
                  topo_rad = 1.0 - cloud_albedo
                  topo_visibility = topo_visibility * (1.0-rlnd_fwsc) &
                                  + topo_rad        * rlnd_fwsc
              endif

              if(sol_alt .gt. twi_alt)then 
!                 Daytime assume topo is lit by sunlight (W/m**2)
!                 Add topo brightness to line of sight sky brightness in radiation space

!                 Using secondary cloud glow as lower bound during twilight
!                 if(sol_alt .gt. 0.)then ! floor of SWI/GTI field (W/m**2)
!                     rindirect = 25.
!                 else
!                     rindirect = 1300. * (10.**glow_secondary_cld / rad_cld_day)
!                 endif

!                 topo_gti_frac = max(topo_swi(i,j),rindirect) / 1300. 
!                 topo_gti_frac =     topo_gti(i,j)            / 1300. 
                  rtopo_red = 2. * gtic(1,i,j) * topo_albedo(1,i,j) * pf_land(1,i,j)
                  rtopo_grn = 2. * gtic(2,i,j) * topo_albedo(2,i,j) * pf_land(2,i,j) 
                  rtopo_blu = 2. * gtic(3,i,j) * topo_albedo(3,i,j) * pf_land(3,i,j) 
                  sky_frac_topo = 1.00               ! hopefully temporary
                  sky_frac_aero = 1.00               ! hopefully temporary
                  red_rad = day_int * rtopo_red*topo_visibility*sky_frac_topo
!                         + counts_to_rad(sky_rgb(0,I,J)) * sky_frac_aero 
                  grn_rad = day_int * rtopo_grn*topo_visibility*sky_frac_topo
!                         + counts_to_rad(sky_rgb(1,I,J)) * sky_frac_aero
                  blu_rad = day_int * rtopo_blu*topo_visibility*sky_frac_topo
!                         + counts_to_rad(sky_rgb(2,I,J)) * sky_frac_aero 
                  if(idebug .eq. 1)then
                    write(6,98)rtopo_grn,topo_gti(i,j),gtic(2,i,j) &
                              ,topo_albedo(2,i,j),topo_solalt(i,j) &
                              ,dist_2_topo(i,j) &
!                             ,nint(sky_rgb(:,i,j)),red_rad,grn_rad,blu_rad
                              ,red_rad,grn_rad,blu_rad,sky_rad(:)
98                  format( &
                        ' rtopo/gti/gtic/alb/tsalt/dst/trad/srad', &
                           f9.3,f9.1,f9.4,f9.3,f9.2,f10.0,2x,3f12.0,2x,3f14.0)
                  endif

                  sky_rad(1) = sky_rad(1) + red_rad
                  sky_rad(2) = sky_rad(2) + grn_rad
                  sky_rad(3) = sky_rad(3) + blu_rad
 
                  call nl_to_RGB(sky_rad(:),glwmid,contrast & 
                          ,128.,0,sky_rgb(0,I,J),sky_rgb(1,I,J),sky_rgb(2,I,J))

              else ! sol_alt < twi_alt
!                 Nighttime assume topo is lit by city lights (nL)
!                 Moonlight can now be added so we have both
!                 emitted and reflected light.
!                 Note: this can be done more in rad space and airglow
!                 can be added if viewpoint is in space using
!                 'clear_rad_c', 'sky_rad', or 'clear_rad_c_nt'

!                 topo_gti_frac = (max(topo_gti(i,j),001.) / 5000.) ** 0.45
!                 Get city + moon colors via 'topo_gtic'?
                  if(.true.)then ! experiment (relative solar to nL)
                    rtopo_red = day_int * (emic(1,i,j) & 
                              + 2. * gtic(1,i,j) * topo_albedo(1,i,j) )
                    rtopo_grn = day_int * (emic(2,i,j) &
                              + 2. * gtic(2,i,j) * topo_albedo(2,i,j) )
                    rtopo_blu = day_int * (emic(3,i,j) &
                              + 2. * gtic(3,i,j) * topo_albedo(3,i,j) )
                  else ! nL (original city lights only)
                    rtopo_red = topo_gti(i,j) * 1.5
                    rtopo_grn = topo_gti(i,j) * 1.0
                    rtopo_blu = topo_gti(i,j) * 0.5
                  endif

                  topo_viseff = topo_visibility * patm ! allow airglow

                  red_rad = rtopo_red*topo_visibility 
                  grn_rad = rtopo_grn*topo_visibility 
                  blu_rad = rtopo_blu*topo_visibility 

                  if(idebug .eq. 1)then
                    write(6,99)rtopo_grn,topo_gti(i,j),gtic(2,i,j) &
                              ,emic(2,i,j) &
                              ,topo_albedo(2,i,j),topo_solalt(i,j) &
                              ,dist_2_topo(i,j) &
!                             ,nint(sky_rgb(:,i,j)),red_rad,grn_rad,blu_rad
                              ,red_rad,grn_rad,blu_rad,sky_rad(:)
99                  format( &
                        ' rtopo/gti/gtic/alb/tsalt/dst/trad/srad', &
                           f9.0,f9.4,2f10.7,f9.3,f9.2,f10.0,2x,3f12.0,2x,3f14.0)
                  endif

                  sky_rad(1) = sky_rad(1)*(1.0-topo_viseff) + red_rad
                  sky_rad(2) = sky_rad(2)*(1.0-topo_viseff) + grn_rad
                  sky_rad(3) = sky_rad(3)*(1.0-topo_viseff) + blu_rad
 
                  call nl_to_RGB(sky_rad(:),glwmid,contrast & 
                          ,128.,0,sky_rgb(0,I,J),sky_rgb(1,I,J),sky_rgb(2,I,J))

              endif ! sol_alt

              rad_sun = 0.

          else ! not looking at terrain
              od_2_topo = 0.

!             Any refraction can be done up front in 'get_glow_obj'.
!             This can theoretically be done here multispectrally 
!             using newly passed in 'od_slant_g_a'.

              if(sol_alt .gt. -3. .or. solalt_limb_true .gt. 0.)then
!             if(.false.)then
!               Add sun to existing sky radiance (blue extinction)
                rad_sun_r = 10.**glow_sun(i,j)  &
                          * trans(cloud_od(i,j) + clr_od(1)) 
                rad_sun_g = 10.**glow_sun(i,j)  &
                          * trans(cloud_od(i,j) + clr_od(2)) 
                rad_sun_b = 10.**glow_sun(i,j)  &
                          * trans(cloud_od(i,j) + clr_od(3)) 

                sky_rad(1) = sky_rad(1) + rad_sun_r
                sky_rad(2) = sky_rad(2) + rad_sun_g
                sky_rad(3) = sky_rad(3) + rad_sun_b    
              endif
                
!             Add moon/stars to existing sky radiance (blue extinction?)
              if(sol_alt .gt. 0. .and. (l_solar_eclipse .eqv. .false.))then
                rad_moon = 10.**glow_moon(i,j)  &
                          * trans(cloud_od(i,j) + aod_slant(i,j)) 
                rad_moon_r = 10.**glow_sun(i,j)  &
                           * trans(cloud_od(i,j) + clr_od(1)) 
                rad_moon_g = 10.**glow_sun(i,j)  &
                           * trans(cloud_od(i,j) + clr_od(2)) 
                rad_moon_b = 10.**glow_sun(i,j)  &
                           * trans(cloud_od(i,j) + clr_od(3)) 

                sky_rad(1) = sky_rad(1) + rad_moon_r
                sky_rad(2) = sky_rad(2) + rad_moon_g
                sky_rad(3) = sky_rad(3) + rad_moon_b
              else ! stars that also contain the moon
!               if(idebug .eq. 1)then
!                 write(6,*)'rgb before stars/moon',sky_rgb(:,I,J)
!               endif

!               Extinction from clouds, gas and aerosols
                do ic = 1,nc
                  if(cloud_od(i,j) .gt. 0.)then
                    transterm(ic) = trans(cloud_od(i,j) + clr_od(ic))
                  else
!                   transterm(ic) = 1.0
                    transterm(ic) = trans(clr_od(ic))
                  endif
                enddo ! ic

                rad_stars_r = 10.**glow_stars(1,i,j) * transterm(1)
                rad_stars_g = 10.**glow_stars(2,i,j) * transterm(2)
                rad_stars_b = 10.**glow_stars(3,i,j) * transterm(3)

                sky_rad(1) = sky_rad(1) + rad_stars_r
                sky_rad(2) = sky_rad(2) + rad_stars_g
                sky_rad(3) = sky_rad(3) + rad_stars_b

                if(idebug .eq. 1)then
                  write(6,100)glow_stars(:,i,j),cloud_od(i,j),aod_slant(i,j) &
                      ,transterm(:),rad_stars_r,rad_stars_g,rad_stars_b
100               format(' glow/cod/aod/trans/rad_stars',3f9.2,5f8.3,3f12.1)
                endif
              endif

              call nl_to_RGB(sky_rad(:),glwmid,contrast & 
                        ,128.,0,sky_rgb(0,I,J),sky_rgb(1,I,J),sky_rgb(2,I,J))

!             if(idebug .eq. 1)then
!               write(6,100)sky_rad,sky_rgb(:,I,J)
!100            format('skyrad/rgb after moon/stars',3f10.0,f9.2)
!             endif

          endif ! looking at terrain

          if(idebug .eq. 1)then
              if(alt_a(i,j) .gt. 1.7 .and. alt_a(i,j) .le. 2.1)then
                call nl_to_RGB(sky_rad(:),glwmid,contrast & 
                          ,128.,0,rtotal,gtotal,btotal)                        
                write(6,105)rtotal,gtotal,btotal,sky_rad(:)
105             format(' total rgb/skyrad should be',3f9.2,94x,3f14.0)
              endif

              rmaglim = b_to_maglim(10.**glow_tot)
              call apply_rel_extinction(rmaglim,alt_a(i,j),od_atm_g+od_atm_a)

              if(i .eq. ni)then
                  write(6,*)' ******** zenith location ************************ od'
                  write(6,111)clear_rad_c(:,i,j),nint(clr_red),nint(clr_grn),nint(clr_blu)
111               format('clrrad/RGB = ',3f12.0,3i4)
              endif
              if(abs(elong_a(i,j) - 90.) .le. 0.5)then
                  write(6,*)' ******** 90 elong location ********************** od'
                  write(6,111)clear_rad_c(:,i,j),nint(clr_red),nint(clr_grn),nint(clr_blu)
              endif
              if(i .eq. isun .and. j .eq. jsun)then
!             if(elong_a(i,j) .le. 0.5)then
                  write(6,*)' ******** solar location ************************* od'
                  write(6,112)glow_sun(i,j),10.**glow_sun(i,j),alt_a(i,j) &
                             ,rad_sun_r,rad_sun_g,rad_sun_b &
                             ,trans(cloud_od(i,j) + aod_slant(i,j)),glwlow 
112               format('glow_sun/alt/rad_sun/trans',f9.4,e16.6,f9.4,3e16.6,f11.8,f9.4)
                  write(6,*)'od_g_slant = ',od_g_slant_a(:,i)
                  write(6,*)'od_o_slant = ',od_o_slant_a(:,i)
                  write(6,*)'od_a_slant = ',od_a_slant_a(:,i)
                  write(6,*)'clr_od     = ',clr_od(:)
              endif
              if(alt_a(i,j) .eq. -90.)then
                  write(6,*)' ******** nadir location ************************* od'
              endif

!             if(sol_alt .ge. -2.0)then      ! daylight
              if(sol_alt .ge. twi_alt)then   ! daylight
                  write(6,116)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat(2,i,j),r_cloud_3d(i,j),cloud_od(i,j),cloud_od_sp(i,j,:),bkscat_alb(i,j) &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                      ,topo_gti(i,j),topo_albedo(1,i,j),aod_ill(i,j),od_2_topo,topo_visibility,cloud_visibility,rintensity_glow &
                      ,nint(sky_rgb(:,i,j)),nint(cld_red),nint(cld_grn),nint(cld_blu)
              elseif(sol_alt .ge. -16.)then ! twilight
                  write(6,117)i,j,alt_a(i,j),azi_a(i,j),elong_s(i,j) & 
                      ,pf_scat(2,i,j),r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb(i,j) &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j) &
                      ,rintensity(1),glow_cld_nt,glow_cld,glow_nt,glow_twi &
                      ,glow_secondary_clr,glow_tot,rmaglim,cloud_visibility &
                      ,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c(:,i,j),nint(clr_red),nint(clr_grn),nint(clr_blu)
              else ! night
                  write(6,118)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat(2,i,j),r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb(i,j) &
!                     ,frac_cloud,airmass_2_cloud(i,j),cloud_rad_c(2,i,j)*1e3,rintensity(1) &
                      ,frac_cloud,airmass_2_cloud(i,j),cloud_rad_c(1,i,j)*1e3,rintensity(1) &
                      ,glow_cld_nt,glow_cld_moon,glow_cld,glow_secondary_clr,rmaglim &
                      ,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c_nt(:,i,j)
              endif
116           format(2i5,3f6.1,f9.3,f11.6,2f7.2,3f6.2,3f8.3,f8.4,f7.1,f9.5,f9.1,f8.3,2f8.5,2f8.3,f9.2,2x,3i4,' cldrgb',1x,3i4)
117           format(2i5,3f9.2,f9.3,f11.6,4f9.3,f9.4,f7.1,f9.3,f8.1,6f7.3,f9.2,2x,3i4,' clrrad',3f10.0,3i4)
118           format(2i5,3f9.2,f9.3,f11.6,4f9.3,f9.6,f7.1,f9.3,f9.3,4f9.3,f9.2,2x,3i4,' clrrad',3f8.2)
          endif

        enddo ! i
        enddo ! j

        write(6,*)' max RGB = ',maxval(sky_rgb(0,:,:)) &
        ,maxval(sky_rgb(1,:,:)),maxval(sky_rgb(2,:,:))

!       Consider max RGB values where elong > 5 degrees

!       Add bounds to rgb values
        do j = 1,nj
        do i = 1,ni
          final_scaling = 1.0
          sky_rgb(:,i,j) = max(min(sky_rgb(:,i,j)*final_scaling,255.),0.)
        enddo ! i
        enddo ! j

!       Convert nl values to spectral radiance
!       sky_sprad = f(sky_cyl_nl)

        return
        end

        subroutine rgb_to_rad(count_r,count_g,count_b &
                             ,glwlow,contrast,offset,rad_r,rad_g,rad_b)

        counts_to_rad(counts) = 10.**(((counts-offset)/contrast)+glwlow)

        if(count_g .le. 0.)then
           rad_r = 0.; rad_g = 0.; rad_b = 0.
           return
        endif

        gamma = 2.2

        rad_g = counts_to_rad(count_g)

        rad_rog = (count_r / count_g) ** gamma
        rad_r = rad_g * rad_rog

        rad_bog = (count_b / count_g) ** gamma
        rad_b = rad_g * rad_bog

        return
        end

        subroutine rad_to_rgb(rad_r,rad_g,rad_b,glwlow,contrast,offset &
                             ,count_r,count_g,count_b)

        rad_to_counts(rad) = (log10(rad)-glwlow)*contrast + offset

        if(rad_g .le. 0.)then
           count_r = 0.; count_g = 0.; count_b = 0.
           return
        endif

        gamma = 2.2

        count_g = rad_to_counts(rad_g)

        count_rog = (rad_r / rad_g) ** (1./gamma)
        count_r = count_g * count_rog

        count_bog = (rad_b / rad_g) ** (1./gamma)
        count_b = count_g * count_bog

        return
        end
