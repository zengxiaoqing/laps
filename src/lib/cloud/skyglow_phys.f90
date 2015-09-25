

        subroutine skyglow_phys(ialt_start,ialt_end,ialt_delt      &! I
                   ,jazi_start,jazi_end,jazi_delt,azi_scale        &! I
                   ,minalt,maxalt,minazi,maxazi,idebug_a           &! I
                   ,sol_alt,sol_azi,view_alt,view_az,twi_0,twi_alt &! I
                   ,isolalt_lo,isolalt_hi,topo_solalt,trace_solalt &! I
                   ,earth_radius,patm,aod_vrt,aod_ray,aod_ray_dir  &! I
                   ,aod_ref,aero_scaleht,dist_2_topo               &! I
                   ,htmsl,redp_lvl,horz_dep                        &! I
                   ,aod_ill,aod_2_topo,aod_tot                     &! I
                   ,l_solar_eclipse,i4time,rlat,rlon               &! I
                   ,clear_radf_c,ag_2d                             &! I
                   ,od_g_slant_a,od_o_slant_a,od_a_slant_a,ext_a   &! O
                   ,clear_rad_c,sky_rad_ave,elong           )       ! O/I

!       Sky glow with solar altitude > 0

        include 'trigd.inc'

        use mem_namelist, ONLY: fcterm, aod_bin, aod_asy, r_missing_data
        use mem_allsky, ONLY: aod_ill_opac,aod_ill_opac_potl            ! I

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
        curvat(hdst,radius) = hdst**2 / (2. * radius)
        expbh(x) = exp(min(x,+80.))

        real linecyl,linecylp
        linecylp(xcos,ycos,x1,r) = (-(2.*xcos*x1) + sqrt((2.*xcos*x1)**2 - 4. * (xcos**2 + ycos**2) * (x1**2 - r**2))) &
                                 /  (2. * (xcos**2 + ycos**2))

        include 'rad.inc'

        cosp(b,x) = (1. + b/2.) * cosd(x/2)**b
        dhg_gterm(pha,g) = 1. / (1. + g**2 - 2.*g*cosd(pha))**1.5
        dhg_fterm(pha,g) = (3. * ((cosd(pha))**2) - 1.) / (2. * (1.+g**2)**1.5)
        dhg(pha,g,f) = (1.-g**2) * (dhg_gterm(pha,g) + f * dhg_fterm(pha,g))
        dhg2(pha,f,p) = (1.-p) * dhg(pha,g1,f) + p * dhg(pha,g2,f)

        real twi_trans_c(nc)           ! transmissivity

        real hg2d(nc), alphav_g, alphav_a ! hg2(nc)
        real srcdir_90(nc),srcdir(nc),clear_int_c(nc)
        real sumi_gc(nc),sumi_ac(nc)
        real sumi_gct(nc),sumi_act(nc)
        real srcdir_a(nc,minazi:maxazi)
        real sumi_g_a(nc,minazi:maxazi),sumi_a_a(nc,minazi:maxazi)
        real ecl_intd(nc),ecl_dir_rat(nc),ecl_scat(nc),sky_rad_ave(nc)

        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)  ! clear sky illumination
        real clear_rad_c0(nc,minalt:maxalt,minazi:maxazi) ! clear sky (no eclipse)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)! integrated
               ! fraction of air molecules illuminated by the sun along
               ! line of sight (consider Earth's shadow + clouds)
               ! topo and multiple scattering is accounted for
        real ag_2d(minalt:maxalt,minazi:maxazi) ! gas airmass (topo/notopo)
        real aod_ill(minalt:maxalt,minazi:maxazi) ! aerosol illuminated
                                      ! optical depth (slant - topo/notopo)
        real aod_tot(minalt:maxalt,minazi:maxazi) ! aerosol illuminated
                                      ! optical depth (slant - in domain)
        real aod_2_topo(minalt:maxalt,minazi:maxazi) ! slant
        real aod_ray(minalt:maxalt,minazi:maxazi) ! aerosol optical 
                                      ! depth (zenithal) may be adjusted
                                      ! for illumination
        real aod_ray_dir(minalt:maxalt,minazi:maxazi) ! aerosol optical 
                                      ! depth (zenithal) may be adjusted
                                      ! for direct illumination / topo
        real patm_ray(minalt:maxalt,minazi:maxazi)! effective patm adjusted
                                      ! for illumination
        real dist_2_topo(minalt:maxalt,minazi:maxazi)
        real topo_solalt(minalt:maxalt,minazi:maxazi)
        real trace_solalt(minalt:maxalt,minazi:maxazi)
        real elong(minalt:maxalt,minazi:maxazi)
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

        real od_g_slant_a(nc,minalt:maxalt) ! use for sun/moon/star attenuation
        real od_o_slant_a(nc,minalt:maxalt)
        real od_a_slant_a(nc,minalt:maxalt)

        parameter (nsc = 3)
        real aod_asy_eff(nsc) ! ,nc (add nc dimension throughout)

        logical l_solar_eclipse, l_dlow, l_dlow2, l_topo_a(minazi:maxazi)

        parameter (nsteps = 100000)
        parameter (nopac = 10)

        real distecl(nopac,nc)
        real tausum_a(nsteps)
        real taumid(nopac),opacmid(nopac),eobsc_a(nopac,nc)
        real eobsc(nc,minazi:maxazi),eobsc_sum(nc),emag_a(minazi:maxazi)
        real sky_rad_scat(nc,minalt:maxalt,minazi:maxazi)            
!       real sky_rad_scata(minazi:maxazi)

!       Airmasses relative to zenith for observer ht
        real ags_a(isolalt_lo:isolalt_hi)
        real aas_a(isolalt_lo:isolalt_hi)
        real aos_a(isolalt_lo:isolalt_hi)

        icd = 1

        eobsc(:,:) = 0. ! initialize
        sky_rad_ave = r_missing_data
        sfc_alb = 0.15  ! pass this in and account for snow cover?
        ssa = 0.90

        aod_bin(1) = .000
        aod_bin(2) = .987
        aod_bin(3) = .013

        aod_asy(1) =  .75
        aod_asy(2) =  .70
        aod_asy(3) = -.65

!       Now read via namelist
!       fcterm = 0.09 ! range from 0.00 to 0.09 (large aerosol population)
!                     ! peak phase function of 20 to 110

        patm_ray = patm

        aero_refht = redp_lvl

        ramp_fc_nt = min(max((-sol_alt/3.0),0.),1.)               
        fcterm2 = fcterm * (1.0 - ramp_fc_nt)
        angstrom_exp_a = 2.4 - (fcterm * 15.)
        angstrom_exp_ha = 2.0

        if(sol_alt .gt. 0. .or. .true.)then
            write(6,*)' skyglow_phys: i4time is ',i4time,l_solar_eclipse
            write(6,*)' aod_bin = ',aod_bin
            write(6,*)' aod_asy = ',aod_asy
            write(6,*)' htmsl / aero_refht = ',htmsl,aero_refht
            write(6,*)' aod_vrt = ',aod_vrt
            write(6,*)' aod_ref = ',aod_ref
            write(6,*)' aod_ray max = ',maxval(aod_ray)
            write(6,*)' fcterm/fcterm2 = ',fcterm,fcterm2
!           write(6,*)' patm_ray max = ',maxval(patm_ray)
            write(6,*)' angstrom_exp_a = ',angstrom_exp_a
        endif

        do ic = 1,nc
            ext_a(ic)  = (wa(ic)/.55)**(-angstrom_exp_a)
            ext_ha(ic) = (wa(ic)/.55)**(-angstrom_exp_ha)
            write(6,*)' ic/wa/ext_a ',ic,wa(ic),ext_a(ic)
        enddo ! ic

!       Obtain reference values of source term
        call get_airmass(90.,htmsl,patm &           ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,1 &          ! I
                         ,ag_90,ao_90,aa_90)        ! O

!       Obtain reference values for sun             
        call get_airmass(sol_alt,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,1 &          ! I
                         ,ag_s,ao_s,aa_s)           ! O

!       Set this for range of values of solar altitude
        if(.false.)then
          ags_a = ag_s/ag_90
          aas_a = aa_s/aa_90
        else
          patm_refht = ztopsa(aero_refht)/1013.25
          do isolalt = isolalt_lo,isolalt_hi
            sol_alt_a = sol_alt + float(isolalt)
            if(sol_alt_a .gt. 0.)then
              call get_airmass(sol_alt_a,0.,1. &                     ! I
                              ,aero_refht,aero_scaleht &             ! I
                              ,earth_radius,1 &                      ! I
                              ,ags_a(isolalt),aos_a(isolalt),aa_dum) ! O

              call get_airmass(sol_alt_a,aero_refht,patm_refht &     ! I
                              ,aero_refht,aero_scaleht &             ! I
                              ,earth_radius,1 &                      ! I
                              ,ag_dum,ao_dum,aas_a(isolalt))         ! O
              write(6,51)isolalt,sol_alt_a,ags_a(isolalt),aos_a(isolalt) &
                                          ,aas_a(isolalt)
51            format('isolalt/solalt/ags_a/aos_a/aas_a',i3,f9.2,3f9.4)

            else ! sun below horizon, elevate to shadow height (eg132396.)
              ht_shadow = curvat(-sol_alt_a*110000.,earth_radius)
              ht_shadow = ht_shadow + (1000.-sol_alt_a*200.) ! buffer
              patm_shadow = ztopsa(ht_shadow) / 1013.25
              patmo3_shadow = patm_o3(ht_shadow)
              call get_airmass(sol_alt_a,ht_shadow,patm_shadow &     ! I
                              ,aero_refht,aero_scaleht &             ! I
                              ,earth_radius,1 &                      ! I
                              ,ag_sz,ao_sz,aa_sz)                    ! O
              ags_a(isolalt) = ag_sz / patm_shadow
              if(patmo3_shadow .gt. 0.)then
                aos_a(isolalt) = ao_sz / patmo3_shadow
              else
                aos_a(isolalt) = 0.
              endif
              aas_a(isolalt) = aa_sz &
                             * expbh((ht_shadow-aero_refht)/aero_scaleht)
              write(6,52)isolalt,sol_alt_a,ht_shadow,patm_shadow &
                ,patmo3_shadow,ags_a(isolalt),aos_a(isolalt),aas_a(isolalt)
52            format('isolalt/solalt/ht/patm/patmo3/ags_a/aos_a/aas_a',i3 &
                                                   ,f9.2,f9.1,2f9.5,3f15.1)
            endif ! solalt > 0

          enddo 
        endif

        if(sol_alt .gt. twi_0)then ! formerly 0.
          write(6,*)' Obtain reference values of source term'
          write(6,*)'ag_90/aa_90 = ',ag_90,aa_90
          write(6,*)'ag_s/aa_s = ',ag_s,aa_s

          do ic = 1,nc
            od_g_vert = ext_g(ic) * patm
            od_o_vert = ext_o(ic) * patm_o3(htmsl)
            od_a_vert = aod_vrt * ext_a(ic)
            if(ic .eq. icd)then
              idebug = 1
            else
              idebug = 0
            endif
            ssa90 = 1.0
            aa_90_o_aa_90 = 1.
            if(aa_90 .gt. 0.)then
                aa_s_o_aa_90 = aa_s/aa_90
            else
                aa_s_o_aa_90 = 0.
            endif
            call get_clr_src_dir(sol_alt,90.,ext_g(ic), &
                od_g_vert,od_a_vert,ext_ha(ic), &
                htmsl,ssa90,ag_90/ag_90,ao_90,aa_90_o_aa_90, &
                aod_ref,aero_refht,aero_scaleht,ag_s/ag_90,aa_s_o_aa_90,ic,idebug, &
                srcdir_90(ic),sumi_gc(ic),sumi_ac(ic),opac_slant, &
                nsteps,ds,tausum_a)
            write(6,*)' Returning with srcdir_90 of ',srcdir_90(ic)
          enddo ! ic
        endif

        do ialt = ialt_start,ialt_end
  
         altray = view_alt(ialt,jazi_start)
         viewalt_limb_true = altray - (-horz_dep)
         if(altray .ge. 0.)then
           htmin_view = htmsl
         else
           htmin_view = htminf(htmsl,altray,earth_radius)
         endif

         call get_airmass(altray,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,0 &          ! I
                         ,ag,ao,aa)                 ! O
         write(6,63)altray,ag,ao,aa
63       format(' returned from get_airmass with alt/ag/ao/aa = ',f9.0,3e12.4)

!        od_g_slant_a(:,ialt) = ext_g(:) * patm * ag/ag90 ! for export
         od_g_slant_a(:,ialt) = ext_g(:) * ag             ! for export
         od_o_slant_a(:,ialt) = (o3_du/300.) * ext_o(:) * ao
!        if(aa_90 .gt. 0.)then
!          od_a_slant_a(:,ialt) = aod_vert * ext_a(:) * aa/aa_90
!        else
!          od_a_slant_a(:,ialt) = 0.
!        endif
         od_a_slant_a(:,ialt) = aod_ref * aa * ext_a(:)

!        do ic = 1,nc
!          od_g_slant(ic,ialt) = ext_g(ic) * patm * ag/ag_90
!        enddo ! ic

!        Set sparse azimuth array for topo calls
         if(altray .ge. 0.)then
           jazi_delt_topo = 1
         elseif(altray .eq. -90.)then
           jazi_delt_topo = maxazi-minazi
         else
           jazi_delt_topo = int(1./cosd(altray))
         endif
         l_topo_a = .false.
         l_topo_a(minazi:maxazi:jazi_delt_topo) = .true.
!        write(6,*)'topo az delt at alt:',altray,jazi_delt_topo

!        Determine src term and ratio with zenith value
         if(sol_alt .gt. 0.)then
           do ic = 1,nc
             od_g_vert = ext_g(ic) * patm
             od_o_vert = ext_o(ic) * patm_o3(htmsl)
             od_a_vert = aod_vrt * ext_a(ic)
!            if((l_solar_eclipse .eqv. .true.) .AND. ic .eq. icd)then
             if((ic .eq. icd .and. altray .eq. nint(altray)) .OR. &
                      altray .eq. 0. .OR. altray .eq. 90.)then
               idebug = 1
               write(6,*)' alt/ag/aa =',altray,ag,aa
               write(6,*)' call get_clr_src_dir for altitude ',altray
             else
               idebug = 0
             endif
             if(aa_90 .gt. 0.)then
                 aa_o_aa_90 = aa/aa_90
                 aa_s_o_aa_90 = aa_s/aa_90
             else
                 aa_o_aa_90 = 0.
                 aa_s_o_aa_90 = 0.
             endif
             call get_clr_src_dir(sol_alt,altray,ext_g(ic), &
                od_g_vert,od_a_vert,ext_ha(ic), &
                htmsl,ssa,ag/ag_90,ao,aa_o_aa_90, &
                aod_ref,aero_refht,aero_scaleht, &
                ag_s/ag_90,aa_s_o_aa_90,ic,idebug, &
                srcdir(ic),sumi_gc(ic),sumi_ac(ic),opac_slant, &
                nsteps,ds,tausum_a)

             if(l_solar_eclipse .eqv. .true.)then
               do iopac = 1,nopac
                 opacmid(iopac) = opac_slant * (float(iopac) - 0.5) / float(nopac)
                 taumid(iopac) = -log(1.0 - opacmid(iopac))              
!                taumid(iopac) = 1.0 ! align with earlier test
                 do i = 1,nsteps
                   sbar = (float(i)-0.5) * ds
                   if(tausum_a(i) .lt. taumid(iopac))distecl(iopac,ic) = sbar
                   if(iopac .eq. nopac .and. ic .eq. icd .and. ialt .eq. ialt_start .and. jazi .eq. jazi_start)then
                     write(6,64)i,tausum_a(i),taumid(iopac),ds,sbar,distecl(iopac,ic)
64                   format('ic/i/tausum/mid/ds/sbar/dst',2i5,3f9.3,2f9.1)
                   endif
                 enddo ! i 
               enddo ! iopac         

               if(ic .eq. ic .AND. (altray .eq. nint(altray) .OR. altray .le. 20.) )then
                write(6,65)ic,altray,srcdir(ic),srcdir(ic)/srcdir_90(ic) &
                          ,opac_slant,opacmid(1),opacmid(nopac),taumid(1),taumid(nopac),distecl(1,ic),distecl(nopac,ic)
65              format( &
                   ' ic/alt/srcdir/ratio/opacsl/opacmid/taumid/distecl:',i3,f7.1,2f9.3,f7.3,2(1x,2f7.3),2f9.1)
               endif

             else ! l_solar_eclipse = F
               if(ic .eq. icd .AND. (altray .eq. nint(altray) .OR. altray .le. 20.) )then
                write(6,66)altray,srcdir(ic),srcdir(ic)/srcdir_90(ic)
66              format(' alt/srcdir/ratio:',3f9.3)
               endif

             endif ! l_solar_eclipse
           enddo ! ic
         endif ! solalt > 0

!        Determine aerosol multiple scattering order
!        altscat = 1.00 * altray + 0.00 * sol_alt
!        altscat = max(altray,sol_alt)
         altscat = sqrt(0.5 * (altray**2 + sol_alt**2))
         call get_airmass(altscat,htmsl,patm &      ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,0 &          ! I
                         ,agdum,aodum,aascat)       ! O
!        scatter_order = 1.0 ! f(altray,sol_alt)
         scatter_order = max(aod_vrt*aascat,1.0)**1.0 ! 0.5-1.5
         scatter_order_t = 1.0
         write(6,*)' altscat/aascat/sco',altscat,aascat,scatter_order

!        Determine effective asymmetry parameter from multiple scattering
!        http://www.arm.gov/publications/proceedings/conf15/extended_abs/sakerin_sm.pdf
         do isc = 1,nsc
           if(aod_asy(isc) .gt. 0.)then
              aod_asy_eff(isc) = aod_asy(isc) ** scatter_order
           elseif(aod_asy(isc) .lt. 0.)then
              aod_asy_eff(isc) = -(abs(aod_asy(isc)) ** scatter_order)
           else
              aod_asy_eff(isc) = 0. 
           endif
         enddo ! isc

!        cosp_frac = max(1.-scatter_order,0.)

!        Get src term at selected azimuths (every 10 degrees)
         dmintopo = minval(dist_2_topo(ialt,:))
         dmaxtopo = maxval(dist_2_topo(ialt,:))
         if(dmintopo .eq. 0. .and. htmin_view .le. 130000.)then 
           l_dlow = .true.  ! incomplete terrain in ring and ray not above atm
         else
           l_dlow = .false. ! complete terrain in ring or ray outside atm
         endif

         if((sol_alt * altray .lt. 100. .or. sol_alt .lt. 0. .or. altray .lt. 0.) .AND. (l_dlow .eqv. .true.) )then
           l_dlow2 = .true.
         else
           l_dlow2 = .false.
         endif

         if(altray .le. 10.0)then
           jazi_d10 = nint(5./azi_scale)
         else
           jazi_d10 = nint(10./azi_scale)
         endif

         write(6,67)altray,sol_alt,sol_alt*altray,dmintopo,htmin_view,l_dlow,l_dlow2
67       format(' altray/sol_alt/prod/dmintopo/htmin_view/l_dlow,l_dlow2',2f8.2,f8.1,2f12.0,2l2)
!        write(6,*)'az range at this alt2',ialt &
!                  ,minval(view_az(ialt,:)),maxval(view_az(ialt,:))
         do jazi = jazi_start,jazi_end,jazi_d10
          altray = view_alt(ialt,jazi)
          view_altitude_deg = altray
          view_azi_deg = view_az(ialt,jazi)
!         write(6,*)'ialt/jazi/view_az',ialt,jazi,view_az(ialt,jazi)

!         if(jazi .eq. jazi_end)then ! test for now
          if(.true.)then ! test for now
            do ic = 1,nc
              od_g_vert = ext_g(ic) * patm
              od_o_vert = ext_o(ic) * patm_o3(htmsl)
              od_a_vert = aod_vrt * ext_a(ic)

!             Needed when sol_alt and/or altray are low
              if((sol_alt * altray .lt. 100. .or. sol_alt .lt. 0. .or. altray .lt. 0.)  &
                               .AND.(l_dlow .eqv. .true.) )then 
!               if(jazi .eq. jazi_start .AND. altray .eq. nint(altray) &
                if(sol_azi .gt. 180)then
                  azi_ref = 270.                 
                else
                  azi_ref = 90.                 
                endif
!               azi_ref = 175. ! volcanic test
                if(view_azi_deg .eq. azi_ref.AND. (altray .eq. nint(altray) .or. abs(viewalt_limb_true) .le. 1.0) & 
!               if(abs(180. - view_azi_deg) .eq. 1.0 .AND. (altray .eq. nint(altray) .or. abs(viewalt_limb_true) .le. 1.0) & 
!                                       .AND. sol_alt .ge. 0. &
                                        .AND. ic .eq. icd)then
                  idebug_clr = 1 ! * idebug_a(ialt,jazi)
                elseif(sol_alt .ge. twi_0 .AND. (ic .eq. icd .or. altray .eq. 90. .or. altray .eq. 0.))then
!               elseif(ialt .eq. ialt_start .AND. ic .eq. icd)then
                  idebug_clr = 1 * idebug_a(ialt,jazi)
                else
                  idebug_clr = 0
                endif

                if(idebug_clr .eq. 1)then
                  write(6,68)ialt,jazi,ic,sol_alt,altray,view_azi_deg &
                            ,dmintopo,dmaxtopo,ao
68                format(/' call   get_clr_src_dir_low',i4,2i5,3f8.2,2f9.0,f9.3)
                endif
                if(aa_90 .gt. 0.)then
                  aa_o_aa_90 = aa/aa_90
                  aa_s_o_aa_90 = aa_s/aa_90
                else
                  aa_o_aa_90 = 0.
                  aa_s_o_aa_90 = 0.
                endif

!               Determine relevant solar altitude along ray
                if(htmsl .gt. 150e3)then ! highalt strategy
                  if(altray .gt. 0.)then
                    refdist_solalt = 0.
                    solalt_ref = sol_alt
                  else
                    refdist_solalt = sqrt( (earth_radius+htmsl)**2 &
                                          - earth_radius**2       )
                    if(trace_solalt(ialt,jazi) .ne. sol_alt)then ! valid trace
                      solalt_ref = trace_solalt(ialt,jazi)
                    else
                      solalt_ref = sol_alt - (altray*cosd(view_azi_deg-sol_azi))
                      solalt_ref = min(solalt_ref,+180.-solalt_ref)
                      solalt_ref = max(solalt_ref,-180.-solalt_ref)
                    endif
                  endif
                else
                  refdist_solalt = 0.
                  solalt_ref = sol_alt
                endif

                od_o_msl = (o3_du/300.) * ext_o(ic)

                call get_clr_src_dir_low(sol_alt,sol_azi, &
                  altray,view_azi_deg, &
                  ext_g(ic),od_g_vert,od_o_msl,od_o_vert, &
                  od_a_vert,ext_ha(ic),htmsl,ssa, &
                  ag/ag_90,ao,aa_o_aa_90, &
                  aod_ref,aero_refht,aero_scaleht, &
                  ag_s/ag_90,aa_s_o_aa_90, &
                  ags_a,aas_a,isolalt_lo,isolalt_hi,ic,idebug_clr, &
                  refdist_solalt,solalt_ref, &
                  srcdir_a(ic,jazi),sumi_g_a(ic,jazi),sumi_a_a(ic,jazi),&
                  opac_slant,nsteps,dsl,tausum_a)
                if(idebug_clr .eq. 1)then
                  write(6,69)ialt,jazi,ic,sol_alt,trace_solalt(ialt,jazi),altray &
                            ,srcdir(ic),srcdir_a(ic,jazi)
69                format(' called get_clr_src_dir_low',i4,2i5,3f8.2,2f9.4/)
                endif

              else ! use generic values for azimuths every 10 degrees
                srcdir_a(ic,jazi) = srcdir(ic)
                sumi_g_a(ic,jazi) = sumi_gc(ic)
                sumi_a_a(ic,jazi) = sumi_ac(ic)

              endif
            enddo ! ic

          endif ! .true.
         enddo ! jazi

         do jazi = jazi_start,jazi_end,2 ! get eclipse parms at selected azimuths
          altray = view_alt(ialt,jazi)
          view_azi_deg = view_az(ialt,jazi)
          if(sol_alt .gt. 0.)then

!           Get eclipse parms if applicable (and effective brightness)
            if(l_solar_eclipse .eqv. .true.)then
              do ic = 1,nc
!               Compute/Interpolate weighted intensity using distecl
                if( ((jazi-1)/2)*2+1 .eq. jazi)then
!               if(.true.)then
                  eobsc_sum(:) = 0.
                  emag_sum = 0.
                  do iopac = 1,nopac
                    call sun_eclipse_parms(i4time,rlat,rlon,htmsl,0 &
                     ,altray,view_azi_deg,distecl(iopac,ic) &
                     ,earth_radius,elgms,emag,eobscf,eobsc_a(iopac,:))
                    eobsc_sum(:) = eobsc_sum(:) + eobsc_a(iopac,:)
                    emag_sum  = emag_sum  + emag
                  enddo

                  eobsc(:,jazi) = eobsc_sum(:) / float(nopac)
                  emag_a(jazi)  = emag_sum  / float(nopac)
!               else ! fill in from previous azimuth
!                 eobsc(:,jazi) = eobsc(:,jazi-1)
!                 emag_a(jazi) = emag_a(jazi-1)
                endif
              enddo ! ic
            endif
          endif ! sun above horizon
         enddo ! jazi

         do jazi = jazi_start,jazi_end ! main azimuth loop
          altray = view_alt(ialt,jazi)
          view_altitude_deg = altray
          view_azi_deg = view_az(ialt,jazi)

!         Interp srcdir(:), sumi_gc, sumi_ac
          if(jazi .lt. jazi_end)then
            jazi_interpl = ((jazi-jazi_start)/jazi_d10) * jazi_d10 + jazi_start
          else
            jazi_interpl = jazi_end - jazi_d10
          endif
          jazi_interph = jazi_interpl + jazi_d10 
          fazi = (float(jazi) - float(jazi_interpl)) / float(jazi_d10)
          srcdir(:) = (1.-fazi) * srcdir_a(:,jazi_interpl) & ! testing
                    +     fazi  * srcdir_a(:,jazi_interph)
          sumi_gc(:) = (1.-fazi) * sumi_g_a(:,jazi_interpl) & ! testing
                     +     fazi  * sumi_g_a(:,jazi_interph)
          sumi_ac(:) = (1.-fazi) * sumi_a_a(:,jazi_interpl) & ! testing
                     +     fazi  * sumi_a_a(:,jazi_interph)

          xs = cosd(sol_alt) * cosd(sol_azi)
          ys = cosd(sol_alt) * sind(sol_azi)
          zs = sind(sol_alt)

          xo = cosd(altray) * cosd(view_azi_deg)
          yo = cosd(altray) * sind(view_azi_deg)
          zo = sind(altray)

!         elong(ialt,jazi) = angleunitvectors(xs,ys,zs,xo,yo,zo)

          idebug = idebug_a(ialt,jazi) * 2

          if(sol_alt .gt. twi_0)then ! formerly 0.

!           Get eclipse parms if applicable (and effective brightness)
            if(l_solar_eclipse .eqv. .true.)then
              do ic = 1,nc
                if(idebug .ge. 1)then
                  idebuge = 2
                else
                  idebuge = 0
                endif

!               Compute/Interpolate weighted intensity using distecl
                if( ((jazi-1)/2)*2+1 .eq. jazi)then
                  continue
                else ! fill in from previous azimuth
                  eobsc(:,jazi) = eobsc(:,jazi-1)
                  emag_a(jazi) = emag_a(jazi-1)
                endif

!               Secondary scattering for each color
!               arg1 = .00004 + float(nc-1) * .00003
                arg1 = .000004 + float(ic-1) * .000002
                arg2 = 5e-7                         

                if(emag_a(jazi) .gt. 1.0)then
                  ecl_scat(ic) = arg1 - (emag_a(jazi)-1.0) * arg2
                else
                  ecl_scat(ic) = arg1 * eobsc(ic,jazi)
                endif

!               Totality darkness fraction
!               Hardwired temporal and spatial variability only used here
                ecl_intd(ic) = (1.0 - eobsc(ic,jazi)) + ecl_scat(ic) * eobsc(ic,jazi)

                ecl_dir_rat(ic) = 1.0 - (ecl_scat(ic) / ecl_intd(ic))

                if(idebuge .ge. 1)then
                  write(6,70)altray,distecl(1,ic),distecl(nopac,ic),elgms &
                      ,emag_a(jazi),eobsc_a(1,ic),eobsc_a(nopac,ic),eobsc(ic,jazi) &
                      ,ecl_scat(ic),ecl_intd(ic),ecl_dir_rat(ic)
70                format('eclipse altray/dist/elg/emag',f9.2,2f9.1,f9.4,f7.4, &
                         ' eobsc',3f7.4,' scat/int/dir',2f10.7,f7.4)
                endif
              enddo ! ic
              if(idebug .ge. 1)then
                  write(6,*)'ecl_scat = ',ecl_scat
              endif
            else ! no eclipse
              ecl_intd(:) = 1.0          
              ecl_dir_rat(:) = 1.0
            endif

!           Fraction of atmosphere illuminated by the sun 
!           clear_radf_c(3,ialt,jazi) = 1. 
            angle_plane = 0.; dist_ray_plane = 0.; ht_ray_plane = 0.

!           Illumination using Rayleigh and HG phase functions
            z = (90. - altray) + refractd_app(altray      ,patm)
            airmass_g = ag
!           sb_corr = 2.0 * (1.0 - (sind(sol_alt)**0.5)) ! magnitudes
            sb_corr = 0.0
!           day_int = 3e9 / 10.**(0.4 * sb_corr) * ecl_intd

!           HG illumination
            if(.false.)then
!             do ic = 1,nc
              hg2 = aod_bin(1) * hg(aod_asy_eff(1),elong(ialt,jazi)) &
                  + aod_bin(2) * hg(aod_asy_eff(2),elong(ialt,jazi)) &
                  + aod_bin(3) * hg(aod_asy_eff(3),elong(ialt,jazi)) 
!             enddo ! ic

              hg2t = hg2

!             non-topo phase function with variable scatter order
              fc = fcterm2 * 0.5**(scatter_order-1.0)
              gc = 2300. / scatter_order**2
              hg2 = (1.-fc) * hg2 + fc * cosp(gc,elong(ialt,jazi))
!             hg2(:) = (1.-fc) * hg2(:) + fc * cosp(gc,elong(ialt,jazi))

!             topo phase function assumes scatter order is 1
              fc = fcterm2
              hg2t = (1.-fc) * hg2t + fc * cosp(2300.,elong(ialt,jazi))

            else
              fb = 0.55**scatter_order
!             g1 = 0.58**scatter_order
              g1 = 0.45**scatter_order
              g2 = 0.962**scatter_order
              hg2 = dhg2(elong(ialt,jazi),fb,fcterm2)

!             topo phase function assumes scatter order is non-topo 
!             value (for now)
              hg2t = hg2

            endif

            if(aod_ray(ialt,jazi) .gt. 0.)then
                aod_dir_rat = aod_ray_dir(ialt,jazi) / aod_ray(ialt,jazi)
            else
                aod_dir_rat = 0.
            endif

!           Consider arg for sideways scattering from downward diffuse
!           is ~hg(90.) or an evaluated integral. Use 0.5 for now.
            if(l_solar_eclipse .eqv. .false.)then
              hg2d = (hg2  * aod_dir_rat) + (0.5 * (1.0 - aod_dir_rat)) 
              hg2t = (hg2t * aod_dir_rat) + (0.5 * (1.0 - aod_dir_rat)) 
            else
              do ic = 1,nc
                arg = aod_dir_rat * ecl_dir_rat(ic)
                hg2d(ic) = (hg2 * arg) + (0.5 * (1.0 - arg)) 
              enddo ! ic
            endif

            if(aod_dir_rat .gt. 1.0001 .OR. aod_dir_rat .lt. 0.0)then
              write(6,*) &
            ' ERROR in skyglow_phys: aod_dir_rat out of bounds',aod_dir_rat
              write(6,*)' ialt/jazi/altray = ',ialt,jazi,altray
              write(6,*)' aod_ray_dir/aod_ray = ',aod_ray_dir(ialt,jazi),aod_ray(ialt,jazi)
              stop
            endif

            rayleigh_pfunc = rayleigh_pf(elong(ialt,jazi))
            rayleigh_gnd = rayleigh_pfunc + sfc_alb * sind(sol_alt)

            if(htmsl .le. 250000. .and. dist_2_topo(ialt,jazi) .eq. 0.)then 
!           if(htmsl .le. 100000e3 .and. dist_2_topo(ialt,jazi) .eq. 0.)then 
              mode_sky = 4 ! simple approach (<150km non-topo case)
              do ic = 1,nc
                day_int = 3e9 * (1.0 - eobsc(ic,jazi))

                if(airmass_g .gt. 0.)then
                  gasf = ag_2d(ialt,jazi) / airmass_g
                else
                  gasf = 1.0
                endif

!               od_a = aod_ref * aa * ext_a(ic)  ! slant od
                od_a_slant = od_a_vert * aa_o_aa_90 
                od_a = od_slant
                if(od_a_slant .gt. 0. .and. sol_alt .gt. twi_alt)then ! normalize with extinction?
!                 aodf = min(aod_ill(ialt,jazi)/od_a_slant,1.0)
                  aodfo = min(aod_ill(ialt,jazi)/aod_tot(ialt,jazi),1.0)
                  aodf = aod_ill_opac(ialt,jazi)/aod_ill_opac_potl(ialt,jazi)
                else
                  aodfo = 1.0
                  aodf = 1.0
                endif

!               We can still add in gasf and aodf
                clear_rad_c(ic,ialt,jazi) = day_int * &
!                   (rayleigh_pfunc * sumi_gc(ic) * gasf + &
                    (rayleigh_gnd   * sumi_gc(ic) * gasf + &
                     hg2d(ic)       * sumi_ac(ic) * aodf   )

                if(idebug .ge. 1 .AND. &
                                    (ic .eq. icd .or. altray .eq. 90.))then
                  write(6,71)ic,day_int,elong(ialt,jazi) &
                      ,sumi_gc(ic),sumi_ac(ic),gasf,aodf,aodfo &
                      ,rayleigh_pfunc,hg2,aod_dir_rat,hg2d(ic) &
                      ,clear_rad_c(ic,ialt,jazi)
71                format('day_int/elong/sumi_g/sumi_a/gasf/aodf/aodfo', &
                         '/rayleigh/hg2/aodr/hg2d/clrrd4', &
                         i2,f12.0,f8.3,2x,2f11.8,2x,3f8.3,2x,4f8.3,f12.0)      
                endif
              enddo ! ic

!           The first option has a discontinuity around -15 degrees above 20km
!           The second option has a brightness jump around 18000m
!           elseif(htmsl .gt. 18000. .and. dist_2_topo(ialt,jazi) .eq. 0.)then ! high vantage point
            elseif(dist_2_topo(ialt,jazi) .eq. 0.)then ! non-topo (high vantage point)
!           elseif(htmsl .gt. 18000.)then ! high vantage point

!             Consider case where 'dist_2_topo' is zero (in a cloud) and
!             'topo_solalt' wasn't set. In this case 'day_int' wouldn't 
!             be properly determined. Here we use info from a new array 
!             called 'trace_solalt'.
              solalt_atmos = trace_solalt(ialt,jazi) 
              if(solalt_atmos .gt. 0.)then
                sol_occ = 1.0 ! visible
              else
                sol_occ = 0.0 ! invisible
              endif

!             Use/restore simpler approach for downward rays from high alt
              do ic = 1,nc
!               Rayleigh illumination
                day_int = 3e9 * sol_occ * (1.0 - eobsc(ic,jazi))
                od_per_am = ext_g(ic)
                rayleigh = brtf(airmass_g,od_per_am) * rayleigh_gnd
                od_a = aod_ref * aa * ext_a(ic)  ! slant od
                brta = opac(od_a) * hg2d(ic)

!               Total illumination
                if(htmin_view .le. 100e3)then
!                 clear_rad_c(ic,ialt,jazi) = day_int * (rayleigh + brta)
                  clear_rad_c(ic,ialt,jazi) = day_int * &
                      (rayleigh       * sumi_gc(ic) + &
                       hg2d(ic)       * sumi_ac(ic)   )
                else
                  clear_rad_c(ic,ialt,jazi) = 0.
                endif

                mode_sky = 3 ! high vantage point >250km (non-topo)

                if(idebug .ge. 1 .and. ic .eq. icd)then
                  write(6,72)day_int,elong(ialt,jazi),airmass_g,brtf(airmass_g,od_per_am) &
                            ,rayleigh,hg2,hg2d(ic),sumi_gc(ic),sumi_ac(ic),clear_rad_c(2,ialt,jazi)
72                format('day_int/elong/ag/brtf/rayleigh', &
                         '/hg2/hg2d/sumi/clrrd3',f12.0,6f8.3,2f10.5,f12.0)      
                endif

               enddo ! ic

               if(idebug .ge. 1)then
                 write(6,73)day_int,elong(ialt,jazi),rayleigh,hg2d(2),clear_rad_c(:,ialt,jazi)
73               format('day_int/elong/hg2/clear_rad :',f12.0,3f8.3,3f12.0)      
               endif

            else ! assume two scattering layers (g+a [lower-1], g [upper-2])
              do ic = 1,nc
                gasfrac = ag_2d(ialt,jazi) / airmass_g ! topo illuminated fraction of gas
                if(dist_2_topo(ialt,jazi) .eq. 0.)then    ! free of topo
                  od_g = ext_g(ic)*airmass_g              ! slant od
!                 od_a = aod_ray(ialt,jazi)*aa*ext_a(ic)  ! slant od
                  od_a = aod_ref           *aa*ext_a(ic)  ! slant od
!                 if(altray .gt. 0.)then
                    nlyr = 2
!                 else ! simplify things
!                   nlyr = 1
!                 endif
                else                                      ! hitting topo
                  ht_to_topo = dist_2_topo(ialt,jazi) * sind(altray)
                  od_g_vert = ext_g(ic) * patm
!                 od_o_vert = ext_o(ic) * patm_o3(htmsl)
                  od_a_vert = aod_vrt * ext_a(ic)
                  alpha_g = od_g_vert / 8000.
                  alpha_a = od_a_vert / aero_scaleht
                  od_g = dist_2_topo(ialt,jazi) * alpha_g ! slant od
                  od_a = dist_2_topo(ialt,jazi) * alpha_a ! slant od
                  if(ht_to_topo .gt. aero_scaleht)then
                    nlyr = 2
                  else
                    nlyr = 1
                  endif
                  if(l_topo_a(jazi) .eqv. .true.)then
                    if( idebug_a(ialt,jazi) .gt. 0 .AND. &
                       (ic .eq. icd .or. sol_alt .lt. 4.)     )then
                      idebug_topo = 1
                    elseif(altray .le. -60. .and. altray .eq. nint(altray)&
                     .and. htmsl .gt. 1000000. .and. ic .eq. icd &
                     .and. jazi .eq. minazi)then
                      idebug_topo = 1
                    else
                      idebug_topo = 0
                    endif
                    if(aa_90 .gt. 0.)then
                      aa_o_aa_90 = aa/aa_90
                    else
                      aa_o_aa_90 = 0.
                    endif

                    if(htmsl .gt. 150e3)then ! highalt strategy
                      refdist_solalt = dist_2_topo(ialt,jazi)
                      solalt_ref = topo_solalt(ialt,jazi)
                    else
                      refdist_solalt = 0.
                      solalt_ref = sol_alt
                    endif

                    od_o_msl = (o3_du/300.) * ext_o(ic)

                    call get_clr_src_dir_topo(sol_alt,sol_azi, &        ! I
                     altray,view_azi_deg, &
                     ext_g(ic),od_g_vert,od_o_msl,od_a_vert, &          ! I
                     htmsl,dist_2_topo(ialt,jazi), &
                     ssa,ag/ag_90,aa_o_aa_90, &
                     aod_ref,aero_refht,aero_scaleht, &                 ! I
                     ags_a,aas_a, &                                     ! I
                     isolalt_lo,isolalt_hi,ic,idebug_topo, &
                     nsteps,refdist_solalt,solalt_ref, &
                     sumi_gct(ic),sumi_act(ic),opac_slant,dst,tausum_a)  

                    if(idebug_topo .gt. 0)then
                      write(6,81)ialt,jazi,ic,sol_alt,altray &
                         ,dist_2_topo(ialt,jazi),sumi_gct(ic),sumi_act(ic)
81                    format(' called get_clr_src_dir_topo',i4,i5,i4,2f8.2,f10.0,2f9.4)
                    endif ! write log info
                  endif ! call get_clr_src_dir_topo

!                 Normalize by extinction?
                  aodfo = aod_ill(ialt,jazi) / aod_2_topo(ialt,jazi) 
                  aodf = aod_ill_opac(ialt,jazi)/aod_ill_opac_potl(ialt,jazi)
                  day_int = 3e9 * (1.0 - eobsc(ic,jazi))
                  clear_rad_topo = day_int * &
                    (sumi_gct(ic) * rayleigh_gnd * clear_radf_c(ic,ialt,jazi) + &
!                    sumi_act(ic) * hg2t         * aodf * ssa) ! take out ssa?
                     sumi_act(ic) * hg2t         * aodf) 
                  if(idebug .ge. 1 .AND. ic .eq. icd)then
                    write(6,*)' computed clear_rad_topo = ',clear_rad_topo,sumi_gct(ic),sumi_act(ic)
                  endif
                  if(l_solar_eclipse .eqv. .false.)then ! experimental
                    clear_rad_c(ic,ialt,jazi) = clear_rad_topo
                  endif
                  elglim = 5.0 ! catch aureole on topo
                  mode_sky = 2 ! hitting topo case

                  if(idebug .ge. 1 .AND. ic .eq. icd)then
                    write(6,91)elong(ialt,jazi),sumi_gct(ic),rayleigh_gnd &
                       ,clear_radf_c(ic,ialt,jazi),sumi_act(ic),hg2t &
                       ,aodf,aodfo,ssa,clear_rad_topo      
91                  format('elg/ig/rayg/radf/ia/hg2t/aodf/aodfo/ssa/clear_rad_topo =',f7.2,3x,3f9.4,3x,4f9.4,f7.2,3x,f12.0)
                  endif

                endif ! hitting topo

                if(dist_2_topo(ialt,jazi) .eq. 0.)then    ! free of topo
                  alphav_g = od_g / 8000.           ! extinction per vert m
                  alphav_a = od_a / aero_scaleht    ! extinction per vert m
                  if(od_a .gt. 0)then  ! lower layer
                    od_g1 = od_a * (alphav_g / alphav_a)
                  else
                    od_g1 = 0.       
                  endif
                  od_g2 = od_g - od_g1 ! upper layer

                  am_sun = airmassf(90. - sol_alt,patm) * (od_g2/od_g)
                  solar_int_g2 = trans(am_sun*ext_g(ic))

!                 Effective Rayleigh Phase Factor considering shadowing
                  rayleigh_pf_eff = rayleigh_gnd * clear_radf_c(ic,ialt,jazi)
                  pf_eff1 = (rayleigh_pf_eff * alphav_g + hg2d(ic) * alphav_a * solar_int_g2) &
                          / (alphav_g + alphav_a * solar_int_g2)
!                 od_1 = od_g1*gasfrac + od_a
                  od_1 = od_g1*gasfrac + aod_ill(ialt,jazi) ! lower layer ill od
                  brt1 = brto(od_1) * pf_eff1 ! lower layer
!                 brt2 = brto(od_g2) * rayleigh_pf_eff
                  brt2 = brto(od_g2*gasfrac) * rayleigh_pf_eff ! upper layer
                  if(brt2 .eq. 0)brt2 = (od_g2*gasfrac) * rayleigh_pf_eff
                  trans1 = trans(od_1)

                  day_int = 3e9 * (1.0 - eobsc(ic,jazi))
                  if(nlyr .eq. 2)then ! free of topo, or high topo case
                    if(l_solar_eclipse .eqv. .false.)then
                      clear_rad_c(ic,ialt,jazi) = day_int * ((1.-trans1)*brt1 + trans1*brt2) * (srcdir(ic)/srcdir_90(ic))
                    else
                      clear_rad_c0(ic,ialt,jazi) = 3e9    * ((1.-trans1)*brt1 + trans1*brt2) * (srcdir(ic)/srcdir_90(ic))
                      clear_rad_c(ic,ialt,jazi) = day_int * ((1.-trans1)*brt1 + trans1*brt2) * (srcdir(ic)/srcdir_90(ic))
                    endif
                  else                ! single layer, hitting (low) topo
                    if(l_solar_eclipse .eqv. .false.)then
                      clear_rad_c(ic,ialt,jazi) = day_int *              brt1                * (srcdir(ic)/srcdir_90(ic))
                    else
                      clear_rad_c0(ic,ialt,jazi) = 3e9    *              brt1                * (srcdir(ic)/srcdir_90(ic))
                      clear_rad_c(ic,ialt,jazi) = day_int *              brt1                * (srcdir(ic)/srcdir_90(ic))
                    endif
                  endif

                  clear_rad_2lyr = clear_rad_c(ic,ialt,jazi)
                  elglim = 0.5
                  mode_sky = 1 ! free of topo (discontinued)
                  write(6,*)' Stop to check if mode_sky=1 is still used'
                  stop
                endif ! free of topo

                if(idebug .ge. 1 .AND. (ic .eq. icd .or. &
                        abs(elong(ialt,jazi)-90.) .le. 0.5 .or. &
                            elong(ialt,jazi)      .le. elglim)        )then
                  jazim1 = max(jazi-1,minazi)
                  write(6,82)day_int,ic,jazi,eobsc(ic,jazim1:jazi),srcdir(ic),srcdir_90(ic)
82                format('day_int/eobsc/srcdir/90',f12.0,2i5,4f8.4)           
                  write(6,83)day_int,elong(ialt,jazi),airmass_g,od_g &
                            ,aod_ray(ialt,jazi),aa &
                            ,od_a,alphav_g*1e3,alphav_a*1e3,od_g1,od_g2 &
                            ,clear_rad_2lyr            
83                format(&
     'day_int/elg/ag/od_g/aod_ray/aa/od_a/alphav_g/alphav_a/od_g1/od_g2/clr_rad :' &
                        ,f12.0,f7.2,2f7.3,f8.3,2f7.3,1x,2f7.3,1x,2f8.4,f12.0)      
                  write(6,84)altray,view_azi_deg,am_sun,solar_int_g2,aascat,scatter_order,hg2,hg2d(ic),gasfrac,aod_ill(ialt,jazi) & 
                            ,aod_ray_dir(ialt,jazi),aod_ray(ialt,jazi),aod_dir_rat
84                format('altaz/amsun/solarintg2/aasc/sco/hg2/hg2d/gasfrac/aodill/dir/ray/rat = ',2f8.2,6f9.4,2f10.6,3f7.3)                  
                  if(mode_sky .eq. 1)then ! can be 1 or 2 here
                    write(6,85)dist_2_topo(ialt,jazi),nlyr,od_1,clear_radf_c(ic,ialt,jazi),rayleigh_pf_eff,pf_eff1,brt1,brt2,trans1
85                  format('dst/nlyr/od_1/radf/pf_eff/pf_eff1/brt1/brt2/trans1',f8.0,i2,8f9.4)
                    write(6,*)'rad recalc',day_int,trans1,brt1,brt2,srcdir(ic),srcdir_90(ic), &
                      day_int * ((1.-trans1)*brt1 + trans1*brt2) * (srcdir(ic)/srcdir_90(ic))
                  endif
                  write(6,86)ag_2d(ialt,jazi),airmass_g,gasfrac
86                format('ag2d/airmass_g/gasfrac',f10.4,f9.3,f10.4)
!                 if(dist_2_topo(ialt,jazi) .gt. 0.)then    ! hit topo
!                 endif
                endif ! idebug

              enddo ! ic

            endif ! mode determination

            if(idebug .ge. 1)then
              rmaglim = b_to_maglim(clear_rad_c(2,ialt,jazi))
            endif

          else ! sun below horizon (twilight / night brightness)
            continue
            mode_sky = 0
          endif ! sun above horizon

          if(idebug .ge. 1)then 
              write(6,111)altray,view_azi_deg,sol_alt,topo_solalt(ialt,jazi) &
                       ,trace_solalt(ialt,jazi),mode_sky,od_a &
                       ,dist_2_topo(ialt,jazi),aodf,aodfo,clear_rad_c(:,ialt,jazi)
111           format('alt/azi/salt-t-t/mode/od_a/dst/aodf/aodfo/clrrd' &
                    ,2f6.1,1x,3f7.2,i3,f12.5,f9.0,2f9.4,3e14.5)

              if(clear_rad_c(1,ialt,jazi) .lt. 0. .or. clear_rad_c(1,ialt,jazi) .gt. 1e20)then
                write(6,*)' ERROR clrrd(1) out of bounds',clear_rad_c(1,ialt,jazi)
!               stop
              endif

              if(idebug .ge. 2)then
                if(altray .ne. 0 .or. view_azi_deg .le. 180.)then
                  write(6,121)rmaglim
121               format('rmaglim = ',f8.3)
                else
                  write(6,122)rmaglim
122               format('rmaglim = ',f8.3,'    horizon')
                endif
              endif
              write(6,*)
              write(6,*)' clearrad1 column:',clear_rad_c(1,minalt:minalt+9,minazi)
          endif ! idebug .eq. 1

         enddo ! jazi
        enddo ! ialt

        I4_elapsed = ishow_timer()

        write(6,*)' clearrad1 column:',clear_rad_c(1,minalt:minalt+9,minazi)

!       Apply during daylight eclipses or during twilight
        if( (sol_alt .ge. 0. .and. (l_solar_eclipse .eqv. .true.)) .OR. &
            (sol_alt .lt. 10. .and. sol_alt .gt. twi_0)               )then

            write(6,*)'Calling get_sky_rad_ave'
            do ic = 1,nc ! secondary scattering in each color
!               Good temporal variability, constant for spatial
                call get_sky_rad_ave(clear_rad_c(ic,:,:) &
                    ,view_alt,view_az,maxalt-minalt+1,maxazi-minazi+1 &
                    ,sol_alt,sol_azi,sky_rad_ave(ic))
                write(6,*)'range of clear_rad_c',minval(clear_rad_c(ic,:,:)),maxval(clear_rad_c(ic,:,:))
                od_g_vert = ext_g(ic) * patm
                od_o_vert = ext_o(ic) * patm_o3(htmsl)
                od_a_vert = aod_vrt * ext_a(ic)
                od_vert = od_g_vert + od_a_vert

!               Assuming non-reflective land for now
                znave = clear_rad_c(ic,maxalt,minazi) / sky_rad_ave(ic)
                if(sol_alt .gt. 0.)then
                  highalt_adjust = max(min(-log10(znave/5.0),3.0),1.0)
                  viewalt_adjust = 1.0
                else
                  highalt_adjust = 1.0 ! + (-sol_alt/6.0)
                  viewalt_adjust = 1.0 ! view_alt(maxalt,minazi)/90. 
                endif
                write(6,191)ialt_end,maxalt,znave,highalt_adjust,viewalt_adjust
191             format(' ialt_end,maxalt,znave,hadj,vadj',2i7,f9.5,2f9.3)

                jazi_dbg = min(540,maxazi)
                if(ic .eq. icd)then
                  write(6,*)' single scat clear rad at ic/azi',ic,view_az(ialt_start,jazi_dbg)
                endif

                clear_rad_ref = clear_rad_c0(ic,maxalt,minazi) ! eclipse
                clear_rad_ref = clear_rad_c(ic,maxalt,minazi)  ! no eclipse

                do ialt = ialt_start,ialt_end
                  altray = view_alt(ialt,jazi_start)
                  arg = sind(min(max(altray,1.5),90.))
                  od_slant = od_vert / arg                      
                  fracg = od_g_vert / od_vert
                  fraca = od_a_vert / od_vert
                  if(ic .eq. icd)then
                    write(6,*)'alt,od_slant,clear_rad_c',altray,od_slant,clear_rad_c(ic,ialt,jazi_dbg)
                  endif

                  sph_rad_ave = sky_rad_ave(ic) * 0.5 ! sphere ave, drk gnd

                  do jazi = jazi_start,jazi_end
                    if(dist_2_topo(ialt,jazi) .eq. 0.)then ! unobstructed by terrain
                      sky_rad_scatg = sph_rad_ave * viewalt_adjust & 
                                    * fracg * opac(od_slant) * highalt_adjust
                      sky_rad_scata = (0.4 * sph_rad_ave * viewalt_adjust + 0.6 * clear_rad_c(ic,ialt,jazi)) & 
                                    * fraca * opac(od_slant) * highalt_adjust
                    else ! topo in path
                      od_slant_a = aod_2_topo(ialt,jazi)
                      sky_rad_scatg = 0. ! neglected term for now
                      sky_rad_scata = sph_rad_ave * opac(od_slant_a) * highalt_adjust * viewalt_adjust
                    endif
                    sky_rad_scat(ic,ialt,jazi) = sky_rad_scatg + sky_rad_scata
                  enddo ! jazi

                  clear_rad_c(ic,ialt,:) = clear_rad_c(ic,ialt,:) + sky_rad_scat(ic,ialt,:)
                  if(ialt .eq. ialt_end)then ! zenith/high point info
                    scatfrac = sky_rad_scat(ic,maxalt,jazi_dbg) / clear_rad_c(ic,maxalt,jazi_dbg)
                    ri_over_i0 = clear_rad_c(ic,maxalt,jazi_dbg) / clear_rad_c0(ic,maxalt,jazi_dbg)
                    ri_over_f0 = clear_rad_c(ic,maxalt,jazi_dbg) / 3e9                            
                    rmaglim = b_to_maglim(clear_rad_c(2,ialt_end,jazi_dbg))
                    write(6,201)ic,altray,view_az(ialt_end,jazi_dbg),sph_rad_ave,od_vert,clear_rad_ref,sky_rad_scat(ic,maxalt,jazi_dbg),clear_rad_c(ic,maxalt,jazi_dbg),scatfrac,ri_over_i0
201                 format(' ic/alt/az/sph_rad_ave/od/c0/scat/rad/rat/scatf/ii0',i2,2f9.2,f11.0,f6.3,f11.0,2f11.0,f7.4,2f10.7)
                    if(ic .eq. icd)then
                      b = clear_rad_c(2,ialt_end,jazi_dbg)
                      write(6,211)b,b_to_v(b),rmaglim
                    endif
211                 format(' b/mag/rmaglim at high point (second scat) =' &
                          ,f10.1,2f9.2)
                  endif ! at zenith 
                enddo ! ialt
            enddo ! ic
        endif

        write(6,*)' clearrad2:',clear_rad_c(2,(minalt+maxalt)/2,minazi:maxazi:10)
        write(6,*)' returning from skyglow_phys'

        return
        end
        
