

        subroutine skyglow_phys(ialt_start,ialt_end,ialt_delt      &! I
                   ,jazi_start,jazi_end,jazi_delt,azi_scale        &! I
                   ,minalt,maxalt,minazi,maxazi,idebug_a           &! I
                   ,sol_alt,sol_azi,view_alt,view_az,twi_0,twi_alt &! I
                   ,sol_lat,sol_lon,solalt_limb_true,del_solalt    &! I
                   ,isolalt_lo,isolalt_hi,topo_solalt,trace_solalt &! I
                   ,earth_radius,patm,sfc_alb_c                    &! I
                   ,aod_vrt,aod_ray,aod_ray_dir                    &! I
                   ,aod_ref,aero_scaleht,dist_2_topo               &! I
                   ,htmsl,redp_lvl,horz_dep,eobsl                  &! I
                   ,aod_ill,aod_2_topo,aod_tot,ext_g               &! I
                   ,l_solar_eclipse,i4time,rlat,rlon               &! I
                   ,clear_radf_c,ag_2d                             &! I
                   ,od_g_slant_a,od_o_slant_a,od_a_slant_a,ext_a   &! O
                   ,sky_rad_scat                                   &! O
                   ,clear_rad_c,sky_rad_ave,elong           )       ! O/I

!       Sky glow with solar altitude > 0

        include 'trigd.inc'

        use mem_namelist, ONLY: fcterm, aod_bin, aod_asy, ssa, r_missing_data
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
        angdist(p1,p2,dlon) = acosd(sind(p1) * sind(p2) + cosd(p1) * cosd(p2) * cosd(dlon))
        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1
        curvat(hdst,radius) = hdst**2 / (2. * radius)
        expbh(x) = exp(min(x,+80.))

        real linecyl,linecylp
        linecylp(xcos,ycos,x1,r) = (-(2.*xcos*x1) + sqrt((2.*xcos*x1)**2 - 4. * (xcos**2 + ycos**2) * (x1**2 - r**2))) &
                                 /  (2. * (xcos**2 + ycos**2))

        include 'rad_nodata.inc'

        cosp(b,x) = (1. + b/2.) * cosd(x/2)**b

        real twi_trans_c(nc)           ! transmissivity

        real hg2d(nc), alphav_g, alphav_a, hg2(nc), hg2t(nc), ssa_eff(nc)
        real srcdir_90(nc),srcdir(nc),clear_int_c(nc)
        real sumi_gc(nc),sumi_ac(nc)
        real sumi_gct(nc),sumi_act(nc)
        real srcdir_a(nc,minazi:maxazi)
        real sumi_g_a(nc,minazi:maxazi),sumi_a_a(nc,minazi:maxazi)
        real ecl_intd(nc),ecl_dir_rat(nc),ecl_scat(nc),sky_rad_ave(nc)

        real clear_rad_c(nc,minalt:maxalt,minazi:maxazi)  ! clear sky illumination
!       real clear_rad_c0(nc,minalt:maxalt,minazi:maxazi) ! clear sky (no eclipse)
        real clear_radf_c(nc,minalt:maxalt,minazi:maxazi)! integrated
               ! fraction of gas illuminated by the sun along line of sight
               ! (consider Earth's shadow + clouds, used when sun is below
               !  the horizon), attenuated from consideration behind clouds
        real ag_2d(minalt:maxalt,minazi:maxazi) ! gas airmass (topo/notopo)
        real aod_ill(minalt:maxalt,minazi:maxazi) ! aerosol illuminated
                                      ! optical depth (slant - topo/notopo)
        real aod_tot(minalt:maxalt,minazi:maxazi) ! aerosol 
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
        real*8 dist_2_topo(minalt:maxalt,minazi:maxazi)
        real topo_solalt(minalt:maxalt,minazi:maxazi)
        real trace_solalt(minalt:maxalt,minazi:maxazi)
        real elong(minalt:maxalt,minazi:maxazi)
        integer idebug_a(minalt:maxalt,minazi:maxazi)

        real view_alt(minalt:maxalt,minazi:maxazi)
        real view_az(minalt:maxalt,minazi:maxazi)

!       If below the limb this probably is set using aerosol reference
!       height under the observer
        real od_g_slant_a(nc,minalt:maxalt)
        real od_o_slant_a(nc,minalt:maxalt)
        real od_a_slant_a(nc,minalt:maxalt)

        logical l_solar_eclipse, l_dlow, l_dlow2, l_topo_a(minazi:maxazi)

        parameter (nsteps = 300000)
!       This larger number could enable high altitude eclipse handling if
!       it doesn't have side effects of slowing down regular runs. We 
!       are declaring just the needed section of the array, based on
!       istart passed into 'get_clr_src_dir'.
        parameter (nsteps_low = 100000)
        parameter (nsteps_topo = 1)
        parameter (nopac = 10)

        real distecl(nopac,nc)
        real tausum_a(nsteps)
        real taumid(nopac),opacmid(nopac),eobsc_a(nopac,nc)
        real eobsc(nc,minazi:maxazi),eobsc_sum(nc),emag_a(minazi:maxazi)
        real sky_rad_scat(nc,minalt:maxalt,minazi:maxazi)            
!       real sky_rad_scata(minazi:maxazi)
        real sfc_alb_c(nc)

!       Airmasses relative to zenith for observer ht
        real ags_a(isolalt_lo:isolalt_hi)
        real aas_a(isolalt_lo:isolalt_hi)
        real aos_a(isolalt_lo:isolalt_hi)

        parameter (n_order = 30)
        real rel_order(1:n_order)  ! P of received photon with this order
        real p_order(0:n_order)    ! P of emitted photon with this order

        icd = 3

        eobsc(:,:) = 0. ! initialize
        sky_rad_ave = r_missing_data
!       sfc_alb_c(:) = 0.15  ! pass this in and account for snow cover?

!       aod_bin(1) = .000
!       aod_bin(2) = .987
!       aod_bin(3) = .013

!       aod_asy(1) =  .75
!       aod_asy(2) =  .70
!       aod_asy(3) = -.65

!       Now read via namelist
!       fcterm = 0.09 ! range from 0.00 to 0.09 (large aerosol population)
!                     ! peak phase function of 20 to 110
        thr_abv_clds = 25e3

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
            write(6,*)' ic/wa/ext_a/aod ',ic,wa(ic),ext_a(ic),aod_vrt*ext_a(ic)
        enddo ! ic

!       Obtain reference values of source term
        call get_airmass(90.,htmsl,patm &             ! I
                         ,aero_refht,aero_scaleht &   ! I
                         ,earth_radius,1 &            ! I
                         ,ag_90,ao_90,aa_90,refr_deg) ! O

!       Obtain reference values for sun             
        call get_airmass(sol_alt,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,1 &          ! I
                         ,ag_s,ao_s,aa_s,refr_deg)  ! O

!       Set this for range of values of solar altitude
        if(.false.)then
          ags_a = ag_s/ag_90
          aas_a = aa_s/aa_90
        else
          patm_refht = ztopsa(aero_refht)/1013.25
          do isolalt = isolalt_lo,isolalt_hi
            sol_alt_a = sol_alt + float(isolalt) * del_solalt
            if(sol_alt_a .gt. 0.)then
              call get_airmass(sol_alt_a,0.,1. &                     ! I
                              ,aero_refht,aero_scaleht &             ! I
                              ,earth_radius,0 &                      ! I
                              ,ags_a(isolalt),aos_a(isolalt),aa_dum &! O
                              ,refr_deg)                             ! O

              call get_airmass(sol_alt_a,aero_refht,patm_refht &     ! I
                              ,aero_refht,aero_scaleht &             ! I
                              ,earth_radius,0 &                      ! I
                              ,ag_dum,ao_dum,aas_a(isolalt),refr_deg)! O
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
                              ,earth_radius,0 &                      ! I
                              ,ag_sz,ao_sz,aa_sz,refr_deg)           ! O
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

        if(solalt_limb_true .gt. twi_0)then ! formerly 0.
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
            istart = 1; iend = nsteps; ds = 25.
            call get_clr_src_dir(sol_alt,90.,ext_g(ic), &
                od_g_vert,od_a_vert,ext_ha(ic), &
                htmsl,ssa90,ag_90/ag_90,ao_90,aa_90_o_aa_90, &
                aod_ref,aero_refht,aero_scaleht, &
                ag_s/ag_90,aa_s_o_aa_90,ic,idebug, &
                istart,iend, &
                srcdir_90(ic),sumi_gc(ic),sumi_ac(ic),opac_slant, &
                nsteps,ds,tausum_a)
            write(6,*)' Returning with srcdir_90 of ',srcdir_90(ic)
          enddo ! ic
        endif

        do ialt = ialt_start,ialt_end
  
         altray = view_alt(ialt,jazi_start)

!        if(altray .eq. -73.)write(6,*)'debug check 2',altray,idebug_a(ialt,:)
         
         viewalt_limb_true = altray - (-horz_dep)
         if(altray .ge. 0.)then
           htmin_view = htmsl
         else
           htmin_view = htminf(htmsl,altray,earth_radius)
         endif

         call get_airmass(altray,htmsl,patm &       ! I
                         ,aero_refht,aero_scaleht & ! I
                         ,earth_radius,0 &          ! I
                         ,ag,ao,aa,refr_deg)        ! O
         write(6,63)altray,ag,ao,aa
63       format(' returned from get_airmass with alt/ag/ao/aa = ',f9.4,3e12.4)

         htsfc = 0.
         call get_ray_info(altray,htmsl,htsfc,earth_radius,emis_ang &
                          ,r_missing_data,dist_to_sfc)

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
           if(htmsl .gt. 100e3)then
             jazi_delt_topo = min(jazi_delt_topo,4)
           else
             jazi_delt_topo = min(jazi_delt_topo,8)
           endif
         endif
         l_topo_a = .false.
         l_topo_a(minazi:maxazi:jazi_delt_topo) = .true.
         write(6,*)'topo az delt at alt:',altray,jazi_delt_topo

!        Determine src term and ratio with zenith value
         if((sol_alt .gt. 0. .and. htmsl .lt. 50000e3) &
                             .OR. l_solar_eclipse .eqv. .true.)then

!          Allow greater range of altitudes
           ds = 25. 
           if(htmsl .gt. 100000.)then
               call get_ray_info(altray,htmsl,100e3,earth_radius,altdum,r_missing_data,dist_toa)
               call get_ray_info(altray,htmsl,  0e3,earth_radius,altdum,r_missing_data,dist_sfc)
               if(dist_sfc .ne. r_missing_data)then
                   istart = nint(dist_toa / ds)                    
                   iend = nint(dist_sfc / ds)      
               else
                   istart = 1
                   iend = nsteps
               endif
           else
               istart = 1
               iend = nsteps
           endif
           ioffset = istart-1

           do ic = 1,nc
             od_g_vert = ext_g(ic) * patm
             od_o_vert = ext_o(ic) * patm_o3(htmsl)
             od_a_vert = aod_vrt * ext_a(ic)

!            if((l_solar_eclipse .eqv. .true.) .AND. ic .eq. icd)then
             if((ic .eq. icd .and. altray/5. .eq. nint(altray/5.)) .OR. &
              altray .eq. 0. .OR. altray .eq. -2. .OR. abs(altray) .eq. 90.)then
               idebug = 1
               if(altray .eq. 90.)write(6,*)' zenith location (skyglow_phys)'
               write(6,*)' alt/ag/aa =',altray,ag,aa
               write(6,631)altray,ic,istart,iend
631            format(' call get_clr_src_dir for altitude ',f9.3,i4,2i9)
               if(altray .eq. 3.0)I4_elapsed = ishow_timer()
             else
               idebug = 0
             endif

             if(aod_vrt .gt. 5.0)then
!               od_a_scat = od_a_slant_a(ic,ialt)
                od_a_scat = aod_ref * ext_a(ic) 
                gscat = aod_asy(2,ic)
                if(idebug .gt. 0)then
                  write(6,*)
                  write(6,*)' Calling mscat_phase 1 for ialt/ic = ',ialt,ic,od_a_scat,gscat
                endif
!               We can still mimic the phase function from 'get_cld_pf' and 'ssa' from 'get_cloud_rad'
                rad = 0.
                call mscat_phase(od_a_scat,ssa(ic),gscat,idebug,rad &
                                ,p_order,rel_order,n_order,gmean,ssa_eff(ic))
                ssa_eff(ic) = ssa_eff(ic) ! * (0.5 + 0.5 * sind(altray))
             else
                ssa_eff(ic) = ssa(ic)
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
                htmsl,ssa_eff(ic),ag/ag_90,ao,aa_o_aa_90, &
                aod_ref,aero_refht,aero_scaleht, &
                ag_s/ag_90,aa_s_o_aa_90,ic,idebug, &
                istart,iend, &
                srcdir(ic),sumi_gc(ic),sumi_ac(ic),opac_slant, &
                nsteps,ds,tausum_a)

             if(l_solar_eclipse .eqv. .true.)then
               tausum_max = maxval(tausum_a)
               opac_max = opac(tausum_max)
               if(ic .eq. icd)then
                 write(6,632)altray,opac_slant,opac_max,istart,iend
632              format(' alt / opac_slant / opac_max = ',f9.3,2f11.6,2i9)
               endif
               do iopac = 1,nopac
!                opacmid(iopac) = opac_slant * (float(iopac) - 0.5) / float(nopac)
                 opacmid(iopac) = opac_max * (float(iopac) - 0.5) / float(nopac)
                 taumid(iopac) = -log(1.0 - opacmid(iopac))              
                 do i = 1,nsteps
                   io = i + ioffset
                   sbar = (float(io)-0.5) * ds
                   if(tausum_a(i) .ge. taumid(iopac))then
                       distecl(iopac,ic) = sbar
                       goto 641
                   endif
                   if(iopac .eq. nopac .and. ic .eq. icd .and. &
                      abs(altray) .le. 5. .and. i .eq. (i/10)*10 )then
                     write(6,64)ic,iopac,i,io,tausum_a(i),taumid(iopac) &
                               ,ds,sbar,distecl(iopac,ic)
64                   format('ic/iopac/i/io/tausum/mid/ds/sbar/dst',2i5,2i8 &
                           ,3f9.3,2f12.1)
                   endif
                 enddo ! i 
641              continue
               enddo ! iopac         

               if(ic .eq. icd .AND. (altray .eq. nint(altray) .OR. altray .le. 20.) )then
                write(6,65)ic,altray,srcdir(ic),srcdir(ic)/srcdir_90(ic) &
                          ,opac_slant,opacmid(1),opacmid(nopac),taumid(1),taumid(nopac),distecl(1,ic),distecl(nopac,ic)
65              format( &
                   ' ic/alt/srcdir/ratio/opacsl/opacmid/taumid/distecl:',i3,f8.3,2f9.3,f7.3,2(1x,2f16.3),2f12.1)
               endif

             else ! l_solar_eclipse = F
               if(ic .eq. icd .AND. (altray .eq. nint(altray) .OR. altray .le. 20.) )then
                write(6,66)altray,srcdir(ic),srcdir(ic)/srcdir_90(ic)
66              format(' alt/srcdir/ratio:',3f9.3)
               endif

             endif ! l_solar_eclipse
           enddo ! ic
         else
           sumi_gc(:) = 0.
           sumi_ac(:) = 0.
         endif ! solalt > 0

!        Determine aerosol multiple scattering order
!        altscat = 1.00 * altray + 0.00 * sol_alt
!        altscat = max(altray,sol_alt)
         altscat = sqrt(0.5 * (altray**2 + sol_alt**2))
         call get_airmass(altscat,htmsl,patm &         ! I
                         ,aero_refht,aero_scaleht &    ! I
                         ,earth_radius,0 &             ! I
                         ,agdum,aodum,aascat,refr_deg) ! O
!        scatter_order = 1.0 ! f(altray,sol_alt)
         scatter_order = max(aod_vrt*aascat,1.0)**1.0 ! 0.5-1.5
         scatter_order_t = 1.0
         write(6,*)' altscat/aascat/sco',altscat,aascat,scatter_order

!        Determine effective asymmetry parameter from multiple scattering
!        http://www.arm.gov/publications/proceedings/conf15/extended_abs/sakerin_sm.pdf

!        cosp_frac = max(1.-scatter_order,0.)

!        Get src term at selected azimuths (every 10 degrees)
!        How should dmintopo be set if the ray hits the earth outside the
!        domain? If htmin_view < 0. then l_dlow should perhaps be False.
         dmintopo = minval(dist_2_topo(ialt,:))
         dmaxtopo = maxval(dist_2_topo(ialt,:))
         if(dmintopo .eq. 0. .and. htmin_view .le. 130000. .and. &
                                   htmin_view .ge. 0.            )then 
!          Incomplete terrain in ring and ray grazes/exits atm
           l_dlow = .true.  
         else
!          Complete terrain in ring or ray doesn't graze/exit atm
           l_dlow = .false. 
         endif

         if((sol_alt * altray .lt. 100. .or. sol_alt .lt. 0. .or. altray .lt. 0.) .AND. (l_dlow .eqv. .true.) )then
           l_dlow2 = .true.
         else
           l_dlow2 = .false.
         endif

         if(altray .le. -10.0)then
           azi_d10 = 1.
         elseif(altray .le. 10.0)then
           azi_d10 = 5.
         else
           azi_d10 = 10.
         endif
         jazi_d10 = nint(azi_d10/azi_scale)

         write(6,67)altray,sol_alt,sol_alt*altray,dmintopo,htmin_view,l_dlow,l_dlow2,emis_ang
67       format(' altray/sol_alt/prod/dmintopo/htmin_view/l_dlow/l_dlow2/em',f8.4,f8.2,f8.1,2f12.0,2l2,f8.3)
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
                azi_ref = nint(sol_azi / azi_d10) * azi_d10

!               azi_ref = 175. ! volcanic test
!               azi_ref = 340. ! 300km test
!               azi_ref = 315. ! limb test  

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
                idebug_clr = idebug_a(ialt,jazi)

                if(idebug_clr .eq. 1)then
                  write(6,68)ialt,jazi,ic,sol_alt,altray,view_azi_deg &
                            ,dmintopo,dmaxtopo,ao
68                format(/' call   get_clr_src_dir_low',i6,2i5,3f8.2,2f9.0,f9.3)
!                 write(6,*)'sol_azi,azi_ref,jazi_d10',sol_azi,azi_ref,jazi_d10
                endif
                if(aa_90 .gt. 0.)then
                  aa_o_aa_90 = aa/aa_90
                  aa_s_o_aa_90 = aa_s/aa_90
                else
                  aa_o_aa_90 = 0.
                  aa_s_o_aa_90 = 0.
                endif

!               Determine relevant solar altitude along ray
                if(htmsl .gt. 100e3)then ! highalt strategy
                  if(altray .gt. 0.)then
                    refdist_solalt = 0.
                    solalt_ref = sol_alt
                  else ! refdist_solalt is hopefully generalized enough
!                   distance to limb
                    refdist_solalt = sqrt( (earth_radius+htmsl)**2 &
                                          - earth_radius**2       )
                    if(trace_solalt(ialt,jazi) .ne. sol_alt)then ! valid trace
                      solalt_ref = trace_solalt(ialt,jazi)
                    else ! angular distance to htmin
                      gcdist_km = (horz_dep * rpd * earth_radius) / 1000.
                      call RAzm_Lat_Lon_GM(rlat,rlon,gcdist_km &
                             ,view_azi_deg,rlat_limb,rlon_limb,istatus)
                      solzen_ref = angdist(rlat_limb,sol_lat &
                                          ,rlon_limb-sol_lon)
                      solalt_ref = 90. - solzen_ref
                    endif
                    if(idebug_clr .eq. 1)then
                      write(6,681)rlat,rlon,sol_lat,sol_lon,rlat_limb,rlon_limb,gcdist_km
681                   format(' lat-lon obs/sol/limb = ',6f8.2,' gc = ',f9.0)
                      write(6,*)'salt/tr/ref-1',sol_alt,trace_solalt(ialt,jazi),solalt_ref
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
                  od_a_vert,ext_ha(ic),htmsl,ssa_eff(ic), &
                  ag/ag_90,ao,aa_o_aa_90, &
                  aod_ref,aero_refht,aero_scaleht, &
                  ag_s/ag_90,aa_s_o_aa_90, &
                  ags_a,aas_a,isolalt_lo,isolalt_hi,del_solalt,ic,idebug_clr, &
                  refdist_solalt,solalt_ref, &
                  srcdir_a(ic,jazi),sumi_g_a(ic,jazi),sumi_a_a(ic,jazi),&
                  opac_slant,nsteps_low,dsl,tausum_a)
                if(idebug_clr .eq. 1)then
                  write(6,69)ialt,jazi,ic,sol_alt,trace_solalt(ialt,jazi),altray &
                            ,srcdir(ic),srcdir_a(ic,jazi)
69                format(' called get_clr_src_dir_low',i7,2i5,3f8.2,2f9.4/)
                endif

              else ! use generic values for azimuths every 10 degrees
                srcdir_a(ic,jazi) = srcdir(ic)
                sumi_g_a(ic,jazi) = sumi_gc(ic)
                sumi_a_a(ic,jazi) = sumi_ac(ic)

              endif
            enddo ! ic

          endif ! .true.
         enddo ! jazi

         if(altray .eq. 90. .and. (l_solar_eclipse .eqv. .true.))then
          I4_elapsed = ishow_timer()
          write(6,*)' get eclipse parms at selected azimuths'
         endif

!        Get eclipse parms at selected azimuths
         jazi_dec = 2
         if(altray .ge. 5.)then
          jazi_dec = max(jazi_dec,nint(1./azi_scale))
         endif

         if(l_solar_eclipse .eqv. .true.)then
          write(6,*)' eclipse parms for alt/jazi_dec',altray,jazi_dec
         endif

         do jazi = jazi_start,jazi_end,jazi_dec
          view_azi_deg = view_az(ialt,jazi)
          if(sol_alt .gt. 0.)then

!           Get eclipse parms if applicable (and effective brightness)
            if(l_solar_eclipse .eqv. .true.)then
              do ic = 1,nc
!               Compute/Interpolate weighted intensity using distecl
!               if( ((jazi-1)/2)*2+1 .eq. jazi)then
                if(mod(jazi,2) .eq. 0)then ! even jazi
                  eobsc_sum(ic) = 0.
                  emag_sum = 0.
                  if(abs(altray) .eq. 90. .and. jazi .eq. jazi_start)then
                    iverbose = 2
                  else
                    iverbose = 0
                  endif
                  do iopac = 1,nopac
                    write(6,*)
                    if(iverbose .eq. 2)write(6,*)' Calling sun_eclipse_parms at zenith, iopac = ',iopac
                    call sun_eclipse_parms(i4time,rlat,rlon,htmsl &
                     ,iverbose,altray,view_azi_deg,distecl(iopac,ic) &
                     ,earth_radius,elgms,emag,eobscf,eobsc_a(iopac,ic))
                    eobsc_sum(ic) = eobsc_sum(ic) + eobsc_a(iopac,ic)
                    emag_sum  = emag_sum  + emag
                  enddo

                  eobsc(ic,jazi) = eobsc_sum(ic) / float(nopac)
                  emag_a(jazi)  = emag_sum  / float(nopac)
                  if(iverbose .eq. 2)write(6,*)' eobsc(ic,jazi) = ',ic,jazi,eobsc(ic,jazi)
!               else ! fill in from previous azimuth
!                 jazim1 = max(jazi-1,jazi_start)
!                 eobsc(ic,jazi) = eobsc(ic,jazim1)
!                 emag_a(jazi) = emag_a(jazim1)
                endif
              enddo ! ic
            endif
          endif ! sun above horizon
         enddo ! jazi

         if(altray .eq. 90. .and. (l_solar_eclipse .eqv. .true.))then
          I4_elapsed = ishow_timer()
         endif

         do jazi = jazi_start,jazi_end ! main azimuth loop
          altray = view_alt(ialt,jazi)
          view_altitude_deg = altray
          view_azi_deg = view_az(ialt,jazi)

!         Interp srcdir(:), sumi_gc(:), sumi_ac(:)
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

          if(solalt_limb_true .gt. twi_0)then ! formerly 0.

!           Get eclipse parms if applicable (and effective brightness)
            if(l_solar_eclipse .eqv. .true.)then

!             Interp eobsc(:), emag_a
              if(jazi .lt. jazi_end)then
                jazi_interpl = ((jazi-jazi_start)/jazi_dec) * jazi_dec + jazi_start
              else
                jazi_interpl = jazi_end - jazi_dec
              endif
              jazi_interph = jazi_interpl + jazi_dec 
              fazi = (float(jazi) - float(jazi_interpl)) / float(jazi_dec)

              eobsc(:,jazi) = (1.-fazi) * eobsc(:,jazi_interpl) &
                            +     fazi  * eobsc(:,jazi_interph)

              emag_a(jazi) = (1.-fazi) * emag_a(jazi_interpl) &
                           +     fazi  * emag_a(jazi_interph)

              do ic = 1,nc
                if(idebug .ge. 1)then
                  idebuge = 2
                else
                  idebuge = 0
                endif

!               Compute/Interpolate weighted intensity using distecl
!               if( (jazi/jazi_dec)*jazi_dec .eq. jazi)then ! even jazi
!                 continue
!               else ! fill in from previous azimuth
!                 jazim1 = max(jazi-1,jazi_start)
!                 eobsc(ic,jazi) = eobsc(ic,jazim1)
!                 emag_a(jazi) = emag_a(jazim1)
!               endif

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
                      ,emag_a(jazi),eobsc_a(1,ic),eobsc_a(nopac,ic) &
                      ,eobsc(ic,jazi),ecl_scat(ic),ecl_intd(ic) &
                      ,ecl_dir_rat(ic)
70                format('eclipse altray/dist/elg/emag',f9.2,2f12.1,f9.4,f8.4, &
                         ' eobsc',3f7.4,' scat/int/dir',f11.7,f10.7,f7.4)
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

!           HG illumination
            if(aod_vrt .lt. 5.0)then
              do ic = 1,nc
!               Check assignments in 'mem_namelist.f90'
                fb = aod_asy(3,ic)**scatter_order
                g1 = aod_asy(2,ic)**scatter_order
                g2 = aod_asy(1,ic)**scatter_order
                hg2(ic) = dhg2(elong(ialt,jazi),fb,fcterm2)

!               topo phase function assumes scatter order is non-topo 
!               value (for now)
                hg2t(ic) = hg2(ic)
              enddo 

              ssa_eff(:) = ssa(:)

            else ! use high OD routine calculations
              do ic = 1,nc
!               Check assignments in 'mem_namelist.f90'
                fb = aod_asy(3,ic)**scatter_order
                g1 = aod_asy(2,ic)**scatter_order
                g2 = aod_asy(1,ic)**scatter_order
                hg2(ic) = dhg2(elong(ialt,jazi),fb,fcterm2)

!               topo phase function assumes scatter order is non-topo 
!               value (for now)
                hg2t(ic) = hg2(ic)
 
!               od_a_scat = od_a_slant_a(ic,ialt)
                od_a_scat = aod_ref * ext_a(ic) 
                gscat = aod_asy(2,ic)
                if(idebug .gt. 0)then
                  write(6,*)
                  write(6,*)' Calling mscat_phase 2 for ialt/ic = ',ialt,ic,od_a_scat,gscat
                endif
                call mscat_phase(od_a_scat,ssa(ic),gscat,idebug,rad &
                                ,p_order,rel_order,n_order,gmean,ssa_eff(ic))
                ssa_eff(ic) = ssa_eff(ic) ! * (0.5 + 0.5 * sind(altray))
              enddo 

            endif

            if(aod_ray(ialt,jazi) .gt. 0.)then
                aod_dir_rat = aod_ray_dir(ialt,jazi) / aod_ray(ialt,jazi)
            else
                aod_dir_rat = 0.
            endif

            if(aod_dir_rat .gt. 1.0001 .OR. aod_dir_rat .lt. 0.0)then
              write(6,*) &
            ' WARNING in skyglow_phys: resetting out of bounds aod_dir_rat',aod_dir_rat
              write(6,*)' ialt/jazi/altazray = ',ialt,jazi,altray,view_azi_deg
              write(6,*)' aod_ray_dir/aod_ray = ',aod_ray_dir(ialt,jazi),aod_ray(ialt,jazi)
!             stop
              aod_dir_rat = 1.0
            endif

!           Consider arg for sideways scattering from downward diffuse
!           is ~hg(90.) or an evaluated integral. Use 0.5 for now.
            if(l_solar_eclipse .eqv. .false.)then
              hg2d = (hg2  * aod_dir_rat) + (0.5 * (1.0 - aod_dir_rat)) 
              hg2t = (hg2t * aod_dir_rat) + (0.5 * (1.0 - aod_dir_rat)) 
            else
              do ic = 1,nc
                arg = aod_dir_rat * ecl_dir_rat(ic)
                hg2d(ic) = (hg2(ic) * arg) + (0.5 * (1.0 - arg)) 
              enddo ! ic
            endif

            rayleigh_pfunc = rayleigh_pf(elong(ialt,jazi))
            rayleigh_gnd = rayleigh_pfunc + sfc_alb_c(2) * sind(sol_alt)

!           Raising this further will reduce the use of mode_sky = 3
            if(htmsl .le. 1000e3 .and. dist_2_topo(ialt,jazi) .eq. 0d0)then 
              mode_sky = 4 ! simple approach (<250km non-topo case)
              do ic = 1,nc
                day_inte = day_int0 * (1.0 - eobsc(ic,jazi))

                if(airmass_g .gt. 0.)then
                  gasf = ag_2d(ialt,jazi) / airmass_g
                else
                  gasf = 1.0
                endif

!               od_a = aod_ref * aa * ext_a(ic)  ! slant od
                od_a_slant = od_a_vert * aa_o_aa_90 
                od_a = od_a_slant
                if(od_a_slant .gt. 0. .and. sol_alt .gt. twi_alt .and. &
                  aod_ill_opac_potl(ialt,jazi) .gt. 0.)then ! normalize with extinction?
!                 aodf = min(aod_ill(ialt,jazi)/od_a_slant,1.0)
                  aodfo = min(aod_ill(ialt,jazi)/aod_tot(ialt,jazi),1.0)
!                 May be large when viewer is in stratosphere
                  aodf = aod_ill_opac(ialt,jazi)/aod_ill_opac_potl(ialt,jazi)
                  aodf = min(aodf,1e8)
                else
                  aodfo = 1.0
                  aodf = 1.0
                endif

                if(solalt_ref .gt. twi_alt)then ! allow crepuscular rays?
!               if(solalt_ref .gt. 0.0)then     ! allow only sun rays
                  radf = clear_radf_c(ic,ialt,jazi)
                else
                  radf = 1.0
                endif

!               We can still add in gasf and aodf
                clear_rad_c(ic,ialt,jazi) = day_inte * &
!                   (rayleigh_pfunc * sumi_gc(ic) * gasf + &
                    (rayleigh_gnd   * sumi_gc(ic) * gasf * radf + &
                     hg2d(ic)       * sumi_ac(ic) * aodf   )

                if(idebug .ge. 1 .OR. clear_rad_c(ic,ialt,jazi) .gt. 1e15)then
!                 write(6,*)'solalt_ref/twi_alt/cradf ',solalt_ref,twi_alt,clear_radf_c(ic,ialt,jazi)
                  write(6,71)ic,day_inte,elong(ialt,jazi) &
                      ,sumi_gc(ic),sumi_ac(ic),gasf,radf,aodf,aodfo &
                      ,rayleigh_pfunc,hg2(ic),aod_dir_rat,hg2d(ic) &
                      ,clear_rad_c(ic,ialt,jazi)
71                format('day_inte/elong/sumi_g/sumi_a/gasf/radf/aodf/aodfo', &
                         '/rayleigh/hg2/aodr/hg2d/clrrd4', &
                         i2,f12.0,f8.3,2x,2f11.8,2x,3f8.3,2x,5f8.3,f12.0)      
                endif
              enddo ! ic

!           The first option has a discontinuity around -15 degrees above 20km
!           The second option has a brightness jump around 18000m
!           elseif(htmsl .gt. 18000. .and. dist_2_topo(ialt,jazi) .eq. 0.)then ! high vantage point
            elseif(dist_2_topo(ialt,jazi) .eq. 0d0)then ! non-topo (high vantage point)
!           elseif(htmsl .gt. 18000.)then ! high vantage point

              mode_sky = 3 ! high vantage point >250km (non-topo)

!             Consider case where 'dist_2_topo' is zero (in a cloud) and
!             'topo_solalt' wasn't set. In this case 'day_inte' wouldn't 
!             be properly determined. Here we use info from a new array 
!             called 'trace_solalt'.
!             solalt_atmos = trace_solalt(ialt,jazi) 
!             if(solalt_atmos .gt. 0.)then
!               sol_occ = 1.0 ! visible
!             else
!               sol_occ = 1.0 ! 0.0 ! invisible
!             endif

!             Use/restore simpler approach for downward rays from high alt
              do ic = 1,nc
!               Rayleigh illumination
                day_inte = day_int0 * (1.0 - eobsc(ic,jazi))
                od_per_am = ext_g(ic)
                rayleigh = brtf(airmass_g,od_per_am) * rayleigh_gnd
                od_a = aod_ref * aa * ext_a(ic)  ! slant od
                brta = opac(od_a) * hg2d(ic)

!               Total illumination
                if(htmin_view .le. 100e3)then
!                 clear_rad_c(ic,ialt,jazi) = day_inte * (rayleigh + brta)
                  clear_rad_c(ic,ialt,jazi) = day_inte * &
                      (rayleigh       * sumi_gc(ic) + &
                       hg2d(ic)       * sumi_ac(ic)   )
                else
                  clear_rad_c(ic,ialt,jazi) = 0.
                endif

                if(idebug .ge. 1 .and. ic .eq. icd)then
                  write(6,72)day_inte,elong(ialt,jazi),airmass_g,brtf(airmass_g,od_per_am) &
                            ,rayleigh,hg2(ic),hg2d(ic),sumi_gc(ic),sumi_ac(ic),clear_rad_c(2,ialt,jazi)
72                format('day_inte/elong/ag/brtf/rayleigh', &
                         '/hg2/hg2d/sumi/clrrd3',f12.0,6f8.3,2f10.5,f12.0)      
                endif

               enddo ! ic

               if(idebug .ge. 1)then
                 write(6,73)day_inte,elong(ialt,jazi),rayleigh,hg2d(2),clear_rad_c(:,ialt,jazi)
73               format('day_inte/elong/hg2/clear_rad :',f12.0,3f8.3,3f12.0)      
               endif

            else ! assume two scattering layers (g+a [lower-1], g [upper-2])
              do ic = 1,nc
                gasfrac = ag_2d(ialt,jazi) / airmass_g ! topo illuminated fraction of gas
                if(dist_2_topo(ialt,jazi) .eq. 0d0)then    ! free of topo
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
!                   if(altray .eq. -73.)write(6,*)'debug check 3',altray,jazi,idebug_a(ialt,jazi),idebug_topo
!                   if(altray .le. -89.07 .and. altray .ge. -89.40 .and. jazi .eq. 2048)then ! limb test
                    if(altray .le. -89.06 .and. altray .ge. -89.40 .and. jazi .eq. 2104)then ! limb test
                      if(ic .eq. icd)idebug_topo = 1
                    endif

                    if(aa_90 .gt. 0.)then
                      aa_o_aa_90 = aa/aa_90
                    else
                      aa_o_aa_90 = 0.
                    endif

                    if(htmsl .gt. 100e3)then ! highalt strategy

!                     Strategy from above dir_low section
                      if(altray .gt. 0.)then
                        refdist_solalt = 0.
                        solalt_ref = sol_alt
                      else ! is refdist_solalt generalized enough?
!                       distance to limb
                        refdist_solalt = sqrt( (earth_radius+htmsl)**2 &
                                              - earth_radius**2       )
                        if(trace_solalt(ialt,jazi) .ne. sol_alt)then ! valid trace
                          solalt_ref = trace_solalt(ialt,jazi)
                        else ! should consider angular distance to surface
                          solalt_ref = sol_alt - (altray*cosd(view_azi_deg-sol_azi))
                          solalt_ref = min(solalt_ref,+180.-solalt_ref)
                          solalt_ref = max(solalt_ref,-180.-solalt_ref)
                        endif
                        if(idebug_clr .eq. 1)then
                          write(6,*)'salt/tr/ref-2',sol_alt,trace_solalt(ialt,jazi),solalt_ref
                        endif
                      endif

!                     Alternative strategy
                      refdist_solalt = dist_2_topo(ialt,jazi)
!                     solalt_ref = topo_solalt(ialt,jazi)

                    else
                      refdist_solalt = 0.
                      solalt_ref = sol_alt
                    endif

                    od_o_msl = (o3_du/300.) * ext_o(ic)

                    call get_clr_src_dir_topo(sol_alt,sol_azi, &        ! I
                     altray,view_azi_deg,emis_ang,r_missing_data, &     ! I
                     ext_g(ic),od_g_vert,od_o_msl,od_a_vert, &          ! I
                     htmsl,dist_2_topo(ialt,jazi), &                    ! I
                     ssa_eff(ic),ag/ag_90,aa_o_aa_90, &                 ! I
                     aod_ref,aero_refht,aero_scaleht, &                 ! I
                     ags_a,aas_a, &                                     ! I
                     isolalt_lo,isolalt_hi,del_solalt,ic,idebug_topo, & ! I
                     nsteps_topo,refdist_solalt,solalt_ref, &           ! I
                     sumi_gct(ic),sumi_act(ic),opac_slant,tausum_a)     ! O
  
                    if(idebug_topo .gt. 0)then
                      write(6,81)ialt,jazi,ic,sol_alt,altray &
                         ,dist_2_topo(ialt,jazi),sumi_gct(ic),sumi_act(ic)
81                    format(' called get_clr_src_dir_topo',i7,i5,i4,2f8.2,f12.0,2f9.4)
                    endif ! write log info
                  endif ! call get_clr_src_dir_topo

!                 Normalize by extinction?
                  if(aod_2_topo(ialt,jazi) .gt. 0.)then
                    aodfo = aod_ill(ialt,jazi) / aod_2_topo(ialt,jazi) 
                    aodfo = min(aodfo,1e8)
                  else
                    aodfo = 1.0
                  endif

                  if(aod_ill_opac_potl(ialt,jazi) .gt. 0.)then
                    aodf = aod_ill_opac(ialt,jazi)/aod_ill_opac_potl(ialt,jazi)
                    aodf = min(aodf,1e8)
                  else
                    aodf = 1.0
                  endif

!                 if(htmsl .gt. 500e3 .and. solalt_ref .lt. 0.)then
                  if(htmsl .gt. thr_abv_clds .and. solalt_ref .lt. twi_alt)then
                    aodf = 1.0
                    radf = 1.0
                  else ! use original values
                    radf = clear_radf_c(ic,ialt,jazi)
                  endif

!                 If we're hitting topo during an eclipse with a low vantage point
!                 and we're looking above the limb this may need to trend closer
!                 to the observer obscuration.
                  if(htmsl .le. 10000. .and. altray .ge. 0.)then  
                    day_inte = day_int0 * (1.0 - eobsl)
                  else ! consider eobsc determined earlier
                    day_inte = day_int0 * (1.0 - eobsc(ic,jazi))
                  endif

                  clear_rad_topo = day_inte * &
                    (sumi_gct(ic) * rayleigh_gnd * radf + &
                     sumi_act(ic) * hg2t(ic)     * aodf) 
                  if(idebug .ge. 1 .AND. ic .eq. icd)then
                    write(6,85)ialt,jazi,ic,clear_rad_topo,sumi_gct(ic),sumi_act(ic)
85                  format(' computed clear_rad_topo =  ',i7,i5,i4,3e15.7)
                  endif
!                 if(l_solar_eclipse .eqv. .false.)then ! old experimental
                    clear_rad_c(ic,ialt,jazi) = clear_rad_topo
!                 endif
                  elglim = 5.0 ! catch aureole on topo
                  mode_sky = 2 ! hitting topo case

                  if(idebug .ge. 1 .AND. ic .eq. icd)then
                    write(6,91)elong(ialt,jazi),sumi_gct(ic),rayleigh_gnd &
                       ,radf,sumi_act(ic),hg2t(ic) &
                       ,aodf,aodfo,ssa_eff(ic),clear_rad_topo      
91                  format('elg/ig/rayg/radf/ia/hg2t/aodf/aodfo/ssa/clear_rad_topo =',f7.2,3x,3f9.4,3x,4f9.4,f7.2,3x,f12.0)
                  endif

                endif ! hitting topo

                if(idebug .ge. 1 .AND. (ic .eq. icd .or. &
                        abs(elong(ialt,jazi)-90.) .le. 0.5 .or. &
                            elong(ialt,jazi)      .le. elglim)        )then
                  jazim1 = max(jazi-1,minazi)
                  write(6,82)ic,jazi,day_inte,eobsc(ic,jazim1:jazi),srcdir(ic),srcdir_90(ic)
82                format('day_inte/eobsc/srcdir/90',2i5,f12.0,4f8.4)           
!                 write(6,83)day_inte,elong(ialt,jazi),airmass_g,od_g &
!                           ,aod_ray(ialt,jazi),aa &
!                           ,od_a,alphav_g*1e3,alphav_a*1e3,od_g1,od_g2 &
!                           ,clear_rad_2lyr            
83                format(&
     'day_inte/elg/ag/od_g/aod_ray/aa/od_a/alphav_g/alphav_a/od_g1/od_g2/clr_rad :' &
                        ,f12.0,f7.2,2f7.3,f8.3,2f7.3,1x,2f7.3,1x,2f8.4,f12.0)      
                  write(6,84)altray,view_azi_deg,am_sun,aascat,scatter_order,hg2(ic),hg2d(ic),gasfrac,aod_ill(ialt,jazi) & 
                            ,aod_ray_dir(ialt,jazi),aod_ray(ialt,jazi),aod_dir_rat
84                format('altaz/amsun/aasc/sco/hg2/hg2d/gasfrac/aodill/dir/ray/rat = ',2f8.2,5f9.4,2f10.6,3f7.3)                  
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
111           format('alt/azi/salt-tp-tr/mode/od_a/dst/aodf/aodfo/clrrd' &
                    ,f8.3,f6.1,1x,3f7.2,i3,f12.5,f12.0,2f9.4,3e14.5)

              if(clear_rad_c(1,ialt,jazi) .lt. 0. .or. clear_rad_c(1,ialt,jazi) .gt. 1e20)then
                write(6,*)' ERROR clrrd(1) out of bounds',clear_rad_c(1,ialt,jazi)
!               stop
              endif

              if(dist_2_topo(ialt,jazi) .le. 0. .and. viewalt_limb_true .le. 0.)then
                write(6,*)' ERROR: dist_2_topo/viewalt_limb =',dist_2_topo(ialt,jazi),viewalt_limb_true
                stop
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

!       Apply during daylight eclipses or during any sun / twilight
        if( ( (sol_alt .ge. 0. .and. (l_solar_eclipse .eqv. .true.)) .OR. &
            (sol_alt .lt. 100. .and. sol_alt .gt. twi_0) ) .AND. & 
                                                  htmsl .le. thr_abv_clds )then
            write(6,*)'Calling get_sky_rad_ave'
            do ic = 1,nc ! secondary scattering in each color
!               Good temporal variability, constant for spatial
                call get_sky_rad_ave(clear_rad_c(ic,:,:) &
                    ,view_alt,view_az,maxalt-minalt+1,maxazi-minazi+1 &
                    ,sky_rad_ave(ic))
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

!               clear_rad_ref = clear_rad_c0(ic,maxalt,minazi) ! eclipse
                clear_rad_ref = clear_rad_c(ic,maxalt,minazi)  ! no eclipse

                do ialt = ialt_start,ialt_end
                  altray = view_alt(ialt,jazi_start)
                  arg = sind(min(max(altray,1.5),90.))
                  od_slant = od_vert / arg                      
                  if(ic .eq. icd)then
                    write(6,*)'alt,od_slant,clear_rad_c',altray,od_slant,clear_rad_c(ic,ialt,jazi_dbg)
                  endif

!                 sphere ave, any gnd
                  sph_rad_ave = sky_rad_ave(ic) * 0.5 * (1. + sfc_alb_c(ic))

                  do jazi = jazi_start,jazi_end
                    if(dist_2_topo(ialt,jazi) .eq. 0d0)then ! unobstructed by terrain
                      od_slant = od_vert / arg                      
                      fracg = od_g_vert / od_vert
                      fraca = od_a_vert / od_vert
!                     sky_rad_scatg = sph_rad_ave * viewalt_adjust & 
!                                   * fracg * opac(od_slant) * highalt_adjust
!                     sky_rad_scata = (0.4 * sph_rad_ave * viewalt_adjust + 0.6 * clear_rad_c(ic,ialt,jazi)) & 
!                                   * fraca * opac(od_slant) * highalt_adjust
                    else ! topo in path
                      od_slant_g = ag_2d(ialt,jazi) * ext_g(ic)
                      od_slant_a = aod_2_topo(ialt,jazi)
                      od_slant = od_slant_g + od_slant_a
!                     sky_rad_scatg = sph_rad_ave * opac(od_slant_g) * highalt_adjust * viewalt_adjust
!                     sky_rad_scata = sph_rad_ave * opac(od_slant_a) * highalt_adjust * viewalt_adjust
                    endif
!                   sky_rad_scat(ic,ialt,jazi) = sky_rad_scatg + sky_rad_scata
                    sky_rad_scat(ic,ialt,jazi) = sph_rad_ave * opac(od_slant) * highalt_adjust * viewalt_adjust
                  enddo ! jazi

                  clear_rad_c(ic,ialt,:) = clear_rad_c(ic,ialt,:) + sky_rad_scat(ic,ialt,:)
                  if(ialt .eq. ialt_end)then ! zenith/high point info
                    scatfrac = sky_rad_scat(ic,maxalt,jazi_dbg) / clear_rad_c(ic,maxalt,jazi_dbg)
!                   ri_over_i0 = clear_rad_c(ic,maxalt,jazi_dbg) / clear_rad_c0(ic,maxalt,jazi_dbg)
!                   ri_over_f0 = clear_rad_c(ic,maxalt,jazi_dbg) / day_int0                            
                    rmaglim = b_to_maglim(clear_rad_c(2,ialt_end,jazi_dbg))
                    write(6,201)ic,altray,view_az(ialt_end,jazi_dbg),sph_rad_ave,od_vert,clear_rad_ref,sky_rad_scat(ic,maxalt,jazi_dbg),clear_rad_c(ic,maxalt,jazi_dbg),scatfrac
201                 format(' ic/alt/az/sph_rad_ave/od/c0/scat/rad/rat/scatf',i2,2f9.2,f11.0,f6.3,f11.0,2f11.0,f7.4,2f10.7)
                    if(ic .eq. icd)then
                      b = clear_rad_c(2,ialt_end,jazi_dbg)
                      write(6,211)b,b_to_v(b),rmaglim
                    endif
211                 format(' b/mag/rmaglim at high point (second scat) =' &
                          ,f10.1,2f9.2)
                  endif ! at zenith 
                enddo ! ialt
            enddo ! ic
        else
            sky_rad_scat(:,:,:) = 0.
        endif

        write(6,*)' clearrad2:',clear_rad_c(2,(minalt+maxalt)/2,minazi:maxazi:10)
        write(6,*)' returning from skyglow_phys'

        return
        end
        

        subroutine get_ray_info(alt,htmsl,htsfc,earth_radius,alt_norm &
                               ,r_missing_data,dist_to_sfc)

        include 'trigd.inc'

        real alt                  ! I elevation angle
        real htmsl                ! I observer height MSL
        real htsfc                ! I sfc (or layer) height 
        real earth_radius         ! I earth radius (meters)
        real r_missing_data       ! I
        real alt_norm             ! O elevation angle rel to sfc normal
        real dist_to_sfc          ! O distance to sfc (meters)

!       Altitude relative to surface normal (emission angle)
        if(alt .ne. 0.)then
          slope = tand(90. - abs(alt))
          htrad = (earth_radius+htmsl) / (earth_radius+htsfc)
          c = htrad * slope
          call line_ellipse(slope,c,0.,0.,1.,1.,r_missing_data,x1,x2,y1,y2)
          if(y1 .ne. r_missing_data)then
            gc = atan2d(+y1,-x1) ! great circle dist to ground pt
            if(alt .gt. 0.)then
              alt_norm = alt - gc
            else
              alt_norm = alt + gc
            endif
            distrad = sqrt((htrad+x1)**2 + y1**2)
            dist_to_sfc = distrad * (earth_radius+htsfc)
          else
            dist_to_sfc = r_missing_data
          endif

!         write(6,1)distrad,htrad,x1,y1
 1        format(' distrad,htrad,x1,y1',4f13.8)
        else
          alt_norm = 0.  
        endif

        return
        end
