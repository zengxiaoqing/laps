
     subroutine get_clr_src_dir(solalt,viewalt,od_g_msl,od_g_vert,od_a_vert,ext_haero,htmsl,ssa,agv,ao,aav,aod_ref,redp_lvl,scale_ht_a,ags,aas,ic,idebug,istart,iend,srcdir,sumi_g,sumi_a,opac_slant,nsteps,ds,tausum_a)

     use mem_namelist, ONLY: earth_radius
     use constants_laps, ONLY: R, GRAV

     include 'trigd.inc'
     include 'rad.inc'

!    Input idebug for select angles

     trans(od) = exp(-min(od,80.))
     opac(od) = 1.0 - trans(od)
     scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5
     real*8 dsdst,dradius_start,daltray,dcurvat2
     dcurvat2(dsdst,dradius_start,daltray) = &
              sqrt(dsdst**2 + dradius_start**2 &
       - (2D0*dsdst*dradius_start*cosd(90D0+daltray))) - dradius_start      

     real tausum_a(istart:iend)
     real*8 sbar,htbar

     real agv      ! gas   relative to zenith value for observer ht/viewalt
     real ao       ! ozone relative to zenith at sea level pressure/viewalt
     real aav      ! aero  relative to zenith value for observer ht/viewalt
     real ags      ! gas   relative to zenith value for observer ht/solalt
     real aos      ! ozone relative to zenith value for observer ht/solalt
     real aas      ! aero  relative to zenith value for observer ht/solalt

     od_o_msl = (o3_du/300.) * ext_o(ic)
     od_o_vert = od_o_msl * patm_o3(htmsl)

!    htmin_view = htminf(htmsl,viewalt,earth_radius)
!    patm_htmsl = ztopsa(htmin_view) / 1013.25
     patm_htmsl = ztopsa(htmsl) / 1013.25
     rearg = 2.0 + 0.667 * patm_htmsl

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * agv        
     od_o_slant = od_o_msl  * ao
     od_a_slant = od_a_vert * aav        
     od_slant = od_g_slant + od_o_slant + od_a_slant
     od_slant = max(od_slant,1e-30) ! QC check
     r_e = 6371.e3

     scale_ht_g = R * ztotsa(htmsl) / GRAV
     alpha_sfc_g = od_g_vert / scale_ht_g

     scale_ht_g = R * ztotsa(0.) / GRAV
     alpha_ref_g = od_g_msl  / scale_ht_g

     alpha_sfc_a = od_a_vert / scale_ht_a
     alpha_ref_a = aod_ref / scale_ht_a

     opac_vert = opac(od_g_vert + od_o_vert + od_a_vert)
     if(opac_vert .eq. 0)then ! very small od
         opac_vert = od_g_vert + od_o_vert + od_a_vert
     endif

     opac_slant = opac(od_slant)
     if(opac_slant .eq. 0)then ! very small od
         opac_slant = od_slant
     endif

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
!    ds = 25.

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'solalt/viewalt = ',solalt,viewalt
         write(6,*)'ags/aas = ',ags,aas
         write(6,*)'agv/aav = ',agv,aav
         write(6,*)'od_g_vert = ',od_g_vert
         write(6,*)'scale_ht_g = ',scale_ht_g
         write(6,*)'od_o_msl/vert = ',od_o_msl,od_o_vert
         write(6,*)'od_o_vert = ',od_o_vert
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_o_slant,ao = ',od_o_slant,ao
!        if(viewalt .ge. 0.)then
           write(6,*)'aorel = ',od_o_slant/od_o_vert
!        endif
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
         write(6,*)'ds/nsteps = ',ds,nsteps
     endif ! i

     tausum = 0.
     tausum_a(:) = 0.
     sumi = 0. ; sumi_g = 0.; sumi_a = 0.
     frac_opac = 0.

     httopill = 100000.

     if(idebug .eq. 1)then
       write(6,*)'           sbar         htbar       dtaug    dtauo   tausum od_sol_a  od_solar   rad       di     sumi_g   sumi_a opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     if(solalt .ge. 0.)then
         frac_iso_a = 1.0
     else
         frac_iso_a = 0.5
     endif

     do i = istart,iend
         sbar = (float(i)-0.5) * ds
         if(htmsl .le. 1000e3)then ! low altitude
           htbar = sbar * sind(viewalt)
           xybar = sbar * cosd(viewalt)
!          htbar = htbar + curvat(xybar,earth_radius*(4./3.)) ! Earth Curvature
           htbar = htbar + xybar**2 / (rearg * r_e) ! Earth Curvature
         else                      ! high altitude
           radius_start_eff = earth_radius + htmsl      
!          Is this showing up negative at nadir?
           htbar = dcurvat2(sbar,dble(radius_start_eff),dble(viewalt)) 
         endif
         htbar_msl = htbar+htmsl

         zapp = 90. - solalt
         od_solar_slant_o3 = od_o_msl * patm_o3(htbar_msl) * airmasso(zapp,htbar_msl)

         patm_bar = ztopsa(htbar_msl) / 1013.25
         scale_ht_g = R * ztotsa(htbar_msl) / GRAV
         alphabar_g = (od_g_msl * patm_bar) / scale_ht_g

!        od_solar_slant_g = od_g_vert * ags * exp(-htbar/scale_ht_g)
         od_solar_slant_g = od_g_msl  * ags * exp(-htbar_msl/scale_ht_g)
!        od_solar_slant_a = od_a_vert * aas * exp(-htbar/scale_ht_a) 
         od_solar_slant_a = aod_ref   * aas * exp(-(htbar_msl-redp_lvl)/scale_ht_a) * frac_iso_a
         od_solar_slant = od_solar_slant_g &
                        + od_solar_slant_o3 + od_solar_slant_a

         alphabar_o = alpha_o3(od_o_msl,htbar_msl)
!        alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         alphabar_a = alpha_ref_a * exp(-(htbar_msl-redp_lvl)/scale_ht_a)

!        AOD stratospheric = .012 (nighttime - layer extends above observer)
         if(htbar_msl .gt. h1_ha .and. htbar_msl .le. h2_ha .and. solalt .lt. 333. .and. aod_ha .gt. 0.)then
             alphabar_a = alphabar_a + alpha_ha * ext_haero
         elseif(htbar_msl .gt. h2_ha .and. htbar_msl .lt. h3_ha)then
             frac_ha = (h3_ha - htbar_msl) / (h3_ha - h2_ha)
             alphabar_a = alphabar_a + alpha_ha * ext_haero * scurve(frac_ha)
         endif

         dtau = ds * (alphabar_g + alphabar_o + alphabar_a)
         tausum = tausum + dtau
!        if(tausum .lt. 1.0)distod1 = sbar
         tausum_a(i) = tausum

         if(htbar .gt. httopill .or. htbar .lt. -(htmsl+500.) .or. frac_opac .gt. 11.4)then
           if(idebug .eq. 1)then
!            write(6,*)' criteria ',htbar,httopill,-(htmsl+500.),frac_opac
             write(6,11)i,sbar,htbar,alphabar_g*ds,alphabar_o*ds,tausum,od_solar_slant_a,od_solar_slant,rad,di,sumi_g,sumi_a,opac_curr,frac_opac,sumi_mean,sumi_extrap
           endif
           goto 900
         endif

         opac_last = opac_curr
         opac_curr = opac(tausum)
         if(opac_curr .eq. 0.)opac_curr = tausum ! very small od
         do = opac_curr - opac_last

!        Experiment with in tandem with 'skyglow_phys'
!        rad = 1.0
         rad = trans(od_solar_slant)
         radg = max(rad,.0001) ! secondary scattering
         rada = max(rad,.00001) ! secondary scattering
         rad = max(rad,.0001) ! secondary scattering

!        di = (dtau * rad) * exp(-tausum)
         if(dtau .gt. 0.)then
             ssa_eff = (alphabar_g + ssa * alphabar_a) / (alphabar_g + alphabar_a)
         else
             ssa_eff = 1.0
         endif

         alphabar = alphabar_g + alphabar_o + alphabar_a
         if(alphabar .gt. 0.)then
             di_g = do * radg       * (alphabar_g / alphabar)
             di_a = do * rada * ssa * (alphabar_a / alphabar)
         else
             di_g = 0.
             di_a = 0.
         endif
         di = di_g + di_a

         sumi_g = sumi_g + di_g
         sumi_a = sumi_a + di_a
         sumi = sumi + di

         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac_curr
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
           if(i .eq. 10 .or. i .eq. 100 .or. i .eq. 1000 .or. i .eq. 10000 .or. i .eq. iend)then
             write(6,11)i,sbar,htbar,alphabar_g*ds,alphabar_o*ds,tausum,od_solar_slant_a,od_solar_slant,rad,di,sumi_g,sumi_a,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i8,f11.0,f13.1,2f10.7,f8.4,f9.6,f9.4,2f9.5,6f9.4)
           endif
         endif
     enddo ! i

900  continue

     srcdir = sumi_extrap

     return
     end

     subroutine get_clr_src_dir_low(solalt,solazi,viewalt,viewazi & ! I
           ,od_g_msl,od_g_vert,od_o_msl,od_o_vert                 & ! I
           ,od_a_vert,ext_haero,htmsl                             & ! I
           ,ssa,agv,ao,aav,aod_ref,redp_lvl,scale_ht_a &            ! I
           ,ags_in,aas_in,ags_a,aas_a,isolalt_lo,isolalt_hi &       ! I
           ,ic,idebug,refdist_solalt,solalt_ref &                   ! I
           ,srcdir,sumi_g,sumi_a,opac_slant,nsteps,ds,tausum_a)

     use mem_namelist, ONLY: earth_radius
     use constants_laps, ONLY: R, GRAV

     include 'trigd.inc'
     include 'rad.inc'

     real*8 opac_last,opac_curr,tausum,trans,opac,od
     real logint

     trans(od) = exp(-min(od,80d0))
     opac(od) = 1d0 - trans(od)
     angdif(X,Y)=MOD(X-Y+540.,360.)-180.
     scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5
     logint(y1,y2,f) = exp(log(y1)*(1.-f) + log(y2)*f)
     curvat2(sdst,radius_start,altray) = &
             sqrt(sdst**2 + radius_start**2 &
               - (2.*sdst*radius_start*cosd(90.+altray))) - radius_start
     real*8 dsdst,dradius_start,daltray,dcurvat2
     dcurvat2(dsdst,dradius_start,daltray) = &
              sqrt(dsdst**2 + dradius_start**2 &
       - (2D0*dsdst*dradius_start*cosd(90D0+daltray))) - dradius_start      
     expbl(x) = exp(max(x,-80.))

     real tausum_a(nsteps)

     real agv      ! gas   relative to zenith value for observer ht/viewalt
     real ao       ! ozone relative to zenith at sea level pressure/viewalt
     real aav      ! aero  relative to zenith value for observer ht/viewalt
     real ah       ! haero relative to zenith value for observer ht/viewalt
     real ags_in   ! gas   relative to zenith value for observer ht/solalt
     real aas_in   ! aero  relative to zenith value for observer ht/solalt
     real ags_a(isolalt_lo:isolalt_hi) ! gas rel to zen val fr obs ht/salt
     real aas_a(isolalt_lo:isolalt_hi) ! aero rel to zen val fr obs ht/salt

!    return ! timing

!    od_o_msl = (o3_du/300.) * ext_o(ic)
     od_o_vert = od_o_msl * patm_o3(htmsl)

!    htmin_view = htminf(htmsl,viewalt,earth_radius)
!    patm_htmsl = ztopsa(htmin_view) / 1013.25
     patm_htmsl = ztopsa(htmsl) / 1013.25
     rearg = 2.0 + 0.667 * patm_htmsl
     radius_earth_eff = earth_radius * (1. + (0.33333 * patm_htmsl))
     radius_start_eff = radius_earth_eff + htmsl

!    High aerosols (just during low sun)
     if(solalt .lt. 333.)then
         od_h_vert = aod_ha * patm_ah(htmsl,h1_ha,h2_ha)
     else
         od_h_vert = 0.
     endif
     if(patm_ah(htmsl,h1_ha,h2_ha) .eq. 1.)then
         ah = am_homo_lyr(90.-viewalt,h1_ha-htmsl,h2_ha-htmsl,earth_radius)
     elseif(patm_ah(htmsl,h1_ha,h2_ha) .eq. 0.)then ! above layer
         ah = 0.
     else ! in layer
         ah = am_homo_lyr(90.-viewalt,0.,h2_ha-htmsl,earth_radius)
     endif

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * agv        
     od_o_slant = od_o_msl  * ao
     od_a_slant = od_a_vert * aav        
     od_h_slant = od_h_vert * ah
     od_slant = od_g_slant + od_o_slant + od_a_slant + od_h_slant
     od_slant = max(od_slant,1e-30) ! QC check
     r_e = 6371.e3

     scale_ht_g = R * ztotsa(htmsl) / GRAV
     alpha_sfc_g = od_g_vert / scale_ht_g

     scale_ht_g = R * ztotsa(0.) / GRAV
     alpha_ref_g = od_g_msl  / scale_ht_g

!    scale_ht_a = 1500. ! Input this?
!    alpha_sfc_a = od_a_vert / scale_ht_a
     alpha_ref_a = aod_ref / scale_ht_a

     opac_vert = opac(dble(od_g_vert + od_o_vert + od_a_vert))
     if(opac_vert .eq. 0)then ! very small od
         opac_vert = od_g_vert + od_o_vert + od_a_vert
     endif

     opac_slant = opac(dble(od_slant))
     if(opac_slant .eq. 0)then ! very small od
         opac_slant = od_slant
     endif

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
     ds = 100. / max(sind(viewalt),0.1)
     istart = max(1,nint((htmsl-100000.)/ds))

     dsolalt_ref = solalt_ref - solalt

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'subroutine get_clr_src_dir_low'
         write(6,*)'solalt/ref/solazi = ',solalt,solalt_ref,solazi
         write(6,*)'viewalt/viewazi = ',viewalt,viewazi
         write(6,*)'ags/aas = ',ags_in,aas_in
         write(6,*)'agv/aav = ',agv,aav
         write(6,*)'od_g_vert/od_g_msl = ',od_g_vert,od_g_msl
         write(6,*)'scale_ht_g = ',scale_ht_g
         write(6,*)'alpha_sfc_g/alpha_ref_g = ',alpha_sfc_g,alpha_ref_g
         write(6,*)'od_o_msl/vert = ',od_o_msl,od_o_vert
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'alpha_ha / ext_haero = ',alpha_ha,ext_haero
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_o_slant/ao = ',od_o_slant,ao
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'od_h_vert/slant = ',od_h_vert,od_h_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
         write(6,*)'istart/ds = ',istart,ds
     endif ! i

!    ags_a = ags_in; aas_a = aas_in ! test
     tausum = 0.
     sumi = 0. ; sumi_g = 0.; sumi_a = 0.
     frac_opac = 0.

     solaltill = max(solalt,0.)
     httopill = max(100000.-sind(solaltill)*60000.,40000.)
     htbotill = -(htmsl + 500.)

     if(idebug .eq. 1)then
       write(6,*)'           sbar    htbar_msl   dtau_g   dtau_a    tausum  dsolalt      ags     aas   od_sol_g  od_solar  rad        di      sumi_g    sumi_a opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     dsolalt_dxy = cosd(angdif(solazi,viewazi)) / 110000.

     if(solalt .ge. 0.)then
         frac_iso_a = 1.0
     else
         frac_iso_a = 0.5
     endif

     iend = max(istart+50000,nsteps)

     do i = istart,iend ! nsteps
         sbar = (float(i)-0.5) * ds

!        htbar = sbar * sind(viewalt)            ! Works at low altitude
!        htbar = htbar + sbar**2 / (rearg * r_e) ! Works at LEO   
!        htbar = curvat2(sbar,radius_start_eff,viewalt) ! Works at GEO
!        Works from moon
         htbar = dcurvat2(dble(sbar),dble(radius_start_eff),dble(viewalt)) 

         htbar_msl = htbar+htmsl
         xybar = (sbar - refdist_solalt) * cosd(viewalt)

         patm_bar = ztopsa(htbar_msl) / 1013.25
         scale_ht_g = R * ztotsa(htbar_msl) / GRAV
         alphabar_g = (od_g_msl * patm_bar) / scale_ht_g

         alphabar_o = alpha_o3(od_o_msl,htbar_msl)
!        alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         alphabar_a = alpha_ref_a * expbl(-(htbar_msl-redp_lvl)/scale_ht_a) 

!        AOD stratospheric = .012 (nighttime - layer extends above observer)
         if(htbar_msl .gt. h1_ha .and. htbar_msl .le. h2_ha .and. solalt .lt. 333. .and. aod_ha .gt. 0.)then
             alphabar_a = alphabar_a + alpha_ha * ext_haero
         elseif(htbar_msl .gt. h2_ha .and. htbar_msl .lt. h3_ha)then
             frac_ha = (h3_ha - htbar_msl) / (h3_ha - h2_ha)
             alphabar_a = alphabar_a + alpha_ha * ext_haero * scurve(frac_ha)
         endif

         dtau = ds * (alphabar_g + alphabar_o + alphabar_a)
         tausum = tausum + dtau
!        if(tausum .lt. 1.0)distod1 = sbar
!        tausum_a(i) = tausum

         if(htbar .gt. httopill .or. htbar .lt. htbotill .or. frac_opac .gt. 11.4 .or. (sbar.gt.refdist_solalt .and. htbar_msl.gt.130e3))then ! efficiency
!          tausum_a(min(i+1,nsteps):nsteps) = tausum
           goto 900
         endif

         opac_last = opac_curr
         opac_curr = opac(tausum)
         do = opac_curr - opac_last

!        Update ags,aas
         dsolalt = dsolalt_ref + dsolalt_dxy * xybar
         solalt_step = solalt + dsolalt
         dsolaltb = &
            max(min(dsolalt,float(isolalt_hi)-.0001),float(isolalt_lo))
         isolaltl = floor(dsolaltb); isolalth = isolaltl + 1
         fsolalt = dsolaltb - float(isolaltl)
!        ags = ags_a(isolaltl) * (1.-fsolalt) + ags_a(isolalth) * fsolalt
!        aas = aas_a(isolaltl) * (1.-fsolalt) + aas_a(isolalth) * fsolalt
         ags = logint(ags_a(isolaltl),ags_a(isolalth),fsolalt)
         aas = logint(aas_a(isolaltl),aas_a(isolalth),fsolalt)

!        Considering smoothing out with solar disk
         horz_dep = horz_depf(htbar_msl,earth_radius)
         if(solalt_step + 0.5 .lt. -horz_dep)then 
            sol_occ = 0.0 ! invisible
         else
            sol_occ = 1.0 ! visible
         endif

         if(solalt_step .ge. 0.)then
           zapp = 90. - solalt_step
           od_solar_slant_o3 = od_o_msl * patm_o3(htbar_msl) * airmasso(zapp,htbar_msl)
         elseif(sol_occ .gt. 0.)then
           zappi = 90. + solalt_step
           patm = exp(-(htbar_msl/scale_ht_g))
           ztruei = zappi + refractd_app(abs(solalt_step) ,patm)
           patmo  = patm_o3(htbar_msl)
           htmin_ray = htminf(htbar_msl,solalt_step,earth_radius)
           patm2o = patm_o3(htmin_ray)
           ztrue0 = 90. + refractd_app(0.,patm2o)
!          alttrue0 = -refractd_app(0.,patm2o)
!          ao1 = airmasso(ztruei,htmin_ray) * patm2o ! ztruei from htmin
           ao2 = airmasso(ztruei,htbar_msl) * patmo  ! upward from observer
           ao3 = airmasso(ztrue0,htmin_ray) * patm2o ! 0 deg from htmin
           ao_calc = (ao3-ao2) + ao3
           od_solar_slant_o3 = od_o_msl * ao_calc * sol_occ
         else
           od_solar_slant_o3 = 999.
         endif

         if(sol_occ .gt. 0.)then

           if(htmsl .le. 150e3 .and. .false.)then
!                        M         L              S
             od_solar_slant = &
                    od_g_vert * ags * exp(-htbar/scale_ht_g)  &
                  + od_solar_slant_o3 &
                  + od_a_vert * aas * exp(-htbar/scale_ht_a) * frac_iso_a
           else ! presently testing
             od_solar_slant_g = od_g_msl * ags * exp(-htbar_msl/scale_ht_g)
             od_solar_slant_a = aod_ref  * aas * exp(-(htbar_msl-redp_lvl)/scale_ht_a) * frac_iso_a
             od_solar_slant = od_solar_slant_g &
                            + od_solar_slant_o3 + od_solar_slant_a
           endif

           rad = trans(dble(od_solar_slant))
           radg = max(rad,.0001) ! secondary scattering
           rada = max(rad,.00001) ! secondary scattering
           rad = max(rad,.0001) ! secondary scattering

!          di = (dtau * rad) * exp(-tausum)
           if(dtau .gt. 0.)then
             ssa_eff = (alphabar_g + ssa * alphabar_a) / (alphabar_g + alphabar_a)
           else
             ssa_eff = 1.0
           endif

           alphabar = alphabar_g + alphabar_o + alphabar_a
           if(alphabar .gt. 0.)then
             di_g = do * radg       * (alphabar_g / alphabar)
             di_a = do * rada * ssa * (alphabar_a / alphabar)
           else
             di_g = 0.
             di_a = 0.
           endif
           di = di_g + di_a

         else ! sol_occ = 0
           od_solar_slant = 999.
           rad = 0. ; di_g = 0. ; di_a = 0. ; di = 0.

         endif ! sol_occ

         sumi_g = sumi_g + di_g
         sumi_a = sumi_a + di_a
         sumi = sumi + di

         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac(tausum)
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
!          if(i .le. 10 .or. ((i .eq. (i/25)*25) .and. i .le. 400) .or. i .eq. 500 .or. i .eq. 1000 .or. i .eq. 1500 .or. i .eq. 2000. .or. i .eq. 2500 .or. i .eq. 10000)then
           if((i .le. 10 .or. i .eq. (i/25)*25) .and. htbar_msl .gt. 0. .and. htbar_msl .le. 110e3)then
             if(sol_occ.gt.0.)then
               write(6,12)horz_dep,solalt_step,htmin_ray,sol_occ ! ,ao2,ao3,ao
             endif
12           format(33x,'   horz_dep/solalt/htmin/sol_occ = ',2f8.2,f9.0,f8.3,3f9.4)
             write(6,11)i,sbar,htbar_msl,alphabar_g*ds,alphabar_a*ds,tausum,dsolalt,ags,aas,od_solar_slant_g,od_solar_slant,rad,di,sumi_g,sumi_a,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i6,f11.0,f12.1,2e10.3,f8.4,f9.4,f10.2,4f9.4,e9.2,2f11.8,2f9.4,2f9.4)
           endif
         endif
     enddo ! i

900  continue

     if(idebug .eq. 1)then
       write(6,11)i,sbar,htbar_msl,alphabar_g*ds,alphabar_a*ds,tausum,dsolalt,ags,aas,od_solar_slant_g,od_solar_slant,rad,di,sumi_g,sumi_a,opac_curr,frac_opac,sumi_mean,sumi_extrap
       if(sol_occ.gt.0.)write(6,12)horz_dep,solalt_step,htmin_ray,sol_occ ! ,ao2,ao3,ao
     endif

     srcdir = sumi_extrap

     return
     end


     subroutine get_clr_src_dir_topo(solalt,solazi,viewalt,viewazi &   ! I
           ,od_g_msl,od_g_vert,od_o_msl,od_a_vert,htmsl,dist_to_topo & ! I
           ,ssa,agv,aav,aod_ref,redp_lvl,scale_ht_a &                  ! I
           ,ags_a,aas_a,isolalt_lo,isolalt_hi &                        ! I
           ,ic,idebug,nsteps,refdist_solalt,solalt_ref &               ! I
           ,sumi_g,sumi_a,opac_slant,tausum_a)                         ! O

     use mem_namelist, ONLY: earth_radius
     use constants_laps, ONLY: R, GRAV
     include 'trigd.inc'
     include 'rad.inc'

     trans(od) = exp(-min(od,80.))
     opac(od) = 1.0 - trans(od)
     angdif(X,Y)=MOD(X-Y+540.,360.)-180.
     curvat(hdst,radius) = hdst**2 / (2. * radius)
     curvat2(sdst,radius_start,altray) = &
             sqrt(sdst**2 + radius_start**2 &
               - (2.*sdst*radius_start*cosd(90.+altray))) - radius_start
     real*8 dsdst,dradius_start,daltray,dcurvat2
     real*8 send,sbar,htbar,dist_to_topo,ds
     dcurvat2(dsdst,dradius_start,daltray) = &
              sqrt(dsdst**2 + dradius_start**2 &
       - (2D0*dsdst*dradius_start*cosd(90D0+daltray))) - dradius_start      
     expbl(x) = exp(max(x,-80.))

     real tausum_a(nsteps)
     real ags_a(isolalt_lo:isolalt_hi) ! gas rel to zen val fr obs ht/salt
     real aas_a(isolalt_lo:isolalt_hi) ! aero rel to zen val fr obs ht/salt

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * agv        
     od_a_slant = od_a_vert * aav        
     od_slant = od_g_slant + od_a_slant
     r_e = 6371.e3
     patm_htmsl = ztopsa(htmsl) / 1013.25
     radius_earth_eff = earth_radius * (1. + (0.33333 * patm_htmsl))
     radius_start_eff = radius_earth_eff + htmsl

     scale_ht_g = R * ztotsa(htmsl) / GRAV
!    alpha_sfc_g = od_g_vert / scale_ht_g

!    scale_ht_g = 8000.
     scale_ht_g = R * ztotsa(0.) / GRAV
     alpha_ref_g = od_g_msl  / scale_ht_g

!    alpha_sfc_a = od_a_vert / scale_ht_a
     alpha_ref_a = aod_ref / scale_ht_a

     opac_vert = opac(od_g_vert + od_a_vert)
     if(opac_vert .eq. 0)then ! very small od
         opac_vert = od_g_vert + od_a_vert
     endif
     opac_slant = opac(od_g_slant + od_a_slant)
     if(opac_slant .eq. 0)then ! very small od
         opac_slant = od_g_slant + od_a_slant
     endif

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
     ds = min(dist_to_topo/10d0,1000d0)           
     if(dist_to_topo .le. 10000d0)then
       nsteps_topo = nint(dist_to_topo/ds) ! max(10,int(htmsl/1000.))
     else
       nsteps_topo = int(dist_to_topo/ds) + 1
     endif
     tausum = 0.
     sumi = 0.; sumi_g = 0.; sumi_a = 0.
     httopill = max(100000.-sind(solalt)*60000.,40000.)
     htbotill = -(htmsl + 500.)
     httopo = 0.
     istart = max(1,nint((htmsl-100000d0)/ds))

     dsolalt_ref = solalt_ref - solalt

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'subroutine get_clr_src_dir_topo'
         write(6,*)'solalt/ref/solazi = ',solalt,solalt_ref,solazi
         write(6,*)'viewalt/viewazi = ',viewalt,viewazi
         write(6,*)'agv/aav = ',agv,aav
         write(6,*)'od_g_vert = ',od_g_vert
         write(6,*)'scale_ht_g = ',scale_ht_g
         write(6,*)'alpha_sfc_g/alpha_ref_g = ',alpha_sfc_g,alpha_ref_g
         write(6,*)'od_o_msl = ',od_o_msl
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'aod_ref/redp_lvl/scale_ht_a',aod_ref,redp_lvl,scale_ht_a
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
         write(6,*)'dist_to_topo = ',dist_to_topo
         write(6,*)'istart/nsteps_topo/ds = ',istart,nsteps_topo,ds
         write(6,*)'htbotill/httopill',htbotill,httopill
     endif ! i

     if(idebug .eq. 1)then ! substitute radg/rada?
       write(6,*)'         sbar      htbar      htmsl     dtau    tausum  dsolalt alphbr_g alphbr_a  od_solar   rad       di     sumi_g   sumi_a  opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     dsolalt_dxy = cosd(angdif(solazi,viewazi)) / 110000.

     iwrite_slant = 0

     do i = istart,nsteps_topo
       sbar = (dble(i)-0.5d0) * ds
       send =  dble(i)        * ds
!      htbar = sbar * sind(viewalt)
!      xybar = sbar * cosd(viewalt)
!      htbar = htbar + curvat(xybar,earth_radius*(4./3.)) ! Earth Curvature
!      htbar = curvat2(sbar,radius_start_eff,viewalt) ! From GEO
       htbar = dcurvat2(sbar,dble(radius_start_eff),dble(viewalt)) !Fm Moon
       htbar_msl_last = htbar_msl
       htbar_msl = htbar+htmsl
       ht_msl_last = ht_msl
       ht_msl = 1.5 * htbar_msl - 0.5 * htbar_msl_last
       if(htbar_msl .lt. httopill)then
         if(send .lt. dist_to_topo)then
           frac_step = 1.0
         else
           frac_step = 1d0 - (send - dist_to_topo)/ds
         endif
         xybar = (sbar - refdist_solalt) * cosd(viewalt)

!        Update ags,aas
         dsolalt = dsolalt_ref + dsolalt_dxy * xybar
         solalt_step = solalt + dsolalt
         dsolaltb = &
            max(min(dsolalt,float(isolalt_hi)-.0001),float(isolalt_lo))
         isolaltl = floor(dsolaltb); isolalth = isolaltl + 1
         fsolalt = dsolaltb - float(isolaltl)
         ags = ags_a(isolaltl) * (1.-fsolalt) + ags_a(isolalth) * fsolalt
         aas = aas_a(isolaltl) * (1.-fsolalt) + aas_a(isolalth) * fsolalt

!        Considering smoothing out with solar disk
         horz_dep = horz_depf(htbar_msl,earth_radius)
         if(solalt_step + 0.5 .lt. -horz_dep)then 
            sol_occ = 0.0 ! invisible
         else
            sol_occ = 1.0 ! visible
         endif

         if(solalt_step .ge. 0.)then
           zapp = 90. - solalt_step
           od_solar_slant_o3 = od_o_msl * patm_o3(htbar_msl) * airmasso(zapp,htbar_msl)
         elseif(sol_occ .gt. 0.)then
           zappi = 90. + solalt_step
           patm = exp(-(htbar_msl/scale_ht_g))
           ztruei = zappi + refractd_app(abs(solalt_step) ,patm)
           patmo  = patm_o3(htbar_msl)
           htmin_ray = htminf(htbar_msl,solalt_step,earth_radius)
           patm2o = patm_o3(htmin_ray)
           ztrue0 = 90. + refractd_app(0.,patm2o)
           ao2 = airmasso(ztruei,htbar_msl) * patmo  ! upward from observer
           ao3 = airmasso(ztrue0,htmin_ray) * patm2o ! 0 deg from htmin
           ao_calc = (ao3-ao2) + ao3
           od_solar_slant_o3 = od_o_msl * ao_calc * sol_occ
         else ! more accurate calculation technique was imported
           od_solar_slant_o3 = 0.
         endif

         if(i .eq. nsteps_topo .and. idebug .eq. 1)then
             write(6,*)'aas terms',isolaltl,aas_a(isolaltl),isolalth,aas_a(isolalth),fsolalt
             write(6,*)'od_solar_slant terms',od_g_vert,ags,htbar,scale_ht_g,exp(-htbar/scale_ht_g),aod_ref,aas,htbar,scale_ht_a,expbl(-(htbar_msl-redp_lvl)/scale_ht_a)
             iwrite_slant = 1
         endif

         if(sol_occ .gt. 0.)then
           od_solar_slant = od_g_msl * ags * expbl(-htbar_msl/scale_ht_g) &
                    + od_solar_slant_o3 &
                    + aod_ref * aas * expbl(-(htbar_msl-redp_lvl)/scale_ht_a)
         else
           od_solar_slant = 999.
         endif

!        alphabar_g = alpha_sfc_g * exp(-htbar/scale_ht_g) 
!        alphabar_g = alpha_ref_g * exp(-htbar_msl/scale_ht_g) 
         patm_bar = ztopsa(htbar_msl) / 1013.25
         scale_ht_g = R * ztotsa(htbar_msl) / GRAV
         alphabar_g = (od_g_msl * patm_bar) / scale_ht_g

!        alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         alphabar_a = alpha_ref_a * expbl(-(htbar_msl-redp_lvl)/scale_ht_a)
         dtau = ds * (alphabar_g + alphabar_a) * frac_step
         tausum = tausum + dtau
!        if(tausum .lt. 1.0)distod1 = sbar
!        tausum_a(i) = tausum

         opac_last = opac_curr
         opac_curr = opac(tausum)
         do = opac_curr - opac_last

!        rad = 1.0
         rad = trans(od_solar_slant)
         radg = rad ! max(rad,.0001) ! secondary scattering
         rada = rad ! max(rad,.00001) ! secondary scattering
!        rad = max(rad,.0001) ! secondary scattering

!        di = (dtau * rad) * exp(-tausum)
         if(dtau .gt. 0.)then
             ssa_eff = (alphabar_g + ssa * alphabar_a) / (alphabar_g + alphabar_a)
         else
             ssa_eff = 1.0
         endif

         di   = do * rad  * ssa_eff
         di_g = do * radg       * alphabar_g / (alphabar_g + alphabar_a)
         di_a = do * rada * ssa * alphabar_a / (alphabar_g + alphabar_a)

         sumi   = sumi   + di
         sumi_g = sumi_g + di_g
         sumi_a = sumi_a + di_a

         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac(tausum)
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
!          if(i .le. 10)then
             write(6,11)i,sbar,htbar,htbar_msl,dtau,tausum,dsolalt,alphabar_g,alphabar_a,od_solar_slant,rad,di,sumi_g,sumi_a,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i6,f11.0,f11.0,f9.0,f9.4,f9.6,f9.4,2f9.6,9f9.4)
!          endif
         endif
       endif ! htbar < httopill

       if(htbar_msl .lt. httopo .or. (htbar_msl .gt. httopill .and. viewalt .gt. 0.))then
         if(idebug .eq. 1)then
           write(6,*)' exit loop:',htbar_msl,ht_msl,httopo,frac_step
         endif
         goto 900
       endif
     enddo ! i

     if(idebug .eq. 1)then
       write(6,*)' end of loop:',htbar_msl,ht_msl,httopo,frac_step
     endif

900  continue

     return
     end

