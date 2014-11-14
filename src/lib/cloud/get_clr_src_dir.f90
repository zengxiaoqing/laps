
     subroutine get_clr_src_dir(solalt,viewalt,od_g_vert,od_a_vert,ag,aa,ags,aas,idebug,srcdir,opac_slant,nsteps,ds,tausum_a)

     include 'trigd.inc'

!    Input idebug for select angles

     trans(od) = exp(-min(od,80.))
     opac(od) = 1.0 - trans(od)

     real tausum_a(nsteps)

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * ag        
     od_a_slant = od_a_vert * aa        
     od_slant = od_g_slant + od_a_slant
     r_e = 6371.e3

     scale_ht_g = 8000.
     alpha_sfc_g = od_g_vert / scale_ht_g

     scale_ht_a = 1500. ! Input this?
     alpha_sfc_a = od_a_vert / scale_ht_a

     opac_vert = opac(od_g_vert + od_a_vert)
     opac_slant = opac(od_g_slant + od_a_slant)

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'solalt/viewalt = ',solalt,viewalt
         write(6,*)'ags/aas = ',ags,aas
         write(6,*)'agv/aav = ',ag,aa
         write(6,*)'od_g_vert = ',od_g_vert
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
     endif ! i

     delta_slant = opac_slant / float(nsteps)

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
     ds = 25.
!    nsteps = 1000
!    nsteps = 10000
     tausum = 0.
     sumi = 0.

     if(idebug .eq. 1)then
       write(6,*)'         sbar     htbar    dtau    tausum  od_solar   rad       di      sumi  opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     do i = 1,nsteps
         sbar = (float(i)-0.5) * ds
         htbar = sbar * sind(viewalt)
         htbar = htbar + sbar**2 / (2.667 * r_e) ! Earth Curvature

         od_solar_slant = od_g_vert * ags * exp(-htbar/scale_ht_g) + od_a_vert * aas * exp(-htbar/scale_ht_a)

         alphabar_g = alpha_sfc_g * exp(-htbar/scale_ht_g) 
         alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         dtau = ds * (alphabar_g + alphabar_a)
         tausum = tausum + dtau
!        if(tausum .lt. 1.0)distod1 = sbar
         tausum_a(i) = tausum

         opac_last = opac_curr
         opac_curr = opac(tausum)
         do = opac_curr - opac_last

!        rad = 1.0
         rad = trans(od_solar_slant)
         rad = max(rad,.03) ! secondary scattering

!        di = (dtau * rad) * exp(-tausum)
         di = do * rad
         sumi = sumi + di
         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac(tausum)
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
           if(i .eq. 10 .or. i .eq. 100 .or. i .eq. 1000 .or. i .eq. 10000 .or. i .eq. nsteps)then
             write(6,11)i,sbar,htbar,dtau,tausum,od_solar_slant,rad,di,sumi,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i6,f9.0,f9.1,10f9.4)
           endif
         endif
     enddo ! i

     srcdir = sumi_extrap

     return
     end

     subroutine get_clr_src_dir_low(solalt,solazi,viewalt,viewazi &
                   ,od_g_vert,od_a_vert,ag,aa,ags_in,aas_in,ags_a,aas_a &
                   ,idebug,srcdir,opac_slant,nsteps,ds,tausum_a)

     include 'trigd.inc'

     trans(od) = exp(-min(od,80.))
     opac(od) = 1.0 - trans(od)
     angdif(X,Y)=MOD(X-Y+540.,360.)-180.

     real tausum_a(nsteps)

     real ags_a(-5:+5), aas_a(-5:+5)

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * ag        
     od_a_slant = od_a_vert * aa        
     od_slant = od_g_slant + od_a_slant
     r_e = 6371.e3

     scale_ht_g = 8000.
     alpha_sfc_g = od_g_vert / scale_ht_g

     scale_ht_a = 1500. ! Input this?
     alpha_sfc_a = od_a_vert / scale_ht_a

     opac_vert = opac(od_g_vert + od_a_vert)
     opac_slant = opac(od_g_slant + od_a_slant)

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'subroutine get_clr_src_dir_low'
         write(6,*)'solalt/solazi = ',solalt,solazi
         write(6,*)'viewalt/viewazi = ',viewalt,viewazi
         write(6,*)'ags/aas = ',ags_in,aas_in
         write(6,*)'agv/aav = ',ag,aa
         write(6,*)'od_g_vert = ',od_g_vert
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
     endif ! i

!    ags_a = ags_in; aas_a = aas_in ! test

     delta_slant = opac_slant / float(nsteps)

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
     ds = 50.
     tausum = 0.
     sumi = 0.
     httopill = max(130000.-sind(solalt)*90000.,40000.)

     if(idebug .eq. 1)then
       write(6,*)'         sbar     htbar    dtau    tausum  dsolalt     ags      aas   od_solar   rad       di      sumi  opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     dsolalt_dxy = cosd(angdif(solazi,viewazi)) / 110000.

     do i = 1,nsteps
         sbar = (float(i)-0.5) * ds
         htbar = sbar * sind(viewalt)
         htbar = htbar + sbar**2 / (2.667 * r_e) ! Earth Curvature
         xybar = sbar * cosd(viewalt)

!        Update ags,aas
         dsolalt = dsolalt_dxy * xybar
         dsolaltb = max(min(dsolalt,4.9999),-5.0)
         isolaltl = floor(dsolaltb); isolalth = isolaltl + 1
         fsolalt = dsolaltb - float(isolaltl)
         ags = ags_a(isolaltl) * (1.-fsolalt) + ags_a(isolalth) * fsolalt
         aas = aas_a(isolaltl) * (1.-fsolalt) + aas_a(isolalth) * fsolalt

         od_solar_slant = od_g_vert * ags * exp(-htbar/scale_ht_g) + od_a_vert * aas * exp(-htbar/scale_ht_a)

         alphabar_g = alpha_sfc_g * exp(-htbar/scale_ht_g) 
         alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         dtau = ds * (alphabar_g + alphabar_a)
         tausum = tausum + dtau
!        if(tausum .lt. 1.0)distod1 = sbar
!        tausum_a(i) = tausum

         if(htbar .gt. httopill)then ! efficiency
!            tausum_a(min(i+1,nsteps):nsteps) = tausum
             goto900
         endif

         opac_last = opac_curr
         opac_curr = opac(tausum)
         do = opac_curr - opac_last

!        rad = 1.0
         rad = trans(od_solar_slant)
         rad = max(rad,.03) ! secondary scattering

!        di = (dtau * rad) * exp(-tausum)
         di = do * rad
         sumi = sumi + di
         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac(tausum)
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
           if(i .eq. 10 .or. i .eq. 100 .or. i .eq. 1000 .or. i .eq. 10000)then
             write(6,11)i,sbar,htbar,dtau,tausum,dsolalt,ags,aas,od_solar_slant,rad,di,sumi,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i6,f9.0,f9.1,13f9.4)
           endif
         endif
     enddo ! i

900  continue
     write(6,11)i,sbar,htbar,dtau,tausum,dsolalt,ags,aas,od_solar_slant,rad,di,sumi,opac_curr,frac_opac,sumi_mean,sumi_extrap

     srcdir = sumi_extrap

     return
     end

