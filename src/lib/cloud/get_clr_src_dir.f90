
     subroutine get_clr_src_dir(solalt,viewalt,od_g_vert,od_a_vert,ag,aa,idebug,srcdir)

!    Input idebug for select angles

     trans(od) = exp(-min(od,80.))
     opac(od) = 1.0 - trans(od)

     nsteps = 2

!    Calculate src term of gas+aerosol along ray normalized by TOA radiance
     od_g_slant = od_g_vert * ag        
     od_a_slant = od_a_vert * aa        
     od_slant = od_g_slant + od_a_slant
     r_e = 6371.e3

     scale_ht_g = 8000.
     alpha_sfc_g = od_g_vert / scale_ht_g

     scale_ht_a = 2000.
     alpha_sfc_a = od_a_vert / scale_ht_a

     opac_vert = opac(od_g_vert + od_a_vert)
     opac_slant = opac(od_g_slant + od_a_slant)

     if(idebug .eq. 1)then
         write(6,*)
         write(6,*)'ag/aa = ',ag,aa
         write(6,*)'od_g_vert = ',od_g_vert
         write(6,*)'od_a_vert = ',od_a_vert
         write(6,*)'opac_vert = ',opac_vert
         write(6,*)'od_g_slant = ',od_g_slant
         write(6,*)'od_a_slant = ',od_a_slant
         write(6,*)'opac_slant = ',opac_slant
         write(6,*)'od_slant = ',od_slant
     endif ! i

     delta_slant = opac_slant / float(nsteps)

     do i = 1,nsteps
         opac_slant_step = (float(i) - 0.5) * delta_slant
         od_slant_step = -log(1.0 - opac_slant_step)
         if(idebug .eq. 1)then
           write(6,*)'i/opac_slant_step/od_slant_step=',i,opac_slant_step,od_slant_step
         endif
     enddo ! i

!    Integral [0 to x] a * exp(-b*x) = a * (1. - exp*(-b*x)) / b
     dx = 25.
     nsteps = 1000
     nsteps = 10000
     tausum = 0.
     sumi = 0.

     if(idebug .eq. 1)then
       write(6,*)'         xbar     htbar    dtau    tausum  od_solar   rad       di      sumi  opac_curr frac_opac sumi_mn sumi_ext'
     endif

     opac_curr = 0.

     do i = 1,nsteps
         xbar = (float(i)-0.5) * dx
         htbar = xbar * sind(viewalt)
         htbar = htbar + xbar**2 / (2.667 * r_e)

         od_solar_vert = od_g_vert * exp(-htbar/scale_ht_g) + od_a_vert * exp(-htbar/scale_ht_a)
         od_solar_slant = od_solar_vert / sind(solalt)

         alphabar_g = alpha_sfc_g * exp(-htbar/scale_ht_g) 
         alphabar_a = alpha_sfc_a * exp(-htbar/scale_ht_a)  
         dtau = dx * (alphabar_g + alphabar_a)
         tausum = tausum + dtau

         opac_last = opac_curr
         opac_curr = opac(tausum)
         do = opac_curr - opac_last

!        rad = 1.0
         rad = trans(od_solar_slant)
         rad = max(rad,.15)

!        di = (dtau * rad) * exp(-tausum)
         di = do * rad
         sumi = sumi + di
         frac_opac = opac_curr / opac_slant

         sumi_mean = sumi / opac(tausum)
!        sumi_extrap = (sumi + (1.-frac_opac) * rad) / opac_slant
         sumi_extrap = (sumi_mean * frac_opac) + (1.-frac_opac) * rad

         if(idebug .eq. 1)then
           if(i .eq. 10 .or. i .eq. 100 .or. i .eq. 1000 .or. i .eq. 10000 .or. i .eq. nsteps)then
             write(6,11)i,xbar,htbar,dtau,tausum,od_solar_slant,rad,di,sumi,opac_curr,frac_opac,sumi_mean,sumi_extrap
11           format(i6,f9.0,f9.1,10f9.4)
           endif
         endif
     enddo ! i

     srcdir = sumi_extrap

     return
     end

