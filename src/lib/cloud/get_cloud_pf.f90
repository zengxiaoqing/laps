
        subroutine get_cld_pf(elong_a,r_cloud_rad,cloud_od,cloud_od_sp,airmass_2_topo,ni,nj & ! I
                             ,pf_scat1,pf_scat2,pf_scat,bkscat_alb) ! O
        include 'trigd.inc'

!       Statement functions
        trans(od) = exp(-od)
        opac(od) = 1.0 - exp(-od)

        include 'rad.inc'
        real elong_a(ni,nj)
        real cloud_od(ni,nj)        ! cloud optical depth
        real cloud_od_sp(ni,nj,nsp) ! cloud species optical depth
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        real pf_scat(ni,nj), pf_scat1(ni,nj), pf_scat2(ni,nj)
        real bkscat_alb(ni,nj)

        do j = 1,nj
        do i = 1,ni

!         a substitute for cloud_rad could be arg2 = (cosd(alt))**3.

!         Phase function that depends on degree of fwd scattering in cloud    
!         pwr controls angular extent of fwd scattering peak of a thin cloud
!         General references: 
!         http://www-evasion.imag.fr/Publications/2008/BNMBC08/clouds.pdf
!         http://rtweb.aer.com/rrtm_frame.html
!         http://wiki.seas.harvard.edu/geos-chem/images/Guide_to_GCRT_Oct2013.pdf
          if(elong_a(i,j) .le. 90.)then ! low elongation phase function
              pwr = 3.0
              ampl = r_cloud_rad(i,j) * 0.7; b = 1.0 + pwr * r_cloud_rad(i,j)
              pf_scat1(i,j) = 0.9 + ampl * (cosd(min(elong_a(i,j),89.99))**b)
              bkscat_alb(i,j) = -99.9 ! flag value for logging
          else                          ! high elongation phase function
!             convert from opacity to albedo
!             cloud_opacity = min(r_cloud_3d(i,j),0.999999)
              bksc_eff_od = cloud_od(i,j) * 0.10 
              cloud_rad_trans = exp(-bksc_eff_od)
              bkscat_alb(i,j) = 1.0 - cloud_rad_trans 
              ampl = 0.15 * bkscat_alb(i,j)
              pf_scat1(i,j) = 0.9 + ampl * (-cosd(elong_a(i,j)))
          endif

!         Separate snow phase function
!         opac_nonsnow  =  opac(cloud_od(i,j) - cloud_od_sp(i,j,4))
          trans_nonsnow = trans(cloud_od(i,j) - cloud_od_sp(i,j,4))
          cloud_od_snow = cloud_od_sp(i,j,4)
          cloud_od_tot  = cloud_od(i,j)

          snow_bin1 = exp(-cloud_od_snow/5.)
          snow_bin2 = 1.0 - snow_bin1

          pf_snow = snow_bin1 * hg(.95,elong_a(i,j)) &
                  + snow_bin2 * hg(0.0,elong_a(i,j))

          if(cloud_od_tot .gt. 0.)then
              snow_factor = trans_nonsnow * (cloud_od_snow / cloud_od_tot)
          else
              snow_factor = 0.
          endif

          pf_scat2(i,j) = pf_snow * snow_factor + pf_scat1(i,j) * (1.0 - snow_factor)

!         Suppress phase function if terrain is close in the light ray
          if(airmass_2_topo(i,j) .gt. 0.)then ! cloud in front of terrain
!             pf_scat(i,j) = pf_scat2(i,j)**r_cloud_rad(i,j)
              pf_scat(i,j) = pf_scat2(i,j)**opac(cloud_od_tot)
          else
              pf_scat(i,j) = pf_scat2(i,j)
          endif

        enddo ! j
        enddo ! i

        return
        end

