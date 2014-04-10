
        subroutine get_cld_pf(elong_a,r_cloud_rad,cloud_od,cloud_od_sp,nsp,airmass_2_topo,ni,nj & ! I
                             ,pf_scat1,pf_scat2,pf_scat,bkscat_alb) ! O
        include 'trigd.inc'

!       Statement functions
        trans(od) = exp(-od)
        opac(od) = 1.0 - exp(-od)
        counts_to_rad(counts) = 10.**((counts/100.)+7.3)
        pf_bk(elong,alb) = 1.0 + cosd(180. - elong)

        include 'rad.inc'
        real elong_a(ni,nj)
        real cloud_od(ni,nj)        ! cloud optical depth (tau)
        real cloud_od_sp(ni,nj,nsp) ! cloud species tau (clwc,cice,rain,snow)
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        real pf_scat(ni,nj), pf_scat1(ni,nj), pf_scat2(ni,nj)
        real bkscat_alb(ni,nj)

        logical l_pf_new /.true./

        do j = 1,nj
        do i = 1,ni

!         a substitute for cloud_rad could be arg2 = (cosd(alt))**3.

!         Phase function that depends on degree of fwd scattering in cloud    
!         pwr controls angular extent of fwd scattering peak of a thin cloud
!         General references: 
!         http://www.uni-leipzig.de/~strahlen/web/research/de_index.php?goto=arctic
!         http://www-evasion.imag.fr/Publications/2008/BNMBC08/clouds.pdf
!         http://rtweb.aer.com/rrtm_frame.html
!         http://wiki.seas.harvard.edu/geos-chem/images/Guide_to_GCRT_Oct2013.pdf

          cloud_od_cice = cloud_od_sp(i,j,2)

          if(l_pf_new)then ! new liquid phase function
            cloud_od_clwc = cloud_od_sp(i,j,1) + cloud_od_sp(i,j,3)

            hgp = max(cloud_od_clwc,1.0)       ! multiple scattering phase func
            sco = hgp

!           Check for being optically thin with direct light source
            clwc_bin1 = exp(-(cloud_od_clwc/10.)**2) ! optically thin
            clwc_bin1 = clwc_bin1 * r_cloud_rad(i,j)**10. ! direct lighting
            clwc_bin1a = 0.80**sco             ! forward peak
            clwc_bin1c = opac(0.00*sco)        ! backscattering
            clwc_bin1b = 1.0 - (clwc_bin1a + clwc_bin1c) ! mid    
            clwc_bin2 = 1.0 - clwc_bin1        ! optically thick 

            pf_clwc = clwc_bin1a * clwc_bin1 * hg(.94**hgp,elong_a(i,j)) & ! corona
                    + clwc_bin1b * clwc_bin1 * hg(.60**hgp,elong_a(i,j)) &
                    + clwc_bin1c * clwc_bin1 * hg(.00     ,elong_a(i,j)) &
                    + clwc_bin2  * hg(-0.3    ,elong_a(i,j))  

            pf_scat1(i,j) = pf_clwc

          else ! old liquid phase function
            if(elong_a(i,j) .le. 90.)then ! low elongation phase function
              pwr = 3.0
              ampl = r_cloud_rad(i,j) * 0.7; b = 1.0 + pwr * r_cloud_rad(i,j)
              pf_scat1(i,j) = 0.9 + ampl * (cosd(min(elong_a(i,j),89.99))**b)
              bkscat_alb(i,j) = -99.9 ! flag value for logging
            else                          ! high elongation phase function
              bksc_eff_od = cloud_od(i,j) * 0.10 
              cloud_rad_trans = exp(-bksc_eff_od)
              bkscat_alb(i,j) = 1.0 - cloud_rad_trans 
              ampl = 0.15 * bkscat_alb(i,j)
              pf_scat1(i,j) = 0.9 + ampl * (-cosd(elong_a(i,j)))
            endif

            pf_scat1(i,j) = counts_to_rad(240. * pf_scat1(i,j)) / counts_to_rad(240.)
          endif

!         Separate snow (+cice) phase function
          cloud_od_snow = cloud_od_sp(i,j,4) + cloud_od_sp(i,j,2)
          trans_nonsnow = trans(cloud_od(i,j) - cloud_od_snow)
          cloud_od_tot  = cloud_od(i,j)

          hgp = max(cloud_od_snow,1.0)       ! multiple scattering phase func
          sco = hgp

!         snow_bin1 = exp(-cloud_od_snow/5.) ! optically thin snow
          snow_bin1a = 0.60**sco             ! forward peak
          snow_bin1c = opac(0.10*sco)        ! backscattering
          snow_bin1b = 1.0 - (snow_bin1a + snow_bin1c) ! mid
!         snow_bin2 = 1.0 - snow_bin1        ! optically thick snow

          pf_snow = snow_bin1a * hg(.98**hgp,elong_a(i,j)) & ! plates
                  + snow_bin1b * hg(.50**hgp,elong_a(i,j)) &
                  + snow_bin1c * hg(.00     ,elong_a(i,j))  
!                 + snow_bin2  * hg(0.0     ,elong_a(i,j))  

          if(cloud_od_tot .gt. 0.)then
              snow_factor = trans_nonsnow * (cloud_od_snow / cloud_od_tot)
          else
              snow_factor = 0.
          endif

          pf_scat2(i,j) = pf_snow * snow_factor + pf_scat1(i,j) * (1.0 - snow_factor)

!         Suppress phase function if terrain is close in the light ray
          if(airmass_2_topo(i,j) .gt. 0.)then ! cloud in front of terrain
              pf_scat(i,j) = pf_scat2(i,j)**(r_cloud_rad(i,j)**2.0)
!             pf_scat(i,j) = pf_scat2(i,j)**opac(cloud_od_tot)
          else
              pf_scat(i,j) = pf_scat2(i,j)
          endif

!         Add simple rainbow (red is 42 degrees radius, blue is 40)
          if(abs(elong_a(i,j)-149.) .lt. 1.0)then
              rain_factor = opac(cloud_od_sp(i,j,3))
              pf_scat(i,j) = pf_scat(i,j) * (1.0 + 3. * rain_factor * r_cloud_rad(i,j)**10.)
          endif

        enddo ! j
        enddo ! i

        return
        end

