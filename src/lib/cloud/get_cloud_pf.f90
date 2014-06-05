
        subroutine get_cld_pf(elong_a,r_cloud_rad,cloud_od,cloud_od_sp,nsp,airmass_2_topo,idebug_a,ni,nj & ! I
                             ,pf_scat1,pf_scat2,pf_scat,pf_thk_a) ! O
        include 'trigd.inc'

!       Statement functions
        trans(od) = exp(-od)
        opac(od) = 1.0 - exp(-od)
        counts_to_rad(counts) = 10.**((counts/100.)+7.3)
!       pf_bk(elong,alb) = 1.0 + cosd(180. - elong)
        alb(bt) = bt / (1.+bt)

        include 'rad.inc'
        real elong_a(ni,nj)
        real cloud_od(ni,nj)        ! cloud optical depth (tau)
        real cloud_od_sp(ni,nj,nsp) ! cloud species tau (clwc,cice,rain,snow)
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        integer idebug_a(ni,nj)
        real pf_scat(nc,ni,nj), pf_scat1(ni,nj), pf_scat2(ni,nj)
        real pf_thk_a(ni,nj)
        real*8 phase_angle_d,phase_corr

        real anglebow1(nc) /41.80,41.20,40.60/
        real anglebow2(nc) /52.50,53.20,53.90/

        logical l_pf_new /.true./

!       Call to phase_func
        dimension v10(0:13,0:70)
        character*10 nm(0:13,0:70)

        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        write(6,*)' max of idebug_a (2) = ',maxval(idebug_a)

        do j = 1,nj
         do i = 1,ni

!         a substitute for cloud_rad could be arg2 = (cosd(alt))**3.

!         Phase function that depends on degree of fwd scattering in cloud    
!         pwr controls angular extent of fwd scattering peak of a thin cloud
!         General references: 
!         http://www.uni-leipzig.de/~strahlen/web/research/de_index.php?goto=arctic
!         http://www-evasion.imag.fr/Publications/2008/BNMBC08/clouds.pdf
!         http://pubs.giss.nasa.gov/docs/1969/1969_Hansen_2.pdf
!         http://rtweb.aer.com/rrtm_frame.html
!         http://wiki.seas.harvard.edu/geos-chem/images/Guide_to_GCRT_Oct2013.pdf

          cloud_od_cice = cloud_od_sp(i,j,2)

          if(.true.)then ! new liquid phase function
            cloud_od_liq = cloud_od_sp(i,j,1) + cloud_od_sp(i,j,3)

            hgp = max(cloud_od_liq,1.0)       ! multiple scattering phase func
            sco = hgp

!           Check for being optically thin with direct light source
            clwc_bin1 = exp(-(cloud_od_liq/10.)**2) ! optically thin
            clwc_bin1 = clwc_bin1 * r_cloud_rad(i,j)**10. ! direct lighting
            clwc_bin1a = 0.80**sco             ! forward peak
!           clwc_bin1c = opac(0.00*sco)        ! backscattering
            clwc_bin1b = 1.0 - (clwc_bin1a)    ! mid    
            clwc_bin2 = 1.0 - clwc_bin1        ! optically thick 

!           Derived using phase function of Venus (a thick cloud analog)
!           Specifically the Venus phase angle magnitude corr
            iplan = 3
            isat = 0
            phase_angle_d = 180. - elong_a(i,j)
            call phase_func(iplan,isat,phase_angle_d,v10,nm,phase_corr)       
            r_ill = (1. + cosd(phase_angle_d)) / 2.

!           Set to 1 or 2 if we're at cloud base?
            pf_thk = (1.94 / (10.**(phase_corr * 0.4))) / r_ill
            pf_thk = min(pf_thk,2.0) ! limit fwd scattering peak

!           This might need expanded support for increased forward scattering 
!           when looking down on a layer of cloud from the air.
!           Eventually sky average r_cloud_rad can help decide the regime for
!           determination of pf_thk? 'scurve' is also available if needed.

!           if(elong_a(i,j) .gt. 90.)then
!               pf_thk = 1.0 +  (elong_a(i,j)-90.)/90.
!           else
!               pf_thk = 1.0 + ((90.-elong_a(i,j))/90.)**2
!           endif
            alb_clwc = alb(0.06*cloud_od_liq)

            radfrac = scurve(r_cloud_rad(i,j)**3) ! high for illuminated clouds
!                     illuminated                unilluminated
            pf_thk = pf_thk*radfrac + hg(-0.0 ,elong_a(i,j)) * (1.-radfrac)
            pf_thk_a(i,j) = pf_thk

            pf_clwc = clwc_bin1a * clwc_bin1 * hg(.94**hgp,elong_a(i,j)) & ! corona
                    + clwc_bin1b * clwc_bin1 * hg(.60**hgp,elong_a(i,j)) &
!                   + clwc_bin1c * clwc_bin1 * hg(.00     ,elong_a(i,j)) &
!                   + clwc_bin2  * hg(-0.3    ,elong_a(i,j))  
!                   + clwc_bin2  * 2.0    ! add albedo term?     
                    + clwc_bin2  * pf_thk ! add albedo term?     

!           clwc_bin1c = .00
            clwc_bin1a = clwc_bin1a - .00
            clwc_bin1b = clwc_bin1b - .00

            pf_rain = clwc_bin1a * clwc_bin1 * hg(.94**hgp,elong_a(i,j)) &         
                    + clwc_bin1b * clwc_bin1 * hg(.60**hgp,elong_a(i,j)) &
!                   + clwc_bin1c * clwc_bin1 * hg(.00     ,elong_a(i,j)) &
!                   + clwc_bin2  * hg(-0.3    ,elong_a(i,j))  
!                   + clwc_bin2  * 2.0    ! add albedo term?     
                    + clwc_bin2  * pf_thk ! add albedo term?     

            if(cloud_od_liq .eq. cloud_od_sp(i,j,1))then ! cloud liquid
                pf_scat1(i,j) = pf_clwc
            else                                         ! rain
                pf_scat1(i,j) = pf_rain
            endif

          endif ! .true.

!         Separate snow (+cice) phase function
          cloud_od_snow = cloud_od_sp(i,j,4) + cloud_od_sp(i,j,2)
          trans_nonsnow = trans(cloud_od(i,j) - cloud_od_snow)
          cloud_od_tot  = cloud_od(i,j)

          hgp = max(cloud_od_snow,1.0)       ! multiple scattering phase func
          sco = hgp

!         snow_bin1 = exp(-cloud_od_snow/5.) ! optically thin snow
          snow_bin1a = 0.60**sco             ! forward peak
          snow_bin1c = opac(0.10*sco)        ! backscattering (thick)
          snow_bin1b = 1.0 - (snow_bin1a + snow_bin1c) ! mid
!         snow_bin2 = 1.0 - snow_bin1        ! optically thick snow

          pf_snow = snow_bin1a * hg(.98**hgp,elong_a(i,j)) & ! plates
                  + snow_bin1b * hg(.50**hgp,elong_a(i,j)) &
                  + snow_bin1c * pf_thk
!                 + snow_bin2  * hg(0.0     ,elong_a(i,j))  

          if(cloud_od_tot .gt. 0.)then
              snow_factor = trans_nonsnow * (cloud_od_snow / cloud_od_tot)
          else
              snow_factor = 0.
          endif

          pf_scat2(i,j) = pf_snow * snow_factor + pf_scat1(i,j) * (1.0 - snow_factor)

!         Suppress phase function if terrain is close in the light ray
          if(airmass_2_topo(i,j) .gt. 0.)then ! cloud in front of terrain
              pf_scat(:,i,j) = pf_scat2(i,j)**(r_cloud_rad(i,j)**2.0)
!             pf_scat(:,i,j) = pf_scat2(i,j)**opac(cloud_od_tot)
          else
              pf_scat(:,i,j) = pf_scat2(i,j)
          endif

!         Add rainbows
!         rain_factor = opac(cloud_od_sp(i,j,3))
!         rain_factor = cloud_od_sp(i,j,3) / cloud_od_tot
          if(cloud_od_sp(i,j,3) .gt. 0.)then
            rain_factor = 1.
          else
            rain_factor = 0.
          endif

          if(rain_factor .gt. 0.)then
            antisol_rad = 180. - elong_a(i,j)
            frac_single = 0.7 + 0.3 * clwc_bin1 ! optically thin
            do ic = 1,nc

!             Reflection inside primary rainbow + supernumerary rainbows
              bow_hw = 0.8
              if(antisol_rad .lt. anglebow1(ic)-bow_hw)then
!                 Calculate super_omega (degrees phase per degree elongation/radius)
                  super_omega = 450. + (float(ic-1) * 75.) 
                  super_phase = ((anglebow1(ic)-bow_hw) - antisol_rad) * super_omega
                  super_amp = 0.1 - (0.1 * cosd(super_phase) * (antisol_rad / (anglebow1(ic)-bow_hw))**15.)
                  bow_int = (2.0 + 1.0 * antisol_rad / (anglebow1(ic)-bow_hw)) * super_amp
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (.053 * bow_int * rain_factor * r_cloud_rad(i,j)**10.)
              endif

!             Primary rainbow (red is 42 degrees radius, blue is 40)
              bow_hw = 0.8
              if(abs(antisol_rad-anglebow1(ic)) .lt. bow_hw)then
                  bow_frac = abs(antisol_rad-anglebow1(ic)) / bow_hw
!                 ratio_bow = max(1.0 / pf_scat(ic,i,j) - 1.0, 0.)
                  bow_int = (1.0 - bow_frac) ** 1.0 
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (1.0 * frac_single * bow_int * rain_factor * r_cloud_rad(i,j)**10.)
              endif

!             Alexander's dark band
              angle_alex = (anglebow1(ic) + anglebow2(ic)) / 2.0     
              alex_hw = 0.5 * (anglebow2(ic)-anglebow1(ic)) - 1.0*bow_hw
              if(abs(antisol_rad-angle_alex) .lt. alex_hw)then
                  frac_alex1 = abs(antisol_rad-angle_alex) / alex_hw
                  frac_alex2 = 1.0 - frac_alex1**8
                  frac_alex = frac_alex2 * rain_factor * r_cloud_rad(i,j)**10.
                  frac_alex = frac_alex * clwc_bin1 ! optically thin
                  pf_alex = .02
                  pf_scat(ic,i,j) = pf_alex * frac_alex + pf_scat(ic,i,j) * (1.0 - frac_alex)
              endif

!             Secondary rainbow (red is 54.5 degrees radius, blue is 52)
              bow_hw = 0.8
              if(abs(antisol_rad-anglebow2(ic)) .lt. bow_hw)then
                  bow_frac = abs(antisol_rad-anglebow2(ic)) / bow_hw 
!                 ratio_bow = max(0.2 / pf_scat(ic,i,j) - 1.0, 0.)
                  bow_int = (1.0 - bow_frac) ** 0.5 
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (0.2 * frac_single * bow_int * rain_factor * r_cloud_rad(i,j)**10.)
              endif
            enddo ! ic

          endif ! rain_factor > 0

          if(idebug_a(i,j) .eq. 1)then
              write(6,101)i,j,pf_clwc,pf_rain,r_cloud_rad(i,j),radfrac,pf_scat1(i,j),trans_nonsnow,snow_factor,rain_factor,pf_scat(2,i,j)
101           format(' clwc/rain/rad/radf/pf1/trans/snow/rain factors = ',2i5,2f9.3,2x,2f9.3,2x,4f9.3,f9.3)
          endif

         enddo ! i (altitude)

         if(j .eq. nj/2)then
             write(6,*)'get_cld_pf (zenith) ',pf_scat(:,ni,j)
         endif

        enddo ! j (azimuth)

        return
        end

