
        subroutine get_scat_pf( &
                elong_a,alt_a,aero_od_src,aero_od_obs,aero_ssa,aero_g      & ! I
               ,r_cloud_rad,cloud_rad_w,cloud_od,cloud_od_sp,nsp           & ! I
               ,cloud_od_sp_w                                              & ! I
               ,emis_ang_a,airmass_2_topo,idebug_a,ni,nj                   & ! I
               ,pf_scat,mtr_msa)                                             ! O

!       subroutine get_cld_pf(elong_a,alt_a,r_cloud_rad,cloud_rad_w,cloud_od,cloud_od_sp &
!                            ,emis_ang_a,nsp,airmass_2_topo,idebug_a,ni,nj & ! I
!                            ,pf_scat1,pf_scat2,pf_scat,pf_thk_a) ! O
        include 'trigd.inc'

        use mem_namelist, ONLY: aod_asy, fcterm

!       Statement functions
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)
        alb(bt) = bt / (1.+bt)
        rad2tau(b,r) = (1.-r)/(r*b)

        include 'rad_nodata.inc'

        real asy_clwc(nc)  /0.945,0.950,0.955/

        pf_thn_clwcf(hgpf,elgf) & 
             = clwc_bin1a * hg(asy_clwc(ic)**hgpf,elgf)& ! corona
             + clwc_bin1b * hg(.60**hgpf         ,elgf)&
             + clwc_bin1c * hg(-.60              ,elgf)&
             + clwc_bin1d * hg(-.00              ,elgf)

        pf_thn_rainf(hgpf,elgf) & 
             = rain_bin1a * hg(0.99**hgpf,elgf)&
             + rain_bin1b * hg(0.75**hgpf,elgf)&
             + rain_bin1c * hg(0.00      ,elgf)&
             + rain_bin1d * hg(-.20      ,elgf)

        pf_thn_snowf(hgpf,elgf) & 
             = arg1 * hg( .999**hgpf,elgf) & ! plates
             + arg2 * hg( .860**hgpf,elgf) &
             + arg3 * hg( .000      ,elgf) & ! isotropic
             + arg4 * hg(-.600      ,elgf)   ! backscat

!       parameter (nsp=4) ! ERROR when passed into 'get_scat_pf'
        
        real elong_a(ni,nj)
        real alt_a(ni,nj)
        real emis_ang_a(ni,nj)        ! 
        real cloud_od(ni,nj)          ! cloud optical depth (tau)
        real cloud_od_sp(ni,nj,nsp)   ! cloud species tau (clwc,cice,rain,snow)
        real cloud_od_sp_w(ni,nj,nsp) ! cloud species tau (clwc,cice,rain,snow)
        real r_cloud_rad(ni,nj)       ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_w(ni,nj)       ! sun to cloud transmissivity (direct+fwd scat) * rad

!       Meteors (mtr) are a composite of hydrometeors and lithometeors (aerosols)
!       From ground the src OD is important near the sun and obs OD less so    
!       From space  the src OD is important near the sun and obs OD less so    
!       From ground the src & obs OD are both important far from the sun
!       From space  the src & obs OD are both important far from the sun
!       Aero ssa is presently applied in get_skyglow, add in 'get_cloud_pf'?   
!       Aero msa is presently applied in get_cloud_rad, add only SSA in 'get_cloud_pf'?

!       Thin clouds
!       High cloud rad, phase angle is driven by single scattering

!       Thick clouds near the sun from ground
!       Low cloud rad with src OD, phase function mostly from altitude

!       Thick clouds near the sun from space
!       High cloud rad with src OD, phase function relates to reflectance and ARF
!            scatter order and MSA from obs OD

!       Thin aerosols
!       High cloud rad, phase angle is driven by single scattering, corrected by SSA
!            scatter order and MSA from obs OD

!       Thick aerosols near the sun from ground
!       Low cloud rad with src OD, phase function mostly from altitude, corrected by MSA
!            scatter order and MSA from src OD (via cloud_rad)

!       Thick aerosols near the sun from space
!       High cloud rad with src OD, phase function relates to reflectance and ARF, corrected by MSA
!            scatter order and MSA from obs OD

        real cloud_od_src_sp(ni,nj,nsp)
!       real cloud_od_obs_sp(ni,nj,nsp) 
        real cloud_ssa_sp(nsp)
        real cloud_asy_sp(nsp)
        real aero_od_src(nc,ni,nj)
        real aero_od_obs(nc,ni,nj)
        real aero_ssa(nc,ni,nj)
        real aero_asy(nc,ni,nj)
        real hg2d(nc), hg2(nc), hg2t(nc)
        real aero_factor(nc)
        real mtr_od_obs(nc,ni,nj)
        real mtr_ssa(nc,ni,nj)
        real mtr_asy(nc,ni,nj)
        real mtr_msa(nc,ni,nj)

        real airmass_2_topo(ni,nj)  ! airmass to topo  
        integer idebug_a(ni,nj)
        real pf_scat(nc,ni,nj), pf_scat1(nc,ni,nj), pf_scat2(nc,ni,nj)
        real pf_clwc(nc),pf_rain(nc)
        real pf_thk_a(ni,nj)
        real*8 phase_angle_d,phase_corr

        real anglebow1(nc) /41.80,41.20,40.60/
        real anglebow2(nc) /52.50,53.20,53.90/

        logical l_pf_new /.true./

!       Call to phase_func
        dimension v10(0:13,0:70)
        character*10 nm(0:13,0:70)

        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        write(6,*)' get_cld_pf: max of idebug_a (2) = ',maxval(idebug_a)
        write(6,*)'                                                                                         alt    elg      cod    pf_thn_c  pf_thk   clwc     rain      rad     radw    radf       pf1     pf2     pfs    trans  sn fctr rn fctr    pf2'
        cloud_ssa_sp(:) = 1.0
        cloud_asy_sp(:) = 0.8

        do j = 1,nj
         do i = 1,ni

          do ic = 1,nc
!           Composite meteor optical properties seen from observer
            mtr_od_obs(ic,i,j) = sum(cloud_od_sp(i,j,:)) + aero_od_obs(ic,i,j)
            mtr_ssa(ic,i,j) = (sum(cloud_od_sp(i,j,:) * cloud_ssa_sp(:)) + aero_od_obs(ic,i,j) * aero_ssa(ic,i,j)) / mtr_od_obs(ic,i,j)
            mtr_asy(ic,i,j) = (sum(cloud_od_sp(i,j,:) * cloud_asy_sp(:)) + aero_od_obs(ic,i,j) * aero_asy(ic,i,j)) / mtr_od_obs(ic,i,j)

            scatter_order = sqrt(1.0 + mtr_od_obs(ic,i,j)**2)
            mtr_msa(ic,i,j) = mtr_ssa(ic,i,j)**scatter_order

!           Aerosol phase function
            fb = aod_asy(3,ic)**scatter_order
            g1 = aod_asy(2,ic)**scatter_order
            g2 = aod_asy(1,ic)**scatter_order
            hg2(ic) = dhg2(elong_a(i,j),fb,fcterm)

            if(mtr_od_obs(ic,i,j) .gt. 0.)then
              aero_factor(ic) = aero_od_obs(ic,i,j)/mtr_od_obs(ic,i,j)
            else
              aero_factor(ic) = 0.
            endif
            
          enddo ! ic
          
!         a substitute for cloud_rad could be arg2 = (cosd(alt))**3.

!         Phase function that depends on degree of fwd scattering in cloud    
!         pwr controls angular extent of fwd scattering peak of a thin cloud
!         General references: 
!         http://www.uni-leipzig.de/~strahlen/web/research/de_index.php?goto=arctic
!         http://www-evasion.imag.fr/Publications/2008/BNMBC08/clouds.pdf
!         http://pubs.giss.nasa.gov/docs/1969/1969_Hansen_2.pdf
!         http://rtweb.aer.com/rrtm_frame.html
!         http://wiki.seas.harvard.edu/geos-chem/images/Guide_to_GCRT_Oct2013.pdf
!         https://www.osapublishing.org/oe/abstract.cfm?uri=oe-23-9-11995

          cloud_od_cice = cloud_od_sp(i,j,2)
          cloud_od_rain = cloud_od_sp(i,j,3)

!         Front weighted species factors
          sumspw = sum(cloud_od_sp_w(i,j,:)) 
          if(sumspw .gt. 0.)then
             clwc_factor_w = cloud_od_sp_w(i,j,1) / sumspw
             cice_factor_w = cloud_od_sp_w(i,j,2) / sumspw
             rain_factor_w = cloud_od_sp_w(i,j,3) / sumspw
             snow_factor_w = cloud_od_sp_w(i,j,4) / sumspw
             aero_factor_w = 0. ! approximate
          else
             clwc_factor_w = 0.
             cice_factor_w = 0.
             rain_factor_w = 0.
             snow_factor_w = 0.
             aero_factor_w = 1. ! approximate
          endif

!         Parameter for looking at normal face of a cloud is set to 1 if
!         we're looking above the horizon, or below the horizon looking
!         straight down. It is 0 if we're looking much below the horizon near
!         the limb.
          if(alt_a(i,j) .le. 0. .and. emis_ang_a(i,j) .gt. 0.)then
              frac_alt = sind(abs(alt_a(i,j)))
              frac_norm = frac_alt * sind(emis_ang_a(i,j)) + (1.-frac_alt) 
          else ! above the horizon
              frac_norm = 1.
          endif

          if(.true.)then ! new liquid phase function
            cloud_od_liq = cloud_od_sp(i,j,1) + cloud_od_sp(i,j,3)
!           cloud_od_liq = cloud_od_liq * frac_norm**2

            hgp = max(cloud_od_liq,1.0)    ! multiple scattering phase func
            sco = hgp

            bf = 0.06

!           Check for being optically thin with direct light source
!           Coefficient of 10. can be lower if needed
            clwc_bin1 = exp(-(cloud_od_liq/10.)**2) ! optically thin
            clwc_albedo = (bf * cloud_od_liq) / (1. + (bf * cloud_od_liq))
            clwc_bin1 = clwc_bin1 * r_cloud_rad(i,j)**10. ! direct lighting
            clwc_bin1a = 0.80**sco                  ! forward peak
            clwc_bin1b =  1.06 - (clwc_bin1a)       ! mid    
            clwc_bin1c =  0.02 ! opac(0.00*sco)     ! backscattering
            clwc_bin1d = -0.08 ! opac(0.00*sco)     ! backscattering
            clwc_bin2 = 1.0 - clwc_bin1             ! optically thick 

!           Derived using phase function of Venus (a thick cloud analog)
!           Specifically the Venus phase angle magnitude corr
            iplan = 3
            isat = 0
            phase_angle_d = 180. - elong_a(i,j)
            call phase_func(iplan,isat,phase_angle_d,v10,nm,phase_corr)       
            r_ill = (1. + cosd(phase_angle_d)) / 2.

!           Set to 1 or 2 if we're at cloud base?
            pf_thk_hr = (1.94 / (10.**(phase_corr * 0.4))) / r_ill
            pf_thk_hr = min(pf_thk_hr,600.0) ! limit fwd scattering peak
!           pf_thk_hr = min(pf_thk_hr,2.0)   ! limit fwd scattering peak

!           Eventually sky average r_cloud_rad can help decide the regime for
!           determination of pf_thk? 'scurve' is also available if needed.

!           if(elong_a(i,j) .gt. 90.)then
!               pf_thk = 1.0 +  (elong_a(i,j)-90.)/90.
!           else
!               pf_thk = 1.0 + ((90.-elong_a(i,j))/90.)**2
!           endif

!           radfrac (high for illuminated clouds) may have multiple
!           purposes including suppressing phase function near terrain
!           and ensuring high pf for backscattering with thick clouds
!           radfrac = scurve(r_cloud_rad(i,j)**3) ! high for illuminated clouds
            radfrac = scurve(scurve(scurve(scurve(r_cloud_rad(i,j)))))
            elgfrac = scurve(elong_a(i,j)/180.) ! * radfrac ! need radfrac?
            alb_clwc = alb(0.06*cloud_od_liq) * elgfrac + (1.-elgfrac)

            pf_thk_lr = (2./3.) * (1. + sind(abs(alt_a(i,j))))
            pf_thk_lr = 1. + 2. * (pf_thk_lr-1.)
!                     illuminated                unilluminated
!           pf_thk = pf_thk*radfrac + hg(-0.,elong_a(i,j)) * (1.-radfrac) &
!                                       * (2./3. * (1. + sind(alt_a(i,j))))
            pf_thk = pf_thk_hr*radfrac*alb_clwc + 2.*pf_thk_lr*(1.-radfrac) 
            pf_thk_a(i,j) = pf_thk

            rain_peak = exp(-rad2tau(.06,r_cloud_rad(i,j))/2.)
            rain_bin1 = exp(-(cloud_od_liq/10.)**2) ! optically thin
            rain_bin1 = rain_bin1 * rain_peak ! direct lighting

            do ic = 1,nc
              pf_thn_clwc = pf_thn_clwcf(hgp,elong_a(i,j)) * (1.-aero_factor(ic)) + hg2(ic) * aero_factor(ic)

              pf_clwc(ic) &
             = clwc_bin1  * pf_thn_clwc &
!            + clwc_bin2  * hg(-0.3    ,elong_a(i,j))  
!            + clwc_bin2  * 2.0    ! add albedo term?     
             + clwc_bin2  * pf_thk ! add albedo term?     

              rain_bin1a =   .10
              rain_bin1b =  1.05
              rain_bin1c = -0.35
              rain_bin1d =   .20

              pf_thn_rain = pf_thn_rainf(hgp,elong_a(i,j)) 

              pf_rain(ic) &
                     = rain_bin1 * pf_thn_rain &
                     + clwc_bin2 * pf_thk ! add albedo term?     

              if(cloud_od_liq .eq. cloud_od_sp(i,j,1))then ! cloud liquid
                  pf_scat1(ic,i,j) = pf_clwc(ic)
              else                                         ! rain
                  pf_scat1(ic,i,j) = pf_rain(ic)
              endif

            enddo ! ic

          endif ! .true.

!         Separate snow (+cice) phase function
          cloud_od_snow = cloud_od_sp(i,j,4) + cloud_od_sp(i,j,2)
          trans_nonsnow = trans(cloud_od(i,j) - cloud_od_snow)
          cloud_od_tot  = cloud_od(i,j)
!         cloud_od_snow = cloud_od_snow * frac_norm**2

          if(cloud_od_snow .gt. 0.)then
              fsnow = cloud_od_sp(i,j,4) / cloud_od_snow
              fcice = cloud_od_sp(i,j,2) / cloud_od_snow
          else
              fsnow = 1.0; fcice = 0.0
          endif

          sco = max(cloud_od_snow,1.0) ! multiple scatter order 
          hgp = sco                    ! multiple scatter order

          snow_bin1 = exp(-cloud_od_snow/5.)      ! optically thin snow
          snow_bin1c = opac(0.10*(sco**1.2 - 1.0))  ! backscattering (thick)
          snow_bin1a = (1. - snow_bin1c)
!         snow_bin1a = 0.50**sco             ! forward peak
!         snow_bin1b = 1.0 - (snow_bin1a + snow_bin1c) ! mid
          snow_bin2 = 1.0 - snow_bin1        ! optically thick snow

          arg1 = (0.50 * fsnow + 0.50 * fcice) * 2.**(-(hgp-1.0))
          arg3 =  0.03 * fsnow + 0.20 * fcice
          arg4 =  0.02 * fsnow + 0.02 * fcice
!         arg2 =  0.45 * fsnow + 0.28 * fcice
          arg2 =  1.0 - (arg1 + arg3 + arg4)

          pf_thn_snow = pf_thn_snowf(hgp,elong_a(i,j)) 

          pf_snow = snow_bin1a * pf_thn_snow &
                  + snow_bin1c * pf_thk
!                 + snow_bin2  * hg(0.0     ,elong_a(i,j))  

          elgfrac = scurve(elong_a(i,j)/180.) * snow_bin2
          alb_snow = alb(0.14*cloud_od_snow) ! * elgfrac + (1.-elgfrac)
          pf_snow = pf_snow * (1.-elgfrac) + 2. * alb_snow * elgfrac

          if(cloud_od_tot .gt. 0.)then
              snow_factor = trans_nonsnow * (cloud_od_snow / cloud_od_tot)
          else
              snow_factor = 0.
          endif

          pf_scat2(:,i,j) = pf_snow * snow_factor + pf_scat1(:,i,j) * (1.0 - snow_factor)

!         Suppress/cap phase function if terrain is close in the light ray
!         Apply when radfrac is low and opacity is low
!         if(airmass_2_topo(i,j) .gt. 0. .and. pf_scat2(2,i,j) .gt. 1.0)then ! cloud in front of terrain
          if(airmass_2_topo(i,j) .gt. 0.)then ! cloud in front of terrain
!             pf_scat(:,i,j) = pf_scat2(:,i,j)**opac(cloud_od_tot)
!             pf_scat(:,i,j) = pf_scat2(:,i,j)**(r_cloud_rad(i,j)**2.0)
!             pf_scat(:,i,j) = pf_scat2(:,i,j) * radfrac + 2. * pf_thk_lr * (1.-radfrac)
              arg_tn = (1. - radfrac) * (1. - opac(cloud_od_tot))
!             pf_scat(:,i,j) = pf_scat2(:,i,j) * radfrac + 2. * pf_thk_lr * (1.-radfrac)
              pf_scat(:,i,j) = pf_scat2(:,i,j) * (1.-arg_tn) + 2. * pf_thk_lr * arg_tn
          else
              pf_scat(:,i,j) = pf_scat2(:,i,j)
          endif

!         Add rainbows
!         Ramping this in between an OD range of .01 to .1
          if(cloud_od_sp(i,j,3) .gt. 0.1)then
              rain_factor1 = 1.0
          elseif(cloud_od_sp(i,j,3) .gt. 0.01)then
              rain_factor1 = log10(cloud_od_sp(i,j,3) + 2.)
          else ! < .01
              rain_factor1 = 0.0
          endif
!         rain_factor = cloud_od_sp(i,j,3) / cloud_od_tot
!         if(cloud_od_sp(i,j,3) .gt. 0.)then
!           rain_factor = 1.
!         else
!           rain_factor = 0.
!         endif

          rain_factor = rain_factor_w
          
          cre = 3.5

!         Call mscat routine here vs in 'get_cloud_rad'?         

          if(rain_factor .gt. 0.)then
            antisol_rad = 180. - elong_a(i,j)
            frac_single = 1.0 ! 0.7 + 0.3 * clwc_bin1 ! optically thin
!           frac_single_hirad = 1.0 - 0.3 * opac(cloud_od) ! hirad
!    1                        + trans(cloud_od - 1.)       ! lowrad
            do ic = 1,nc

!             Reflection inside primary rainbow + supernumerary rainbows
              bow_hw = 0.8
              if(antisol_rad .lt. anglebow1(ic))then
!                 Calculate super_omega (degrees phase per degree elongation/radius)
                  super_omega = 450. + (float(ic-1) * 75.) 
                  super_phase = (anglebow1(ic) - antisol_rad) * super_omega
                  super_amp = 0.4 * (antisol_rad / (anglebow1(ic)))**60.
                  bow_int = (2.0 + 1.0 * antisol_rad / (anglebow1(ic)-bow_hw)) * super_amp
                  horn1dist = anglebow1(ic) - antisol_rad
                  horn1 = 1.0 * exp(-horn1dist/5.)
                  horn1 = horn1 * (1. + super_amp * (cosd(super_phase)-1.)/2.)
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (3.0 * horn1 * rain_factor * cloud_rad_w(i,j)**cre)
                  if(ic .eq. 2 .and. idebug_a(i,j) .eq. 1)then
                      write(6,91)antisol_rad,horn1,rain_factor,cloud_rad_w(i,j)**3.,super_amp
91                    format(' antisol_rad,horn1,rainfact,cldrdw3,sup',f9.3,f9.6,3f9.3)
                  endif
              endif

!             Primary rainbow (red is 42 degrees radius, blue is 40)
              bow_hw = 0.8
              if(abs(antisol_rad-anglebow1(ic)) .le. bow_hw .AND. antisol_rad .ge. anglebow1(ic))then
                  bow_frac = abs(antisol_rad-anglebow1(ic)) / bow_hw
!                 ratio_bow = max(1.0 / pf_scat(ic,i,j) - 1.0, 0.)
                  bow_int = (1.0 - bow_frac) ** 1.0 
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (3.0 * frac_single * bow_int * rain_factor * cloud_rad_w(i,j)**cre)
              endif

!             Alexander's dark band
              angle_alex = (anglebow1(ic) + anglebow2(ic)) / 2.0     
              alex_hw = 0.5 * (anglebow2(ic)-anglebow1(ic)) - 1.0*bow_hw
              if(abs(antisol_rad-angle_alex) .lt. alex_hw)then
                  frac_alex1 = abs(antisol_rad-angle_alex) / alex_hw
                  frac_alex2 = 1.0 - frac_alex1**8
                  frac_alex = frac_alex2 * rain_factor * cloud_rad_w(i,j)**cre ! direct lighting
!                 frac_alex = frac_alex2 * rain_factor * rain_peak            ! direct lighting
!                 frac_alex = frac_alex * rain_bin1 ! optically thin
                  pf_alex = .02
                  pf_scat(ic,i,j) = pf_alex         * frac_alex &
                                  + pf_scat(ic,i,j) * (1.0 - frac_alex)
                  if(ic .eq. 2 .and. idebug_a(i,j) .eq. 1)then
                      write(6,93)rain_peak,rain_bin1,frac_alex2,frac_alex
93                    format(' rain_peak/rain_bin1/alex2/alex',2f9.6,f9.3,f9.6)
                  endif
              endif

!             Secondary rainbow (red is 54.5 degrees radius, blue is 52)
              bow_hw = 0.8
              if(abs(antisol_rad-anglebow2(ic)) .le. bow_hw .AND. antisol_rad .le. anglebow2(ic))then
                  bow_frac = abs(antisol_rad-anglebow2(ic)) / bow_hw 
                  bow_int = (1.0 - bow_frac) ** 0.5 
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (0.56* frac_single * bow_int * rain_factor * cloud_rad_w(i,j)**cre)
              endif

!             Reflection outside secondary rainbow + supernumerary rainbows
              bow_hw = 0.8
              if(antisol_rad .gt. anglebow2(ic) .AND. antisol_rad .lt. anglebow2(ic) + 30.)then
!                 Calculate super_omega (degrees phase per degree elongation/radius)
                  super_omega = 450. + (float(ic-1) * 75.) 
                  super_phase = (antisol_rad - anglebow2(ic)) * super_omega
                  super_amp = 0.4 * (1.0 - (antisol_rad - anglebow2(ic))/30.)**15.
                  super_amp = max(super_amp,0.0)
!                 super_amp = (1.0 - (antisol_rad - anglebow2(ic))/30.)**15.
                  bow_int = (2.0 + 1.0 * antisol_rad / (anglebow1(ic)-bow_hw)) * super_amp
                  horn2dist = antisol_rad - anglebow2(ic)
                  horn2 = 1.0 * exp(-horn2dist/2.0)
                  horn2 = horn2 * (1. + super_amp * (cosd(super_phase)-1.)/2.)
                  pf_scat(ic,i,j) = pf_scat(ic,i,j) + (0.56 * horn2 * rain_factor * cloud_rad_w(i,j)**cre)
                  if(ic .eq. 2 .and. idebug_a(i,j) .eq. 1)then
                      write(6,94)antisol_rad,horn2
94                    format(' antisol_rad,horn2',f9.3,f9.6)
                  endif
              endif

            enddo ! ic

          endif ! rain_factor > 0

          if(idebug_a(i,j) .eq. 1)then
!             write(6,101)i,j,alt_a(i,j),elong_a(i,j),cloud_od_tot,pf_thk_hr,pf_thk,pf_clwc(2),pf_rain(2),r_cloud_rad(i,j),cloud_rad_w(i,j),radfrac,pf_scat1(2,i,j),pf_scat2(2,i,j),pf_scat(2,i,j),trans_nonsnow,snow_factor,rain_factor,pf_scat(2,i,j)
              write(6,101)i,j,alt_a(i,j),elong_a(i,j),cloud_od_tot,pf_thn_clwc,pf_thk,pf_clwc(2),pf_rain(2),r_cloud_rad(i,j),cloud_rad_w(i,j),radfrac,pf_scat1(2,i,j),pf_scat2(2,i,j),pf_scat(2,i,j),trans_nonsnow,snow_factor,rain_factor,pf_scat(2,i,j)
101           format(' alt/elg/cod/thnc/thk/clwc/rain/rad/radw/radf/pf1/pf2/pfs/trans/sn/rn fctrs = ',i4,i5,f6.1,f8.2,5f9.3,2x,3f8.4,2x,6f8.3,f9.3)
              fb = aod_asy(3,2)**scatter_order
              g1 = aod_asy(2,2)**scatter_order
              g2 = aod_asy(1,2)**scatter_order
              write(6,102)bf,cloud_od_liq,clwc_bin2,alb_clwc,frac_norm,clwc_factor_w,cice_factor_w,rain_factor,snow_factor_w,aero_factor(2),scatter_order,fb,g1,g2,hg2(2)
102           format(' bf/od/clwc_bin2/alb_clwc/fnrm/',f9.3,f9.4,2f12.6,f9.4,' factorw ',4f9.4,' aerof/sct/hg2 ',6f9.2)
          endif

         enddo ! i (altitude)

         if(j .eq. nj/2)then
             write(6,*)' get_cld_pf (zenith/nadir) ',pf_scat(:,ni,j)
         endif

        enddo ! j (azimuth)

        return
        end

