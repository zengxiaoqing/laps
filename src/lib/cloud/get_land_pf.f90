
        subroutine get_lnd_pf(elong_a,alt_a,azi_a &                     ! I
                             ,topo_gti,topo_albedo &                    ! I
                             ,topo_lf,topo_sc,transm_obs &              ! I
                             ,gtic,dtic,btic &                          ! I
                             ,dist_2_topo,topo_solalt,topo_solazi,azi_scale &       ! I
                             ,sol_alt,sol_azi,nsp,airmass_2_topo &      ! I
                             ,idebug_pf,ni,nj,i4time,rlat,rlon,htmsl &  ! I
                             ,topo_lat,topo_lon &                       ! I
                             ,pf_land,cld_brdf,emis_ang_a)              ! O

        use mem_namelist, ONLY: r_missing_data,earth_radius
        use cloud_rad, ONLY: ghi_zen_toa, zen_kt
        include 'trigd.inc'
        include 'rad.inc'

        parameter (sph_sqdg = 4. * 180.**2 / pi)

!       Statement functions
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)
        alb(bt) = bt / (1.+bt)
        rad2tau(b,r) = (1.-r)/(r*b)
        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.
        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1
        ph_exp(ampl1,azidiff1) = exp(ampl1 * cosd(azidiff1)) / (1. + abs(ampl1)*.16)
        hg_cyl(g,pha) = hg(g,pha) / (1. + 1.3*g**2)   ! integrate to ~1
        fslope(x,s) = s*x - (s-1.)*x**((s+1.)/s)
        angdist(p1,p2,dlon) = acosd(sind(p1) * sind(p2) + cosd(p1) * cosd(p2) * cosd(dlon))

        real elong_a(ni,nj)
        real alt_a(ni,nj)           ! 90. - u  (replace with emis angle?)
        real azi_a(ni,nj)
        real topo_gti(ni,nj)        ! topo normal global irradiance
        real topo_albedo(nc,ni,nj)  ! topo albedo
        real topo_solalt(ni,nj)     ! 90. - u0 (solar altitude)
        real topo_solazi(ni,nj)     ! 90. - u0 (solar altitude)
        real topo_lat(ni,nj)        ! 
        real topo_lon(ni,nj)        ! 
        real topo_lf(ni,nj)         ! topo land fraction
        real topo_sc(ni,nj)         ! topo snow/ice
        real rlat_a(ni,nj)          ! 
        real rlon_a(ni,nj)          ! 
        real Phase_angle_d(ni,nj)   ! 
        real Specular_ref_angle_d(ni,nj) 
        real emis_ang_a(ni,nj)      ! 
        real azi_fm_lnd_a(ni,nj)    ! 
        real gtic(nc,ni,nj)         ! spectral terrain GNI
        real dtic(nc,ni,nj)         ! spectral terrain diffuse NI 
        real btic(nc,ni,nj)         ! spectral terrain beam (direct) NI 
        real*8 dist_2_topo(ni,nj),sinarc     
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        integer idebug_pf(ni,nj)
        real pf_land(nc,ni,nj)      ! anisotropic reflectance factor (ARF)
                                    ! (weighted by direct/diffuse illumination)
        real cld_brdf(nc,ni,nj)

        real nonspot
        real foam(nc) /.026,.026,.026/

!       topo_gti may be too low when sun is near horizon
!       solar elev    topo_gti
!          0            12
!          1            34
!          2            56

        write(6,*)' subroutine get_lnd_pf... solazi = ',sol_azi
        iwrite = 0

        write(6,*)' sum of idebug_pf = ',sum(idebug_pf)

        write(6,*)' call satgeom ',range_m
        range_m = earth_radius + htmsl ! 42155680.00
        rlat_a = rlat
        rlon_a = rlon
        call satgeom(i4time,topo_lat,topo_lon,ni,nj,rlat_a,rlon_a & ! I
                    ,range_m,r_missing_data,Phase_angle_d &         ! I/O
                    ,Specular_ref_angle_d,emis_ang_a &              ! O
                    ,azi_fm_lnd_a,istatus)                          ! O

        write(6,11)
11      format('  i    j  ic   alt_a    azi_a  azitolnd  tsolazi azidiffg', &
               '  ampl_l    fland    fsnow   fwater   phland   phsnow   ', &
               'phwater    ph1    radfrac   dst2topo   gndarc  toposalt emis_ang  specang    alb    cldbrdf') 

        do j = 1,nj
         do i = 1,ni

            sol_clr = (ghi_zen_toa * zen_kt) * sind(max(sol_alt,1.5))
!           tfrac = topo_gti(i,j) / sol_clr                   
            azidiff = angdif(azi_a(i,j),sol_azi)
 
!           Approximate specular reflection angle
            sinarc = sind(90d0+dble(alt_a(i,j)))*dist_2_topo(i,j)/dble(earth_radius)
            gnd_arc = asind(min(sinarc,1d0))
            gnd_arc2 = gnd_arc * 2.

!           'azidiffg' is computed from a ground reference point
            if(topo_solazi(i,j) .ne. r_missing_data)then 
                azidiffg = angdif(azi_fm_lnd_a(i,j),topo_solazi(i,j))
            else
                azidiffg = 180.
            endif
            c = 180. - (gnd_arc + (90.+alt_a(i,j)))
            emis_ang = c - 90.
!           topo_salt = sol_alt + gnd_arc
            specangvert = abs(emis_ang - topo_solalt(i,j))
!           specangvert2 = specangvert * sind(emis_ang)
            alt_mean = 0.5 * (emis_ang + topo_solalt(i,j))
            azidiffsp = azidiffg / sind(max(alt_mean,.001))
            azidiffsp = min(max(azidiffsp,-180.),+180.)
            if(.false.)then
              specang = angdist(emis_ang,topo_solalt(i,j),azidiffsp)
            else
!             'specang' is calculated to be the zenith angle of the tilt required
!             for the water to reflect between the light source at altitude 
!             'topo_solalt', the observer at altitude 'emis_ang', and separated 
!             in azimuth by 'azidiffg'. Coordinates are relative to the surface
!             of the water at the point of reflection. This is a midpoint between
!             vectors approach.
              x1 = cosd(topo_solalt(i,j)) + cosd(emis_ang) * cosd(azidiffg)
              y1 =                          cosd(emis_ang) * sind(azidiffg)
              z1 = sind(topo_solalt(i,j)) + sind(emis_ang)
              specang = acosd(z1/sqrt(x1**2+y1**2+z1**2))
            endif

!           Land surface type
!           fland = scurve((1. - topo_albedo(2,i,j))**2)
!           fsnow = 1.0 - fland ! 0 is land, 1 is snow

            fsnow = (topo_albedo(2,i,j) - 0.15) / (1.0 - 0.15)
            fsnow = max(fsnow,0.)
            fsnow = topo_sc(i,j)

            if(topo_lf(i,j) .lt. 1.0)then
                fwater = 1. - topo_lf(i,j)
                fland = topo_lf(i,j)
            elseif(topo_lf(i,j) .eq. 0.0)then
                fwater = 1.0
                fland = 0.
            else ! land fraction = 1
                fwater = 0.0
                fland = 1.0
            endif

            fland  = fland  * (1.0 - fsnow)
            fwater = fwater * (1.0 - fsnow)

            fice = 0.0

            do ic = 1,nc
              tfrac = transm_obs      
              alt_thresh = 22. * ext_g(ic) 
              if(sol_alt .lt. alt_thresh .and. sol_alt .ge. 0.0)then
                tfrac = tfrac * (sol_alt/alt_thresh)**2
              endif

!             radfrac - 1 means relatively high direct / global horizontal ratio
!                       0 means zero direct / global ratio
!                       calculate from (gtic - dtic) / gtic  
              if(gtic(ic,i,j) .gt. 0. .and. dist_2_topo(i,j) .gt. 0.)then
                if(btic(ic,i,j) .le. gtic(ic,i,j))then
                  radfrac = btic(ic,i,j)/gtic(ic,i,j)         
                  iwarn = 0
                else
                  radfrac = 1.0
                  iwarn = 1
                endif
              else
                radfrac = 0.
                iwarn = 0
              endif

!             Land
!             Should look brighter opposite direct sun in low sun case
              spot = 0.005 ! fraction of energy in the spot 
              nonspot = (1. - spot)
!             First arg is max amplitude for low sun on the horizon
!             Second arg reduces amplitude for mid sun
              ampl_l = -2.0 * cosd(sol_alt)**4 * cosd(alt_a(i,j))**5 

!             Elliptical opposition effect
              aspect = 3. ! horizontal contraction of spot
              hemisphere_int = 2. / aspect
              alt_antisolar = abs(alt_a(i,j)+sol_alt) ! alt diff from antisolar point
!             azi_antisolar = 180. - (abs(azidiff)*cosd(sol_alt)) * aspect
              azi_antisolar = 180. - abs(azidiff)     ! azi diff from antisolar point
              azi_antisolar_eff = fslope(azi_antisolar*cosd(sol_alt)/180.,aspect) * 180.
              elong_antisolar = 180. - sqrt(alt_antisolar**2 + azi_antisolar_eff**2)
              elong_antisolar = max(elong_antisolar,0.)
              elong_eff = elong_antisolar * cosd(sol_alt)**2 + elong_a(i,j) * sind(sol_alt)**2

              arf_bl = nonspot * ph_exp(ampl_l,azidiff) &
!                    + spot    * hg(-.90,elong_a(i,j))               
                     + spot    * hg(-.90,elong_eff) / hemisphere_int
              arf_dl = 1.0
              phland = arf_bl * radfrac + arf_dl * (1. - radfrac)  

!             Snow
!             http://www.sciencedirect.com/science/article/pii/S002240739900028X
              ampl_s = cosd(alt_a(i,j))
              fszen = sind(topo_solalt(i,j))
              g2 = 0.45 * ampl_s 
              hg_2param = 0.25 * 1.0 + 0.75 * hg_cyl(g2,azidiff)
              brdf_szen = 0.9 + 0.1 * (2. * sind(emis_ang))
              arf_bs = hg_2param * (1.-fszen) + brdf_szen * fszen
              arf_ds = 1.0
              phsnow = arf_bs * radfrac + arf_ds * (1. - radfrac)  
              cld_brdf(ic,i,j) = 1.0 ! arf_bs

!             Water
              g = 0.6 * radfrac
              ampl_w = cosd(alt_a(i,j))
              turbidity = 1.0
              frac_aero = turbidity / (1. + turbidity)
              power_rad = 1.0 ! 0.5 * frac_aero + 1.0 * (1.-frac_aero)
              radfracw = radfrac**power_rad ! aerosols are effectively 

!             Use approximate specular reflection angle
              glint_radius = 16.0
              tsolalt_eff = max(topo_solalt(i,j),3.)
  
!             Compute ARFs relative to amount of incident light reflected 
              if(.false.)then ! original method
                glint_sr = pi * (glint_radius * rpd)**2 * sind(tsolalt_eff)
                arf_bw = 2. / glint_sr
              elseif(.true.)then ! http://www.atmos-meas-tech.net/3/813/2010/amt-3-813-2010.pdf
                sigmax = 0.5 * glint_radius * rpd
                sigmay = 0.5 * glint_radius * rpd
                pslope = 1. / (2. * pi * sigmax * sigmay) 
                cos2th = sind(tsolalt_eff) * sind(emis_ang) + cosd(tsolalt_eff) * cosd(emis_ang) * cosd(azidiffg)
                cosb = (sind(tsolalt_eff) + sind(emis_ang)) / sqrt(2. + 2. * cos2th)
!               Ratio of reflectance to Fresnel coefficient (max is 3 and should be 6?)
                arf_bw = pi * pslope / (4. * sind(tsolalt_eff) * sind(emis_ang) * cosb**4.)
              else           ! old possible method
                sun_sqdg = pi * (.533239/2.)**2
                arf_bw = ( (0.5 * sph_sqdg) / sun_sqdg) / (glint_radius / ((0.5*.533239)**2) )
                arf_bw = arf_bw / sind(tsolalt_eff)
              endif
           
!             argexp = min(specang/glint_radius,8.)
              argexp = min(specang/(glint_radius/2.),8.) ! for new specang
              specamp = exp(-(argexp)**2)

!             phwater = 0.2 + 3.0 * specamp ! ampl_w * hg(g,azidiff) / (1. + 1.3*g**2)

!             Add sunlight/skylight reflecting off the water surface
!             arg_glint = opac(phwater * fwater)
              phi_b = (180. - elong_a(i,j)) * rpd * 0.5
!             alb_glint = .045 / sind(alt_eff)**.75
              alb_glint = fresnel(1.33,phi_b)
              fresnel1 = fresnel(1.33,phi_b)
!             fresnel2 = foam(ic) * (1. - fresnel1) + (1.-foam(ic)) * fresnel1
              fresnel2 = foam(ic) * 0.5             + (1.-foam(ic)) * fresnel1 ! black sky albedo

              phi_d = (90. - emis_ang) * rpd
              fresnel_d = fresnel(1.33,phi_d)                                  ! white sky albedo
              arf_dw = 1.0 ! fresnel_d / .08

!             smooth_water_albedo = radfracw * fresnel2 + (1. - radfracw) * fresnel_d
              smooth_water_albedo = .08
              water_albedo = foam(ic) * 0.5 + (1.-foam(ic)) * smooth_water_albedo
              fresnel_mean = fresnel2 * radfracw + fresnel_d * (1.-radfracw)     ! "gray" sky albedo
              fresnel_arf = fresnel_mean / smooth_water_albedo

!             Override geog albedo values over water, based on land fraction
!             topo_albedo(ic,i,j) = smooth_water_albedo * fwater + topo_albedo(ic,i,j) * (1. - fwater)

!             phwater = arf_bw * specamp
!             phwater = arf_bw * specamp + f(radfracw,arf_dw)
              phwater = (arf_bw * specamp * radfracw + arf_dw * (1. - radfracw)) * fresnel_arf

!             Ice
              phice = 1.0

!             Check for valid scenario
              if(alt_a(i,j) .gt. 0.)then
                phwater = 1. ! default value
              endif

              if(topo_solalt(i,j) .ge. 0.)then ! light source above horizon or land normal
                ph1 = phland * fland + phsnow * fsnow + phwater * fwater
              else
                ph1 = 1.     ! default value
                phwater = 1. ! default value
              endif

!             if((i .eq. ni-100 .and. j .eq. (j/40)*40) .OR.  &
              if(ic .eq. 2)then

                call check_nan(ph1,istat_nan)
                if(istat_nan .ne. 1)then
                  write(6,*)' ERROR: NaN in pf_land ',i,j,sinarc,phi_b/rpd,phi_d/rpd
                endif

!               if((i .eq. 64 .and. j .eq. (j/40)*40) .OR. (ph1 .lt. 0. .and. dist_2_topo(i,j) .gt. 0.) .OR.&
!                ( (abs(azidiff) .lt. azi_scale/2. .or. abs(azidiff) .gt. (180.-azi_scale/2.)) &
!                         .and. i .eq. (i/5)*5 .and. alt_a(i,j) .lt. 5.) .OR. &
!               if( (i .eq. (i/40)*40 .AND. j .eq. nj/2) .OR. &
                if( (idebug_pf(i,j) .eq. 1) .OR. &
                       (alt_a(i,j) .eq. -90. .and. j .eq. 2161) .OR. (istat_nan .ne. 1) )then ! nadir
                  write(6,1)i,j,ic,alt_a(i,j),azi_a(i,j),azi_fm_lnd_a(i,j)+180.,topo_solazi(i,j),azidiffg,ampl_l,fland,fsnow,fwater,phland,phsnow,phwater,ph1,radfrac,dist_2_topo(i,j),gnd_arc,topo_solalt(i,j),emis_ang,specang,topo_albedo(2,i,j),cld_brdf(2,i,j)
1                 format(/i4,i5,i2,f9.4,4f9.2,9f9.4,f12.0,4f9.2,2f9.3)
                  write(6,111)alt_antisolar,azi_antisolar,azi_antisolar_eff,elong_antisolar,elong_a(i,j),elong_eff
111               format('   antisolar alt/azi/azeff/elg/elg_a/eff',6f9.2)
!                 if(abs(emis_ang - emis_ang_a(i,j)) .gt. 0.1 .AND. emis_ang_a(i,j) .gt. 0.)then
                  if(emis_ang .gt. 0.)then
                    write(6,112)emis_ang,emis_ang_a(i,j),azi_fm_lnd_a(i,j),topo_lat(i,j),topo_lon(i,j),topo_solalt(i,j),topo_solazi(i,j),topo_lf(i,j),topo_sc(i,j)
112                 format('   EMISANG INFO:',9f9.3)
                  endif
                  if(fsnow .gt. 0.01)then
                    write(6,113)ampl_s,fszen,g2,hg_2param,brdf_szen,arf_bs 
113                 format('   SNOW BRDF:   ',6f9.3)
                  elseif(fwater .gt. 0.01)then
                    write(6,114)fresnel_mean,fresnel_arf,arf_bw,arf_dw,radfracw,specamp,phwater
114                 format('   WATER BRDF:  ',7f9.3)
                  endif
                endif
              endif

              pf_land(ic,i,j) = ph1

            enddo ! ic

            if(idebug_pf(i,j) .eq. 1 .and. fwater .gt. 0.5)then
                if(iwarn .eq. 1)write(6,*)' alt_a = ',i,j,alt_a(i,j)
                write(6,2)elong_a(i,j),topo_gti(i,j),sol_alt,transm_obs,radfrac,fland,fsnow,fwater,spot
2               format('   elg/tgti/solalt/trnm/radf/fland/fsnow/fwater/fspot',f9.3,8f9.4)
                write(6,3)gtic(:,i,j),dtic(:,i,j),btic(:,i,j),iwarn
!               format('   gtic',3f11.6,'  dtic',3f11.6,'  btic',3f11.6,i3)
3               format('   gtic',3e12.4,'  dtic',3e12.4,'  btic',3e12.4,i3)
                iwrite = iwrite + 1
            endif

            emis_ang_a(i,j) = emis_ang ! more accurate than from 'sat_geom'

         enddo ! i (altitude)

        enddo ! j (azimuth)

        write(6,*)' get_lnd_pf: range of lf = ',minval(topo_lf(:,:)) &
                                               ,maxval(topo_lf(:,:))

        write(6,*)' get_lnd_pf: range of sc = ',minval(topo_sc(:,:)) &
                                               ,maxval(topo_sc(:,:))

        write(6,*)' get_lnd_pf: range of pf = ',minval(pf_land(2,:,:)) &
                                               ,maxval(pf_land(2,:,:))

        return
        end


        real function Fresnel(nWat,phi)

        implicit none
        real nWat,phi, aRef,aDif,aSum,Rpar,Rper,foam,phi90

        parameter (phi90 = ASIN(1.))

        foam = 0. ! .026
        
!       Phi is incidence angle in radians
        if (phi .gt. phi90) then ! exceeds 90 degrees       
           Fresnel = 1.
        elseif (phi .ne. 0.) then
           aRef = ASIN(SIN(phi)/nWat)
           aDif = phi-aRef
           aSum = phi+aRef
           Rpar = TAN(aDif)/TAN(aSum)
           Rper = SIN(aDif)/SIN(aSum)
           Fresnel = 0.5*(Rpar*Rpar+Rper*Rper)
        else
           aSum = (nwat-1.)/(nWat+1.)
           Fresnel = aSum*aSum
        end if

        fresnel = foam * (1. - fresnel) + (1.-foam) * fresnel

        return
        end 


