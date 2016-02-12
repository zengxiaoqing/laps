
        subroutine get_lnd_pf(elong_a,alt_a,azi_a &                     ! I
                             ,topo_gti,topo_albedo,transm_obs &         ! I
                             ,gtic,dtic,btic &                          ! I
                             ,dist_2_topo,topo_solalt,azi_scale &       ! I
                             ,sol_alt,sol_azi,nsp,airmass_2_topo,idebug_a,ni,nj & ! I
                             ,pf_land) ! O

        use mem_namelist, ONLY: r_missing_data,earth_radius
        use cloud_rad, ONLY: ghi_zen_toa, zen_kt
        include 'trigd.inc'
        include 'rad.inc'

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

        real elong_a(ni,nj)
        real alt_a(ni,nj)           ! 90. - u
        real azi_a(ni,nj)
        real topo_gti(ni,nj)        ! topo normal global irradiance
        real topo_albedo(nc,ni,nj)  ! topo albedo
        real topo_solalt(ni,nj)     ! 90. - u0 (solar altitude)
        real gtic(nc,ni,nj)         ! spectral terrain GNI
        real dtic(nc,ni,nj)         ! spectral terrain diffuse NI 
        real btic(nc,ni,nj)         ! spectral terrain beam (direct) NI 
        real dist_2_topo(ni,nj)     
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        integer idebug_a(ni,nj)
        real pf_land(nc,ni,nj)      ! anisotropic reflectance factor (ARF)
                                    ! (weighted by direct/diffuse illumination)

        real nonspot

!       topo_gti may be too low when sun is near horizon
!       solar elev    topo_gti
!          0            12
!          1            34
!          2            56

        write(6,*)' subroutine get_lnd_pf...'
        iwrite = 0

        write(6,11)
11      format('  i    j  ic   alt_a    azi_a   sol_azi azidiff  ', &
               ' ampl_l    fland    fsnow   fwater   phland   phsnow   phwater    ph1    radfrac  dst2topo gndarc  toposalt emis_ang  specang    alb') 

        do j = 1,nj
         do i = 1,ni

            sol_clr = (ghi_zen_toa * zen_kt) * sind(max(sol_alt,1.5))
!           tfrac = topo_gti(i,j) / sol_clr                   
            azidiff = angdif(azi_a(i,j),sol_azi)
 
!           Approximate specular reflection angle
            gnd_arc = asind(sind(90.+alt_a(i,j))*dist_2_topo(i,j)/earth_radius)
            c = 180. - (gnd_arc + (90.+alt_a(i,j)))
            emis_ang = c - 90.
!           topo_salt = sol_alt + gnd_arc
            specangvert = abs(emis_ang - topo_solalt(i,j))
            specangvert2 = specangvert * sind(emis_ang)
            azidiff2 = azidiff * cosd(emis_ang)
            specang = sqrt(specangvert**2 + azidiff2**2)

            do ic = 1,nc
              tfrac = transm_obs      
              alt_thresh = 22. * ext_g(ic) 
              if(sol_alt .lt. alt_thresh .and. sol_alt .ge. 0.0)then
                tfrac = tfrac * (sol_alt/alt_thresh)**2
              endif

!             radfrac - 1 means relatively high direct / global horizontal ratio
!                       0 means zero direct / global ratio
!                       calculate from (gtic - dtic) / gtic  
              if(gtic(ic,i,j) .gt. 0.)then
                if(btic(ic,i,j) .le. gtic(ic,i,j))then
                  radfrac = btic(ic,i,j)/gtic(ic,i,j)         
                  iwarn = 0
                else
                  radfrac = 1.0
                  iwarn = 1
                endif
              else
                radfrac = 0.
              endif

!             fland = scurve((1. - topo_albedo(2,i,j))**2)
!             fsnow = 1.0 - fland ! 0 is land, 1 is snow

              fsnow = (topo_albedo(2,i,j) - 0.15) / (1.0 - 0.15)
              fsnow = max(fsnow,0.)
              fland = 1.0 - fsnow

              if(topo_albedo(1,i,j) .eq. .08 .and. topo_albedo(2,i,j) .eq. .08 .and. topo_albedo(3,i,j) .eq. .08)then
                  fwater = 1.0
                  fland = 0.
                  fsnow = 0.
              else
                  fwater = 0.0
              endif

              fice = 0.0

!             Land
!             Should look brighter opposite direct sun in low sun case
              spot = 0.005 ! fraction of energy in the spot 
              nonspot = (1. - spot)
              ampl_l = -1.0 * cosd(sol_alt)**2 * cosd(alt_a(i,j))**5 

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

              arf_b = nonspot * ph_exp(ampl_l,azidiff) &
!                   + spot    * hg(-.90,elong_a(i,j))               
                    + spot    * hg(-.90,elong_eff) / hemisphere_int
              arf_d = 1.0
              phland = arf_b * radfrac + arf_d * (1. - radfrac)  

!             Snow
              ampl_s = cosd(alt_a(i,j))
              fszen = sind(topo_solalt(i,j))
              g2 = 0.45 * ampl_s ! previously 0.50 * ampl_s
              hg_2param = 0.25 * 1.0 + 0.75 * hg_cyl(g2,azidiff)
              brdf_szen = 0.7 + 0.3 * (2. * sind(alt_a(i,j)))
              arf_b = hg_2param * (1.-fszen) + brdf_szen * fszen
              arf_d = 1.0
              phsnow = arf_b * radfrac + arf_d * (1. - radfrac)  

!             Water
              g = 0.6 * radfrac
              ampl_w = cosd(alt_a(i,j))

!             Use approximate specular reflection angle
              argexp = min(specang/5.,8.)
              specamp = exp(-(argexp)**2)

              phwater = 0.2 + 3.0 * specamp ! ampl_w * hg(g,azidiff) / (1. + 1.3*g**2)

!             Ice
              phice = 1.0

              ph1 = phland * fland + phsnow * fsnow + phwater * fwater

!             if((i .eq. ni-100 .and. j .eq. (j/40)*40) .OR.  &
              if(ic .eq. 2)then
                if((i .eq. 1 .and. j .eq. (j/40)*40) .OR. ph1 .lt. 0. .OR.&
                 ( (abs(azidiff) .lt. azi_scale/2. .or. abs(azidiff) .gt. (180.-azi_scale/2.)) &
                          .and. i .eq. (i/5)*5 .and. alt_a(i,j) .lt. 5.) )then
                  write(6,1)i,j,ic,alt_a(i,j),azi_a(i,j),sol_azi,azidiff,ampl_l,fland,fsnow,fwater,phland,phsnow,phwater,ph1,radfrac,dist_2_topo(i,j),gnd_arc,topo_solalt(i,j),emis_ang,specang,topo_albedo(2,i,j)
1                 format(i4,i5,i2,4f9.2,9f9.4,f9.0,5f9.2)
                  write(6,111)alt_antisolar,azi_antisolar,azi_antisolar_eff,elong_antisolar,elong_a(i,j),elong_eff
111               format('   antisolar alt/azi/azeff/elg/elg_a/eff',6f9.2)
                endif
              endif

              if((elong_a(i,j) .gt. 179.8 .and. iwrite .eq. 0) .OR. iwarn .eq. 1)then
                write(6,2)elong_a(i,j),topo_gti(i,j),sol_alt,transm_obs,radfrac,fland,fsnow,fwater,spot
2               format(' elg/tgti/solalt/trnm/radf/fland/fsnow/fwater/fspot',f9.3,8f9.4)
                write(6,3)gtic(:,i,j),dtic(:,i,j),btic(:,i,j),iwarn
!               format(' gtic',3f11.6,'  dtic',3f11.6,'  btic',3f11.6,i3)
3               format(' gtic',3e12.4,'  dtic',3e12.4,'  btic',3e12.4,i3)
                iwrite = iwrite + 1
              endif

              pf_land(ic,i,j) = ph1

            enddo ! ic

         enddo ! i (altitude)

        enddo ! j (azimuth)

        write(6,*)' get_lnd_pf: range of pf = ',minval(pf_land(2,:,:)) &
                                               ,maxval(pf_land(2,:,:))

        return
        end

