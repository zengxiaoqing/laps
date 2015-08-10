
        subroutine get_lnd_pf(elong_a,alt_a,azi_a &                     ! I
                             ,topo_swi,topo_albedo,transm_obs &         ! I
                             ,dist_2_topo,topo_solalt &                 ! I
                             ,sol_alt,sol_azi,nsp,airmass_2_topo,idebug_a,ni,nj & ! I
                             ,pf_land) ! O

        use mem_namelist, ONLY: r_missing_data,earth_radius
        include 'trigd.inc'

!       Statement functions
        trans(od) = exp(-min(od,80.))
        opac(od) = 1.0 - trans(od)
        alb(bt) = bt / (1.+bt)
        rad2tau(b,r) = (1.-r)/(r*b)
        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.

        include 'rad.inc'
        real elong_a(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real topo_swi(ni,nj)        ! topo normal global irradiance
        real topo_albedo(nc,ni,nj)  ! topo albedo
        real topo_solalt(ni,nj)     ! solar altitude
        real dist_2_topo(ni,nj)     
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        integer idebug_a(ni,nj)
        real pf_land(nc,ni,nj)

!       topo_swi may be too low when sun is near horizon
!       solar elev    topo_swi
!          0            12
!          1            34
!          2            56

        scurve(x) = (-0.5 * cos(x*3.14159265)) + 0.5  ! range of x/scurve is 0 to 1

        write(6,*)' subroutine get_lnd_pf...'
        iwrite = 0

        write(6,11)
11      format('  i    j  ic   alt_a    azi_a   sol_azi azidiff  ', &
               '  ampl     fland    fsnow   fwater   phland   phsnow   phwater    ph1    radfrac  dst2topo gndarc  toposalt  specang    alb') 

        do j = 1,nj
         do i = 1,ni

            sol_clr = 1109.46 * sind(max(sol_alt,1.5))
!           tfrac = topo_swi(i,j) / sol_clr                   
            do ic = 1,nc
              tfrac = transm_obs      
              alt_thresh = 22. * ext_g(ic) 
              if(sol_alt .lt. alt_thresh .and. sol_alt .ge. 0.0)then
                tfrac = tfrac * (sol_alt/alt_thresh)**2
              endif
              radfrac = scurve(tfrac**3) ! high for illuminated clouds
!                       illuminated                unilluminated

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

!             Land
              spot = radfrac * fland
              arg2 = spot * 0.010 ; arg1 = (1. - arg2)
              azidiff = angdif(azi_a(i,j),sol_azi)
              ampl = -0.7 * cosd(sol_alt)**2 * cosd(alt_a(i,j))**5 * radfrac

              phland = exp(ampl * cosd(azidiff))
!             phland = 1.0 + (ampl * cosd(azidiff))

!             Snow
!             Should approximately integrate to 1 over the "cylinder"
              g = 0.3 * radfrac
              ampl_s = cosd(alt_a(i,j))
              phsnow = ampl_s * hg(g,azidiff) / (1. + 1.3*g**2)

!             Water
              g = 0.6 * radfrac
              ampl_w = cosd(alt_a(i,j))

              gnd_arc = asind(sind(90.+alt_a(i,j))*dist_2_topo(i,j)/earth_radius)
              c = 180. - (gnd_arc + (90.+alt_a(i,j)))
              emis_ang = c - 90.
!             topo_salt = sol_alt + gnd_arc
              specangvert = abs(emis_ang - topo_solalt(i,j))
              specangvert2 = specangvert * sind(emis_ang)
              azidiff2 = azidiff * cosd(emis_ang)
              specang = sqrt(specangvert2**2 + azidiff2**2)
              argexp = min(specang/5.,8.)
              specamp = exp(-(argexp)**2)

              phwater = 0.2 + 3.0 * specamp ! ampl_w * hg(g,azidiff) / (1. + 1.3*g**2)

              ph1 = phland * fland + phsnow * fsnow + phwater * fwater

!             if((i .eq. ni-100 .and. j .eq. (j/40)*40) .OR.  &
              if(ic .eq. 2)then
                if((i .eq. 1 .and. j .eq. (j/40)*40) .OR.  &
                 (abs(azidiff) .lt. 0.10 .and. i .eq. (i/5)*5 .and. alt_a(i,j) .lt. 5.) )then
                  write(6,1)i,j,ic,alt_a(i,j),azi_a(i,j),sol_azi,azidiff,ampl,fland,fsnow,fwater,phland,phsnow,phwater,ph1,radfrac,dist_2_topo(i,j),gnd_arc,topo_solalt(i,j),specang,topo_albedo(2,i,j)
1                 format(i4,i5,i2,4f9.2,9f9.4,f9.0,4f9.2)
                endif
              endif

              if(elong_a(i,j) .gt. 179.8 .and. iwrite .eq. 0)then
                write(6,2)elong_a(i,j),topo_swi(i,j),sol_alt,transm_obs,radfrac,fland,fsnow,fwater,spot
2               format(' elg/tswi/solalt/trnm/radf/albt/spot/ph-lsw1',f9.3,7f9.4)
                iwrite = iwrite + 1
              endif

              pf_land(ic,i,j) & 
                      = arg1 * ph1 & 
                      + arg2 * hg(-.90,elong_a(i,j)) 

            enddo ! ic

         enddo ! i (altitude)

        enddo ! j (azimuth)

        write(6,*)' get_lnd_pf: range of pf = ',minval(pf_land(2,:,:)) &
                                               ,maxval(pf_land(2,:,:))

        return
        end

