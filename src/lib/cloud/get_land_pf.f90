
        subroutine get_lnd_pf(elong_a,alt_a,azi_a,topo_swi,topo_albedo,cloud_od,cloud_od_sp,transm_obs,sol_alt,sol_azi,nsp,airmass_2_topo,idebug_a,ni,nj & ! I
                             ,pf_land) ! O
        include 'trigd.inc'

!       Statement functions
        trans(od) = exp(-od)
        opac(od) = 1.0 - exp(-od)
        alb(bt) = bt / (1.+bt)
        rad2tau(b,r) = (1.-r)/(r*b)
        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.

        include 'rad.inc'
        real elong_a(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real cloud_od(ni,nj)        ! cloud optical depth (tau)
        real cloud_od_sp(ni,nj,nsp) ! cloud species tau (clwc,cice,rain,snow)
        real topo_swi(ni,nj)        ! topo normal global irradiance
        real topo_albedo(nc,ni,nj)  ! topo albedo
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
11      format('  i    j     alt_a   azi_a   sol_azi  azidiff  ', &
               '  ampl   snw_fac   phland   phsnow    ph1    radfrac     alb') 

        do j = 1,nj
         do i = 1,ni

            sol_clr = 1109.46 * sind(max(sol_alt,1.5))
!           tfrac = topo_swi(i,j) / sol_clr                   
            tfrac = transm_obs      
            if(sol_alt .lt. 2.0 .and. sol_alt .ge. 0.0)then
                tfrac = tfrac * (sol_alt/2.0)**2
            endif
            radfrac = scurve(tfrac**3) ! high for illuminated clouds
!                     illuminated                unilluminated

            alb_term = scurve((1. - topo_albedo(2,i,j))**2)
            snow_factor = 1.0 - alb_term ! 0 is land, 1 is snow

            spot = radfrac * alb_term
            arg2 = spot * 0.010 ; arg1 = (1. - arg2)
            azidiff = angdif(azi_a(i,j),sol_azi)
            ampl = 0.5 * cosd(sol_alt)**5 * cosd(alt_a(i,j))**5 * radfrac

            phland = exp(ampl * cosd(azidiff))

!           Should approximately integrate to 1 over the "cylinder"
            g = 0.3 * radfrac
            phsnow = hg(g,azidiff) / (1. + 1.3*g**2)

            ph1 = phland * (1.-snow_factor) + phsnow * snow_factor

            if((i .eq. ni-100 .and. j .eq. (j/40)*40) .OR.  &
               (abs(azidiff) .lt. 0.10 .and. i .eq. (i/10)*10 .and. alt_a(i,j) .lt. 0.) )then
                write(6,1)i,j,alt_a(i,j),azi_a(i,j),sol_azi,azidiff,ampl,snow_factor,phland,phsnow,ph1,radfrac,topo_albedo(2,i,j)
1               format(i4,i5,4f9.2,6f9.4,f9.2)
            endif

            if(elong_a(i,j) .gt. 179.8 .and. iwrite .eq. 0)then
                write(6,2)elong_a(i,j),topo_swi(i,j),sol_alt,transm_obs,radfrac,alb_term,spot
2               format(' elg/tswi/solalt/trnm/radf/albt/spot/ph-ls1',f9.3,6f9.4)
                iwrite = iwrite + 1
            endif

            do ic = 1,nc
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

