

        subroutine get_sky_rgb(r_cloud_3d,r_cloud_rad,cloud_rad_c, &
                               clear_rad_c, &
                               glow,od_atm_a, &
                               airmass_2_cloud,airmass_2_topo, &
                               topo_swi,topo_albedo, & 
                               alt_a,azi_a,elong_a,ni,nj,sol_alt,sol_az, &
                               sky_rgb)                                     ! O

        use mem_namelist, ONLY: r_missing_data
        include 'trigd.inc'

        parameter (nc = 3)

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,ni,nj)  ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,ni,nj)  ! integrated fraction of air illuminated by the sun along line of sight                        
        real glow(ni,nj)            ! skyglow
        real airmass_2_cloud(ni,nj) ! airmass to cloud 
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        real topo_swi(ni,nj)        ! terrain illumination
        real topo_albedo(nc,ni,nj)  ! terrain albedo
        real aod_2_cloud(ni,nj)     ! future use
        real aod_2_topo(ni,nj)      ! future use
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real elong_a(ni,nj)
        real rintensity(nc)

        real sky_rgb(0:2,ni,nj)

        write(6,*)' get_sky_rgb: sol_alt = ',sol_alt

        if(ni .eq. nj)then ! polar
            write(6,*)' slice from SW to NE through midpoint'
        else
            write(6,*)' slices at 45,225 degrees azimuth'
        endif

        if(ni .eq. nj)write(6,11)
11      format('    i   j      alt      azi     elong   pf_scat   opac       od      alb     cloud  airmass   rad    rinten   airtopo  switopo  topoalb  topovis   cld_visb  glow          skyrgb')

        do j = 1,nj
        do i = 1,ni

         if(alt_a(i,j) .eq. r_missing_data)then
          sky_rgb(:,i,j) = 0.          
         else
          idebug = 0
          if(ni .eq. nj)then ! polar
!             if(i .eq. ni/2 .AND. j .eq. (j/5) * 5)then
              if(i .eq.    j .AND. j .eq. (j/5) * 5)then ! SW/NE
                  idebug = 1
              endif
          else ! cyl
              if(j .eq. 226 .or. j .eq. 46)then ! constant azimuth
                  idebug = 1
                  if(i .eq. 1)write(6,11)   
              endif
          endif

!         a substitute for cloud_rad could be arg2 = (cosd(alt))**3.

!         using 'rill' is a substitute for considering the slant path
!         in 'get_cloud_rad'
          rill = (1. - cosd(elong_a(i,j)))/2.
          whiteness_thk = r_cloud_rad(i,j) ! * rill ! default is 0.
!         rint_coeff =  380. * (1. - whiteness_thk) ! default is 380.
!         rint_coeff =  900. * (1. - whiteness_thk) ! default is 380.

!         rintensity = 250. - ((abs(r_cloud_3d(i,j)-0.6)**2.0) * rint_coeff)

!         Phase function that depends on degree of forward scattering in cloud    
!         pwr controls angular extent of forward scattering peak of a thin cloud
!         Useful reference: http://www-evasion.imag.fr/Publications/2008/BNMBC08/clouds.pdf
          if(elong_a(i,j) .le. 90.)then
              pwr = 3.0
              ampl = r_cloud_rad(i,j) * 0.7; b = 1.0 + pwr * r_cloud_rad(i,j)
              pf_scat = 0.9 + ampl * (cosd(min(elong_a(i,j),89.99))**b)
              cloud_od   = -99.9 ! flag value for logging
              bkscat_alb = -99.9 ! flag value for logging
          else
!             convert from opacity to albedo
!             bkscat_alb = r_cloud_3d(i,j) ** 2.0 ! approx opacity to albedo
              cloud_opacity = min(r_cloud_3d(i,j),0.999999)
              cloud_od = -log(1.0 - cloud_opacity)
              bksc_eff_od = cloud_od * 0.12 ! > .10 due to machine epsilon
              cloud_rad_trans = exp(-bksc_eff_od)
              bkscat_alb = 1.0 - cloud_rad_trans 
              ampl = 0.15 * bkscat_alb
              pf_scat = 0.9 + ampl * (-cosd(elong_a(i,j)))
          endif

!         Potential intensity of cloud if it is opaque 
!           (240. is nominal intensity of a white cloud far from the sun)
!           (0.25 is dark cloud base value)                                  
          rintensity(:) = 240. * ( (0.25 + 0.75 * cloud_rad_c(:,i,j)) * pf_scat)

!         rintensity = min(rintensity,255.)
          rintensity = max(rintensity,0.)

          cld_red = nint(rintensity(1))                
          cld_grn = nint(rintensity(2))                        
          cld_blu = nint(rintensity(3))                         

          if(sol_alt .ge. 0.)then ! Daylight from skyglow routine
              rintensity_glow = min(((glow(i,j)-7.) * 100.),255.)
              clr_red = rintensity_glow * rintensity_glow / 255.
              clr_grn = rintensity_glow * rintensity_glow / 255.
              clr_blu = rintensity_glow  
          else ! Twilight / Night from clear_rad_c array
!             rintensity_glow = 200. * clear_rad_c(1,i,j)**0.3
!             rintensity_glow = min(((8.0      -7.) * 100.),255.)
              hue = clear_rad_c(1,i,j)
              sat = clear_rad_c(2,i,j)
              arg = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15             
              rintensity_glow = max(min(((arg     -7.) * 100.),255.),70.)
              call hsi_to_rgb(hue,sat,rintensity_glow,clr_red,clr_grn,clr_blu)
!             clr_red = rintensity_glow * clear_rad_c(2,i,j)
!             clr_blu = rintensity_glow * clear_rad_c(3,i,j)
!             clr_grn = 0.5 * (clr_red + clr_blu)                               
          endif

          sol_alt_red_thr = 1.0 + (od_atm_a * 40.)

          if(sol_alt .le. sol_alt_red_thr)then
              redness = min((sol_alt_red_thr - sol_alt) / sol_alt_red_thr,1.0)
          else
              redness = 0.
          endif

          elong_red = 12.0
          if(elong_a(i,j) .le. elong_red)then
!         if(.false.)then                               
              red_elong = (elong_red - elong_a(i,j)) / elong_red
              clr_red = clr_red * 1.0
              clr_grn = clr_grn * (1. - redness * red_elong)**0.3 
              clr_blu = clr_blu * (1. - redness * red_elong)
          endif

          od_atm_g = 0.28
          od_2_cloud = (od_atm_g + od_atm_a) * airmass_2_cloud(i,j)

!         cloud_visibility = exp(-0.28*airmass_2_cloud(i,j))
!         cloud_visibility = exp(-0.14*airmass_2_cloud(i,j))
          cloud_visibility = exp(-0.43*od_2_cloud)

!         Use clear sky values if cloud cover is less than 0.5
          frac_cloud = r_cloud_3d(i,j)
!         frac_cloud = 0.0 ; Testing
          frac_cloud = frac_cloud * cloud_visibility  
          frac_clr = 1.0 - frac_cloud                         
          sky_rgb(0,I,J) = clr_red * frac_clr + cld_red * frac_cloud
          sky_rgb(1,I,J) = clr_grn * frac_clr + cld_grn * frac_cloud
          sky_rgb(2,I,J) = clr_blu * frac_clr + cld_blu * frac_cloud
          sky_rgb(:,I,J) = min(sky_rgb(:,I,J),255.)

!         Use topo value if airmass to topo > 0
          if(airmass_2_topo(i,j) .gt. 0.)then
              od_2_topo = (od_atm_g + od_atm_a) * airmass_2_topo(i,j)
              topo_visibility = exp(-0.43*od_2_topo)                    
              if(airmass_2_cloud(i,j) .gt. 0. .AND. airmass_2_cloud(i,j) .lt. airmass_2_topo(i,j)) then
                  continue
              else ! note 0.45 is display gamma correction  
                  topo_swi_frac = (max(topo_swi(i,j),001.) / 1000.) ** 0.45
                  rtopo_red = 150. * topo_swi_frac * (topo_albedo(1,i,j)/.15)**0.45
                  rtopo_grn = 150. * topo_swi_frac * (topo_albedo(2,i,j)/.15)**0.45
                  rtopo_blu = 150. * topo_swi_frac * (topo_albedo(3,i,j)/.15)**0.45

                  sky_rgb(0,I,J) = nint(rtopo_red*topo_visibility + sky_rgb(0,I,J)*(1.0-topo_visibility) )
                  sky_rgb(1,I,J) = nint(rtopo_grn*topo_visibility + sky_rgb(1,I,J)*(1.0-topo_visibility) )
                  sky_rgb(2,I,J) = nint(rtopo_blu*topo_visibility + sky_rgb(2,I,J)*(1.0-topo_visibility) )
              endif
          endif

          if(idebug .eq. 1)then
              if(sol_alt .ge. 0.)then
                  write(6,102)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j),pf_scat,r_cloud_3d(i,j) &
                         ,cloud_od,bkscat_alb &
                         ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                         ,topo_swi(i,j),topo_albedo(1,i,j),topo_visibility,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j))
              else  
                  write(6,103)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j),pf_scat,r_cloud_3d(i,j) &
                         ,cloud_od,bkscat_alb &
                         ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                         ,topo_swi(i,j),topo_albedo(1,i,j),topo_visibility,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c(:,i,j)
              endif
102           format(2i5,3f9.2,6f9.3,f7.4,f9.1,f9.3,f9.1,3f9.3,f9.2,2x,3i6)
103           format(2i5,3f9.2,6f9.3,f7.4,f9.1,f9.3,f9.1,3f9.3,f9.2,2x,3i6,' clrrad',3f10.6)
          endif

         endif ! missing data tested via altitude

        enddo ! i
        enddo ! j

        return
        end
