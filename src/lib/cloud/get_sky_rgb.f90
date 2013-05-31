

        subroutine get_sky_rgb(r_cloud_3d,r_cloud_rad,glow, &
                               airmass_2_cloud,airmass_2_topo, & 
                               alt_a,azi_a,elong_a,ni,nj,sol_alt,sol_az, &
                               sky_rgb)                                     ! O

        use mem_namelist, ONLY: r_missing_data
        include 'trigd.inc'

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real glow(ni,nj)            ! skyglow
        real airmass_2_cloud(ni,nj) ! airmass to cloud 
        real airmass_2_topo(ni,nj)  ! airmass to topo  
        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real elong_a(ni,nj)

        real sky_rgb(0:2,ni,nj)

        write(6,*)' get_sky_rgb'

        if(ni .eq. nj)then ! polar
            write(6,*)' slice from SW to NE through midpoint'
        else
            write(6,*)' slices at 45,225 degrees azimuth'
        endif

        if(ni .eq. nj)write(6,11)
11      format('    i   j      alt      azi     elong    opac     cloud  airmass  rad     rinten   airtopo  cld_visb       skyrgb')

        itopo_red = 150; itopo_grn = 150 ; itopo_blu = 0

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
          whiteness_thk = r_cloud_rad(i,j) ! * rill  ! default is 0.
!         rint_coeff =  380. * (1. - whiteness_thk) ! default is 380.
          rint_coeff =  900. * (1. - whiteness_thk) ! default is 380.

          rintensity = 250. - ((abs(r_cloud_3d(i,j)-0.6)**2.0) * rint_coeff)

          rintensity = min(rintensity,255.)
          rintensity = max(rintensity,0.)

          grayness = r_cloud_3d(i,j) ** 0.4
          grayness = 0.4 + (0.6 * grayness)
          grayness = 1.0                         

          if(sol_alt .le. 3.0 .and. sol_alt .gt. 0.)then
              redness = (3.0 - sol_alt) / 3.0
          else
              redness = 0.
          endif

          cld_red = nint(rintensity * grayness)                
          cld_grn = nint(rintensity * grayness * (1. - redness)**0.3)                   
          cld_blu = nint(rintensity            * (1. - redness))

          rintensity_glow = min(((glow(i,j)-7.) * 100.),255.)
          clr_red = rintensity_glow * rintensity_glow / 255.
          clr_grn = rintensity_glow * rintensity_glow / 255.
          clr_blu = rintensity_glow  

          if(elong_a(i,j) .le. 8.0)then
              red_elong = (8.0 - sol_alt) / 8.0
              clr_red = clr_red * 1.0
              clr_grn = clr_grn * (1. - redness * red_elong)**0.3 
              clr_blu = clr_blu * (1. - redness * red_elong)
          endif

!         cloud_visibility = exp(-0.28*airmass_2_cloud(i,j))
          cloud_visibility = exp(-0.14*airmass_2_cloud(i,j))

!         Use clear sky values if cloud cover is less than 0.5
          frac_cloud = r_cloud_3d(i,j)
!         frac_cloud = 0.0 ; Testing
          frac_cloud = frac_cloud * cloud_visibility  
          frac_clr = 1.0 - frac_cloud                         
          sky_rgb(0,I,J) = clr_red * frac_clr + cld_red * frac_cloud
          sky_rgb(1,I,J) = clr_grn * frac_clr + cld_grn * frac_cloud
          sky_rgb(2,I,J) = clr_blu * frac_clr + cld_blu * frac_cloud

!         Use topo value if airmass to topo > 0
          if(airmass_2_topo(i,j) .gt. 0.)then
              if(airmass_2_cloud(i,j) .gt. 0. .AND. airmass_2_cloud(i,j) .le. airmass_2_topo(i,j)) then
                  continue
              else
                  sky_rgb(0,I,J) = itopo_red
                  sky_rgb(1,I,J) = itopo_grn
                  sky_rgb(2,I,J) = itopo_blu
              endif
          endif

          if(idebug .eq. 1)then
              write(6,101)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j),r_cloud_3d(i,j) &
                         ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity,airmass_2_topo(i,j) &
                         ,cloud_visibility,nint(sky_rgb(:,i,j))
101           format(2i5,3f9.2,3f9.3,f7.4,f9.1,2f9.3,2x,3i6)
          endif

         endif ! missing data tested via altitude

        enddo ! i
        enddo ! j

        return
        end
