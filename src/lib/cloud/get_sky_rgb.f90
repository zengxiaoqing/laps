

        subroutine get_sky_rgb(r_cloud_3d,glow,airmass_2_cloud & 
                              ,alt_a,azi_a,ni,nj,sky_rgb)   

        use mem_namelist, ONLY: r_missing_data

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real glow(ni,nj)            ! skyglow
        real airmass_2_cloud(ni,nj) ! airmass to cloud 
        real alt_a(ni,nj)
        real azi_a(ni,nj)

        real sky_rgb(0:2,ni,nj)

        write(6,*)' get_sky_rgb'

        write(6,*)' slice from left to right through midpoint'
        write(6,*)'   i   j   alt   cloud    skyrgb'

        do i = 1,ni
        do j = 1,nj

         if(alt_a(i,j) .eq. r_missing_data)then
          sky_rgb(:,i,j) = 0.          
         else
          if(j .eq. nj/2 .AND. i .eq. (i/5) * 5)then
              idebug = 1
          else
              idebug = 0
          endif

          rintensity = 250. - ((abs(r_cloud_3d(i,j)-0.6)**2.0) * 380.)

          rintensity = min(rintensity,255.)
          rintensity = max(rintensity,0.)

          grayness = r_cloud_3d(i,j) ** 0.4
          grayness = 0.4 + (0.6 * grayness)

          sky_rgb(0,I,J) = nint(rintensity * grayness)                
          sky_rgb(1,I,J) = nint(rintensity * grayness)                   
          sky_rgb(2,I,J) = nint(rintensity)

          rintensity_glow = min(((glow(i,j)-7.) * 100.),255.)
          clr_red = rintensity_glow * rintensity_glow / 255.
          clr_grn = rintensity_glow * rintensity_glow / 255.
          clr_blu = rintensity_glow                                    

!         Use clear sky values if cloud cover is less than 0.5
          frac_cloud = r_cloud_3d(i,j)
!         frac_cloud = 0.0 ; Testing
          frac_clr = 1.0 - frac_cloud                         
          sky_rgb(0,I,J) = clr_red * frac_clr + sky_rgb(0,I,J) * frac_cloud
          sky_rgb(1,I,J) = clr_grn * frac_clr + sky_rgb(1,I,J) * frac_cloud
          sky_rgb(2,I,J) = clr_blu * frac_clr + sky_rgb(2,I,J) * frac_cloud

          if(idebug .eq. 1)then
              write(6,101)i,j,alt_a(i,j),frac_cloud,nint(sky_rgb(:,i,j))
101           format(2i5,f9.2,f9.3,2x,3i6)
          endif

         endif ! missing data tested via altitude

        enddo ! j
        enddo ! i

        return
        end
