
        subroutine get_airglow(alt_a,ni,nj,obs_glow_zen & ! I
                              ,patm,htmsl,horz_dep &      ! I
                              ,airmass_2_topo,frac_lp &   ! I
                              ,clear_rad_c_nt)            ! O

!       Calculate sky glow due to airglow. This takes into
!       account the limb when determining airglow. 

        use mem_namelist, ONLY: earth_radius
        include 'trigd.inc'
        include 'rad.inc'

        real alt_a(ni,nj)
        real clear_rad_c_nt(nc,ni,nj)      ! night sky brightness
                                          ! 3 color radiance (Nanolamberts)

        parameter (nlyr = 2)

        real glow_alt(nc),airglow(nc),airglow_zen(nc),airglow_sum(nc)
        real ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_abv,flyr_blw

        do ialt = 1,ni ! Process all azimuths at once for this altitude

          alt = alt_a(ialt,1)

          glow_alt(:) = 0.
          airglow_sum(:) = 0.

          do ilyr = 1,nlyr

            if(ilyr .eq. 1)then    ! Molecular Oxygen + Sodium
              airglow_zen(1) = 75. ! nL (range from ~60-90 with solar cycle)
              airglow_zen(2) = 75. 
              airglow_zen(3) = 35. 
              ht_lyr = 85000.
              thk_lyr = 10000.
            else                   ! Atomic Oxygen
              airglow_zen(1) = 50. ! nL 
              airglow_zen(2) = 0. 
              airglow_zen(3) = 0. 
              ht_lyr = 225000.
              thk_lyr=  75000.
            endif

            bot_lyr = ht_lyr - thk_lyr/2.
            top_lyr = ht_lyr + thk_lyr/2.

            flyr_abv = max( min( (top_lyr-htmsl)/thklyr, 1.) ,0.)
            flyr_blw = max( min( (htmsl-top_lyr)/thklyr, 1.) ,0.)

            if(ialt .eq. 1)then
              write(6,*)' get_airglow: abv/blw',flyr_abv,flyr_blw
              write(6,*)' get_airglow: htmsl',htmsl
              write(6,*)' get_airglow: patm/obs_glow_zen',patm,obs_glow_zen                    
            endif

            call get_airglow_lyr(htmsl,alt,airglow_zen,ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_blw,flyr_abv,airglow)
            airglow_sum(:) = airglow_sum(:) + airglow(:)

          enddo ! ilyr

!         Fill output array for airglow                         
          do jazi = 1,nj                 
            glow_alt(:) = airglow_sum(:)                           ! (nL)
            clear_rad_c_nt(:,ialt,jazi) = glow_alt(:)              ! (nL)
          enddo ! jazi

!         if(ialt .eq. ni)then
          if(mod(alt,5.) .eq. 2.)then
!           write(6,5)alt,h,fracair,obs_glow_zen,airglow(2),glow_alt(2)
!5          format(' get_airglow: alt/h/fair/obsg/airg/glow alt',f9.2,f10.0,f9.3,3f10.0)
            write(6,5)alt,obs_glow_zen,airglow(2),glow_alt(2)
5           format(' get_airglow: alt/obsg/airg/glow alt',f9.2,3f10.0)
!           if(htmsl .ge. top_lyr)then ! above airglow
!             write(6,*)'   horz_dep_airglow = ',horz_dep_airglow
!           endif
          endif

        enddo ! ialt

        return
        end
       

        subroutine get_airglow_lyr(htmsl,alt,airglow_zen,ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_blw,flyr_abv,airglow)

        use mem_namelist, ONLY: earth_radius
        include 'trigd.inc'
        include 'rad.inc'

!       http://aas.aanda.org/articles/aas/pdf/1998/01/ds1449.pdf (Section 6)
        am_thsh(z,h,r) = 1./ sqrt( 1.-( (r/(r+h))**2 * (sind(z))**2) )

!       http://www.wolframalpha.com/input/?i=integrate+1%2F%28sqrt%281-%28r%2F%28r%2Bh%29%29%5E2+sin%5E2z%29%29+dh+
        am_thsh_indefint(z,h,r) = (2.*h**2 + 4.*h*r + r**2 * cosd(2.*z) + r**2) / &
            (2.*(h+r)*sqrt(1.-r**2*sind(z)**2/(h+r)**2))

        am_thsh_defint(z,h1,h2,r) = am_thsh_indefint(z,h2,r) - am_thsh_indefint(z,h1,r) 

        real airglow(nc),airglow_zen(nc)

!         https://farm7.staticflickr.com/6116/6258799449_17eb754b08_o.jpg
!         http://spaceflight.nasa.gov/gallery/images/station/crew-30/hires/iss030e007397.jpg
          if(htmsl .le. bot_lyr)then ! below airglow
            z = 90. - alt
            h = ht_lyr - htmsl
            thickness = thk_lyr

!           Case when shell is entirely higher than observer,
!           looking above or below horizon
!           am1 = am_thsh(z,h,earth_radius)
 
            h1 = h - thickness/2.
            h2 = h + thickness/2.
            if(h1 .ne. 0. .OR. z .ne. 90.)then
              am2 = (1./thickness) * am_thsh_defint(z,h1,h2,earth_radius)
            else ! indefint at h1 is near zero or can blow up
              am2 = (1./thickness) * am_thsh_indefint(z,h2,earth_radius)
            endif

            airglow(:) = airglow_zen(:) * am2

          elseif(htmsl .ge. top_lyr)then ! above airglow
            horz_dep_airglow = horz_depf(htmsl-top_lyr,earth_radius)
            if(alt .lt. -horz_dep_airglow)then
              z = 90. ! 90. + alt
              htmin_ray = htminf(htmsl,alt,earth_radius) ! msl
              h = top_lyr - htmin_ray
              thickness = 10000.
              h = max(h,thickness) ! approximation for partial thickness

              erad_eff = earth_radius + htmin_ray

!             Case when shell is entirely higher than observer,
!             looking above or below horizon
!             am1 = am_thsh(z,h,earth_radius)
 
              h1 = h - thickness/2.
              h2 = h + thickness/2.
              if(h1 .ne. 0. .OR. z .ne. 90.)then
                am2 = (1./thickness) * am_thsh_defint(z,h1,h2,erad_eff)
              else ! indefint at h1 is near zero or can blow up
                am2 = (1./thickness) * am_thsh_indefint(z,h2,erad_eff)
              endif

              if(alt .gt. -horz_dep)then ! pass twice through shell
                airglow(:) = airglow_zen(:) * am2 * 2.
              else                       ! pass once through shell
                airglow(:) = airglow_zen(:) * am2
              endif

            else
              h = 0. ! flag value
              airglow(:) = 0.

            endif

          else
            airglow(:) = 0.
          endif

        return
        end
