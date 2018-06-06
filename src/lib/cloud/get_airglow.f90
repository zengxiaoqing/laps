
        subroutine get_airglow(alt_a,ni,nj,nc,obs_glow_zen,i4time & ! I
                              ,patm,htmsl,horz_dep &                ! I
                              ,airmass_2_topo,frac_lp &             ! I
                              ,clear_rad_c_nt)                      ! O

!       Calculate sky glow due to airglow. This takes into
!       account the limb when determining airglow. 

!       http://adsbit.harvard.edu//full/1997PASP..109.1181K/0001184.000.html

        use mem_namelist, ONLY: earth_radius
        include 'trigd.inc'
        include 'rad_nodata.inc'

        real alt_a(ni,nj) 
        real clear_rad_c_nt(nc,ni,nj)     ! night sky brightness
                                          ! 3 color radiance (Nanolamberts)

        parameter (nlyr = 2)

        real glow_alt(nc),airglow(nc),airglow_zen(nc),airglow_sum(nc)
        double precision jd
!       real ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_abv,flyr_blw

!       Compute Astronomical Julian Date
        call i4time_to_jd(i4time,jd,istatus)

!       Compute Besselian Year
        by = 1900.0 + (jd - 2415020.31352) / 365.242198781

!       Sunspot Cycle Phase (minimum = 0 degrees)
        sunspot_cycle_deg = modulo(((by - 2019.0) / 11.0) * 360.,360.)
        airglow_zen_nl_log = log10(23.334) - log10(49.5/23.334) * cosd(sunspot_cycle_deg)
        airglow_zen_nl = 10. ** airglow_zen_nl_log ! ranges from 11.0-49.5 with solar cycle
!       airglow_zen_nl = 0.1 ! test

        write(6,*)' get_airglow: Besselian Year is      ',by,i4time
        write(6,*)' get_airglow: sunspot cycle (deg) is ',sunspot_cycle_deg
        write(6,*)' get_airglow: airglow_zen_nl (nL) is ',airglow_zen_nl

        do ialt = 1,ni ! Process all azimuths at once for this altitude

          alt = alt_a(ialt,1)

          glow_alt(:) = 0.
          airglow_sum(:) = 0.

          do ilyr = 1,nlyr

            if(ilyr .eq. 1)then     ! Molecular Oxygen + Sodium
              airglow_zen(1) = airglow_zen_nl 
              airglow_zen(2) = airglow_zen_nl 
              airglow_zen(3) = airglow_zen_nl * .467
              ht_lyr = 85000.
              thk_lyr = 10000.
            else                    ! Atomic Oxygen
              airglow_zen(1) = airglow_zen_nl * .667 ! nL 
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

            call get_airglow_lyr(htmsl,alt,airglow_zen,nc,ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_blw,flyr_abv,airglow)
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
            glow_s10 = nl_to_s10(glow_alt(2))
            glow_mag = s10_to_magsecsq(glow_s10)
            write(6,5)alt,obs_glow_zen,airglow(2),glow_alt(2),glow_s10,glow_mag
5           format(' get_airglow: alt/obsg/airg/glow nl-s10-mag',f9.2,4f10.0,f9.2)
!           if(htmsl .ge. top_lyr)then ! above airglow
!             write(6,*)'   horz_dep_airglow = ',horz_dep_airglow
!           endif
          endif

        enddo ! ialt

        return
        end
       

        subroutine get_airglow_lyr(htmsl,alt,airglow_zen,nc,ht_lyr,thk_lyr,bot_lyr,top_lyr,flyr_blw,flyr_abv,airglow)

        use mem_namelist, ONLY: earth_radius
        include 'trigd.inc'
        include 'rad_nodata.inc'

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

        subroutine get_aurora(alt_a,ni,nj,nc,obs_glow_zen & ! I
                             ,patm,htmsl,horz_dep &         ! I
                             ,airmass_2_topo,frac_lp &      ! I
                             ,clear_rad_c_nt)               ! O

        parameter (nsfc = 2)
        parameter (nseg = 2)

        real alt                      ! I elevation angle
        real htmsl                    ! I observer height MSL
        real htsfc                    ! I sfc (or layer) height 
        real earth_radius             ! I earth radius (meters)
        real r_missing_data           ! I
        real alt_norm(nseg,nsfc)      ! O elevation angle rel to sfc normal
        real dist_to_sfc(nseg,nsfc)   ! O distance to sfc (meters)
        real segl(nseg),segh(nseg)

!       Calculate 

!       Aurora intensity a function of magnetic latitude, longitude, Kp, alt

!       Trace line of sight through sphere at 10000m intervals through region
!       between 100km and 400km altitude. Determine bounds of ray between
!       these altitudes (1st/last interection of 400km altitude). There can          
!       be up to two ray segments through this layer.

        do ialt = -10,-10,-20

          do ihtmsl = 0,450,50

            write(6,*)
            write(6,*)' htmsl-km    alt   htsfc-km    dist2sfc1    dist2sfc2       dist2gnd      altn1   altn2     seglen'

            htmsl = float(ihtmsl) * 1000.
            alt = float(ialt)

            iverbose = 1

            call get_aurora_segments(htmsl,alt,earth_radius,iverbose,segl,segh)

          enddo        
        enddo        
          
        return
        end

        subroutine get_aurora_segments(htmsl,alt,earth_radius,iverbose,segl,segh)

        include 'trigd.inc'

        parameter (nsfc = 2)
        parameter (nseg = 2)

        real alt                      ! I elevation angle
        real htmsl                    ! I observer height MSL
        real htsfc                    ! I sfc (or layer) height 
        real earth_radius             ! I earth radius (meters)
        real r_missing_data           ! I
        real alt_norm(nseg,nsfc)      ! O elevation angle rel to sfc normal
        real dist_to_sfc(nseg,nsfc)   ! O distance to sfc (meters)
        real segl(nseg),segh(nseg)

        horz_depf(htmsl,erad) = acosd(erad/(erad+max(htmsl,0.)))

        rlarge = 1e10
        r_missing_data = 1e37
        earth_radius = 6370e3

        htsfcl = 100e3
        htsfch = 400e3
        htgnd = 0e3

            if(htmsl .le. htsfch .and. htmsl .ge. htsfcl)then
              inlayer = 1
            else
              inlayer = 0
            endif

!           Determine horz_dep
            horz_dep = horz_depf(htmsl,earth_radius)
             
            if(iverbose .eq. 1)write(6,*)

            if(alt .lt. (-horz_dep))then
              nseg_pot_ray = 1
            else
              nseg_pot_ray = 2
            endif

            dist_to_gnd = rlarge
            segl(:) = r_missing_data
            segh(:) = r_missing_data

!           Test for condition with ground
            if(htmsl .gt. htgnd .and. alt .le. 0.)then
                call get_ray_info(alt,htmsl,htgnd,earth_radius,alt_norm &
                                 ,r_missing_data,dist_to_gnd)
                if(dist_to_gnd .eq. r_missing_data)dist_to_gnd = rlarge
            elseif(htmsl .le. htgnd .and. alt .le. 0.)then
                dist_to_gnd = 0.
            elseif(htmsl .ge. htgnd .and. alt .ge. 0.)then
                dist_to_gnd = rlarge
            endif

            do isfc = 1,2
              ihtsfc = 100 + (isfc-1) * 300
              htsfc = float(ihtsfc) * 1000.

              dist_to_sfc(:,isfc) = 0.
              alt_norm(:,isfc) = r_missing_data

!             Test for condition with surface
!             rlarge means the surface isn't intersected, zero means ray starts at that surface
              if(htmsl .gt. htsfc .and. alt .le. 0.)then
                  call get_ray_info2(alt,htmsl,htsfc,earth_radius,alt_norm(:,isfc) &
                                    ,r_missing_data,dist_to_sfc(:,isfc))
              elseif(htmsl .le. htsfc)then
                  call get_ray_info2(alt,htmsl,htsfc,earth_radius,alt_norm(:,isfc) &
                                    ,r_missing_data,dist_to_sfc(:,isfc))
              endif

              if(iverbose .eq. 1)write(6,1)ihtmsl,alt,ihtsfc,dist_to_sfc(:,isfc),dist_to_gnd,alt_norm(:,isfc)
1             format(i8,f10.3,i8,4x,2f13.0,2x,f13.0,2x,2f8.3,2x,f13.0)
           enddo        

!          Segment applies before or after the htmin value for downward rays, or segment 1 for upward
           do iseg = 1,nseg_pot_ray
              if(iverbose .eq. 1)write(6,*)' potential segment ',iseg,inlayer,dist_to_sfc(iseg,:)
              if(htmsl .gt. htsfch .and. iseg .eq. 1)then                  ! above upper layer
                 horz_depl = horz_depf(htmsl-htsfcl,earth_radius+htsfcl)
                 horz_deph = horz_depf(htmsl-htsfch,earth_radius+htsfch)
                 if(iverbose .eq. 1)write(6,*)' above upper layer ',-horz_depl,-horz_deph
                 if(alt .lt. (-horz_deph) .and. alt .gt. (-horz_depl))then ! minimum height is within layer
                    if(iverbose .eq. 1)write(6,*)' min height is within layer ',iseg,inlayer,dist_to_sfc(iseg,:)
                    segl(iseg) = dist_to_sfc(1,2)
                    segh(iseg) = dist_to_sfc(2,2)
                 endif
              elseif(inlayer .eq. 0)then                                   ! outside layer
                 if(dist_to_sfc(iseg,1) .le. rlarge .and. dist_to_sfc(iseg,2) .le. rlarge)then
                    segl(iseg) = dist_to_sfc(iseg,1)
                    segh(iseg) = dist_to_sfc(iseg,2)
                 endif
              elseif(iseg .eq. 1)then                                      ! in layer and segment = 1
                 horz_depl = horz_depf(htmsl-htsfcl,earth_radius+htsfcl)
                 if(iverbose .eq. 1)write(6,*)' within layer horz_dep,horz_depl is',-horz_dep,-horz_depl
                 if(alt .ge. 0.)then
                    segl(iseg) = 0.
                    segh(iseg) = dist_to_sfc(iseg,2)
                 else
                    segl(iseg) = 0.
                    if(alt .lt. 0. .and. alt .gt. (-horz_depl))then        ! minimum height is within layer
                      segh(iseg) = dist_to_sfc(2,2)
                      if(iverbose .eq. 1)write(6,*)' minimum height within layer',iseg,segl(iseg),segh(iseg)
                    else
                      segh(iseg) = dist_to_sfc(iseg,1)
                    endif
                 endif
              endif
           enddo

           if(iverbose .eq. 1)write(6,*)'         segments',ihtmsl,alt,segl(1),segh(1),segl(2),segh(2)

        return
        end

        subroutine get_ray_info2(alt,htmsl,htsfc,earth_radius,alt_norm &
                                ,r_missing_data,dist_to_sfc)

        include 'trigd.inc'

        real alt                  ! I elevation angle
        real htmsl                ! I observer height MSL
        real htsfc                ! I sfc (or layer) height 
        real earth_radius         ! I earth radius (meters)
        real r_missing_data       ! I
        real alt_norm(2)          ! O elevation angle rel to sfc normal
        real dist_to_sfc(2)       ! O distance to sfc (meters)

!       Altitude relative to surface normal (emission angle)
        if(alt .ne. 0.)then
          slope = tand(90. + alt)
          htrad = (earth_radius+htmsl) / (earth_radius+htsfc)
          c = htrad * slope
          call line_ellipse(slope,c,0.,0.,1.,1.,r_missing_data,x1,x2,y1,y2)

          if(y1 .ne. r_missing_data)then
            gc1 = atan2d(+y1,-x1) ! great circle dist forward to ground pt
            if(alt .gt. 0.)then
              alt_norm(1) = alt - gc1
            else
              alt_norm(1) = alt + gc1
            endif
            distrad = sqrt((htrad+x1)**2 + y1**2)
            dist_to_sfc(1) = distrad * (earth_radius+htsfc)
            if(gc1 .lt. 0.)then
              dist_to_sfc(1) = r_missing_data
            endif
          else
            dist_to_sfc(1) = r_missing_data
          endif

          if(y2 .ne. r_missing_data)then
            gc2 = atan2d(+y2,-x2) ! great circle dist forward to ground pt
            if(alt .gt. 0.)then
              alt_norm(2) = alt - gc2
            else
              alt_norm(2) = alt + gc2
            endif
            distrad = sqrt((htrad+x2)**2 + y2**2)
            dist_to_sfc(2) = distrad * (earth_radius+htsfc)
            if(gc2 .lt. 0.)then
              dist_to_sfc(2) = r_missing_data
            endif
          else
            dist_to_sfc(2) = r_missing_data
          endif

          write(6,1)dist_to_sfc,htrad,x1,y1,x2,y2,gc1,gc2
 1        format(' d2sfc,htrad,x1,y1,x2,y2,gc1,gc2 ',2f10.0,5f10.4,2f8.2)
        else
          alt_norm(:) = 0.  
        endif

        return
        end
