
        subroutine get_sky_rgb(r_cloud_3d,cloud_od,r_cloud_rad, &
                               cloud_rad_c, &
                               clear_rad_c, &
                               glow,glow_stars,od_atm_a, &
                               airmass_2_cloud,airmass_2_topo, &
                               topo_swi,topo_albedo, & 
                               aod_2_cloud,aod_2_topo, &
                               alt_a,azi_a,elong_a,ni,nj,sol_alt,sol_az, &
                               moon_alt,moon_az,moon_mag, &
                               sky_rgb)                                 ! O

        use mem_namelist, ONLY: r_missing_data
        include 'trigd.inc'

        addlogs(x,y) = log10(10.**x + 10.**y)

        parameter (nc = 3)

        real r_cloud_3d(ni,nj)      ! cloud opacity
        real cloud_od(ni,nj)        ! cloud optical depth
        real r_cloud_rad(ni,nj)     ! sun to cloud transmissivity (direct+fwd scat)
        real cloud_rad_c(nc,ni,nj)  ! sun to cloud transmissivity (direct+fwd scat) * solar color/int
        real clear_rad_c(nc,ni,nj)  ! integrated fraction of air illuminated by the sun along line of sight                 
                                    ! (accounting for Earth shadow + clouds)
        real clear_rad_c_nt(3)      ! HSV night sky brightness
        real glow(ni,nj)            ! skyglow (log b in nanolamberts)
        real glow_stars(ni,nj)      ! starglow (log b in nanolamberts)
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
        real moon_alt,moon_az,moon_mag

        write(6,*)' get_sky_rgb: sol_alt = ',sol_alt
        write(6,*)' moon alt/az/mag = ',moon_alt,moon_az,moon_mag

!       Brighten resulting image at night
        if(sol_alt .lt. -16.)then
            ramp_night = 1.0
        else
            ramp_night = 1.0
        endif

        if(ni .eq. nj)then ! polar
            write(6,*)' slice from SW to NE through midpoint'
        else
            write(6,*)' slices at 45,225 degrees azimuth'
        endif

        if(ni .eq. nj)write(6,11)
11      format('    i   j      alt      azi     elong   pf_scat   opac       od      alb     cloud  airmass   rad    rinten   airtopo  switopo  topoalb   topood  topovis  cld_visb  glow          skyrgb')

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
              cloud_odl  = -99.9 ! flag value for logging
              bkscat_alb = -99.9 ! flag value for logging
          else
!             convert from opacity to albedo
!             bkscat_alb = r_cloud_3d(i,j) ** 2.0 ! approx opacity to albedo
              cloud_opacity = min(r_cloud_3d(i,j),0.999999)
              cloud_odl = -log(1.0 - cloud_opacity)
!             bksc_eff_od = cloud_odl     * 0.12 ! > .10 due to machine epsilon
              bksc_eff_od = cloud_od(i,j) * 0.10 
              cloud_rad_trans = exp(-bksc_eff_od)
              bkscat_alb = 1.0 - cloud_rad_trans 
              ampl = 0.15 * bkscat_alb
              pf_scat = 0.9 + ampl * (-cosd(elong_a(i,j)))
          endif

          if(sol_alt .ge. -16.)then ! Day/twilight from clear_rad_c array
!             Potential intensity of cloud if it is opaque 
!               (240. is nominal intensity of a white cloud far from the sun)
!               (0.25 is dark cloud base value)                                  
              rintensity(:) = 240. * ( (0.25 + 0.75 * cloud_rad_c(:,i,j)) * pf_scat)
!             rintensity = min(rintensity,255.)
              rintensity = max(rintensity,0.)
          else ! nighttime from default surface lighting
              rintensity(:) = 120.
          endif

          cld_red = nint(rintensity(1))                
          cld_grn = nint(rintensity(2))                        
          cld_blu = nint(rintensity(3))                         

          if(sol_alt .ge. 0.)then ! Daylight from skyglow routine
              arg = glow(i,j) + log10(clear_rad_c(3,i,j))
!             rintensity_glow = min(((glow(i,j)-7.) * 100.),255.)
              rintensity_glow = min(((arg      -7.) * 100.),255.)
              clr_red = rintensity_glow *  rintensity_glow / 255.
              clr_grn = rintensity_glow * (rintensity_glow / 255.)**0.80
              clr_blu = rintensity_glow  
          elseif(sol_alt .ge. -16.)then ! Twilight from clear_rad_c array
              hue = clear_rad_c(1,i,j)
              sat = clear_rad_c(2,i,j)
              if(clear_rad_c(3,i,j) .gt. 1.0)then
                  arg2 = log10(clear_rad_c(3,i,j))
                  arg = glow(i,j) + (arg2 + 5.0*arg2**1.5) * 0.15                                       
              else
                  arg = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15             
              endif
              rintensity_floor = 65. + sol_alt
              rintensity_glow = max(min(((arg     -7.) * 100.),255.),rintensity_floor)
              call hsi_to_rgb(hue,sat,rintensity_glow,clr_red,clr_grn,clr_blu)
!             if(idebug .eq. 1)then
!                 write(6,*)'Clr RGB = ',nint(clr_red),nint(clr_grn),nint(clr_blu)
!             endif
!             clr_red = rintensity_glow * clear_rad_c(2,i,j)
!             clr_blu = rintensity_glow * clear_rad_c(3,i,j)
!             clr_grn = 0.5 * (clr_red + clr_blu)                               
          else ! Night from clear_rad_c array (default flat value of glow)
              call get_clr_rad_nt(alt_a(i,j),azi_a(i,j),clear_rad_c_nt)
              hue = clear_rad_c_nt(1)
              sat = clear_rad_c_nt(2)
              glow_nt = log10(clear_rad_c_nt(3)) ! log nL           

              if(moon_alt .gt. 0.)then ! add moon mag condition
!                 Glow from Rayleigh, no clear_rad crepuscular rays yet
!                 argm = glow(i,j) + log10(clear_rad_c(3,i,j)) * 0.15
                  glow_moon = glow(i,j)          ! log nL                 
                  glow_tot = addlogs(glow_nt,glow_moon)
              else
                  glow_tot = glow_nt
                  glow_moon = 0.
              endif

!             Add in stars. Stars have a background glow of 1.0
              glow_tot = addlogs(glow_tot,glow_stars(i,j))

!             if((idebug .eq. 1 .and. moon_alt .gt. 0.) .OR. glow_stars(i,j) .gt. 1.0)then
              if((idebug .eq. 1) .OR. glow_stars(i,j) .gt. 1.0)then
                  write(6,91)glow_nt,glow_moon,glow_stars(i,j),glow_tot
91                format(' glow: nt/moon/stars/tot = ',4f9.3)
                  idebug = 1 ! for subsequent writing at this grid point
              endif

              rintensity_glow = max(min(((glow_tot - 2.5) * 100.),255.),20.)
              call hsi_to_rgb(hue,sat,rintensity_glow,clr_red,clr_grn,clr_blu)
          endif

          sol_alt_red_thr = 1.0 + (od_atm_a * 40.)

          if(sol_alt .le. sol_alt_red_thr .and. sol_alt .gt. -16.0)then
              redness = min((sol_alt_red_thr - sol_alt) / sol_alt_red_thr,1.0)
          else
              redness = 0.
          endif

          elong_red = 12.0
          if(elong_a(i,j) .le. elong_red)then
!         if(.false.)then                        
              red_elong = (elong_red - elong_a(i,j)) / elong_red
!             write(6,*)' alt/elong/redelong: ',alt_a(i,j),elong_a(i,j),red_elong
              clr_red = clr_red * 1.0
              clr_grn = clr_grn * (1. - redness * red_elong)**0.3 
              clr_blu = clr_blu * (1. - redness * red_elong)
          endif

!                     Rayleigh  Ozone   Mag per optical depth            
          od_atm_g = (0.1451  + .016) / 1.086
          od_2_cloud = (od_atm_g + od_atm_a) * airmass_2_cloud(i,j)

!         Empirical correction to account for bright clouds being visible
!         cloud_visibility = exp(-0.28*airmass_2_cloud(i,j))
!         cloud_visibility = exp(-0.14*airmass_2_cloud(i,j))
!         cloud_visibility = exp(-0.43*od_2_cloud)
          cloud_visibility = exp(-0.71*od_2_cloud) 

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
              topo_visibility = exp(-1.00*od_2_topo)                    

              if(airmass_2_cloud(i,j) .gt. 0. .AND. airmass_2_cloud(i,j) .lt. airmass_2_topo(i,j)) then
                  topo_visibility = topo_visibility * (1.0 - r_cloud_3d(i,j))
              endif

              topo_swi_frac = (max(topo_swi(i,j),001.) / 1000.) ** 0.45
              rtopo_red = 150. * topo_swi_frac * (topo_albedo(1,i,j)/.15)**0.45
              rtopo_grn = 150. * topo_swi_frac * (topo_albedo(2,i,j)/.15)**0.45
              rtopo_blu = 150. * topo_swi_frac * (topo_albedo(3,i,j)/.15)**0.45

              sky_rgb(0,I,J) = nint(rtopo_red*topo_visibility + sky_rgb(0,I,J)*(1.0-topo_visibility) )
              sky_rgb(1,I,J) = nint(rtopo_grn*topo_visibility + sky_rgb(1,I,J)*(1.0-topo_visibility) )
              sky_rgb(2,I,J) = nint(rtopo_blu*topo_visibility + sky_rgb(2,I,J)*(1.0-topo_visibility) )
          else
              od_2_topo = 0.
          endif

          sky_rgb(:,i,j) = min(sky_rgb(:,i,j) * ramp_night,255.)

          if(idebug .eq. 1)then
              if(sol_alt .ge. 0.)then
                  write(6,102)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat,r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                      ,topo_swi(i,j),topo_albedo(1,i,j),od_2_topo,topo_visibility,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j))
              elseif(sol_alt .ge. -16.)then  
                  write(6,103)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat,r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                      ,topo_swi(i,j),topo_albedo(1,i,j),od_2_topo,topo_visibility,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c(:,i,j),nint(clr_red),nint(clr_grn),nint(clr_blu)                              
              else ! night
                  write(6,103)i,j,alt_a(i,j),azi_a(i,j),elong_a(i,j) & 
                      ,pf_scat,r_cloud_3d(i,j),cloud_od(i,j),bkscat_alb &
                      ,frac_cloud,airmass_2_cloud(i,j),r_cloud_rad(i,j),rintensity(1),airmass_2_topo(i,j) &
                      ,topo_swi(i,j),topo_albedo(1,i,j),od_2_topo,topo_visibility,cloud_visibility,rintensity_glow,nint(sky_rgb(:,i,j)),clear_rad_c_nt(:)
              endif
102           format(2i5,3f9.2,6f9.3,f7.4,f9.1,f9.3,f9.1,4f9.3,f9.2,2x,3i6)
103           format(2i5,3f9.2,6f9.3,f7.4,f9.1,f9.3,f9.1,4f9.3,f9.2,2x,3i6,' clrrad',3f10.6,3i6)
          endif

         endif ! missing data tested via altitude

        enddo ! i
        enddo ! j

        return
        end

        subroutine get_clr_rad_nt(alt,azi,clear_rad_c_nt)

        real clear_rad_c_nt(3)      ! HSV night sky brightness
                                    ! Nanolamberts

        z = 90. - alt        
        airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))

        airmass_lit = 0.
        sat_ramp = 1.0
        sat_twi_ramp = 0.6

        rint_alt_ramp = sqrt(airmass)

        glow_lp = 500. ! from city lights (nL)
        glow_alt = glow_lp * rint_alt_ramp

!       HSI
        hue = exp(-airmass_lit*0.2) ! 0:R 1:B
        clear_rad_c_nt(1) = hue                             ! Hue
        clear_rad_c_nt(2) = abs(hue-0.5) * 0.8 * sat_ramp   &                        
                                             * sat_twi_ramp ! Sat
        clear_rad_c_nt(3) = glow_alt                        ! Int

        return
        end
       
        subroutine get_starglow(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi,rlat,rlon,alt_scale,azi_scale,glow_stars)

        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real glow_stars(minalt:maxalt,minazi:maxazi) ! log nL

        parameter (nstars = 320)
        real dec_d(nstars),ra_d(nstars),mag_stars(nstars)
        real alt_stars(nstars),azi_stars(nstars),ext_mag(nstars),lst_deg
        real*8 angdif,jed,r8lon,lst,has(nstars),phi,als,azs,ras,x,y,decr

        character*20 starnames(nstars)

        ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)-3.1415926535897932D0

        write(6,*)' subroutine get_starglow...'
        write(6,*)' lat/lon = ',rlat,rlon             

        rpd = 3.14159265 / 180.
        phi = rlat * rpd
        r8lon = rlon          

        call i4time_to_jd(i4time,jed,istatus)
        call sidereal_time(jed,r8lon,lst)

        lst_deg = lst / rpd
        write(6,*)' sidereal time (deg) = ',lst_deg               

!       Obtain stars data
        call read_stars(nstars,ns,dec_d,ra_d,mag_stars,starnames)

        do is = 1,ns        
          RAS = ra_d(is) * rpd
          DECR=dec_d(is) * rpd

          HAS(is)=angdif(lst,RAS)
          call equ_to_altaz_r(DECR,HAS(is),PHI,ALS,AZS)

          alt_stars(is) = als / rpd
          azi_stars(is) = azs / rpd
        enddo ! is

        I4_elapsed = ishow_timer()

!       Obtain planets data
        do iobj = 2,6
          if(.true.)then
            ns = ns + 1
            write(6,*)' Call sun_planet for ',iobj
            call sun_planet(i4time,iobj,rlat,rlon,dec_d(ns),ra_d(ns),alt_stars(ns),azi_stars(ns),elgms_r4,mag_stars(ns))
            starnames(ns) = 'planet'
          endif
        enddo ! iobj

        I4_elapsed = ishow_timer()

        do is = 1,ns        
          patm = 1.0     

          if(alt_stars(is) .gt. 0.0)then
            call calc_extinction(90.          ,patm,airmass,zenext)
            call calc_extinction(alt_stars(is),patm,airmass,totexto)
            ext_mag(is) = totexto - zenext 
          else
            ext_mag(is) = 0.0
          endif

          if(is .le. 100 .OR. (is .gt. 100 .and. is .ge. ns-4))then
            write(6,7)is,starnames(is),dec_d(is),ra_d(is),has(is)/rpd,alt_stars(is),azi_stars(is),mag_stars(is),ext_mag(is)
7           format('dec/ra/ha/al/az/mg/ext',i4,1x,a20,1x,f9.1,4f9.3,2f9.1)
          endif
        enddo ! is

        glow_stars = 1.0 ! cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

            size_glow_sqdg = 1.0  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size
            size_glow_sqdg = 0.3  ! empirical middle ground 

            alt = alt_a(ialt,jazi)
            azi = azi_a(ialt,jazi)

            if(alt .ge. -2.)then

              alt_cos = min(alt,89.)
              alt_dist = alt_scale / 2.0
              azi_dist = alt_dist / cosd(alt_cos)         

              do is = 1,ns        
                if(abs(alt_stars(is)-alt) .le. alt_dist .AND. abs(azi_stars(is)-azi) .le. azi_dist)then
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    rmag_per_sqarcsec = mag_stars(is) + ext_mag(is) + delta_mag                  

!                   Convert to nanolamberts
!                   glow_stars(ialt,jazi) = 5.0 - (mag_stars(is)+ext_mag(is))*0.4

                    glow_nl = v_to_b(rmag_per_sqarcsec)
                    glow_stars(ialt,jazi) = log10(glow_nl)             

                    if(is .le. 50 .AND. abs(azi_stars(is)-azi) .le. 0.5)then
                        write(6,91)is,rmag_per_sqarcsec,delta_mag,glow_nl,glow_stars(ialt,jazi)
91                      format(' rmag_per_sqarcsec/dmag/glow_nl/glow_stars = ',i4,4f10.3)             
                    endif
                endif ! within star kernel
              enddo ! is
            endif ! alt > -2.
        enddo ! jazi
        enddo ! ialt

        return
        end 

        subroutine read_stars(nstars,ns,dec_d,ra_d,mag_stars,starnames)

        real dec_d(nstars),ra_d(nstars),mag_stars(nstars)
        character*20 starnames(nstars)

        character*150 static_dir,filename
        character*120 cline

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/stars.dat'          

        is = 0
        open(51,file=filename,status='old')
1       read(51,2,err=3,end=9)cline
2       format(a)
3       is = is + 1
!       write(6,*)cline

        if(is .le. 0)then
            do ic = 1,80
                write(6,*)' char ',ic,cline(ic:ic)
            enddo ! ic
        endif

        read(cline,4,err=5)ih,im,dec_d(is),mag_stars(is)
4       format(49x,i2,1x,i2,1x,f5.0,28x,f5.0)
5       continue
!       read(cline(66:70),*,err=6)mag_stars(is)
6       continue

!       Convert coordinates
        ra_d(is) = float(ih)*15. + float(im)/4.

        starnames(is) = cline(30:49) 
!       write(6,*)' name is: ',starnames(is)
!       write(6,*)' mag string is: ',cline(91:95)
!       write(6,*)'ih/im',ih,im,dec_d(is),ra_d(is),mag_stars(is)

        goto 1

9       close(51)

        ns = is

        write(6,*)' read_stars completing with # stars of ',ns

        return
        end
