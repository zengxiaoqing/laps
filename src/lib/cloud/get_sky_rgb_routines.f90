
        subroutine get_clr_rad_nt_2d(alt_a,ni,nj,obs_glow_zen & ! I
                                    ,patm,htmsl,horz_dep &      ! I
                                    ,clear_rad_c_nt)            ! O

        use mem_namelist, ONLY: earth_radius
        include 'trigd.inc'
        include 'rad.inc'

        am_thsh(z,h,r) = 1./ sqrt( 1.-( (r/(r+h))**2 * (sind(z))**2) )

        am_thsh_indefint(z,h,r) = (2.*h**2 + 4.*h*r + r**2 * cosd(2.*z) + r**2) / &
            (2.*(h+r)*sqrt(1.-r**2*sind(z)**2/(h+r)**2))

        am_thsh_defint(z,h1,h2,r) = am_thsh_indefint(z,h2,r) - am_thsh_indefint(z,h1,r) 

        real alt_a(ni,nj)
        real azi_a(ni,nj)
        real clear_rad_c_nt(3,ni,nj)      ! night sky brightness
                                          ! 3 color radiance (Nanolamberts)

        real glow_alt(3),airglow(3),airglow_zen(3)

        airglow_zen(1) = 75. ! nL (range from ~60-90 with solar cycle)
        airglow_zen(2) = 75. 
        airglow_zen(3) = 35. 

        ht_lyr = 85000.
        thk_lyr = 10000.
        bot_lyr = ht_lyr - thk_lyr/2.
        top_lyr = ht_lyr + thk_lyr/2.

        flyr_abv = max( min( (top_lyr-htmsl)/thklyr, 1.) ,0.)
        flyr_blw = max( min( (htmsl-top_lyr)/thklyr, 1.) ,0.)

        write(6,*)' get_clr_rad_nt_2d: abv/blw ',flyr_abv,flyr_blw

        do ialt = 1,ni ! Process all azimuths at once for this altitude

          alt = alt_a(ialt,1)
          z = min(90. - alt,91.)        
          airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))

          airmass_lit = 0.
          sat_ramp = 1.0
          sat_twi_ramp = 0.4

!         https://farm7.staticflickr.com/6116/6258799449_17eb754b08_o.jpg
!         http://spaceflight.nasa.gov/gallery/images/station/crew-30/hires/iss030e007397.jpg
          if(htmsl .le. bot_lyr)then ! below airglow
            rint_alt_ramp = sqrt(airmass)

!           airglow(:) = airglow_zen(:) * rint_alt_ramp

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
            horz_dep_airglow = horz_depf(htmsl-90000.,earth_radius)
            if(alt .lt. -horz_dep_airglow)then
              z = 90. ! 90. + alt
              htmin_ray = htminf(htmsl,alt,earth_radius) ! msl
              h = 90000. - htmin_ray
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
              rint_alt_ramp = 1.

            endif

          else
            airglow(:) = 0.
            rint_alt_ramp = 1.
          endif

!         glow_lp = 500. ! from city lights (nL)
          glow_lp = (obs_glow_zen * patm * rint_alt_ramp) ! from city lights (nL)
          glow_alt(:) = glow_lp + airglow(:) ! from city lights + airglow (nL)

!         HSI
!         hue = exp(-airmass_lit*0.2) ! 0:R 1:B
!         clear_rad_c_nt(1) = hue                             ! Hue
!         clear_rad_c_nt(2) = abs(hue-0.5) * 0.8 * sat_ramp   &                        
!                                              * sat_twi_ramp ! Sat

!         Now returned as radiance (nL)
          do ic = 1,3
            clear_rad_c_nt(ic,ialt,:) = glow_alt(ic)          ! Int
          enddo

!         if(ialt .eq. ni)then
          if(ialt .eq. ni .or. alt .eq. nint(alt))then
            write(6,1)alt,h,am2,obs_glow_zen,airglow(2),glow_lp,glow_alt(2)
1           format(' get_clr_rad_nt_2d: alt,h,am2,obsg,airg,glow lp/alt',f9.2,f10.0,f9.3,4f10.0)
!           if(htmsl .ge. 90000.)then ! above airglow
!             write(6,*)'   horz_dep_airglow = ',horz_dep_airglow
!           endif
          endif

        enddo ! ialt

        return
        end
       
        subroutine get_starglow(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi,rlat,rlon,alt_scale,azi_scale,horz_dep,glow_stars)

        include 'rad.inc'

!       http://arxiv.org/pdf/astro-ph/9706111.pdf
!       http://adsabs.harvard.edu/abs/1986PASP...98..364G
!       1 nL = 1.31e6  photons (550nm) sec**-1 cm**-2 sterad**-1
!       1 nL = 1.31e10 photons (550nm) sec**-1  m**-2 sterad**-1
!       1 Watt = 2.76877e18 photons (550nm) / sec
!       1 nL = 4.7313e-9 watts m**-2 sterad**-1 (550nm)
!       wm2sr_to_nl(x) = 2.113e8 * x     ! (550nm)
!       wm2sr_to_nl(x) = 2145.707e9 * x
        wm2srum_to_s10(x) = x / 1.261e-8 ! 550nm
        s10_to_wm2srum(x) = x * 1.261e-8 ! 550nm
        wm2_to_wm2nm(x)     = 0.00136 * x ! solar 550nm approx
        wm2sr_to_wm2srnm(x) = 0.00136 * x ! solar 550nm approx
        wm2sr_to_wm2srum(x) = x * 1.36 ! solar 550nm approx
        wm2srum_to_wm2sr(x) = x / 1.36 ! solar 550nm approx
!       wm2srum_to_nl(x) = wm2sr_to_nl(wm2srum_to_wm2sr(x)) ! sol 550nm
!       s10_to_nl(x) = wm2srum_to_nl(s10_to_wm2srum(x)) ! solar 550nm


        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real glow_stars(3,minalt:maxalt,minazi:maxazi) ! log nL
        real rad_stars(3,minalt:maxalt,minazi:maxazi) ! nL

!       Local arrays
        real dec_a(minalt:maxalt,minazi:maxazi)
        real ha_a(minalt:maxalt,minazi:maxazi)
        real ra_a(minalt:maxalt,minazi:maxazi)
        real gallat_a(minalt:maxalt,minazi:maxazi)
        real gallon_a(minalt:maxalt,minazi:maxazi)
        real helioeclipticlat_a(minalt:maxalt,minazi:maxazi)
        real helioeclipticlon_a(minalt:maxalt,minazi:maxazi)

        character*20 c_conv

        parameter (nstars = 2000)
        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        real alt_stars(nstars),azi_stars(nstars),ext_mag(nstars),lst_deg
        real angdif
        real*8 dangdif,jed,r8lon,lst,has(nstars),phi,als,azs,ras,decr
        real*8 dalt,dazi,dha,ddec,dra,sol_meanlon,sol_meananom

        character*20 starnames(nstars)

        DANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)-3.1415926535897932D0
        ANGDIFD(X,Y)=MOD(X-Y+540.,360.)-180.
!       addmags(a,b)=log10(10.**(-a*0.4) + 10.**(-b*0.4)) * (-2.5)
        addlogs(x,y) = log10(10.**x + 10.**y)

        write(6,*)' subroutine get_starglow...'
        write(6,*)' lat/lon = ',rlat,rlon             

        phi = rlat * rpd
        r8lon = rlon          

        call i4time_to_jd(i4time,jed,istatus)
        call sidereal_time(jed,r8lon,lst)

        lst_deg = lst / rpd
        write(6,*)' sidereal time (deg) = ',lst_deg               

!       Obtain stars data
        call read_bsc(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames)

        do is = 1,ns        
          RAS = ra_d(is) * rpd
          DECR=dec_d(is) * rpd

          HAS(is)=dangdif(lst,RAS)
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
            bmv(ns) = 0.
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
            write(6,7)is,starnames(is),dec_d(is),ra_d(is),has(is)/rpd,alt_stars(is),azi_stars(is),mag_stars(is),bmv(is),ext_mag(is)
7           format('dec/ra/ha/al/az/mg/bmv/ext',i5,1x,a20,1x,f9.1,4f9.3,3f7.1)
          endif
        enddo ! is

!       http://www.ing.iac.es/Astronomy/observing/conditions/skybr/skybr.html
!       http://aas.aanda.org/articles/aas/pdf/1998/01/ds1449.pdf
!       http://download.springer.com/static/pdf/80/art%253A10.1186%252FBF03352136.pdf?auth66=1427760212_d8fc6d84b219c957974dd2b221bc226a&ext=.pdf
        glow_stars = 1.0 ! initialize to cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

!       Coordinate conversions
        nalt = (maxalt-minalt) + 1
        nazi = (maxazi-minazi) + 1
        c_conv = 'altaz2decra'
        call angcoord_2d(c_conv,nalt,nazi,lst_deg,arg2,alt_a,azi_a,rlat,dec_a,ra_a)

        c_conv = 'decra2galactic'
        call angcoord_2d(c_conv,nalt,nazi,arg1,arg2,dec_a,ra_a,argphi,gallat_a,gallon_a)

!       http://en.wikipedia.org/wiki/Position_of_the_Sun
        sol_meanlon  = 280.460 + 0.9856474 * (jed-2451545.0D0)
        sol_meananom = 357.528 + 0.9856003 * (jed-2451545.0D0)
        sollon = sol_meanlon + 1.915 * sind(sol_meananom) + .02 * sind(2. * sol_meananom)
        sollon = mod(sollon,360.)
        write(6,*)' sollon (computed) ',sollon
        c_conv = 'decra2helioecliptic'
        call angcoord_2d(c_conv,nalt,nazi,sollon,arg2,dec_a,ra_a,argphi,helioeclipticlat_a,helioeclipticlon_a)

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

            altg = alt_a(ialt,jazi)
            azig = azi_a(ialt,jazi)

            size_glow_sqdg = 0.3  ! empirical middle ground 
            size_glow_sqdg = 0.5  ! empirical middle ground 
            size_glow_sqdg = 0.0625  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size?

            if(altg .ge. -horz_dep)then

!             Milky Way (mainly from stars from 6-16 magnitude)
              explatterm = (abs(gallat_a(ialt,jazi)/20.))**1.5
              explonterm = (angdifd(gallon_a(ialt,jazi),0.)/140.)**2
              s10gal = 25. + 250. * exp(-explatterm) * exp(-explonterm) 
              if(s10gal .gt. 0.)then
                vgal = s10_to_magsecsq(s10gal)
                glow_gal = log10(v_to_b(vgal))
              else
                glow_gal = 0.
              endif

!             Zodiacal Light S10 units 140 - 90 sin(beta) ecliptic latitude
              elat = helioeclipticlat_a(ialt,jazi)
              elon = helioeclipticlon_a(ialt,jazi)
              if(elon .gt. 180.)elon = 360. - elon
              s10zod_pole = 60. ! high ecliptic latitude
              s10geg_ecliptic =   40. * cosd((elon-180.)/2.)**150.0
              s10zod_ecliptic =  900. * cosd((elon)/2.0001)**7.7 + 140. + s10geg_ecliptic
!             s10zod_ecliptic = 4000. * cosd((elon)/2.)**16.0 + 150. +
!             s10geg_ecliptic
!             Note the gegenschein of s10 = 40 can be added near
!             helioecliptic 180.
              eclp = 1.9 - 0.9*abs(elon/180.)
              eclp = 1.0  ! control ecliptic sharpness
              ecllatterm = abs(sind(elat))**eclp
              s10zod = 10.**(log10(s10zod_pole) * ecllatterm + &
                             log10(s10zod_ecliptic) * (1.-ecllatterm))
!             s10zod =            (s10zod_pole) * ecllatterm + &
!                                 (s10zod_ecliptic) * (1.-ecllatterm) 
!             s10zod = s10zod*10. ! debug for visibility

!             F corona section (inner)
!             http://arxiv.org/pdf/0909.1722.pdf (far corona in Fig 1.)
!             "Interplanetary Dust" edited by B. Grun
              sqarcsec_per_sqdeg = 3600.**2
              size_glow_sqdg = 0.2    ! sun/moon area           
              delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
              distr = sqrt(elat**2 + elon**2)
              distr2 = max(distr,0.40) ! account for pixel size
              srad = distr2/0.25
              spowerk =  5.7 +       srad  * 2.0
              spowerf = 10.0 + log10(srad) * 2.0
              spowerf2 = 8.3 +       srad  * 0.16
              smag_eff = -26.7 + 2.5 * min(spowerk,spowerf)
              rmag_per_sqarcsec = smag_eff + delta_mag
              elong = min(sqrt(elon**2 + elat**2),180.)
              coeff = 3.5 - 2.5*cosd(elong/2)**180.
              if(elong .gt. 0.)then
!               pcoeff = 1.0 - elon**2/(elon**2 + elat**2)
!               pcoeff = abs(elat) / sqrt(elon**2 + elat**2)
!               The last exponent controls ecliptic sharpness
                pcoeff = (abs(elat) / sqrt(elon**2 + elat**2))**1.3 
              else
                pcoeff = 0.
              endif
              fexp = -2.5*(1.-pcoeff) + -2.8*pcoeff
              f_wm2srum = coeff * srad**fexp
              f_s10 = wm2srum_to_s10(f_wm2srum)
!             f_nl  = wm2srum_to_nl(f_wm2srum)

              s10zod2 = max(s10zod,f_s10)

              sbu = s10zod2 * 1.25
              vzod = s10_to_magsecsq(s10zod2)
              glow_zod = log10(v_to_b(vzod))

              do ic = 1,3
                glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi) + 10.**glow_zod + 10.**glow_gal)             
!               glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi))             
              enddo ! ic

              if((azig .eq. 0.0 .or. azig .eq. 180.) .AND. &
                  (altg .eq. 0.0 .or. altg .eq. 10. .or. altg .eq. 20. .or. altg .eq. 30. .or. altg .eq. 40. .or. altg .eq. 60. .or. altg .eq. 90.) &
                                                                                                                            )then
                write(6,71)altg,azig,dec_a(ialt,jazi),ra_a(ialt,jazi) &
                          ,gallat_a(ialt,jazi),gallon_a(ialt,jazi) &
                          ,helioeclipticlat_a(ialt,jazi),helioeclipticlon_a(ialt,jazi)
71              format(' alt-azi/dec-ra/gallat-lon/helioecllat-lon',4(2x,2f8.2))
                write(6,*)'angdif/explonterm/s10gal',angdifd(gallon_a(ialt,jazi),0.),explonterm,s10gal
                write(6,72)elat,elon,srad,spowerk,spowerf,spowerf2,f_wm2srum,s10zod,f_s10,s10zod2,f_wm2srum*1e8,coeff,fexp
72              format(' elat/elon',6f9.2,f10.6,4f10.0,2f9.2)

                write(6,*)'s10geg_ecl/s10zod_ecl/s10zod',s10geg_ecliptic,s10zod_ecliptic,s10zod
                write(6,*)'s10zod/s10gal/glowstars',s10zod2,s10gal,glow_stars(2,ialt,jazi)
              endif

              alt_cos = min(altg,89.)
              alt_dist = alt_scale / 2.0          * 1.3
              azi_dist = alt_dist / cosd(alt_cos) * 1.3   

              do is = 1,ns        
!               Place star in center of grid box
                alts = nint(alt_stars(is)/alt_scale) * alt_scale
                azis = nint(azi_stars(is)/azi_scale) * azi_scale
                if(abs(alts-altg) .le. alt_dist .AND. abs(azis-azig) .le. azi_dist)then
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    rmag_per_sqarcsec = mag_stars(is) + ext_mag(is) + delta_mag                  

!                   Convert to log nanolamberts
                    glow_nl = log10(v_to_b(rmag_per_sqarcsec))
                    do ic = 1,3
                      if(ic .eq. 1)then
                        glow_delt = (bmv(is)-0.6) * 0.7 * 0.4
                      elseif(ic .eq. 2)then
                        glow_delt = 0.
                      elseif(ic .eq. 3)then
                        glow_delt = -(bmv(is)-0.6) * 0.4
                      endif

                      glow_nl = glow_nl + glow_delt

                      glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi) + 10.**glow_nl)             

                      if(is .le. 50 .AND. abs(azis-azig) .le. azi_scale)then
                        if(ic .eq. 2)then
                          write(6,91)is,mag_stars(is),bmv(is),rmag_per_sqarcsec,delta_mag,glow_nl,glow_stars(ic,ialt,jazi)
91                        format( &
                          ' mag/bmv/rmag_per_sqsec/dmag/glow_nl/glow_stars = ',i4,2f8.2,4f10.3)             
                        endif
                      endif

                    enddo ! ic
                endif ! within star kernel
              enddo ! is
            endif ! alt > -20.
        enddo ! jazi
        enddo ! ialt

        rad_stars(:,:,:) = 10.**glow_stars(:,:,:)
 
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

        subroutine read_bsc(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames)

!       http://tdc-www.harvard.edu/catalogs/bsc5.html

        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        character*20 starnames(nstars)

        character*150 static_dir,filename
        character*120 cline
        character*1 c1_dec

        call get_directory('static',static_dir,len_dir)
        filename = static_dir(1:len_dir)//'/bsc5.dat'          

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

        read(cline,4,err=5)ih,im,c1_dec,id,idm,mag_stars(is),bmv(is)
4       format(75x,i2,i2,4x,a1,i2,i2,14x,f5.0,2x,f5.2)
5       continue
!       read(cline(66:70),*,err=6)mag_stars(is)
6       continue

!       Convert coordinates
        if(c1_dec .eq. '+')then
            rsign = 1.0
        else
            rsign = -1.0
        endif
        dec_d(is) = rsign * (float(id) + float(idm)/60.)
        ra_d(is) = float(ih)*15. + float(im)/4.

        starnames(is) = cline(5:14) 

        if(dec_d(is) .eq. 0. .and. ra_d(is) .eq. 0)then ! QC
            is = is - 1
            iqc = 0
        else
            iqc = 1
        endif

        if(mag_stars(is) .lt. 2.0 .and. iqc .eq. 1)then ! magnitude limit
!           write(6,*)' name is: ',
!           write(6,*)' mag string is: ',cline(91:95)
            write(6,13)is,starnames(is),ih,im,dec_d(is),ra_d(is),mag_stars(is),bmv(is)
13          format(i4,1x,a10,' ih/im',2i4,2f7.2,' mag/bmv',2f7.1)
        endif

        if(mag_stars(is) .gt. 5.0)then ! magnitude limit
            is = is - 1
        endif

        goto 1

9       close(51)

        ns = is

        write(6,*)' read_bsc completing with # stars of ',ns

        return
        end
       
        subroutine get_glow_obj(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi &
                               ,alt_scale,azi_scale &
                               ,alt_obj_in,azi_obj_in,mag_obj &
                               ,diam_deg,horz_dep,glow_obj)

        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real glow_obj(minalt:maxalt,minazi:maxazi) ! log nL

        real ext_mag,lst_deg,mag_obj
        real*8 angdif,jed,lst,has,als,azs,ras,x,y,decr

        character*20 starnames

        ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)-3.1415926535897932D0

        write(6,*)' subroutine get_glow_obj...'
        write(6,21)alt_obj_in,azi_obj_in,mag_obj,diam_deg
21      format(' alt/az/mag/diam = ',4f9.3)

        r_missing_data = 1e37

!       Place object in center of grid box
        alt_obj = nint(alt_obj_in/alt_scale) * alt_scale
        azi_obj = nint(azi_obj_in/azi_scale) * azi_scale

        rpd = 3.14159265 / 180.

        I4_elapsed = ishow_timer()

!       Calculate extinction
        patm = 1.0     

        if(alt_obj .gt. 0.0)then
            call calc_extinction(90.    ,patm,airmass,zenext)
            call calc_extinction(alt_obj,patm,airmass,totexto)
            ext_mag = totexto - zenext 
        else
            ext_mag = 0.0
        endif
        ext_mag = 0. ! test of attenuation in 'get_sky_rgb'

        write(6,7)dec_d,ra_d,has/rpd,alt_obj,azi_obj,mag_obj,ext_mag
7       format('dec/ra/ha/al/az/mg/ext',1x,f9.1,4f9.3,2f9.1)

        glow_obj = 1.0 ! cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

!       if(mag_obj .lt. -5.)then ! sun or moon (0.5 degree diameter)
        if(diam_deg .ge. 0.25)then ! sun or moon (0.5 degree diameter)
            alt_dist = diam_deg / 2.0
        else                     ! star or planet (point object)
            alt_dist = alt_scale / 2.0
        endif

        radius_deg = diam_deg / 2.0
    
        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

            size_glow_sqdg = 1.0  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size
            size_glow_sqdg = 0.3  ! empirical middle ground 

            altg = alt_a(ialt,jazi)
            azig = azi_a(ialt,jazi)

            if(alt .ge. -2. .or. alt .ge. -horz_dep)then

                alt_cos = min(altg,89.)
                azi_dist = alt_dist / cosd(alt_cos)         

!               Calculate distance in object/grid radii in cyl projection
                distr = sqrt(((alt_obj-altg)/alt_dist)**2 + ((azi_obj-azig)/azi_dist)**2)

                frac_lit = 1.0
                grid_frac_obj = 1.0

                if(diam_deg .ge. 0.75 .and. distr .le. radius_deg)then  ! solar corona
                    size_glow_sqdg = 0.2    ! sun/moon area           
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    distr = sqrt(((alt_obj-alt))**2 + ((azi_obj-azi)*cosd(alt_cos))**2)
                    if(distr .ge. 0.25)then ! corona
                        distr2 = max(distr,0.40) ! account for pixel size 
                        spowerk =  5.7 +      (distr2/0.25) * 2.0
                        spowerf = 10.0 + log10(distr2/0.25) * 2.0
                        smag_eff = -26.7 + 2.5 * min(spowerk,spowerf)
                        rmag_per_sqarcsec = smag_eff + ext_mag + delta_mag 
                    else                    ! dark side of moon
                        rmag_per_sqarcsec = +2.5 + ext_mag + delta_mag
                    endif
!                   if(abs(alt_obj-alt) .le. 0.15)then                  
                    if(.true.)then                                      
                        write(6,81)ialt,jazi,altg,alt_obj,azig,azi_obj,distr,smag_eff
81                      format(' altg,alt_obj,azig,azi_obj,distr (deg) =',2i4,5f9.3,f7.2)
                    endif
                elseif(diam_deg .ge. 0.25 .and. distr .le. 2.0)then ! regular sun or moon
                    size_glow_sqdg = 0.2    ! sun/moon area           

!                   Calculate fraction of grid box illuminated by object
                    if(.false.)then ! very simple anti-aliasing 
                      if(distr .le. 0.5)then ! object radii
                        frac_lit = 1.0
                      elseif(distr .gt. 1.5)then
                        frac_lit = 0.0
                      else
                        frac_lit = 1.5 - distr
                      endif
                    else           ! more accurate anti-aliasing scheme
                      radius_pix = radius_deg / alt_scale
                      ricen = (azi_obj-azig)/azi_dist
                      rjcen = (alt_obj-altg)/alt_dist
                      aspect_ratio = 1. / cosd(min(abs(alt_obj),89.))
                      call antialias_ellipse(radius_pix,ricen,rjcen,aspect_ratio,frac_lit,0,0,1)
                    endif

                    write(6,89)ialt,jazi,altg,alt_obj,azig,azi_obj,distr,frac_lit
89                  format(' altg/alt_obj/azig/azi_obj/distr (deg)/frac_lit =',2i4,5f9.3,f7.2)

                    if(frac_lit .gt. 0.)then
                        delta_mag = log10((size_glow_sqdg*sqarcsec_per_sqdeg)/frac_lit)*2.5
                        rmag_per_sqarcsec = mag_obj + ext_mag + delta_mag                  
                    else
                        rmag_per_sqarcsec = r_missing_data
                    endif

                elseif(distr .le. radius_deg)then ! star
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    rmag_per_sqarcsec = mag_obj + ext_mag + delta_mag                  
                else
                    rmag_per_sqarcsec = r_missing_data
                endif

                if(rmag_per_sqarcsec .ne. r_missing_data)then
!                 Convert to nanolamberts
                  glow_nl = v_to_b(rmag_per_sqarcsec)
                  glow_obj(ialt,jazi) = log10(glow_nl)             

                  if(.true.)then                          
!                 if(abs(alt_obj-alt) .le. 0.15)then                  
                        if(glow_nl .lt. 1e5)then
                            write(6,91)rmag_per_sqarcsec,delta_mag,glow_nl,glow_obj(ialt,jazi)
91                          format(' rmag_per_sqarcsec/dmag/glow_nl/glow_obj = ',2f10.3,e12.4,f11.2)
                        else
                            write(6,92)rmag_per_sqarcsec,delta_mag,glow_nl,glow_obj(ialt,jazi)
92                          format(' rmag_per_sqarcsec/dmag/glow_nl/glow_obj = ',2f10.3,e12.4,f11.2,' ***GLOW***')
                        endif
                  endif
                endif ! within star kernel
            endif ! alt > -2.
        enddo ! jazi
        enddo ! ialt

        return
        end 

        subroutine get_twi_glow_ave(glow,alt_a,azi_a,ni,nj &
                                   ,sol_alt,sol_az,twi_glow_ave)

        ANGDIF(X,Y)=MOD(X-Y+540.,360.)-180.

        real glow(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)

        cnt = 0.

!       Compare with -22. - sol_alt (as an integrated magnitude)
!       Variable would be 'twi_mag'
        do i = 1,ni
        do j = 1,nj
            diffaz = angdif(azi_a(i,j),sol_az)
            if( abs(diffaz) .le. 90. .AND. &
               (alt_a(i,j)  .le. 20. .and. alt_a(i,j) .ge. 0.) )then            
                cnt = cnt + 1.0
                sum = sum + 10.**glow(i,j)
            endif
        enddo ! j
        enddo ! i

        if(cnt .gt. 0.)then
            twi_glow_ave = log10(sum/cnt)
        else ! tested at -6.2 degrees
            write(6,*)' WARNING in get_twi_glow_ave - using empirical function'
            twi_glow_ave = 8.8 + (sol_alt*0.3)
            write(6,*)' sol_alt/sol_az/twi_glow_ave = ' &
                       ,sol_alt,sol_az,twi_glow_ave
        endif

        return
        end

        subroutine get_sky_rad_ave(rad,alt_a,azi_a,ni,nj &
                                  ,sol_alt,sol_az,sky_rad_ave)

        real rad(ni,nj)
        real alt_a(ni,nj)
        real azi_a(ni,nj)

        cnt = 0.
        sum = 0.

        do i = 1,ni
          cosi = cosd(alt_a(i,1))
          if(alt_a(i,1) .ge. 0.)then ! above horizontal
            do j = 1,nj
                cnt = cnt + (1.0      * cosi)
                sum = sum + (rad(i,j) * cosi)     
            enddo ! j
          endif
        enddo ! i

        if(cnt .gt. 0.)then
            sky_rad_ave = sum/cnt
        else
            write(6,*)' ERROR in get_sky_rad_ave'
            write(6,*)' sol_az = ',sol_az
            sky_rad_ave = 10. 
        endif

        return
        end

        subroutine apply_rel_extinction(rmaglim,alt,od_zen)

        include 'trigd.inc'

        z = 90. - alt
        airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))
        od_slant = od_zen * (airmass - 1.0)
        delta_mag = od_slant * 1.086
        rmaglim = rmaglim - delta_mag

!       write(6,1)airmass,od_zen,od_slant,delta_mag,rmaglim 
1       format(6x,'airmass,od_zen,od_slant,delta_mag,rmaglim ',5f9.2)

        return
        end

        subroutine nl_to_rad(nl,nc,wa,rad)

!       Partly based on equations by Brad Schaffer,  Sky and Telescope, 
!                                                    Februrary 1992

        real nl(nc),wa(nc),rad(nc)

        dist_au = 1.0
        ts = 5700 ! Kelvin

        do ic = 1,nc
            bb = (.0000374*wa(nc)**-5) / (exp(1.43/(wa(ic)*ts))-1.0)
            solar_rad = (1./dist_au**2) * bb * 1e-8 ! erg/cm2/sr/A
            rad(ic) = (nl(ic)/3e9) * solar_rad
        enddo ! ic

        continue

        return
        end

        subroutine get_auroraglow(i4time,alt_a,azi_a,minalt,maxalt,minazi,maxazi,rlat,rlon,alt_scale,azi_scale,glow_aurora)

        real glow_aurora(3,minalt:maxalt,minazi:maxazi) ! log nL
        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)

        return
        end
