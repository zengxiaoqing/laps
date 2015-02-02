
        subroutine get_clr_rad_nt(alt,azi,clear_rad_c_nt)

        include 'trigd.inc'

        real clear_rad_c_nt(3)      ! HSV night sky brightness
                                    ! Nanolamberts

        z = min(90. - alt,91.)        
        airmass = 1. / (cosd(z) + 0.025 * exp(-11 * cosd(z)))

        airmass_lit = 0.
        sat_ramp = 1.0
        sat_twi_ramp = 0.4

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
        real glow_stars(3,minalt:maxalt,minazi:maxazi) ! log nL

        parameter (nstars = 2000)
        real dec_d(nstars),ra_d(nstars),mag_stars(nstars),bmv(nstars)
        real alt_stars(nstars),azi_stars(nstars),ext_mag(nstars),lst_deg
        real*8 angdif,jed,r8lon,lst,has(nstars),phi,als,azs,ras,x,y,decr

        character*20 starnames(nstars)

        ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)-3.1415926535897932D0
!       addmags(a,b)=log10(10.**(-a*0.4) + 10.**(-b*0.4)) * (-2.5)
        addlogs(x,y) = log10(10.**x + 10.**y)

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
        call read_bsc(nstars,ns,dec_d,ra_d,mag_stars,bmv,starnames)

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

        glow_stars = 1.0 ! cosmic background level

        sqarcsec_per_sqdeg = 3600.**2

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

            size_glow_sqdg = 0.3  ! empirical middle ground 
            size_glow_sqdg = 0.5  ! empirical middle ground 
            size_glow_sqdg = 0.0625  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size?

            alt = alt_a(ialt,jazi)
            azi = azi_a(ialt,jazi)

            if(alt .ge. -2.)then

              alt_cos = min(alt,89.)
              alt_dist = alt_scale / 2.0          * 1.3
              azi_dist = alt_dist / cosd(alt_cos) * 1.3   

              do is = 1,ns        
                if(abs(alt_stars(is)-alt) .le. alt_dist .AND. abs(azi_stars(is)-azi) .le. azi_dist)then
                    delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                    rmag_per_sqarcsec = mag_stars(is) + ext_mag(is) + delta_mag                  

!                   Convert to nanolamberts
!                   glow_stars(:,ialt,jazi) = 5.0 - (mag_stars(is)+ext_mag(is))*0.4

                    glow_nl = log10(v_to_b(rmag_per_sqarcsec))
                    do ic = 1,3
!                     glow_stars(:,ialt,jazi) = log10(glow_nl)             
!                     glow_stars(:,ialt,jazi) = addlogs(glow_stars(ic,ialt,jazi),glow_nl)             
                      glow_stars(ic,ialt,jazi) = log10(10.**glow_stars(ic,ialt,jazi) + 10.**glow_nl)             
                    enddo ! ic

                    if(is .le. 50 .AND. abs(azi_stars(is)-azi) .le. 0.5)then
                        write(6,91)is,rmag_per_sqarcsec,delta_mag,glow_nl,glow_stars(2,ialt,jazi)
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
                               ,alt_obj,azi_obj,mag_obj,diam_deg,glow_obj)

        real alt_a(minalt:maxalt,minazi:maxazi)
        real azi_a(minalt:maxalt,minazi:maxazi)
        real glow_obj(minalt:maxalt,minazi:maxazi) ! log nL

        real ext_mag,lst_deg,mag_obj
        real*8 angdif,jed,lst,has,als,azs,ras,x,y,decr

        character*20 starnames

        ANGDIF(X,Y)=DMOD(X-Y+9.4247779607694D0,6.2831853071796D0)-3.1415926535897932D0

        write(6,*)' subroutine get_glow_obj...'
        write(6,*)' alt/az/mag = ',alt_obj,azi_obj,mag_obj

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

        do ialt = minalt,maxalt
        do jazi = minazi,maxazi

            size_glow_sqdg = 1.0  ! alt/az grid
            size_glow_sqdg = 0.1  ! final polar kernel size
            size_glow_sqdg = 0.3  ! empirical middle ground 

            alt = alt_a(ialt,jazi)
            azi = azi_a(ialt,jazi)

            if(alt .ge. -2.)then

              alt_cos = min(alt,89.)
              azi_dist = alt_dist / cosd(alt_cos)         

                distr = sqrt(((alt_obj-alt)/alt_dist)**2 + ((azi_obj-azi)/azi_dist)**2)

                if(abs(alt_obj-alt) .le. alt_dist .AND. abs(azi_obj-azi) .le. azi_dist)then
                    if(diam_deg .ge. 0.75)then  ! solar corona
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
!                       if(abs(alt_obj-alt) .le. 0.15)then                  
                        if(.true.)then                                      
                            write(6,81)ialt,jazi,alt,alt_obj,azi,azi_obj,distr,smag_eff
81                          format(' alt,alt_obj,azi,azi_obj,distr (deg) =',2i4,5f9.3,f7.2)
                        endif
                    elseif(diam_deg .ge. 0.25)then ! regular sun or moon
                        size_glow_sqdg = 0.2    ! sun/moon area           
                        delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                        rmag_per_sqarcsec = mag_obj + ext_mag + delta_mag                  
                    else                        ! star
                        delta_mag = log10(size_glow_sqdg*sqarcsec_per_sqdeg)*2.5
                        rmag_per_sqarcsec = mag_obj + ext_mag + delta_mag                  
                    endif

!                   Convert to nanolamberts
                    glow_nl = v_to_b(rmag_per_sqarcsec)
                    glow_obj(ialt,jazi) = log10(glow_nl)             

                    if(.true.)then                          
!                   if(abs(alt_obj-alt) .le. 0.15)then                  
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
            do j = 1,nj
                cnt = cnt + (1.0      * cosi)
                sum = sum + (rad(i,j) * cosi)     
            enddo ! j
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
