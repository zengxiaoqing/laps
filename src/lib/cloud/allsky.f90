
       subroutine cyl_to_polar(cyl,polar,minalt,maxalt,minazi,maxazi,alt_scale,azi_scale,polat,pomag,pox,poy,alt_a,azi_a,rotew,rotz &
                              ,iplo,iphi,jplo,jphi,ni_polar,nj_polar)

       include 'trigd.inc'
       use mem_namelist, ONLY: r_missing_data

       real cyl(minalt:maxalt,minazi:maxazi)
       real polar(iplo:iphi,jplo:jphi)
       real alt_a(iplo:iphi,jplo:jphi)
       real azi_a(iplo:iphi,jplo:jphi)
       real delti,deltj,atan3d

       real ri_a(iplo:iphi,jplo:jphi)
       real rj_a(iplo:iphi,jplo:jphi)

       logical l_stereo /.false./

       rstereo(alt) = tand((90.-alt)/2.) ! 0 at zenith, 1 at horizon

       posign = polat / 90.

       altmin = float(minalt) * alt_scale
       azimin = float(minazi) * azi_scale
       azimax = float(maxazi) * azi_scale

       ri_polar_mid = ((ni_polar - 1.0) / 2.0) + 1.0 ! - float(ni_polar-1) * poy
       rj_polar_mid = ((nj_polar - 1.0) / 2.0) + 1.0 ! + float(nj_polar-1) * pox

       write(6,*)'cyl_to_polar ',ni_polar,nj_polar,ri_polar_mid,rj_polar_mid,posign
       write(6,*)'minalt,maxalt,minazi,maxazi',minalt,maxalt,minazi,maxazi
       write(6,*)'alt_scale,azi_scale',alt_scale,azi_scale
       write(6,*)'altmin,azimin',altmin,azimin
       write(6,*)'polat,pomag,rotew,rotz',polat,pomag,rotew,rotz
       write(6,*)'pox,poy',pox,poy
       write(6,*)'iplo,iphi,jplo,jphi',iplo,iphi,jplo,jphi

!      do iaz = minazi,maxazi,20
!          write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz,cyl(maxalt/2,iaz)
!      enddo ! iaz

       write(6,*)' cyl_to_polar: shape of cyl is ',shape(cyl)
       write(6,*)' iplo,iphi,jplo,jphi = ',iplo,iphi,jplo,jphi

       alt_a = r_missing_data
       azi_a = r_missing_data
       polar = r_missing_data

!      Determine cylindrical indices 'ri_a' and 'rj_a' for each polar grid point
       do ip = iplo,iphi
       do jp = jplo,jphi
           rip = ip; rjp = jp
           delti = rip - ri_polar_mid + float(ni_polar-1) * poy
           deltj = rjp - rj_polar_mid - float(nj_polar-1) * pox
           zenith_angle = sqrt(delti**2 + deltj**2) * (90. / ri_polar_mid)
           zenith_angle = zenith_angle / pomag
           alt = (90. - zenith_angle) * posign
           azi_last = azi
           azi = atan3d(-deltj*posign,delti) ! minus sign on deltj flips left/right

!          Determine sign convention of azimuth range
           if(azimin .lt. 0. .and. azi .gt. 180. .and. azi .gt. azimax)then
               azi2 = azi - 360. ! azi2 is between -180. and 0.
           else
               azi2 = azi
           endif

!          Perform coordinate rotation (around E-W horizon axis)
           if(rotew .ne. 0. .or. rotz .ne. 0.)then
             if(azi2 .ne. azi)then
                 write(6,*)' software error, negative azi',azi,azi2
                 stop
             endif

             alt_orig = alt
             azi_orig = azi ! between 0 and 360

             sindec = SIND(alt)
             cosdec = COSD(alt)
             sinphi = SIND(90.-rotew)
             cosphi = COSD(90.-rotew)
             cosha  = COSD(azi)

             alt=ASIND (sinphi*sindec+cosphi*cosdec*cosha)
             cosarg = (cosphi*sindec-sinphi*cosdec*cosha)/cosd(alt)
             cosarg = min(max(cosarg,-1.),+1.)
             azi=180. - ACOSD(cosarg)

             if(modulo(azi_orig,360.) .gt. 180.)then
                azi = 360.0 - azi
             endif

!            Apply azimuth rotation at the end
             azi = modulo(azi + rotz,360.)
             if(abs(alt) .eq. 90.)azi = azi_orig

             if(ip .eq. ni_polar/2 .and. jp .eq. (jp/50)*50)then
               write(6,*)'rot',alt_orig,azi_orig,alt,azi
             endif

             azi2 = azi
           endif
           
           if(alt*posign .ge. 0.)then
!              i_cyl = nint(alt)    
!              j_cyl = nint(azi) 
!              polar(ip,jp) = cyl(i_cyl,j_cyl)
               ri_a(ip,jp) = (alt -altmin) / alt_scale + 1.0 ! real I index in CYL array, offset to start at 1
               rj_a(ip,jp) = (azi2-azimin) / azi_scale + 1.0 ! real J index in CYL array, offset to start at 1

               if(ip .eq. ni_polar/2 .and. jp .eq. (jp/10)*10)then
                   iprint = 1
               else
                   iprint = 0
               endif
!              if(jp .eq. nint(rj_polar_mid))then
!              if((azi-50.) * (azi_last-50.) .lt. 0. .OR. ip .eq. 257)then ! cross 50.
!              if((azi-50.) * (azi_last-50.) .lt. 0. .AND. jp .gt. 265)then ! cross 50.
               if(iprint .eq. 1)then
                 icyl = nint(ri_a(ip,jp)) + minalt
                 jcyl = nint(rj_a(ip,jp)) + minazi
                 if(icyl .ge. minalt .and. icyl .le. maxalt .and. jcyl .ge. minazi .and. jcyl .le. maxazi)then
                   write(6,1)ip,jp,alt,azi,ri_a(ip,jp),rj_a(ip,jp),icyl,jcyl,cyl(icyl,jcyl)
                 else
                   write(6,1)ip,jp,alt,azi,ri_a(ip,jp),rj_a(ip,jp),icyl,jcyl
                 endif
 1               format('ip/jp/alt/azi/ri/rj/icyl/jcyl/cyl',2i6,2f9.3,2f9.2,2i7,f9.1)
               endif

               alt_a(ip,jp) = alt
               azi_a(ip,jp) = azi 
           endif
       enddo ! jp
       enddo ! ip

       imax = maxalt-minalt+1
       jmax = maxazi-minazi+1
       write(6,*)' imax/jmax = ',imax,jmax

       gamma = 2.2
       where(cyl .ne. r_missing_data)cyl = cyl**gamma

!      Bilinear interpolation is most suited when the output grid spacing
!      less than 1.5 times the input grid spacing (determine this ratio)
       frame_width_deg = 180. / pomag
       polar_sp = frame_width_deg / float(ni_polar)
       polar_over_cyl_sp = polar_sp / alt_scale 
       write(6,*)'polar_sp = ',polar_sp
       write(6,*)'alt_scale = ',alt_scale
       write(6,*)'polar_over_cyl_sp = ',polar_over_cyl_sp
       
       nip_crop = iphi-iplo+1
       njp_crop = jphi-jplo+1

       if(polar_over_cyl_sp .lt. 1.5)then
           write(6,*)' Using bilinear interpolation'
           call bilinear_laps_2d(ri_a(iplo:iphi,jplo:jphi) &
                                ,rj_a(iplo:iphi,jplo:jphi) &
                                ,imax,jmax,nip_crop,njp_crop &
                                ,cyl,polar(iplo:iphi,jplo:jphi))

       else ! allow for antialiasing approach
!          Loop through polar indices           
           write(6,*)' Using antialiasing interpolation method'
           do ip = iplo,iphi
           do jp = jplo,jphi

               if(ip .ge. 440 .and. ip .le. 450 .and. jp .eq. 1153)then
                   iprint = 1
               else
                   iprint = 0
               endif

!              Cylindrical indices are given by 'ri_a' and 'rj_a'

!              Azimuth direction
               rjcyl = rj_a(ip,jp)
               j1 = int(rjcyl)
               j2 = j1+1
               fracj = mod(rjcyl,1.0)
               fracj1 = 1. - fracj
               fracj2 = fracj

!              Find overlap of each circumferential stripe in circularized pixel          
               ricyl = ri_a(ip,jp)
               icen = (nint(ricyl)-1) + minalt

               sumij = 0.
               do jcyl = j1,j2 ! Azimuth
                   sumi = 0.
                   do icylrel = -2,+2 ! Altitude
                       icyl = icen + icylrel
                       relpolarcen = (float(icylrel) + mod(ricyl,1.)) / polar_over_cyl_sp
                       relpolarl = relpolarcen - 0.5 / polar_over_cyl_sp
                       relpolarh = relpolarcen + 0.5 / polar_over_cyl_sp
                       call get_overlap(relpolarl,relpolarh,-0.5,+0.5,r_missing_data,icond,overlap,x1o,x2o)
                       if(icyl .gt. minalt .and. icyl .le. maxalt .and. jcyl .gt. minazi .and. jcyl .le. maxazi)then
                           sumi = sumi + cyl(icyl,jcyl) * overlap 
                           if(iprint .eq. 1)then
                               write(6,21)jcyl,icyl,ricyl,relpolarcen,relpolarl,relpolarh,cyl(icyl,jcyl),overlap,sumi
21                             format(i5,i8,f9.3,' relp',3f9.3,' ovr ',3f9.3)
                           endif
                       endif
                   enddo ! icyl
                   if(jcyl .eq. j1)then
                       sumij = sumij + sumi * fracj1
                   else
                       sumij = sumij + sumi * fracj2
                   endif
               enddo ! jcyl

               polar(ip,jp) = sumij

               if(iprint .eq. 1)then
                   write(6,*)' ip/jp/sumij-polar',ip,jp,sumij,sumij**(1./gamma)
               endif

           enddo ! jp
           enddo ! ip           
       endif

       where(polar .ne. r_missing_data)polar = polar**(1./gamma)

!      write(6,*)' polar array',polar(ni_polar/2,1:nj_polar)

       write(6,*)' alt(iplo,jplo) = ',iplo,jplo,alt_a(iplo,jplo)       
       write(6,*)' alt(iphi,jplo) = ',iphi,jplo,alt_a(iphi,jplo)       
       write(6,*)' alt(iplo,jphi) = ',iplo,jphi,alt_a(iplo,jphi)       
       write(6,*)' alt(iphi,jphi) = ',iphi,jphi,alt_a(iphi,jphi)       

!      Fill in corners of polar array (below horizon) with zero values
       where(alt_a(:,:) .eq. r_missing_data)polar(:,:) = 0.
      
       return
       end
