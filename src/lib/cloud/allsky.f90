
       subroutine cyl_to_polar(cyl,polar,minalt,maxalt,minazi,maxazi,alt_scale,azi_scale,polat,pomag,alt_a,azi_a,ni_polar,nj_polar)

       use mem_namelist, ONLY: r_missing_data

       real cyl(minalt:maxalt,minazi:maxazi)
       real polar(ni_polar,nj_polar)
       real alt_a(ni_polar,nj_polar)
       real azi_a(ni_polar,nj_polar)
       real delti,deltj,atan3d

       real ri_a(ni_polar,nj_polar)
       real rj_a(ni_polar,nj_polar)

       posign = polat / 90.

       altmin = float(minalt) * alt_scale
       azimin = float(minazi) * azi_scale

       ri_polar_mid = ((ni_polar - 1.0) / 2.0) + 1.0
       rj_polar_mid = ((nj_polar - 1.0) / 2.0) + 1.0

       write(6,*)'cyl_to_polar ',ni_polar,nj_polar,ri_polar_mid,rj_polar_mid,posign
       write(6,*)'minalt,maxalt,minazi,maxazi',minalt,maxalt,minazi,maxazi
       write(6,*)'alt_scale,azi_scale',alt_scale,azi_scale
       write(6,*)'altmin,azimin',altmin,azimin
       write(6,*)'polat,pomag',polat,pomag

!      do iaz = minazi,maxazi,20
!          write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz,cyl(maxalt/2,iaz)
!      enddo ! iaz

       write(6,*)' cyl_to_polar: shape of cyl is ',shape(cyl)

       alt_a = r_missing_data
       azi_a = r_missing_data
       polar = r_missing_data

       do ip = 1,ni_polar
       do jp = 1,nj_polar
           rip = ip; rjp = jp
           delti = rip-ri_polar_mid
           deltj = rjp-rj_polar_mid
           zenith_angle = sqrt(delti**2 + deltj**2) * (90. / ri_polar_mid)
           zenith_angle = zenith_angle / pomag
           alt = (90. - zenith_angle) * posign
           if(alt*posign .ge. 0.)then
               azi_last = azi
               azi = atan3d(-deltj*posign,delti) ! minus sign on deltj flips left/right
!              i_cyl = nint(alt)    
!              j_cyl = nint(azi) 
!              polar(ip,jp) = cyl(i_cyl,j_cyl)
               ri_a(ip,jp) = (alt-altmin) / alt_scale + 1.0 ! real I index in CYL array, offset to start at 1
               rj_a(ip,jp) = (azi-azimin) / azi_scale + 1.0 ! real J index in CYL array, offset to start at 1
!              if(jp .eq. nint(rj_polar_mid))then
!              if((azi-50.) * (azi_last-50.) .lt. 0. .OR. ip .eq. 257)then ! cross 50.
               if((azi-50.) * (azi_last-50.) .lt. 0. .AND. jp .gt. 265)then ! cross 50.
                 icyl = nint(ri_a(ip,jp)) + minalt
                 jcyl = nint(rj_a(ip,jp)) + minazi
!                if(icyl .ge. minalt .and. icyl .le. maxalt .and. jcyl .ge. minazi .and. jcyl .le. maxazi)then
!                  write(6,1)ip,jp,alt,azi,ri_a(ip,jp),rj_a(ip,jp),icyl,jcyl,cyl(icyl,jcyl)
!                else
!                  write(6,1)ip,jp,alt,azi,ri_a(ip,jp),rj_a(ip,jp),icyl,jcyl
!                endif
 1               format('ip/jp/alt/azi/ri/rj/icyl/jcyl/cyl',2i5,2f9.2,2f9.2,2i5,f9.1)
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

       call bilinear_laps_2d(ri_a,rj_a,imax,jmax,ni_polar,nj_polar,cyl,polar)

       where(polar .ne. r_missing_data)polar = polar**(1./gamma)

!      write(6,*)' polar array',polar(ni_polar/2,1:nj_polar)

       write(6,*)' alt(1,1) = ',alt_a(1,1)       

!      Fill in corners of polar array (below horizon) with zero values
       where(alt_a(:,:) .eq. r_missing_data)polar(:,:) = 0.
      
       return
       end
