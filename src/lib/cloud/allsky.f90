
       subroutine cyl_to_polar(cyl,polar,minalt,maxalt,nj_cyl,alt_scale,azi_scale,alt_a,azi_a,ni_polar,nj_polar)

       use mem_namelist, ONLY: r_missing_data

       real cyl(minalt:maxalt,0:nj_cyl)
       real polar(ni_polar,nj_polar)
       real alt_a(ni_polar,nj_polar)
       real azi_a(ni_polar,nj_polar)
       real delti,deltj,atan3d

       real ri_a(ni_polar,nj_polar)
       real rj_a(ni_polar,nj_polar)

       altmin = float(minalt) * alt_scale

       ri_polar_mid = ((ni_polar - 1.0) / 2.0) + 1.0
       rj_polar_mid = ((nj_polar - 1.0) / 2.0) + 1.0

       write(6,*)'cyl_to_polar ',ni_polar,nj_polar,ri_polar_mid,rj_polar_mid
       write(6,*)'minalt,maxalt,nj_cyl',minalt,maxalt,nj_cyl
       write(6,*)'alt_scale,azi_scale',alt_scale,azi_scale

       do iaz = 0,nj_cyl,20
           write(6,*)'iaz,cyl(maxalt/2,iaz)',iaz,cyl(maxalt/2,iaz)
       enddo ! iaz

       alt_a = r_missing_data
       azi_a = r_missing_data
       polar = r_missing_data

       do ip = 1,ni_polar
       do jp = 1,nj_polar
           rip = ip; rjp = jp
           delti = rip-ri_polar_mid
           deltj = rjp-rj_polar_mid
           zenith_angle = sqrt(delti**2 + deltj**2) * (90. / ri_polar_mid)
           alt = 90. - zenith_angle
           if(alt .ge. 0.)then
               azi = atan3d(+deltj,delti) ! minus sign on deltj flips left/right
!              i_cyl = nint(alt)    
!              j_cyl = nint(azi) 
!              polar(ip,jp) = cyl(i_cyl,j_cyl)
               ri_a(ip,jp) = (alt-altmin) / alt_scale + 1.0 ! real I index in CYL array, offset to start at 1
               rj_a(ip,jp) =  azi         / azi_scale + 1.0 ! real J index in CYL array, offset to start at 1
!              if(jp .eq. nint(rj_polar_mid))then
!                  write(6,*)'ip,jp,ri,rj',ip,jp,ri_a(ip,jp),rj_a(ip,jp)
!              endif
               alt_a(ip,jp) = alt
               azi_a(ip,jp) = azi 
           endif
       enddo ! jp
       enddo ! ip

       imax = maxalt-minalt+1
       jmax = nj_cyl+1

       call bilinear_laps_2d(ri_a,rj_a,imax,jmax,ni_polar,nj_polar,cyl,polar)

       write(6,*)' alt(1,1) = ',alt_a(1,1)       

!      Fill in corners of polar array (below horizon) with zero values
       where(alt_a(:,:) .eq. r_missing_data)polar(:,:) = 0.
      
       return
       end
