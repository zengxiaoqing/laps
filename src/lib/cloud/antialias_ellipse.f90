
       subroutine antialias_ellipse(radius,ricen,rjcen,aspect_ratio,array,ni,nj,iverbose)

!      radius        radius in pixels of ellipse vertical axis        I
!      ricen         location in pixels of ellipse center             I
!      ni,nj         half size of pixel box to evaluate               I
!      array         fractional area of pixels inside ellipse         O

       real array(-ni:ni,-nj:nj)
       real m,k

       a = radius * aspect_ratio 
       b = radius
       area_sum = 0.

       r_missing_data = 1e37

       do i = -ni,+ni
       do j = -nj,+nj

!        Subdivide grid box into horizontal lines
         sum = 0.
         xlo = float(i) - 0.5
         xhi = float(i) + 0.5
         
         do jsub = 0,9
           y = float(j) + float(jsub)/10. - 0.45
           c = y
           m = 0.
           h = ricen
           k = rjcen
           call line_ellipse(m,c,h,k,a,b,r_missing_data,x1,x2,y1,y2)
           call get_overlap(x1,x2,xlo,xhi,r_missing_data,icond,overlap,x1o,x2o)
           if(iverbose .ge. 2)write(6,1)i,j,y,xlo,xhi,x1,x2,overlap,icond,radius,aspect_ratio
1          format(2i4,f9.2,2f9.2,2f9.4,f9.4,i3,2f9.2)
           sum = sum + overlap
         enddo ! isub

         area = sum / float(10)

         if(iverbose .ge. 1)write(6,2)i,j,ricen,rjcen,area
2        format(' i/j/ricen/rjcen/area (sq pix)',2x,2i3,2x,2f7.2,f10.4)

         array(i,j) = area

         area_sum = area_sum + area

       enddo ! j
       enddo ! i

       if(iverbose .ge. 1)then
         pi = 3.14159
         area_theo = (pi * radius**2) * aspect_ratio
         write(6,3)radius,aspect_ratio,area_sum,area_theo
3        format(' sum of illuminated area (sq pix): rad/asp/area_sum/theo = ',4f9.5)
       endif
       
       return
       end

       subroutine get_overlap(x1,x2,xlo,xhi,r_missing_data,icond,overlap,x1o,x2o)

       if(x1 .le. xlo .and. x2 .ge. xhi)then
         icond = 1
         overlap = 1.0
         x1o = xlo
         x2o = xhi
       elseif(x2 .le. xlo .or. x1 .ge. xhi)then
         icond = 2
         overlap = 0.0
         x1o = r_missing_data
         x2o = r_missing_data
       elseif(x1 .le. xlo .and. x2 .ge. xlo)then
         icond = 3
         overlap = x2 - xlo
         x1o = xlo             
         x2o = x2            
       elseif(x1 .le. xhi .and. x2 .ge. xhi)then
         icond = 4
         overlap = xhi - x1
         x1o = x1              
         x2o = xhi           
       elseif(x1 .eq. r_missing_data .or. x2 .eq. r_missing_data)then
         icond = 5
         overlap = 0.0
         x1o = r_missing_data
         x2o = r_missing_data
       elseif(x1 .ge. xlo .and. x2 .le. xhi)then
         icond = 6
         overlap = x2-x1
         x1o = x1
         x2o = x2
       else
         icond = 7
         overlap = r_missing_data
         x1o = r_missing_data
         x2o = r_missing_data
       endif

       return
       end

       subroutine line_ellipse(m,c,h,k,a,b,r_missing_data,x1,x2,y1,y2)

!      http://www.ambrsoft.com/TrigoCalc/Circles2/Ellipse/EllipseLine.htm

!      Find intersection of line and circle
!      Line is defined as y = mx + c

       real m,k

       e = c-k
       d = c + m*h

       disc = a**2*m**2+b**2-d**2-k**2+2*d*k

       if(disc .lt. 0.)then
         x1 = r_missing_data
         x2 = r_missing_data
         y1 = r_missing_data
         y2 = r_missing_data
       else
         denom = a**2*m**2+b**2
         x1 = (b**2*h - m*a**2*e    - a*b*sqrt(disc))   / denom
         x2 = (b**2*h - m*a**2*e    + a*b*sqrt(disc))   / denom
         y1 = (b**2*d + k*a**2*m**2 - a*b*m*sqrt(disc)) / denom
         y2 = (b**2*d + k*a**2*m**2 + a*b*m*sqrt(disc)) / denom
       endif

       return
       end

       subroutine antialias_phase(radius,ricen,rjcen,aspect_ratio,alt_scale,azi_scale,va,rill,array,ni,nj,iverbose)

!      radius        radius in pixels of ellipse vertical axis        I
!      ricen         location in pixels of ellipse center             I
!      ni,nj         half size of pixel box to evaluate               I
!      aspect_ratio  circle projected onto alt/az grid                I
!      va            vertex angle of center of illuminated limb       I
!                    (0 is up, 90 is right) 
!      array         fractional area of pixels inside ellipse         O

       parameter (rpd = 3.14159/180.)
         
!      radius of ellipse with theta measured from a (major) axis
       rell(a,b,theta) = (a * b) / sqrt( (b*cos(theta))**2 + (a*sin(theta))**2 )
!      rell(a,b,theta) = sqrt( (b*cos(theta))**2 + (a*sin(theta))**2 )

       real array(-ni:ni,-nj:nj)
       real m,k

       xpix = ricen
       ypix = rjcen

       xdeg = xpix * azi_scale / aspect_ratio
       ydeg = ypix * alt_scale

       radpix = radius
       raddeg = sqrt(xdeg**2+ydeg**2)

       bb = radius
       area_sum = 0.

!      Ratio of terminator minor axis / major axis
       aspect2 = abs((2.*rill)-1.)

       r_missing_data = 1e37

       do i = -ni,+ni ! azimuth direction
       do j = -nj,+nj ! altitude direction

!        Subdivide grid box into horizontal lines
         sum = 0.
         if(iverbose .ge. 2)write(6,*)'  i   j    xpix     ypix     rdeg    th_up   th_minor th_major   rell1    rell2  aspect_rat rinc'
         do isub = 0,9 
         do jsub = 0,9 
           xpix = -ricen + float(i) + float(isub)/10. - 0.45 ! azimuth direction relative to object
           ypix = -rjcen + float(j) + float(jsub)/10. - 0.45 ! altitude direction relative to object

           xdeg = xpix * azi_scale / aspect_ratio
           ydeg = ypix * alt_scale

           rdeg = sqrt(xdeg**2+ydeg**2)
           theta_up = atan3(xdeg,ydeg)     ! theta relative to up direction
           theta_minor = theta_up - va*rpd ! theta relative to minor axis (illuminated side)
           theta_major=theta_minor+90.*rpd ! theta relative to major axis

!          Determine r_el1, r_el2 and compare
           rell1 = 0.25
           rell2 = rell(rell1,rell1*aspect2,theta_major)
           if(rdeg .lt. rell1)then
             if(rill .ge. 0.5)then ! gibbous
               if(cos(theta_minor).ge.0.)then ! illuminated half
                 rinc = 1.      
               elseif(rdeg .lt. rell2)then
                 rinc = 1.      
               else
                 rinc = 0.
               endif
             else
               if(cos(theta_minor).lt.0.)then ! unilluminated half
                 rinc = 0.      
               elseif(rdeg .lt. rell2)then
                 rinc = 0.      
               else
                 rinc = 1.
               endif
             endif
           else
             rinc = 0.
           endif
           if(iverbose .ge. 2)write(6,1)i,j,xpix,ypix,rdeg,theta_up/rpd,theta_minor/rpd,theta_major/rpd,rell1,rell2,aspect2,rinc
1          format(2i4,10f9.3)
           sum = sum + rinc    
         enddo ! jsub
         enddo ! isub

         area = sum / float(100)

         if(iverbose .ge. 1)write(6,2)i,j,ricen,rjcen,rill,raddeg,area
2        format(' i/j/ricen/rjcen/rill/raddeg/area (sq pix)',2x,2i3,2x,4f7.2,f10.4)

         array(i,j) = area

         area_sum = area_sum + area

       enddo ! j
       enddo ! i

       if(iverbose .ge. 1)then
         pi = 3.14159
         area_theo = (pi * radius**2) * aspect_ratio
         write(6,3)area_sum,area_theo
3        format(' sum of illuminated area (sq pix): area_sum/theo = ',2f9.5)
       endif

       write(6,*)'array(0,0)=',array(0,0)

       return
       end

       
