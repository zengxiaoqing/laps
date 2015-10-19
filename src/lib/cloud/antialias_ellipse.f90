
       subroutine antialias_ellipse(radius,ricen,rjcen,aspect_ratio,array,ni,nj,iverbose)

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
         do jsub = 0,9
           y = float(j) + float(jsub)/10. - 0.45
           xlo = float(i) - 0.5
           xhi = float(i) + 0.5
           c = y
           m = 0.
           h = ricen
           k = rjcen
           call line_ellipse(m,c,h,k,a,b,r_missing_data,x1,x2,y1,y2)
           if(x1 .le. xlo .and. x2 .ge. xhi)then
             icond = 1
             overlap = 1.0
           elseif(x2 .le. xlo .or. x1 .ge. xhi)then
             icond = 2
             overlap = 0.0
           elseif(x1 .le. xlo .and. x2 .ge. xlo)then
             icond = 3
             overlap = x2 - xlo
           elseif(x1 .le. xhi .and. x2 .ge. xhi)then
             icond = 4
             overlap = xhi - x1
           elseif(x1 .eq. r_missing_data .or. x2 .eq. r_missing_data)then
             icond = 5
             overlap = 0.0
           elseif(x1 .ge. xlo .and. x2 .le. xhi)then
             icond = 6
             overlap = x2-x1
           else
             icond = 7
             overlap = r_missing_data
           endif
           if(iverbose .ge. 2)write(6,1)i,j,y,xlo,xhi,x1,x2,overlap,icond
1          format(2i4,f9.2,2f9.2,2f9.4,f9.4,i3)
           sum = sum + overlap
         enddo ! isub

         area = sum / float(10)

         if(iverbose .ge. 1)write(6,2)i,j,ricen,rjcen,area
2        format(' i/j/ricen/rjcen/area',2x,2i3,2x,2f7.2,f10.4)

         array(i,j) = area

         area_sum = area_sum + area

       enddo ! j
       enddo ! i

       if(iverbose .ge. 1)then
         pi = 3.14159
         area_theo = (pi * radius**2) * aspect_ratio
         write(6,*)'area_sum/theo = ',area_sum,area_theo
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
