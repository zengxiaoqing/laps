

       subroutine contour_settings(a,ni,nj,clow,chigh,cint,zoom,scale)       

       real*4 a(ni,nj)

       call get_r_missing_data(r_missing_data,istatus)

       call array_range(a,ni,nj,rmin,rmax,r_missing_data)

       rmax = rmax / scale
       rmin = rmin / scale

       range = (rmax-rmin) / sqrt(zoom)

       if(range .gt. 2000)then
           cint = 400.
       elseif(range .gt. 600)then
           cint = 100.
       elseif(range .gt. 200)then
           cint = 50.
       elseif(range .gt. 60)then ! From 60-200, cint = 10 (6  - 20 contours)
           cint = 10.
       elseif(range .gt. 30)then ! From  30-60, cint = 5  (6  - 12 contours)
           cint = 5.
       elseif(range .gt. 6)then  ! From   6-30, cint = 2  (3  - 15 contours)
           cint = 2.
       else ! range < 6          ! From   0- 6, cint = 1  (0  -  6 contours)
           cint = 1.
       endif

       cint_2 = cint * 2.

       clow  = int(rmin)/int(cint_2) * int(cint_2)
       chigh = int(rmax)/int(cint) * int(cint) + int(cint)
 
       write(6,*)' Subroutine contour_settings....',range,zoom,cint

       return
       end


       subroutine array_range(a,ni,nj,rmin,rmax,r_missing_data)

       real*4 a(ni,nj)

       rmin =  abs(r_missing_data)
       rmax = -abs(r_missing_data)

       do i = 1,ni
       do j = 1,nj
           if(a(i,j) .ne. r_missing_data)then
               rmin = min(rmin,a(i,j))
               rmax = max(rmax,a(i,j))
           endif
       enddo ! j
       enddo ! i

       return
       end
