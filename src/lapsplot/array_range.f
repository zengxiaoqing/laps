

       subroutine contour_settings(a,ni,nj,clow,chigh,cint)

       real*4 a(ni,nj)

       call get_r_missing_data(r_missing_data,istatus)

       call array_range(a,ni,nj,rmin,rmax,r_missing_data)

       range = rmax-rmin

       if(range .gt. 2000)then
           cint = 400.
       elseif(range .gt. 600)then
           cint = 100.
       elseif(range .gt. 200)then
           cint = 50.
       elseif(range .gt. 60)then
           cint = 10.
       elseif(range .gt. 20)then
           cint = 5.
       elseif(range .gt. 5)then
           cint = 2.
       else ! range < 5
           cint = 1.
       endif

       cint_2 = cint * 2.

       clow  = int(rmin)/int(cint_2) * int(cint_2)
       chigh = int(rmax)/int(cint) * int(cint) + int(cint)

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
