       subroutine check(data,r_missing_data,istatus,nx_l,ny_l)
c
       real*4 data(nx_l,ny_l)
c
       icnt = 0
       istatus = 0       ! error return
c
       do j=1,ny_l
       do i=1,nx_l
         if(data(i,j).ne.r_missing_data)goto 10
c        if(data(i,j).gt.0. .and. data(i,j).lt.bad) go to 10
         icnt = icnt - 1
10       enddo !i
       enddo !j
c
       if(icnt .lt. 0) then
         istatus = icnt
       else
         istatus = 1
       endif
c
       return
       end
