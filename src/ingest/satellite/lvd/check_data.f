       subroutine check(data,r_missing_data,istatus,nx_l,ny_l)
c
       real data(nx_l,ny_l)
c
       icnt = 0
       istatus = 0       ! error return
c
       write(6,*)'Subroutine check - r_missing_data is ',r_missing_data

       do j=1,ny_l
       do i=1,nx_l
         if(data(i,j).ne.r_missing_data)goto 10
c        if(data(i,j).gt.0. .and. data(i,j).lt.bad) go to 10
         if(icnt .ge. -10)then
            write(6,*)' check missing data on model grid at'
     1               ,i,j,data(i,j)
         endif
         icnt = icnt - 1
10     enddo !i
       enddo !j
c
       if(icnt .lt. 0) then
         istatus = icnt
         write(6,*)' check routine missing data fraction is '
     1             ,-float(icnt)/float(nx_l*ny_l)          
       else
         istatus = 1
       endif
c
       return
       end

c------------------------------------------------------------------

      subroutine check_field_ave(nx,ny,data,thresh,istatus)

      implicit none

      integer nx, ny
      integer i,j,nn,istatus

      real data(nx,ny)
      real rmxbtemp,rmnbtemp,btempsum,ave 
      real r_missing_data
      real thresh

      print*
      print*,'Computing field ave: '

      call get_r_missing_data(r_missing_data,istatus)

      rmxbtemp=0.0
      rmnbtemp=999.
      nn=0
      btempsum=0.0
      do j=1,ny
      do i=1,nx
         if(data(i,j).ne.r_missing_data)then
            nn=nn+1
            btempsum=btempsum+data(i,j)
            rmxbtemp=max(data(i,j),rmxbtemp)
            rmnbtemp=min(data(i,j),rmnbtemp)
         endif
      enddo
      enddo
      if(nn.gt.0)then
         ave=btempsum/nn
         print*,'------------------------------------'
         print*,'Field max:  ',rmxbtemp,' (K)'
         print*,'Field min:  ',rmnbtemp,' (K)'
         print*,'Field ave:  ',ave, ' (K)'
c            print*,'Field sdev: ',sdev,' (K)'
c            print*,'Field adev: ',adev,' (K)'
         print*,'------------------------------------'
      else
         print*,'No Stats computed: missing data'
      endif
      print*

      if(ave.lt.thresh)then
         print*,'field avg < threshold. no lvd output.'
         istatus = 0
      else
         istatus = 1
      endif

      return
      end
