      subroutine count_range(ndimx,ndimy,ndimch,imax,jmax,nch,
     &isndrdata,imaximum,iminimum,istatus)
c
c routine determines the range for count values for each sounding
c data channel
c
      implicit none

      integer*4 i,j,k
      integer*4 ndimch,ndimy,ndimx(jmax)
      integer*4 imax,jmax,nch
      integer*4 istatus
      integer*4 isndrdata(imax,jmax,nch)
      integer*4 imaximum(nch)
      integer*4 iminimum(nch)
      integer*4 i2_missing_data
      integer*4 maxthresh
      integer*4 icntm(nch)
      integer*4 icnteth(nch)

      call get_i2_missing_data(i2_missing_data,istatus)
      if(istatus.ne.1)goto 999

      maxthresh=65535

      do k=1,ndimch

         imaximum(k)=0
         iminimum(k)=maxthresh
         icntm(k)=0
         icnteth(k)=0

         do j=1,ndimy

           do i=1,ndimx(j)

             if(isndrdata(i,j,k).ne.i2_missing_data.and.
     &          isndrdata(i,j,k).lt.maxthresh)then

               if(isndrdata(i,j,k).gt.imaximum(k))then
                  imaximum(k)=isndrdata(i,j,k)
               endif

               if(isndrdata(i,j,k).lt.iminimum(k))then
                  iminimum(k)=isndrdata(i,j,k)
               endif

             elseif(isndrdata(i,j,k).eq.i2_missing_data)then
               icntm(k)=icntm(k)+1
             elseif(isndrdata(i,j,k).gt.maxthresh)then
               icnteth(k)=icnteth(k)+1
             endif

            enddo
         enddo
      enddo
c
      do k=1,ndimch
       write(6,*)'Ch/I2Missing/ExceedThresh: ',k,icntm(k),icnteth(k)
      enddo

      goto 1000
999   write(6,*)'Error getting i2_missing_data'
 
1000  return
      end
