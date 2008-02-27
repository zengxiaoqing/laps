      Subroutine check_radar_data(imax,jmax,baddata,data,istatus)
c
c routine checks for remapped data points that
c are out of bounds, 
c
      implicit none
      integer    i,j
      integer    imax,jmax
      integer    istatus
      real     data(imax,jmax)
      real     baddata
      real     r_missing_data
c ============================================
      istatus=0

      call get_r_missing_data(r_missing_data,istatus)
      if(istatus.ne.1)then
         write(6,*)'Error returned: get_r_missing_data'
         goto 990
      endif

      do j=1,jmax
      do i=1,imax
         if(data(i,j).lt.baddata.or.
     &      data(i,j).eq.r_missing_data)goto 100
         data(i,j)=r_missing_data
         istatus=istatus-1
100   enddo
      enddo
      goto 1000

990   write(6,*)'Error in check_radar_data'

1000  return
      end
