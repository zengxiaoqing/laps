      subroutine ldarray_i2r_gwc(image_data,n_elem_max,n_lines_max,
     &istart,jstart,iend,jend,r_image,nelem,nlines,istatus)
c
      implicit none

      integer n_elem_max
      integer n_lines_max
      integer nelem
      integer nlines
      integer istart
      integer jstart
      integer iend
      integer jend
      integer istatus
      integer i,j,ii,jj

      integer image_data(n_elem_max,n_lines_max)

      real    r_image(nelem,nlines)

      istatus = 1
      if(iend-istart+1 .gt. nelem)then
         write(6,*)'Error: i array bounds - ldarray_i2r_gwc'
         write(6,*)'Terminating: ',iend-istart+1, nelem
         istatus=-1
         goto 999
      endif

      if(jend-jstart+1 .gt. nlines)then
         write(6,*)'Error: j array bounds - ldarray_i2r_gwc'
         write(6,*)'Terminating: ',jend-jstart+1, nlines
         istatus=-1
         goto 999
      endif

      jj = 0
      do j=jstart,jend
         jj = jj + 1
         ii = 0
         do i=istart,iend
            ii = ii + 1
            r_image(ii,jj) = float(image_data(i,j))*4.0  ! Convert to 10 bit data
         enddo
      enddo

999   return
      end
