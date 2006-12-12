      subroutine image_compare(n_lines,n_elems,r_missing_data,
     &                   image_main,image_prev,
     &                   iqstatus)
c
c This routine compares pixel by pixel image_main and
c image_prev to fill missing data points in image_main
c
      implicit none

      integer n_elems
      integer n_lines
      integer iqstatus
      integer i,j
      real    r_missing_data
      real    image_main(n_elems,n_lines)
      real    image_prev(n_elems,n_lines)
c
      iqstatus=0
      do j=1,n_lines
      do i=1,n_elems

         if(image_main(i,j).eq.r_missing_data)then
            if(image_prev(i,j).ne.r_missing_data)then

               image_main(i,j)=image_prev(i,j)
               iqstatus=iqstatus+1

            endif
         endif
      enddo
      enddo
      return
      end
