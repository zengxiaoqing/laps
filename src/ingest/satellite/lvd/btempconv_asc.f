       subroutine btemp_convert_asc(nlines,nelems,
     &          r_missing_data,image,istatus)
c
       implicit none

       integer i,j
       integer istatus
       integer nelems
       integer nlines

       real*4  image(nelems,nlines)
       real*4  r_missing_data

       do j=1,nlines
       do i=1,nelems

          if(image(i,j).ne.r_missing_data)then
             image(i,j)=image(i,j)/10.0
          endif

       enddo
       enddo

       return
       end
