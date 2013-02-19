       subroutine btemp_convert_asc(nlines,nelems,
     &          r_missing_data,image,scale,istatus)
c
       implicit none

       integer i,j
       integer istatus
       integer nelems
       integer nlines

       real  image(nelems,nlines)
       real  r_missing_data
       real  scale

       do j=1,nlines
       do i=1,nelems

          if(image(i,j).ne.r_missing_data)then
             image(i,j)=image(i,j)*scale
          endif

       enddo
       enddo

       write(6,*)'btemp_convert_asc range: ',minval(image),maxval(image)

       return
       end
