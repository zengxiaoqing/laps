       subroutine btemp_convert(nelems,nlines,
     &                        rcnt_to_btemp_lut,
     &                        r_missing_data,
     &                        image)
c
c
c
       implicit none

       integer i,j
       integer nelems
       integer nlines

       real  image(nelems,nlines)
       real  rcnt_to_btemp_lut(0:4095)
       real  r_missing_data,rmin,rmax

       call array_range(image,nelems,nlines,rmin,rmax,r_missing_data)
       write(6,*)' non-missing range before btemp conversion ',rmin,rmax

       do j=1,nlines
       do i=1,nelems

          if(image(i,j).ne.r_missing_data)then
             image(i,j)=rcnt_to_btemp_lut(int(image(i,j)))
          endif

       enddo
       enddo

       call array_range(image,nelems,nlines,rmin,rmax,r_missing_data)
       write(6,*)' non-missing range after  btemp conversion ',rmin,rmax

       return
       end
