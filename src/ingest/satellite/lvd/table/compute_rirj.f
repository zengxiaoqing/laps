
      Subroutine compute_rirj(nx,ny,n_elems_tot,n_lines_tot,
     &nelems,nlines,rlat,rlon,lat,lon,ri,rj,istatus)
c
c routine generates the ri/rj arrays for mapping laps grid points
c to satellite grid points. 
c
      implicit none

      integer i,j,k,l
      integer nx,ny
      integer lstart, kstart
      integer n_lines_tot
      integer n_elems_tot
      integer nlines
      integer nelems
      integer istatus

      real*4    rlat(n_elems_tot,n_lines_tot)
      real*4    rlon(n_elems_tot,n_lines_tot)
      real*4    ri(nx,ny)
      real*4    rj(nx,ny)
      real*4    iline,jline
      real*4    idiff,jdiff

      real*4    fraclat
      real*4    fraclon
      real*4    rlat_diff_deg
      real*4    rlon_diff_deg

      real*4    lat(nx,ny)
      real*4    lon(nx,ny)

      logical   found_line
      logical   found_elem
      logical   first_time

      integer ishow_timer
      integer init_timer
      integer itstatus

      write(6,*)'This could take awhile, maybe hours.'
      write(6,*)'Some output is coming, and timing stats'
      write(6,*)

      itstatus=init_timer()
      itstatus=ishow_timer()

      lstart=nlines
      kstart=1
      do j = 1,ny
         jline = float(j)/10.
         jdiff = jline - int(jline)
c        first_time = .true.
         l=lstart
      do i = 1,nx
         iline = float(i)/10.
         idiff = iline - int(iline)

         if(idiff.eq.0.00 .and. jdiff.eq.0.00)then
            write(6,29)i,j
29          format(1x,'LAPS(i,j) ',2i6)
         end if

         l=lstart
         found_line=.false.
         found_elem=.false.
         do while(.not.found_line)
            if(l.gt.1)then
            k=kstart
            found_elem=.false.
         do while(.not.found_elem)
            if(k.lt.nelems)then

            if(rlat(i,j).ne.0.0.and.rlon(i,j).ne.0.0)then

            if( (lat(i,j).gt.rlat(k,l)       .and.
     &           lat(i,j).le.rlat(k,l-1) )   .and.
     &           (lon(i,j).gt.rlon(k,l)  .and.
     &            lon(i,j).le.rlon(k+1,l))    )then

               rlat_diff_deg=rlat(k,l-1)-rlat(k,l)
               rlon_diff_deg=rlon(k+1,l)-rlon(k,l)
               if(rlat_diff_deg.eq.0.0000000)then
                  write(6,*)i,j,k,l,'rlat_diff=0.0'
               endif
               if(rlon_diff_deg.eq.0.0000000)then
                  write(6,*)i,j,k,l,'rlon_diff=0.0'
               endif

               fraclat=(lat(i,j)-rlat(k,l))/rlat_diff_deg
               fraclon=(lon(i,j)-rlon(k,l))/rlon_diff_deg

               ri(i,j) = float(k)+fraclon
               rj(i,j) = float(l)-fraclat
               found_line=.true.
               found_elem=.true.
c              if(first_time)then
c                 lstart=l
c              endif

            endif
            k=k+1

            endif

            else
               found_elem=.true.
            endif

         enddo

         else
            found_line=.true.
         endif

         l=l-1

         enddo

      enddo
         if(jdiff.eq.0.00)then
           itstatus=ishow_timer()
         end if
      enddo
c ------------------------------------------------------------------------------
900   return
      end
