      subroutine interp_to_sfc(sfc_2d,field_3d,heights_3d,ni,nj,nk,
     &                         badflag,interp_2d)
c
c=================================================================
c
c     Routine to interpolate a 3-d field to a 2-d surface.  The
c     2-d surface can be the actual "surface" (the topography),
c     or any other 2-d surface within the grid.
c
c     Based on code in the LAPS wind analysis, Steve Albers, FSL.
c
c     Original: 04-09-99  Peter A. Stamus, NOAA/FSL
c     Changes:
c
c=================================================================
c
      real sfc_2d(ni,nj), field_3d(ni,nj,nk), heights_3d(ni,nj,nk)
      real interp_2d(ni,nj)
c
      write(6,*)' Interpolating 3-d field to 2-d surface.'

c
c..... Interpolate from the 3-d grid to the 2-d surface at each point.
c
      do j=1,nj
      do i=1,ni
c
         zlow = height_to_zcoord2(sfc_2d(i,j),heights_3d,ni,nj,nk,
     &                                                  i,j,istatus)
         if(istatus .ne. 1)then
            write(6,*) ' Error in height_to_zcoord2 in interp_to_sfc',
     &                 istatus
            write(6,*) i,j,zlow,sfc_2d(i,j),(heights_3d(i,j,k),k=1,nk)
            return
         endif
c
         klow = max(zlow, 1.)
         khigh = klow + 1
         fraclow = float(khigh) - zlow
         frachigh = 1.0 - fraclow
c
         if( field_3d(i,j,klow)  .eq. badflag .or.
     &       field_3d(i,j,khigh) .eq. badflag) then

            write(6,3333)i,j
 3333       format(' Warning: cannot interpolate to sfc at ',2i5)
            interp_2d(i,j) = badflag

         else

            interp_2d(i,j) = field_3d(i,j,klow ) * fraclow  +
     &                       field_3d(i,j,khigh) * frachigh

         endif
c
      enddo !i
      enddo !j
c
c..... That's all.
c
      return
      end
