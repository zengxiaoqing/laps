      subroutine filter_2dx(field,ix,iy,iz,smth)
c
c *** Subprogram:  smooth - Smooth a meteorological field.
c     Author:  Stan Benjamin 
c     Date  :  90-06-15
c
c *** Abstract:  Shapiro smoother. 
c 
c *** Program history log: 
c        85-12-09  S. Benjamin - Original version
c        96-06-16  J. Snook    - Modified to do 3d RAMS fields
c                              - hold array is dynamically allocated
c 
c *** Usage:  call smooth(field,ix,iy,iz,smth) 
c
c *** Input argument list: 
c        field    - real array  field(ix,iy,iz)
c                               Meteorological field
c        ix       - integer     x coordinates of field
c        iy       - integer     y coordinates of field
c        iz       - integer     z coordinates of field
c        smth     - real      
c
c *** Output argument list:   
c        field    - real array  field(ix,iy,iz)
c                               Smoothed meteorological field
c 
c *** Remarks:  Reference:  Shapiro, 1970: "Smoothing, filtering, and
c        boundary effects", Rev. Geophys. Sp. Phys., 359-387.
c
c     This filter is of the type 
c        z(i) = (1-s)z(i) + s(z(i+1)+z(i-1))/2
c     for a filter which is supposed to damp 2dx waves completely
c     but leave 4dx and longer with little damping,
c     it should be run with 2 passes using smth (or s) of 0.5
c     and -0.5.
c_______________________________________________________________________________
c   
      
c
      integer ix,iy,iz,i,j,k,i1,i2,it
c
      real field(ix,iy,iz),
     .     hold(ix,2),
     .     smth,smth1,smth2,smth3,smth4,smth5,
     .     sum1,sum2
c_______________________________________________________________________________
c
      smth1=0.25*smth*smth
      smth2=0.50*smth*(1.-smth)
      smth3=(1.-smth)*(1.-smth)
      smth4=(1.-smth)
      smth5=0.5*smth
c
      do k=1,iz
c
         do j=1,2
         do i=1,ix
            hold(i,j)=0.
         enddo
         enddo
c
         i1=2
         i2=1
         do j=2,iy-1
            it=i1
            i1=i2
            i2=it
            do i=2,ix-1
               sum1=field(i-1,j+1,k)+field(i-1,j-1,k)
     .             +field(i+1,j+1,k)+field(i+1,j-1,k)
               sum2=field(i  ,j+1,k)+field(i+1,j  ,k)
     .             +field(i  ,j-1,k)+field(i-1,j  ,k)
               hold(i,i1)=smth1*sum1+smth2*sum2+smth3*field(i,j,k)
            enddo
            if (j .eq. 2) goto 200
            do i=2,ix-1
               field(i,j-1,k)=hold(i,i2)
            enddo
200         continue
         enddo
c
         do i=2,ix-1
            field(i,iy-1,k)=hold(i,i1)
         enddo
c
         do i=2,ix-1
            field(i,1,k)=smth4*field(i,1,k) 
     .                  +smth5*(field(i-1,1,k)+field(i+1,1,k))
            field(i,iy,k)=smth4*field(i,iy,k) 
     .                   +smth5*(field(i-1,iy,k)+field(i+1,iy,k))
         enddo
c
         do j=2,iy-1
            field(1,j,k)=smth4*field(1,j,k) 
     .                  +smth5*(field(1,j-1,k)+field(1,j+1,k))
            field(ix,j,k)=smth4*field(ix,j,k) 
     .                   +smth5*(field(ix,j-1,k)+field(ix,j+1,k))
         enddo
c
      enddo
c
      return
      end

