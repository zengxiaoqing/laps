      Subroutine balstagger(u,v,om,phi,t,nx,ny,nz,wr1,wr2,wr3,p,idstag)
c  This routine takes a standard LAPS field with all variables at
c  each grid point and produces the E-stagger appropriate for applying
c  the dynamic balancing in qbalpe.f  idstag > 0 staggers (LAPS -> stagger)
c  idstag < 0 destaggers (stagger -> LAPS)

      real u(nx,ny,nz),v(nx,ny,nz),om(nx,ny,nz),t(nx,ny,nz), 
     & phi(nx,ny,nz),wr1(nx,ny,nz),wr2(nx,ny,nz),wr3(nx,ny,nz)
     & ,p(nz)

c  gas const
     r=287.04

      if(idstag.eq.0) return 'no staggering accomplished'
      if(idstag.gt.0) then

c set vertical stagger first
c first wind level for balcon is second level in LAPS
c omega is shifted one-half in vertical
c t is shifted likewise
c phi is shifted upward one level like winds
c level nz will be the same as laps for all fields

       do k=1,nz-1
        do j=1,ny
         do i=1,nx
          u(i,j,k)=u(i,j,k+1)
          v(i,j,k)=v(i,j,k+1)
          t(i,j,k)=(t(i,j,k)+t(i,j,k+1))*.5
          om(i,j,k)=(om(i,j,k)+om(i,j,k+1))*.5
         enddo
        enddo
       enddo

c horzizontal stagger

       do k=1,nz
        do j=1,ny-1
         do i=1,nx-1
          u(i,j,k)=(u(i+1,j+1,k)+u(i,j+1,k))*.5
          v(i,j,k)=(v(i+1,j+1,k)+v(i+1,j,k))*.5
          t(i,j,k)=(t(i,j,k)+t(i+1,j+1,k)+t(i,j+1,k)+t(i+1,j,k))*.25
c omega is already on the standard horizontal mesh
         enddo
        enddo
       enddo
c re integrate phi hydrostatically from staggered t using phi 1 from LAPS
       do k=1,nz-1
        do j=1,nj
         do i=1,nx
          phi(i,j,k)=phi(i,j,k)+r*t(i,j,k)*alog(p(k)/p(k+1))
         enddo
        enddo
       enddo
c for level nz 
       do j=nj
        do i=1,nx
          phi(i,j,nz)=phi(i,j,nz)+r*t(i,j,nz)*alog(p(nz-1)/p(nz))
        enddo
       enddo
       
       return

      else! de-stagger

c we begin by writing level 1 staggered into the laps level 1; use work arrays
c winds first
       do j=1,ny
        do i=1,nx
         wr1(i,j,1)=u(i,j,1)
         wr2(i,j,1)=v(i,j,1)
         wr3(i,j,1)=om(i,j,1)
        enddo
       enddo
c horizontal destagger

       do k=2,nz
        do j=2,ny
         do i=2,nx
          wr1(i,j,k)=(u(i,j-1,k-1)+u(i-1,j-1,k-1))*.5
          wr2(i,j,k)=(v(i-1,j,k-1)+v(i-1,j-1,k-1))*.5
          wr3(i,j,k)=(om(i,j,k)+om(i,j,k-1))*.5
         enddo
        enddo
       enddo

c now some of the boundary rows and columns

       do k=2,nz
        wr1(1,1,k)=u(1,1,k-1)
        wr2(1,1,k)=v(1,1,k-1)
        wr3(1,1,k)=(om(1,1,k-1)+om(1,1,k))*.5
       do i=2,nx
        wr1(i,1,k)=wr1(i,2,k)
        wr2(i,1,k)=v(i-1,1,k-1)
        wr3(i,1,k)=(om(i,1,k)+om(i,1,k-1))*.5
       enddo
       do j=2,ny
        wr1(1,j,k)=u(1,j-1,k-1)
        wr2(1,j,k)=wr2(2,j,k)
        wr3(1,j,k)=(om(1,j,k)+om(1,j,k-1))*.5
       enddo
      enddo ! on k

c transfer over to u,v,om arrays

      do k=1,nz
       do j=1,ny
        do i=1,nx
         u(i,j,k)=wr1(i,j,k)
         v(i,j,k)=wr2(i,j,k)
         om(i,j,k)=wr3(i,j,k)
        enddo
       enddo
      enddo
 
c now geopotential

      do k=2,nz
       do j=2,ny
        do i=2,nx
         wr3(i,j,k)=(phi(i-1,j-1,k-1)+phi(i,j-1,k-1)+phi(i-1,j,k-1)
     &                +phi(i,j,k-1))*.25
         wr2(i,j,k)=(t(i-1,j-1,k-1)+t(i,j-1,k-1)+t(i-1,j,k-1)+
     &                +t(i,j,k))*.25
        enddo
       enddo
     enddo
     do k=2,nz
      wr3(1,1,k)=phi(1,1,k-1)
      wr2(1,1,k)=t(1,1,k-1)
     enddo
     do k=2,nz
      do j=2,ny
       wr3(1,j,k)=(phi(1,j,k-1)+phi(1,j-1,k-1_*.5
       wr2(1,j,k)=(t(i,1,k-1)+t(i,1,k-1))*.5
      enddo
      do i=2,nx
       wr3(i,1,k)=(phi(i,1,k-1)+phi(i,1,k-1))*.5
       wr2(i,1,k)=(t(i,1,k-1)+t(i,1,k-1))*.5
      enddo
     enddo
     do j=1,ny
      do i=1,nx
       wr3(i,j,1)=wr3(i,j,2)+r*wr2(i,j,2)*alog(p(2)/p(1))
       wr2(i,j,1)=wr2(i,j,2)
      enddo
     enddo
c call phig to compute temps from heights on the native laps grid
c temps will sill be staggered and must be destaggered below
     call phig(wr3,wr2,nx,ny,nz,0,p,-1)
     do k=2,nz
      do j=1,ny
       do i=1,nx
        wr2(i,j,k)=(wr2(i,j,k)+wr2(i,j,k-1))*.5
       enddo
      enddo
     enddo

     do k=1,nz
      do j=1,ny
       do i=1,nx
        phi(i,j,k)=wr3(i,j,k)
        t(i,j,k)=wr2(i,j,k)
       enddo
      enddo
     enddo

     return
     end


c
      subroutine phig(phi,t,nx,ny,nz,itshif,p,ittop)
c
c *** phig computes heights from temps or vice versa for
c        ittop equal to 1 and -1 respectively.
c
      implicit none
c
      integer*4 nx,ny,nz
     .         ,itshif,ittop
     .         ,i,j,k
c
      real*4 phi(nx,ny,nz),t(nx,ny,nz)
     .      ,p(nz)
     .      ,z(50),r,g,rog,gor,ddz,ddzi
c
      parameter (r=287.053,g=9.80665,rog=r/g,gor=g/r)
c_______________________________________________________________________________
c
      do k=1,nz
         z(k)=alog(p(1)/p(k))
      enddo
c
c *** Only for data that has no temps at lvl 1 (1025mb).
c 
      if (itshif .ne. 0) then
         do j=1,ny
         do i=1,nx
            t(i,j,1)=t(i,j,2)*2.-t(i,j,3)
         enddo
         enddo
      endif
c
c *** Scheme assumes 1000mb ht is correct.
c
      if (ittop .eq. 1) then
         do k=2,nz
            ddz=z(k)-z(k-1)
            do j=1,ny
            do i=1,nx
               phi(i,j,k)=(t(i,j,k))*rog*ddz+phi(i,j,k-1)
            enddo
            enddo
         enddo
      elseif (ittop .eq. -1) then
         do k=2,nz
            ddzi=1./(z(k)-z(k-1))
            do j=1,ny
            do i=1,nx
               t(i,j,k)=gor*(phi(i,j,k)-phi(i,j,k-1))*ddzi
            enddo
            enddo
         enddo
      endif
c
      return
      end
c
     
