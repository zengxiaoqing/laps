      subroutine compute_latlon(nnxp,nnyp,n_staggers
     + ,deltax,xtn_ret,ytn_ret,lats,lons,istatus)


      implicit none

      integer nnxp,nnyp
      integer n_staggers
      integer i,j,k

      real    deltax,deltay
      real    mdlat,mdlon
      real    stagger_ew,stagger_ns
      real    erad

c these use internally
      real    xtn(nnxp)
      real    ytn(nnyp)
c these arrays returned for geodat subroutine
      real    xtn_ret(nnxp)
      real    ytn_ret(nnyp)

      real    xmn1,ymn1

      real    lat(nnxp,nnyp)
      real    lon(nnxp,nnyp)

c A-c staggers contained within these arrays.
      real    lats(nnxp,nnyp,n_staggers)
      real    lons(nnxp,nnyp,n_staggers)

      integer istatus

      print*,'calculate lat/lon at stagger grid points.'

      deltay=deltax

      call get_grid_center(mdlat,mdlon,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error returned: get_grid_center'
         return
      endif

      call get_earth_radius(erad,istatus)
      if(istatus .ne. 1)then
         write(6,*)' Error calling get_earth_radius'
         return
      endif

      do k=1,n_staggers

         if(k.eq.1)then

c           Get X/Y for lower left corner
            CALL POLAR_GP(mdlat,mdlon,XMN1,YMN1,DELTAX,DELTAY,
     +  NNXP,NNYP)
            stagger_ns=0
            stagger_ew=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

            do i=1,nnxp
               xtn_ret(i)=xtn(i)
            enddo
            do j=1,nnyp
               ytn_ret(j)=ytn(j)
            enddo

         elseif(k.eq.2)then    !this is the B stagger (E-W stagger)

            stagger_ew=0.5*deltax
            stagger_ns=0
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

         elseif(k.eq.3)then !this is the C stagger (N-S stagger)

            call get_grid_center(mdlat,mdlon,istatus)
            call get_xytn(deltax,deltay,nnxp,nnyp,0,0,xmn1,ymn1
     +,xtn,ytn)
            stagger_ew=0
            stagger_ns=0.5*deltay 
            call get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn) 

         endif
C
C*****************************************************************
C*  Convert it to lat/lon using the library routines.            *

         Do J = 1,nnyp
         Do I = 1,nnxp

            call xy_to_latlon(xtn(i),ytn(j),erad ! ,90.,std_lon
     1                                          ,lat(I,J),lon(I,J))

c             print *,'i,j,xtn,ytn,pla,lplo=',i,j,xtn,ytn,pla,plo
         enddo
         enddo

c this is the A grid
         do j=1,nnyp
         do i=1,nnxp

            lats(i,j,k)=lat(i,j)
            lons(i,j,k)=lon(i,j)

         enddo
         enddo

      enddo

c     open(40,file='latlon_stagger.dat',form='unformatted')
c     write(40)lats,lons
c     close(40)

      istatus = 1

      return
      end

c ----------------------------------------------

      subroutine get_xytn(deltax,deltay,nnxp,nnyp,stagger_ew
     +,stagger_ns,xmn1,ymn1,xtn,ytn)

      implicit none

      integer  nnxp,nnyp
      integer  i,j

      real    deltax,deltay
      real    stagger_ew,stagger_ns
      real    xmn1,ymn1
      real    xmn(nnxp)
      real    ymn(nnyp)
      real    xtn(nnxp)
      real    ytn(nnyp)

      xmn(1)=xmn1+stagger_ew
      ymn(1)=ymn1+stagger_ns

      DO 600 I=2,NNXP
         XMN(I)=XMN(I-1)+DELTAX
 600  CONTINUE
      XMN(NNXP)=2*XMN(NNXP-1)-XMN(NNXP-2)
 
      DO 610 J=2,NNYP
         YMN(J)=YMN(J-1)+DELTAY
 610  CONTINUE
      YMN(NNYP)=2*YMN(NNYP-1)-YMN(NNYP-2)
 
      DO 650 I=2,NNXP
         XTN(I)=.5*(XMN(I)+XMN(I-1))
 650  CONTINUE
      XTN(1)=1.5*XMN(1)-.5*XMN(2)

      DO 660 J=2,NNYP
         YTN(J)=.5*(YMN(J)+YMN(J-1))
 660  CONTINUE
      YTN(1)=1.5*YMN(1)-.5*YMN(2)

      return
      end
