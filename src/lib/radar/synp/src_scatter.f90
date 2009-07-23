      module src_scatter
      
      use src_hilfe
      use src_global
      contains
        subroutine rayleighsph(aovrb)


    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

!
!-----------------------------------------------------------------------
!     This routine calculates the scattering amplitudes for a spheriod
!     (oblate/prolate/spherical) samll compared with wavelength using
!     the Rayleigh's approximation.
!-----------------------------------------------------------------------
!
!      common/glhydbdy/glbdy(12),glcad(2,6)
!
!      complex*16  pvz,phy
!      common/glraylei/pvz,phy
!
      complex*16 hat
      real*8 vlght,eccn,eccnsq,lz, ffreq, ddeq, ddmx
    real*8, parameter :: ppi=3.14159265358979d0
      data vlght/3.d8/
!
!===> Calculate lz (L3) for oblate/prolate spheroid.
!
      ffreq  = freq*1e9
      refrc = cmplx(refre,refim)
      ddeq   = deq *0.001
      ddmx   = dmx *0.001
 !     aovrb = glbdy(9)
     
    
      eps=refrc**2
      if(ddeq.le.0.0d0) ddeq=ddmx*min(1.d0,aovrb)/aovrb**(2.d0/3.d0)    !


!
      if(aovrb.gt.1.0d0) then

!
!-->  prolate spheroid:
!
      eccnsq=1.0d0-1.0d0/(aovrb*aovrb)

      eccn=sqrt(eccnsq)
   lz=(0.5d0*dlog((1.0d0+eccn)/(1.0d0-eccn))/eccn-1.0d0)*(1.0d0-eccnsq)/eccnsq
!
      else if(aovrb.lt.1.0d0) then

!
!-->  oblate spheroid:
!
      eccnsq=1.0d0/(aovrb*aovrb)-1.0d0!
    
      eccn=sqrt(eccnsq)

      lz=(1.0d0+eccnsq)*(1.0d0-datan(eccn)/eccn)/eccnsq
!
      else
!
!-->  sphere:
!
      lz=1.d0/3.d0
!
      endif
!
!===> Calculate polarizibility and scattering amplitude components.
!
      hat=(ppi*ffreq/vlght)**2*ddeq**3/6.0d0*(eps-1.0d0)    !

      phy=hat/((eps-1.0d0)*(1.0d0-lz)/2.0d0+1.0d0)    !
      pvz=hat/((eps-1.0d0)*lz+1.0d0)    !
   
!
!
      return
      end subroutine


!
!=======================================================================
!
      subroutine miesphere()


      Implicit none

!
!-----------------------------------------------------------------------
!     This routine calculates the scattering amplitudes for a sphere of
!     arbitrary size using the Mie theory algorithm described by van de
!     Hulst, 1981: "Light Scattering by Small Particles", Chapter 9.
!     This routine uses SI units.  by C. Tang, cxt@lance.colostate.edu
!-----------------------------------------------------------------------
!

      Include 'parameter.incf'
      Include 'variablen.incf'
    
!      common/glhydbdy/glbdy(12),glcad(2,6)
!      common/ensemmtx/stks

      
      real  sem(7), nst
      complex*16  as(mrank),bs(mrank),samp(4,mds)
      real*8  x,wvnm, rpsix(mrank),tx(2),uu,pia,pib,pic,tau,&
            chix0,chix1,chix2,psix0,psix1
      complex*16  dpsiz(mrank),tz(2),ztmp,xix0,xix1,sumh,sumv
      
      integer :: i, ids, idr, n,ia1
      integer, parameter :: nds = 2
      real*8, parameter :: degrad = 1.745329251994329d-2
      real, parameter :: phsct(mds) =(0.d0,1.8d2),thsct(mds)=(0d0,1.8d2)
      real :: freq0
      real*8 :: ddeq
      include 'constants.incf'
    
      
!
!===> Series terminated after n terms
!
     
      refrc = cmplx(refre,refim)
      ddeq   = deq *0.001    !equivolume diameter of the particle in mm


      freq0 = freq*1.e9
      wvnm=2*pi*freq0/vlght
      

      x=wvnm*ddeq/2
    
      n=INT(x+4*x**(1.0d0/3.0d0)+20)



!===> Calculate logarithmic derivatives of Riccati-Bessel functions psiz
!     (of complex argument refrc*x), and ratio of psix (of real argument
!     x, rpsix(i)=psix(i-1)/psix(i)).
!
!-->  the tx and tz are used to recur down to the highest order function
!     that are needed.  first set the nozero starting order values for
!     downward recursion, then recur downward to obtain order n.
!

      nst = n+abs(refrc)*x+15
      tx(1) = (nst+1) / x
      tz(1)= (nst+1) / (x * refrc)
 
      do 10 i=nst-1,n,-1
      ia1=i+1
      tx(2)=(i+ia1)/x-1.0d0/tx(1)
      tx(1)=tx(2)
      ztmp = ia1 / (x*refrc)
      tz(2) = ztmp - 1.0d0 / (tz(1) + ztmp)
      tz(1)= tz(2)
 10   continue
      dpsiz(n)=tz(2)
      rpsix(n)=tx(2)
!
!-->  continue downward recursion to order 1.
!
      do 20 i=n-1,1,-1
      ia1=i+1
      rpsix(i)=(i+ia1)/x-1.0d0/rpsix(ia1)
      ztmp=ia1/x/refrc
      dpsiz(i)=ztmp-1.0d0/(dpsiz(ia1)+ztmp)
 20   continue
!
!===> Recur upward to calculate Riccati-Bessel functions psix and chix
!     all of real argument x.  Then compute the a and b arrays.
!
      psix1=sin(x)
      chix0=cos(x)
      chix1=chix0/x+psix1

      do 30 i=1,n
      ia1=i+1
      psix0=psix1
      psix1=psix1/rpsix(i)

      ztmp=(i+ia1)/(i*ia1*dcmplx(0.0d0,wvnm))
      xix0=dcmplx(psix0,chix0)
      xix1=dcmplx(psix1,chix1)
      tz(1)=dpsiz(i)/refrc+i/x
      tz(2)=dpsiz(i)*refrc+i/x
      as(i)=ztmp*(tz(1)*psix1-psix0)/(tz(1)*xix1-xix0)
      bs(i)=ztmp*(tz(2)*psix1-psix0)/(tz(2)*xix1-xix0)

      chix2=(i+ia1)*chix1/x-chix0
      chix0=chix1
      chix1=chix2
 30   continue
!
!===> Starting with order 1 and 2, recur upward to calculate associated
!     Legendre functions divided by sin(theta) (argument cos(theta), at
!     azimuthal mode m=1).  Then compute scattering amplitudes.
!
      do 50 ids=1,nds

      uu=cos(thsct(ids)*degrad)
      pia=1.0d0
      pib=3.0d0*uu
      sumv=bs(1)+as(1)*uu + bs(2)*pib+as(2)*(6.0d0*uu*uu-3.0d0)
      sumh=as(1)+bs(1)*uu + as(2)*pib+bs(2)*(6.0d0*uu*uu-3.0d0)
      do 40 i=3,n
      pic=((2*i-1)*uu*pib-i*pia)/(i-1)
      tau=i*uu*pic-(i+1)*pib
      pia=pib
      pib=pic
      sumv=sumv + bs(i)*pic+as(i)*tau
      sumh=sumh + as(i)*pic+bs(i)*tau
 40   continue
    
      samp(1,ids)= sumv*cos(phsct(ids)*degrad)
!      samp(2,ids)=-sumh*sin(phsct(ids)*degrad)
!      samp(3,ids)= sumv*sin(phsct(ids)*degrad)
!      samp(4,ids)= sumh*cos(phsct(ids)*degrad)

 50   continue
!
!===> Elements in extinction and backscattering Mueller matrices.
!
    
      sem(1)=-2*real(samp(1,1)*cmplx(0.0,-1.0))
    
      sem(4)= abs(samp(1,2))**2
      sem(6)=-sem(4)

      do 60 idr=1,ndr
        
      stks( 1,idr)= sem(1)
      stks( 2,idr)= 0.0
      stks( 3,idr)= 0.0
      stks( 4,idr)= 0.0
      stks( 5,idr)= 0.0

      stks( 6,idr)= sem(4)
      stks( 7,idr)= 0.0
      stks( 8,idr)= 0.0
      stks( 9,idr)= 0.0
      stks(10,idr)= sem(4)
      stks(11,idr)= 0.0
      stks(12,idr)=-sem(4)
      stks(13,idr)= 0.0
      stks(14,idr)=-sem(4)

 60   continue
!
!

      return
      end subroutine
      
      end module src_scatter
