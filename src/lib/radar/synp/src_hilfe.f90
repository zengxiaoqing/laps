      module src_hilfe
      
      use src_global
      contains

          function gammln(xx)


      real*4 cof(6),stp,half,one,fpf,x,tmp,ser, xx
      real gammln

      data cof,stp/76.18009173d0, -86.50532033d0, 24.01409822d0,&
         -1.231739516d0,.120858033d-2,-.536382d-5,2.50662827465d0/
      data half,one,fpf/0.5d0,1.0d0,5.5d0/

      x=xx-one
      tmp=x+fpf
      tmp=(x+half)*log(tmp)-tmp
      ser=one
      do 11 j=1,6
         x=x+one
         ser=ser+cof(j)/x
 11      continue
      gammln=exp(tmp+log(stp*ser))
      return
      end function gammln

      subroutine scatamp(thinc,sem,sctmodel, nrank)

    implicit none
      integer :: nrank
      character sctmodel*3
      real sem(7),thinc
!       write(*,*) sctmodel
      if(sctmodel.eq.'tmt') call tmtaddprc(thinc,sem, nrank)
      if(sctmodel.eq.'ray') call rayaddprc(thinc,sem)


      return
      end subroutine

            function cadfunc(cang,cdtyp,calow,caupp,caavr,cadev)
!***********************************************************************
!-----------------------------------------------------------------------
!     This routine defines the functional forms of commonly used canting
!     angle distributions.
!-----------------------------------------------------------------------
include 'variablen.incf'
      integer idtyp
      real ::  ca, cadfunc, tmp1, cang

!===> Check canting angle distribution type.

      idtyp=cdtyp+0.01

      goto (10,20,30) idtyp
!
!===> Uniform orientation probability:
!
 10   cadfunc = 1.0
      goto 100




!      
!===> Gaussian distribution:
!
 20   cadfunc = exp( -0.5*((cang-caavr)/cadev)**2 )
      goto 100
!
!===> Simple harmonic oscillator orientation probability:
!
 30   tmp1=2*(cang-caavr)/(caupp-calow) !double-sided
      if(abs(caavr-calow).lt.1.e-8) tmp1=tmp1/2 !single-sided
      if( abs(tmp1).ge.1 ) tmp1=1-1.0e-6
      cadfunc = 1.0/sqrt(1.0-tmp1*tmp1)
      goto 100


 100  return
      end function

            function erf(x)


!     returns the error function
      real x, erf
      if(x .lt. 0)then
         erf=-gammp(.5,x**2)
      else
         erf=gammp(.5,x**2)
      endif
      return
      end function

       function gammp(a,x)

!    returns the incomplete gamma function
      real x,a,gammp,gammcf, gln
      if((x .lt. 0) .or. (a .le. 0))then
         print*,'problem in the gammp function'
         stop
      endif
      if(x .lt. a+1.)then
         call gser(gammp,a,x,gln)
      else
         call gcf(gammcf,a,x,gln)
         gammp=1.-gammcf
      endif
      return
      end function


       subroutine gser(gamser,a,x,gln)
!     returns the incomplete gamma function P(a,x) evaluated by its series
!     representationas GAMSER
      integer n
      real x,ap,a,sum,del,gamser,gln
      integer, parameter ::itmax=100
      real, parameter :: eps=3.e-7
      gln=log(gammln(a))
      if(x .le. 0)then
         if(x .lt. 0)print*, 'x is less than zero in gser routine'
         gamser=0
         return
      endif
      ap=a
      sum=1./a
      del=sum
      do 11 n=1,itmax
         ap=ap+1.
         del=del*x/ap
         sum=sum+del
         if(abs(del) .lt. abs(sum)*eps)goto 1
 11      continue
      print*, 'A too large, itmax too small'
 1    gamser=sum*exp(-x+a*log(x)-gln)
      return
      end subroutine

      subroutine gcf(gammcf,a,x,gln)
!     returns the incomplete gamma function Q(a,x) evaluated by its 
!     continued fraction representation as GAMMCF
      real gold,a0,a1,b0,b1,fac,an,ana,anf,gammcf,a,x,gln,g
      integer, parameter :: itmax=100
      real, parameter :: eps=3e-7
      integer :: n
      gln=log(gammln(a))
      gold=0.
      a0=1.
      a1=x
      b0=0.
      b1=1.
      fac=1.
      do 11 n=1,itmax
         an=REAL(n)
         ana=an-a
         a0=(a1+a0*ana)*fac
         b0=(b1+b0*ana)*fac
         anf=an*fac
         a1=x*a0+anf*a1
         b1=x*b0+anf*b1
         if(a1 .ne. 0)then
            fac=1./a1
            g=b1*fac
            if(abs((g-gold)/g) .lt. eps)goto 1
            gold=g
         endif
 11      continue
         print*, 'A too large, ITMAX too small'
 1       gammcf=exp(-x+a*alog(x)-gln)*g
         return
         end subroutine


            function fdb(x)
!
!===> Function fbd changes varible x into dB scale.
!
        
        real :: fdb, x
        fdb=10.0*alog10( amax1(1.0e-37,x) )
        
       end function

       subroutine watereps (temp,epswtr,epsice)
!***************************************************************************
implicit none

INCLUDE 'parameter.incf'
include 'fields.incf'

!
!-----------------------------------------------------------------------
!     This routine calculates the complex relative dielectric constant
!     of liquid and ice phase water using an empirical model developed
!     by P. S. Ray, 1972, Applied Optics, 11(8):1836-1844
!     Applicable range: wvln: 1 to 600 mm;  freq: 0.5 to 300 GHz
!     Temperarure (Celsius degree):  ice: -20 to 0, water: -20 to 50
!     exp(+jwt) convention
!-----------------------------------------------------------------------
!


!***********************************************************************

COMPLEX, INTENT(OUT) :: epswtr, epsice

!Locale Variablen
real*4 :: temp, tt, tempx
real :: x, wvln
REAL*8 :: y, sigma2, freqx, xx1, yy1 
REAL   :: einf, einf1, alpha, alpha1, estamf, estamf1
real :: ewu, ew0, tptw, xx, yy, tptwfq
REAL*8, Parameter  :: sigma1=12.5664e8
real           :: ri_real(2), ri_im(2), f, m_real, b_real, m_im, b_im
real           :: T_grid(4), lam_grid(21), lam_low, lam_high, lam
real           :: ri_real_final, ri_im_final, T_high, T_low, testre, testtt, wo1, testim
complex        :: ri(4,21)
integer :: kk
character*255 :: static_dir
integer :: len_dir
include 'constants.incf'

!*************************************************************************
!write(*,*) freq, 'freq'

freqx = freq *1.e9

wvln =1000.* vlght/freqx


 lam = 1e6*2.997925*1e8/freqx ! microns


!***********************************************************************
!Berechnung der dielectrischen Konstanten für Wasser      
!***********************************************************************

einf1=5.27137+(0.0216474-0.131198e-2*temp)*temp      !(7a)

alpha1=-16.8129/(temp+273.0)+0.609265e-1    !(7b)
  
tt=temp    -25.0
     
estamf1=78.54*(1.0-4.579e-3*tt+1.19e-5*tt*tt-2.8e-8*tt**3)-einf1


xx1=((0.33836e-2*exp(2513.98/(temp+273.0)))/wvln)**(1-alpha1)


yy1=1.0+2.0*xx1*sin(alpha1*pi/2.0)+xx1*xx1
  
epswtr=cmplx(estamf1*(1.0+xx1*sin(alpha1*pi/2.0))/yy1+einf1, -estamf1*xx1*cos(alpha1*pi/2.0)/yy1+sigma1*wvln/18.8496e10)

!**************************************************************************
!Berechnung der dielectrischen Konstanten für Eis
!**************************************************************************

tempx = temp
if (temp .gt. -1.) tempx = -1.5
call get_directory('static',static_dir,len_dir)
open(unit=88, file=static_dir(1:len_dir)//'/warren84.tab', status='old')
 read(88,*) T_grid
 T_grid = T_grid 
!print*, T_grid
 read(88,*) lam_grid
! print*, lam_grid
 read(88,*) ri
!print*, ri
!print*, temp, 'temp'
!stop

 if (tempx.lt.T_grid(4).or.tempx.gt.T_grid(1)) then
!  print*,'temperature out of range'
  tempx = tempx+10.
 ! stop
 endif
! write(*,*) lam, lam_grid(21), lam_grid(1)
 if (lam.gt.lam_grid(21).or.lam.lt.lam_grid(1)) then
  print*,'frequency out of range'
  stop
 endif

!interpoliere in T & lam

 do i = 1, 3
  if (tempx.lt.T_grid(i).and.tempx.ge.T_grid(i+1)) then
   do k = 1, 2
    kk = i+k-1
    do j = 1, 20
     if (lam.ge.lam_grid(j).and.lam.lt.lam_grid(j+1)) then
      lam_low = lam_grid(j)
      lam_high = lam_grid(j+1)
      m_real = (REAL(ri(kk,j+1))-REAL(ri(kk,j)))/(LOG(lam_high)-LOG(lam_low))
      b_real = REAL(ri(kk,j))-m_real*LOG(lam_low)
      m_im = (IMAG(ri(kk,j+1))-IMAG(ri(kk,j)))/(LOG(lam_high)-LOG(lam_low))
      b_im = IMAG(ri(kk,j))-m_im*LOG(lam_low)
      ri_real(k) = m_real*LOG(lam) + b_real
      ri_im(k) = m_im*LOG(lam) + b_im
     endif
    enddo   
   enddo
   T_low = T_grid(i+1)
   T_high = T_grid(i)
   m_real = (ri_real(1)-ri_real(2))/(T_high-T_low)
   b_real = ri_real(2)-m_real*T_low
   m_im = (ri_im(1)-ri_im(2))/(T_high-T_low)
   b_im = ri_im(2)-m_im*T_low
   ri_real_final = m_real*tempx + b_real
   ri_im_final = m_im*tempx + b_im
  endif
 enddo 


testre =  ri_real_final**2 -  ri_im_final**2
testim = 2 *  ri_real_final*ri_im_final
!write(*,*) testre, testim
!write(*,*) epsice
!stop
epsice=cmplx(testre, testim)
close(88)


End subroutine watereps

!
!=======================================================================
!
      subroutine refeffect(refre,refim,scmix,epsmat,epsinc)


!vorsicht noch nicht getestet, weil die Übergabe von Tmatrix mit 6 variablen
! nicht funktioniert!
!
!-----------------------------------------------------------------------
!     This routine computes the effective refractive index following the
!     from Maxwell-Garnet.
!
!     epsmat : dielectric constant of the matrix material
!     epsinc : dielectric constant of spherical inclusions
!     scmix  : volume fraction of spherical inclusions 
!-----------------------------------------------------------------------
!
include 'variablen.incf'
!Lokale Parameter: 
complex    :: alpha,beta
real       :: refre, refim, scmix
     
!       write(*,*) scmix, 'scmix'
!       write(*,*) epsmat, 'epsmat'
!       write(*,*) epsinc, 'epsinc'
      ctmp=scmix*(epsinc-epsmat)/(epsinc+2*epsmat)
!    write(*,*) ctmp
      alpha=3*ctmp
      beta=1-ctmp
      epseff=epsmat*(1+(alpha/beta))
      refre=real(csqrt(epseff))
      refim=aimag(csqrt(epseff))
     
  
      return

end subroutine refeffect




      subroutine gauslegquads(x,w,x1,x2,mn,n)
!
!-----------------------------------------------------------------------
!     calculates the abscissas (x) and weights (w) for n-point Gaussian-
!     Legendre quadrature over an integration interval of (x1,x2).
!-----------------------------------------------------------------------
!
      real  w(mn),x(mn),x1,x2, m
      real  crt,p1,p2,p3,pi,pp,srt3,xl,xm,z,z1
      integer j,k,i, mn, n
      data pi,srt3,crt/3.14159265358979e0,5.773502691896256e-1,1.0e-14/

      m=(n+1)/2
      xm=0.5*(x2+x1)
      xl=0.5*(x2-x1)

      if(n.eq.1) then
      x(1)=xm+srt3*xl
      w(1)=1.0
      return
      endif

      do 30 i=1,m

      z=cos(pi*(i-0.25)/(n+0.5))
      k=0

 10   continue
      p1=1.0
      p2=0.0
      do 20 j=1,n
      p3=p2
      p2=p1
      p1=((2*j-1)*z*p2-(j-1)*p3)/j
 20   continue
      pp=n*(z*p1-p2)/(z*z-1)
      z1=z
      z=z1-p1/pp

      k=k+1
      if(abs(z-z1).gt.crt.and.k.le.50)go to 10
      x(i)=xm-xl*z
      x(n+1-i)=xm+xl*z
      w(i)=2.0*xl/((1.0-z*z)*pp*pp)
      w(n+1-i)=w(i)

 30   continue


      return
      end subroutine


      subroutine drvcubspln(x,y,nmax,n,yp1,ypn,y2)
!
!-----------------------------------------------------------------------
!     Given array x and y of length n containing a tabulated function,
!     i.e. y(i)=f(x(i)), with x(1)<x(2)<...<x(n), and given values yp1
!     and ypn for the first derivative at point 1 and n, respectively,
!     this routine returns an array y2 of length n which contains the
!     second derivatives at the tabulated points.
!     if yp1 and/or ypn exceed 1.0e30, the natural boundary values of
!     yp1 and ypn are adapted, i.e., yp1=ypn=0.0
!-----------------------------------------------------------------------
!
      real  un, qn, p, sig, u(99), y2, &
    yp1, ypn, x, y
      integer i, n, nmax
      dimension x(nmax),y(nmax),y2(nmax)

      if ( n.eq.1 ) goto 50
      if (nmax.gt.99) goto 30

      if(yp1.gt.0.999e30) then
        y2(1)=0.0
        u(1)=0.0
      else
        y2(1)=-0.5
        u(1)=3.0/(x(2)-x(1))*((y(2)-y(1))/(x(2)-x(1))-yp1)
      endif

      do 10 i=2,n-1
      sig=(x(i)-x(i-1))/(x(i+1)-x(i))
      p=sig*y2(i-1)+2.0
      y2(i)=(sig-1.0)/p
      u(i)=(6.0*((y(i+1)-y(i))/(x(i+1)-x(i))-(y(i)-y(i-1))/(x(i)-x(i-1)))/(x(i+1)-x(i-1))-sig*u(i-1))/p
 10   continue

      if(ypn.gt.0.999e30) then
        qn=0.0
        un=0.0
      else
        qn=0.5
        un=(3.0/(x(n)-x(n-1)))*(ypn-(y(n)-y(n-1))/(x(n)-x(n-1)))
      endif

      y2(n)=(un-qn*u(n-1))/(qn*y2(n-1)+1.0)
      do 20 i=n-1,1,-1
      y2(i)=y2(i)*y2(i+1)+u(i)
 20   continue
      return
 30   write(6,40)
 40   format(' Subroutine DRVCUBSPLN can not handle array with size'/'larger than 99, change u(99) into a larger array')
      return
 50   y2(1)=0.0


      return
      end subroutine

            subroutine rayaddprc(thinc,sem)

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'
    

      real thinc,sem(7)

!      complex*16 pvz,phy
!      common/glraylei/pvz,phy

      complex*16  sfb(2)
      real*8  safb(4),sabb(4)
      equivalence (safb,sfb)
!
!===> Scattering amplitudes.
!
!      write(*,*) 'rayleigh', thinc
!      write(*,*) pvz, phy
      sfb(1)= phy*cos(thinc)**2+pvz*sin(thinc)**2
      sfb(2)= phy

      sabb(1)= safb(1)
      sabb(2)= safb(2)
      sabb(3)=-safb(3)
      sabb(4)=-safb(4)
!
!===> Elements in extinction and backscattering Mueller matrices.
!
      sem(1)=-safb(2)-safb(4)
      sem(2)=-safb(2)+safb(4)
      sem(3)= safb(1)-safb(3)

      sem(4)= (sabb(1)*sabb(1)+sabb(2)*sabb(2))/2
      sem(5)=-(sabb(3)*sabb(3)+sabb(4)*sabb(4))/2+sem(4)
      sem(4)= 2*sem(4)-sem(5)
      sem(6)=  sabb(1)*sabb(3)+sabb(2)*sabb(4)
      sem(7)=  sabb(2)*sabb(3)-sabb(1)*sabb(4)


      return
      end subroutine

            subroutine tmtaddprc(thinc,sem, nrank)
    
    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

!
      real thinc,sem(7)
!
!      integer, parameter :: m2rank=2*mrank1
    integer :: imd,i1,i2,j1,ia1,mdp,&
        nshft,nsnsft,n1,nmode, nrank
      
   !   complex*16  refrc,eps, tmat
   !   common/partemwv/refrc,eps,wvnm,kxysym
   !   common/gltmatrx/nrank,nmode,tmat(m2rank,m2rank,mrank1)
!
      real*8  degrad,uu,vv,pnmllg(mrank1+1)
      complex*16  c1,cj,abi(m2rank),dnrm(m2rank),fgh(m2rank),fgv(m2rank),&
         sfb(2),sbb(2)
      real*8  safb(4),sabb(4)
      equivalence (safb,sfb),(sabb,sbb)
      data degrad/1.745329251994329d-2/
!
!===> Initialize scattering amplitudes to zero to accumulate with mode.
!
!      write(*,*) wvnm, 'muh'
!stop
      cj=cmplx(0.0,-1.0)
      sfb(1)=0.0
      sfb(2)=0.0
      sbb(1)=0.0
      sbb(2)=0.0
!
!===> Loop for each azimuthal (phi) mode mdp=0, 1, 2, ... ,nrank-1).
!     set shifting indix (nshft) for matrix compression when n<m, m!=0.
!
!      write(*,*) '1', n2mode

      do 70 imd=1,n2mode
!       write(*,*) imd
      mdp=imd-1
      if(n2rank.eq.1) mdp=1
      nshft=mdp-1
      if(nshft.lt.0) nshft=0
      nsnsft=n2rank-nshft
!
!===> normalization factor (=(-)^n*Dmn/wvnm) for each mode.
!
      n1=2*nshft+2
      dnrm(1)=4*(-1)**(nshft+1)*float(n1+1)/float(n1*(nshft+2))/wvnm
      do 10 i=2,n1
      dnrm(1)=dnrm(1)/i
 10   continue
      do 20 i=2,nsnsft
      j=i+nshft
      dnrm(i)=-dnrm(i-1)*((2*j+1)*(j-1)*(j-mdp))/((j+1)*(j+mdp)*(2*j-1))
 20   continue
!
!===> Calculate the incident and scattering field coefficients
!
!-->  evaluate the associated Legendre functions for the incident wave,
!     evaluate incident field coefficients [abi], and compute scattering
!     field coefficients fgv and fgh (=-[T] times [abi]).
!     Here fgv and fgh include normalization factor: (-)^n*Dmn/wvnm.
!
      uu=cos(thinc)
      vv=sin(thinc)
!      write(*,*) 'ah'
!stop
      call tmtgenlgp(uu,vv,pnmllg,mdp,mrank1, n2rank)
!
!      write(*,*) pnmllg, nrank
!stop
      c1=1.
      do 30 i=1,n2rank
      c1=c1*cj
      if(i.le.nshft) goto 30
      ia1=i+1
      i1=i-nshft
      i2=i1+nsnsft
      abi(i1)=(mdp*pnmllg(ia1))*c1
      abi(i2)=(i*uu*pnmllg(ia1)-(i+mdp)*pnmllg(i))*c1
 30   continue
!
      do 50 i=1,nsnsft
      i1=i+nsnsft
      fgv(i )=0.0
      fgv(i1)=0.0
      fgh(i )=0.0
      fgh(i1)=0.0
      do 40 j=1,nsnsft
      j1=j+nsnsft
      fgv(i )=fgv(i )-tmt(i ,j,imd)*abi(j )+cj*tmt(i ,j1,imd)*abi(j1)
      fgv(i1)=fgv(i1)-tmt(i1,j,imd)*abi(j )+cj*tmt(i1,j1,imd)*abi(j1)
      fgh(i )=fgh(i )+tmt(i ,j,imd)*abi(j1)-cj*tmt(i ,j1,imd)*abi(j )
      fgh(i1)=fgh(i1)+tmt(i1,j,imd)*abi(j1)-cj*tmt(i1,j1,imd)*abi(j )
 40   continue
      fgv(i )= fgv(i )*dnrm(i)
      fgv(i1)= fgv(i1)*dnrm(i)
      fgh(i )= fgh(i )*dnrm(i)
      fgh(i1)=-fgh(i1)*dnrm(i)
 50   continue
!
!===> forward and backward scattering amplitudes
!
      i2=(-1)**nshft
      do 60 i=1,nsnsft
      i1=i+nsnsft
      i2=-i2
      sfb(1)=sfb(1)+(-cj*abi(i )*fgv(i)+abi(i1)*fgv(i1))
      sfb(2)=sfb(2)+( cj*abi(i1)*fgh(i)+abi(i )*fgh(i1))
      sbb(1)=sbb(1)+(-cj*abi(i )*fgv(i)-abi(i1)*fgv(i1))*i2
      sbb(2)=sbb(2)+(-cj*abi(i1)*fgh(i)+abi(i )*fgh(i1))*i2
 60   continue

 70   continue
!
!===> Independent elements in the extinction and backscattering Mueller
!     matrices.
!
      sem(1)=-safb(2)-safb(4)
      sem(2)=-safb(2)+safb(4)
      sem(3)= safb(1)-safb(3)

      sem(4)= (sabb(1)*sabb(1)+sabb(2)*sabb(2))/2
      sem(5)=-(sabb(3)*sabb(3)+sabb(4)*sabb(4))/2+sem(4)
      sem(4)= 2*sem(4)-sem(5)
      sem(6)=  sabb(1)*sabb(3)+sabb(2)*sabb(4)
      sem(7)=  sabb(2)*sabb(3)-sabb(1)*sabb(4)


      return
      end subroutine

         subroutine tmtgenlgp(u,v,pnmllg,m,mn, nrank)

    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'


!
!-----------------------------------------------------------------------
!     This routine generates associated Legendre functions (argument u=
!     cos(theta)) of the order from 0 to n divided by v=sin(theta) for a
!     given azimuthal mode m.  Element pnmllg(i)=(i-1)th order function.
!-----------------------------------------------------------------------
!

      integer mn, nrank
      real*8  u,v,pnmllg(mn+1)
      real*8  cnn,pla,plb,plc
    integer :: j1,j2,ibeg,j3,n1,nn

        nn = nrank
      n1=nn+1
    
!
!===> Cawrite(*,*) npnt, n0rankse A): abs(u)=1 or v=0 (theta = 0 or 180 degrees) and m != 1
!     all functions are zero, where those of m=0 are forced to be zero.
!
!      write(*,*) abs(u), m
!      write(*,*) pnmllg
!      stop
      if(abs(u).eq.1.0 .and. m.ne.1) then
      do 10 i=1,n1
      pnmllg(i)=0.0
10    continue
      return
      endif
 !     write(*,*) pnmllg
 !     stop
!
!===> Case B): abs(u)=1 or v=0 (theta = 0 or 180 degrees) and m=1
!     starts with order 0,1,2
!
      if(abs(u).eq.1.0 .and. m.eq.1) then
      pnmllg(1)=0.0
      pla=1.0
      plb=3.0*u
      pnmllg(2)=pla
      pnmllg(3)=plb
      ibeg=4
      goto 30
      endif
!
!===> Case C): general case, i.e., abs(u)<1 or v!=0
!     starts with nn=m since pnmllg(nn<m)=0
!
      cnn=1.0
      do 20 i=1,m
      pnmllg(i)=0.0
      cnn=cnn*(m+i)
 20   continue
     
      pla=cnn*(v/2.0)**m/v
      plb=(2*m+1)*u*pla
      pnmllg(m+1)=pla
      pnmllg(m+2)=plb
      ibeg=m+3
!
!===> Recur upward to obtain all remaining orders of nn<m
!
 30   if(ibeg.gt.n1) return
      j1=2*ibeg-3
      j2=ibeg+m-2
      j3=ibeg-m-1
!      write(*,*) ibeg, n1
      do 40 i=ibeg,n1
      plc=(j1*u*plb-j2*pla)/j3
      pnmllg(i)=plc
      pla=plb
      plb=plc
      j1=j1+2
      j2=j2+1
      j3=j3+1
 40   continue
     
      return
      end subroutine


!-----------------------------------------------------------------

!    

      end module src_hilfe


 MODULE MODI_GAMMA !5
!########################
!
INTERFACE
!

FUNCTION GAMMA(PX)  RESULT(PGAMMA) !10
REAL, INTENT(IN)                                  :: PX
REAL                                              :: PGAMMA
END FUNCTION GAMMA
!
END INTERFACE
!
END MODULE MODI_GAMMA
!     ##################################

      FUNCTION GAMMA(PX)  RESULT(PGAMMA) !10
!     ##################################
!     
!
!!****  *GAMMA * -  Gamma  function  
!!                   
!!
!!    PURPOSE
!!    -------
!       The purpose of this function is to compute the Generalized gamma
!    function of its argument.
!    
!
!!**  METHOD
!!    ------
!!
!!    EXTERNAL
!!    --------
!!      NONE
!!
!!    IMPLICIT ARGUMENTS
!!    ------------------
!!      None
!!
!!    REFERENCE
!!    ---------
!!      Press, Teukolsky, Vetterling and Flannery: Numerical Recipes, 206-207
!!
!!
!!    AUTHOR
!!    ------
!!        Jean-Pierre Pinty *LA/OMP*
!!
!!    MODIFICATIONS
!!    -------------
!!      Original     7/11/95
!
!*       0. DECLARATIONS
!           ------------
!
IMPLICIT NONE
!
!*       0.1 declarations of arguments and result
!
REAL, INTENT(IN)                     :: PX
REAL                                 :: PGAMMA
!
!*       0.2 declarations of local variables
!
INTEGER                              :: JJ ! Loop index
REAL                                 :: ZSER,ZSTP,ZTMP,ZX,ZY,ZCOEF(6)
!
ZCOEF(1) = 76.18009172947146
ZCOEF(2) =-86.50532032941677
ZCOEF(3) = 24.01409824083091
ZCOEF(4) = -1.231739572450155
ZCOEF(5) =  0.1208650973866179E-2
ZCOEF(6) = -0.5395239384953E-5
ZSTP     =  2.5066282746310005
!
ZX = PX
ZY = ZX
ZTMP =  ZX + 5.5
ZTMP = (ZX + 0.5)*ALOG(ZTMP) - ZTMP
ZSER = 1.000000000190015
!
DO JJ = 1 , 6
  ZY = ZY + 1.0
  ZSER = ZSER + ZCOEF(JJ)/ZY
END DO
!
PGAMMA = EXP( ZTMP + ALOG( ZSTP*ZSER/ZX ) )
RETURN
!
END FUNCTION GAMMA
