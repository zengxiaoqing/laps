      module src_tmatrix

      use src_hilfe
      use src_global
      contains



       subroutine tmtbessel(nn,npnt)

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'
!
!-----------------------------------------------------------------------
!     This routine generates the ratios and products of spherical Bessel
!     (bj), Neumann (by), and first kind Hankel (bh) functions of the
!     order from 0 to n, of real argument x (in bjx, byx, bhx) and/or
!     complex*16 argument z (in bjz).  The ratios and products are defined
!     as: r_i=b_(i-1)/b_i and p_ij=(b1_i)*(b2_j), i,j=1,2,...,n.
!     Currently, rbhx,rbjx,rbjz,pjzhx,pjzjx are of particular interest.
!     there are npnt number of points at which functions are evaluated.
!-----------------------------------------------------------------------
!
!      integer, parameter::mpnt=32
!
    integer :: nn, ia1,is1,js1,ipnt,i1, nst, npnt
 !      integer :: nrank

!      common/surfgrid/u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
!      common/bessfunc/pjzhx(mrank1,mrank1,mpnt),pjzjx(mrank1,mrank1,mpnt),&
!           rbjz(mrank1,mpnt),rbhx(mrank1,mpnt),rbjx(mrank1,mpnt)
!      common/partemwv/refrc,eps,wvnm,kxysym
!
      real*8  bjx0,bjx1,byx0,byx1,tx(2),x,xhi,xlo
      complex*16  bjz0,bjz1,tz(2),z
!      logical*1 lhi(8),llo(8)
      real*8 lhi(8),llo(8)
      equivalence (xhi,lhi),(xlo,llo)
      data llo,lhi/0,0,0,1,0,0,0,0, 127,127,255,0,0,0,0,0/
!
!===> Loop for each integration points
!
      do 60 ipnt=1,npnt            !

!         print*,xlo,'hhhh'

      x=xkr(ipnt)    !
      z=refrc*x    !


!
!===> rbjx and rbjz array calculations.
!
!-->  the t array is used to recur down to the highest order function
!     that are needed.  first set the nozero starting order values for
!     downward recursion.
!
      nst=nn+10*(1+abs(z))        !nn=nrank laufindex in tmtnrkconv.
!    write(*,*) nn
      tx(1)=dmax1((2*nst+1)/x,xlo)
      tz(1)=dmax1((2*nst+1)/abs(z),xlo)
!
!-->  recur downward to obtain orders nn.
!
      do 10 i=nst-1,nn,-1
      i1=2*i+1
      tx(2)=i1/x-1/tx(1)
      tx(1)=tx(2)
      tz(2)=i1/z-1/tz(1)
      tz(1)=tz(2)
 10   continue
      rbjx(nn,ipnt)=tx(2)
      rbjz(nn,ipnt)=tz(2)

!
!-->  continue downward recursion to order 1.
!
      do 20 i=nn-1,1,-1
      ia1=i+1
      i1=ia1+i
      rbjx(i,ipnt)=i1/x-1/rbjx(ia1,ipnt)
      rbjz(i,ipnt)=i1/z-1/rbjz(ia1,ipnt)
 20   continue
!
!===> rbhx, pjzhx and pjzjx array calculations.
!
!-->  obtain the zeroth and first order functions, then construct the
!     first order ratio and product functions.
!
      bjx0= sin(x)/x
      byx0=-cos(x)/x
      bjz0= sin(z)/z
      bjx1=bjx0/x+byx0
      byx1=byx0/x-bjx0
      bjz1=(bjz0-cos(z))/z

!     if(abs(rbjx(1,ipnt)/(bjx0/bjx1)-1).ge.1.e-9) write(6,30) 'rbjx'
!     if(abs(rbjz(1,ipnt)/(bjz0/bjz1)-1).ge.1.e-9) write(6,30) 'rbjz'
!30   format(2x, 'error occured while evaluating ',a4)
      rbhx(1,ipnt)=dcmplx(bjx0,-byx0)/dcmplx(bjx1,-byx1)
      pjzhx(1,1,ipnt)=bjz1*dcmplx(bjx1,-byx1)
   !   write(*,*) pjzhx(1,1,ipnt), bjz1, bjz0, z
    !  write(*,*)  x in ordnung
!      write(*,*) refrc
      pjzjx(1,1,ipnt)=bjz1*bjx1

!
!-->  recur upward to order nn.
!
      do 50 i=2,nn

      is1=i-1
!    write(*,*) rbhx(i,ipnt)
      rbhx(i,ipnt)=1/((2*i-1)/x-rbhx(is1,ipnt))
!    write(*,*) pjzhx(is1,1,ipnt)        ! in ordnung

      pjzhx(i,1,ipnt)=pjzhx(is1,1,ipnt)/rbjz(i,ipnt) !IN ORDNUGN

      pjzhx(1,i,ipnt)=pjzhx(1,is1,ipnt)/rbhx(i,ipnt)
      pjzjx(i,1,ipnt)=pjzjx(is1,1,ipnt)/rbjz(i,ipnt)
      pjzjx(1,i,ipnt)=pjzjx(1,is1,ipnt)/rbjx(i,ipnt)
      do 40 j=2,nn
      js1=j-1
      pjzhx(j,i,ipnt)=pjzhx(js1,i,ipnt)/rbjz(j,ipnt)    !nicht in ordnung
                            !rbjz in ordnung

      pjzjx(j,i,ipnt)=pjzjx(js1,i,ipnt)/rbjz(j,ipnt)
 40   continue
 50   continue

 60   continue


      return
      end subroutine

!***************************************************
        subroutine tmtbuildup(npnt,aovrb)
!***************************************************

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

!      integer, parameter::mpnt=32
!      real*8 dkr,u,v,wt,xkr,
      real :: wcnzc2
!      complex*16  refrc,eps
    integer :: npnt, ipnt
 !     common/surfgrid/u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
 !     common/partemwv/refrc,eps,wvnm,kxysym
!      common/tmatbody/scshp,deq,dmx,wcnxl,wcnzc,wcnlm
!      common/mist/aovrb
      real*8  degrad,dsa,r1,r2,r3,r4,ulow,uupp
      data degrad/1.745329251994329d-2/
    include 'constants.incf'
!**************************************************
!
!===> Set execution parameters into SI units.
!
!***********************************************!


      iscshp=scshp+0.1    !=1

      if(iscshp.eq.1) kxysym=1
      if(iscshp.eq.2) kxysym=0        !kxysym = 1

      if(kxysym.eq.1) npnt=(npnt+1)/2    !npnt = 16
!*************************************************
!
!===> Set the quadrature integration points over the surface of the
!     object and the weighting values.  u=cos(theta) is the integration
!     variable.  if particle is symmytric about xy plane, only the upper
!     half space is needed.
!
      ulow=kxysym-1
      uupp=cos(0.0)        !u, wt, ulow, uupp, mpnt, npnt

!************************************************
      call tmtgauslegquad(u,wt,ulow,uupp,mpnt,npnt)
!************************************************


      do ipnt=1,npnt
      v(ipnt)=sqrt(1.0-u(ipnt)**2)    !v(ipnt)
      if(kxysym.eq.1) wt(ipnt)=2*wt(ipnt)    !wt(ipnt)
      enddo
!*************************************************

!
!===> Prepare radius-related-only quantities that are used repeatedly:
!     kr=wvnm*radius and its derivative of theta times sin(theta)
!     (denoted as xkr and dkr) at each integration point.
!     Note: here dsa==particle dimension on z-axis.
!
!=>   spheroids.  a=b: sphere, a<b:oblate, a>b: prolate, a: z-dimension.
!


      if(iscshp.eq.1) then
      if(deq.gt.0.0) then
      dsa=deq*aovrb**(2.d0/3.d0)    !

      else
      dsa=dmx*min(1.d0,aovrb)    !
      endif
     !  write(*,*) dsa

      r1=wvnm*dsa/2    !
!**********************************************************************
      do 10 ipnt=1,npnt
      r2=sqrt(u(ipnt)**2+(v(ipnt)*aovrb)**2)    !

      xkr(ipnt)=r1/r2    !

      dkr(ipnt)=r1*u(ipnt)*v(ipnt)**2*(1.0-aovrb**2)/r2**3    !

 10   continue
      endif        !r2
!*******************************************************************
!
!=>   P. K. Wang's conical particle.
!     wcnzc=z-dimension, wcnxl=x-dimension, wcnlm=shape factor>1.
!    schleife wird nicht betreten fuer regen, aber fuer graupel


      if(iscshp.eq.2) then
      wcnzc2=wcnzc/2.0

      r1=-0.9
      do 24 ipnt=1,npnt
      r2=pi*wcnzc2*v(ipnt)/u(ipnt)/wcnxl

      do 20 i=1,100
      r3=sqrt(1.0d0-r1*r1)
      r4=r1*acos(r1/wcnlm)+r3*r3/sqrt(wcnlm*wcnlm-r1*r1)+r2*r3
      xkr(ipnt)=r1+(r3*acos(r1/wcnlm)-r1*r2)*r3/r4

      if(abs(1.d0-xkr(ipnt)/r1).le.1.d-9) goto 22
      r1=xkr(ipnt)

      if(abs(r1).gt.(1.0d0-1.d-14)) r1=(1.d0-1.d-14)*r1/abs(r1)
 20   continue
 22   dkr(ipnt)=wvnm*wcnzc2*r1*(v(ipnt)/u(ipnt))**2*(1.0d0-r2*r3/(v(ipnt)**2*r4))
      xkr(ipnt)=wvnm*wcnzc2*xkr(ipnt)/u(ipnt)

 24   continue
      r1=1.d0/wcnlm/wcnlm
      deq=(6*(wcnxl/pi)**2*wcnzc2*(3.2889+(0.2667+(0.0382+(0.0113+(0.0046+&
         0.0024*r1)*r1)*r1)*r1)*r1))**(1.0d0/3.0d0)
      r2=-0.4/wcnlm
      do 26 i=1,100
      r3=sqrt(1.0d0-r2*r2)
      r4=r2-(r2*acos(r2)*r3+r1-r2*r2)/(cos(r2)*(1-2*r2*r2)/r3-3*r2)
      if(abs(1.d0-r4/r2).le.1.d-9) goto 28
      r2=r4
 26   continue
 28   aovrb=wcnzc2/(wcnxl/pi*sqrt(1-(r4*wcnlm)**2)*acos(r4))
      dmx=wcnzc
      if(aovrb.lt.1.0) dmx=dmx/aovrb
      endif


      return
      end subroutine

        subroutine tmtdriver(n0rank,n0pnt,crtnm, aovrb, nrank)

    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'


    integer ::nmode, imd, n0rank, n0pnt,nshft, mdp,&
          kmdcon, n2snsft, npnt, nrank
!
!      common/glhydbdy/glbdy(12),glcad(2,6)
!
    real crtnm
!      integer, parameter :: mpnt=32
!      real*8  wvnm
!      complex*16 tmt
!      common/partemwv/refrc,eps,wvnm,kxysym
!      common/abtmtrxs/tmat(m2rank,m2rank),a(m2rank,m2rank)
!      common/gltmatrx/n2rank,n2mode,tmt(m2rank,m2rank,mrank1)
!      common/tmatbody/scshp,deq,dmx,wcnxl,wcnzc,wcnlm
!
    real  freq0
      real*8  ssmd1
      complex*16 svvnm,shhnm!, tmat(m2rank, m2rank)
      real*8  degrad,pi,vlght
      data vlght,pi,degrad/3.d8,3.14159265358979d0,1.745329251994329d-2/
!
!===> Set the quadrature integration points over the surface of the
!     object and the weighting values.  u=cos(theta) is the integration
!     variable.  if particle is symmytric about xy plane, only the upper
!     half space is needed.  Note: change units to SI units.
!
!       write(*,*) 'test'
      refrc = cmplx(refre,refim)    !
!      scshp = glbdy(6)    !
      deq   = deq *0.001    !
      dmx   = dmx *0.001    !
 !     aovrb = glbdy(9)    !
!    write(*,*) aovrb
      wcnxl = wcnxl *0.001    !
      wcnzc = wcnzc *0.001    !
!      wcnlm = glbdy(12)    !
    freq0=freq*1.e9    !

      eps=refrc**2    !
      wvnm=2*pi*freq0/vlght    !
!**********************************************************************

      npnt=min0(n0pnt,mpnt)    !starting tmatrix rank, automatically increases                 !with size, bis jetzt konstant 5, VORSICHT!!

      n0rank=max0(n0rank,5)
      crtnm=amin1(crtnm,1.0e-2)

      !passt alles

!***********************************************************************
      call tmtbuildup(npnt,aovrb)    !getestest fuer regen
!***********************************************************************
!
!===> Decide convergent nrank value using special property at thinc=0.
!     the ratios and products of spherical Bessel functions and Hankel
!     functions at each integration point, up to the order of nrank,
!     and multiplied with proper weights are computed here.
!
     nrank =0
      call tmtnrkconvrg(n0rank,nrank, npnt,crtnm)!, tmat)
!**********************************************************************
!         write(*,*) 'hallo'
!
!===> Initialize scattering amplitudes to zero to accumulate with mode.
!
      svvnm=0
      shhnm=0
!
!===> Loop for each azimuthal (phi) mode mdp=0, 1, 2, ... ,nrank-1).
!
      ssmd1=0.0
      kmdcon=0
      do 70 imd=1,nrank
!         write(*,*) imd, nrank
      mdp=imd-1
      if(nrank.eq.1) mdp=1
!
!-->  integrate matrices A and B over the surface (integration variable
!     u=cos(theta)) and calculate the Transition-matrix T=B*A^(-1).

      call tmtmtrxabt(nrank,mdp,npnt,tmat,aaa)



!
!-->  layout T-matrix
!
      call tmtnmdconvrg(nrank, mdp, svvnm,shhnm,ssmd1,crtnm,kmdcon)!, tmat)
!      write(*,*) kmdcon, 'kmdcon'
      if(kmdcon.eq.1) goto 80

      nshft=mdp-1
      if(nshft.lt.0) nshft=0
      n2snsft=2*(nrank-nshft)
      do 65 j=1,n2snsft
      do 65 i=1,n2snsft
      tmt(i,j,imd)=tmat(i,j)
 65   continue

 70   continue

 80   nmode=imd-1

      n2rank=nrank
      n2mode=nmode
!      write(*,*) n2mode, nmode, imd
!stop
      return
      end subroutine

      subroutine tmtgauslegquad(x,w,xx1,xx2,mn,nn)

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'
!
!-----------------------------------------------------------------------
!     calculates the abscissas (x) and weights (w) for n-point Gaussian-
!     Legendre quadrature over an integration interval of (xx1,xx2).
!-----------------------------------------------------------------------
!
      integer mn
      real*8  w(mn),x(mn),xx1,xx2
      real*8  crt,p1,p2,p3,pi,pp,srt3,xl,xm,z,zz1
    integer :: nn
      data pi,srt3,crt/3.14159265358979d0,5.773502691896256d-1,1.0d-14/
  !     write(*,*) 'hallo'
      m=(nn+1)/2
      xm=0.5*(xx2+xx1)
      xl=0.5*(xx2-xx1)

      if(nn.eq.1) then
      x(1)=xm+srt3*xl
       w(1)=1.0
      return
      endif
!************************************************************************
      do 30 i=1,m

      z=cos(pi*(i-0.25)/(nn+0.5))
      k=0

 10   continue
      p1=1.0
      p2=0.0
      do 20 j=1,nn
      p3=p2
      p2=p1
      p1=((2*j-1)*z*p2-(j-1)*p3)/j
 20   continue
      pp=nn*(z*p1-p2)/(z*z-1)
      zz1=z
      z=zz1-p1/pp

      k=k+1
      if(abs(z-zz1).gt.crt.and.k.le.50)go to 10
      x(i)=xm-xl*z
      x(nn+1-i)=xm+xl*z
      w(i)=2.0*xl/((1.0-z*z)*pp*pp)
      w(nn+1-i)=w(i)

 30   continue


      return
      end subroutine


      subroutine tmtmtrxabt(nrank,mdp,npnt,bb,aa)

    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

!      integer, parameter :: mpnt=32
!      integer, parameter :: m2rank=2*mrank1
    integer :: mdp, npnt, nshft, msq,lqm, nsnsft,n2sft

    integer :: jam,nrank,  jtm, ja1, j2, j1, itj, jtja1, &
        iaj, i2, i1, ia1,itia1, iam, itm, ipnt

      complex*16 aa(m2rank,m2rank),bb(m2rank,m2rank)
  !    common/surfgrid/u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
  !    common/bessfunc/pjzhx(mrank1,mrank1,mpnt),pjzjx(mrank1,mrank1,mpnt),&
  !        rbjz(mrank1,mpnt),rbhx(mrank1,mpnt),rbjx(mrank1,mpnt)
!      common/partemwv/refrc,eps,wvnm,kxysym
    !  common/abtmtrxs/b(m2rank,m2rank),a(m2rank,m2rank)

      real*8  bil,bjk,bjxr,pnmllg(mrank1+1),pp00,pp01,pp10,pp11,qtusq,&
           qty,r01,r1,r2,uu,vv,x,y
     complex*16  ail1,ail2,ajk1,ajk2,akj1,akj2,ali1,ali2,bhxjz,bhxr,&
           bjxjz,bjzr,c01,c1,c2,cil1,cil2,cjk,ckj

!
!===> Set shifting indix (nshft) for matrix compression when n<m, m!=0.
!     and initialize matrices A and B to zero for surface integration.
!

     nsnsft = 0
!       write(*,*) 'test', nrank
      nshft=mdp-1    !0
      if(nshft.lt.0) nshft=0
      nsnsft=nrank-nshft    !5

      n2sft=2*nsnsft    !10


      aa = 0.0
      bb = 0.0

!
!===> Integrating over the surface (integration variable: u=cos(theta))
!     To reduce computer time, frequently used variable combinations are
!     sorted and calculated a priori.   a:+, s:-, t:*, d:/
!
     lqm=1
     if(mdp.eq.0) lqm=2    !1

     do 50 ipnt=1,npnt        !1=>16

!
!-->  calculate associated Legendre functions at each integration point.
!
      uu=u(ipnt)    !!
      vv=v(ipnt)    !!

      call tmtgenlgp(uu,vv,pnmllg,mdp,mrank1, nrank)
!        write(*,*) 'xxx', nrank
    msq=mdp*mdp    !1
      x=xkr(ipnt)    !

      y=dkr(ipnt)    !

      qty=lqm*y    !

      qtusq=lqm*uu*uu    !

!
!-->  loop for each row (index i) of matrices A(row,col) and B(row,col).
!
     do 40 i=1,nrank
!       write(*,*) i, nrank
     if(i.le.nshft) goto 40
     iam=i+mdp    !
     itm=i*mdp    !
     ia1=i+1    !

     itia1=i*ia1    !
     i1=i-nshft        !
     i2=i1+nsnsft    !

     bhxr=rbhx(i,ipnt)    !
     bjxr=rbjx(i,ipnt)    !
     c01=x*bhxr-i    !
     r01=x*bjxr-i    !

!
!-->  loop for each column (index j) of matrices A(i,j) and B(i,j).
!
      do 30 j=1,nrank
!       write(*,*) j, nrank
     if(j.le.nshft) goto 30
     jam=j+mdp        !
     jtm=j*mdp        !
     ja1=j+1        !
     jtja1=j*ja1    !

     iaj=i+j        !
     itj=i*j        !
     j1=j-nshft        !
     j2=j1+nsnsft    !

     bjzr=rbjz(j,ipnt)    !
     bhxjz=pjzhx(j,i,ipnt)    !
!    write(*,*) pjzhx(j,i,ipnt)
     bjxjz=pjzjx(j,i,ipnt)    !
!     write(*,*) pjzjx(j,i,ipnt)
!
     pp00=pnmllg(i  )*pnmllg(j  )    !
     pp01=pnmllg(i  )*pnmllg(ja1)    !
     pp10=pnmllg(ia1)*pnmllg(j  )    !
     pp11=pnmllg(ia1)*pnmllg(ja1)    !

!-->  parameter bil (IML,LMI,IML',LMI'), and bjk (JMK,KMJ,JMK',KMJ')
     bil=x*( iaj*uu*pp11-iam*pp01-jam*pp10 )    !
     bjk=x*(lqm*(iam*jam*pp00-((itj+itm)*pp10+(itj+jtm)*pp01)*uu)+(msq+itj*qtusq)*pp11 )    !


!
!-->  calculate the [A] and [B] submatrices.  Using symmetry property of
!     associated Legendre function for simplification to xy-plane mirror
!     symmetric particle: I=L=0 if (i+j) is even, J=K=0 if (i+j) is odd.
!
!-->  increment the -(I+refrc*L) and (L+refrc*I) submatrices.
!     when mdp(or m)=0, Iomnemn=Lomnemn=0.
!
      if(mdp.eq.0) goto 20

      if(kxysym.eq.1 .and. mod(iaj,2).eq.0) goto 20

!
!-->  cil1 (IML, LMI), cil2 (IML',LMI'), and ail1,ali1,ail2,ali2
     r1=pp11*y
      c1=itia1*bjzr-itj*(ia1+ja1)/x
      cil1=r1*(c1+jtja1*bhxr)
      cil2=r1*(c1+jtja1*bjxr)
      c1=x*bjzr-j
      c2=x-i*bjzr+itj/x
      ail1=bhxr*c1+c2
      ali1=ail1+(eps-1)*x
      ail2=bjxr*c1+c2
      ali2=ail2+(eps-1)*x
      aa(i2,j1)=aa(i2,j1)+mdp*bhxjz*(ail1*bil+cil1)
      aa(i1,j2)=aa(i1,j2)-mdp*bhxjz*(ali1*bil+cil1)/refrc

      bb(i2,j1)=bb(i2,j1)+mdp*bjxjz*(ail2*bil+cil2)
      bb(i1,j2)=bb(i1,j2)-mdp*bjxjz*(ali2*bil+cil2)/refrc
!       write(*,*) i2, j1, bb(i2,j1)
!       write(*,*) i1,j2, bb(i1,j2)
!-->  increment the (K+refrc*J) and (J+refrc*K) submatrices.
!
20   if(kxysym.eq.1 .and. mod(iaj,2).eq.1) goto 30

!
!-->  cjk (JMK,JMK'), ckj (KMJ,KMJ'), and akj1,ajk1,akj2,ajk2
      r1=jtja1*(i*uu*pp11-iam*pp01)    !

      r2=itia1*(j*uu*pp11-jam*pp10)    !
      cjk=(r1-eps*r2)*qty        !
      ckj=(r1-r2)*qty            !
      c1=x*bjzr-j            !

      ajk1=c1-eps*c01            !
      ajk2=c1-eps*r01            !
      akj1=c1-c01            !
      akj2=c1-r01            !

      aa(i1,j1)=aa(i1,j1)+bhxjz*(akj1*bjk+ckj)        !
      aa(i2,j2)=aa(i2,j2)+bhxjz*(ajk1*bjk+cjk)/refrc    !

      bb(i1,j1)=bb(i1,j1)+bjxjz*(akj2*bjk+ckj)        !
      bb(i2,j2)=bb(i2,j2)+bjxjz*(ajk2*bjk+cjk)/refrc    !
!      write(*,*) i1,j1, bb(i1,j1)
!      write(*,*) i2,j2, bb(i2,j2)
!      write(*,*) bb(i1,j1)

 30   continue

 40   continue

 50   continue
!      write(*,*) 'tst'
!
!===> Calculate the Transition-matrix T=B*A^(-1) (stored in B).
!write(*,*) bb(1,1), bb(3,2),bb(5,9)

      call tmtprcssm(n2sft,aa,bb)
!write(*,*) 'aus'
!stop

      return
      end subroutine


       subroutine tmtnmdconvrg(nrank, mdp,svvnm,shhnm,ssmd1,crtnm,kmdcon)!, tmat)

    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

      real*8 ssmd1
      complex*16 svvnm,shhnm

      integer, parameter:: mpnt=32
    integer :: j1,i1,i2, mdp, ia1, n1,nshft, nsnsft,&
        kmdcon, nrank
    real :: crtnm
!     real*8  wvnm
!      complex*16  refrc,eps
!      common/partemwv/refrc,eps,wvnm,kxysym
!!      common/abtmtrxs/tmat(m2rank,m2rank),a(m2rank,m2rank)
!   complex*16 :: tmat(m2rank, m2rank)

     real*8  ssmd2,uu,vv,pnmllg(mrank1+1)
      complex*16  c1,cj,abi(m2rank),dnrm(m2rank),fgh(m2rank),fgv(m2rank)
!
!===> Set shifting indix (nshft) for matrix compression when n<m, m!=0.
!
!    write(*,*) mdp, nrank

      cj=cmplx(0.0,-1.0)
      nshft=mdp-1
      if(nshft.lt.0) nshft=0
      nsnsft=nrank-nshft
!      write(*,*) 'nsnsft, nrank, nshft'
!       write(*,*) nsnsft, nrank, nshft
!
!===> normalization factor (=(-)^n*Dmn/wvnm) for each mode.
!
      n1=2*nshft+2
      dnrm(1)=4*(-1)**(nshft+1)*float(n1+1)/float(n1*(nshft+2))/wvnm
!     write(*,*) dnrm(1), nshft, n1, wvnm
!top
      do i=2,n1
      dnrm(1)=dnrm(1)/i
      enddo
!    write(*,*) dnrm(1)        !
      do i=2,nsnsft
!    write(*,*) nsnsft    !falsch
      j=i+nshft
      dnrm(i)=-dnrm(i-1)*((2*j+1)*(j-1)*(j-mdp))/((j+1)*(j+mdp)*(2*j-1))
!     write(*,*) dnrm(i)
    enddo
!top
!
!===> Determin the largest nmode value needed, using back-scattering
!     alone x axis
!
      uu=0.0
      vv=1.0
      call tmtgenlgp(uu,vv,pnmllg,mdp,mrank1, nrank)

      c1=1.
      do 33 i=1,nrank
      c1=c1*cj
      if(i.le.nshft) goto 33
     ia1=i+1
      i1=i-nshft
      i2=i1+nsnsft
      abi(i1)=(mdp*pnmllg(ia1))*c1
      abi(i2)=(-(i+mdp)*pnmllg(i))*c1
 33   continue

      do 34 i=1,nsnsft
      i1=i+nsnsft
      fgv(i )=0.0
      fgv(i1)=0.0
      fgh(i )=0.0
      fgh(i1)=0.0
      do j=1,nsnsft
      j1=j+nsnsft
      fgv(i )=fgv(i )-tmat(i ,j)*abi(j )+cj*tmat(i ,j1)*abi(j1)
      fgv(i1)=fgv(i1)-tmat(i1,j)*abi(j )+cj*tmat(i1,j1)*abi(j1)
      fgh(i )=fgh(i )+tmat(i ,j)*abi(j1)-cj*tmat(i ,j1)*abi(j )
      fgh(i1)=fgh(i1)+tmat(i1,j)*abi(j1)-cj*tmat(i1,j1)*abi(j )
      enddo
      fgv(i )= fgv(i )*dnrm(i)
      fgv(i1)= fgv(i1)*dnrm(i)
      fgh(i )= fgh(i )*dnrm(i)
      fgh(i1)=-fgh(i1)*dnrm(i)
 34   continue
      i2=(-1)**nshft
      do i=1,nsnsft
      i1=i+nsnsft
      i2=-i2
!      write(*,*) fgv(i1), dnrm(i)
      svvnm=svvnm+(-cj*abi(i )*fgv(i)-abi(i1)*fgv(i1))*i2
      shhnm=shhnm+(-cj*abi(i1)*fgh(i)+abi(i )*fgh(i1))*i2
      enddo
      ssmd2=abs(svvnm)+abs(shhnm)
!     write(*,*) ssmd2, svvnm, shhnm
!top
      if(abs(1.0-ssmd1/ssmd2).le.crtnm) then
      kmdcon=1
      goto 100
      endif
      ssmd1=ssmd2


 100  return
      end subroutine


       subroutine tmtnrkconvrg(n0rank,nrank,npnt,crtnm)!, tmat)

    implicit none

    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'


!
!-----------------------------------------------------------------------
!     This routine determines the largest nrank value needed by current
!     environment for convergent solution, defined by critirion crtnm.
!     It uses special property of back-scattering along symmetry axis.
!-----------------------------------------------------------------------
!
!      integer, parameter :: mpnt=32
!
    integer :: i1,j1,ipnt,mdp, n0rank, npnt, nrank
!      real*8  dkr,xkr, rbjx, wvnm
    complex*16 svvnm!,tmat(m2rank,m2rank)
!      complex*16 pjzhx,pjzjx,rbjz,rbhx
!      common/surfgrid/u(mpnt),v(mpnt),wt(mpnt),xkr(mpnt),dkr(mpnt)
!      common/bessfunc/pjzhx(mrank1,mrank1,mpnt),pjzjx(mrank1,mrank1,mpnt),&
!           rbjz(mrank1,mpnt),rbhx(mrank1,mpnt),rbjx(mrank1,mpnt)
!      common/partemwv/refrc,eps,wvnm,kxysym
!      common/abtmtrxs/tmat(m2rank,m2rank),a(m2rank,m2rank)
    real :: CRTNM
      real*8  r1,r2,ssnr1,ssnr2,x
      complex*16 c1,c2,cj
!       write(*,*) 'hallo'
      cj=cmplx(0.0,-1.0)
      mdp=1
      ssnr1=0.
      ssnr2=0.
      nrank = n0rank
      r2 = 0
      c2 = 0
      tmat = 0.
!      write(*,*) 'hallo', nrank
!**************************************************************
    x=1.


    mschleife: do while (x>crtnm)


    ssnr1=ssnr2
!        write(*,*) '1'
    call tmtbessel(nrank,npnt)
!         write(*,*) '2'
        do 10 ipnt=1,npnt
        do 10 i=1,nrank
        rbjz(i,ipnt)=rbjz(i,ipnt)*refrc

        do 10 j=1,nrank

        pjzhx(j,i,ipnt)=pjzhx(j,i,ipnt)*wt(ipnt)    !
!    write(*,*) pjzhx(j,i,ipnt) !nochmal geprueft
    pjzjx(j,i,ipnt)=pjzjx(j,i,ipnt)*wt(ipnt)    !
 10     continue
!    write(*,*) nrank, mdp, npnt, mrank1
      call tmtmtrxabt(nrank,mdp,npnt,tmat,aaa)
!     write(*,*) tmat(1,1), tmat(10,6)
!stop
!        write(*,*) '4'
!
       svvnm=0

      c1=0.5/wvnm    !

   do 20 i=1,nrank
      i1=i+nrank    !
      r1=(2*i+1.0)/(i*(i+1))    !
      c1=c1*cj    !
      c2=c1

      do 20 j=1,nrank
!    write(*,*) j
      j1=j+nrank
      r2=j*(j+1)*r1
      c2=c2*cj

    ! tmat in ordnung !!!
      svvnm=svvnm-r2*c2*(tmat(i1,j)+tmat(i,j1)-cj*(tmat(i1,j1)-tmat(i,j)))


 20   continue




       ssnr2=abs(svvnm)

      x=abs(1-ssnr1/ssnr2)
!      write(*,*) nrank
!       write(*,*) tmat(i1,j), i1, j
!      write(*,*) x, ssnr1, ssnr2, svvnm, crtnm
      nrank = nrank+1
      end do mschleife
   !   if(abs(1-ssnr1/ssnr2).le.crtnm) goto 40
   !   ssnr1=ssnr2

!30 continue
     !**************************************************
 40   if(nrank.ge.mrank1) write(6,*) 'Warning: '&
     &'T-matrix may not be convergent, try again with larger mrank1.'

      nrank=nrank-2    !geändert!!!
      n0rank=nrank-1

      return
      end subroutine


            subroutine tmtprcssm(nn,as,bs)
    implicit none


    INCLUDE 'parameter.incf'
    Include 'fields.incf'
    include 'variablen.incf'

!
!-----------------------------------------------------------------------
!     This routine calculates matrix production T=B*A^(-1) by means of
!     Gauss-Jordan column reduction which is equivalent to row reduction
!     applied to the transpose T'=A'^(-1)*B'.
!     The output matrix T(row,col) is stored in the final B(row,col).
!-----------------------------------------------------------------------
!
!      integer, parameter::mpnt=32
!
    integer :: nn, imax, jj
      complex*16 bs(m2rank,m2rank),as(m2rank,m2rank)
!      common/abtmtrxs/b(m2rank,m2rank),a(m2rank,m2rank)
      complex*16 ajimax,arat
      integer :: ls(m2rank)
!
!===> Start column reduction of A, or equivalently, row reduction of A'.
!     indix notations:  i=row of A and B; j and k=column of A and B.

      ls = 0
!      write(*,*) ls

      j = 0
       i =0
      imax = 0
      ajimax = 0

      do 50 j=1,nn   !nn=10
!    write(*,*) j, 'j'
!        write(*,*) ls
!
!-->  search for the element of maximum magnitude in the jth column of A
!
      ajimax=as(1,j)

      imax=1

      do 10 i=2,nn
      if(abs(as(i,j)).gt.abs(ajimax)) then
      ajimax=as(i,j)
      imax=i
      endif

 10   continue

!  write(*,*) 'test'
!  write(*,*) ls
!
!-->  normalize the jth column of A and B by ajimax.
!
      do 20 i=1,nn
      as(i,j)=as(i,j)/ajimax
      bs(i,j)=bs(i,j)/ajimax


 20   continue

!
!-->  use column transformations to obtain zeros for the imax-th row but
!     element of the jth column of A.  the same is applied to B.
!
      k = 0
      do 40 k=1,nn

      if(k.ne.j) then

      arat=-as(imax,k)
      do 30 i=1,nn
      if(abs(as(i,j)).gt.0.0) as(i,k)=arat*as(i,j)+as(i,k)
      if(abs(bs(i,j)).gt.0.0) bs(i,k)=arat*bs(i,j)+bs(i,k)
 30   continue
      as(imax,k)=0.0

      endif
 !  write(*,*) bs(i,k)
 40   continue
!
!-->  store column counter of A in array ls such that ls contains the
!     location of the pivot element of each row (after reduction)
!
!      write(*,*) ls
      ls(imax) = 0
      ls(imax)= j

!       write(*,*) ls(imax), j, imax
!       write(*,*) ls
 50   continue
!      write(*,*) ls, 'ls'
!      write(*,*) nn


!
!===> perform column interchanges on A as indicated in array ls.
!     put k=ls(k), if k<i, that column has already been involved in an
!     interchange so iterate k=ls(k) until k>i (corresponding to a
!     column stored above the jth column) is obtained.
!
!      write(*,*) nn, 'aus'
!      write(*,*) ls
      do 80 jj=1,nn
!         write(*,*) jj
      k=jj
!     write(*,*) k
 60   k=ls(k)
!   !  write(*,*) k
!      write(*,*) 'k', k,'ls',  ls(k), 'j', jj
     if(k .lt. jj) goto 60
!        write(*,*) 'k', k,'ls',  ls(k), 'j', jj
      if(k.gt.0) then
      do 70 i=1,nn
      arat=bs(i,jj)
      bs(i,jj)=bs(i,k)
      bs(i,k)=arat
 !     write(*,*) bs(i,jj), bs(i,k)
 70   continue
      endif!

 80   continue
!write(*,*) bs
!stop
      return
      end subroutine

      end module src_tmatrix
