c
c
      subroutine kalman(dta,F,Y,p,it,w,v,x,xt,imax,m,atime)
c
c*********************************************************************
c
c     Routine to apply Kalman Filter to a set of obs for QC purposes.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c
c     Notes:
c          Variables "m" and "im,jm" must be the same size
c
c*********************************************************************
c
      parameter (badflag=-99.9)
      parameter (im=500,jm=500)
      Real K(im,im),P(m,m),F(m,m),II(im,im)
      Real H(im,im),HT(im,im),FT(im,im),G(im)
      Real X(m),Y(m),XT(m)
      Real PT(im,im),ZZ(im,im)
      Real A(im,im),B(im),C(im),Z(im),D(im,im),E(im,im)
      Real w(m,m),v(m,m)
      real dta(m),UU(im,im),VV(im,im)
      integer sca(2,2),scb(2,2),scf(2,2),on,off
c
c Initialize arrays
c
      jmax=imax
      on=1
      off=0
      iiii=20998
c
c     fill initial matrix values
c
      do j=1,jmax
      do i=1,imax
         II(i,j) = 0.
         H(i,j) = 0.
         ZZ(i,j) = 0.
         PT(i,j) = 0. 
      enddo !i
      enddo !j
      do i=1,2
      do j=1,2
         sca(i,j) = 0
         scf(i,j) = 0 
         scb(i,j) = 0
      enddo !j
      enddo !i
c
c writeout parameter settings
c
      Do i=1,imax
         Z(i) = 0.
         II(I,I) = 1.
         H(i,i) = 1.
      enddo !i
c
c     Set obs 
c first guess - initial
c
      Do i=1,imax
         x(i) = dta(i)  
      enddo !i
c
c XT=FX
c
      Call mvmult(F,X,XT,imax,imax,1,m)
      call trans(F,FT,imax,imax,m)
c     call writev(F,imax,imax,m,'   F        ',atime,off,0.)
c     call writev(XT,imax,1,m,'   XT       ',atime,off,0.)
c
c PT=FPFT+T
c
      call mvmult(F,P,A,imax,imax,imax,m)
      call mvmult(A,FT,PT,imax,imax,imax,m)
      call addmv(PT,W,PT,1.,imax,imax,m)
c     call writev(PT,imax,imax,m,'   PT       ',atime,off,0.)
c
c K=PTH/(HPTHT+V)
c
      call mvmult(H,PT,A,imax,imax,imax,m)
      call trans(H,HT,imax,imax,m)
      Call mvmult(A,HT,E,imax,imax,imax,m)
      call addmv(E,V,ZZ,1.,imax,imax,m)
c     call writev(ZZ,imax,imax,m,'HPTHT+V 2INV',atime,off,0.)
      call matrixanal(ZZ,imax,imax,m, ' A ')
      call trans(ZZ,A,imax,imax,m)
      call replace(A,UU,imax,imax,m,m)
      call svdcmp(UU,imax,imax,m,m,B,VV)
c     call writev(B ,imax,1,m,'DIAG WJ     ',atime,on,0.)
c     call writev(UU,imax,imax,m,'  UU svdcmp ',atime,off,0.)
c     call writev(VV,imax,imax,m,'  VV svdcmp ',atime,off,0.)
      wmax=0.
      do j=1,imax
         if(b(j) .gt. wmax) wmax = b(j)
         g(j) = b(j) 
      enddo !j
      wmin=wmax*1.e-6
      do j=1,imax
         if(b(j) .lt. wmin) b(j) = 0.  
      enddo !j
      call mvmult(PT,H,A,imax,imax,imax,m)
      call trans(A,D ,imax,imax,m)
      do j=1,imax
         do i=1,imax
            c(i) = d(i,j)
         enddo !i 
         call svbksb(UU,b,VV,imax,imax,m,m,c,z)
         do i=1,imax
            zz(i,j) = z(i)
         enddo !i
      enddo!on j
      call trans(ZZ,K,imax,imax,m)
c
c     call invert(A,imax,m,E)  
c     
c this is single value decomposition solution for A
c
      call trans(UU,zz,imax,imax,m)
      call zero(E,m,m)     
      call mvmult(E,ZZ,UU,imax,imax,imax,m)
      call mvmult(VV,UU,ZZ,imax,imax,imax,m) 
      call trans(ZZ,A ,imax,imax,m)
c     call writev(D,imax,imax,m,'PT TRANS    ',atime,off,0.)
c     call writev(A,imax,imax,m,'AT INVERTED ',atime,off,0.)
      call writev(K,imax,imax,m,'KALMAN GAIN ',atime,off,0.)
c
c     est obs loop
c
c X=XT+K(Y-HXT)      
c
      call mvmult(H,XT,C,imax,imax,1,m)
      call addmv(Y,C,B,-1.,imax,1,m)
      call mvmult(K,B,C,imax,imax,1,m)
      call addmv(XT,C,X,1.,imax,1,m)
      call mvmult(k,H,E,imax,imax,imax,m)
      call addmv(II,E,A,-1.,imax,imax,m)
      call mvmult(A,PT,P,imax,imax,imax,m)
      do i=1,imax
         write(*,*) 'i,x,xt,y,k,w,v ',i,x(i),xt(i),
     &              y(i),k(i,i),w(i,i),v(i,i)
      enddo !i
      call writev(P,imax,imax,m,'ANAL COV ERR',atime,off,0.)
c
      return
      end
