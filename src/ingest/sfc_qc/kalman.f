c
      subroutine kalman(F,dta,Y,p,it,w,v,x,xt,imax,m,atime,stn)
c
c*********************************************************************
c
c     Routine to apply Kalman Filter to a set of obs for QC purposes.
c     
c     Original: John McGinley, NOAA/FSL  Spring 1998
c     Changes:
c       21 Aug 1998  Peter Stamus, NOAA/FSL
c          Make code dynamic, housekeeping changes, for use in LAPS.
c       09 Dec 1999  John McGinley and Peter Stamus, NOAA/FSL
c          New version; additional housekeeping changes too.
c
c*********************************************************************
c
      Real K(m,m),P(m,m),F(m,m),II(m,m)
      Real H(m,m),HT(m,m),FT(m,m),G(m)
      Real X(m),Y(m),XT(m)
      Real PT(m,m),ZZ(m,m)
      Real A(m,m),B(m),C(m),Z(m),D(m,m),E(m,m)
      Real w(m,m),v(m,m)
      real dta(m),UU(m,m),VV(m,m)
      integer sca(2,2),scb(2,2),scf(2,2),on,off
      character atime*(*),stn(m)*5
c
c Initialize arrays
c
      on=1
      off=0
      iiii=20998
c
c     fill initial matrix values
c
      call zero(II, m,m)
      call zero( H, m,m)
      call zero(ZZ, m,m)
      call zero(PT, m,m)
      call zero(A, m,m)
      call zero(E, m,m)
      call zero(HT,m,m)
      
c
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
         II(i,i) = 1.
         H(i,i) = 1.
      enddo !i
c
c     Set obs 
c first guess - initial
c
      Do i=1,imax
         X(i) = dta(i)  
      enddo !i
c
c XT=FX
c
      Call mvmult(F,X,XT,imax,imax,1,m)
      call trans(F,FT,imax,imax,m)
c     call writev(F,imax,imax,m,'   F        ',atime,on ,1.0)
c     call writev(X,imax,1,m,'   X       ',atime,on,10000.)
c     call writev(XT,imax,1,m,'   XT       ',atime,on,10000.)
c
c PT=FPFT+T
c
      call mvmult(F,P,A,imax,imax,imax,m)
      call mvmult(A,FT,PT,imax,imax,imax,m)
      call addmv(PT,W,PT,imax,imax,m)
      call writev(PT,imax,imax,m,'   PT       ',atime,off,0.)
c
c K=PTH/(HPTHT+V)
c
      call mvmult(H,PT,A,imax,imax,imax,m)
      call trans(H,HT,imax,imax,m)
      Call mvmult(A,HT,E,imax,imax,imax,m)
      call addmv(E,V,ZZ,imax,imax,m)
      call writev(ZZ,imax,imax,m,'HPTHT+V 2INV',atime,off ,0.)
      idiag=0
      call matrixanal(ZZ,imax,imax,m,idiag, ' HPTHT+V  ')
      if(idiag.eq.1) then 
       call fastinv(ZZ,imax,imax,m)
       call mvmult(PT,H,A,imax,imax,imax,m)
       call mvmult(A,ZZ,K,imax,imax,imax,m)
       go to 34
      endif
      call trans(ZZ,A,imax,imax,m)
      call replace(A,UU,imax,imax,m,m)
      call svdcmp(UU,imax,imax,m,m,B,VV,m)
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
         call svbksb(UU,b,VV,imax,imax,m,m,c,z,m)
         do i=1,imax
            zz(i,j) = z(i)
         enddo !i
      enddo!on j
      call trans(ZZ,K,imax,imax,m)
c     call invert(A,imax,m,E,m)  
c     
c.....  This is single value decomposition solution for A
c
      call trans(UU,zz,imax,imax,m)
      call zero(E,m,m)     
      call mvmult(E,ZZ,UU,imax,imax,imax,m)
      call mvmult(VV,UU,ZZ,imax,imax,imax,m) 
      call trans(ZZ,A ,imax,imax,m)
c     call writev(D,imax,imax,m,'PT TRANS    ',atime,off,0.)
c     call writev(A,imax,imax,m,'AT INVERTED ',atime,off,0.)
 34   call writev(K,imax,imax,m,'KALMAN GAIN ',atime,off,0.)
c
c.....  Estimate obs loop
c
c X=XT+K(Y-HXT)      
c
      call mvmult(H,XT,C,imax,imax,1,m)
      call submv(Y,C,B,imax,1,m)
      call mvmult(K,B,C,imax,imax,1,m)
      call addmv(XT,C,X,imax,1,m)
c
c P=(I-K)PT
c
      call mvmult(K,H,E,imax,imax,imax,m)
      call submv(II,E,A,imax,imax,m)
      call mvmult(A,PT,P,imax,imax,imax,m)
c
      sum = 0.
      write(6,2000) 
 2000 format(1x,' Stn Indx',' Kalman X ',' Forecast ',' Observatn'
     &,' KalmGn','     W    ','     V    ')
      do i=1,imax
         write(6,1098) stn(i),i,X(i)-10000.,XT(i)-10000.,
     &              Y(i)-10000.,K(i,i),w(i,i),v(i,i)
 1098 format(1x,a5,i4,3f10.3,f7.4,2f10.3)
         sum=sum+K(i,i)
      enddo !i
      print*, 'MEAN KALMAN ',sum/float(imax)
      write(6,*) 'MEAN KALMAN ',sum/float(imax)
      call writev(P,imax,imax,m,'ANAL COV ERR',atime,off,0.)
c
      return
      end
c
c
      Subroutine kalmod(F,yta,byta,dta,ta,wmt,wot,wbt,offset,
     &                  imax,mwt,m)
c
c*********************************************************************
c
c     Kalman Filter tool.
c     
c     Original: John McGinley, NOAA/FSL  December 1999
c     Changes:
c
c       09 Dec 1999  Peter Stamus, NOAA/FSL
c          Housekeeping changes.
c
c*********************************************************************
c
      real yta(m),byta(m),dta(m),ta(m),wmt(m),wot(m),wbt(m)
      real mwt(m,m),F(m,m),a,b,c
c
      do i=1,imax
         sum=wmt(i)+wot(i)+wbt(i)
         a=0.5*(wmt(i)+wbt(i))/sum
         b=0.5*(wot(i)+wmt(i))/sum
         c=0.5*(wbt(i)+wot(i))/sum
         sum=0.
         sum1=0.
         if(mwt(i,i).eq.1.) then
             print*,'Station ',i,' is isolated: set buddy trend to 0'
             sum=0.
          else
            do j = 1,imax
             if(i.eq.j) go to 1
             sum=mwt(i,j)/(1.-mwt(i,i))*yta(j)+sum
             F(i,j)=0.
 1           continue
            enddo
         endif
         byta(i)=sum
         F(i,i)=a*(1.+yta(i)/(ta(i)+offset)) + 
     &          c*(1.+dta(i)/(ta(i)+offset)) +
     &          b*(1+byta(i)/(ta(i)+offset))
      enddo
c
      return
      end

