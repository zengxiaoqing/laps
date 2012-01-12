       subroutine multiwind_noZ(u,v,rms,u0,v0,x,y,height
     1                         ,n,xx,yy,ht,vr,rmsmax,ier)
c***********************************************************************
c description      : To derive horizontal wind from multi-radar radial
c                    velocity at constant height without reflectivity data.
c                                                                       
c I/O parameters   :                                                    
c  I/O/W   name,            type,      description
c    O     u,v              real       the horizontal wind component of x,y coordinate
c                                      in unit of m/s. 
c    O     rms              real       root mean square error relative to initial guess wind 
c                                      in unit of m/s.
c                                      rms=sqrt( (u-u0)**2 + (v-v0)**2 )
c    I     u0,v0            real     initial guess (or model wind) horizontal wind
c                                      in unit of m/s.
c    I     x,y              real     the (x,y) coordinate at the wind position.
c                                      in unit of meter.
c    I     height           real     the height at the wind position 
c                                      in unit of meter.
c    I     n                integer    number of input radar data. (maximum number, n=4 )
c    I     xx(n)          real array the x-coordinate of radar center in unit of meter.
c    I     yy(n)          real array the y-coordinate of radar center in unit of meter.
c    I     ht(n)          real array the height of radar antena in unit of meter.
c    I     vr(n)          real array the radial velocity observed by radar in unit of m/s.
c    I     rmsmax           real       the maximum value of root mean square error relative
c                                      to initial guess wind in unit of m/s.
c                                      If the rms of multi-radar (n>2) wind is greater than 
c                                      rmsmax, the final computed horizontal wind is obtained 
c                                      by the single radar method.
c    O     ier              integer    =0, success message.
c                                      =1, n < 1, or n > 4.
c                                      =2, ht() > height or ht() > 30000 meters.
c                                      =3, sqrt( (x-xx())**2+(y-yy())**2 ) < 1. meter.
c                                      =4, single Doppler method used due to ambiguities
c
c Date :
c   May. 14, 2004 (S.-M. Deng)
c***********************************************************************

      real x,y,height
      real xx(n),yy(n),ht(n),vr(n)
      real u,v,rms,u0,v0,rmsmax
      integer n,ier,i,m

      real xa(4),yb(4),d(4),vh(4)
      real uu(6),vv(6),srms(6),vcos(6),wu(6),wv(6)
      real rcos,su,sv
      real vt

c-----------------------------------------------------------------------
c*  To check input data.

      ier=0
      if( (n.lt.1).or.(n.gt.4) )then
          ier=1
          return
      endif
      rmsmax=10.

      do i=1,n
         if( (ht(i).gt.height).or.(ht(i).gt.30000.) )then
             ier=2
             return
         endif
      enddo
      do i=1,n
         xa(i)=x-xx(i)
         yb(i)=y-yy(i)
         d(i)=sqrt( xa(i)**2+yb(i)**2 )
         if( d(i).lt.1.0 )then
             ier=3
             return
         endif
      enddo 

c-----------------------------------------------------------------------
c*  To compute the horizontal radial velocity.

      do i=1,n
         rcos=d(i)/sqrt( d(i)**2+(ht(i)-height)**2 )
         vh(i)=vr(i)*rcos
      enddo 

c-----------------------------------------------------------------------
c*  For single radar (n=1), to reserve tangent component of initial wind,
c                and combine observed horizontal radial velocity. 

      if( n.eq.1 )then
          vt=( u0*yb(1)-v0*xa(1) )/d(1)
          u=( vh(1)*xa(1)+vt*yb(1) )/d(1)
          v=( vh(1)*yb(1)-vt*xa(1) )/d(1)
          rms=sqrt( (u-u0)**2 + (v-v0)**2 )
          return
      endif

c-----------------------------------------------------------------------
c*  For dual-radar (n=2), to compute horizontal wind by Dual-Doppler method.

      if( n.eq.2 )then
          call duwind(u,v,rms,vcos(1),u0,v0,xa(1),xa(2)
     1               ,yb(1),yb(2),d(1),d(2),vh(1),vh(2))
          if( rms.gt.rmsmax )go to 100 
          return
      endif 

c-----------------------------------------------------------------------
c*  For multi-radar (n=3), to compute horizontal wind.

      if( n.eq.3 )then
          call duwind(wu(1),wv(1),srms(1),vcos(1),u0,v0,xa(1),xa(2)
     1               ,yb(1),yb(2),d(1),d(2),vh(1),vh(2))
          call duwind(wu(2),wv(2),srms(2),vcos(2),u0,v0,xa(1),xa(3)
     1               ,yb(1),yb(3),d(1),d(3),vh(1),vh(3))
          call duwind(wu(3),wv(3),srms(3),vcos(3),u0,v0,xa(2),xa(3)
     1               ,yb(2),yb(3),d(2),d(3),vh(2),vh(3))
          m=0
          do i=1,3
             if( (srms(i).lt.rmsmax).and.(vcos(i).lt.0.95) )then
                  m=m+1
                  uu(m)=wu(i)
                  vv(m)=wv(i)
             endif
          enddo
          if( m.eq.0 )go to 100
          su=0.
          sv=0.
          do i=1,m
             su=su+uu(i)
             sv=sv+vv(i)
          enddo
          u=su/float(m)
          v=sv/float(m)
          rms=sqrt( (u-u0)**2 + (v-v0)**2 )
          return
      endif

c-----------------------------------------------------------------------
c*  For multi-radar (n=4), to compute horizontal wind.

      call duwind(wu(1),wv(1),srms(1),vcos(1),u0,v0,xa(1),xa(2)
     1           ,yb(1),yb(2),d(1),d(2),vh(1),vh(2))
      call duwind(wu(2),wv(2),srms(2),vcos(2),u0,v0,xa(1),xa(3)
     1           ,yb(1),yb(3),d(1),d(3),vh(1),vh(3))
      call duwind(wu(3),wv(3),srms(3),vcos(3),u0,v0,xa(1),xa(4)
     1           ,yb(1),yb(4),d(1),d(4),vh(1),vh(4))
      call duwind(wu(4),wv(4),srms(4),vcos(4),u0,v0,xa(2),xa(3)
     1           ,yb(2),yb(3),d(2),d(3),vh(2),vh(3))
      call duwind(wu(5),wv(5),srms(5),vcos(5),u0,v0,xa(2),xa(4)
     1           ,yb(2),yb(4),d(2),d(4),vh(2),vh(4))
      call duwind(wu(6),wv(6),srms(6),vcos(6),u0,v0,xa(3),xa(4)
     1           ,yb(3),yb(4),d(3),d(4),vh(3),vh(4))
      m=0
      do i=1,6
         if( (srms(i).lt.rmsmax).and.(vcos(i).lt.0.95) )then
              m=m+1
              uu(m)=wu(i)
              vv(m)=wv(i)
         endif
      enddo
      if( m.eq.0 )go to 100
      su=0.
      sv=0.
      do i=1,m
         su=su+uu(i)
         sv=sv+vv(i)
      enddo
      u=su/float(m)
      v=sv/float(m)
      rms=sqrt( (u-u0)**2 + (v-v0)**2 )
      return

 100  continue ! single Doppler method with best fit radar
      rms=rmsmax
      do i=1,n
         vt=( u0*yb(i)-v0*xa(i) )/d(i)
         uu(i)=( vh(i)*xa(i)+vt*yb(i) )/d(i)
         vv(i)=( vh(i)*yb(i)-vt*xa(i) )/d(i)
         srms(i)=sqrt( (uu(i)-u0)**2 + (vv(i)-v0)**2 )
         if( rms.gt.srms(i) )then
             rms=srms(i)
             u=uu(i)
             v=vv(i) 
         endif
      enddo

      ier = 4

      return
      end

      subroutine duwind(u,v,rms,vcos,u0,v0,x1,x2,y1,y2,d1,d2,v1,v2)
      real u,v,rms,vcos,u0,v0,x1,x2,y1,y2,d1,d2,v1,v2
      real det
      det=x1*y2-x2*y1
      if( abs(det).lt.2. )then
          u=-999.0
          v=-999.0
          rms=999.0
          vcos=1.0
          return
      endif
      u=( v1*d1*y2-v2*d2*y1 )/det
      v=( v2*d2*x1-v1*d1*x2 )/det
      rms=sqrt( (u-u0)**2 + (v-v0)**2 )
      vcos=abs( x1*x2+y1*y2 )/(d1*d2)
      return
      end
