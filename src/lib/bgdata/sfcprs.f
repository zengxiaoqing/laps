      SUBROUTINE SFCPRS (T,Q,HEIGHT,tsfc,tdsfc,TER,P,IMX,JMX,KX,
     &                   PSFC)

c inputs are from laps analyzed 2- 3d fields
C     INPUT       T        ANALYZED TEMPERATURE     3D
C                 Q        Analyzed MIXXING RATIO   3D
C                 HEIGHT   ANALYZED HEIGHT          3D
C                 TER      TERRAIN                  2D
C                 IMX      DIMENSION E-W
C                 JMX      DIMENSION N-S
C                 KX       NUMBER OF VERTICAL LEVELS
C                 P        LAPS Pressure levels     1D
C                 tdsfc    LAPS Dew Point Temp      2D

c
c J. Smart    09-22-98:	Original Version: This is used to compute
c                       sfc p for lgb when using NOGAPS1.0 deg since
c                       this field currently does not come with the
c                       model grids at AFWA.
c    "        02-01-99: Recompute tdsfc with new psfc and tsfc for
c                       consistency.
c
      IMPLICIT NONE

      REAL       G
      REAL       R
      parameter (G         = 9.8,
     &           R         = 287.04)

      INTEGER    I
      INTEGER    IMX
      INTEGER    J
      INTEGER    JMX
      INTEGER    K
      INTEGER    KX
      INTEGER    it

      LOGICAL    lfndz

      REAL       HEIGHT      ( IMX , JMX, KX )
      REAL       P           ( KX )
      REAL       PSFC        ( IMX , JMX )
      REAL       TSFC        ( IMX , JMX )
      REAL       Q           ( IMX , JMX , KX )
      REAL       tdsfc       ( imx,  jmx )
      REAL       T           ( IMX , JMX , KX )
      REAL       tbarv,dz,tvsfc,tvk,tbar
      REAL       qsfc,xe,psfc1
      REAL       esat
      REAL       r_log_p,t_sfc,p_sfc
      REAL       TER         ( IMX , JMX )
      REAL       rterm
      REAL       ssh2
      REAL       make_td

      do j=1,jmx
      do i=1,imx

         k=0
         lfndz=.false.
         do while(.not.lfndz)
            k=k+1
            if(height(i,j,k).gt.ter(i,j))then

               lfndz=.true.

c first guess psfc without moisture consideration
               tbar=(tsfc(i,j)+t(i,j,k))*0.5
               dz=height(i,j,k)-ter(i,j)
               p_sfc=p(k)*exp(G/(R*tbar)*dz)

c recompute psfc with moisture consideration
c              it=tdsfc(i,j)*100.
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              qsfc=0.622*xe/(p_sfc-xe)
c              qsfc=qsfc/(1.+qsfc)

               qsfc=ssh2(p_sfc,tsfc(i,j)-273.15,tdsfc(i,j)-273.15,0.)
     &              *.001  !kg/kg

               tvsfc=tsfc(i,j)*(1.+0.608*qsfc)
               tvk=t(i,j,k)*(1.+0.608*q(i,j,k))
               tbarv=(tvsfc+tvk)*.5

               p_sfc=(p(k)*exp(G/(R*tbarv)*dz))
c
c recompute Td sfc using recomputed p and Tv
               tdsfc(i,j)=make_td(p_sfc,tvsfc-273.15,qsfc*1000.,0.)
     &                    +273.15
c
c recompute sfc T with new sfc p. "devirtualize" and return dry T.
c

c                 r_log_p = log(p_sfc/p(k))
c                 rterm   = (2.*G*dz)/(R*r_log_p)
c                 tvsfc   = rterm - tvk

               t_sfc   = tvsfc/(1.+0.608*qsfc)

c                 r_log_p = log(p(k)/p_sfc)
c                 tbarv   = -(G*dz/R)/r_log_p
c                 tvsfc = tbarv*2.-tvk
c                 t_sfc   = tvsfc/(1.+0.608*qsfc)

               tsfc(i,j)=t_sfc               
               psfc(i,j)=p_sfc*100.      !return p units = pascals

            endif
            if(k.eq.kx)lfndz=.true.
         enddo
      enddo
      enddo 
      
      return
      end
