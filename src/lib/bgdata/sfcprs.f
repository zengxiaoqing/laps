      SUBROUTINE SFCPRS (bgmodel,t,q,height,tsfc,tdsfc,ter,p,
     &IMX,JMX,KX,psfc)

c inputs are from laps analyzed 2- 3d fields
C     INPUT       t        ANALYZED TEMPERATURE     3D
C                 q        Analyzed MIXXING RATIO   3D
C                 height   ANALYZED HEIGHT          3D
C                 ter      TERRAIN                  2D
C                 IMX      DIMENSION E-W
C                 JMX      DIMENSION N-S
C                 KX       NUMBER OF VERTICAL LEVELS
C                 p        LAPS Pressure levels     1D
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
      INTEGER    bgmodel

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
      REAL       r_log_p,t_sfc,p_sfc,td_sfc
      REAL       TER         ( IMX , JMX )
      REAL       rterm
      REAL       ssh2
      REAL       make_td
      REAL       make_ssh
      REAL       dzp,dtdz
c     REAL       dpp,dp,pfrac,zfrac
      REAL       plo,phi,pla,slope
      REAL       td1,td2
c     REAL       td_sfc_logp,td_sfc_lapse
      REAL       td_sfc_slope
c
c if bgmodel = 6 or 8 then tdsfc is indeed td sfc
c if bgmodel = 4      then tdsfc is rh
c otherwise tdsfc is input with q
c
      do j=1,jmx
      do i=1,imx

         k=0
         lfndz=.false.
         do while(.not.lfndz)
            k=k+1
            if(height(i,j,k).gt.ter(i,j))then
               lfndz=.true.
c
c first guess psfc without moisture consideration
c
               tbar=(tsfc(i,j)+t(i,j,k))*0.5
               dz=height(i,j,k)-ter(i,j)
               p_sfc=p(k)*exp(G/(R*tbar)*dz)

c
c recompute psfc with moisture consideration. snook's orig code.
c              it=tdsfc(i,j)*100.
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              qsfc=0.622*xe/(p_sfc-xe)
c              qsfc=qsfc/(1.+qsfc)
c --------------------------------------------------

               if(bgmodel.eq.6.or.bgmodel.eq.8)then
                  qsfc=ssh2(p_sfc,tsfc(i,j)-273.15
     &                          ,tdsfc(i,j)-273.15,-47.)*.001  !kg/kg
               elseif(bgmodel.eq.4)then
                  qsfc=make_ssh(p_sfc,tsfc(i,j)-273.15
     &                          ,tdsfc(i,j)/100.,-47.)*.001  !kg/kg
               else 
                  qsfc=tdsfc(i,j)
               endif

               tvsfc=tsfc(i,j)*(1.+0.608*qsfc)
               tvk=t(i,j,k)*(1.+0.608*q(i,j,k))
               tbarv=(tvsfc+tvk)*.5
               p_sfc=(p(k)*exp(G/(R*tbarv)*dz))

               if(k.gt.1)then
c
                  plo=alog(p(k-1))
                  phi=alog(p(k))
                  pla=alog(p_sfc)
                  slope=(plo-pla)/(plo-phi)
                  dzp=height(i,j,k)-height(i,j,k-1)
c
c currently using slope in Td calc as this gives an interpolated value
c between the other two. However, the calc using pfrac is not far off,
c and perhaps a bit cooler than the others.
c
c                 zfrac=dz/dzp
c                 dpp=alog(p(k))-alog(p(k-1))
c                 dp=alog(p(k))-alog(p_sfc)
c                 pfrac=dp/dpp

                  dtdz=(t(i,j,k-1)-t(i,j,k))/dzp
                  t_sfc=t(i,j,k)+dtdz*dz

c                 t_sfc=t(i,j,k)+(t(i,j,k)-t(i,j,k-1))*zfrac

                  td2=make_td(p(k),t(i,j,k)-273.15,q(i,j,k)*1000.
     &,-47.)
                  td1=make_td(p(k-1),t(i,j,k-1)-273.15,q(i,j,k-1)*1000.
     &,-47.)
                  td_sfc_slope=td2-slope*(td1-td2)

c                 td_sfc_logp=td2+(td2-td1)*pfrac
c                 td_sfc_lapse=td2+((td1-td2)/dzp)*dz

                  td_sfc=td_sfc_slope
               endif
c 
c recompute Td sfc using recomputed p and Tv
c
c              tdsfc(i,j)=make_td(p_sfc,tvsfc-273.15,qsfc*1000.,-47.)
c    &                    +273.15
c
c recompute sfc T with new sfc p. "devirtualize" and return dry T.
c
c                 r_log_p = log(p_sfc/p(k))
c                 rterm   = (2.*G*dz)/(R*r_log_p)
c                 tvsfc   = rterm - tvk

c              t_sfc   = tvsfc/(1.+0.608*qsfc)

c                 r_log_p = log(p(k)/p_sfc)
c                 tbarv   = -(G*dz/R)/r_log_p
c                 tvsfc = tbarv*2.-tvk
c                 t_sfc   = tvsfc/(1.+0.608*qsfc)

               tsfc(i,j)=t_sfc               
               tdsfc(i,j)=td_sfc+273.15
               psfc(i,j)=p_sfc*100.      !return p units = pascals

            endif
            if(k.eq.kx)lfndz=.true.
         enddo
      enddo
      enddo 
      
      return
      end
