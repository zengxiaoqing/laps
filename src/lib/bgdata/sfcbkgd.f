      SUBROUTINE SFCBKGD(bgmodel,t,q,height,tsfc,tdsfc,ter,p,
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
c    "        11-18-99: Put psfc,tsfc, and tdsfc comps into subroutine
c
      IMPLICIT NONE

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
      REAL       rh3d        ( IMX , JMX , KX )
      REAL       rh2d        ( IMX , JMX )
      REAL       tdsfc       ( imx,  jmx )
      REAL       T           ( IMX , JMX , KX )
      REAL       esat
      REAL       TER         ( IMX , JMX )
      REAL       make_td
      REAL       make_rh
      REAL       make_ssh
      REAL       t_ref,badflag
c
c if bgmodel = 6 or 8 then tdsfc is indeed td sfc
c if bgmodel = 4      then tdsfc is rh (WFO - RUC)
c if bgmodel = 9      then no surface fields input. Compute all from 3d
c                     fields. q3d used. (NOS - ETA)
c otherwise tdsfc is input with q
c 
      t_ref=-47.0
      if(bgmodel.eq.3.or.bgmodel.eq.9)then
         do k=1,kx
            do j=1,jmx
            do i=1,imx
               rh3d(i,j,k)=make_rh(p(k),t(i,j,k)-273.15,q(i,j,k)*1000.
     +,t_ref)*100.
            enddo
            enddo
         enddo

         badflag=0.
         call interp_to_sfc(ter,rh3d,height,imx,jmx,kx,
     &                      badflag,tdsfc)
         call interp_to_sfc(ter,t,height,imx,jmx,kx,badflag,
     &                      tsfc)
      endif

      do j=1,jmx
      do i=1,imx

         k=0
         lfndz=.false.
         do while(.not.lfndz)
            k=k+1
            if(height(i,j,k).gt.ter(i,j))then
               lfndz=.true.
c
c
c recompute psfc with moisture consideration. snook's orig code.
c              it=tdsfc(i,j)*100.
c              it=min(45000,max(15000,it))
c              xe=esat(it)
c              qsfc=0.622*xe/(p_sfc-xe)
c              qsfc=qsfc/(1.+qsfc)
c --------------------------------------------------

               if(k.gt.1)then

                  call compute_sfc_bgfields(bgmodel,imx,jmx,kx,i,j,k
     &,ter(i,j),height,t,p,q,t_ref,psfc(i,j),tsfc(i,j),tdsfc(i,j))

               else
                  tsfc(i,j)=t(i,j,k)
                  tdsfc(i,j)=make_td(p(k),t(i,j,k)-273.15,q(i,j,k)*1000.
     &,t_ref)+273.15
               endif

            endif
            if(k.eq.kx)lfndz=.true.
         enddo
      enddo
      enddo 
      
      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine compute_sfc_bgfields(bgm,nx,ny,nz,i,j,k,ter,height
     &,t,p,q,t_ref,psfc,tsfc,tdsfc)
c
c J. Smart 11/18 put existing code in subroutine for other process use in laps
c
      implicit none

      integer nx,ny,nz
      integer i,j,k            !I, i,j,k coordinate of point for calculation
      integer bgm              !I, model type {if = 0, then tdsfc input = qsfc}

      real*4 p(nz)             !I, pressure of levels
      real*4 ter               !I, Terrain height at i,j
      real*4 t_ref             !I, reference temp for library moisture conv routines

      real*4 height(nx,ny,nz)  !I, heights of pressure levels 3d
      real*4 t(nx,ny,nz)       !I, temperatures 3d 
      real*4 q(nx,ny,nz)       !I, specific humidity 3d
      real*4 psfc              !O, output surface pressure, pa
      real*4 tsfc              !I/O input model T, output recomputed T
      real*4 tdsfc             !I/O input model Td, output  recomputed Td

      real*4 qsfc              ! surface spec hum, Input as q or computed internally

      real*4 tbar
      real*4 td1,td2
      real*4 G,R
      real*4 ssh2,make_ssh,make_td
      real*4 dz,dzp,dtdz
      real*4 tvsfc,tvk,tbarv

      parameter (G         = 9.8,
     &           R         = 287.04)
c
c first guess psfc without moisture consideration
c
       tbar=(tsfc+t(i,j,k))*0.5
       dz=height(i,j,k)-ter
       psfc=p(k)*exp(G/(R*tbar)*dz)

       if(bgm.eq.6.or.bgm.eq.8)then      !Td
          qsfc=ssh2(psfc,tsfc-273.15
     &                  ,tdsfc-273.15,t_ref)*.001  !kg/kg

       elseif(bgm.eq.3.or.bgm.eq.4.or.bgm.eq.9)then  !RH
          qsfc=make_ssh(psfc,tsfc-273.15
     &                      ,tdsfc/100.,t_ref)*.001  !kg/kg

       else                              !q
          qsfc=tdsfc
       endif

c pressure
       tvsfc=tsfc*(1.+0.608*qsfc)
       tvk=t(i,j,k)*(1.+0.608*q(i,j,k))
       tbarv=(tvsfc+tvk)*.5
       psfc=(p(k)*exp(G/(R*tbarv)*dz))*100.  !return units = pa


       dzp=height(i,j,k)-height(i,j,k-1)
c temp
       dtdz=(t(i,j,k-1)-t(i,j,k))/dzp
       tsfc=t(i,j,k)+dtdz*dz

c dew point temp
       td2=make_td(p(k),t(i,j,k)-273.15,q(i,j,k)*1000.,t_ref)
       td1=make_td(p(k-1),t(i,j,k-1)-273.15,q(i,j,k-1)*1000.,t_ref)
       tdsfc=td2+((td1-td2)/dzp)*dz+273.15
c
       return
       end
