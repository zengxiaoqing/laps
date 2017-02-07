      SUBROUTINE SFCBKGD(bgmodel,t,q,height,tsfc,qsfc_i,
     &                   tdsfc_i,tdsfc_o,ter,p,IMX,JMX,KX,psfc,
     &                   nx_pr,ny_pr,istatus)

c inputs are from laps analyzed 2- 3d fields
C     INPUT       t        ANALYZED TEMPERATURE     3D
C                 q        Analyzed MIXING RATIO    3D
C                 height   ANALYZED HEIGHT          3D
C                 qsfc_i   Q Surface                2D
C                 ter      TERRAIN                  2D
C                 IMX      DIMENSION E-W
C                 JMX      DIMENSION N-S
C                 KX       NUMBER OF VERTICAL LEVELS
C                 p        LAPS Pressure levels     1D/3D (Pa)
c
c Compute surface variables on hi-res terrain using first 3D model level
c above the LAPS terrain. This routine is also called for reduced pressure
c where 'ter' is passed in as a constant array.

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
      INTEGER    nx_pr,ny_pr,ip,jp
      INTEGER    it, idebug, istatus, ishow_timer
      INTEGER    istat_qsfc_i, istat_tdsfc_i
      INTEGER    bgmodel

      LOGICAL    lfndz

      REAL       HEIGHT      ( IMX , JMX, KX )
      REAL       P           ( nx_pr,ny_pr,KX ) ! Pa
      REAL       PSFC        ( IMX , JMX )
      REAL       TSFC        ( IMX , JMX )
      REAL       Q           ( IMX , JMX , KX )
      REAL       rh3d        ( IMX , JMX , KX )
      REAL       rh2d        ( IMX , JMX )
      REAL       qsfc_i      ( imx,  jmx )
      REAL       tdsfc_i     ( imx,  jmx )
      REAL       tdsfc_o     ( imx,  jmx )
      REAL       T           ( IMX , JMX , KX )
      REAL       esat
      REAL       TER         ( IMX , JMX )
      REAL       make_td
      REAL       make_rh
      REAL       make_ssh
      REAL       t_ref,badflag
      REAL       tdsfc_o_min, tdsfc_o_max
      REAL       qsfc_i_min, qsfc_i_max
      REAL       p_mb
      REAL       r_missing_data
c
c if bgmodel = 6 or 8 then tdsfc_i is used         
c if bgmodel = 4      then qsfc_i is rh (WFO - RUC)
c if bgmodel = 3      then qsfc_i is rh and q is rh
c if bgmodel = 9      then no surface fields input. Compute all from 3d
c                     fields. q3d used. (NOS - ETA)
c otherwise qsfc_i is used directly with a test for g/kg or dimensionless
c 
      call get_r_missing_data(r_missing_data,istatus)

      write(6,*)
      write(6,*)' Subroutine sfcbkgd, bgmodel = ',bgmodel

      qsfc_i_min = minval(qsfc_i)
      qsfc_i_max = maxval(qsfc_i)
      write(6,*)' qsfc_i range = ',qsfc_i_min,qsfc_i_max
      if(qsfc_i_min .eq. r_missing_data .AND.
     1   qsfc_i_max .eq. r_missing_data       )then
         write(6,*)' NOTE: qsfc_i has missing data values'
         istat_qsfc_i = 0
      else
         istat_qsfc_i = 1
      endif

      write(6,*)' tdsfc_i range = ',minval(tdsfc_i),maxval(tdsfc_i)
      if(minval(tdsfc_i) .eq. r_missing_data .AND.
     &   maxval(tdsfc_i) .eq. r_missing_data       )then
         write(6,*)' NOTE: tdsfc_i has missing data values'
         istat_tdsfc_i = 0
      else
         istat_tdsfc_i = 1
      endif

      if(istat_qsfc_i .eq. 0 .and. istat_tdsfc_i .eq. 0)then
         write(6,*)
     1           ' WARNING: both qsfc_i and tdsfc_i have missing values'
      endif

      write(6,*)' ter range = ',minval(ter),maxval(ter)

      t_ref=-132.0
      if(bgmodel.eq.3.or.bgmodel.eq.9)then
         write(6,*)' bgmodel is ',bgmodel,' convert 3D q to rh'
         do k=1,kx
            do j=1,jmx
            do i=1,imx
               ip = min(i,nx_pr)
               jp = min(j,ny_pr)
               p_mb = p(ip,jp,k) / 100.
               rh3d(i,j,k)=make_rh(p_mb,t(i,j,k)-273.15,q(i,j,k)*1000.
     +,t_ref)*100.
            enddo
            enddo
         enddo

         badflag=0.
         write(6,*)
     1  ' Interp 3D T and RH to hi-res terrain, set tdsfc_o array to RH' 
         write(6,*)' height bottom level range = '
     1             ,minval(height(:,:,1)),maxval(height(:,:,1))
         call interp_to_sfc(ter,rh3d,height,imx,jmx,kx,
     &                      badflag,tdsfc_o) ! Here tdsfc_o temporarily is RH
         call interp_to_sfc(ter,t,height,imx,jmx,kx,badflag,
     &                      tsfc)
      endif

      istatus = ishow_timer()

      write(6,*)' Compute sfc fields using 3D model data'

      idebug = 1

      do j=1,jmx
      do i=1,imx
         k=0
         lfndz=.false.
         do while(.not.lfndz)
            k=k+1
            if(height(i,j,k).gt.ter(i,j))then
               lfndz=.true.

               if(k.gt.0)then

                  call compute_sfc_bgfields(bgmodel,imx,jmx,kx,i,j,k
     &                ,ter(i,j),height,t,p,q,t_ref,psfc(i,j),tsfc(i,j)
     &                ,qsfc_i(i,j),qsfc_i_min,qsfc_i_max,istat_qsfc_i
     &                ,tdsfc_i(i,j),istat_tdsfc_i
     &                ,tdsfc_o(i,j),idebug,ip,jp,nx_pr,ny_pr)      

                  idebug = 0

               endif

            endif
            if(k.eq.kx)lfndz=.true.
         enddo
      enddo
      enddo 

      istatus = 1
      
      tdsfc_o_min = minval(tdsfc_o)
      tdsfc_o_max = maxval(tdsfc_o)
      if(tdsfc_o_max .gt. 1000.)then
          write(6,*)' ERROR: tdsfc_o is out of bounds'
          istatus = 0
      endif

      write(6,*)' tdsfc_o range = ',tdsfc_o_min,tdsfc_o_max
      write(6,*)' psfc range = ',minval(psfc),maxval(psfc)

      if(minval(psfc) .lt. 50000.)then
          write(6,*)' ERROR: psfc is out of bounds'
          istatus = 0
      endif

      write(6,*)' returning from sfcbkgd...'
      write(6,*)

      return
      end
c
c----------------------------------------------------------------------------
c
      subroutine compute_sfc_bgfields(bgm,nx,ny,nz,i,j,k,ter,height
     &,t,p,q,t_ref,psfc_pa,tsfc,qsfc,qsfc_i_min,qsfc_i_max,istat_qsfc_i
     &,tdsfc_i,istat_tdsfc_i,tdsfc
     &,idebug,ip,jp,nx_pr,ny_pr)
c
c J. Smart 11/18/99 put existing code in subroutine for other process use in laps
c
      implicit none

      integer nx,ny,nz
      integer i,j,k            !I, i,j,k coordinate of point for calculation
      integer bgm              !I, model type {if = 0, then tdsfc input = qsfc}
      integer nx_pr,ny_pr,ip,jp
      integer istatus, idebug, init_td                          
      integer istat_qsfc_i,istat_tdsfc_i
      data init_td/0/
      save init_td

      real   p(nx_pr,ny_pr,nz) !I, pressure of levels (Pa)
      real   ter               !I, Terrain height at i,j
      real   t_ref             !I, reference temp for library moisture conv routines

      real   height(nx,ny,nz)  !I, heights of pressure levels 3d
      real   t(nx,ny,nz)       !I, temperatures 3d 
      real   q(nx,ny,nz)       !I, specific humidity 3d
      real   psfc_pa           !O, output surface pressure, pa
      real   tsfc              !I/O input sfc T, output recomputed T
      real   qsfc              !I   input sfc Q (g/kg)           
      real   tdsfc_i           !I   input sfc Td
      real   tdsfc             !O   hi-res Td (K)           

      real   qsfc_l            !L   surface spec hum, Input as q or computed internally (dimensionless)
      real   psfc_mb           !L 

      real   tbar,tsfc_c
      real   td1,td2,tdsfc_c
      real   p_mb,p_mb_p1,p_mb_m1                                 
      real   G,R
      real   ssh2,make_ssh,make_td
      real   dz,dzp,dtdz
      real   tvsfc,tvk,tbarv
      real   r_missing_data
      real   qsfc_i_min, qsfc_i_max

      parameter (G         = 9.8,
     &           R         = 287.04)
c
c if first guess values are missing data then return missing data
c
      call get_r_missing_data(r_missing_data,istatus)

      ip = min(i,nx_pr)
      jp = min(j,ny_pr)
      p_mb = p(ip,jp,k) / 100.

      if(tsfc.lt.500.0.and.t(i,j,k).lt.500.0) then 
c
c first guess psfc_mb without moisture consideration
c
          tbar=(tsfc+t(i,j,k))*0.5
          dz=height(i,j,k)-ter
          psfc_mb=p_mb*exp(G/(R*tbar)*dz)

!         Calculate qsfc_l according to model (bgm=0 for reduced P)
!         This is used below for virtual temperature and sfc P reduction
          if(istat_qsfc_i .eq. 1)then                  
           if(bgm.eq.0.or.bgm.eq.6.or.
     &        bgm.eq.8.or.bgm.eq.12)then                   ! qsfc is Td
             tsfc_c  = tsfc-273.15
             tdsfc_c = qsfc-273.15
             if(tdsfc_c .lt. -200.)then ! a la ssh2 error check
                 if(init_td .le. 100)write(6,*)
     1           ' WARNING: setting qsfc_l to zero, tdsfc_c = ',tdsfc_c
                 init_td = init_td + 1
!                tdsfc_c = tsfc_c - 30.                                 
                 qsfc_l = 0.0
             else
                 if(idebug .eq. 1)write(6,*)' qsfc is Td: calling ssh2'
                 qsfc_l=ssh2(psfc_mb,tsfc_c,tdsfc_c,t_ref)*.001 ! kg/kg
             endif

           elseif(bgm.eq.3.or.bgm.eq.4.or.bgm.eq.9)then    ! qsfc is RH
!    1                                .or.bgm.eq.13)then   ! qsfc is RH
             if(idebug .eq. 1)write(6,*)' qsfc is RH: calling make_ssh'
             qsfc_l=make_ssh(psfc_mb,tsfc-273.15
     &                      ,qsfc/100.,t_ref)*.001 !kg/kg

           else                                            ! qsfc is qsfc
             if(qsfc_i_max .gt. .050)then                 ! likely g/kg
                 qsfc_l=qsfc*.001                         ! make dimensionless
                 if(idebug .eq. 1)then
                     write(6,*)
     1                      ' compute_sfc_bgfields: qsfc_l = qsfc*.001 '       
     1                        ,qsfc_l,qsfc
                     write(6,*)' thus qsfc_i is assumed g/kg'
                 endif
             else
                 qsfc_l=qsfc                              ! keep dimensionless
                 if(idebug .eq. 1)then
                     write(6,*)' compute_sfc_bgfields: qsfc_l = qsfc '
     1                        ,qsfc_l,qsfc
                     write(6,*)' thus qsfc_i is assumed dimensionless'       
                 endif
             endif
           endif
    
           if(qsfc_l .lt. 0.0 .or. qsfc_l .gt. .050)then
             write(6,*)' WARNING in compute_sfc_bgfields, qsfc_l = '
     1                 ,qsfc_l      
             psfc_mb =r_missing_data
             tsfc =r_missing_data
             tdsfc=r_missing_data
           endif

          elseif(istat_tdsfc_i .eq. 1)then ! q is missing, so use td more directly
           tsfc_c  = tsfc-273.15
           tdsfc_c = tdsfc_i-273.15
           qsfc_l=ssh2(psfc_mb,tsfc_c,tdsfc_c,t_ref)*.001 ! kg/kg

          else ! both q and td are missing
           qsfc_l = 0.

          endif ! valid range for qsfc_i

c pressure
          tvsfc=tsfc*(1.+0.608*qsfc_l)
          tvk=t(i,j,k)*(1.+0.608*q(i,j,k))
          tbarv=(tvsfc+tvk)*.5
          psfc_mb=(p_mb*exp(G/(R*tbarv)*dz))

!         if(j .eq. ny/2)then
!             write(6,101)i,tvsfc,tvk,tbarv,height(i,j,k),dz,p_mb,psfc_pa
!    1                   ,(psfc_pa/100.-p_mb) / dz
!101          format(' i,tvsfc,tvk,tbarv,ht,dz,pk,psfc',i4,6f8.2,f9.2
!    1                                                 ,e15.6)
!         endif

          if(k.gt.1)then ! calculate tsfc,qsfc

             dzp=height(i,j,k)-height(i,j,k-1)
c temp
             dtdz=(t(i,j,k-1)-t(i,j,k))/dzp
             tsfc=t(i,j,k)+dtdz*dz

c dew point temp
             td2=make_td(p_mb,t(i,j,k)-273.15,q(i,j,k)*1000.,t_ref)
             p_mb_m1 = p(ip,jp,k-1) / 100.
             td1=make_td(p_mb_m1,t(i,j,k-1)-273.15,q(i,j,k-1)*1000.
     .,t_ref)
             tdsfc=td2+((td1-td2)/dzp)*dz+273.15

             if(idebug .eq. 1)then
                 write(6,*)' k/tdsfc (C) = ',k,tdsfc               
             endif
c
          else ! k=1: calculate tsfc,qsfc and psfc_mb

             dzp=height(i,j,k+1)-height(i,j,k)
             p_mb_p1 = p(ip,jp,k+1) / 100.
             td2=make_td(p_mb,t(i,j,k)-273.15
     &,q(i,j,k)*1000.,t_ref)+273.15
             td1=make_td(p_mb_p1,t(i,j,k+1)-273.15
     &,q(i,j,k+1)*1000.,t_ref)+273.15
             dtdz=(t(i,j,k)-t(i,j,k+1))/dzp
             tsfc=t(i,j,k)+dtdz*dz
             tbar=(tsfc+t(i,j,k))*0.5
             psfc_mb=p_mb*exp(G/(R*tbar)*dz)
             tdsfc=td2+((td1-td2)/dzp)*dz

             if(idebug .eq. 1)then
                 write(6,*)' k/tdsfc = ',k,tdsfc                                    
             endif

          endif

          if(psfc_mb .ne. r_missing_data)then
             psfc_pa = psfc_mb * 100.
          else
             psfc_pa = r_missing_data
          endif

       else
          psfc_pa =r_missing_data
          tsfc =r_missing_data
          tdsfc=r_missing_data
       endif

       return
 
       end
