cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS
cdis
cdis    This software and its documentation are in the public domain and
cdis    are furnished "as is."  The United States government, its
cdis    instrumentalities, officers, employees, and agents make no
cdis    warranty, express or implied, as to the usefulness of the software
cdis    and documentation for any purpose.  They assume no responsibility
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making
cdis    the modifications.  If significant modifications or enhancements
cdis    are made to this software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis
cdis
cdis
cdis
cdis
cdis
cdis

      subroutine plot_windob(dir,spd,ri,rj,lat,lon,imax,jmax,relsize)

      include 'lapsparms.cmn'

      real*4 lat(imax,jmax),lon(imax,jmax)

      call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!     write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
 1234 format(1x,4i5,4e12.4,i4)

!     This keeps about the size of barbs relative to the domain
!     du=(imax)/200. * relsize   

!     This tries to keep the same size of barbs relative to the grid points
      du = relsize

      call get_border(imax,jmax,x_1,x_2,y_1,y_2)
      call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))

      rot = (standard_longitude - lon(nint(ri),nint(rj))) / 57.295

!     Convert ri and rj to x1 and y1 (U and V)
!     call supcon(alat,alon,x1,y1)
!     x1 = umin + (umax - umin) * (ri-1.) / float(imax-1)
!     y1 = vmin + (vmax - vmin) * (rj-1.) / float(jmax-1)

      if(dir .gt. -400.)then
          call barbs(spd,dir,ri,rj,du,rot,-1e10,+1e10,-1e10,+1e10)
      endif

    1 continue
      return
      end


      subroutine plot_vr(i,j,vel,imax,jmax,c1_plottype
     1                ,alat_radar,alon_radar)

      character*1 c1_plottype

!     call getset(mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype)
!     write(6,1234) mxa,mxb,mya,myb,umin,umax,vmin,vmax,ltype
!1234 format(1x,4i5,4e12.4,i4)
!     du=(umax-umin)/200.

!     call get_border(imax,jmax,x_1,x_2,y_1,y_2)
!     call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))

!     du = (umax-umin)/imax
!     dv = (vmax-vmin)/jmax

      du = 1. / float(imax-1)
      dv = 1. / float(jmax-1)

      u = float(i-1) / float(imax-1)
      v = float(j-1) / float(jmax-1)

      if(.true.)then

!     call get_border(imax,jmax,x_1,x_2,y_1,y_2)
!     call set(x_1,x_2,y_1,y_2,1.,float(imax),1.,float(jmax))
!     u = i
!     v = j
!     du = 1.
!     dv = 1.

      endif

      if(c1_plottype .eq. 'y')then

          icol = min(max(  110. + vel / 2.0   ,  100.  )  ,120.)
          call setusv_dum(2hIN,icol)
!         call setusv_dum(2hIN,40)

          do uu = u-du/2.,u+du/2.,du/25.
              call line(uu,v-dv/2.,uu,v+dv/2.)
          enddo

      elseif(c1_plottype .eq. 'a')then ! Plot arrows


      endif

      return
      end


      subroutine map(iproj,selat,selon,nwlat,nwlon)
      real nwlat,nwlon
      cenlat=90. ! (selat+nwlat)/2.
      cenlon=-105. ! (selon+nwlon)/2.
      rot=0.
      igrid=0
      idetail=4
      idot_pat=0
      call supmap(iproj,cenlat,cenlon,0.0,selat,selon,nwlat,nwlon,-2,
     1igrid,idetail,idot_pat,ier)
!     call savesup_gp('laps_sup.parms',Istatus)
      return
      end


c
        SUBROUTINE BARBS(SPD,DIR,U,V,DU,PROJROT,umin,umax,vmin,vmax)
C
C---Wind barb plotter from barbs_gp.
C
c       INCLUDE 'SYSDISK:[GUDOC]EDFVAXBOX.FOR'
C       Paul Schultz    30-JUN-1982     Original version
C       HOAGENSON       02-SEP-1982     Changed name from BARBS
C       Schultz         25-OCT-1984     Changed zero-wind symbol
C
C---DIR,SPD     WIND DIRECTION AND SPEED
C---U,V         NCAR COORDINATES OF OBS
C---DU          ARBITRARILY CHOSEN INCREMENT OF SCREEN SPACE
C---PROJROT     ROTATION OF BARB DUE TO MAP PROJECTION
C
        DATA DEGRAD/.01745329/
C
C---RETURN IF MISSING SPEED OR DIRECTION.
C
        IF (DIR .LT. 0. .OR. SPD .LT. 0.) RETURN
        IF (DIR .GT. 360. .OR. SPD .GT. 200.) RETURN
C
C---DIRECTIONS:
C
        DR=DIR*DEGRAD+PROJROT
        DR1=(DIR+60.)*DEGRAD+PROJROT
        SIND=SIN(DR)
        SIND1=SIN(DR1)
        COSD=COS(DR)
        COSD1=COS(DR1)
C
C---LENGTHS:
C
        STAFF=DU*4.                     !  MULTIPLIER ARBITRARILY CHOSEN
        BARB=STAFF*.5
        ADD=STAFF*.3
C
C---SPEED AND COUNTERS:
C
        N50=0
        N10=0
C
!       if(u .le. umin .or. u .ge. umax .or. v .le. vmin .or. v .ge. vmax)
!       1                                                       return

        IF (SPD .LT. 1.0) THEN
        CALL PWRIT (U,V,'O',1,6,0,0)
        RETURN
        END IF
        IF (SPD .LT. 2.5) THEN
        X1=U
        Y1=V
        X2=X1+SIND*STAFF
        Y2=Y1+COSD*STAFF

        if(x2 .lt. umin .or. x2 .gt. umax .or. y2 .lt. vmin .or. y2 .gt.
     1 vmax)
     1                                                  return

        CALL LINE(X1,Y1,X2,Y2)
        return
        endif
C
        SP=SPD+2.5
C
        DO WHILE (SP .GE. 50.)
        N50=N50+1
        SP=SP-50.
        END DO
C
        DO WHILE (SP .GE. 10.)
        N10=N10+1
        SP=SP-10.
        END DO
C
C---DRAW STAFF
C
        X1=U
        Y1=V
        X2=X1+SIND*STAFF
        Y2=Y1+COSD*STAFF

        if(x2 .lt. umin .or. x2 .gt. umax .or. y2 .lt. vmin .or. y2 .gt.
     1 vmax)
     1                                                  return

        CALL LINE(X1,Y1,X2,Y2)
C
C---PLOT HALF-BARB, IF NECESSARY
C
        IF (SP .GE. 5.) THEN
        X1=X2+SIND1*ADD
        Y1=Y2+COSD1*ADD

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        CALL LINE(X1,Y1,X2,Y2)
        IF (N50 .NE. 0 .OR. N10 .NE. 0) GO TO 40
        X1=X2+SIND*ADD
        Y1=Y2+COSD*ADD

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        CALL LINE(X1,Y1,X2,Y2)
        RETURN
        END IF
C
40      X1=X2
        Y1=Y2
C
C---PLOT BARBS, IF NECESSARY
C
        DO 50 I=1,N10
        X2=X1+SIND*ADD
        Y2=Y1+COSD*ADD
        X3=X2+SIND1*BARB
        Y3=Y2+COSD1*BARB

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        CALL FRSTPT(X1,Y1)

!       if(x2 .le. umin .or. x2 .ge. umax .or. y2 .le. vmin .or. y2 .ge. vmax)
!       1                                                       return

        CALL VECTOR(X2,Y2)

!       if(x3 .le. umin .or. x3 .ge. umax .or. y3 .le. vmin .or. y3 .ge. vmax)
!       1                                                       return

        CALL VECTOR(X3,Y3)
        X1=X2
        Y1=Y2
50      CONTINUE
C
C---PLOT FLAGS, IF NECESSARY
C
        DO 60 I=1,N50
        X2=X1+SIND*ADD
        Y2=Y1+COSD*ADD
        X3=X2+SIND1*BARB
        Y3=Y2+COSD1*BARB

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        CALL FRSTPT(X1,Y1)

!       if(x2 .le. umin .or. x2 .ge. umax .or. y2 .le. vmin .or. y2 .ge. vmax)
!       1                                                       return

        CALL VECTOR(X2,Y2)

!       if(x3 .le. umin .or. x3 .ge. umax .or. y3 .le. vmin .or. y3 .ge. vmax)
!       1                                                       return

        CALL VECTOR(X3,Y3)

!       if(x1 .le. umin .or. x1 .ge. umax .or. y1 .le. vmin .or. y1 .ge. vmax)
!       1                                                       return

        CALL VECTOR(X1,Y1)
        X1=X2
        Y1=Y2
60      CONTINUE
C
        RETURN
        END
c
