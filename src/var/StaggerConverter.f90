!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis
!dis
!dis
!dis

SUBROUTINE StaggerZ(vin,nx,ny,nz,ptop,psfc,znw,iorder,logP,vout)

!==========================================================
!  This routine computes a grid function on a vertically
!  stagger grid using a given uniform grid function with
!  options of 2nd/4th order accuracy and LogP transfer.
!
!  Input:
!	ptop: 	Top pressure;
!	psfc:	Sfc pressure;
!	znw:	Grid spacing ratio relative to psfc (nz);
!	vin:	Input grid function on a uniform grid;
!		(nx*ny*nz)
!	nx:	Number of uniform grid point in X;
!	ny:	Number of uniform grid point in Y;
!	nz:	Number of uniform vertical points;
!	logP:	0 -- standard interpolation;
!		1 -- log(P) interpolation;
!       iorder: 2 or 4 for second or fourth order;
!
!  Output:
!	vout:	Output grid funcion on a stagger grid.
!		(nx*ny*(nz-1))
!
!  HISTORY: JUL. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,logP,iorder
  REAL, INTENT(IN) ::  ptop,psfc(nx,ny),znw(nz)
  REAL, INTENT(IN) ::    vin(nx,ny,nz)
  REAL, INTENT(OUT) ::   vout(nx,ny,nz-1)

  ! Local variables:
  INTEGER :: i,j,k,l,m
  REAL :: p,x(4)

  ! Lower boundary:
  DO j=1,ny
    DO i=1,nx

      DO l=1,iorder
        x(l) = znw(l)*(psfc(i,j)-ptop)+ptop
      ENDDO

      ! First stagger grid point:
      p = 0.5*(x(1)+x(2))

      ! Log(P):
      IF (logP .EQ. 1) THEN
        x = LOG(x)
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
        CALL Intplt2(p,x,vin(i,j,1:iorder),vout(i,j,1))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,1:iorder),vout(i,j,1))

    ENDDO
  ENDDO

  ! Interior:
  DO k=2,nz-2

    m = k-iorder/2

    DO j=1,ny
      DO i=1,nx

        DO l=1,iorder
          x(l) = znw(m+l)*(psfc(i,j)-ptop)+ptop
        ENDDO

        ! First stagger grid point:
        p = 0.5*(x(iorder/2)+x(iorder/2+1))

        ! Log(P):
        IF (logP .EQ. 1) THEN
          x = LOG(x)
          p = LOG(p)
        ENDIF

        ! Interpolation:
        IF (iorder .EQ. 2) &
          CALL Intplt2(p,x,vin(i,j,k:k+1),vout(i,j,k))
        IF (iorder .EQ. 4) &
          CALL Intplt4(p,x,vin(i,j,k-1:k+2),vout(i,j,k))

      ENDDO
    ENDDO

  ENDDO

  ! Upper boundary:
  m = nz-iorder
  DO j=1,ny
    DO i=1,nx

      DO l=1,iorder
        x(l) = znw(m+l)*(psfc(i,j)-ptop)+ptop
      ENDDO

      ! First stagger grid point:
      p = 0.5*(x(iorder-1)+x(iorder))

      ! Log(P):
      IF (logP .EQ. 1) THEN
        x = LOG(x)
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
        CALL Intplt2(p,x,vin(i,j,nz-1:nz),vout(i,j,nz-1))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,nz-3:nz),vout(i,j,nz-1))

    ENDDO
  ENDDO

END SUBROUTINE StaggerZ

SUBROUTINE UnStaggerZ(vin,nx,ny,nz,ptop,psfc,znu,iorder,logP,vout)

!==========================================================
!  This routine computes a grid function on a uniform grid
!  using a given stagger grid function with options of a
!  second/fourth order accuracy and Log P transfer.
!
!  Input:
!	ptop: 	Top pressure;
!	psfc:	Sfc pressure;
!	znw:	Eta value for a uniform grid;
!	vin:	Input grid function on a stagger grid;
!		(nx*ny*(nz-1))
!	nx:	Number of uniform grid point in X;
!	ny:	Number of uniform grid point in Y;
!	nz:	Number of uniform vertical points;
!	logP:	0 -- standard interpolation;
!		1 -- log(P) interpolation;
!	iorder: 2 or 4 order accuracy;
!
!  Output:
!	vout:	Output grid funcion on a uniform grid.
!		(nx*ny*nz)
!
!  HISTORY: JUL. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz,logP,iorder
  REAL, INTENT(IN) ::  ptop,psfc(nx,ny),znu(nz-1)
  REAL, INTENT(IN) ::    vin(nx,ny,nz-1)
  REAL, INTENT(OUT) ::   vout(nx,ny,nz)

  ! Local variables:
  INTEGER :: i,j,k,l,m
  REAL :: p,x(4)

  ! Lower boundary:
  DO j=1,ny
    DO i=1,nx

      DO l=1,iorder
        x(l) = znu(l)*(psfc(i,j)-ptop)+ptop
      ENDDO

      ! First stagger grid point:
      p = 1.5*x(1)-0.5*x(2)

      ! Log(P):
      IF (logP .EQ. 1) THEN
        x = LOG(x)
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
        CALL Intplt2(p,x,vin(i,j,1:2),vout(i,j,1))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,1:4),vout(i,j,1))

      ! Second stagger grid point:
      p = 0.5*(x(1)+x(2))

      ! Log(P):
      IF (logP .EQ. 1) THEN
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
        CALL Intplt2(p,x,vin(i,j,1:2),vout(i,j,2))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,1:4),vout(i,j,2))

    ENDDO
  ENDDO

  ! Interior:
  DO k=3,nz-2

    m = k-iorder/2-1

    DO j=1,ny
      DO i=1,nx

        DO l=1,iorder
          x(l) = znu(m+l)*(psfc(i,j)-ptop)+ptop
        ENDDO

        ! First stagger grid point:
        p = 0.5*(x(iorder/2)+x(iorder/2+1))

        ! Log(P):
        IF (logP .EQ. 1) THEN
          x = LOG(x)
          p = LOG(p)
        ENDIF

        ! Interpolation:
        IF (iorder .EQ. 2) &
          CALL Intplt2(p,x,vin(i,j,k-1:k),vout(i,j,k))
        IF (iorder .EQ. 4) &
          CALL Intplt4(p,x,vin(i,j,k-2:k+1),vout(i,j,k))

      ENDDO
    ENDDO

  ENDDO

  ! Upper boundary:
  m = nz-iorder-1
  DO j=1,ny
    DO i=1,nx

      DO l=1,iorder
        x(l) = znu(m+l)*(psfc(i,j)-ptop)+ptop
      ENDDO

      ! First stagger grid point from the top:
      p = 1.5*x(iorder)-0.5*x(iorder-1)

      ! Log(P):
      IF (logP .EQ. 1) THEN
        x = LOG(x)
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
	CALL Intplt2(p,x,vin(i,j,nz-2:nz-1),vout(i,j,nz))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,nz-4:nz-1),vout(i,j,nz))

      ! Second stagger grid point:
      p = 0.5*(x(iorder-1)+x(iorder))

      ! Log(P):
      IF (logP .EQ. 1) THEN
        p = LOG(p)
      ENDIF

      ! Interpolation:
      IF (iorder .EQ. 2) &
	CALL Intplt2(p,x,vin(i,j,nz-2:nz-1),vout(i,j,nz-1))
      IF (iorder .EQ. 4) &
        CALL Intplt4(p,x,vin(i,j,nz-4:nz-1),vout(i,j,nz-1))

    ENDDO
  ENDDO

END SUBROUTINE UnStaggerZ
