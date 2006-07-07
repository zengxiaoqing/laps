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

SUBROUTINE StaggerLogP(ptop,psfc,znw,vin,nx,ny,nz,vout)

!==========================================================
!  This routine computes a grid function on a vertically
!  stagger grid using a given uniform grid function over a 
!  Log(p) grid.
!
!  Input:
!	ptop: 	Top pressure;
!	psfc:	Sfc pressure;
!	znw:	Eta value for a uniform grid;
!	vin:	Input grid function on a uniform grid;
!		(nx*ny*nz)
!	nx:	Number of uniform grid point in X;
!	ny:	Number of uniform grid point in Y;
!	nz:	Number of uniform vertical points;
!
!  Output:
!	vout:	Output grid funcion on a stagger grid.
!		(nx*ny*(nz-1))
!
!  HISTORY: APR. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL*4, INTENT(IN) ::  ptop,psfc(nx,ny),znw(nz)
  REAL, INTENT(IN) ::    vin(nx,ny,nz)
  REAL, INTENT(OUT) ::   vout(nx,ny,nz-1)

  ! Local variables:
  INTEGER :: i,j,k
  REAL :: a,b,p,p1,p2

  ! Linear interpolation:
  DO k=1,nz-1

    DO j=1,ny
      DO i=1,nx

	p1 = znw(k)*(psfc(i,j)-ptop)+ptop
	p2 = znw(k+1)*(psfc(i,j)-ptop)+ptop
	! Pressure value at the stagger grid:
	p = 0.5*(p1+p2)
    	a = LOG(p1/p)
    	b = LOG(p/p2)

    	vout(i,j,k) = (b*vin(i,j,k  )+ &
                       a*vin(i,j,k+1))/(a+b)

      ENDDO
    ENDDO

  ENDDO

END SUBROUTINE StaggerLogP


SUBROUTINE UnStaggerLogP(ptop,psfc,znu,vin,nx,ny,nz,vout)

!==========================================================
!  This routine computes a grid function on a uniform grid
!  using a given stagger grid function over a Log(p) grid.
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
!
!  Output:
!	vout:	Output grid funcion on a uniform grid.
!		(nx*ny*nz)
!
!  HISTORY: APR. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: nx,ny,nz
  REAL*4, INTENT(IN) ::  ptop,psfc(nx,ny),znu(nz-1)
  REAL, INTENT(IN) ::    vin(nx,ny,nz-1)
  REAL, INTENT(OUT) ::   vout(nx,ny,nz)

  ! Local variables:
  INTEGER :: i,j,k
  REAL :: a,b,p,p1,p2

  ! Linear interpolation:

  ! Bottom boundary:
  DO j=1,ny
    DO i=1,nx
      p1 = znu(1)*(psfc(i,j)-ptop)+ptop
      p2 = znu(2)*(psfc(i,j)-ptop)+ptop
      ! Pressure value at the stagger grid:
      p = p1+0.5*(p1-p2)
      a = LOG(p1/p)
      b = LOG(p/p2)

      vout(i,j,1) = (b*vin(i,j,1)+ &
                     a*vin(i,j,2))/(a+b)

    ENDDO
  ENDDO

  ! Interior:
  DO k=2,nz-1

    DO j=1,ny
      DO i=1,nx

	p1 = znu(k-1)*(psfc(i,j)-ptop)+ptop
	p2 = znu(k)*(psfc(i,j)-ptop)+ptop
	! Pressure value at the stagger grid:
	p = 0.5*(p1+p2)
    	a = LOG(p1/p)
    	b = LOG(p/p2)

    	vout(i,j,k) = (b*vin(i,j,k-1)+ &
                       a*vin(i,j,k  ))/(a+b)

      ENDDO
    ENDDO

  ENDDO

  ! Top boundary:
  DO j=1,ny
    DO i=1,nx
      p1 = znu(nz-1)*(psfc(i,j)-ptop)+ptop
      p2 = znu(nz-2)*(psfc(i,j)-ptop)+ptop
      ! Pressure value at the stagger grid:
      p = p1+0.5*(p1-p2)
      a = LOG(p1/p)
      b = LOG(p/p2)

      vout(i,j,nz) = (b*vin(i,j,nz-1)+ &
                      a*vin(i,j,nz-2))/(a+b)

    ENDDO
  ENDDO

END SUBROUTINE UnStaggerLogP
