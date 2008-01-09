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

SUBROUTINE Intplt2(p,x,vin,vout)

!==========================================================
!  This routine interpolates a value using a second order
!  scheme.
!
!  Input:
!	p: position where to interpolate;
!	x: positions of 2 grid points;
!	vin: input grid values;
!  Output:
!	vout: interpolated value;
!
!  Reference: intpl2.doc on MAC notebook.
!
!  HISTORY:
!	JUL. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  REAL, INTENT(IN) :: p,x(2),vin(2)
  REAL, INTENT(OUT) :: vout

  INTEGER :: i,j
  REAL :: e(2),r,b,c

  ! Check the 4 grid points:
  DO i=1,2
    DO j=1,2
      IF ((i .NE. j) .AND. &
        (ABS(x(i)-x(j)) .LT. 1.0e-5*x(i) ) ) THEN
        PRINT*,'Intplt2: There are overlap gridpoints, STOP'
	STOP
      ENDIF
    ENDDO
  ENDDO

  DO i=1,2
    IF (ABS(p-x(i)) .LT. 1.0e-5*MAXVAL(x) ) THEN
      e = 0.0
      e(i) = 1.0
      GOTO 1
    ENDIF
  ENDDO

  ! e2:
  e(2) = (p-x(1))/(x(2)-x(1))

  ! e1:
  e(1) = 1.0-e(2)

1 CONTINUE

  ! Interpolate:
  vout = 0.0
  DO i=1,2
    vout = vout+e(i)*vin(i)
  ENDDO

END SUBROUTINE Intplt2


SUBROUTINE Intplt4(p,x,vin,vout)

!==========================================================
!  This routine interpolates a value using a fourth order
!  scheme.
!
!  Input:
!	p: position where to interpolate;
!	x: positions of 4 grid points;
!	vin: input grid values;
!  Output:
!	vout: interpolated value;
!
!  Reference: intpl2.doc on MAC notebook.
!
!  HISTORY:
!	JUL. 2006 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  REAL, INTENT(IN) :: p,x(4),vin(4)
  REAL, INTENT(OUT) :: vout

  INTEGER :: i,j
  REAL :: e(4),r,b,c

  ! Check the 4 grid points:
  DO i=1,4
    DO j=1,4
      IF ((i .NE. j) .AND. &
        (ABS(x(i)-x(j)) .LT. 1.0e-5*x(i) ) ) THEN
        PRINT*,'Intplt4: There are overlap gridpoints, STOP'
	STOP
      ENDIF
    ENDDO
  ENDDO

  DO i=1,4
    IF (ABS(p-x(i)) .LT. 1.0e-5*MAXVAL(x) ) THEN
      e = 0.0
      e(i) = 1.0
      GOTO 1
    ENDIF
  ENDDO

  ! e4:
  r = 1.0
  c = 1.0
  DO i=1,3
    r = r*(p   -x(i))
    c = c*(x(4)-x(i))
  ENDDO
  e(4) = r/c

  ! e3:
  r = 1.0
  b = 1.0
  c = 1.0
  DO i=1,2
    r = r*(p   -x(i))
    c = c*(x(4)-x(i))
    b = b*(x(3)-x(i))
  ENDDO
  e(3) = (r-c*e(4))/b

  ! e2:
  e(2) = (p-x(1)-(x(4)-x(1))*e(4)-(x(3)-x(1))*e(3))/(x(2)-x(1))

  ! e1:
  e(1) = 1.0-e(2)-e(3)-e(4)

1 CONTINUE

  ! Interpolate:
  vout = 0.0
  DO i=1,4
    vout = vout+e(i)*vin(i)
  ENDDO

END SUBROUTINE Intplt4
