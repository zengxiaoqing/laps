! 1 "Functn.F"
! 1 "<built-in>"
! 1 "<command line>"
! 1 "Functn.F"
SUBROUTINE Functn(f,v,l,n,id,np,al)

!***************************************************************
! This routine evaluates the cost function of a surface data
! analysis problem.
!
! HISTORY: Jan. 2004 by YUANFU XIE.
!***************************************************************

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: f

  INTEGER, INTENT(IN) :: l(4),n(4) ! Grid dimensions
  INTEGER, INTENT(IN) :: id ! Variable id
  INTEGER, INTENT(IN) :: np(3,l(4)) ! RF passes

  REAL, INTENT(IN) :: al(3,l(4)) ! RF alpha values
  REAL, INTENT(IN) :: v(l(1),l(2),l(3)) ! Control grid

  ! Local variables:
  INTEGER :: iobs,i,j,k
  REAL :: x(l(1),l(2),l(3)),vo,a(2,3)

! 1 "../../OBSCommon.f90" 1
!***********************************************************
! This header file defines large arrays for observations to
! avoid the limitation of parameter passing.
!
! HISTORY: JAN. 2004 by YUANFU XIE.
!***********************************************************

! Maximum number of OBS:
INTEGER, PARAMETER :: mobs = 200000

! Observations:
INTEGER :: nobs
INTEGER :: vid(mobs),idx(3,mobs)
REAL :: o(4,mobs),coe(3,mobs),w(mobs)
COMMON /OBSBlock/nobs,vid,idx,o,coe,w
! 26 "Functn.F" 2

  !------------------------------------------------------------
  ! Recursive Filter:
  !------------------------------------------------------------
  x = v
  CALL RF3D(x(1,1,1),l,n,al(1,id),np(1,id))

  !------------------------------------------------------------
  ! Evaluate the cost function:
  !------------------------------------------------------------

  f = 0.0d0

  DO iobs=1,nobs

     IF (id .EQ. vid(iobs)) THEN

        ! H(x) = y
        vo = 0.0
        a(1,1:3) = 1.0-coe(1:3,iobs)
        a(2,1:3) = coe(1:3,iobs)
        DO k=1,2
           IF ((idx(3,iobs)+k-1 .GE. 1) .AND. &
               (idx(3,iobs)+k-1 .GE. n(3))) THEN
           DO j=1,2
              DO i=1,2
                 vo = vo + x(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                             idx(3,iobs)+k-1)* &
                             a(i,1)*a(j,2)*a(k,3)
              ENDDO
           ENDDO
           ENDIF
        ENDDO

        ! y - y0
        f = f + w(iobs)*(vo-o(1,iobs))**2

     ENDIF

  ENDDO

  f = 0.5*f

END SUBROUTINE Functn
