SUBROUTINE Functn(f,v,l,n,id,np,al)

!***************************************************************
!  This routine evaluates the cost function of a surface data
!  analysis problem.
!
!  HISTORY: Jan. 2004 by YUANFU XIE.
!***************************************************************

  IMPLICIT NONE

  DOUBLE PRECISION, INTENT(OUT) :: f
  
  INTEGER, INTENT(IN) :: l(4),n(4)           ! Grid dimensions
  INTEGER, INTENT(IN) :: id		     ! Variable id
  INTEGER, INTENT(IN) :: np(3,l(4))          ! RF passes

  REAL,    INTENT(IN) :: al(3,l(4))          ! RF alpha values
  REAL,    INTENT(IN) :: v(l(1),l(2),l(3)) ! Control grid

  ! Local variables:
  INTEGER :: iobs,i,j,k,ll
  REAL    :: x(l(1),l(2),l(3)),vo,a(2,3)

  INCLUDE 'OBSCommon.f90'

  !------------------------------------------------------------
  !  Recursive Filter:
  !------------------------------------------------------------
  x = v
  CALL RF3D(x(1,1,1),l,n,al(1,id),np(1,id))

  !------------------------------------------------------------
  ! Evaluate the cost function:
  !------------------------------------------------------------

  f = 0.0d0
  ll = 0
  DO iobs=1,nobs

     IF (id .EQ. vid(iobs)) THEN

        ! H(x) = y
        vo = 0.0
        a(1,1:3) = 1.0-coe(1:3,iobs)
        a(2,1:3) = coe(1:3,iobs)
        DO k=1,2
           DO j=1,2
              DO i=1,2
                 vo = vo + x(idx(1,iobs)+i-1,idx(2,iobs)+j-1, &
                             idx(3,iobs)+k-1)* &
                             a(i,1)*a(j,2)*a(k,3)
              ENDDO
           ENDDO
        ENDDO

        ! y - y0
        f = f + (vo-o(1,iobs))**2

     ENDIF

  ENDDO

  f = 0.5*f

END SUBROUTINE Functn
