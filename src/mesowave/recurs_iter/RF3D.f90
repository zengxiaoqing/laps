SUBROUTINE RF3D(u,l,n,a,np)

!****************************************************
!  This routine applies the one-D recursive filter in
!  x, y, and z direction np times with alpha value of
!  a(1),a(2),a(3) to a three-D array u.
!
!  HISTORY: APR. 2003 by YUANFU XIE.
!****************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: l(3),n(3),np(3)
  REAL,    INTENT(IN) :: a(3)
  REAL                :: u(l(1),l(2),l(3))

  ! Local variables:
  INTEGER :: j,k

  ! X direction:
  DO k=1,n(3)
     DO j=1,n(2)
        CALL RF1D(u(1:n(1),j,k),n(1),a(1),np(1))
     ENDDO
  ENDDO

  ! Y direction:
  DO k=1,n(3)
     DO j=1,n(1)
        CALL RF1D(u(j,1:n(2),k),n(2),a(2),np(2))
     ENDDO
  ENDDO

  ! Z direction:
  DO k=1,n(2)
     DO j=1,n(1)
        CALL RF1D(u(j,k,1:n(3)),n(3),a(3),np(3))
     ENDDO
  ENDDO

END SUBROUTINE RF3D
