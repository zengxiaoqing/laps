SUBROUTINE Intplt3d(x,l,n,d,dm,idx,coe,ier)

!*************************************************
!  This routine computes an interpolation indices
!  and coefficients using bi-linear interpolation.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: l(3),n(3)
  REAL,    INTENT(IN) :: x(3),d(3),dm(2,3)

  INTEGER, INTENT(OUT) :: idx(3),ier
  REAL,    INTENT(OUT) :: coe(3)

  ! Local variables:
  INTEGER :: i

  ier = 0
  ! Check in box?
  DO i=1,3
     IF (x(i) .LT. dm(1,i)) ier = -i
     IF (x(i) .GT. dm(2,i)) ier = i
  ENDDO

  ! Indices:
  idx = (x-dm(1,1:3))/d

  ! Coefficients:
  coe = (x-idx*d-dm(1,1:3))/d

  idx = idx+1

END SUBROUTINE Intplt3d
