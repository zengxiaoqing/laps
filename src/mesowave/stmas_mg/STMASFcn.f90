SUBROUTINE Functn3D(fctn,grid,ngrd,obsv,nobs,indx,coef,wght)

!==========================================================
!  This routine evaluates the STMAS cost function.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3),nobs,indx(6,nobs)
  REAL, INTENT(IN) :: grid(ngrd(1),ngrd(2),ngrd(3))
  REAL, INTENT(IN) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  REAL, INTENT(OUT) :: fctn

  ! Local variables:
  INTEGER :: i,j,k,io
  REAL    :: gvl,pnl

  ! Initial:
  fctn = 0.0

  ! Jo term:
  DO io=1,nobs

    ! Grid value at obs position:
    gvl = 0.0
    DO k=3,6,3
      DO j=2,5,3
	DO i=1,4,3
	  gvl = gvl+grid(indx(i,io), &
			 indx(j,io), &
			 indx(k,io))* &
		coef(i,io)*coef(j,io)*coef(k,io)
	ENDDO
      ENDDO
    ENDDO

    ! Residual:
    fctn = fctn+(obsv(1,io)-gvl)**2

  ENDDO

  fctn = fctn*0.5

  ! Penalty term: Penalize the laplacian:

  ! X:
  pnl = 0.0
  DO k=1,ngrd(3)
    DO j=1,ngrd(2)
      DO i=2,ngrd(1)-1
        pnl = pnl + &
          (grid(i-1,j,k)-2.0*grid(i,j,k)+grid(i+1,j,k))**2
      ENDDO
    ENDDO
  ENDDO

  ! Y:
  DO k=1,ngrd(3)
    DO j=2,ngrd(2)-1
      DO i=1,ngrd(1)
        pnl = pnl + &
          (grid(i,j-1,k)-2.0*grid(i,j,k)+grid(i,j+1,k))**2
      ENDDO
    ENDDO
  ENDDO

  ! T:
  DO k=2,ngrd(3)-1
    DO j=1,ngrd(2)
      DO i=1,ngrd(1)
        pnl = pnl + &
          (grid(i,j,k-1)-2.0*grid(i,j,k)+grid(i,j,k+1))**2
      ENDDO
    ENDDO
  ENDDO

  fctn = fctn+penalt*pnl

END SUBROUTINE Functn3D
