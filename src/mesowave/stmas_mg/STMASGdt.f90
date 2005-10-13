SUBROUTINE Gradnt3D(grdt,grid,ngrd,obsv,nobs,indx,coef,wght)

!==========================================================
!  This routine evaluates gradients of STMAS cost function.
!
!  HISTORY:
!	Creation: JUN. 2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3),nobs,indx(6,nobs)
  REAL, INTENT(IN) :: grid(ngrd(1),ngrd(2),ngrd(3))
  REAL, INTENT(IN) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  REAL, INTENT(OUT) :: grdt(ngrd(1),ngrd(2),ngrd(3))

  ! Local variables:
  INTEGER :: i,j,k,io
  REAL :: gvl,pnl
  REAL :: lpl(ngrd(1),ngrd(2),ngrd(3))
  REAL :: gll(ngrd(1),ngrd(2),ngrd(3))

  ! Initial:
  grdt = 0.0

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

    ! Gradient:
    DO k=3,6,3
      DO j=2,5,3
        DO i=1,4,3
          grdt(indx(i,io),indx(j,io),indx(k,io)) = &
            grdt(indx(i,io),indx(j,io),indx(k,io))- &
            (obsv(1,io)-gvl)*coef(i,io)*coef(j,io)*coef(k,io)
        ENDDO
      ENDDO
    ENDDO

  ENDDO

  ! Penalty:

  gll = 0.0
  ! X:
  lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
         grid(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) &
    -2.0*grid(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) &
        +grid(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3))
  gll(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) = &
  gll(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3)) + &
        lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))
  gll(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3)) = &
  gll(3:ngrd(1)  ,1:ngrd(2),1:ngrd(3)) + &
        lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))
  gll(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
  gll(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) - &
    2.0*lpl(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3))

  ! Y:
  lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
         grid(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) &
    -2.0*grid(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) &
        +grid(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3))
  gll(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) = &
  gll(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3)) + &
        lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))
  gll(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3)) = &
  gll(1:ngrd(1),3:ngrd(2)  ,1:ngrd(3)) + &
        lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))
  gll(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
  gll(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) - &
    2.0*lpl(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3))

  ! T:
  lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
         grid(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) &
    -2.0*grid(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) &
        +grid(1:ngrd(1),1:ngrd(2),3:ngrd(3)  )
  gll(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) = &
  gll(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2) + &
        lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)
  gll(1:ngrd(1),1:ngrd(2),3:ngrd(3)  ) = &
  gll(1:ngrd(1),1:ngrd(2),3:ngrd(3)  ) + &
        lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)
  gll(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
  gll(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) - &
    2.0*lpl(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1)

  ! Add gll to grdt:
  grdt(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = &
    grdt(1:ngrd(1),1:ngrd(2),1:ngrd(3)) + &
    2.0*penalt*gll(1:ngrd(1),1:ngrd(2),1:ngrd(3))

END SUBROUTINE Gradnt3D
