SUBROUTINE Grid2Obs(indx,coef,obsv,nobs,wght,ngrd,dxyt,domn)

!==========================================================
!  This routine finds the indices and coefficients for an
!  interpolation scheme from a grid function to observation
!  sites.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)
  REAL, INTENT(IN) :: dxyt(3),domn(2,3)

  INTEGER, INTENT(OUT) :: indx(6,nobs)
  INTEGER, INTENT(INOUT) :: nobs
  REAL, INTENT(OUT) :: coef(6,nobs)
  REAL, INTENT(INOUT) :: obsv(4,nobs),wght(nobs)

  ! Local variables:
  INTEGER :: i,ier
  INTEGER :: nib		! Number of obs in box

  ! Count obs in box:
  nib = 0
  ! Interpolation for each obs:
  DO i=1,nobs
    CALL Intplt3D(obsv(2:4,i),ngrd,dxyt,domn, &
		  indx(1,i),coef(1,i),ier)

    ! Check:
    IF (ier .EQ. 0) THEN
      nib = nib+1

      ! Save the obs and its weight:
      obsv(1:4,nib) = obsv(1:4,i)
      wght(nib) = wght(i)

      indx(1:6,nib) = indx(1:6,i)
      coef(1:6,nib) = coef(1:6,i)

    ENDIF
  ENDDO

  ! Count inbox obs:
  nobs = nib

END SUBROUTINE Grid2Obs

SUBROUTINE Intplt3D(pstn,ngrd,gspc,domn,indx,coef,ierr)

!==========================================================
!  This routine returns interpolation coefficients and
!  indices of a given location from a grid in 3D space.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)	! Numbers gridpoint
  REAL,    INTENT(IN) :: pstn(3)	! Interpolate point
  REAL,    INTENT(IN) :: gspc(3)	! Grid spacing
  REAL,    INTENT(IN) :: domn(2,3)	! Domain

  INTEGER, INTENT(OUT) :: indx(6)	! Indices
  INTEGER, INTENT(OUT) :: ierr		! Error flag
  REAL,    INTENT(OUT) :: coef(6)	! Coefficients

  ! Local variables:
  INTEGER :: i

  ierr = 0
  ! Check if it is in box?
  DO i=1,3
     IF (pstn(i) .LT. domn(1,i)) ierr = -i
     IF (pstn(i) .GT. domn(2,i)) ierr = i
  ENDDO

  ! Indices:
  indx(1:3) = (pstn-domn(1,1:3))/gspc

  ! Coefficients:
  coef(1:3) = (pstn-indx(1:3)*gspc-domn(1,1:3))/gspc

  indx(1:3) = indx(1:3)+1

  indx(4:6) = MIN(indx(1:3)+1,ngrd)
  coef(4:6) = coef(1:3)
  coef(1:3) = 1.0-coef(4:6)

END SUBROUTINE Intplt3D
