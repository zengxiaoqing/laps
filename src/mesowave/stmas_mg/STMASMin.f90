SUBROUTINE Minimize(sltn,mlds,ngrd,obsv,nobs,indx,coef,wght)

!==========================================================
!  This routine minimizes the STMAS cost function.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mlds(3),ngrd(3),nobs,indx(6,nobs)
  REAL, INTENT(IN) :: obsv(4,nobs),coef(6,nobs),wght(nobs)
  REAL, INTENT(INOUT) :: sltn(mlds(1),mlds(2),mlds(3))

  !** LBFGS_B variables:
  INTEGER, PARAMETER :: msave=7		! Max iter save

  CHARACTER*60 :: ctask,csave		! Evaluation flag

  REAL :: wkspc(ngrd(1)*ngrd(2)*ngrd(3)*(2*msave+4)+ &
	        12*msave*msave+12*msave)! Working space
  REAL :: bdlow(ngrd(1)*ngrd(2)*ngrd(3))! Lower bounds
  REAL :: bdupp(ngrd(1)*ngrd(2)*ngrd(3))! Upper bounds
  REAL :: factr,dsave(29)

  INTEGER :: iprnt,isbmn,isave(44)
  INTEGER :: nbund(ngrd(1)*ngrd(2)*ngrd(3)) ! Bound flags
  INTEGER :: iwrka(3*ngrd(1)*ngrd(2)*ngrd(3))

  LOGICAL :: lsave(4)

  !** End of LBFGS_B declarations.

  ! Local variables:
  INTEGER :: itr,sts,nvr,i,j,k
  REAL :: fcn,fnp,fnm,ggg,eps
  REAL :: gdt(ngrd(1),ngrd(2),ngrd(3))	! Gradients
  REAL :: grd(ngrd(1),ngrd(2),ngrd(3))  ! Grid function
  REAL :: egd(ngrd(1),ngrd(2),ngrd(3))  ! Grid function

  ! Start LBFGS_B:
  ctask = 'START'
  nvr = ngrd(1)*ngrd(2)*ngrd(3)		! Number of controls
  factr = 1.0d2
  iprnt = 1
  isbmn = 1

  ! No bounds for controls:
  nbund = 0

  ! Initial:
  itr = 0
  grd = sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3))

  ! Iterations:
1 CONTINUE

  CALL LBFGSB(nvr,msave,grd,bdlow,bdupp,nbund,fcn,gdt,factr, &
	wkspc,iwrka,ctask,iprnt,isbmn,csave,lsave,isave,dsave)

  ! Exit iteration if succeed:
  IF (ctask(1:11) .EQ. 'CONVERGENCE') GOTO 2

  ! Function and gradient values are needed:
  IF (ctask(1:2) .EQ. 'FG') THEN

    ! Function value:
    CALL Functn3D(fcn,grd,ngrd,obsv,nobs,indx,coef,wght)

    ! Gradient values:
    CALL Gradnt3D(gdt,grd,ngrd,obsv,nobs,indx,coef,wght)

    ! Check gradients:
    ! i=int(0.5*ngrd(1))
    ! j=int(0.3*ngrd(2))
    ! k=3
    ! eps = 1.0e-1
    ! egd = grd
    ! egd(i,j,k) = egd(i,j,k)+eps
    ! CALL Functn3D(fnp,egd,ngrd,obsv,nobs,indx,coef,wght)
    ! egd(i,j,k) = egd(i,j,k)-2.0*eps
    ! CALL Functn3D(fnm,egd,ngrd,obsv,nobs,indx,coef,wght)
    ! ggg = (fnp-fnm)*0.5/eps
    ! WRITE(*,9) ggg,gdt(i,j,k),i,j,k

  ENDIF
9 FORMAT('STMAS>STMASMin: Gradient check: ',2E12.4,3I6)

  ! Exit if irregularity is encountered:
  IF ((ctask(1:2) .NE. 'FG') .AND. (ctask(1:5) .NE. 'NEW_X')) THEN
    WRITE(*,*) 'STMAS>STMASMin: Irregularity termination of LBFGSB'
    GOTO 2
  ENDIF

  ! A new iteration point:
  IF (ctask(1:5) .EQ. 'NEW_X') itr = itr+1
  IF (itr .LT. maxitr) GOTO 1

  ! Exit of LBFGSB iteration:
2 CONTINUE

  ! Save the solution:
  sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = grd

END SUBROUTINE Minimize
