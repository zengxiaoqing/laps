SUBROUTINE Minimize(id)

!*********************************************************
!  This routine minimizes the cost function derived from
!  a surface analysis using recursive filter.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*********************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id

  ! LBFGS_B variables:
  INCLUDE 'LBFGSB.f90'

  ! Local variables:
  INTEGER          :: itr,nvar
  ! Adjoint variables:
  DOUBLE PRECISION :: f,adf,dv(mvar),g(mvar)
  REAL             :: ada(mx,my,mt),rv(mvar)
  

  ! Start LBFGS_B
  ctask = 'START'
  factr = 1.0d+2
  iprnt = 1
  isbmn = 1
  nbund = 0

  ! Initial:
  itr = 0
  nvar = n(1)*n(2)*n(3)
  rv(1:nvar) = RESHAPE(a(1:n(1),1:n(2),1:n(3),id),(/nvar/))
  dv(1:nvar) = rv(1:nvar)

  ! Looping:
1 CONTINUE
  CALL LBFGSB(nvar,msave,dv,bdlow,bdupp,nbund,f,g,factr,wk, &
              iwrka,ctask,iprnt,isbmn,csave,lsave,isave,dsave)

  ! Transfer a vector to a grid:
  rv(1:nvar) = dv(1:nvar)
  a(1:n(1),1:n(2),1:n(3),id) = RESHAPE(rv(1:nvar), &
                               (/n(1),n(2),n(3)/))
  ! Exit if succeed:
  IF (ctask(1:11) .EQ. 'CONVERGENCE') GOTO 2

  ! Function and gradient are needed:
  IF (ctask(1:2) .EQ. 'FG') THEN

     ! Function value:
     CALL functn(f,a(1,1,1,id),l,n,id,np,al)

     ! Gradient value:
     adf = 1.0d0
     ada = 0.0
     CALL adfunctn( a(1,1,1,id), l, n, id, np, al, adf, ada )

     rv(1:nvar) = RESHAPE(ada(1:n(1),1:n(2),1:n(3)), (/ nvar /))
     g(1:nvar) = rv(1:nvar)

  ENDIF

  ! Exit if irregularity is encountered:
  IF ((ctask(1:2) .NE. 'FG') .AND. (ctask(1:5) .NE. 'NEW_X')) THEN
     PRINT*,'Error in LBFGS_B'
     GOTO 2
  ENDIF

  ! Count number of iteration and compute relative error:
  IF (ctask(1:5) .EQ. 'NEW_X') THEN
     itr = itr+1
     ! print*,''
  ENDIF

  ! If number of iterations of LBFGS_B exceeds the limit:
  IF (itr .LT. maxitr) GOTO 1

  ! When an error in LBFGS_B occurs: exit
2 CONTINUE

  ! Convert to analysis:
  CALL RF3D(a(1,1,1,id),l,n,al(1,id),np(1,id))
  
END SUBROUTINE Minimize
