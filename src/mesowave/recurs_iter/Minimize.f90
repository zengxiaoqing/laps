SUBROUTINE Minimize(id,ds)

!*********************************************************
!  This routine minimizes the cost function derived from
!  a surface analysis using recursive filter.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!           Sep. 2004 by YUANFU XIE penalizing div.
!*********************************************************

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: id
  REAL,    INTENT(IN) :: ds(3)

  ! LBFGS_B variables:
  INCLUDE 'LBFGSB.f90'

  ! Local variables:
  INTEGER          :: itr,nvar,idp,istatus
  ! Adjoint variables:
  DOUBLE PRECISION :: f,adf,dv(mvar),g(mvar)
  REAL             :: ada(mx,my,mt,2),rv(mvar)
  
  ! Unified u/v analysis:
  idp = id
  IF (id .EQ. 201) idp = id+1

  ! Start LBFGS_B
  ctask = 'START'
  factr = 1.0d+2
  iprnt = 1
  isbmn = 1

  ! Allocate memory:
  ALLOCATE(bdlow(mvar),bdupp(mvar),nbund(mvar),iwrka(3*mvar), &
	wk(mvar*(2*msave+4)+12*msave*msave+12*msave), &
	STAT=istatus)
  IF (istatus .NE. 0) THEN
     PRINT*,'Minimize: no space for LBFGSB workspace'
     STOP
  ENDIF

  nbund = 0

  ! Initial:
  itr = 0
  nvar = n(1)*n(2)*n(3)*(idp-id+1)
  rv(1:nvar) = RESHAPE(a(1:n(1),1:n(2),1:n(3),id:idp),(/nvar/))
  dv(1:nvar) = rv(1:nvar)

  ! Looping:
1 CONTINUE
  CALL LBFGSB(nvar,msave,dv,bdlow,bdupp,nbund,f,g,factr,wk, &
              iwrka,ctask,iprnt,isbmn,csave,lsave,isave,dsave)

  ! Transfer a vector to a grid:
  rv(1:nvar) = dv(1:nvar)
  a(1:n(1),1:n(2),1:n(3),id:idp) = RESHAPE(rv(1:nvar), &
                               (/n(1),n(2),n(3),idp-id+1/))
  ! Exit if succeed:
  IF (ctask(1:11) .EQ. 'CONVERGENCE') GOTO 2

  ! Function and gradient are needed:
  IF (ctask(1:2) .EQ. 'FG') THEN

     ! Function value:
     IF ((id .NE. 201) .AND. (id .NE. 301)) THEN
        CALL functn(f,a(1,1,1,id),l,n,id,np,al)
     ELSE
        ! CALL functndiv(f,a(1,1,1,id),l,n,ds,id,np,al)
     ENDIF

     ! Gradient value:
     adf = 1.0d0
     ada = 0.0
     IF (id .NE. 201) THEN
        CALL adfunctn( a(1,1,1,id), l, n, id, np, al, adf, ada )
     ELSE
	! CALL adfunctndiv( a(1,1,1,id), l, n, ds, id, np, al, adf, ada )
     ENDIF

     rv(1:nvar) = RESHAPE(ada(1:n(1),1:n(2),1:n(3),1:idp-id+1), &
	(/ nvar /))
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
  IF (id .EQ. 201) &
     CALL RF3D(a(1,1,1,idp),l,n,al(1,idp),np(1,idp))

  ! Deallocate memory:
  DEALLOCATE(bdlow,bdupp,nbund,iwrka,wk,STAT=istatus)
  
END SUBROUTINE Minimize
