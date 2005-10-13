SUBROUTINE STMASAna(anal,ngrd,dxyt,domn,bkgd,nfrm, &
		    obsv,nobs,wght,ospc,indx,coef, &
		    ipar,rpar)
		    

!==========================================================
!  This routine reads into necessary namelists for STMAS
!  analysis through a multigrid technique.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)	! Grid numbers
  INTEGER, INTENT(IN) :: nfrm		! Bkgd time frames
  INTEGER, INTENT(INOUT) :: nobs		! Number of obs
  INTEGER, INTENT(IN) :: indx(6,nobs)	! Indices (intepolate)
  INTEGER, INTENT(IN) :: ipar(1)	! Integer parameter
  REAL, INTENT(IN) :: dxyt(3)		! Grid spacing
  REAL, INTENT(IN) :: domn(2,3)		! Domain
  REAL, INTENT(IN) :: bkgd(ngrd(1),ngrd(2),nfrm)
  REAL, INTENT(INOUT) :: obsv(4,nobs)
  REAL, INTENT(INOUT) :: wght(nobs)	! Obs weightings
  REAL, INTENT(IN) :: ospc(3)		! Obs spacing
  REAL, INTENT(IN) :: coef(6,nobs)	! Coeffients
  REAL, INTENT(IN) :: rpar(1)		! Real parameter

  REAL, INTENT(INOUT) :: anal(ngrd(1),ngrd(2),ngrd(3))

  ! Local variables:
  INTEGER :: lvl(3)			! Number of multigrid levels
  INTEGER :: lvc(3)			! Account of the levels
  INTEGER :: mgd(3)			! Multigrid points
  INTEGER :: mld(3)			! Leading dimensions
  INTEGER :: inc(3)			! Grid change increment
  INTEGER :: mcl			! Number of multigrid cycles
  INTEGER :: mlv			! Maximum number levels
  INTEGER :: idx(6,nobs)		! Indices at a multigrid
  INTEGER :: i,j,k,l,ier
  REAL :: dis,rsz			! Distance and resizes
  REAL :: dgd(3)			! Grid spacing at a multigrid
  REAL :: coe(6,nobs)			! Coefficients at a multigrid
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: sln				! Multigrid solutions

  ! Number of multigrid V-cycles:
  mcl = 0

  DO i=1,3
    lvl(i) = 0
    dis = 0.5*(domn(2,i)-domn(1,i))	! Start from 3 gridpoints
1   CONTINUE
    IF (dis .GT. ospc(i)) THEN
      lvl(i) = lvl(i)+1
      dis = dis*0.5
      GOTO 1
    ENDIF
    lvl(i) = MIN0(lvl(i),INT(ALOG(FLOAT(ngrd(i)-1))/ALOG(2.0))-2)
    lvl(i) = MAX0(lvl(i),0)
    IF (verbal .EQ. 1) WRITE(*,2) lvl(i)
  ENDDO
2 FORMAT('STMASAna: Number of levels: ',I3)

  ! Leading dimensions:
  mld = 2**(lvl+1)+1

  ! Allocate memory for mutligrid solutions:
  ALLOCATE(sln(mld(1),mld(2),mld(3)), STAT=ier)

  ! Start multigrid analysis:
  mgd = 2
  mlv = MAXVAL(lvl(1:3))

  ! Multigrid cycles:
  DO l=1,mcl*2+1

    ! Down/Up cycle:
    IF (MOD(l,2) .EQ. 1) THEN
      rsz = 2.0
    ELSE
      rsz = 0.5
    ENDIF

    ! Through all possible levels:
    lvc = 0
    DO k=1,mlv

      ! Redefine number of gridpoints:
      DO j=1,3
        IF (lvc(j) .LE. lvl(j)) THEN
	  mgd(j) = (mgd(j)-1)*rsz+1
	  inc(j) = 2			! grid change
	  lvc(j) = lvc(j)+1		! count levels
	ELSE
	  inc(j) = 1			! grid unchange
	ENDIF
      ENDDO

      ! Interpolate a multigrid to observations:
      dgd = (domn(2,1:3)-domn(1,1:3))/FLOAT(mgd-1)
      CALL Grid2Obs(idx,coe,obsv,nobs,wght,mgd,dgd,domn)

      ! Initial guesses:
      IF ((l .EQ. 1) .AND. (k .EQ. 1)) &
	sln = 0.5*(MAXVAL(obsv(1,1:nobs))-MINVAL(obsv(1,1:nobs)))

      ! Down and up cycle:
      IF (MOD(l,2) .EQ. 0) THEN
	CALL Projectn(sln,mld,mgd,inc)	! Up cycle
      ELSE
	CALL Interpln(sln,mld,mgd,inc)	! Down cycle
      ENDIF

      ! Analyzing:
      CALL Minimize(sln,mld,mgd,obsv,nobs,idx,coe,wght)

    ENDDO
  ENDDO

  ! Map the multigrid solution to the grid requested:
  CALL Mul2Grid(sln,mld,mgd,anal,ngrd,domn)

  ! Deallocate multigrid:
  DEALLOCATE(sln,STAT=ier)

END SUBROUTINE STMASAna
