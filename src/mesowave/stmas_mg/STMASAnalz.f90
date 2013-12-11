!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

MODULE STMASAnalz

!==========================================================
!  This module contains STMAS variational analysis based
!  a multigrid fitting technique.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modification:
!                 10-2008 by min-ken hsieh
!                 11-2008 by min-ken hsieh
!                 STMASAna
!                 Functn3D
!                 Gradnt3D
!==========================================================

  USE Definition

CONTAINS

SUBROUTINE STMASAna(anal,ngrd,dxyt,domn,bkgd,nfrm, &
		    obsv,nobs,wght,stna,ospc,indx, &
		    coef, bund,ipar,rpar,vnam,pnlt,&
                    slvl,ucvr,diag)
		    

!==========================================================
!  This routine reads into necessary namelists for STMAS
!  analysis through a multigrid technique.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modification:
!                 08-2008 by min-ken hsieh
!                 add parameter stna to map each obs its stn name for STMASVer
!       Modification:
!                 10-2008 by min-ken hsieh
!                 bound option only apply to last level
!                 pass in penalty for each var
!                 pass in slevel for each var
!       Modification:
!                 11-2008 by min-ken hsieh
!                 pass in uncovr for each var
!                 call IntplBkg to interpolate bkg to multigrid
!                 and pass interpolated bkg and uncover (hbg/huc) to Minimize
!       Modification:
!                 12-2008 by Yuanfu Xie
!                 pass in diag(nol) array for J_b term.
!       Modification:
!                 01-2009 by Yuanfu Xie
!                 calculate number of multigrid levels using slvl
!                 assuming level 1 with gridpoint 3 3 3.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3)	! Grid numbers
  INTEGER, INTENT(IN) :: nfrm		! Bkgd time frames
  INTEGER, INTENT(INOUT) :: nobs	! Number of obs
  INTEGER, INTENT(IN) :: indx(6,nobs)	! Indices (intepolate)
  INTEGER, INTENT(IN) :: bund		! Bound constraints
  INTEGER, INTENT(IN) :: ipar(1)	! Integer parameter
  INTEGER, INTENT(IN) :: slvl           ! Level to start analysis. added by min-ken
  REAL, INTENT(IN) :: dxyt(3)		! Grid spacing
  REAL, INTENT(IN) :: domn(2,3)		! Domain
  REAL, INTENT(IN) :: bkgd(ngrd(1),ngrd(2),nfrm)
  REAL, INTENT(INOUT) :: obsv(4,nobs)
  REAL, INTENT(INOUT) :: wght(nobs)	! Obs weightings
  REAL, INTENT(IN) :: pnlt		! penalty for each variable
  CHARACTER*20, INTENT(INOUT):: stna(nobs)
					! Obs station name by min-ken hsieh
  CHARACTER*4, INTENT(IN):: vnam        ! variable name(used in lvl addition and output each level result)
  REAL, INTENT(IN) :: ospc(3)		! Obs spacing
  REAL, INTENT(IN) :: coef(6,nobs)	! Coeffients
  REAL, INTENT(IN) :: rpar(1)		! Real parameter
  REAL, INTENT(IN) :: diag(ngrd(1),ngrd(2))! Diagnol array for J_b

  REAL, INTENT(INOUT) :: anal(ngrd(1),ngrd(2),ngrd(3))

  LOGICAL, INTENT(IN) :: ucvr(ngrd(1),ngrd(2),ngrd(3)) !uncovered array by min-ken

  ! Local variables:
  INTEGER :: lvl(3)			! Number of multigrid levels
  INTEGER :: lvc(3)			! Account of the levels
  INTEGER :: mgd(3)			! Multigrid points
  INTEGER :: mld(3)			! Leading dimensions
  INTEGER :: inc(3)			! Grid change increment
  INTEGER :: mcl			! Number of multigrid cycles
  INTEGER :: mlv			! Maximum number levels
  INTEGER :: idx(6,nobs)		! Indices at a multigrid
  INTEGER :: i,j,k,l,ier,ii,ij,ik
  REAL :: dis,rsz			! Distance and resizes
  REAL :: dgd(3)			! Grid spacing at a multigrid
  REAL :: coe(6,nobs)			! Coefficients at a multigrid
  REAL, ALLOCATABLE, DIMENSION(:,:,:) &
       :: sln				! Multigrid solutions
  INTEGER, ALLOCATABLE, DIMENSION(:,:,:) &
       :: huc				! Interpolated uncover by min-ken


  ! Number of multigrid V-cycles:
  mcl = 0

  ! Calculate total multigrid levels needed:
  DO i=1,3
    lvl(i) = 0
    ! Initial grid distance 3 gridpoints:
    dis = (domn(2,i)-domn(1,i))/2.0
1   CONTINUE
    IF (dis .GT. ospc(i)) THEN
      lvl(i) = lvl(i)+1
      dis = dis*0.5
      GOTO 1
    ENDIF
    ! LVL counts addition multigrid levels except the current 3x3x3 level
    lvl(i) = MIN0(lvl(i),INT(ALOG(FLOAT(ngrd(i)-1))/ALOG(2.0))-1)
  ENDDO
  IF (verbal .EQ. 1) WRITE(*,2) lvl(1:3)
2 FORMAT('STMASAna: Number of levels: ',3I3)

  ! Leading dimensions according to the total levels:
  mld = 2**(lvl+1)+1

  ! Allocate memory for mutligrid solutions:
  ALLOCATE(sln(mld(1),mld(2),mld(3)), STAT=ier)
  
  ! Allocate interpolated bkg and uncover
  ALLOCATE(huc(mld(1),mld(2),mld(3)), STAT=ier)
  
  ! Start multigrid analysis:

  ! modified by min-ken hsieh
  ! mgd should be determinated by slvl,
  ! but it cannot be bigger than mld
  DO i=1,3
    mgd(i) = MIN0(2**(slvl-1)+1,mld(i))
  ENDDO
  ! Maximum multigrid levels:
  mlv = MAXVAL(lvl(1:3))-slvl+1

  ! Multigrid cycles:
  DO l=1,mcl*2+1

    ! Down/Up cycle:
    IF (MOD(l,2) .EQ. 1) THEN
      rsz = 2.0
    ELSE
      rsz = 0.5
    ENDIF

    ! Through all possible levels:
    ! modified by min-ken hsieh
    ! because we start from slvl
    ! lvc are no longer starting from 0
    !lvc = 0
    lvc = slvl - 1

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
      CALL Grid2Obs(idx,coe,obsv,nobs,wght,stna,mgd,dgd,domn)

      ! added by min-ken hsieh
      ! Interpolate a background to multigrid
      CALL IntplBkg(ucvr,ngrd,huc,mgd,mld,dgd,dxyt,domn)

      ! Initial guesses:
      IF ((l .EQ. 1) .AND. (k .EQ. 1)) &
	!sln = 0.5*(MAXVAL(obsv(1,1:nobs))-MINVAL(obsv(1,1:nobs)))
	sln = 0.0

      ! Down and up cycle:
      IF (MOD(l,2) .EQ. 0) THEN
	CALL Projectn(sln,mld,mgd,inc)	! Up cycle
      ELSE
	CALL Interpln(sln,mld,mgd,inc)	! Down cycle
      ENDIF

      ! added by min-ken hsieh
      ! Smooth:
      !CALL Smoother(sln,mld,mgd)

      ! Analyzing:
      !bound the pcp analysis on finest level
      IF (k.EQ.mlv) THEN
        CALL Minimize(sln,mld,mgd,obsv,nobs,idx,coe,wght,bund,pnlt,huc)
      ELSE
        CALL Minimize(sln,mld,mgd,obsv,nobs,idx,coe,wght,0,pnlt,huc)
      ENDIF

    ENDDO
  ENDDO

  ! Map the multigrid solution to the grid requested:
  CALL Mul2Grid(sln,mld,mgd,anal,ngrd,domn)

  ! Deallocate multigrid:
  DEALLOCATE(sln,STAT=ier)
  DEALLOCATE(huc,STAT=ier)

END SUBROUTINE STMASAna

SUBROUTINE Smoother(sltn,mlds,ngrd)

!==========================================================
!  This routine smooth a grid with its neighbor values 
!
!  HISTORY:
!	Creation: 9-2008 by min-ken hsieh
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mlds(3),ngrd(3)
  REAL, INTENT(INOUT) :: sltn(mlds(1),mlds(2),mlds(3))


  ! smooth:
  ! x-direction
  sltn(2:ngrd(1)-1,1:ngrd(2),1:ngrd(3)) = &
    0.5*sltn(1:ngrd(1)-2,1:ngrd(2),1:ngrd(3))+ &
    0.5*sltn(3:ngrd(1),1:ngrd(2),1:ngrd(3)) 
  ! y-direction
  sltn(1:ngrd(1),2:ngrd(2)-1,1:ngrd(3)) = &
    0.5*sltn(1:ngrd(1),1:ngrd(2)-2,1:ngrd(3))+ &
    0.5*sltn(1:ngrd(1),3:ngrd(2),1:ngrd(3)) 
  ! t-direction
  sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1) = &
    0.5*sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2)+ &
    0.5*sltn(1:ngrd(1),1:ngrd(2),3:ngrd(3)) 

END SUBROUTINE Smoother

SUBROUTINE Projectn(sltn,mlds,ngrd,incr)

!==========================================================
!  This routine projects a grid value function to a coarse
!  resolution whose numbers of gridpoints are ngrd.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mlds(3),incr(3),ngrd(3)
  REAL, INTENT(INOUT) :: sltn(mlds(1),mlds(2),mlds(3))

  ! Local variables:
  INTEGER :: mgd(3)

  mgd = incr*(ngrd-1)+1

  ! Projection based grid increments:
  sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)) = &
    sltn(1:mgd(1):incr(1),1:mgd(2):incr(2),1:mgd(3):incr(3))

END SUBROUTINE Projectn

SUBROUTINE Interpln(sltn,mlds,ngrd,incr)

!==========================================================
!  This routine interpolates a grid value function to a 
!  fine resolution with ngrd gridpoints.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mlds(3),ngrd(3),incr(3)
  REAL, INTENT(INOUT) :: sltn(mlds(1),mlds(2),mlds(3))

  ! Local variables:
  INTEGER :: mgd(3)

  mgd = (ngrd-1)/incr+1

  ! Projection based grid increments:
  sltn(1:ngrd(1):incr(1),1:ngrd(2):incr(2),1:ngrd(3):incr(3)) = &
    sltn(1:mgd(1),1:mgd(2),1:mgd(3))

  ! Interpolate:
  
  ! X direction:
  IF (incr(1) .EQ. 2) &
    sltn(2:ngrd(1)-1:2,1:ngrd(2):incr(2),1:ngrd(3):incr(3)) = &
      0.5*( &
      sltn(1:ngrd(1)-2:2,1:ngrd(2):incr(2),1:ngrd(3):incr(3))+ &
      sltn(3:ngrd(1)  :2,1:ngrd(2):incr(2),1:ngrd(3):incr(3)) )

  ! Y direction:
  IF (incr(2) .EQ. 2) &
    sltn(1:ngrd(1),2:ngrd(2)-1:2,1:ngrd(3):incr(3)) = 0.5*( &
      sltn(1:ngrd(1),1:ngrd(2)-2:2,1:ngrd(3):incr(3))+ &
      sltn(1:ngrd(1),3:ngrd(2)  :2,1:ngrd(3):incr(3)) )

  ! T direction:
  IF (incr(3) .EQ. 2) &
    sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3)-1:2) = 0.5*( &
      sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2:2)+ &
      sltn(1:ngrd(1),1:ngrd(2),3:ngrd(3)  :2) )

END SUBROUTINE Interpln

SUBROUTINE Mul2Grid(sltn,mled,mgrd,anal,ngrd,domn)

!==========================================================
!  This routine projects a grid value function to a coarse
!  resolution.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mled(3),mgrd(3),ngrd(3)
  REAL, INTENT(IN) :: sltn(mled(1),mled(2),mled(3)),domn(2,3)
  REAL, INTENT(OUT) :: anal(ngrd(1),ngrd(2),ngrd(3))

  ! Local variables:
  INTEGER :: i,j,k,idx(6),ix,iy,it,ier
  REAL :: xyt(3),dsm(3),dsn(3),coe(6)

  ! Grid spacings:
  dsm = (domn(2,1:3)-domn(1,1:3))/FLOAT(mgrd-1)
  dsn = (domn(2,1:3)-domn(1,1:3))/FLOAT(ngrd-1)

  ! Interpolate:
  DO k=1,ngrd(3)
    xyt(3) = domn(1,3)+(k-1)*dsn(3)
    DO j=1,ngrd(2)
      xyt(2) = domn(1,2)+(j-1)*dsn(2)
      DO i=1,ngrd(1)
	xyt(1) = domn(1,1)+(i-1)*dsn(1)

	CALL Intplt3d(xyt,mgrd,dsm,domn,idx,coe,ier)

	! Evaluate:
	anal(i,j,k) = 0.0
	DO it=3,6,3
	  DO iy=2,5,3
	    DO ix=1,4,3
	      anal(i,j,k) = anal(i,j,k)+ &
		sltn(idx(ix),idx(iy),idx(it))* &
		coe(ix)*coe(iy)*coe(it)
	    ENDDO
	  ENDDO
	ENDDO

      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE Mul2Grid

SUBROUTINE Minimize(sltn,mlds,ngrd,obsv,nobs,indx,coef, &
		    wght,bund,pnlt,huc)

!==========================================================
!  This routine minimizes the STMAS cost function.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       Modified: 11-2008 by min-ken hsieh.
!                 pass in hbg,huc to calcualte Jb
!       Modified: 12-2013 by YUANFU XIE:
!                 New LBFGSB 3.0 fixes machine eps calculation
!                 Yuanfu adjust LBFGSB 3.0 for solving a super
!                 large minimization.
!                 a) Switch to LBFGSB.3.0 with Yuanfu's
!                    modification for super large minimization
!                 b) Adjust wkspc dimension for new LBFGSB and
!                    add a new parameter pg for new LBFGSB
!                 c) Change automatic arrays of wkspc, wrka, 
!                    bdlow, bdupp, nbund to allocatables
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: mlds(3),ngrd(3),nobs,indx(6,nobs)
  INTEGER, INTENT(IN) :: bund
  REAL, INTENT(IN) :: obsv(4,nobs),coef(6,nobs),wght(nobs)
  REAL, INTENT(INOUT) :: sltn(mlds(1),mlds(2),mlds(3))
  REAL, INTENT(IN) :: pnlt	!penalty for each variable

  INTEGER, INTENT(IN) :: huc(mlds(1),mlds(2),mlds(3))	! interpolated uncover array by min-ken

  !** LBFGS_B variables:
  INTEGER, PARAMETER :: msave=7		! Max iter save

  CHARACTER*60 :: ctask,csave		! Evaluation flag

  REAL, ALLOCATABLE :: wkspc(:),bdlow(:),bdupp(:)
  REAL :: factr,pg,dsave(29)            ! Yuanfu added pg for using setulb

  INTEGER :: iprnt,isbmn,isave(44)
  INTEGER, ALLOCATABLE :: nbund(:) ! Bound flags
  INTEGER, ALLOCATABLE :: iwrka(:)

  LOGICAL :: lsave(4)

  !** End of LBFGS_B declarations.

  ! Local variables:
  INTEGER :: itr,sts,nvr,i,j,k,istatus
  REAL :: fcn,fnp,fnm,ggg,eps
  REAL :: gdt(ngrd(1),ngrd(2),ngrd(3))	! Gradients
  REAL :: grd(ngrd(1),ngrd(2),ngrd(3))  ! Grid function
  REAL :: egd(ngrd(1),ngrd(2),ngrd(3))  ! Grid function
  
  ! Dec 2013: Yuanfu changed the wkspc dimension from (2*msave+4)
  ! to (2*msave+5) for both old and new versions of LBFGSB:
  allocate(wkspc(ngrd(1)*ngrd(2)*ngrd(3)*(2*msave+5)+ &
                12*msave*msave+12*msave), &
           iwrka(3*ngrd(1)*ngrd(2)*ngrd(3)), STAT=istatus)
  allocate(bdlow(ngrd(1)*ngrd(2)*ngrd(3)), &
           bdupp(ngrd(1)*ngrd(2)*ngrd(3)), &
           nbund(ngrd(1)*ngrd(2)*ngrd(3)), STAT=istatus)
           

  ! Start LBFGS_B:
  ctask = 'START'
  nvr = ngrd(1)*ngrd(2)*ngrd(3)		! Number of controls
  factr = 1.0d2
  pg = 1.0e-4        ! Yuanfu added for using setulb
  iprnt = 1
  isbmn = 1

  ! Simple bound constraints for controls (only low bound used here):
  nbund = bund
  if (bund .EQ. 1) bdlow = 0.0

  ! Initial:
  itr = 0
  grd = sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3))

  ! Iterations:
1 CONTINUE

  CALL SETULB(nvr,msave,grd,bdlow,bdupp,nbund,fcn,gdt,factr, &
	      pg,wkspc,iwrka,ctask,iprnt,csave,lsave,isave,dsave)

  ! Exit iteration if succeed:
  IF (ctask(1:11) .EQ. 'CONVERGENCE') GOTO 2

  ! Function and gradient values are needed:
  IF (ctask(1:2) .EQ. 'FG') THEN

    ! Function value:
    CALL Functn3D(fcn,grd,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

    ! Gradient values:
    CALL Gradnt3D(gdt,grd,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

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

SUBROUTINE Functn3D(fctn,grid,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

!==========================================================
!  This routine evaluates the STMAS cost function.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!       Modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       Modified: 11-2008 by min-ken hsieh.
!                 pass in bkg params to calcualte Jb
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3),mlds(3),nobs,indx(6,nobs)
  REAL, INTENT(IN) :: grid(ngrd(1),ngrd(2),ngrd(3))
  REAL, INTENT(IN) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  REAL, INTENT(IN) :: pnlt
  REAL, INTENT(OUT) :: fctn

  INTEGER, INTENT(IN) :: huc(mlds(1),mlds(2),mlds(3))	! uncover by min-ken

  ! Local variables:
  INTEGER :: i,j,k,io
  REAL    :: gvl,pnl

  ! Initial:
  fctn = 0.0
  penalt = pnlt

  ! Jb term:
  DO k=1,ngrd(3)
    DO j=1,ngrd(2)
      DO i=1,ngrd(1)
        fctn = fctn+huc(i,j,k)*(grid(i,j,k)**2)
      ENDDO
    ENDDO
  ENDDO

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

SUBROUTINE Gradnt3D(grdt,grid,ngrd,obsv,nobs,indx,coef,wght,pnlt,huc,mlds)

!==========================================================
!  This routine evaluates gradients of STMAS cost function.
!
!  HISTORY:
!	Creation: JUN. 2005 by YUANFU XIE.
!       Modified: 10-2008 by min-ken hsieh.
!                 pass in penalty for each var
!       Modified: 11-2008 by min-ken hsieh.
!                 pass in hbg, huc to calcualte Jb
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3),mlds(3),nobs,indx(6,nobs)
  REAL, INTENT(IN) :: grid(ngrd(1),ngrd(2),ngrd(3))
  REAL, INTENT(IN) :: coef(6,nobs),obsv(4,nobs),wght(nobs)
  REAL, INTENT(IN) :: pnlt
  REAL, INTENT(OUT) :: grdt(ngrd(1),ngrd(2),ngrd(3))

  INTEGER, INTENT(IN) :: huc(mlds(1),mlds(2),mlds(3))	! uncover array by min-ken

  ! Local variables:
  INTEGER :: i,j,k,io
  REAL :: gvl,pnl
  REAL :: lpl(ngrd(1),ngrd(2),ngrd(3))
  REAL :: gll(ngrd(1),ngrd(2),ngrd(3))

  ! Initial:
  grdt = 0.0
  penalt = pnlt

  ! Jb term:
  
  DO k=1,ngrd(3)
    DO j=1,ngrd(2)
      DO i=1,ngrd(1)
        grdt(i,j,k) = huc(i,j,k)*(grid(i,j,k))
      ENDDO
    ENDDO
  ENDDO


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

SUBROUTINE STMASInc

!==========================================================
!  This routine adds the analysis increments to background
!  fields using the LAPS land/water factors.
!
!  HISTORY:
!	Creation: 9-2005 by YUANFU XIE.
!==========================================================

  IMPLICIT NONE

  ! Local variables:
  INTEGER i,j,k,l

  DO l=1,numvar
    IF (needbk(l) .EQ. 1) THEN
      DO k=1,numtmf

        ! Add increment to the background:
	analys(1:numgrd(1),1:numgrd(2),k,l) = &
	  analys(1:numgrd(1),1:numgrd(2),k,l)+ &
	  bkgrnd(1:numgrd(1),1:numgrd(2),k,l)

	! Treat land factor:
!	analys(1:numgrd(1),1:numgrd(2),k,l) = &
!	  lndfac(1:numgrd(1),1:numgrd(2))* &
!	  analys(1:numgrd(1),1:numgrd(2),k,l)+ &
!	  (1.0-lndfac(1:numgrd(1),1:numgrd(2)))* &
!	  bkgrnd(1:numgrd(1),1:numgrd(2),k,l)
      ENDDO
    ENDIF
  ENDDO

END SUBROUTINE STMASInc

SUBROUTINE IntplBkg(ucvr,ngrd,huc,mgd,mld,dgd,dxyt,domn)

!==========================================================
!  This routine interpolate background field and uncover to
!  multigrid.
!
!  HISTORY:
!	Creation: 11-2008 by min-ken hsieh
!
!==========================================================

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: ngrd(3),mgd(3),mld(3)
  REAL, INTENT(IN) :: dgd(3),dxyt(3),domn(2,3)
  LOGICAL, INTENT(IN) :: ucvr(ngrd(1),ngrd(2),ngrd(3))

  INTEGER, INTENT(OUT) :: huc(mld(1),mld(2),mld(3))

  ! Local variables:
  INTEGER :: ier,i,j,k,ix,iy,it
  REAL :: xyt(3)
  INTEGER :: idx(6,mgd(1),mgd(2),mgd(3))
  REAL :: coe(6,mgd(1),mgd(2),mgd(3))
  LOGICAL :: lhuc(mld(1),mld(2),mld(3))

  ! Interpolation for each background grid:
  DO k=1,mgd(3)
    xyt(3) = domn(1,3)+(k-1)*dgd(3)
    DO j=1,mgd(2)
      xyt(2) = domn(1,2)+(j-1)*dgd(2)
      DO i=1,mgd(1)
        xyt(1) = domn(1,1)+(i-1)*dgd(1)
        CALL Intplt3D(xyt(1:3),ngrd,dxyt,domn, &
	    	      idx(1,i,j,k),coe(1,i,j,k),ier)

        ! Check:
        IF (ier .NE. 0) THEN
          PRINT*, 'BACKGROUND INTERPOLATION ERROR!!',i,j,k,xyt(1),xyt(2),xyt(3)
        ELSE
          ! Evaluate:
          huc(i,j,k) = 0
	  lhuc(i,j,k) = .TRUE.
	  DO it=3,6,3
	    DO iy=2,5,3
	      DO ix=1,4,3
		! see if this grid is uncovered
		lhuc(i,j,k) = lhuc(i,j,k) .AND. &
		  ucvr(idx(ix,i,j,k),idx(iy,i,j,k),idx(it,i,j,k))
		   
	      ENDDO
	    ENDDO
	  ENDDO

          IF(lhuc(i,j,k)) THEN
            huc(i,j,k) = 1
          ENDIF

        ENDIF
        
      ENDDO
    ENDDO
  ENDDO

END SUBROUTINE IntplBkg

END MODULE STMASAnalz
