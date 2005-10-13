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
    sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3):2) = 0.5*( &
      sltn(1:ngrd(1),1:ngrd(2),1:ngrd(3)-2:2)+ &
      sltn(1:ngrd(1),1:ngrd(2),2:ngrd(3)  :2) )

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
