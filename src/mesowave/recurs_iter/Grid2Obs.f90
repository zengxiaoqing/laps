SUBROUTINE Grid2Obs

!*************************************************
!  This routine maps grid functions to observation
!  sites.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*************************************************

  IMPLICIT NONE

  INTEGER :: iobs,ier

  DO iobs=1,nobs

     CALL Intplt3d(o(2,iobs),l,n,d,dm, &
                   idx(1,iobs),coe(1,iobs),ier)

     IF (ier .NE. 0) THEN
        w(iobs) = 0.0
        idx(1:3,iobs) = 1
        coe(1:3,iobs) = 0.0
     ENDIF

  ENDDO

END SUBROUTINE Grid2Obs
