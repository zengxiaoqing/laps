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
	analys(1:numgrd(1),1:numgrd(2),k,l) = &
	  lndfac(1:numgrd(1),1:numgrd(2))* &
	  analys(1:numgrd(1),1:numgrd(2),k,l)+ &
	  (1.0-lndfac(1:numgrd(1),1:numgrd(2)))* &
	  bkgrnd(1:numgrd(1),1:numgrd(2),k,l)
      ENDDO
    ENDIF
  ENDDO

END SUBROUTINE STMASInc
