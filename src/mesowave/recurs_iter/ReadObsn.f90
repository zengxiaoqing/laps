SUBROUTINE ReadObsn

!*********************************************************
!  This routine reads in all observations avaiable on this
!  surface analysis.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*********************************************************

  IMPLICIT NONE

  INTEGER :: m

  OPEN(unit=10,file=datafile(1:namelens),status='old')

  ! Domain x time interval:
  READ(10,*) dm(1:2,1:3)
  d(1:3) = (dm(2,1:3)-dm(1,1:3))/FLOAT(n(1:3)-1)

  nobs = 1		! Count observations
  m = 0			! Count variables to be analyzed

1 FORMAT(i2,5e14.6)
2 CONTINUE
  READ(10,1,END=9) vid(nobs),o(1:4,nobs),w(nobs)
  IF (m .LT. vid(nobs)) m = vid(nobs)
  nobs = nobs+1
  GOTO 2

9 CONTINUE

  ! Close open file:
  CLOSE(10)

  nobs = nobs-1
  
  ! Check dimension for variables:
  IF (m .GT. mv) THEN
     PRINT*,'ReadObsn: Too much variables'
     STOP
  ENDIF
  IF (m .LE. 0) THEN
     PRINT*,'ReadObsn: No variables to analyze'
     STOP
  ENDIF
  ! Check if number of obs is fit into the array:
  IF (nobs .GT. mobs) THEN
     PRINT*,'ReadObsn: Too many obs'
     STOP
  ENDIF
  IF (nobs .LE. 0) THEN
     PRINT*,'ReadObsn: no obs to analyze'
     STOP
  ENDIF

END SUBROUTINE ReadObsn
