SUBROUTINE Namelist

!*********************************************************
!  This routine reads in namelist for the surface analysis
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*********************************************************

  USE Definition

  INTEGER :: i

  OPEN(unit=10,file='run/Namelist.txt',status='old')
  READ(10,*) ! skip headline

  l(1:4) = (/mx,my,mt,mv/)

  ! Grid dimensions:
  READ(10,*) n(1),n(2),n(3),n(4)

  ! Check:
  IF ((n(1) .GT. mx) .OR. (n(2) .GT. my) .OR. &
      (n(3) .GT. mt) .OR. (n(4) .GT. mv)) THEN
     PRINT*,'Namelist: Analysis array is too small!'
     STOP
  ENDIF

  ! Data filename:
  READ(10,*) namelens,datafile

  ! Recursive filters:
  READ(10,*) ! Skip a specification line
  DO i=1,n(4)
     READ(10,*) al(1:3,i),np(1:3,i)
  ENDDO

  ! Number of minimization iterations:
  READ(10,*) maxitr

  ! Number of recursive filter iterations:
  READ(10,*) nrf(1:n(4))

  CLOSE(10)

END SUBROUTINE Namelist
