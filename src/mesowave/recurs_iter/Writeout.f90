SUBROUTINE Writeout

!*************************************************
!  This routine writes out the solution to a file.
!
!  HISTORY: FEB. 2004 by YUANFU XIE.
!*************************************************

  USE Definition

  IMPLICIT NONE

  OPEN(unit=10,file='analysis.dat')

  ! write(10,*) n

  ! WRITE(10,*) s(1:n(1),1:n(2),1:n(3),1:n(4))

  write(10,*) n(1),n(2),1,n(4)

  WRITE(10,*) s(1:n(1),1:n(2),n(3),1:n(4))

  CLOSE(10)

  !OPEN(unit=10,file='../dat/analysis.bin',form='unformatted')

  !write(10) n

  !WRITE(10) s(1:n(1),1:n(2),1:n(3),1:n(4))

  !CLOSE(10)

END SUBROUTINE Writeout
