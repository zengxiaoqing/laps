PROGRAM MAIN

!*********************************************************
!  This program analyzes a set of surface in time using a
!  recursive filter.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*********************************************************

  USE Definition
  USE Initialize
  USE Minimizatn
  USE ConfigLaps

  IMPLICIT NONE

  ! Local variable:
  INTEGER :: id

  CALL LapsInfo
  CALL LSO_Data

  ! CALL Namelist

  ! CALL ReadObsn

  CALL Grid2Obs

  !CALL Minimize
  DO id=1,n(4)
     CALL Iterates(id,bkgd,ldf,nx,ny,ncycles,nvlaps,nfic)
     PRINT*,'Variable ',id,' has been analyzed'
  ENDDO

  ! CALL Writeout
  CALL Writeout
  CALL WriteAnalysis(s(1:n(1),1:n(2),1:n(3),1:n(4)),n)

END PROGRAM MAIN
