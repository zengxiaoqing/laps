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
  INTEGER :: id,ierr
  REAL    :: ds(3)

  CALL LapsInfo
  CALL LSO_Data

  ! CALL Namelist

  ! CALL ReadObsn

  CALL Grid2Obs

  !CALL Minimize
  ds(1) = grid_spacingx
  ds(2) = grid_spacingy
  ds(3) = d(3)
  DO id=1,n(4)
     IF (id .NE. 4) THEN  ! do not analyze station pressure
        CALL Iterates(id,bkgd,ldf,nx,ny,ds,ncycles,nvlaps,nfic)
        PRINT*,'Variable ',id,' has been analyzed'
     ENDIF
  ENDDO

  ! Release memory of ldf:
  DEALLOCATE(bkgd,ldf,STAT=ierr)

  ! CALL Writeout
  ! CALL Writeout
  CALL WriteAnalysis(s(1:n(1),1:n(2),1:n(3),1:n(4)),n)

END PROGRAM MAIN
