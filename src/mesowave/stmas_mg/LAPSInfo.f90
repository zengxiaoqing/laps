SUBROUTINE LAPSInfo

!==========================================================
!  This routine configures the STMAS for its analyses.
!  Tasks:
!	1. 
!
!  NOTE: three letter variables are local; six global;
!
!  HISTORY: 
! 	Creation: YUANFU XIE	6-2005
!==========================================================

  IMPLICIT NONE

  ! Local variables:
  CHARACTER*9 :: fnm
  INTEGER :: err	! Error indicator
  INTEGER :: i

  !*********************
  ! LAPS configuration:
  !*********************

  ! Get number of gridpoints:
  CALL GET_GRID_DIM_XY(numgrd(1),numgrd(2),err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting numgrd'

  ! Get LAPS cycle time:
  CALL GET_LAPS_CYCLE_TIME(lapsdt,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting cycle time'

  ! Get current system time:
  CALL GET_SYSTIME(i4time,fnm,err)

  ! Get a flag for missing data:
  CALL GET_R_MISSING_DATA(mising,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting mising flag'

  ! Get a flag for bad surface data:
  CALL GET_SFC_BADFLAG(badsfc,err)
  IF (err .ne. 1) PRINT*, 'STMAS>LAPSInfo: Error getting badsfc flag'

  CALL GET_MAXSTNS(mxstts,err)
  IF (err .ne. 1) print*, 'STMAS>LAPSInfo: Error getting maxstations'

  ! Check:
  IF (verbal .EQ. 1) THEN
    WRITE(*,1) numgrd(1:3),lapsdt,mxstts,mising,badsfc
  ENDIF
1 FORMAT('STMAS>LAPSInfo: Num  gridpoints: ',3I6,/, &
	 'STMAS>LAPSInfo: LAPS cycle time: ',I6,/, &
	 'STMAS>LAPSInfo: Maxnumber sites: ',I6,/, &
	 'STMAS>LAPSInfo: Missing/bad ids: ',e16.8,f16.5)

END SUBROUTINE LAPSInfo
