SUBROUTINE LAPS_DIVIDER

!==============================================================================
!  THIS ROUTINE PRINTS A LINE DIVIDER TO FORMAT AN OUTPUT.
!
!  HISTORY:
!	CREATION:	YUANFU XIE	JUN 2007
!==============================================================================

  CHARACTER :: SYMBOL
  INTEGER   :: LOOPVR,NCOUNT

  SYMBOL = '='
  NCOUNT = 70

  WRITE(6,*) (SYMBOL,LOOPVR=1,NCOUNT)

END SUBROUTINE LAPS_DIVIDER
