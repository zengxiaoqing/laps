
!====================================================
!  This header file declares necessary variables for
!  the optimization package: LBFGSB.
!
!  Author: YUANFU XIE   NOAA/FSL
!  Date:   AUGUST 2000
!
!====================================================
      
!.....Variables required by LBFGS_B:

INTEGER, PARAMETER :: msave=7
INTEGER, PARAMETER :: mvar = mx*my*mt
CHARACTER*60       :: ctask,csave
DOUBLE PRECISION,ALLOCATABLE,DIMENSION(:) :: &
	wk,bdlow,bdupp
DOUBLE PRECISION   :: factr,dsave(29)
INTEGER            :: iprnt,isbmn,isave(44)
INTEGER,ALLOCATABLE,DIMENSION(:) :: nbund,iwrka
LOGICAL            :: lsave(4)

!.....End of LBFGS_B declarations.
