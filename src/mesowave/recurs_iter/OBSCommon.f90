!***********************************************************
!  This header file defines large arrays for observations to
!  avoid the limitation of parameter passing.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!***********************************************************

! Maximum number of OBS:
INTEGER, PARAMETER :: mobs = 200000

! Observations:
INTEGER :: nobs
INTEGER :: vid(mobs),idx(3,mobs)
REAL    :: o(4,mobs),coe(3,mobs),w(mobs)
COMMON /OBSBlock/nobs,vid,idx,o,coe,w
