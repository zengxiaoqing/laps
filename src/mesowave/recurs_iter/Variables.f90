!****************************************************
!  This part of definition module defines all global
!  variables in this data assimilation.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!****************************************************

!===================
!  File names:
!===================

CHARACTER*60 :: datafile
INTEGER      :: namelens

!===================
!  Parameters:
!===================

INTEGER      :: maxitr

!===================
!  Arrays:
!===================

INTEGER      :: l(4),n(4),np(1:3,1:mv),nrf(mv)
REAL         :: a(mx,my,mt,mv),dm(2,3),d(3),al(1:3,1:mv)
REAL         :: s(mx,my,mt,mv),qc_cons(mv)
