MODULE Initialize

!****************************************************
!  This module initializes the GPS data assimilation
!  package.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!****************************************************

  USE Definition
  USE Util_Tools

CONTAINS

  INCLUDE 'Namelist.f90'
  INCLUDE 'ReadObsn.f90'
  INCLUDE 'Grid2Obs.f90'

END MODULE Initialize
