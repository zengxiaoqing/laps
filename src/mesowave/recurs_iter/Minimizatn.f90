MODULE Minimizatn

!*********************************************************
!  This module defines a minimization procedure analyzing
!  surface observations in time sequence.
!
!  HISTORY: JAN. 2004 by YUANFU XIE.
!*********************************************************

  USE Definition

CONTAINS
  
  INCLUDE 'Functn.f90'
  INCLUDE 'Functn_ad.f90'
  ! INCLUDE 'FunctnDiv.f90'
  ! INCLUDE 'FunctnDiv_ad.f90'
  INCLUDE 'Iterates.f90'
  INCLUDE 'Minimize.f90'
  INCLUDE 'RF1D.f90'
  INCLUDE 'RF1D_ad.f90'
  INCLUDE 'RF3D.f90'
  INCLUDE 'RF3D_ad.f90'

END MODULE Minimizatn
