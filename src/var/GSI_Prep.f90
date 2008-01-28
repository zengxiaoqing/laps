!dis    Forecast Systems Laboratory
!dis    NOAA/OAR/ERL/FSL
!dis    325 Broadway
!dis    Boulder, CO     80303
!dis
!dis    Forecast Research Division
!dis    Local Analysis and Prediction Branch
!dis    LAPS
!dis
!dis    This software and its documentation are in the public domain and
!dis    are furnished "as is."  The United States government, its
!dis    instrumentalities, officers, employees, and agents make no
!dis    warranty, express or implied, as to the usefulness of the software
!dis    and documentation for any purpose.  They assume no responsibility
!dis    (1) for the use of the software and documentation; or (2) to provide
!dis    technical support to users.
!dis
!dis    Permission to use, copy, modify, and distribute this software is
!dis    hereby granted, provided that the entire disclaimer notice appears
!dis    in all copies.  All modifications to this software must be clearly
!dis    documented, and are solely the responsibility of the agent making
!dis    the modifications.  If significant modifications or enhancements
!dis    are made to this software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis

PROGRAM GSI_Prep

!==========================================================
!  This program is a three-dimensional variational analysis
!  for LAPS. It can be either GSI or STMAS depending on the
!  variational namelist parameter.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!==========================================================

  IMPLICIT NONE

  ! LAPS configuration:
  CALL LAPS_Conf
  
  ! LAPS PREP:
  CALL LAPS_Bkgd

  ! LAPS OBS:
  CALL LAPS_Obsv

  ! Conversion to GSI:
  CALL GSI_BkObs

  ! Deallocate dynamic memory:
  CALL LAPS_Remv

  ! Radar Level III BUFR data:
  CALL GSI_radarbufr

  ! NOAA polar satellite 1D BUFR data: AMSU-A, AMSU-B, and HIRS
  ! CALL GSI_noaa1d2bufr

END PROGRAM GSI_Prep
