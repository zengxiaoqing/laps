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

SUBROUTINE LAPS_Bkgd

!==========================================================
!  This routine retrieves background fields by LAPS ingest.
!  The 3 D background fields are:
!	height, pressure, temperature, specific humidity,
!       u wind, v wind.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!==========================================================

  USE LAPS_Parm 

  IMPLICIT NONE
  
  ! Local variables:
  CHARACTER*3 :: varname
  INTEGER :: status

  ! Height:
  varname = 'HT'
  CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
		height3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving height'
    CALL LAPS_Remv
    STOP
  ENDIF

  ! Pressure:
  CALL get_pres_3d(i4time,n(1),n(2),n(3),pressr3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving pressure'
    CALL LAPS_Remv
    STOP
  ENDIF
  CALL get_pres_1d(i4time,n(3),pressr1d,status)

  ! Temperature:
  varname = 'T3'
  CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
		temptr3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving temperature'
    CALL LAPS_Remv
    STOP
  ENDIF

  ! Specific:
  ! varname = 'RH'
  ! CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
  !		sphumd3d,status)
  ! IF (status .NE. 1) THEN
  !  WRITE(6,*) 'LAPS_Pars: error retrieving relative humidity'
  !  CALL LAPS_Remv
  !  STOP
  !ENDIF
  varname = 'SH'
  CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
		sphumd3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving specific humidity'
    CALL LAPS_Remv
    STOP
  ENDIF

  ! U wind:
  varname = 'U3'
  CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
		u_wind3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving U wind'
    CALL LAPS_Remv
    STOP
  ENDIF

  ! V wind:
  varname = 'V3'
  CALL get_modelfg_3d(i4time,varname,n(1),n(2),n(3), &
		v_wind3d,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Pars: error retrieving V wind'
    CALL LAPS_Remv
    STOP
  ENDIF

END SUBROUTINE LAPS_Bkgd
