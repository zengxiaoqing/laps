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

SUBROUTINE LAPS_Allc

!==========================================================
!  This routine allocates memory for LAPS dynamic arrays.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!==========================================================

  USE LAPS_Parm
  USE MEM_NAMELIST

  IMPLICIT NONE

  INTEGER :: status

  ! One dimension arrays:
  ALLOCATE(pressr1d(n(3)),height1d(n(3)),STAT=status)
  IF (status .NE. 0) THEN
    WRITE(6,*) 'LAPS_Allc: error to allocate memory for one-D arrays'

    ! Deallocate all previously allocated arrays:
    CALL LAPS_Remv
    STOP
  ENDIF
  ALLOCATE(rlat_radar(max_radars),rlon_radar(max_radars), &
	   rhgt_radar(max_radars),n_grid_vel(max_radars))

  ALLOCATE(grid_radar_vel(n(1),n(2),n(3),max_radars), STAT=status)
  IF (status .NE. 0) THEN
    WRITE(6,*) 'LAPS_Allc: error to allocate memory for 4-D arrays'

    ! Deallocate all previously allocated arrays:
    CALL LAPS_Remv
    STOP
  ENDIF

  ! Lat/lon and topography:
  ALLOCATE(lat(n(1),n(2)),lon(n(1),n(2)),topo(n(1),n(2)), &
	STAT=status)
  IF (status .NE. 0) THEN
    WRITE(6,*) 'LAPS_Allc: error to allocate memory for lat/lon/topo'

    ! Deallocate all previously allocated arrays:
    CALL LAPS_Remv
    STOP
  ENDIF

  ! 3D variables:
  ALLOCATE(height3d(n(1),n(2),n(3)),pressr3d(n(1),n(2),n(3)), &
	   temptr3d(n(1),n(2),n(3)),sphumd3d(n(1),n(2),n(3)), &
           u_wind3d(n(1),n(2),n(3)),v_wind3d(n(1),n(2),n(3)), &
	STAT=status)
  IF (status .NE. 0) THEN
    WRITE(6,*) 'LAPS_Allc: error to allocate memory for 3D variables'

    ! Deallocate all previously allocated arrays:
    CALL LAPS_Remv
    STOP
  ENDIF

END SUBROUTINE LAPS_Allc


SUBROUTINE LAPS_Remv

!==========================================================
!  This routine deallocates memory for LAPS dynamic arrays.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!==========================================================

  USE LAPS_Parm

  IMPLICIT NONE

  INTEGER :: status

  ! One dimension arrays:
  DEALLOCATE(pressr1d,height1d, STAT=status)

  DEALLOCATE(rlat_radar,rlon_radar,rhgt_radar,n_grid_vel, &
	     grid_radar_vel, STAT=status)

  ! Lat/lon and topography:
  DEALLOCATE(lat,lon,topo, STAT=status)

  ! 3D variables:
  DEALLOCATE(height3d,pressr3d,temptr3d,sphumd3d, &
    	     u_wind3d,v_wind3d, STAT=status)

END SUBROUTINE LAPS_Remv
