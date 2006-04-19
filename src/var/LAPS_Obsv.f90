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

SUBROUTINE LAPS_Obsv

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
  INTEGER :: useable_radar,istatus_remap_pro,status
  INTEGER :: INIT_TIMER,ntmin,ntmax
  REAL*4 :: u(n(1),n(2),n(3)),v(n(1),n(2),n(3))
  REAL*4 :: ulaps(n(1),n(2),n(3)),vlaps(n(1),n(2),n(3))
  REAL*4 :: wt(n(1),n(2),n(3)),weight_prof


  ! Profiler Data:
  ntmin = -1
  ntmax = +1

  status = INIT_TIMER()

  CALL get_wind_3d_obs(n(1),n(2),n(3),rmissing,imissing, &
	i4time,height3d,height1d,max_pr,max_pr_lvls,     &
	weight_prof,l_raob,l_cdw,n_sao,n_pirep,lat,lon,  &
	ntmin,ntmax,u,v,ulaps,vlaps,wt,maxobs,obs_point, &
	nobs_point,rlat_radar,rlon_radar,rhgt_radar,     &
	useable_radar,n_grid_vel,grid_radar_vel,         &
	istatus_remap_pro,status)

END SUBROUTINE LAPS_Obsv
