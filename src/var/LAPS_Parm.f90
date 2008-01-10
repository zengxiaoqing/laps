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

MODULE LAPS_Parm

!==========================================================
!  This module defines all necessary constant and variables
!  for LAPS parameters.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!========================================================== 

  IMPLICIT NONE

  INTEGER, PARAMETER :: maxxobs = 100000	! Maximum obs limit

  include 'barnesob.inc'

  CHARACTER*16 :: asctime

  INTEGER :: n(3)		! Spatial dimensions
  INTEGER :: timelen		! Time interval
  INTEGER :: i4time		! System time
  INTEGER :: imissing		! Integer missing data
  INTEGER :: nobs		! Number of observations
  INTEGER :: n_tobs		! Number of temp observations
  INTEGER :: nobs_point		! Number of point observations
  INTEGER :: n_meso		! Number of mesonet data
  INTEGER :: n_sao		! Number of sfc aviation obs
  ! INTEGER :: n_pirep		! Number of pilot report obs
  ! INTEGER :: max_pr		! Max profiles
  INTEGER :: max_pr_lvls	! Max levels
  ! INTEGER :: max_radars		! Maximum number of radars
  INTEGER :: maxtobs		! Maximum number of temp obs
  INTEGER :: istat_radar_vel

  LOGICAL :: l_raob, l_cdw, l_radial
  ! LOGICAL :: l_use_raob

  REAL :: rmissing		! Real missing data
  REAL :: dxy

  ! Dynamic arrays:
  REAL,ALLOCATABLE,DIMENSION(:) :: pressr1d,height1d
  REAL,ALLOCATABLE,DIMENSION(:) :: rlat_radar,rlon_radar
  REAL,ALLOCATABLE,DIMENSION(:) :: rhgt_radar,n_grid_vel
  REAL,ALLOCATABLE,DIMENSION(:,:) :: lat,lon,topo
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: height3d,pressr3d
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: temptr3d,sphumd3d
  REAL,ALLOCATABLE,DIMENSION(:,:,:) :: u_wind3d,v_wind3d
  REAL,ALLOCATABLE,DIMENSION(:,:,:,:) :: grid_radar_vel

  REAL :: obs_temp(maxxobs,12)
  type (barnesob) :: obs_point(maxxobs)

END MODULE LAPS_Parm
