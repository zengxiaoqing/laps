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

SUBROUTINE LAPS_Conf

!==========================================================
!  This routine initializes LAPS configuration parameters:
!  i4time, dimensions and so on.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006
!========================================================== 

  USE LAPS_Parm
  USE MEM_NAMELIST

  IMPLICIT NONE

  ! Local variables:
  CHARACTER :: a9time*9,fnm*9,hr*2,mins*2,jday*5
  INTEGER :: status,i4time_sys
!  INTEGER ::thresh_2_radarobs_lvl_unfltrd, &
!              thresh_4_radarobs_lvl_unfltrd, &
!              thresh_9_radarobs_lvl_unfltrd, i4time_sys
!  REAL :: weight_bkg_const_wind,weight_radar,rms_thresh_wind
!  INTEGER :: max_obs

  ! Wind parameters:
  CHARACTER*150 :: static_dir, filename
  INTEGER       :: len_dir

!  NAMELIST /wind_nl/ l_raob, l_cdw, l_radial, &
!                     thresh_2_radarobs_lvl_unfltrd, &
!                     thresh_4_radarobs_lvl_unfltrd, &
!                     thresh_9_radarobs_lvl_unfltrd, &
!                     weight_bkg_const_wind, &
!                     weight_radar, &
!                     rms_thresh_wind, &
!                     max_pr,max_pr_lvls,max_obs

  ! Spatial dimensions:
  CALL get_grid_dim_xy(n(1),n(2),status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error in horizontal dimensions'
    STOP
  ENDIF
  CALL get_laps_dimensions(n(3),status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error in vertical dimension'
    STOP
  ENDIF
PRINT*,'N = ',n

  ! System time:
  CALL get_systime(i4time,a9time,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error in system times'
    STOP
  ELSE
    CALL get_directory('log',filename,len_dir)
    filename = filename(1:len_dir)//'i4time.txt'
    ! OPEN(10,file='i4time.txt')
    OPEN(10,file=filename(1:len_dir+10))
    WRITE(10,*) i4time
    CLOSE(10)
  ENDIF
  CALL get_systime_all(i4time_sys,fnm,hr,mins,asctime,jday,status)
  IF (i4time .NE. i4time_sys) THEN
    PRINT*,'LAPS_Conf: error in reading background at wrong time'
    STOP
  ENDIF
  CALL get_laps_cycle_time(timelen,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error in LAPS cycle time'
    STOP
  ENDIF
PRINT*,'Time: ',i4time,a9time,timelen

  ! Missing data:
  CALL get_r_missing_data(rmissing,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error obtaining real missing_data'
    STOP
  ENDIF
  CALL get_i2_missing_data(imissing,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error obtaining integer missing_data'
    STOP
  ENDIF

  ! Get wind parameters:
  CALL GET_DIRECTORY('static',static_dir,len_dir)
  filename = static_dir(1:len_dir)//'wind.nl'
  CALL READ_NAMELIST_LAPS ('wind',filename)
  ! Get temperature parameters:
  filename = static_dir(1:len_dir)//'temp.nl'
  CALL READ_NAMELIST_LAPS ('temp_anal',filename)

  ! Radar info:
  CALL get_max_radars(max_radars,status)
    IF (status .NE. 1) THEN
      WRITE(6,*) 'Error obtaining max_radars'
      STOP
    ENDIF

  ! Meso, SAO and PIREP:
  CALL get_meso_sao_pirep(N_MESO,N_SAO,N_PIREP,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'Error obtaining N_MESO, N_SAO, and N_PIREP'
    STOP
  ENDIF

  ! Allocate LAPS dynamic arrays:
  CALL LAPS_Allc

  ! Read lat/lon and topo data:
  CALL read_static_grid(n(1),n(2),'LAT',lat,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error get LAPS LAT'
    STOP
  ENDIF
  CALL read_static_grid(n(1),n(2),'LON',lon,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error get LAPS LON'
    STOP
  ENDIF
  CALL read_static_grid(n(1),n(2),'AVG',topo,status)
  IF (status .NE. 1) THEN
    WRITE(6,*) 'LAPS_Conf: error get LAPS topography'
    STOP
  ENDIF
PRINT*,'LAT/LON/AVG: ',lat(1,1),lon(1,1),topo(1,1)

  ! Grid spacing:
  CALL get_grid_spacing_actual(lat(n(1)/2+1,n(2)/2+1), &
		               lon(n(1)/2+1,n(2)/2+1), &
			       dxy,status)
PRINT*,'Grid spacing: ',dxy

END SUBROUTINE LAPS_Conf
