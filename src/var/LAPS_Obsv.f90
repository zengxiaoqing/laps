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
!  This routine retrieves observation through LAPS ingest.
!
!  HISTORY: 
! 	Creation: YUANFU XIE	3-2006	Wind obs;
!	Modified: YUANFU XIE	5-2006	add temp obs;
!==========================================================

  USE LAPS_Parm
  USE MEM_NAMELIST

  IMPLICIT NONE

  ! Local variables:
  INTEGER :: useable_radar,istatus_remap_pro,status
  INTEGER :: INIT_TIMER,ntmin,ntmax
  REAL :: u(n(1),n(2),n(3)),v(n(1),n(2),n(3))
  REAL :: ulaps(n(1),n(2),n(3)),vlaps(n(1),n(2),n(3))
  REAL :: wt(n(1),n(2),n(3)),weight_prof

  INTEGER :: max_snd_grd,max_snd_lvl
  LOGICAL :: l_adj_hgt
  REAL :: bg_weight

  ! Profiler Data:
  ntmin = -1
  ntmax = +1

  status = INIT_TIMER()

  CALL get_wind_3d_obs(n(1),n(2),n(3),rmissing,imissing, &
	i4time,height3d,height1d,max_pr,max_pr_lvls,     &
	weight_prof,l_raob,l_cdw,n_sao,n_pirep,lat,lon,  &
	ntmin,ntmax,u,v,ulaps,vlaps,wt,maxxobs,obs_point, &
	nobs_point,rlat_radar,rlon_radar,rhgt_radar,     &
	useable_radar,n_grid_vel,istatus_remap_pro,status)

  ! Temperature parameters:
  ! CALL get_temp_parms(l_raob,l_use_raob,l_adj_hgt,       &
  !                    bg_weight,rms_thresh,              &
  !                    pres_mix_thresh,max_snd_grd,max_snd_lvl,maxtobs,   &
  !                    status)

  ! Assign temperature parameters using the new temperature module:
  l_raob = l_read_raob_t	
  maxtobs = max_obs
  CALL get_temp_3d_obs(max_snd_grid,max_snd_levels,max_obs)

END SUBROUTINE LAPS_Obsv

SUBROUTINE get_temp_3d_obs(max_snd_grd,max_snd_lvl,max_obs)

!==========================================================
!  This routine retrieves 3D temperature observation data.
!
!  HISTORY:
!	Creation: YUANFU XIE	5-2006.
!==========================================================

  USE LAPS_Parm

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: max_snd_grd,max_snd_lvl,max_obs

  INTEGER :: i_obstype
  INTEGER :: status

  include 'tempobs.inc'

  ! Sonde temperature data:
  CALL get_temp_snd(max_snd_grd,max_snd_lvl,temp_obs,max_obs,status)

  ! ACAR temperature data:
  CALL get_temp_acar(temp_obs,max_obs,status)

  ! Copy to the output array:
  obs_temp(1:max_obs,1:12) = temp_obs(1:max_obs,1:12)

END SUBROUTINE get_temp_3d_obs

SUBROUTINE get_temp_snd(max_snd_grd,max_snd_lvl,temp_obs,maxaobs,error)

!==========================================================
!  This routine retrieves temperature obs from sonde data.
!
!  HISTORY:
!	Adapted from insert_tobs: YUANFU XIE	5-2006.
!==========================================================

  USE LAPS_Parm

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: max_snd_grd,max_snd_lvl,maxaobs
  INTEGER, INTENT(OUT) :: error
  REAL, INTENT(OUT) :: temp_obs(maxaobs,12)	! tempobs.inc
  
  INTEGER :: i4_window_raob
  INTEGER :: status,n_rass,n_snde,n_tsnd
  REAL :: lattsnd(max_snd_grd),lontsnd(max_snd_grd)
  REAL :: tsnd(max_snd_grd,n(3))
  REAL :: inst_err_tsnd(max_snd_grd)
  REAL :: bias_tsnd(max_snd_grd,n(3)),bias_htlow(max_snd_grd)
  CHARACTER*5 c5name(max_snd_grd)
  CHARACTER*8 c8obstype(max_snd_grd)
  LOGICAL :: l_struct

  INTEGER :: k,isnd,n_qc_snd
  INTEGER :: igrid(max_snd_grd),jgrid(max_snd_grd)
  REAL :: ri,rj,p_pa,sh,tvir,tamb,devirt_sh
  LOGICAL :: l_string_contains,l_qc

  i4_window_raob = 0 		! Same as insert_tobs
  l_struct = .true.

  ! Read sonde temperature:
   CALL read_tsnd(i4time,height3d,temptr3d,sphumd3d, &
		 pressr3d,lattsnd,lontsnd,lat,lon, &
                 max_snd_grd,max_snd_lvl,tsnd,inst_err_tsnd,c5name, &
                 c8obstype,l_raob,l_struct, &
                 i4_window_raob,bias_htlow, &
                 n_rass,n_snde,n_tsnd,timelen,n(1), &
                 n(2),n(3),rmissing,status)

  ! QC temperature observation data:
  n_qc_snd = 0
  n_tobs = 0
  DO isnd=1,n_tsnd

    ! Location on the grid:
    CALL latlon_to_rlapsgrid(lattsnd(isnd),lontsnd(isnd), &
                             lat,lon,n(1),n(2),ri,rj,status)

    ! Count obs within the domain:
    IF (status .EQ. 1) THEN
      igrid(isnd) = nint(ri)
      jgrid(isnd) = nint(rj)

      ! Find the sonde data:
      DO k=1,n(3)
        IF (tsnd(isnd,k) .NE. rmissing) THEN
	  ! RASS observes the virtue temperature:
          IF(l_string_contains(c8obstype(isnd),'RASS',status)) THEN       
            ! Convert from virtual temperature to temperature
            tvir = tsnd(isnd,k)
            sh = sphumd3d(igrid(isnd),jgrid(isnd),k)       
            p_pa = pressr3d(igrid(isnd),jgrid(isnd),k)    
            tamb = devirt_sh(tvir,sh,p_pa)
          ELSE
            sh = 0.
            tamb = tsnd(isnd,k)
          END IF
	  ! Save bias:
          bias_tsnd(isnd,k) =  tamb - &
            temptr3d(igrid(isnd),jgrid(isnd),k)

	  ! HARD QC: 
          IF (ABS(bias_tsnd(isnd,k)) .GT. 10.) THEN
            l_qc = .true.
            WRITE(6,*) ' ABS(Temp - FIRST GUESS) > 10., Temp NOT USED'       
          ENDIF
        ENDIF
      ENDDO

      ! Good sonde obs:
      IF (.NOT. l_qc) THEN
        n_qc_snd = n_qc_snd+1
        ! Count the observations:
        DO k=1,n(3)
          IF (bias_tsnd(isnd,k) .NE. rmissing) THEN
            n_tobs = n_tobs+1

            IF (n_tobs .GT. maxaobs) THEN
              WRITE(6,*) 'Too many temperature sonde obs'
              error = 0
              RETURN
            ENDIF

            ! Pass the obs to temperature obs array:
	    ! Note: Referring to tempobs.inc under include
            temp_obs(n_tobs,1) = igrid(isnd)
            temp_obs(n_tobs,2) = jgrid(isnd)
            temp_obs(n_tobs,3) = k
            temp_obs(n_tobs,4) = igrid(isnd)
            temp_obs(n_tobs,5) = jgrid(isnd)
            temp_obs(n_tobs,6) = k
            temp_obs(n_tobs,8) = &
              temptr3d(igrid(isnd),jgrid(isnd),k) &    
              + bias_tsnd(isnd,k)
            temp_obs(n_tobs,9) = bias_tsnd(isnd,k)
            temp_obs(n_tobs,10) = &
              1.0 / inst_err_tsnd(isnd)**2
            temp_obs(n_tobs,11) = inst_err_tsnd(isnd)
          ENDIF
        ENDDO
      ENDIF
    ENDIF

  ENDDO

  WRITE(6,*) '# of sonde temperature obs: ',n_tobs,n_qc_snd

END SUBROUTINE get_temp_snd

SUBROUTINE get_temp_acar(temp_obs,maxaobs,status)

!==========================================================
!  This routine reads in the temperature obs from ACAR.
!
!  HISTORY:
!	Adapted from insert_tobs: YUANFU XIE	5-2006.
!==========================================================

  USE LAPS_Parm
  USE MEM_NAMELIST

  IMPLICIT NONE

  INTEGER, INTENT(IN) :: maxaobs
  INTEGER, INTENT(OUT) :: status
  REAL, INTENT(INOUT) :: temp_obs(maxaobs,12)

  INTEGER :: n_good_acars

  CALL rd_acars_t(i4time,height3d,temptr3d &
                 ,n_pirep,n_good_acars,'pin' &
                 ,n(1),n(2),n(3),lat,lon,rmissing &
                 ,temp_obs,maxaobs,n_tobs,status)

END SUBROUTINE get_temp_acar


