!dis
!dis    Open Source License/Disclaimer, Forecast Systems Laboratory
!dis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
!dis
!dis    This software is distributed under the Open Source Definition,
!dis    which may be found at http://www.opensource.org/osd.html.
!dis
!dis    In particular, redistribution and use in source and binary forms,
!dis    with or without modification, are permitted provided that the
!dis    following conditions are met:
!dis
!dis    - Redistributions of source code must retain this notice, this
!dis    list of conditions and the following disclaimer.
!dis
!dis    - Redistributions in binary form must provide access to this
!dis    notice, this list of conditions and the following disclaimer, and
!dis    the underlying source code.
!dis
!dis    - All modifications to this software must be clearly documented,
!dis    and are solely the responsibility of the agent making the
!dis    modifications.
!dis
!dis    - If significant modifications or enhancements are made to this
!dis    software, the FSL Software Policy Manager
!dis    (softwaremgr@fsl.noaa.gov) should be notified.
!dis
!dis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
!dis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
!dis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
!dis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
!dis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
!dis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
!dis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
!dis
!dis                                                                   

PROGRAM wfoprep

  ! This is a new version of the wfoprep program.  Unlike the previous
  ! version of this programs distributed with LAPS, this version is
  ! focused on acquiring data for lateral boundary conditions to support
  ! a mesoscale NWP model run.  For initialization, the new program
  ! "lapsprep" builds the files necessary for the initial condition, assuming
  ! you wish to initialize with LAPS.  Furthermore, the old version was
  ! set up to support only the SFM, which we are no longer supporting.  This
  ! version can output files for ingest by MM5v3/REGRIDDER, RAMS 4.3 in RALPH2
  ! format, or WRFSI/hinterp input.   This is controlled by the output_type
  ! entry(ies) in the wfoprep.nl.  Finally, this new version is cast in
  ! free-form F90 code, including the use of modules.
  !
  ! History
  ! -------
  ! Sep 2001:  New f90 version for MM5v3/RAMS4.3/WRFSI support.  
  !            B. Shaw, NOAA/CIRA/FSL FRD/LAPB
  !
  !
  USE wfoprep_setup 
  USE wfo_models
  USE map_utils
  USE wfoprep_mm5
  USE wfoprep_wrf
 
  IMPLICIT NONE  
  REAL, PARAMETER           :: rmissingval = -9999.
  INTEGER,PARAMETER         :: maxl=75 
  INTEGER,PARAMETER         :: maxt=50
  LOGICAL                   :: filefound
  INTEGER                   :: fcstsec(maxt)
  INTEGER                   :: goodlevs
  REAL                      :: goodpct
  INTEGER                   :: i,ii
  INTEGER                   :: i4time_offset
  INTEGER                   :: i4time_last
  INTEGER                   :: i4time_valid,i4time_valid1,i4time_valid2
  INTEGER                   :: i4time_cycle
  INTEGER, ALLOCATABLE      :: i4times_avail_max(:)
  INTEGER, ALLOCATABLE      :: i4times_avail(:)
  INTEGER                   :: istatus
  INTEGER                   :: k,kk
  CHARACTER (LEN=13)        :: last_cycle_processed
  CHARACTER (LEN=13)        :: latest_cycle_wfo

  INTEGER                   :: m
  CHARACTER (LEN=256)       :: modelfile_wfo
  INTEGER                   :: nfid
  INTEGER                   :: ntimes,ntimes_needed
  INTEGER                   :: outfreq_sec
  CHARACTER (LEN=256)       :: proclog
  INTEGER                   :: ta
  INTEGER                   :: t_id,ht_id,u_id,v_id,rh_id,msl_id
  INTEGER                   :: nz_t,nz_ht,nz_u,nz_v,nz_rh,nz_msl
  INTEGER                   :: np_t,np_ht,np_u,np_v,np_rh,np_rh_raw,np_msl
  LOGICAL                   :: havesfc_t,havesfc_ht,havesfc_u,havesfc_msl
  LOGICAL                   :: havesfc_v,havesfc_rh
  CHARACTER (LEN=10)        :: t_levels_c(maxl)
  CHARACTER (LEN=10)        :: ht_levels_c(maxl)
  CHARACTER (LEN=10)        :: rh_levels_c(maxl)
  CHARACTER (LEN=10)        :: u_levels_c(maxl)
  CHARACTER (LEN=10)        :: v_levels_c(maxl)
  CHARACTER (LEN=10)        :: msl_levels_c(maxl)
  REAL                      :: t_plevels(maxl)
  REAL                      :: ht_plevels(maxl)
  REAL                      :: rh_plevels(maxl),rh_plevels_raw(maxl)
  REAL                      :: u_plevels(maxl)
  REAL                      :: v_plevels(maxl)
  REAL                      :: msl_plevels(maxl)
  INTEGER                   :: t_kbotp,t_ktopp, t_ksfc
  INTEGER                   :: rh_kbotp,rh_ktopp,rh_ksfc
  INTEGER                   :: ht_kbotp,ht_ktopp,ht_ksfc
  INTEGER                   :: u_kbotp,u_ktopp,u_ksfc
  INTEGER                   :: v_kbotp,v_ktopp,v_ksfc
  INTEGER                   :: msl_kbotp,msl_ktopp,msl_ksfc
  CHARACTER (LEN=13)        :: wfofname
  CHARACTER (LEN=13),EXTERNAL :: cvt_i4time_wfo_fname13
  LOGICAL, ALLOCATABLE      :: t_inv(:,:)
  LOGICAL, ALLOCATABLE      :: u_inv(:,:)
  LOGICAL, ALLOCATABLE      :: v_inv(:,:)
  LOGICAL, ALLOCATABLE      :: rh_inv(:,:)
  LOGICAL, ALLOCATABLE      :: rh_inv_raw(:,:)
  LOGICAL, ALLOCATABLE      :: ht_inv(:,:)
  LOGICAL, ALLOCATABLE      :: msl_inv(:,:)
  LOGICAL, ALLOCATABLE      :: goodtime_flag (:)
  INTEGER                   :: n_goodtimes
  CHARACTER(LEN=10)         :: mslname
  INTEGER, ALLOCATABLE      :: time_index(:) 
  INTEGER                   :: nfssttimes
  TYPE(proj_info)           :: proj
  REAL                      :: weight1
  ! Data arrays
  REAL, ALLOCATABLE         :: z3d(:,:,:),z3d1(:,:,:),z3d2(:,:,:)
  REAL, ALLOCATABLE         :: zsf(:,:),zsf1(:,:),zsf2(:,:)
  REAL, ALLOCATABLE         :: t3d(:,:,:),t3d1(:,:,:),t3d2(:,:,:)
  REAL, ALLOCATABLE         :: tsf(:,:),tsf1(:,:),tsf2(:,:)
  REAL, ALLOCATABLE         :: rh3d(:,:,:),rh3d1(:,:,:),rh3d2(:,:,:)
  REAL, ALLOCATABLE         :: rhsf(:,:),rhsf1(:,:),rhsf2(:,:) 
  REAL, ALLOCATABLE         :: u3d(:,:,:),u3d1(:,:,:),u3d2(:,:,:)
  REAL, ALLOCATABLE         :: usf(:,:),usf1(:,:),usf2(:,:) 
  REAL, ALLOCATABLE         :: v3d(:,:,:),v3d1(:,:,:),v3d2(:,:,:)
  REAL, ALLOCATABLE         :: vsf(:,:),vsf1(:,:),vsf2(:,:)
  REAL, ALLOCATABLE         :: slp(:,:),slp1(:,:),slp2(:,:)
  REAL, ALLOCATABLE         :: fsstsum(:,:), fsst(:,:)
  REAL, ALLOCATABLE         :: topo(:,:)
  REAL, ALLOCATABLE         :: data3d_temp(:,:,:)
  REAL                      :: tdsf_c, psf_mb
  REAL, EXTERNAL            :: dwpt_laps, twet_fast
  INTEGER                   :: ix,jy
  LOGICAL                   :: fixvar1, fixvar2

  istatus = 1
  PRINT '(A)', 'Calling setup routine...'
  CALL setup_wfoprep(istatus)
  IF (istatus .ne. 1) THEN
    PRINT '(A)', 'DPREP: Problem with setup configuration.'
    STOP
  ENDIF

  ! At this point, the wfoprep_setup module contains a list of models
  ! and information relating to their availability (e.g., frequency
  ! of model runs, delay time before a cycle is available, etc.).  We
  ! will loop through each valid model and process as requested and
  ! available.

  model_loop: DO m = 1, num_models

    outfreq_sec = output_freq(m) * 3600
    PRINT '(2A)', '*** Processing model ', TRIM(model_name(m))
    ! Based on the namelist entries that tell us how frequently this
    ! model is run and how many hours after cycle time it becomes available,
    ! we can do some computations to see which cycle we should be trying to
    ! process.

    ! First, compute the i4time_offset value, which is the current time
    ! minus the amount of time between this model's cycle time and when
    ! the entire run is available.

    i4time_offset = i4time_now - model_delay(m)*3600
    
    ! Round i4time_offset to nearest hour
    i4time_offset = i4time_offset - MOD(i4time_offset,3600)

    ! Now, compute the closest cycle time for this model that is equal 
    ! to or earlier than the i4time_offset time.

    i4time_cycle = (i4time_offset / (model_run_freq(m)*3600)) * &
                       model_run_freq(m)*3600    
    i4time_last = i4time_cycle + max_fcst_len(m)*3600

    ! What is the cycle time closest?
    latest_cycle_wfo = cvt_i4time_wfo_fname13(i4time_cycle)
    print *, '    Latest possible cycle for this model = ', latest_cycle_wfo

    ! Check the {model}.last file in the drep output directory
    ! to see if this cycle has been processed yet.
 
    PRINT *, '    Checking log of previously processed cycles.' 
    proclog = TRIM(ext_data_path) // TRIM(output_name(m)) // '.last'
    INQUIRE (FILE=proclog, EXIST=filefound)
    IF (filefound) THEN
      OPEN(FILE=proclog,UNIT=11,STATUS='OLD',FORM='FORMATTED', &
           ACCESS='SEQUENTIAL')
      READ(11,'(A13)') last_cycle_processed
      CLOSE(11)
      PRINT '(2A)', '       Last cycle processed for this model:', &
           last_cycle_processed
      IF (last_cycle_processed .GE. latest_cycle_wfo) THEN
        PRINT '(A)', '       Current cycle <= last cycle processed. Skipping...'
        CYCLE model_loop
      ENDIF
    ELSE  
      PRINT '(A)', '    No record of processing for any previous cycles.'
    ENDIF
  
    ! See if the file required for this cycle is available
    modelfile_wfo = TRIM(model_path(m)) // '/' // latest_cycle_wfo

    PRINT '(2A)', '    Searching for data file: ',TRIM(modelfile_wfo)
    INQUIRE(FILE=modelfile_wfo, EXIST=filefound)
    IF (.NOT. filefound) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#############################'
      PRINT '(A)', 'SKIPPING MODEL: NOT AVAILABLE'
      PRINT '(A)', '#############################' 
      PRINT '(A)', ' ' 
      CYCLE model_loop
    ENDIF

    ! Open the file
    CALL open_wfofile(modelfile_wfo,nfid,istatus)
    IF (istatus .NE. 1)THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#############################'
      PRINT '(A)', 'SKIPPING MODEL: I/O ERROR'
      PRINT '(A)', '#############################'
      PRINT '(A)', ' '
      CYCLE model_loop                          
    ENDIF
   
    ! Get grid/projection info
    CALL get_wfomodel_proj(nfid,model_name(m),proj,istatus)
    IF (istatus .NE. 1) THEN
            PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: PROJ INFO PROBLEM'
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '     
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF

    ! Process time info
    CALL get_wfomodel_fcsttimes(nfid,ntimes,fcstsec,istatus)
    IF (istatus.NE.1) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: TIME INFO PROBLEM'
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '     
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF
    ALLOCATE(i4times_avail_max(ntimes))
    i4times_avail_max(1:ntimes) = i4time_cycle + fcstsec(1:ntimes)
    IF (i4time_last .GT. i4times_avail_max(ntimes)) THEN
      PRINT *, 'INFO: This source does not support max_fcst length of ',max_fcst_len(m)
      PRINT *, 'INFO: Resetting to ', (i4times_avail_max(ntimes)-i4time_cycle)/36
      i4time_last = i4times_avail_max(ntimes)
    ELSE  
      ! Compute the ntimes_needed value, which reprenents the number
      ! of time periods actually needed to get us equal to or just past
      ! i4time_last
      ntimes_needed = 1
      find_ntimes_needed: DO i = 2, ntimes
        ntimes_needed = ntimes_needed + 1
        IF (i4times_avail_max(i) .GE. i4time_last) THEN
          EXIT find_ntimes_needed
        ENDIF
      ENDDO find_ntimes_needed
        
    ENDIF

    ! Get level information for state variables (t,ht,u,v,rh) and then
    ! get the inventory  

    ! Temperature
    CALL get_wfomodel_var_levels(nfid,'t         ',t_id, nz_t, t_levels_c, np_t, &
                  t_plevels,t_kbotp,t_ktopp, havesfc_t,t_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: NO T DATA        '
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '                  
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF
  
    ! Height is only mandatory for 3D data sets
    IF (model_code(m).GT.1) THEN
      CALL get_wfomodel_var_levels(nfid,'gh        ',ht_id, nz_ht, ht_levels_c, np_ht, &
                 ht_plevels,ht_kbotp,ht_ktopp, havesfc_ht, ht_ksfc,istatus)
      IF (istatus.NE.1)THEN
        PRINT '(A)', ' '
        PRINT '(A)', '#################################'
        PRINT '(A)', 'SKIPPING MODEL: NO Z DATA        '
        PRINT '(A)', '#################################'
        PRINT '(A)', ' '  
        CALL close_wfofile(nfid,istatus)
        CYCLE model_loop
      ENDIF                           
    ENDIF 
     
    ! Relative humidity
    CALL get_wfomodel_var_levels(nfid,'rh        ',rh_id, nz_rh, rh_levels_c, np_rh_raw, &
                  rh_plevels_raw,rh_kbotp,rh_ktopp, havesfc_rh,rh_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: NO RH DATA        '
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '                    
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF                      

    ! U-wind component
    CALL get_wfomodel_var_levels(nfid,'uw        ',u_id, nz_u, u_levels_c, np_u, &
                  u_plevels, u_kbotp,u_ktopp,havesfc_u,u_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: NO U DATA        '
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '        
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF                  

    ! V-wind component
    CALL get_wfomodel_var_levels(nfid,'vw        ',v_id, nz_v, v_levels_c, np_v, &
                  v_plevels,v_kbotp,v_ktopp, havesfc_v,v_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT '(A)', ' '
      PRINT '(A)', '#################################'
      PRINT '(A)', 'SKIPPING MODEL: NO V DATA        '
      PRINT '(A)', '#################################'
      PRINT '(A)', ' '   
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF     

    ! MSLP is mandatory for 3D data sets.  Check for Eta MSLP first.  If non-existent,
    ! then check for standard mslp
  
    IF (model_code(m).GT.1) THEN 
      CALL get_wfomodel_var_levels(nfid,'emsp      ',msl_id, nz_msl, &
                 msl_levels_c, np_msl, &
                 msl_plevels,msl_kbotp,msl_ktopp, havesfc_msl, msl_ksfc,istatus)
      IF (istatus.EQ.1)THEN
        PRINT *, 'Using ETA MSLP field'
        mslname = 'emsp      '
      ELSE
        CALL get_wfomodel_var_levels(nfid,'pmsl      ',msl_id, nz_msl, &
                 msl_levels_c, np_msl, &
                 msl_plevels,msl_kbotp,msl_ktopp, havesfc_msl, msl_ksfc,istatus) 
        IF (istatus .EQ. 1) THEN
          PRINT *, 'Using PMSL field'
          mslname = 'pmsl      '
        ELSE
          CALL get_wfomodel_var_levels(nfid,'mmsp      ',msl_id, nz_msl, &
                 msl_levels_c, np_msl, &
                 msl_plevels,msl_kbotp,msl_ktopp, havesfc_msl, msl_ksfc,istatus)
          IF (istatus .EQ. 1) THEN
            PRINT *, 'Using MMSP field (MAPS Sea-level Pressure)'
            mslname = 'mmsp      '
          ELSE
            PRINT '(A)', ' '
            PRINT '(A)', '#################################'
            PRINT '(A)', 'SKIPPING MODEL: NO MSLP DATA        '
            PRINT '(A)', '#################################'
            PRINT '(A)', ' '  
            CALL close_wfofile(nfid,istatus)
            CYCLE model_loop
          ENDIF
        ENDIF
      ENDIF
    ENDIF    
         
    ! Now, check to make sure we have T,ht, u, and v for all of the same
    ! pressure levels if this is a 3D data set.  This seems to be the case
    ! for all supported datasets.  However, RH is a bit trickier.  In some
    ! cases, RH is not provided above 300 mb or below 975 mb.  Deal with this
    ! later.

    IF ( model_code(m) .GT. 1) THEN
      IF ( (np_t .NE. np_ht) .OR. (np_u .NE. np_ht) .OR. &
           (np_u .NE. np_ht) ) THEN
         PRINT '(A)', ' '
          PRINT '(A)', '####################################'
          PRINT '(A)', 'SKIPPING MODEL: DIMENSION MISMATCHES'
          PRINT '(A)', '####################################'
          PRINT '(A)', ' '              
         DEALLOCATE(i4times_avail_max)
         CALL close_wfofile(nfid,istatus)
         CYCLE model_loop
      ENDIF
    ENDIF

    ! We have made it this far, so we must have variables we need, but we
    ! need to check their inventory to see which time periods are complete.
    ! If this is a 3D data source (model code > 1), then we need to have at
    ! least the top and bottom level values for each variable for the time
    ! to be included in the process.  Also must have the first and last
    ! time periods.  
 
    n_goodtimes = 0 
    IF (ALLOCATED(goodtime_flag))DEALLOCATE(goodtime_flag)
    ALLOCATE(goodtime_flag(ntimes_needed))
    goodtime_flag(:) = .true.

    IF (ALLOCATED(t_inv)) DEALLOCATE(t_inv)
    ALLOCATE(t_inv(nz_t,ntimes))
    CALL get_wfomodel_var_inv(nfid,'t         ',nz_t,ntimes,t_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '####################################'
      PRINT '(A)', 'SKIPPING MODEL: NO T INVENTORY      '
      PRINT '(A)', '####################################'
      PRINT '(A)', ' ' 
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF

    IF (ALLOCATED(u_inv)) DEALLOCATE(u_inv)
    ALLOCATE(u_inv(nz_u,ntimes))
    CALL get_wfomodel_var_inv(nfid,'uw        ',nz_u,ntimes,u_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '####################################'
      PRINT '(A)', 'SKIPPING MODEL: NO U INVENTORY      '
      PRINT '(A)', '####################################'
      PRINT '(A)', ' '    
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF      

    IF (ALLOCATED(v_inv)) DEALLOCATE(v_inv)
    ALLOCATE(v_inv(nz_v,ntimes))
    CALL get_wfomodel_var_inv(nfid,'vw        ',nz_v,ntimes,v_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '####################################'
      PRINT '(A)', 'SKIPPING MODEL: NO V INVENTORY      '
      PRINT '(A)', '####################################'
      PRINT '(A)', ' '    
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF       

    IF (ALLOCATED(rh_inv_raw)) DEALLOCATE(rh_inv_raw)
    ALLOCATE(rh_inv_raw(nz_rh,ntimes))
    ! Also allocate an array that matches the temperature array
    IF (ALLOCATED(rh_inv)) DEALLOCATE(rh_inv)
    ALLOCATE(rh_inv(nz_t,ntimes))
    rh_inv(:,:) = .false.
    IF (nz_rh .EQ. nz_t) rh_inv = rh_inv_raw
    CALL get_wfomodel_var_inv(nfid,'rh        ',nz_rh,ntimes,rh_inv_raw,istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', ' '
      PRINT '(A)', '####################################'
      PRINT '(A)', 'SKIPPING MODEL: NO RH INVENTORY      '
      PRINT '(A)', '####################################'
      PRINT '(A)', ' '    
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF   

    IF (ALLOCATED(ht_inv)) DEALLOCATE(ht_inv)
    IF (model_code(m).GT.1) THEN
      ALLOCATE (ht_inv(nz_ht,ntimes))
      CALL get_wfomodel_var_inv(nfid,'gh        ',nz_ht,ntimes,ht_inv,istatus)
      IF (istatus .NE. 1) THEN
        PRINT '(A)', ' '
        PRINT '(A)', '####################################'
        PRINT '(A)', 'SKIPPING MODEL: NO Z INVENTORY      '
        PRINT '(A)', '####################################'
        PRINT '(A)', ' '    
        CALL close_wfofile(nfid,istatus)
        CYCLE model_loop
      ENDIF   
    ENDIF  
 
    IF (ALLOCATED(msl_inv)) DEALLOCATE(msl_inv)
    IF (model_code(m).GT.1) THEN
      ALLOCATE (msl_inv(nz_msl,ntimes))
      CALL get_wfomodel_var_inv(nfid,mslname,nz_msl,ntimes,msl_inv,istatus)
      IF (istatus .NE. 1) THEN
        PRINT '(A)', ' '
        PRINT '(A)', '####################################'
        PRINT '(A)', 'SKIPPING MODEL: NO MSLP INVENTORY      '
        PRINT '(A)', '####################################'
        PRINT '(A)', ' '   
        CALL close_wfofile(nfid,istatus)
        CYCLE model_loop
      ENDIF
    ENDIF                                                                

    ! Here is the actual check for minimum requirements.
    ! Minimum variable requirements:
    !  ht, T, U, V, RH on pressure levels + MSLP
    ! For now, require that 100% of the levels be present for each
    ! variable, and that both the top and bottom level must be
    ! present.  In the future, we can add vertical interpolatoin
    ! to allow for missing levels, thus allowing us to reduce the
    ! 100% level threshold

    invloop: DO i = 1,ntimes_needed
 
      IF (model_code(m) .GT. 1) THEN

        ! Check MSL pressure
        IF (.NOT.msl_inv(msl_ksfc,i)) THEN
          PRINT *, 'WARNING: Missing MSL for this time:',i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF

        ! Check height field, requiring most levels to be present!
        goodlevs = 0
        htinvloop: DO k = ht_kbotp, ht_ktopp
          IF (ht_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  htinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_ht)
        IF ((goodpct .LT. min_vert_frac).OR.(.NOT.ht_inv(ht_kbotp,i)).OR.&
            (.NOT.ht_inv(ht_ktopp,i))) THEN
          PRINT *, 'WARNING: Height inventory failed vertical check:', &
              goodpct, i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF

        ! Check temperatures
        goodlevs = 0
        tinvloop: DO k = t_kbotp, t_ktopp
          IF (t_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  tinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_t)
        IF ((goodpct .LT. min_vert_frac).OR.(.NOT.t_inv(t_kbotp,i)).OR.&
            (.NOT.t_inv(t_ktopp,i))) THEN
          PRINT *, 'WARNING: Temperature inventory failed vertical check:', &
              goodpct, i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF               

        ! Check RH
        goodlevs = 0
        rhinvloop: DO k = rh_kbotp, rh_ktopp
          IF (rh_inv_raw(k,i)) goodlevs = goodlevs + 1
        ENDDO  rhinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_rh_raw)
        IF ((goodpct .LT. min_vert_frac).OR.(.NOT.rh_inv_raw(rh_kbotp,i)).OR.&
            (.NOT.rh_inv_raw(rh_ktopp,i))) THEN
          PRINT *, 'WARNING: RH inventory failed vertical check:', &
              goodpct, i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF                   

        ! Check U 
        goodlevs = 0
        uinvloop: DO k = u_kbotp, u_ktopp
          IF (u_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  uinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_u)
        IF ((goodpct .LT. min_vert_frac).OR.(.NOT.u_inv(u_kbotp,i)).OR.&
            (.NOT.u_inv(u_ktopp,i))) THEN
          PRINT *, 'WARNING: U inventory failed vertical check:', &
              goodpct, i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF                

        ! Check V
        goodlevs = 0
        vinvloop: DO k = v_kbotp, v_ktopp
          IF (v_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  vinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_v)
        IF ((goodpct .LT. min_vert_frac).OR.(.NOT.v_inv(v_kbotp,i)).OR.&
            (.NOT.v_inv(v_ktopp,i))) THEN
          PRINT *, 'WARNING: V inventory failed vertical check:', &
              goodpct,i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF 

      ENDIF
      IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3))THEN
        
        ! We need surface values

        IF (.NOT. t_inv(t_ksfc,i)) THEN
          PRINT *, 'WARNING: Missing surface T', i
          goodtime_flag(i) = .false.
          CYCLE invloop 
        ENDIF

        IF (.NOT. u_inv(u_ksfc,i)) THEN
          PRINT *, 'WARNING: Missing surface U', i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF    

        IF (.NOT. v_inv(u_ksfc,i)) THEN
          PRINT *, 'WARNING: Missing surface V',i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF     

        IF (.NOT. rh_inv_raw(rh_ksfc,i)) THEN
          PRINT *, 'WARNING: Missing surface RH',i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF     

      ENDIF  
      ! If we made it this far without cycling the loop, then
      ! this is a good level
      n_goodtimes = n_goodtimes + 1
      goodtime_flag(i) = .true.
    ENDDO invloop 
           
    ! OK, now lets make sure that enough time periods
    ! passed the above check and that the first and last
    ! times are available
    goodlevs = 0
    DO i = 1,ntimes_needed
      IF (goodtime_flag(i)) goodlevs = goodlevs + 1
    ENDDO
    goodpct = FLOAT(goodlevs)/FLOAT(ntimes_needed)
    IF ( (goodpct .LT. min_time_frac).OR.(.NOT.goodtime_flag(1)).OR.&
         (.NOT.goodtime_flag(ntimes_needed)))THEN
      PRINT *, ' '
      PRINT *, '###################################################'
      PRINT *, 'SKIPPING MODEL:  TIME INVENTORY CHECK FAILED'
      PRINT *, 'goodpct/min_time_frac = ', goodpct,min_time_frac
      PRINT *, 'goodtime_Flag(1) = ', goodtime_Flag(1)
      PRINT *, 'ntimes_needed = ', ntimes_needed
      PRINT *, 'goodtime_flag(ntimes_needed) =',goodtime_flag(ntimes_needed)
      PRINT *, '###################################################'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF
         
    ! Allocate the i4times_avail array to the actual number
    ! of usable time periods and populate it
    IF (ALLOCATED(i4times_avail)) DEALLOCATE(i4times_avail) 
    ALLOCATE(i4times_avail(n_goodtimes))
    IF (ALLOCATED(time_index)) DEALLOCATE(time_index)
    ALLOCATE(time_index(n_goodtimes))
    ii = 1
    DO i = 1, ntimes_needed
      IF (goodtime_flag(i)) THEN
        i4times_avail(ii) = i4times_avail_max(i)
        time_index(ii) = i
        ii = ii + 1
      ENDIF
    ENDDO
    DEALLOCATE(i4times_avail_max)
 

    ! Allocate 3 arrays for each state variable.  Two will be needed for the
    ! data at two times bounding the time of interest.  The third is where
    ! the time interpolated values will be saved. 
   
    IF (model_code(m) .GT. 1) THEN   
      ALLOCATE (z3d     (proj%nx,proj%ny,np_ht) )
      ALLOCATE (z3d1    (proj%nx,proj%ny,np_ht) )
      ALLOCATE (z3d2    (proj%nx,proj%ny,np_ht) )
      ALLOCATE (t3d     (proj%nx,proj%ny,np_t) )
      ALLOCATE (t3d1    (proj%nx,proj%ny,np_t) )
      ALLOCATE (t3d2    (proj%nx,proj%ny,np_t) )
      ALLOCATE (rh3d     (proj%nx,proj%ny,np_t) )
      ALLOCATE (rh3d1    (proj%nx,proj%ny,np_rh_raw) )
      ALLOCATE (rh3d2    (proj%nx,proj%ny,np_rh_raw) ) 
      ALLOCATE (u3d     (proj%nx,proj%ny,np_u) )
      ALLOCATE (u3d1    (proj%nx,proj%ny,np_u) )
      ALLOCATE (u3d2    (proj%nx,proj%ny,np_u) )  
      ALLOCATE (v3d     (proj%nx,proj%ny,np_v) )
      ALLOCATE (v3d1    (proj%nx,proj%ny,np_v) )
      ALLOCATE (v3d2    (proj%nx,proj%ny,np_v) )  
      ALLOCATE (slp     (proj%nx,proj%ny) )
      ALLOCATE (slp1    (proj%nx,proj%ny) )  
      ALLOCATE (slp2    (proj%nx,proj%ny) ) 
    ENDIF 
    IF ( (model_code(m).EQ.1).OR.(model_code(m).EQ.3) )THEN
      ALLOCATE (tsf     (proj%nx,proj%ny) )
      ALLOCATE (tsf1    (proj%nx,proj%ny) )
      ALLOCATE (tsf2    (proj%nx,proj%ny) )    
      ALLOCATE (rhsf    (proj%nx,proj%ny) )
      ALLOCATE (rhsf1   (proj%nx,proj%ny) )
      ALLOCATE (rhsf2   (proj%nx,proj%ny) )    
      ALLOCATE (usf     (proj%nx,proj%ny) )
      ALLOCATE (usf1    (proj%nx,proj%ny) )
      ALLOCATE (usf2    (proj%nx,proj%ny) )    
      ALLOCATE (vsf     (proj%nx,proj%ny) )
      ALLOCATE (vsf1    (proj%nx,proj%ny) )
      ALLOCATE (vsf2    (proj%nx,proj%ny) ) 
      ALLOCATE (fsstsum (proj%nx,proj%ny) )
      ALLOCATE (fsst    (proj%nx,proj%ny) )
      fsstsum = 0.
      nfssttimes = 0
    ENDIF    


    ! Initialize the i4time_valid1 and i4time_valid2, which are the two available
    ! times which bound the desired output time
    ta = 1
    i4time_valid1 = i4times_avail(ta)
    i4time_valid2 = i4times_avail(ta+1)

    ! Read in the data for these two time periods
    IF (model_code(m) .GT. 1) THEN
      ! Get the pressure level data arrays
  
      ! MSLP
      CALL read_wfomodel_data(nfid,msl_id,proj,time_index(ta), &
                              msl_ksfc,msl_ksfc,slp1,istatus)
      CALL read_wfomodel_data(nfid,msl_id,proj,time_index(ta+1), &
                              msl_ksfc,msl_ksfc,slp2,istatus) 
      ! Height
      CALL read_wfomodel_data(nfid,ht_id,proj,time_index(ta), &
                              ht_kbotp,ht_ktopp,z3d1,istatus)
      CALL read_wfomodel_data(nfid,ht_id,proj,time_index(ta+1), &
                              ht_kbotp,ht_ktopp,z3d2,istatus)    
      ! Temperature 
      CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta), &
                              t_kbotp,t_ktopp,t3d1,istatus)
      CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta+1), &
                              t_kbotp,t_ktopp,t3d2,istatus)  

      ! RH
      CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta), &
                              rh_kbotp,rh_ktopp,rh3d1,istatus)
      CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta+1), &
                              rh_kbotp,rh_ktopp,rh3d2,istatus)   
  
      ! U
      CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta), &
                              u_kbotp,u_ktopp,u3d1,istatus)
      CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta+1), &
                              u_kbotp,u_ktopp,u3d2,istatus)  

      ! V
      CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta), &
                              v_kbotp,v_ktopp,v3d1,istatus)
      CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta+1), &
                              v_kbotp,v_ktopp,v3d2,istatus) 
    ENDIF
    IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3))THEN 
      CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta), &
                              t_ksfc,t_ksfc,tsf1,istatus)
      CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta+1), &
                              t_ksfc,t_ksfc,tsf2,istatus) 
      CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta), &
                              u_ksfc,u_ksfc,usf1,istatus)
      CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta+1), &
                              u_ksfc,u_ksfc,usf2,istatus)   
      CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta), &
                              v_ksfc,v_ksfc,vsf1,istatus)
      CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta+1), &
                              v_ksfc,v_ksfc,vsf2,istatus)    
      CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta), &
                              rh_ksfc,rh_ksfc,rhsf1,istatus)
      CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta+1), &
                              rh_ksfc,rh_ksfc,rhsf2,istatus)       

    ENDIF
    ALLOCATE(topo(proj%nx,proj%ny))
    CALL get_wfomodel_topo(nfid,proj,topo,istatus) 

    ! Main loop over all desired output times
    output_time_loop: DO i4time_valid = i4time_cycle,i4time_last,outfreq_sec
     ! PRINT *, 'Getting data for i4time = ', i4time_valid 
      IF (i4time_valid .GT. i4time_valid2) THEN 
        ! We need to get new bounding data times
        read_avail_loop: DO WHILE(i4time_valid .GT. i4time_valid2)
          ta = ta + 1
          i4time_valid1 = i4times_avail(ta)
          i4time_valid2 = i4times_avail(ta+1)
        ENDDO read_avail_loop
        ! Read the data for these 2 times
        IF (model_code(m) .GT. 1) THEN
          ! Get the pressure level data arrays

          ! MSLP
          CALL read_wfomodel_data(nfid,msl_id,proj,time_index(ta), &
                              msl_ksfc,msl_ksfc,slp1,istatus)
          CALL read_wfomodel_data(nfid,msl_id,proj,time_index(ta+1), &
                              msl_ksfc,msl_ksfc,slp2,istatus)
          ! Height
          CALL read_wfomodel_data(nfid,ht_id,proj,time_index(ta), &
                              ht_kbotp,ht_ktopp,z3d1,istatus)
          CALL read_wfomodel_data(nfid,ht_id,proj,time_index(ta+1), &
                              ht_kbotp,ht_ktopp,z3d2,istatus)
          ! Temperature
          CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta), &
                              t_kbotp,t_ktopp,t3d1,istatus)
          CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta+1), &
                              t_kbotp,t_ktopp,t3d2,istatus)

          ! RH
          CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta), &
                              rh_kbotp,rh_ktopp,rh3d1,istatus)
          CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta+1), &
                              rh_kbotp,rh_ktopp,rh3d2,istatus)

          ! U
          CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta), &
                              u_kbotp,u_ktopp,u3d1,istatus)
          CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta+1), &
                              u_kbotp,u_ktopp,u3d2,istatus)

          ! V
          CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta), &
                              v_kbotp,v_ktopp,v3d1,istatus)
          CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta+1), &
                              v_kbotp,v_ktopp,v3d2,istatus)
        ENDIF  
        IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3))THEN
          CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta), &
                              t_ksfc,t_ksfc,tsf1,istatus)
          CALL read_wfomodel_data(nfid,t_id,proj,time_index(ta+1), &
                              t_ksfc,t_ksfc,tsf2,istatus)
          CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta), &
                              u_ksfc,u_ksfc,usf1,istatus)
          CALL read_wfomodel_data(nfid,u_id,proj,time_index(ta+1), &
                              u_ksfc,u_ksfc,usf2,istatus)
          CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta), &
                              v_ksfc,v_ksfc,vsf1,istatus)
          CALL read_wfomodel_data(nfid,v_id,proj,time_index(ta+1), &
                              v_ksfc,v_ksfc,vsf2,istatus)
          CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta), &
                              rh_ksfc,rh_ksfc,rhsf1,istatus)
          CALL read_wfomodel_data(nfid,rh_id,proj,time_index(ta+1), &
                              rh_ksfc,rh_ksfc,rhsf2,istatus)
        ENDIF                                                                                                                              
      ENDIF

      ! Handle case where RH has fewer levels than
      ! the other state variables (only applies if
      ! we are obtaining upper air data (model_code >1)

      ! If we are processing upper air data from this model,
      ! then run through some cleanup to account for missing
      ! levels. 
      IF (model_code(m) .GT. 1) THEN

        ! Handle case where RH has fewer levels than
        ! the other state variables (only applies if
        ! we are obtaining upper air data
        np_rh = np_t
        IF (np_rh_raw .LT. np_t) THEN
          PRINT *, 'WARNING: Expanding RH'
          ALLOCATE(data3d_temp(proj%nx,proj%ny,np_rh)) 
          ! Do the rh3d1 array
          data3d_temp(:,:,:) = rmissingval
          rh1_zloop: DO k = 1,np_rh_raw
            t1_zloop: DO kk = 1, np_rh
              IF  (rh_plevels_raw(k).EQ.t_plevels(kk)) THEN
                rh_inv(kk+t_kbotp-1,time_index(ta))= &
                  rh_inv_raw(k+rh_kbotp-1,time_index(ta))
                IF (rh_inv(kk+t_kbotp-1,time_index(ta))) THEN
                  data3d_temp(:,:,kk) = rh3d1(:,:,k)
                ENDIF
                EXIT t1_zloop
              ENDIF
            ENDDO t1_zloop
          ENDDO rh1_zloop

          DEALLOCATE(rh3d1)
          ALLOCATE(rh3d1(proj%nx,proj%ny,np_rh))
          rh3d1(:,:,:) = data3d_temp(:,:,:)
          
          ! Fill in top level if not already filled
          IF (MAXVAL(rh3d1(:,:,np_rh)) .EQ. rmissingval) THEN
             print *, 'Setting top level RH to 5%'
             rh3d1(:,:,np_rh) = 5.0  ! Very dry
             rh_inv(t_ktopp,time_index(ta)) = .true. 
          ENDIF
          ! Fill in bottom level if not already filled (the 
          ! MesoEta does not have the 1000mb value!)
          IF (MAXVAL(rh3d1(:,:,1)) .EQ. rmissingval) THEN
             rh3d1(:,:,1) = rh3d1(:,:,2)
             rh_inv(t_kbotp,time_index(ta)) = .true.
          ENDIF


          ! Repeat for the rh3d2 array
          data3d_temp(:,:,:) = rmissingval
          rh2_zloop: DO k = 1,np_rh_raw
            t2_zloop: DO kk = 1, np_rh
              IF (rh_plevels_raw(k).EQ.t_plevels(kk)) THEN
                rh_inv(kk+t_kbotp-1,time_index(ta+1))= &
                  rh_inv_raw(k+rh_kbotp-1,time_index(ta+1))
                IF (rh_inv(kk+t_kbotp-1,time_index(ta+1))) THEN
                  data3d_temp(:,:,kk) = rh3d2(:,:,k)
                ENDIF          
                EXIT t2_zloop
              ENDIF
            ENDDO t2_zloop
          ENDDO rh2_zloop
          DEALLOCATE(rh3d2)
          ALLOCATE(rh3d2(proj%nx,proj%ny,np_rh))
          rh3d2(:,:,:) = data3d_temp(:,:,:)

          ! Fill in top level if not already filled
          IF (MAXVAL(rh3d2(:,:,np_rh)) .EQ. rmissingval) THEN
             rh3d2(:,:,np_rh) = 5.0  ! Very dry
             rh_inv(t_ktopp,time_index(ta+1)) = .true.
          ENDIF        
          ! Fill in bottom level if not already filled
          IF (MAXVAL(rh3d2(:,:,1)) .EQ. rmissingval) THEN
             rh3d2(:,:,1) = rh3d2(:,:,2)
             rh_inv(t_kbotp,time_index(ta+1)) = .true.
          ENDIF

          DEALLOCATE(data3d_temp)
          rh_plevels(1:np_rh) = t_plevels(1:np_t) 
        ELSE
          rh_plevels = rh_plevels_raw
        ENDIF

        ! Clean up all 3D arrays to fill in any missing levels

        ! Z3d
        fixvar1 = .false.
        fixvar2 = .false.
        DO k = 1,np_ht
          IF (.NOT. ht_inv(k-1+ht_kbotp,time_index(ta))) THEN
            z3d1(:,:,k) = rmissingval
            fixvar1 = .true.
          ENDIF
          IF (.NOT. ht_inv(k-1+ht_kbotp,time_index(ta+1))) THEN
            z3d2(:,:,k) = rmissingval
            fixvar2 = .true.
          ENDIF  
          IF (fixvar1) THEN
            PRINT *, 'WARNING: Filling in missing levels for z3d1'
            ! Call the fill routine to fill in for z3d1
            CALL fill_missing_levs(proj%nx,proj%ny,np_ht,ht_plevels, &
                  z3d1,rmissingval,2)
          ENDIF
          IF (fixvar2) THEN
            ! Call the vinterp routine to fill in for z3d2
            PRINT *, 'WARNING: Filling in missing levels for z3d2'  
            CALL fill_missing_levs(proj%nx,proj%ny,np_ht,ht_plevels, &
                  z3d2,rmissingval,2) 
          ENDIF
        ENDDO

        ! t3d
        fixvar1 = .false.
        fixvar2 = .false.
        DO k = 1,np_t
          IF (.NOT. t_inv(k-1+t_kbotp,time_index(ta))) THEN
            t3d1(:,:,k) = rmissingval
            fixvar1 = .true.
          ENDIF
          IF (.NOT. t_inv(k-1+t_kbotp,time_index(ta+1))) THEN
            t3d2(:,:,k) = rmissingval
            fixvar2 = .true.
          ENDIF
          IF (fixvar1) THEN
            ! Call the vinterp routine to fill in for t3d1
            PRINT *, 'WARNING: Filling in missing levels for t3d1' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_t,t_plevels, &
                  t3d1,rmissingval,2) 
          ENDIF
          IF (fixvar2) THEN
            ! Call the vinterp routine to fill in for t3d2
            PRINT *, 'WARNING: Filling in missing levels for t3d2' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_t,t_plevels, &
                  t3d2,rmissingval,2) 
          ENDIF
        ENDDO 

        ! u3d
        fixvar1 = .false.
        fixvar2 = .false.
        DO k = 1,np_u
          IF (.NOT. u_inv(k-1+u_kbotp,time_index(ta))) THEN
            u3d1(:,:,k) = rmissingval
            fixvar1 = .true.
          ENDIF
          IF (.NOT. u_inv(k-1+u_kbotp,time_index(ta+1))) THEN
            u3d2(:,:,k) = rmissingval
            fixvar2 = .true.
          ENDIF
          IF (fixvar1) THEN
            ! Call the vinterp routine to fill in for u3d1
            PRINT *, 'WARNING: Filling in missing levels for u3d1' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_u,u_plevels, &
                  u3d1,rmissingval,1)
          ENDIF
          IF (fixvar2) THEN
            ! Call the vinterp routine to fill in for u3d2
            PRINT *, 'WARNING: Filling in missing levels for u3d2' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_u,u_plevels, &
                  u3d2,rmissingval,1)  
          ENDIF
        ENDDO                         

        ! v3d
        fixvar1 = .false.
        fixvar2 = .false.
        DO k = 1,np_v
          IF (.NOT. v_inv(k-1+v_kbotp,time_index(ta))) THEN
            v3d1(:,:,k) = rmissingval
            fixvar1 = .true.
          ENDIF
          IF (.NOT. v_inv(k-1+v_kbotp,time_index(ta+1))) THEN
            v3d2(:,:,k) = rmissingval
            fixvar2 = .true.
          ENDIF
          IF (fixvar1) THEN
            ! Call the vinterp routine to fill in for v3d1  
            PRINT *, 'WARNING: Filling in missing levels for v3d1' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_v,v_plevels, &
                  v3d1,rmissingval,1)  
          ENDIF
          IF (fixvar2) THEN
            ! Call the vinterp routine to fill in for v3d2
            PRINT *, 'WARNING: Filling in missing levels for z3d2' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_v,v_plevels, &
                  v3d2,rmissingval,1) 
          ENDIF                   
        ENDDO

        ! rh3d
        fixvar1 = .false.
        fixvar2 = .false.
        DO k = 1,np_rh
          IF (.NOT. rh_inv(k-1+t_kbotp,time_index(ta))) THEN
            rh3d1(:,:,k) = rmissingval
            fixvar1 = .true.
          ENDIF
          IF (.NOT. rh_inv(k-1+t_kbotp,time_index(ta+1))) THEN
            rh3d2(:,:,k) = rmissingval
            fixvar2 = .true.
          ENDIF
          IF (fixvar1) THEN
            ! Call the vinterp routine to fill in for rh3d1
            PRINT *, 'WARNING: Filling in missing levels for rh3d1' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_rh,rh_plevels, &
                  rh3d1,rmissingval,1) 
          ENDIF
          IF (fixvar2) THEN
            ! Call the vinterp routine to fill in for rh3d2
            PRINT *, 'WARNING: Filling in missing levels for rh3d2' 
            CALL fill_missing_levs(proj%nx,proj%ny,np_rh,rh_plevels, &
                  rh3d2,rmissingval,1) 
          ENDIF
        ENDDO           
      ENDIF                  

      ! Time to do time interpolation

      ! PRINT *, 'Using bounding i4times of ', i4time_valid1,i4time_valid2
      ! At this point, we have all of the data we need for two
      ! bounding times.  So interpolate to the desired time,
      ! then output the data
      IF (i4time_valid .EQ. i4time_valid1) THEN
        IF (model_code(m).GT.1) THEN
          z3d = z3d1
          t3d = t3d1
          u3d = u3d1
          v3d = v3d1
          rh3d = rh3d1
          slp = slp1
        ENDIF
        IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3)) THEN
          tsf = tsf1
          rhsf = rhsf1
          usf = usf1
          vsf = vsf1
        ENDIF
  
      ELSE IF(i4time_valid .EQ. i4time_valid2) THEN
        IF (model_code(m).GT.1) THEN
          z3d = z3d2
          t3d = t3d2
          u3d = u3d2
          v3d = v3d2
          rh3d = rh3d2
          slp = slp2
        ENDIF
        IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3)) THEN
          tsf = tsf2
          rhsf = rhsf2
          usf = usf2
          vsf = vsf2
        ENDIF       
        
      ELSE
        ! time interpolate
        weight1 = FLOAT(i4time_valid2-i4time_valid)/FLOAT(i4time_valid2-i4time_valid1)
        IF (model_code(m).GT.1) THEN
          z3d = weight1*z3d1 + (1.-weight1)*z3d2
          t3d = weight1*t3d1 + (1.-weight1)*t3d2  
          u3d = weight1*u3d1 + (1.-weight1)*u3d2  
          v3d = weight1*v3d1 + (1.-weight1)*v3d2  
          rh3d = weight1*rh3d1 + (1.-weight1)*rh3d2  
          slp = weight1*slp1 + (1.-weight1)*slp2  
        ENDIF
        IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3)) THEN
          tsf = weight1*tsf1 + (1.-weight1)*tsf2
          rhsf = weight1*rhsf1 + (1.-weight1)*rhsf2
          usf = weight1*usf1 + (1.-weight1)*usf2
          vsf = weight1*vsf1 + (1.-weight1)*vsf2
        ENDIF                         

      ENDIF

      ! Output the data for this time
      IF (output_type(1:3).EQ. 'mm5') THEN
        IF (model_code(m) .GT. 1) THEN
          PRINT *, 'Writing MM5v3 3D state variables for i4time =',i4time_valid
          CALL output_mm5v3_basic(i4time_cycle,i4time_valid, proj, &
               np_ht,np_t,np_u,np_v,np_rh, &
               ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
               z3d, t3d, u3d, v3d, rh3d, slp, topo, &
               ext_data_path, output_name(m),MM5MODE_NEW,istatus)
          IF (model_code(m) .EQ. 3) THEN
            PRINT *, 'Writing MM5v3 2D state variables for i4time =',i4time_valid 
            CALL output_mm5v3_sfc(i4time_cycle,i4time_valid, proj, &
               tsf, usf, vsf, rhsf, &
               ext_data_path, output_name(m),MM5MODE_APPEND,istatus)
          ENDIF
        ELSE 
          PRINT *, 'Writing MM5v3 2D state variables for i4time =',i4time_valid 
           CALL output_mm5v3_sfc(i4time_cycle,i4time_valid, proj, &
               tsf, usf, vsf, rhsf,  &
               ext_data_path, output_name(m), MM5MODE_NEW,istatus) 
        ENDIF
      ELSE IF (output_type(1:3).EQ.'wrf') THEN
        IF (model_code(m) .GT. 1) THEN 
          PRINT *, 'Writing WRF 3D state variables for i4time =',i4time_valid
          CALL output_wrf_basic(i4time_cycle,i4time_valid, proj, &
               np_ht,np_t,np_u,np_v,np_rh, &
               ht_plevels, t_plevels, u_plevels, v_plevels, rh_plevels, &
               z3d, t3d, u3d, v3d, rh3d, slp, topo, &
               ext_data_path, output_name(m),WRFMODE_NEW,istatus)
          IF (model_code(m) .EQ. 3) THEN
            PRINT *, 'Writing WRF 2D state variables for i4time =',i4time_valid 
            CALL output_wrf_sfc(i4time_cycle,i4time_valid, proj, &
               tsf, usf, vsf, rhsf, &
               ext_data_path, output_name(m),WRFMODE_APPEND,istatus)
          ENDIF
        ELSE
          PRINT *, 'Writing WRF 2D state variables for i4time =',i4time_valid 
          CALL output_wrf_sfc(i4time_cycle,i4time_valid, proj, &
               tsf, usf, vsf, rhsf,  &
               ext_data_path, output_name(m), WRFMODE_NEW,istatus)
        ENDIF                                                                   
      ELSE IF (output_type(1:4).EQ.'rams') THEN
        print *, 'rams output coming soon'
        stop 
      ELSE
        PRINT *,'Unrecognized output format requested: ', output_type
        STOP
      ENDIF
   
      ! If we have surface data, compute a fake SST field
      ! from the wet bulb temperature for this time.  However,
      ! we want to average wet bulb T over all time periods and write
      ! one SST field per run
      IF ( (model_code(m) .EQ. 1) .OR. (model_code(m).EQ.3))THEN 
        fsst(:,:) = 300. ! initialize to something not too bad
        DO jy = 1,proj%ny
          DO ix = 1,proj%nx
            ! Compute tdsf (deg C)
            tdsf_c = dwpt_laps(tsf(ix,jy)-273.15,rhsf(ix,jy))
           
            ! Estimate psf from topo 
            psf_mb = 1013.25*(1.-(topo(ix,jy)*2.257e-5))**(5.259)
            
            ! Compute fsst (wet bulb T) from tsf/rhsf/psf 
            fsst(ix,jy) = twet_fast(tsf(ix,jy)-273.15,tdsf_c,psf_mb)+273.15
          ENDDO
        ENDDO 
        ! Add fsst to fsstsum and increment nfssttimes
        fsstsum(:,:) = fsstsum(:,:) + fsst(:,:) 
        nfssttimes = nfssttimes + 1 
      ENDIF
    ENDDO output_time_loop 
 
    ! If nfssttimes > 0, compute mean fsst and output   
    IF (nfssttimes .GT. 0) THEN

      fsst(:,:) = fsstsum(:,:)/FLOAT(nfssttimes)
      print *, 'Min/Max Fake SST: ', minval(fsst),maxval(fsst)

      ! Call output routine for specific model
      IF (output_type(1:3) .EQ. 'mm5') THEN
        CALL output_mm5v3_fsst(i4time_cycle,i4time_cycle, proj, &
               fsst, ext_data_path, output_name(m), MM5MODE_NEW,istatus)

      ENDIF
    ENDIF
    ! Deallocate 3 arrays for each state variable.  Two will be needed for the
    ! data at two times bounding the time of interest.  The third is where
    ! the time interpolated values will be saved.

    IF (model_code(m) .GT. 1) THEN
      DEALLOCATE (z3d)
      DEALLOCATE (z3d1)
      DEALLOCATE (z3d2)
      DEALLOCATE (t3d)
      DEALLOCATE (t3d1)
      DEALLOCATE (t3d2)
      DEALLOCATE (rh3d)
      DEALLOCATE (rh3d1)
      DEALLOCATE (rh3d2)
      DEALLOCATE (u3d)
      DEALLOCATE (u3d1)
      DEALLOCATE (u3d2)
      DEALLOCATE (v3d)
      DEALLOCATE (v3d1)
      DEALLOCATE (v3d2)
      DEALLOCATE (slp)
      DEALLOCATE (slp1)
      DEALLOCATE (slp2)
    ENDIF
    IF ( (model_code(m).EQ.1).OR.(model_code(m).EQ.3)) THEN
      DEALLOCATE (tsf)
      DEALLOCATE (tsf1)
      DEALLOCATE (tsf2)
      DEALLOCATE (rhsf)
      DEALLOCATE (rhsf1)
      DEALLOCATE (rhsf2)
      DEALLOCATE (usf)
      DEALLOCATE (usf1)
      DEALLOCATE (usf2)
      DEALLOCATE (vsf)
      DEALLOCATE (vsf1)
      DEALLOCATE (vsf2)
      DEALLOCATE (fsst)
      DEALLOCATE (fsstsum)
    ENDIF                                                            
    ! Show this cycle for this model as having been processed.
    PRINT '(2A)', '    Updating process log:',TRIM(proclog)
    OPEN(FILE=proclog,UNIT=11,FORM='FORMATTED', &
           ACCESS='SEQUENTIAL',STATUS='REPLACE') 
    WRITE(11, '(A13)') latest_cycle_wfo
    CLOSE(11)   
     
     ! Deallocate arrays specific to this model
     IF (ALLOCATED(topo)) DEALLOCATE(topo) 
     IF (ALLOCATED(i4times_avail)) DEALLOCATE(i4times_avail)
     IF (ALLOCATED(i4times_avail_max)) DEALLOCATE(i4times_avail_max)
     CALL close_wfofile(nfid,istatus)
  ENDDO model_loop
  
END PROGRAM wfoprep
  
