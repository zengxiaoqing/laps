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
  ! format, or WRFSI/hinterp input.   This is controlled by the output_format
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

  INTEGER,PARAMETER         :: maxl=100
  INTEGER,PARAMETER         :: maxt=100
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
  INTEGER                   :: k
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
  INTEGER                   :: np_t,np_ht,np_u,np_v,np_rh,np_msl
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
  REAL                      :: rh_plevels(maxl)
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
  LOGICAL, ALLOCATABLE      :: ht_inv(:,:)
  LOGICAL, ALLOCATABLE      :: msl_inv(:,:)
  LOGICAL, ALLOCATABLE      :: goodtime_flag (:)
  INTEGER                   :: n_goodtimes
  CHARACTER(LEN=10)         :: mslname
  INTEGER, ALLOCATABLE      :: time_index(:) 
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
  REAL, ALLOCATABLE         :: topo(:,:)
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
      PRINT '(A)', '       NOT AVAILABLE...SKIPPING.'
      CYCLE model_loop
    ENDIF

    ! Open the file
    CALL open_wfofile(modelfile_wfo,nfid,istatus)
    IF (istatus .NE. 1)THEN
      PRINT *, 'Skipping due to problem opening file: ', TRIM(modelfile_wfo)
      CYCLE model_loop
    ENDIF
   
    ! Get grid/projection info
    CALL get_wfomodel_proj(nfid,model_name(m),proj,istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Skipping due to problem getting projection info.'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF

    ! Process time info
    CALL get_wfomodel_fcsttimes(nfid,ntimes,fcstsec,istatus)
    IF (istatus.NE.1) THEN
      PRINT *, 'Problem getting time information.  Skipping.'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF
    ALLOCATE(i4times_avail_max(ntimes))
    i4times_avail_max(1:ntimes) = i4time_cycle + fcstsec(1:ntimes)
    IF (i4time_last .GT. i4times_avail_max(ntimes)) THEN
      PRINT *, 'This source does not support max_fcst length of ',max_fcst_len(m)
      PRINT *, 'Resetting to ', (i4times_avail_max(ntimes)-i4time_cycle)/36
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
      PRINT *,'No temperature data found.  Skipping this model.'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF
  
    ! Height is only mandatory for 3D data sets
    IF (model_code(m).GT.1) THEN
      CALL get_wfomodel_var_levels(nfid,'gh        ',ht_id, nz_ht, ht_levels_c, np_ht, &
                 ht_plevels,ht_kbotp,ht_ktopp, havesfc_ht, ht_ksfc,istatus)
      IF (istatus.NE.1)THEN
        PRINT *,'No height data found for this 3d source. Skipping this model.'
        CALL close_wfofile(nfid,istatus)
        CYCLE model_loop
      ENDIF                           
    ENDIF 
     
    ! Relative humidity
    CALL get_wfomodel_var_levels(nfid,'rh        ',rh_id, nz_rh, rh_levels_c, np_rh, &
                  rh_plevels,rh_kbotp,rh_ktopp, havesfc_rh,rh_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT *,'No relative humidity data found.  Skipping this model.'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF                      

    ! U-wind component
    CALL get_wfomodel_var_levels(nfid,'uw        ',u_id, nz_u, u_levels_c, np_u, &
                  u_plevels, u_kbotp,u_ktopp,havesfc_u,u_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT *,'No u-component wind data found.  Skipping this model.'
      CALL close_wfofile(nfid,istatus)
      CYCLE model_loop
    ENDIF                  

    ! V-wind component
    CALL get_wfomodel_var_levels(nfid,'vw        ',v_id, nz_v, v_levels_c, np_v, &
                  v_plevels,v_kbotp,v_ktopp, havesfc_v,v_ksfc,istatus)
    IF (istatus.NE.1)THEN
      PRINT *,'No v-component wind data found.  Skipping this model.'
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
          PRINT *, 'No mean sea level pressure found!'
          CALL close_wfofile(nfid,istatus)
          CYCLE model_loop
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
         PRINT *, 'Mismatch in pressure levels between T,HT,U, and V.'
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
      PRINT *, 'Problem getting T inventory.'
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF

    IF (ALLOCATED(u_inv)) DEALLOCATE(u_inv)
    ALLOCATE(u_inv(nz_u,ntimes))
    CALL get_wfomodel_var_inv(nfid,'uw        ',nz_u,ntimes,u_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Problem getting U inventory.'
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF      

    IF (ALLOCATED(v_inv)) DEALLOCATE(v_inv)
    ALLOCATE(v_inv(nz_v,ntimes))
    CALL get_wfomodel_var_inv(nfid,'vw        ',nz_v,ntimes,v_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Problem getting V inventory.'
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF       

    IF (ALLOCATED(rh_inv)) DEALLOCATE(rh_inv)
    ALLOCATE(rh_inv(nz_rh,ntimes))
    CALL get_wfomodel_var_inv(nfid,'rh        ',nz_rh,ntimes,rh_inv,istatus)
    IF (istatus .NE. 1) THEN
      PRINT *, 'Problem getting RH inventory.'
      CALL close_wfofile(nfid,istatus)  
      CYCLE model_loop
    ENDIF   
             
    IF (ALLOCATED(ht_inv)) DEALLOCATE(ht_inv)
    IF (model_code(m).GT.1) THEN
      ALLOCATE (ht_inv(nz_ht,ntimes))
      CALL get_wfomodel_var_inv(nfid,'gh        ',nz_ht,ntimes,ht_inv,istatus)
      IF (istatus .NE. 1) THEN
        PRINT *, 'Problem getting HT inventory.'
        CALL close_wfofile(nfid,istatus)
        CYCLE model_loop
      ENDIF   
    ENDIF  
 
    IF (ALLOCATED(msl_inv)) DEALLOCATE(msl_inv)
    IF (model_code(m).GT.1) THEN
      ALLOCATE (msl_inv(nz_msl,ntimes))
      CALL get_wfomodel_var_inv(nfid,mslname,nz_msl,ntimes,msl_inv,istatus)
      IF (istatus .NE. 1) THEN
        PRINT *, 'Problem getting MSLP inventory.'
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
          PRINT *, 'Missing MSL for this time.',msl_ksfc,i
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF

        ! Check height field, requiring all levels to be present!
        goodlevs = 0
        htinvloop: DO k = ht_kbotp, ht_ktopp
          IF (ht_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  htinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_ht)
        IF ((goodpct .LT. 1.0).OR.(.NOT.ht_inv(ht_kbotp,i)).OR.&
            (.NOT.ht_inv(ht_ktopp,i))) THEN
          PRINT *, 'Height inventory failed check:', &
              ht_inv(ht_kbotp:ht_ktopp,i),ht_kbotp,ht_ktopp
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF

        ! Check temperatures
        goodlevs = 0
        tinvloop: DO k = t_kbotp, t_ktopp
          IF (t_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  tinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_t)
        IF ((goodpct .LT. 1.0).OR.(.NOT.t_inv(t_kbotp,i)).OR.&
            (.NOT.t_inv(t_ktopp,i))) THEN
          PRINT *, 'Temperature inventory failed check:', &
              t_inv(t_kbotp:t_ktopp,i)
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF               

        ! Check RH
        goodlevs = 0
        rhinvloop: DO k = rh_kbotp, rh_ktopp
          IF (rh_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  rhinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_rh)
        IF ((goodpct .LT. 1.0).OR.(.NOT.rh_inv(rh_kbotp,i)).OR.&
            (.NOT.rh_inv(rh_ktopp,i))) THEN
          PRINT *, 'RH inventory failed check:', &
              rh_inv(rh_kbotp:rh_ktopp,i)
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF                   

        ! Check U 
        goodlevs = 0
        uinvloop: DO k = u_kbotp, u_ktopp
          IF (u_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  uinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_u)
        IF ((goodpct .LT. 1.0).OR.(.NOT.u_inv(u_kbotp,i)).OR.&
            (.NOT.u_inv(u_ktopp,i))) THEN
          PRINT *, 'U inventory failed check:', &
              u_inv(u_kbotp:u_ktopp,i)
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF                

        ! Check V
        goodlevs = 0
        vinvloop: DO k = v_kbotp, v_ktopp
          IF (v_inv(k,i)) goodlevs = goodlevs + 1
        ENDDO  vinvloop
        goodpct = FLOAT(goodlevs)/FLOAT(np_v)
        IF ((goodpct .LT. 1.0).OR.(.NOT.v_inv(v_kbotp,i)).OR.&
            (.NOT.v_inv(v_ktopp,i))) THEN
          PRINT *, 'V inventory failed check:', &
              v_inv(v_kbotp:v_ktopp,i)
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF 

      ENDIF
      IF ((model_code(m).EQ.1).OR.(model_code(m).EQ.3))THEN
        
        ! We need surface values

        IF (.NOT. t_inv(t_ksfc,i)) THEN
          PRINT *, 'Missing surface T'
          goodtime_flag(i) = .false.
          CYCLE invloop 
        ENDIF

        IF (.NOT. u_inv(u_ksfc,i)) THEN
          PRINT *, 'Missing surface U'
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF    

        IF (.NOT. v_inv(u_ksfc,i)) THEN
          PRINT *, 'Missing surface V'
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF     

        IF (.NOT. rh_inv(rh_ksfc,i)) THEN
          PRINT *, 'Missing surface RH'
          goodtime_flag(i) = .false.
          CYCLE invloop
        ENDIF     

      ENDIF  
      ! If we made it this far without cycling the loop, then
      ! this is a good level
      n_goodtimes = n_goodtimes + 1
      goodtime_flag(i) = .true.
    ENDDO invloop 
           
    ! OK, now lets make sure that 80% ofthe time periods
    ! passed the above check and that the first and last
    ! times are avaible
    goodlevs = 0
    DO i = 1,ntimes_needed
      IF (goodtime_flag(i)) goodlevs = goodlevs + 1
    ENDDO
    goodpct = FLOAT(goodlevs)/FLOAT(ntimes_needed)
    IF ( (goodpct .LT. 0.80).OR.(.NOT.goodtime_flag(1)).OR.&
         (.NOT.goodtime_flag(ntimes_needed)))THEN
      PRINT *, 'Model run failed time inventory checks.  Skipping.'
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
      ALLOCATE (rh3d     (proj%nx,proj%ny,np_rh) )
      ALLOCATE (rh3d1    (proj%nx,proj%ny,np_rh) )
      ALLOCATE (rh3d2    (proj%nx,proj%ny,np_rh) ) 
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
      IF (output_format(1:3).EQ. 'mm5') THEN
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
      ELSE IF (output_format(1:3).EQ.'wrf') THEN
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
      ELSE IF (output_format(1:4).EQ.'rams') THEN
        print *, 'rams output coming soon'
        stop 
      ELSE
        PRINT *,'Unrecognized output format requested: ', output_format
        STOP
      ENDIF
  
    ENDDO output_time_loop 

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
  
