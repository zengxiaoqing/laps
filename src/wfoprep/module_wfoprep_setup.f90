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

MODULE wfoprep_setup

  ! This module contains data and subroutines related to the setup of
  ! the wfoprep program. 
  !  
  ! History
  ! -------
  ! Sep 2001 - Original version..............B. Shaw, NOAA/CIRA/FSL

  IMPLICIT NONE

  ! Set default declarations to private 
  PRIVATE

  ! Private variables available only within module
  INTEGER, PARAMETER           :: max_sources = 10
  ! Public variables shared within and outside of this module 
 
  INTEGER, PUBLIC              :: num_models 
  INTEGER, PUBLIC              :: num_formats
  INTEGER, PUBLIC              :: model_code(max_sources)
  INTEGER, PUBLIC              :: model_run_freq(max_sources)
  INTEGER, PUBLIC              :: model_delay(max_sources)
  INTEGER, PUBLIC              :: max_fcst_len(max_sources)
  INTEGER, PUBLIC              :: output_freq(max_sources)
  REAL, PUBLIC                 :: min_time_frac
  REAL, PUBLIC                 :: min_vert_frac
  CHARACTER (LEN=256)          :: fxa_data
  CHARACTER (LEN=132),PUBLIC   :: model_name(max_sources)
  CHARACTER (LEN=256),PUBLIC   :: model_path(max_sources)
  CHARACTER (LEN=10), PUBLIC   :: output_type
  CHARACTER (LEN=32), PUBLIC   :: output_name(max_sources)
  CHARACTER (LEN=256),PUBLIC   :: laps_data_root
  CHARACTER (LEN=256),PUBLIC   :: ext_data_path

  ! Time stuff 
  INTEGER, PUBLIC              :: i4time_now
  INTEGER, PUBLIC              :: hour_now
  INTEGER, PUBLIC              :: year4_now
  INTEGER, PUBLIC              :: year2_now
  INTEGER, PUBLIC              :: day3_now
  CHARACTER(LEN=9),PUBLIC      :: a9time_now
   
  ! Routines publicly availble to all routines using this module
  PUBLIC  setup_wfoprep

CONTAINS
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_wfoprep(istatus)
    ! Calls all of the wfoprep initialization routines

    IMPLICIT NONE
    INTEGER :: istatus

    istatus = 1

    CALL get_data_roots(istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', 'SETUP_DPREP: Problem getting data roots.'
      istatus = 0
      RETURN
    ENDIF

    CALL read_wfoprep_nl(istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', 'SETUP_DPREP: Problem processing namelist.'
      istatus = 0
      RETURN
    ENDIF   

    CALL get_systime(i4time_now,a9time_now,istatus)
    IF (istatus .NE. 1) THEN
      PRINT '(A)', 'SETUP_DPREP: Problem getting systime.'
      istatus = 0
      RETURN
    ELSE
      PRINT '(A,I12,1x,A9)','LAPS i4time/a9time = ', i4time_now,a9time_now
    ENDIF  
    READ (a9time_now(1:2), '(I2)') year2_now
    IF (year2_now .LT. 90) THEN
      year4_now = year2_now + 2000
    ELSE
      year4_now = year2_now + 1900
    ENDIF
    READ (a9time_now(3:5), '(I3)') day3_now
    READ (a9time_now(6:7), '(I2)') hour_now
    RETURN

  END SUBROUTINE setup_wfoprep
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE get_data_roots(istatus)

    ! Sets up the LAPS_DATA_ROOT and EXT_DATA_ROOT.  
    IMPLICIT NONE
    INTEGER                    :: istatus
    CHARACTER(LEN=256)         :: ext_data_root


    PRINT '(A)', 'Setting up data roots for wfoprep.'  

    CALL GETENV("LAPS_DATA_ROOT", laps_data_root)
    IF (laps_data_root(1:10) .EQ. '          ') THEN
      PRINT '(A)', 'LAPS_DATA_ROOT not set.'
      istatus = 0
      RETURN
    ENDIF

    ! See if EXT_DATA_ROOT is set.  This allows the flexibility for
    ! users of the MM5-LAPS distribution to share an extdataroot between
    ! pregrid and wfoprep.  If not set, then set it to the appropriate
    ! LAPS directory.  
    CALL GETENV("EXT_DATA_ROOT", ext_data_root)
    IF (ext_data_root(1:10) .EQ. '          ') THEN
      ext_data_path = TRIM(laps_data_root) // '/lapsprd/dprep/'
    ELSE
      ext_data_path = TRIM(ext_data_root) // '/extprd/'
    ENDIF
    PRINT '(2A)', 'SBNPREP output will be written in ', TRIM(ext_data_path)

    PRINT '(A)', 'Successfully set up data roots.'
    istatus = 1
    RETURN

  END SUBROUTINE get_data_roots
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  SUBROUTINE read_wfoprep_nl(istatus)
  
    ! Subroutine to read the sbnprep.nl namelist file 

    IMPLICIT none

    INTEGER, PARAMETER          :: nlunit = 22
    INTEGER                     :: i
    INTEGER                     :: istatus
    INTEGER                     :: ioerr
    CHARACTER(LEN=256)          :: nlfile 
    LOGICAL                     :: nlfound
    LOGICAL                     :: done_proc
    CHARACTER(LEN=132)          :: model_name_tmp 
    NAMELIST / wfoprep_nl / model_code, model_name, &
                            model_run_freq, model_delay, max_fcst_len, &
                            output_freq,output_type, fxa_data, &
                            output_name,min_time_frac,min_vert_frac

    PRINT '(A)', 'Processing wfoprep.nl namelist information.'
    ! Initialize the namelist values
    model_code(:) = -1
    model_name(:) = '                                                '
    model_path(:) = '                                                '
    model_run_freq(:) = 0 
    model_delay(:) = 0
    max_fcst_len(:) = 0
    output_freq(:) = 0
    output_type = '          '
    min_time_frac = 0.8
    min_vert_frac = 0.5
    ! Build the file name of the namelist and make sure it exists
    nlfile = TRIM(laps_data_root) // '/static/wfoprep.nl'
    INQUIRE(file=nlfile, EXIST=nlfound)
    IF (.NOT. nlfound) THEN
      PRINT '(2a)', 'Namelist file not found: ', trim(nlfile)
      istatus = 0
      RETURN
    ENDIF 
    
    ! Open the file, rewind for safety, then read the namelist data.
    OPEN (FILE=nlfile, UNIT=nlunit, FORM='FORMATTED',STATUS='OLD', &
          IOSTAT=ioerr)
    IF (ioerr .NE. 0) THEN
      PRINT '(A,I4)', 'Error opening namelist file: ', ioerr
      istatus = 0
      RETURN
    ENDIF

    REWIND(nlunit)
    READ (nlunit, NML=wfoprep_nl,IOSTAT=ioerr)
    IF (ioerr .NE. 0) THEN
      PRINT '(A,I4)', 'I/O error reading namelist: ', ioerr
      istatus = 0
      CLOSE(nlunit)
      RETURN
    ELSE
      CLOSE(nlunit)
    ENDIF

    ! Process the data
    num_models = 1
    done_proc = .false.

    ! Get FXA_DATA path
    IF (fxa_data(1:10).EQ. '          ') THEN
      CALL GETENV('FXA_DATA', fxa_data)
      IF (fxa_data(1:10).EQ. '          ') THEN 
        PRINT *, 'Error, no FXA_DATA path set in namelist or environment.'
        STOP
      ENDIF
    ENDIF

    DO WHILE ((.NOT. done_proc).AND.(num_models .LT. max_sources))

      IF (model_code(num_models) .GT. 0)  THEN
        model_name_tmp = model_name(num_models)
        IF (model_name_tmp(1:5).NE.'     ') THEN
          IF (model_run_freq(num_models) .GT. 0) THEN
            PRINT '(A,I3)', 'Will process code', model_code(num_models)
            ! Build path for this model
            model_path(num_models) = TRIM(fxa_data) // '/Grid/SBN/netCDF/' &
                // TRIM(model_name(num_models))
          ELSE
            PRINT '(2A)', 'No model_run_freq set for ', &
                           TRIM(model_name(num_models))
            done_proc = .true.
          ENDIF
        ELSE
          PRINT '(A,I2)', 'No name provided for model code: ', &
                          model_code(num_models)
          done_proc = .true.
        ENDIF
      ELSE 
        PRINT '(A)', 'No more valid models.'
        done_proc = .true.
      ENDIF
      IF (.NOT. done_proc) THEN
        num_models = num_models + 1
      ELSE
        num_models = num_models - 1
      ENDIF
    ENDDO
    IF (num_models .LT. 1) THEN
      PRINT '(A)', 'No valid model sources have been specified!'
      PRINT '(A)', 'Check your wfoprep.nl settings.'
      istatus = 0
      RETURN
    ELSE
      PRINT '(A)', 'wfoprep will process the following sources...'
      PRINT '(A)', 'CODE   NAME                   PATH'
      PRINT '(A)', '----   --------------------   --------------------------------------'
      DO i = 1, num_models
        PRINT '(I4,3x,A20,3x,A)',model_code(i),model_name(i)(1:20), &
              TRIM(model_path(i))
      ENDDO
    ENDIF

    ! Check output format for validity
    IF ( (output_type .NE. 'mm5       ') .AND. &
         (output_type .NE. 'rams      ') .AND. &
         (output_type .NE. 'wrf       ') ) THEN
      PRINT '(2A)', 'Invalid output format requested: ', &
              output_type
      PRINT *, 'Must be one of:  mm5, rams, wrf'
      istatus = 0
      RETURN
    ELSE
      PRINT '(A)', 'wfoprep will output the following format:',output_type
    ENDIF

    PRINT '(A)', 'Successfully processed namelist.'
    
    istatus = 1
    RETURN
  END SUBROUTINE read_wfoprep_nl
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE wfoprep_setup
