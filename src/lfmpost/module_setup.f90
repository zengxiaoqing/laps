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

MODULE setup

  USE mm5v3_io
  USE wrfsi_static
  USE wrf_netcdf
  USE time_utils
  USE map_utils
  ! Contains routines for reading model setup

  IMPLICIT NONE
  
  INTEGER, PARAMETER            :: max_domains = 10
  ! File/path names
  CHARACTER(LEN=1)              :: domain_num_str
  CHARACTER(LEN=2)              :: domain_num_str2
  CHARACTER(LEN=255)            :: laps_data_root
  CHARACTER(LEN=255)            :: mm5_data_root
  CHARACTER(LEN=255)            :: moad_dataroot
  CHARACTER(LEN=255)            :: lfm_data_root
  CHARACTER(LEN=255)            :: lfmprd_dir
  CHARACTER(LEN=255)            :: data_file
  CHARACTER(LEN=255)            :: terrain_file
  CHARACTER(LEN=10)             :: mtype 
  ! Run configuration
  INTEGER                       :: num_domains
  INTEGER                       :: domain_num
  INTEGER                       :: kprs
  INTEGER                       :: domarg
  REAL, ALLOCATABLE             :: prslvl(:)
  REAL                          :: redp_lvl
  LOGICAL                       :: keep_fdda
  LOGICAL                       :: split_output
  LOGICAL                       :: realtime
  LOGICAL                       :: make_donefile
  INTEGER                       :: max_wait_sec
  LOGICAL                       :: proc_by_file_num
  INTEGER                       :: start_file_num
  INTEGER                       :: stop_file_num
  INTEGER                       :: file_num_inc
  LOGICAL                       :: file_num3
  LOGICAL                       :: make_laps(max_domains)
  LOGICAL                       :: write_to_lapsdir(max_domains)
  LOGICAL                       :: make_v5d(max_domains)
  LOGICAL                       :: make_points(max_domains)
  INTEGER                       :: v5d_compress
  CHARACTER(LEN=32)             :: model_name
  LOGICAL                       :: do_smoothing
  LOGICAL                       :: use_model_pbl 
  LOGICAL                       :: gribsfc(max_domains)
  LOGICAL                       :: gribua(max_domains)
  INTEGER                       :: table_version
  INTEGER                       :: center_id
  INTEGER                       :: subcenter_id
  INTEGER                       :: process_id(max_domains)

  ! Time information
  INTEGER                       :: num_times_avail
  INTEGER                       :: num_times_to_proc
  INTEGER                       :: start_time_index
  INTEGER                       :: stop_time_index
  INTEGER                       :: time_index_inc
  REAL                          :: output_freq_min
  CHARACTER(LEN=24)             :: cycle_date
  CHARACTER(LEN=24),ALLOCATABLE :: times_to_proc(:)
  INTEGER,ALLOCATABLE           :: sim_tstep(:)
  CHARACTER(LEN=24), PARAMETER  :: static_date='1900-01-01_00:00:00.0000'
 
  ! New stuff related to formatting of point forecast files

  CHARACTER(LEN=24),ALLOCATABLE :: point_times(:)
  CHARACTER (LEN=3)             :: point_tz_label
  INTEGER                       :: point_tz_utcoffset
  CHARACTER (LEN=3)             :: point_windspd_units
  CHARACTER (LEN=1)             :: point_temp_units 
  CHARACTER (LEN=5)             :: point_vent_units

  ! Model domain configuration info
  INTEGER                       :: nx
  INTEGER                       :: ny
  INTEGER                       :: ksigh
  INTEGER                       :: ksigf
  TYPE(proj_info)               :: proj
  REAL                          :: grid_spacing
  REAL, ALLOCATABLE             :: terdot(:,:)
  REAL, ALLOCATABLE             :: latdot(:,:)
  REAL, ALLOCATABLE             :: londot(:,:)
  REAL, ALLOCATABLE             :: mapfac_d(:,:)
  REAL, ALLOCATABLE             :: coriolis(:,:)
  REAL                          :: proj_cent_lat
  REAL                          :: proj_cent_lon
  REAL                          :: truelat1
  REAL                          :: truelat2
  REAL                          :: cone_factor
  REAL                          :: pole_point
  REAL                          :: corner_lats(4)
  REAL                          :: corner_lons(4)
  REAL, ALLOCATABLE             :: sigmah ( : )
  REAL, ALLOCATABLE             :: sigmaf ( : )
  CHARACTER(LEN=32)             :: projection
  LOGICAL                       :: clwflag
  LOGICAL                       :: iceflag
  LOGICAL                       :: graupelflag
  REAL                          :: Ptop
  REAL                          :: PmslBase
  REAL                          :: TmslBase
  REAL                          :: dTdlnPBase
  REAL                          :: TisoBase

  ! Stuff for point output
  TYPE point_struct
    CHARACTER(LEN=10)           :: id
    REAL                        :: lat
    REAL                        :: lon
    REAL                        :: i
    REAL                        :: j
    INTEGER                     :: elevation
    REAL                        :: hi_temp
    CHARACTER(LEN=16)           :: hi_temp_time
    REAL                        :: lo_temp
    CHARACTER(LEN=16)           :: lo_temp_time
    REAL                        :: avg_temp
    REAL                        :: avg_dewpt
    REAL                        :: total_pcp
    REAL                        :: total_snow
    INTEGER                     :: output_unit
    CHARACTER(LEN=80)           :: customer
  END TYPE point_struct
  TYPE(point_struct),ALLOCATABLE :: point_rec(:)
  INTEGER                        :: num_points 
 
CONTAINS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE setup_lfmpost
   
    IMPLICIT NONE
    INTEGER  :: lun_data, lun_terrain, status
    LOGICAL   :: file_ready
    CHARACTER(LEN=2) :: domarg_str
    ! Main lfmpost.nl namelist is in LAPS_DATA_ROOT, so we must
    ! have one set
   
    CALL GETENV("LAPS_DATA_ROOT",laps_data_root)

    ! Get the 1st argument which tells us which model type is being
    ! run.  (mtype)

    CALL GETARG(1,mtype)
    IF (mtype .EQ. 'lfmpost.ex') THEN
      !Must be an HP!
      CALL GETARG(2,mtype)
      CALL GETARG(3,domarg_str)
    ELSE  
      CALL GETARG(2, domarg_str)
    ENDIF
   
    READ (domarg_str, '(I2)') domarg
    PRINT *, 'Domain number provided as argument: ',domarg
    IF(domarg .GT. 0) THEN
       domain_num = domarg
    ELSE 
       domain_num = 1
    ENDIF
    IF ((mtype .EQ. 'WRF       ').OR.(mtype .EQ. 'wrf       ')) THEN
      mtype = 'wrf       '
    ELSEIF ((mtype .EQ. 'WRF2      ').OR.(mtype .EQ. 'wrf2      ')) THEN
      mtype = 'wrf2      '
    ELSEIF ((mtype .EQ. 'MM5       ').OR.(mtype .EQ. 'mm5       ')) THEN
      mtype = 'mm5       '
    ELSE
      print *, 'Unrecognized model type provided as first arg: ', mtype
      print *, 'WRF and MM5 are supported.'
      STOP
    ENDIF
    
    IF (mtype(1:3) .eq. "mm5" ) THEN
      CALL GETENV("MM5_DATA_ROOT", mm5_data_root)
      IF (mm5_data_root(1:3).EQ."   ") THEN
        PRINT *, 'SETUP: Need to set MM5_DATA_ROOT environment variable.'
        STOP
      ENDIF

      lfmprd_dir = TRIM(mm5_data_root)//"/mm5prd"
      lfm_data_root = mm5_data_root
    ELSEIF (mtype(1:3) .eq. "wrf") THEN
      CALL  GETENV("MOAD_DATAROOT", moad_dataroot)
      IF (moad_dataroot(1:3).EQ."   ") THEN
        PRINT *, 'SETUP: Need to set MOAD_DATAROOT environment variable.'
        STOP
      ENDIF
      
      lfmprd_dir = TRIM(moad_dataroot)//"/wrfprd"
      lfm_data_root = moad_dataroot
    ENDIF

    ! Lets make sure we get it right
    PRINT *, 'Model prd directory:  ', trim(lfmprd_dir)

    CALL read_namelist

    WRITE(domain_num_str, '(I1)') domain_num
    WRITE(domain_num_str2, '(I2.2)') domain_num

    IF (mtype .eq. 'mm5') THEN
      IF (proc_by_file_num) THEN
        CALL make_data_file_name(mm5_data_root, domain_num_str, split_output,&
                             start_file_num, data_file, file_num3)
      ELSE
        CALL make_data_file_name(mm5_data_root, domain_num_str, split_output,&
                             0,data_file,file_num3)
      ENDIF
      print *, 'Initial data file=',TRIM(data_file)
      CALL make_terrain_file_name(mm5_data_root, domain_num_str,terrain_file)
      print *, 'Terrain file =', TRIM(terrain_file)
      CALL open_mm5v3(terrain_file, lun_terrain, status)
      IF (split_output) THEN
        INQUIRE(FILE=data_file, EXIST=file_ready)
        IF (.NOT.file_ready) CALL io_wait(data_file,max_wait_sec)
      ENDIF
      CALL open_mm5v3(data_file, lun_data, status)
      CALL time_setup(lun_data)
      CALL model_setup(lun_data, lun_terrain)
      CLOSE (lun_data)
      CLOSE (lun_terrain)
    
    ELSEIF(mtype(1:3) .eq. 'wrf') THEN
      ! Set up WRF model
      print *, 'Setting up times for ', mtype
      IF (mtype .EQ. 'wrf') THEN
        CALL make_wrf_file_name(lfmprd_dir, domain_num,0,data_file)
        CALL wrf_time_setup
      ELSEIF(mtype .EQ. 'wrf2') THEN
        CALL wrf2_time_setup
        CALL make_wrf2_file_name(lfmprd_dir, domain_num, &
            times_to_proc(1), data_file)
      ENDIF
      print *, 'Initial WRF file: ', TRIM(data_file)
      INQUIRE(FILE=data_file,EXIST=file_ready)
      IF (.NOT.file_ready) THEN
        CALL wrfio_wait (data_file,max_wait_sec)
      ELSE
        IF (realtime) CALL sleep(60)
      ENDIF
      CALL open_wrfnc(data_file,lun_data,status)
      CALL model_setup(lun_data,lun_terrain)
      CLOSE(lun_data)

    ENDIF
    ! If we want to make points, then set this up now
    IF (make_points(domain_num)) THEN
      CALL init_points(status)
      IF (status .NE. 0) THEN
        PRINT *, 'Point output requested, but cannot be fulfilled.'
        PRINT *, 'Ensure MM5_DATA_ROOT/static/lfmpost_points.txt file exists.'
        make_points(domain_num) = .false.
      ELSE 
        PRINT '(A,I3)', 'Points initialized: ', num_points
      ENDIF
    ENDIF
    print *, 'Ended setup_lfmpost'
  END SUBROUTINE setup_lfmpost
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE read_namelist

    IMPLICIT NONE
    INTEGER, PARAMETER          :: max_levels = 200
    INTEGER                     :: k,unit, nml_unit
    INTEGER                     :: status
    REAL                        :: levels_mb(max_levels)
    LOGICAL                     :: used
    CHARACTER(LEN=255)          :: namelist_file    
    CHARACTER(LEN=32)           :: lfm_name(max_domains)

    NAMELIST /lfmpost_nl/ num_domains, keep_fdda, split_output, levels_mb,&
        redp_lvl,lfm_name , proc_by_file_num, start_file_num, stop_file_num, &
             file_num_inc, file_num3, make_laps,realtime, write_to_lapsdir, &
             make_donefile,make_v5d, v5d_compress,max_wait_sec, do_smoothing, &
             gribsfc,gribua,table_version,center_id,subcenter_id,process_id, &
             make_points,point_tz_utcoffset,point_tz_label,point_windspd_units,&
             point_temp_units, point_vent_units,use_model_pbl

    IF (lfmprd_dir(1:3).EQ. "   ") THEN
       namelist_file = "lfmpost.nl"
    ELSE
       namelist_file = TRIM(lfmprd_dir) // '/../static/lfmpost.nl'
    ENDIF

    ! Get a unit number
    CALL get_file_unit(nml_unit)
    OPEN (UNIT=nml_unit, FILE=TRIM(namelist_file), FORM='FORMATTED', &
          STATUS='OLD',IOSTAT=status)
    IF (status.NE.0) THEN
      PRINT *, 'Error opening namelist file (',TRIM(namelist_file),'):',&
         status
      CALL ABORT
    ENDIF

    ! Initialize some values 
    num_domains = 1
    model_name = '                                '
    kprs = 1
    levels_mb(:)=-1.
    redp_lvl = 0. 
    keep_fdda = .true.
    split_output = .true.
    max_wait_sec = 3600
    proc_by_file_num = .false.
    start_file_num = 0
    stop_file_num = 999
    file_num_inc = 1
    file_num3 = .false.
    make_v5d(:) = .false.
    make_laps(:) = .true.
    make_points(:) = .false.
    v5d_compress = 2
    realtime = .true.
    make_donefile = .true.
    do_smoothing = .true.
    ! Default GRIB settings
    gribsfc(:) = .false.
    gribua(:) = .false.
    table_version = 2
    center_id = 59   ! FSL
    subcenter_id = 2 ! LAPB
    process_id(:) = 255
    write_to_lapsdir(:) = .false.
    point_tz_utcoffset = 0
    point_tz_label = 'UTC'
    point_temp_units = 'F'
    point_windspd_units = 'KTS'
    point_vent_units = 'M^2/S'
    use_model_pbl = .true.
   READ(UNIT=nml_unit, NML=lfmpost_nl)
    CLOSE (nml_unit)

    ! See if we have to have LAPS_DATA_ROOT
    IF (make_laps(domain_num)) THEN
      IF (laps_data_root(1:3) .EQ. "   ") THEN
         PRINT *, 'No LAPS_DATA_ROOT environment variable set!'
         PRINT *, 'Either set LAPS_DATA_ROOT or set MAKE_LAPS = .FALSE.'
         STOP
      ENDIF
    ENDIF
  
    ! Check domain number
    IF (domain_num .GT. num_domains) THEN
      PRINT *, 'You requested domain number: ', domain_num
      PRINT *, 'But lfmpost.nl indicates only ', num_domains, &
                ' domains are configured.'
      STOP
    ENDIF

    ! Count up number of levels requested.  They must be in monotonically 
    ! decreasing order (bottom to top atmospherically)

    check_levels: DO k = 2,max_levels
      IF ( (levels_mb(k).GT.0).AND.(levels_mb(k).LT.levels_mb(k-1)))THEN
        kprs = kprs + 1
      ELSE
        EXIT check_levels
      ENDIF
    ENDDO check_levels

    IF ((kprs .LT. 2).AND.(levels_mb(1).LT.0)) THEN
      PRINT *, 'No valid pressure levels set in levels_mb!'
      CALL ABORT
    ENDIF
    ALLOCATE(prslvl(kprs))
    prslvl(1:kprs) = levels_mb(1:kprs)*100.  ! Convert to Pascals!
    PRINT *, 'Setup Information'
    PRINT *, '-----------------'
    PRINT '(A,I4,A)', 'Interpolating to ', kprs, ' pressure levels.'
    DO k = 1,kprs
      PRINT '(A,F9.1,A)', 'Level: ', prslvl(k), 'Pa'
    ENDDO
    model_name =lfm_name(domain_num)

    ! Check point settings 
    IF (make_points(domain_num)) THEN
      ! Timezone stuff
      IF( (point_tz_utcoffset .LT. -12).OR.&
          (point_tz_utcoffset .GT.  12)) THEN
        PRINT *, 'Bad POINT_TZ_UTCOFFSET specified for points: ', &
            point_tz_utcoffset
        STOP
      ENDIF

      ! Temp Units
      IF ((point_temp_units .EQ. 'F').OR.(point_temp_units.EQ.'f'))THEN
        point_temp_units = 'F'
      ELSEIF((point_temp_units .EQ. 'C').OR.(point_temp_units.EQ.'c'))THEN
        point_temp_units = 'C'
      ELSEIF((point_temp_units .EQ. 'K').OR.(point_temp_units.EQ.'k'))THEN
        point_temp_units = 'K'
      ELSE
        PRINT *, 'Bad POINT_TEMP_UNITS:', point_temp_units
        STOP
      ENDIF
  
      ! Wind units
      IF ((point_windspd_units(1:1).EQ. 'K').OR.&
          (point_windspd_units(1:1).EQ. 'k') ) THEN
         point_windspd_units = 'KTS'
      ELSEIF ((point_windspd_units(1:3).EQ. 'MPH').OR.&
          (point_windspd_units(1:3).EQ. 'mph') ) THEN
         point_windspd_units = 'MPH'
      ELSEIF ((point_windspd_units(1:3).EQ. 'M/S').OR.&
         (point_windspd_units(1:3).EQ. 'm/s') ) THEN
         point_windspd_units = 'M/S'
      ELSE
        PRINT *, 'Bad POINT_WINDSPD_UNITS:', point_windspd_units
        STOP
      ENDIF
   
      ! Ventilation index units
      IF ((point_vent_units(1:5).EQ. 'm^2/s').OR.&
          (point_vent_units(1:5).EQ. 'M^2/S') ) THEN
         point_vent_units = 'm^2/s'
      ELSEIF ((point_vent_units(1:5).EQ. 'kt-ft').OR.&
         (point_vent_units(1:5).EQ. 'KT-FT') ) THEN
         point_vent_units = 'KT-FT'
      ELSE
        PRINT *, 'Bad POINT_VENT_UNITS:', point_windspd_units
        STOP
      ENDIF
    ENDIF

  END SUBROUTINE read_namelist
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE time_setup(lun_data)

    IMPLICIT NONE
    INTEGER,INTENT(IN)          :: lun_data
    INTEGER                     :: min_to_add
    INTEGER                     :: num_to_skip
    INTEGER                     :: status
    INTEGER                     :: sim_start_year
    INTEGER                     :: sim_start_month
    INTEGER                     :: sim_start_day
    INTEGER                     :: sim_start_hour
    INTEGER                     :: sim_start_min
    INTEGER                     :: sim_start_sec
    INTEGER                     :: sim_start_frac
    INTEGER                     :: t
    REAL                        :: sim_stop_min
    REAL                        :: dom_start_min
    REAL                        :: dom_stop_min
    REAL                        :: fdda_start_min
    REAL                        :: fdda_stop_min
    REAL                        :: tapfrq
    REAL                        :: buffrq
    LOGICAL                     :: fdda_on
 
    CHARACTER(LEN=24)           :: initial_date,new_date
    CHARACTER(LEN=255)          :: wrfnl

    CALL get_mm5_time_info(lun_data,sim_start_year, sim_start_month, &
                         sim_start_day, sim_start_hour, &
                         sim_start_min, sim_start_sec, sim_start_frac, &
                         sim_stop_min, dom_start_min, dom_stop_min, &
                         fdda_on, fdda_start_min, fdda_stop_min, &
                         tapfrq, buffrq, status)

    ! Compute total number of times available in this domain
    num_times_avail = NINT(dom_stop_min-dom_start_min)/NINT(tapfrq) + 1    
    print *, 'num_times_avail = ', num_times_avail 
    WRITE(initial_date, &
       '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2,".",I4.4)') &
      sim_start_year,sim_start_month,sim_start_day,sim_start_hour, &
      sim_start_min, sim_start_sec, sim_start_frac
    PRINT '(2A)', 'Simulation start time: ',initial_date
    cycle_date = initial_date
    IF ( (fdda_on).AND.(.NOT.keep_fdda))THEN
      min_to_add = fdda_stop_min - dom_start_min
      CALL geth_newdate(new_date, cycle_date(1:16),min_to_add)
      cycle_date = new_date
    ENDIF
    IF (dom_start_min .GT.0.) THEN
      CALL geth_newdate(new_date, initial_date(1:16), NINT(dom_start_min))
      cycle_date = new_date
    ENDIF
    PRINT '(2A)', 'Model cycle time (excluding FDDA): ',cycle_date
    PRINT '(2A)', 'Domain initial time: ', initial_date
    IF ((.NOT.proc_by_file_num).AND.(split_output)) THEN
      IF ( (fdda_on).AND.(.NOT.keep_fdda).AND.&
           (dom_start_min .LT. fdda_stop_min) ) THEN
        num_to_skip = NINT(fdda_stop_min-dom_start_min)/NINT(tapfrq)
        num_times_to_proc = num_times_avail - num_to_skip
        start_time_index = num_to_skip + 1
      ELSE
        num_to_skip = 0
        num_times_to_proc = num_times_avail
        start_time_index = 1
      ENDIF
      stop_time_index = num_times_avail
      time_index_inc = 1
    ELSE
      num_to_skip = start_file_num 
      stop_file_num = MIN(stop_file_num,num_times_avail-1)
      stop_time_index = stop_file_num + 1
      num_times_to_proc = stop_file_num - start_file_num + 1
      start_time_index = start_file_num + 1
      time_index_inc = file_num_inc
    ENDIF
    ALLOCATE(times_to_proc(num_times_to_proc))
    PRINT *, 'Total number of output times: ', num_times_to_proc
    DO t = 1,num_times_avail, time_index_inc
      min_to_add = (t-1)*NINT(tapfrq)
      IF((proc_by_file_num).AND.(split_output))min_to_add=(t-1)*NINT(buffrq)
       CALL geth_newdate(new_date, initial_date(1:16),min_to_add)
      IF ((t .ge. start_time_index).AND.(t.le.stop_time_index)) THEN
        PRINT *, 'Will Process Date: ', new_date
        times_to_proc(t-num_to_skip) = new_date
      ENDIF
    ENDDO
    PRINT *, ' '

  END SUBROUTINE time_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE model_setup(lun_data,lun_terrain)

    IMPLICIT NONE
    INTEGER              :: k, km1
    INTEGER, INTENT(IN)  :: lun_data
    INTEGER, INTENT(IN)  :: lun_terrain
    INTEGER              :: status
    REAL, ALLOCATABLE    :: tempsigma ( : )
    REAL                 :: lat1,lon1,stdlon,dx_m,dy_m
    CHARACTER(LEN=255)   :: wrfinitfile
    INTEGER              :: luninit
    LOGICAL              :: fileready

    IF (mtype .eq. 'mm5') THEN
      CALL get_mm5_map(lun_terrain,projection, proj_cent_lat, proj_cent_lon, &
                     truelat1, truelat2, cone_factor, pole_point, &
                     grid_spacing, nx, ny, status)
      ALLOCATE (latdot (nx,ny))
      ALLOCATE (londot (nx,ny))
      ALLOCATE (terdot (nx,ny))
      CALL get_mm5_2d(lun_terrain,'LATITDOT ',static_date,latdot,'D   ',status)
      CALL get_mm5_2d(lun_terrain,'LONGIDOT ',static_date,londot,'D   ',status)
      CALL get_mm5_2d(lun_terrain,'TERRAIN  ',static_date,terdot,'D   ',status)
      lat1 = latdot(1,1)
      lon1 = londot(1,1) 
      stdlon = proj_cent_lon
    ELSEIF (mtype(1:3) .eq. 'wrf') THEN
     
     ! Use routines in module wrfsi_static to read from static.wrfsi
     CALL get_wrfsi_static_dims(moad_dataroot,domain_num, nx,ny)
     CALL get_wrfsi_static_proj(moad_dataroot,domain_num, &
          projection, lat1,lon1,dx_m,dy_m, stdlon, truelat1,truelat2)
     IF (dx_m .ne. dy_m) THEN
       PRINT *, 'WRF dx != dy...not ready for this!',dx_m,dy_m
       STOP
     ELSE 
       grid_spacing = dx_m
     ENDIF
     ALLOCATE (latdot (nx,ny))
     ALLOCATE (londot (nx,ny))
     ALLOCATE (terdot (nx,ny))
 
     ! Get the lat/lon for the non-staggered grid
     CALL get_wrfsi_static_latlon(moad_dataroot,domain_num,'N',&
            latdot,londot)
    
     ! Get the terrain height
!     CALL get_wrfsi_static_2d(moad_dataroot, 'avg', terdot)
     wrfinitfile = TRIM(moad_dataroot)//'/wrfprd/wrfinput_d'//&
         domain_num_str2
     INQUIRE(FILE=wrfinitfile,EXIST=fileready)
     IF (.NOT.fileready) THEN
       CALL wrfio_wait(wrfinitfile,max_wait_sec)
     ELSE
       IF (realtime) CALL sleep(60)
     ENDIF
     CALL open_wrfnc(wrfinitfile,luninit,status)
     IF (status .NE. 0) THEN
       PRINT *, 'WRF Input File not found: ', TRIM(wrfinitfile)
       PRINT *, 'Needed for terrain heights!'
       STOP
     ENDIF
     CALL get_wrfnc_2d(luninit, "HGT","A",nx,ny,1,terdot,status)  
     CALL close_wrfnc(luninit)
    ENDIF

    ! Use the map_set routine to set up the projection information structure
    SELECT CASE(projection)
      CASE ('LAMBERT CONFORMAL               ')
        CALL map_set(PROJ_LC,lat1,lon1,grid_spacing, &
                     stdlon,truelat1,truelat2,nx,ny,proj)
      CASE ('POLAR STEREOGRAPHIC             ')
        CALL map_set(PROJ_PS,lat1,lon1,grid_spacing, &
                     stdlon,truelat1,truelat2,nx,ny,proj) 
      CASE ('MERCATOR                        ')
        CALL map_set(PROJ_MERC,lat1,lon1,grid_spacing, &
                     stdlon,truelat1,truelat2,nx,ny,proj) 
    END SELECT

    IF (make_v5d(domain_num)) THEN  
      ALLOCATE (mapfac_d(nx,ny)) 
      ALLOCATE (coriolis(nx,ny))
      IF (mtype .eq. 'mm5') THEN
        CALL get_mm5_2d(lun_terrain,'MAPFACDT ',static_date, &
                     mapfac_d,'D   ',status)
        CALL get_mm5_2d(lun_terrain,'CORIOLIS ',static_date, &
                     coriolis,'D   ',status)
      ELSEIF (mtype(1:3) .eq. 'wrf') THEN
        CALL get_wrfsi_static_2d(moad_dataroot,domain_num, 'mfl', mapfac_d)
        CALL get_wrfsi_static_2d(moad_dataroot,domain_num, 'cph', coriolis)
      ENDIF

    ENDIF
    PRINT *, ' '
    PRINT '(A,2F10.2)', 'Min/Max value of terrain: ', minval(terdot), &
                         maxval(terdot)
    PRINT * , ' '
    corner_lats(1) = latdot(1,1)
    corner_lons(1) = londot(1,1)
    corner_lats(2) = latdot(1,ny)
    corner_lons(2) = londot(1,ny)
    corner_lats(3) = latdot(nx,ny)
    corner_lons(3) = londot(nx,ny)
    corner_lats(4) = latdot(nx,1)
    corner_lons(4) = londot(nx,1)
    PRINT *, ' '
    PRINT *, 'Corner points from dot point lat/lon arrays:'
    PRINT *, '============================================'
    PRINT *, ' '
    PRINT '(F8.3,1x,F8.3,10x,F8.3,1x,F8.3)', &
            corner_lats(2),corner_lons(2), corner_lats(3), corner_lons(3)
    PRINT *, '      (NW)----------------------(NE)'
    PRINT *, '        |                        |'
    PRINT *, '        |                        |'
    PRINT *, '      (SW)----------------------(SE)'                 
    PRINT '(F8.3,1x,F8.3,10x,F8.3,1x,F8.3)', &
            corner_lats(1),corner_lons(1), corner_lats(4), corner_lons(4)
    PRINT *, ' '
    !DEALLOCATE(latdot)
    !DEALLOCATE(londot)
    
    IF (mtype .eq. 'mm5') THEN
      CALL get_mm5_misc(lun_data, ksigh, ptop, PmslBase, TmslBase, &
                      dTdlnPBase, TisoBase, &
                      clwflag, iceflag, graupelflag)
      ksigf = ksigh+1
      ALLOCATE (tempsigma(ksigh))
      ALLOCATE (sigmah (ksigh))
      ALLOCATE (sigmaf (ksigf))
      CALL get_mm5_1d(lun_data, 'SIGMAH   ', static_date, tempsigma, status)
      IF (status .NE. 0) THEN
        print *, 'SIGMAH not in output file!'
        stop
      ENDIF
      ! Re-order sigma levels to be from ground up and compute sigmaf
      DO k = 1, ksigh
        sigmah(k) = tempsigma(ksigh + 1 - k)
      ENDDO
      sigmaf(1) = 1.0
      sigmaf(ksigf) = 0.0
      DO k = 2, ksigh
        km1 = k - 1
        sigmaf(k) = sigmah(km1) - (sigmaf(km1) - sigmah(km1))
      ENDDO    
      DEALLOCATE(tempsigma)
      PRINT '(A,I4)', 'Number of half sigma levels: ', ksigh
      PRINT '(A)', 'LEVEL    SIGFULL   SIGHALF'
      DO k = 1, ksigh
        PRINT '(I4,4x,F8.6,2x,F8.6)', k, sigmaf(k), sigmah(k)
      ENDDO
      PRINT '(I4,4x,F8.6)', ksigf, sigmaf(ksigf)
      PRINT *, ' ' 
    ELSEIF (mtype(1:3) .eq. 'wrf') THEN
      IF (mtype .EQ. 'wrf') THEN
        CALL get_wrf_misc(lun_data,ksigh,ksigf,ptop, clwflag, iceflag, &
                        graupelflag)
      ELSEIF(mtype .EQ. 'wrf2') THEN
        CALL get_wrf2_misc(lun_data,ksigh,ksigf,ptop,clwflag, iceflag, &
                        graupelflag)
      ENDIF
      ALLOCATE(sigmah(ksigh))
      ALLOCATE(sigmaf(ksigf))
      CALL get_wrf_1d(lun_data, 'ZNU',sigmah,status)
      CALL get_wrf_1d(lun_data, 'ZNW',sigmaf,status)
      
      ! Some things that WRF does not provide but are initialized
      ! here in a hardcoded way to maintain compatibility with
      ! original MM5POST core code.
      PmslBase = 100000.  !NOT CURRENTLY USED
      TmslBase = 275.     !NOT CURRENTLY USED
      dTdlnPBase = 50.    ! Only used for downward T extrapolation
      TisoBase =  0.      !NOT CURRENTLY USED
       
    ENDIF

    PRINT '(A,I4)', 'Number of half sigma levels: ', ksigh
    PRINT '(A)', 'LEVEL    SIGFULL   SIGHALF'
    DO k = 1, ksigh
      PRINT '(I4,4x,F8.6,2x,F8.6)', k, sigmaf(k), sigmah(k)
    ENDDO
    PRINT '(I4,4x,F8.6)', ksigf, sigmaf(ksigf)
    PRINT *, ' '
   
    RETURN
  END SUBROUTINE model_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE init_points(status)

    ! Reads the lfmpost_points.txt file, checks validity of each requested
    ! point, initializes their output files, etc.  Returns non-zero status
    ! if point initialization cannot be completed.  
    USE map_utils
    IMPLICIT NONE
    INTEGER, INTENT(OUT)        :: status
    INTEGER, PARAMETER       :: max_points = 1000
    TYPE(point_struct),ALLOCATABLE :: points_temp(:)
    TYPE(point_struct)       :: point
    CHARACTER(LEN=200)       :: pointfile
    INTEGER                  :: pointunit
    LOGICAL                  :: pointfileexists
    LOGICAL                  :: lunused
    INTEGER                  :: ic, jc
    INTEGER                  :: year,year2,month,day,hour,minute,sec,jjj
    CHARACTER(LEN=9 )        :: cyclestr
    CHARACTER(LEN=200)       :: outfile
    INTEGER                  :: outunit
    CHARACTER (LEN=10)       :: id
    REAL                     :: lat
    REAL                     :: lon
    CHARACTER(LEN=80)        :: customer
    REAL                     :: elevation
    REAL, EXTERNAL           :: bint
    CHARACTER(LEN=3)         :: domnum_str
    INTEGER                  :: nestpt, i
    CHARACTER(LEN=13)        :: odate,ndate 
    CHARACTER(LEN=24)        :: tempdate
    CHARACTER(LEN=3)         :: model_title
    status = 0
    num_points = 0
    pointfile = TRIM(lfmprd_dir) // '/../static/lfmpost_points.txt'
    INQUIRE(FILE=pointfile, EXIST=pointfileexists)
    IF (.NOT. pointfileexists) THEN
      PRINT *, 'Point file missing: ', pointfile
      status = 1
    ELSE
      ! Get a logical unit number to use
      CALL get_file_unit(pointunit)
      OPEN(FILE=pointfile, UNIT=pointunit, STATUS='OLD',FORM='FORMATTED', &
           ACCESS='SEQUENTIAL')
      ALLOCATE(points_temp(max_points))
      DO WHILE (num_points .LT. max_points)
        READ(pointunit, '(I2,1x,a10,1x,F8.4,1x,f9.4,1x,a)',END=99) &
         nestpt,id,lat,lon,customer
 
        IF (nestpt .EQ. domain_num) THEN
          point%id = id
          point%lat = lat
          point%lon = lon
          point%customer = customer
          ! Compute the I/J  for this point
        
          CALL latlon_to_ij(proj,point%lat,point%lon,point%i,point%j)
          IF (ABS(point%i-1.).LT..001) point%i = 1.
          IF (ABS(point%i-nx).LT..001) point%i = nx
          IF (ABS(point%j-1.).LT..001) point%j = 1.
          IF (ABS(point%j-ny).LT..001) point%j = ny

          IF ((point%i .GE. 1.).AND.(point%i .LE. nx).AND.&
              (point%j .GE. 1.).AND.(point%j .LE. ny) ) THEN
            PRINT *, 'Initializing point location ', point%id
            ic = NINT(point%i)
            jc = NINT(point%j) 
         
            elevation = bint(point%i,point%j,terdot,nx,ny) 
            point%elevation = NINT(elevation*3.2808)
            point%hi_temp = -999.9
            point%lo_temp = 999.9
            point%hi_temp_time = '00/00/0000 00:00'
            point%lo_temp_time = '00/00/0000 00:00'
            point%avg_temp = 0.
            point%avg_dewpt = 0.
            point%total_pcp = 0.
            point%total_snow = 0.
            ! Create output file name, find a logical unit number, open the
            ! file, and write out the header
            CALL split_date_char(times_to_proc(1),year, &
                                 month,day,hour,minute,sec)
            jjj = COMPUTE_DAY_OF_YEAR(year,month,day)
            year2 = MOD(year,1000)
            WRITE(cyclestr,'(I2.2,I3.3,I2.2,I2.2)') year2,jjj,hour,minute 
            WRITE (domnum_str, '("d",I2.2)') domain_num
            outfile = TRIM(lfmprd_dir) // '/' // domnum_str // &
                      '/points/' // TRIM(point%id) // '_' // cyclestr //  &
                      '_fcst.txt'
            CALL get_file_unit(outunit)
            point%output_unit = outunit
            num_points = num_points + 1
            points_temp(num_points) = point 
            OPEN(FILE=outfile, UNIT=outunit, FORM='FORMATTED', &
               ACCESS='SEQUENTIAL')
            WRITE(outunit,'("****************************************", &
                           &"******************************************************")')
            WRITE(outunit,&
              '("LOCATION: ",A,2x,"LAT: ",F8.4,2x,"LON: ",F9.4,2x, &
             &"I: ",F7.2,2x,"J: ",F7.2)') point%id, point%lat, &
                   point%lon, point%i,point%j
            IF (mtype(1:3) .EQ. 'mm5') THEN 
              model_title = 'MM5'
            ELSEIF (mtype(1:3) .EQ. 'wrf') THEN
              model_title = 'WRF'
            ELSE
              model_title = 'UNK'
            ENDIF
            WRITE(outunit,  &
               '(A3,1x,F4.1," km    FORECAST CYCLE: ",A,2x,"DOM: ",&
               & I2,2x,"MODEL ELEVATION: ",I4)') &
                 model_title,grid_spacing/1000.,cyclestr, domain_num, point%elevation
            WRITE(outunit,'("****************************************", &
                           &"******************************************************")')
            WRITE(outunit, &
     '("DATE       TIME  TMP DPT RH  WIND   CEI VIS  WEATHER  PRECP SNOW VENT   MIXHT PBLWND HM HH Fbg")')
            WRITE(outunit, &
     '(A3,8x,A3,3x,A1,3x,A1,3x,"%",3x,"Dg@",A3,1x,"hft mile",10x,"in",4x,"in",3x,A5,2x,"FtAGL",1x,"Dg@",A3)') & 
         point_tz_label, point_tz_label,point_temp_units,&
         point_temp_units, point_windspd_units,point_vent_units,point_windspd_units

            WRITE(outunit, &
     '("---------- ----- --- --- --- ------ --- ---- -------- ----- ---- ------ ----- ------ -- -- ---")')

          ELSE
            PRINT *, 'Point location ',TRIM(point%id),' outside of domain!'
          ENDIF
        ENDIF
      ENDDO
 99   IF (num_points .EQ. 0.) THEN
        status = 1
        PRINT *, 'No valid points found to process'
      ELSE 
        ALLOCATE(point_rec(num_points))
        point_rec = points_temp(1:num_points)
        DEALLOCATE(points_temp)
        status = 0
      ENDIF
    ENDIF 

    ! Set up point times time zone offset
    ALLOCATE(point_times(num_times_to_proc))
    IF (point_tz_utcoffset .EQ. 0) THEN
      point_times = times_to_proc
    ELSE
      DO i = 1, num_times_to_proc
        tempdate = times_to_proc(i)
        odate = tempdate(1:13)
        CALL geth_newdate(ndate,odate,point_tz_utcoffset)
        tempdate(1:13) = ndate
        point_times(i) = tempdate
      ENDDO
    ENDIF
    RETURN
  END SUBROUTINE init_points
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrf_time_setup

    ! Sets up the list of expected times and files to find/process for
    ! a run of the WRF model by reading the namelist.input file.

    IMPLICIT NONE

    ! Here are the variables declared for being able to read the WRF 1.2
    ! namelist sections required
    INTEGER  :: time_step_max, max_dom, dyn_opt, rk_ord, diff_opt
    INTEGER  :: km_opt, damp_opt, isfflx, ifsnow, icloud, num_soil_layers
    INTEGER  :: spec_bdy_width, spec_zone, relax_zone, tile_sz_x, tile_sz_y
    INTEGER  :: numtiles, debug_level
    NAMELIST /namelist_01/ time_step_max,max_dom,dyn_opt,rk_ord,diff_opt,&
              km_opt, damp_opt, isfflx, ifsnow, icloud, num_soil_layers, &
              spec_bdy_width, spec_zone, relax_zone, tile_sz_x, tile_sz_y, &
              numtiles, debug_level
   
    INTEGER  :: grid_id, level, s_we, e_we, s_sn, e_sn, s_vert, e_vert
    INTEGER  :: time_step_count_output, frames_per_outfile
    INTEGER  :: time_step_count_restart, time_step_begin_restart
    INTEGER  :: time_step_sound
    NAMELIST /namelist_02/ grid_id, level,s_we,e_we,s_sn,e_sn, s_vert,e_vert,&
                           time_step_count_output, frames_per_outfile,&
                           time_step_count_restart, time_step_begin_restart,&
                           time_step_sound


    REAL     :: dx,dy,dt,ztop,zdamp,dampcoef
    LOGICAL  :: non_hydrostatic
    REAL     :: smdiv, emdiv, epssm, khdif, kvdif, mix_cr_len
    REAL     :: radt, bldt, cudt, gmt
    INTEGER  :: julyr, julday
    NAMELIST /namelist_03/ dx,dy,dt,ztop,zdamp,dampcoef,&
                           non_hydrostatic,&
                           smdiv, emdiv, epssm, khdif, kvdif, mix_cr_len,&
                           radt, bldt, cudt, gmt, julyr, julday


    INTEGER  :: start_year, start_month, start_day, start_hour, start_second
    INTEGER  :: start_minute, end_year, end_month, end_day, end_hour, &
                end_minute, end_second
    INTEGER  :: interval_seconds, real_data_init_type
    NAMELIST /namelist_05/start_year,start_month,start_day,start_hour, &
                          start_minute,start_second,end_year, end_month,&
                          end_day, end_hour, &
                          end_minute, end_second, &
                          interval_seconds, real_data_init_type

    CHARACTER(LEN=255)  :: wrfnl
    LOGICAL             :: found_wrfnl,used
    INTEGER             :: nml_unit,unit
    INTEGER             :: t, sec_to_add,status
    CHARACTER(LEN=24)   :: new_date

    found_wrfnl = .false.
    ! First, look for "namelist.input" in lfmprd_dir
    wrfnl = TRIM(lfmprd_dir) // '/namelist.input'
    INQUIRE(FILE=wrfnl, EXIST=found_wrfnl)

    ! If not found, look in lfmprd_dir/../static for wrf.nl
    IF (.NOT. found_wrfnl) THEN
      wrfnl = TRIM(lfmprd_dir) // '/../static/wrf.nl'
      INQUIRE(FILE=wrfnl,EXIST=found_wrfnl)
    ENDIF

    IF (.NOT. found_wrfnl) THEN
      wrfnl = 'namelist.input'
      INQUIRE(FILE=wrfnl,EXIST=found_wrfnl)
    ENDIF
    IF (.NOT. found_wrfnl) THEN
      PRINT *, 'Unable to find valid WRF namelist.input'
      STOP
    ENDIF
 
    ! Presumably, we have a valid file, so open and read the 
    ! four appropriate namelist sections.  

    ! Get a unit number
    CALL get_file_unit(nml_unit)
    OPEN (UNIT=nml_unit, FILE=TRIM(wrfnl), FORM='FORMATTED', &
          STATUS='OLD',IOSTAT=status)
    IF (status.NE.0) THEN
      PRINT *, 'Error opening namelist file (',TRIM(wrfnl),'):',&
         status
      CALL ABORT
    ENDIF
    REWIND(nml_unit)
    READ(nml_unit, NML=namelist_01)
    READ(nml_unit, NML=namelist_02)
    READ(nml_unit, NML=namelist_03)
    READ(nml_unit, NML=namelist_05)
    CLOSE(nml_unit)

    num_times_avail = time_step_max/time_step_count_output + 1
    
    ! For now, hardcode to process all available times
    num_times_to_proc = num_times_avail
    start_time_index = 1
    stop_time_index = num_times_avail
    time_index_inc = 1
 
    ! Fill cycle date
    WRITE(cycle_date, &
     '(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2,".0000")') &
      start_year, start_month, start_day,start_hour, &
      start_minute, start_second

    ALLOCATE (times_to_proc (num_times_to_proc))
    ALLOCATE (sim_tstep   (num_times_to_proc))
    PRINT *, 'Total number of output times: ', num_times_to_proc
    DO t = 1,num_times_avail, time_index_inc
      sec_to_add = (t-1)*NINT(dt)*time_step_count_output
      sim_tstep(t) = NINT(sec_to_add/dt)
      CALL geth_newdate(new_date, cycle_date(1:19),sec_to_add)
      PRINT *, 'Will Process Date: ', new_date,sim_tstep(t)
      times_to_proc(t) = new_date
    ENDDO
    PRINT *, ' '
    output_freq_min = dt * FLOAT(time_step_count_output) / 60.
    

  END SUBROUTINE wrf_time_setup
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  SUBROUTINE wrf2_time_setup
 
    ! This subroutine sets up list of times to process for WRFV2.  It
    ! requires a special "wrftimes.nl" file to be in the
    ! MOAD_DATAROOT/wrfprd area.  This namelist is a subset of the
    ! time control items from the standard WRFv2 namelist.input file
    
    IMPLICIT NONE
    
    ! Namelist entries
    INTEGER,DIMENSION(max_domains) :: start_year,   end_year
    INTEGER,DIMENSION(max_domains) :: start_month,  end_month
    INTEGER,DIMENSION(max_domains) :: start_day,   end_day
    INTEGER,DIMENSION(max_domains) :: start_hour,   end_hour
    INTEGER,DIMENSION(max_domains) :: start_minute, end_minute
    INTEGER,DIMENSION(max_domains) :: start_second, end_second
    INTEGER,DIMENSION(max_domains) :: history_interval
    INTEGER,DIMENSION(max_domains) :: frames_per_outfile

    ! Others
    INTEGER,PARAMETER  :: max_times = 500
    CHARACTER(LEN=24), ALLOCATABLE  :: possible_times(:)
    CHARACTER(LEN=24) :: hdate, hdate_beg, hdate_end,hdate_new
    CHARACTER(LEN=255)              :: wrftimesnl
    INTEGER  :: first_time_index, unit ,status

    NAMELIST /wrftimes/ &
      start_year, start_month, start_day, &
      start_hour, start_minute, start_second, &
      end_year, end_month, end_day, end_hour, end_minute, &
      end_second, history_interval, frames_per_outfile

    ! Initialization of some namelist variables that might
    ! not be present.
    start_minute(:) = 0
    start_second(:) = 0
    end_minute(:)   = 0
    end_second(:)   = 0
    history_interval(:) = -1
    frames_per_outfile(:) = 1

    wrftimesnl = TRIM(moad_dataroot)//'/wrfprd/wrftimes.nl'

    ! Get a unit number
    CALL get_file_unit(unit)
    OPEN (UNIT=unit, FILE=TRIM(wrftimesnl), FORM='FORMATTED', &
          STATUS='OLD',IOSTAT=status)
    IF (status.NE.0) THEN
      PRINT *, 'Error opening namelist file (',TRIM(wrftimesnl),'):',&
         status
      CALL ABORT
    ENDIF
    READ(unit, NML=wrftimes)

    IF (frames_per_outfile(domain_num) .NE. 1) THEN
       print *, 'FRAMES_PER_OUTFILE MUST = 1 for LFMPOST!'
       STOP
    ENDIF

    ALLOCATE(possible_times(max_times))
    num_times_to_proc = 1
    WRITE(hdate,50) &
        start_year(domain_num),start_month(domain_num),start_day(domain_num),&
        start_hour(domain_num),start_minute(domain_num), &
        start_second(domain_num) 
50  FORMAT(I4.4,"-",I2.2,"-",I2.2,"_",I2.2,":",I2.2,":",I2.2,".0000")
    possible_times(1) = hdate
    WRITE(hdate_end,50) end_year(domain_num), end_month(domain_num), &
      end_day(domain_num), end_hour(domain_num), end_minute(domain_num), &
      end_second(domain_num) 
    DO WHILE (hdate .LE. hdate_end)
      CALL geth_newdate(hdate,hdate(1:16),history_interval(domain_num))
      IF (hdate .LE. hdate_end) THEN
        num_times_to_proc = num_times_to_proc+1
        IF (num_times_to_proc .GT. max_times) THEN 
          print *, 'Exceeded maximum input times!'
          STOP 
        ENDIF
        possible_times(num_times_to_proc) = hdate
      ENDIF
    ENDDO
    ALLOCATE(times_to_proc(num_times_to_proc))
    times_to_proc(:) = possible_times(1:num_times_to_proc)
    DEALLOCATE(possible_times)

    ! For now, hardcode to process all available times
    num_times_avail = num_times_to_proc
    start_time_index = 1
    stop_time_index = num_times_to_proc
    time_index_inc = 1
    cycle_date = times_to_proc(1)
    RETURN
  END SUBROUTINE wrf2_time_setup    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
END MODULE setup
