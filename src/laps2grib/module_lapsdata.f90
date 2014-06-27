MODULE lapsdata
!!
!! Configures the list of LAPS variables and their associated GRIB-2 
!! table entries

  IMPLICIT NONE

  ! 3D Isobaric variable metadata
  TYPE var_iso3d
    INTEGER              :: qbal ! 0=use standard, 1=use balance output
    CHARACTER(LEN=3)     :: ext  ! subdirectory and extension
    CHARACTER(LEN=3)     :: var  ! netcdf variable name
    REAL                 :: pmin ! minimum pressure level to output
    REAL                 :: pmax ! maximum pressure level to output
    REAL                 :: conv_fac ! Used to convert units to GRIB standard
    INTEGER              :: scale_fac ! Scale factor when gribing (power of 10)
    INTEGER              :: discipline 
    INTEGER              :: cat_id
    INTEGER              :: parm_id
  END TYPE var_iso3d

  TYPE var_2d
    CHARACTER(LEN=3)     :: ext
    CHARACTER(LEN=3)     :: var
    REAL                 :: conv_fac
    INTEGER              :: scale_fac
    INTEGER              :: lev1_type
    INTEGER              :: lev1_scale
    INTEGER              :: lev1_value
    INTEGER              :: lev2_type
    INTEGER              :: lev2_scale
    INTEGER              :: lev2_value
    INTEGER              :: discipline
    INTEGER              :: cat_id
    INTEGER              :: parm_id
  END TYPE var_2d

  TYPE var_a2d
    CHARACTER(LEN=3)     :: ext
    CHARACTER(LEN=3)     :: var
    REAL                 :: conv_fac
    INTEGER              :: scale_fac
    INTEGER              :: lev1_type
    INTEGER              :: lev1_scale
    INTEGER              :: lev1_value
    INTEGER              :: lev2_type
    INTEGER              :: lev2_scale
    INTEGER              :: lev2_value
    INTEGER              :: discipline
    INTEGER              :: cat_id
    INTEGER              :: parm_id

    INTEGER              :: eyear
    INTEGER              :: emon
    INTEGER              :: eday
    INTEGER              :: ehour
    INTEGER              :: emin
    INTEGER              :: esec
    INTEGER              :: ntimes
    INTEGER              :: ntimes_miss
    INTEGER              :: stattype
    INTEGER              :: periodtype
    INTEGER              :: etime_unit
    INTEGER              :: etime_value
  END TYPE var_a2d

  INTEGER, PARAMETER     ::  max_iso = 40
  INTEGER, PARAMETER     ::  max_2d  = 50
  INTEGER, PARAMETER     ::  max_a2d = 10

  TYPE(var_iso3d)        :: meta_iso3d(max_iso)
  TYPE(var_2d)           :: meta_2d(max_2d)
  TYPE(var_a2d)          :: meta_accum2d(max_a2d)
  INTEGER                :: n_iso, n_2d, n_accum2d

CONTAINS

  SUBROUTINE get_data_config(laps_data_root,vname)

    IMPLICIT NONE
    CHARACTER(LEN=*),INTENT(IN)    :: laps_data_root
    CHARACTER(LEN=*),INTENT(IN)    :: vname
    CHARACTER(LEN=256)             :: vtab   
    LOGICAL                        :: file_exists
    INTEGER,PARAMETER              :: lun = 10
    CHARACTER(LEN=80)              :: fline
    INTEGER                        :: fstat,line_len
    INTEGER  :: qbal,scale_fac,discipline,cat_id,parm_id
    INTEGER  :: lev1_type,lev1_scale,lev1_value
    INTEGER  :: lev2_type,lev2_scale,lev2_value
    REAL     :: conv_fac, pmin,pmax
    CHARACTER(LEN=3) :: var, ext
    vtab = TRIM(laps_data_root)//'/static/'//TRIM(vname)
    INQUIRE (FILE=vtab,EXIST=file_exists) 

    IF (.NOT. file_exists) THEN 
      PRINT *, "-- Config file not found:  ",TRIM(vtab)
      STOP 'no config file'
    ENDIF

    n_iso = 0
    n_2d = 0
    n_accum2d = 0
    OPEN (FILE=vtab,UNIT=lun,ACCESS='SEQUENTIAL',FORM='FORMATTED',STATUS='OLD')
    readfile: DO 
      READ(lun,'(A)',IOSTAT=fstat) fline
      IF (fstat .EQ. 0) THEN
        line_len = LEN_TRIM(fline)
        IF (fline(1:2) .EQ. '3d') THEN
          READ(fline(3:line_len),*,IOSTAT=fstat) qbal,ext,var,pmin,pmax,conv_fac,scale_fac,&
               discipline,cat_id,parm_id
          IF (fstat .EQ. 0) THEN
            n_iso = n_iso + 1
            IF (n_iso .GT. max_iso) THEN
              PRINT *, "Number of 3D variables requested exceeds max_iso parameter!"
              STOP 'module_lapsdata.f90'
            ENDIF
            meta_iso3d(n_iso)%qbal = qbal
            meta_iso3d(n_iso)%ext = ext
            meta_iso3d(n_iso)%var = var
            meta_iso3d(n_iso)%pmin = pmin
            meta_iso3d(n_iso)%pmax = pmax
            meta_iso3d(n_iso)%conv_fac = conv_fac
            meta_iso3d(n_iso)%scale_fac = scale_fac
            meta_iso3d(n_iso)%discipline = discipline
            meta_iso3d(n_iso)%cat_id = cat_id
            meta_iso3d(n_iso)%parm_id = parm_id
          ELSE
            PRINT *, "-- Format problem with 3d variable line: ", TRIM(fline)
          ENDIF

        ELSEIF(fline(1:2) .EQ. '2d')THEN
          READ(fline(3:line_len),*,IOSTAT=fstat) ext,var,conv_fac,scale_fac,&
                lev1_type,lev1_scale,lev1_value,lev2_type,lev2_scale,lev2_value, &
                discipline,cat_id,parm_id
          IF (fstat .EQ. 0) THEN
            n_2d = n_2d + 1
            IF (n_2d .GT. max_2d) THEN 
              PRINT *, "Number of 2D variables requested exceeds max_2d parameter!"
              STOP 'module_lapsdata.f90'
            ENDIF
            meta_2d(n_2d)%ext = ext
            meta_2d(n_2d)%var = var
            meta_2d(n_2d)%conv_fac = conv_fac
            meta_2d(n_2d)%scale_fac = scale_fac
            meta_2d(n_2d)%lev1_type = lev1_type
            meta_2d(n_2d)%lev1_scale = lev1_scale
            meta_2d(n_2d)%lev1_value = lev1_value
            meta_2d(n_2d)%lev2_type = lev2_type
            meta_2d(n_2d)%lev2_scale = lev2_scale
            meta_2d(n_2d)%lev2_value = lev2_value 
            meta_2d(n_2d)%discipline = discipline
            meta_2d(n_2d)%cat_id = cat_id
            meta_2d(n_2d)%parm_id = parm_id

          ELSE
            PRINT *, "-- Format problem with 2d variable line: ", TRIM(fline) 
          ENDIF

        ELSEIF(fline(1:3) .EQ. 'a2d')THEN
          READ(fline(4:line_len),*,IOSTAT=fstat) ext,var,conv_fac,scale_fac,&
                lev1_type,lev1_scale,lev1_value,lev2_type,lev2_scale,lev2_value, &
                discipline,cat_id,parm_id
          IF (fstat .EQ. 0) THEN
            n_accum2d = n_accum2d + 1

            IF (n_accum2d .GT. max_a2d) THEN 
              PRINT *, "Number of Accumulated 2D variables requested exceeds max_a2d parameter!"
              STOP 'module_lapsdata.f90'
            ENDIF

            meta_accum2d(n_accum2d)%ext = ext
            meta_accum2d(n_accum2d)%var = var
            meta_accum2d(n_accum2d)%conv_fac = conv_fac
            meta_accum2d(n_accum2d)%scale_fac = scale_fac
            meta_accum2d(n_accum2d)%lev1_type = lev1_type
            meta_accum2d(n_accum2d)%lev1_scale = lev1_scale
            meta_accum2d(n_accum2d)%lev1_value = lev1_value
            meta_accum2d(n_accum2d)%lev2_type = lev2_type
            meta_accum2d(n_accum2d)%lev2_scale = lev2_scale
            meta_accum2d(n_accum2d)%lev2_value = lev2_value 
            meta_accum2d(n_accum2d)%discipline = discipline
            meta_accum2d(n_accum2d)%cat_id = cat_id
            meta_accum2d(n_accum2d)%parm_id = parm_id
            ! Octet 35-40, reset to current runtime in laps2grib program
            meta_accum2d(n_accum2d)%eyear = 0
            meta_accum2d(n_accum2d)%emon  = 0
            meta_accum2d(n_accum2d)%eday  = 0
            meta_accum2d(n_accum2d)%ehour = 0
            meta_accum2d(n_accum2d)%emin  = 0
            ! Octet 41-43, 47-49
            meta_accum2d(n_accum2d)%esec  = 0
            meta_accum2d(n_accum2d)%ntimes= 1      !42 indicator of time unit (see octet 49)
            meta_accum2d(n_accum2d)%ntimes_miss= 0 !43 missing precip causes accumulation to reset
            meta_accum2d(n_accum2d)%stattype   = 1 !47 1=Accumulation, 0=Avg
            meta_accum2d(n_accum2d)%periodtype = 2 !48 2=Succesive runs have same start time of accum
            meta_accum2d(n_accum2d)%etime_unit = 0 !49 0=min, 1=hour, 2=day
            ! Octet 50, reset to current runtime in laps2grib program
            meta_accum2d(n_accum2d)%etime_value= 60 !50 length of time range (see octet 49)

          ELSE
            PRINT *, "-- Format problem with Accumulated 2d variable line: ", TRIM(fline) 
          ENDIF

        ENDIF

      ELSEIF (fstat .LT. 0) THEN
        EXIT readfile
      ENDIF

    ENDDO readfile

    IF ( (n_iso .LE. 0).AND.(n_2d .LE. 0)) THEN
      PRINT *, "-- No variables configured in config file. Using default values."
      CALL config_data_static
    ENDIF 
    RETURN
  END SUBROUTINE get_data_config

  SUBROUTINE config_data_static

    ! This sets up the list of variables to
    ! process.  Eventually, we should make this some sort
    ! of user-editable table that gets read in.

    IMPLICIT NONE
   
    ! Initialize the counter for the 3d variables 
    n_iso = 0


    ! Fill in the metadata for each variable, then increment the counter

    ! Temperature
    n_iso = n_iso + 1 
    meta_iso3d(n_iso)%qbal = 0
    meta_iso3d(n_iso)%ext = 'lt1'
    meta_iso3d(n_iso)%var = 't3 '
    meta_iso3d(n_iso)%pmin = 1000.
    meta_iso3d(n_iso)%pmax = 105000.
    meta_iso3d(n_iso)%conv_fac = 1.
    meta_iso3d(n_iso)%scale_fac = 1
    meta_iso3d(n_iso)%discipline = 0
    meta_iso3d(n_iso)%cat_id = 0
    meta_iso3d(n_iso)%parm_id = 0

    ! Height
    n_iso = n_iso + 1
    meta_iso3d(n_iso)%qbal = 0 
    meta_iso3d(n_iso)%ext = 'lt1'
    meta_iso3d(n_iso)%var = 'ht '
    meta_iso3d(n_iso)%pmin = 1000.
    meta_iso3d(n_iso)%pmax = 105000.
    meta_iso3d(n_iso)%conv_fac = 1.
    meta_iso3d(n_iso)%scale_fac = 0
    meta_iso3d(n_iso)%discipline = 0
    meta_iso3d(n_iso)%cat_id = 3
    meta_iso3d(n_iso)%parm_id = 5

    ! RH wrt Liquid
    n_iso = n_iso + 1
    meta_iso3d(n_iso)%qbal = 0
    meta_iso3d(n_iso)%ext = 'lh3'
    meta_iso3d(n_iso)%var = 'rhl'
    meta_iso3d(n_iso)%pmin = 1000.
    meta_iso3d(n_iso)%pmax = 105000.
    meta_iso3d(n_iso)%conv_fac = 1.
    meta_iso3d(n_iso)%scale_fac = 1
    meta_iso3d(n_iso)%discipline = 0
    meta_iso3d(n_iso)%cat_id = 1
    meta_iso3d(n_iso)%parm_id = 1

    ! U wind
    n_iso = n_iso + 1
    meta_iso3d(n_iso)%qbal = 0
    meta_iso3d(n_iso)%ext = 'lw3'
    meta_iso3d(n_iso)%var = 'u3 '
    meta_iso3d(n_iso)%pmin = 1000.
    meta_iso3d(n_iso)%pmax = 105000.
    meta_iso3d(n_iso)%conv_fac = 1.
    meta_iso3d(n_iso)%scale_fac = 1
    meta_iso3d(n_iso)%discipline = 0
    meta_iso3d(n_iso)%cat_id = 2
    meta_iso3d(n_iso)%parm_id = 2

    ! V wind
    n_iso = n_iso + 1
    meta_iso3d(n_iso)%qbal = 0
    meta_iso3d(n_iso)%ext = 'lw3'
    meta_iso3d(n_iso)%var = 'v3 '
    meta_iso3d(n_iso)%pmin = 1000.
    meta_iso3d(n_iso)%pmax = 105000.
    meta_iso3d(n_iso)%conv_fac = 1.
    meta_iso3d(n_iso)%scale_fac = 1
    meta_iso3d(n_iso)%discipline = 0
    meta_iso3d(n_iso)%cat_id = 2
    meta_iso3d(n_iso)%parm_id = 3

 
    ! 2 D variables
    n_2d = 0

    ! MSLP
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'msl'
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 0
    meta_2d(n_2d)%lev1_type = 1
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 0
    meta_2d(n_2d)%lev2_type = 255 
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 3
    meta_2d(n_2d)%parm_id = 1

    ! PSFC
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'ps '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 0
    meta_2d(n_2d)%lev1_type = 1
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 0
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 3
    meta_2d(n_2d)%parm_id = 0

    ! 2m Temp
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 't  '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 1
    meta_2d(n_2d)%lev1_type = 103
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 2
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 0
    meta_2d(n_2d)%parm_id = 0

    ! 2m Dewpoint
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'td '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 1
    meta_2d(n_2d)%lev1_type = 103
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 2
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 0
    meta_2d(n_2d)%parm_id = 6

    ! 2m RH
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'rh '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 0
    meta_2d(n_2d)%lev1_type = 103
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 2
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 1
    meta_2d(n_2d)%parm_id = 1

    ! 10m U
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'u  '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 0
    meta_2d(n_2d)%lev1_type = 103
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 10
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 2
    meta_2d(n_2d)%parm_id = 2

    ! 10m V
    n_2d = n_2d + 1
    meta_2d(n_2d)%ext = 'lsx'
    meta_2d(n_2d)%var = 'v  '
    meta_2d(n_2d)%conv_fac = 1.
    meta_2d(n_2d)%scale_fac = 0
    meta_2d(n_2d)%lev1_type = 103
    meta_2d(n_2d)%lev1_scale = 0
    meta_2d(n_2d)%lev1_value = 10
    meta_2d(n_2d)%lev2_type = 255
    meta_2d(n_2d)%lev2_scale = 255
    meta_2d(n_2d)%lev2_value = 255
    meta_2d(n_2d)%discipline = 0
    meta_2d(n_2d)%cat_id = 2
    meta_2d(n_2d)%parm_id = 3

    RETURN
  END SUBROUTINE config_data_static

END MODULE lapsdata
