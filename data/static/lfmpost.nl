 &lfmpost_nl
 KEEP_FDDA = .FALSE., 
 SPLIT_OUTPUT = .TRUE.,
 MAX_WAIT_SEC = 900,
 REALTIME = .true.,
 MAKE_DONEFILE = .true.,
 LEVELS_MB = 1100, 1050, 1000, 950, 900, 850, 800, 750, 700, 650, 600, 550, 500, 450, 400, 350, 300, 250, 200, 150, 100, 
 REDP_LVL = 0.,
 PROC_BY_FILE_NUM = .FALSE.,
 START_FILE_NUM = 0,
 STOP_FILE_NUM = 999,
 FILE_NUM_INC = 1,
 FILE_NUM3 = .FALSE.,
 V5D_COMPRESS = 1,
 DO_SMOOTHING = .false.,
 USE_MODEL_PBL = .TRUE.,
 TABLE_VERSION = 2,
 CENTER_ID = 59,
 SUBCENTER_ID = 2, 
 POINT_TZ_UTCOFFSET = -7,
 POINT_TZ_LABEL = 'MST',
 POINT_TEMP_UNITS = 'F',
 POINT_WINDSPD_UNITS = 'MPH'
 POINT_VENT_UNITS = 'KT-FT'

 NUM_DOMAINS = 1,
 LFM_NAME = 'mm5', 'mm5', 'mm5' ,
 MAKE_LAPS = .TRUE. , .TRUE. , .TRUE.  ,
 WRITE_TO_LAPSDIR = .FALSE.,.FALSE.,.FALSE.,
 MAKE_V5D = .FALSE.,.FALSE.,.FALSE.,
 GRIBSFC = .TRUE.,.TRUE.,.TRUE.,
 GRIBUA = .TRUE.,.TRUE.,.TRUE.
 PROCESS_ID = 1,2,3,
 MAKE_POINTS = .TRUE.,.FALSE.,.FALSE.,

&END
c
c  The following entries apply to all domains:
c 
c  KEEP_FDDA:  Only used for MM5 when model runs in FDDA mode with
c    a pre-forecast nudging period.  If set to false, the lfmpost
c    program will ignore the pre-forecast period.
c  SPLIT_OUTPUT:  For MM5 only.  If set to true, code assumes
c     you have split the MM5 output into multiple files.  If
c     false, it assumes all data is in MMOUT_DOMAINx
c  MAX_WAIT_SEC:  Number of seconds after which the code will
c    abort if it doesn't find the expected output file.
c  REALTIME:  If true, MAX_WAIT_SEC is used to wait for files
c    to be created. 
c  MAKE_DONEFILE: Logical to control whether or not .fsf.done
c    and .fua.done files get created.  Only has an effect 
c    if realtime = .true.
c  LEVELS_MB:  Pressure levels (in descending pressure order,
c    specified in mb) to which data will be interpolated.
c  REDP_LVL:  Height in meters to which the LAPS reduced
c             pressure will use. (If using fsf/fua for LAPS
c             input, this should be set to be the same 
c             as in LAPS_DATA_ROOT/static/surface_analysis.nl
c  PROC_BY_FILE_NUM:  For MM5 split output files, setting to true
c   allows you to process selected files from the run.
c  START_FILE_NUM/STOP_FILE_NUM/FILE_NUM_INC:  For MM5 split output
c   option when PROC_BY_FILE_NUM is true.
c  FILE_NUM3:  Set to true if using MM5 split output that has a 
c    3-digit number sequence.
c  VIS5D_COMPRESS:  Set to 1, 2, or 4 for maximum to no compression
c   for Vis5D files.
c  DO_SMOOTHING:  True causes most fields to be smoothed in output. 
c   Exceptions are precip and cloud fields.
c  USE_MODEL_PBL: Logical.  If true, lfmpost will attempt to use the
c    models PBL height if available and valid.  Otherwise, it will
c    generate a pbl height internally.  For now, WRF always uses an 
c    internally generated value.
c  TABLE_VERSION/CENTER_ID/SUBCENTER_ID:  Used to 
c    control these values in the PDS section of the output GRIB
c    files.  NOTE:  GRIB output uses FSL/LAPB (Center 59/Subcenter 2)
c    custom GRIB table.  Contact shaw@fsl.noaa.gov if you need a copy
c    of this GRIB table.  
c  POINT_TZ_UTCOFFSET:  Integer value to set the offset from UTC to 
c  use when displaying valid times in the point forecast files. For
c  example, if you wanted mountain standard time, you would set this
c  to -7.  For UTC time, leave it set to 0.
c  POINT_TZ_LABEL:  Label for time zone in valid time column.  This
c   is a 3-character string.  Default is "UTC".  For example, if you
c   want Mountain Standard Time for the valid time column, set 
c   this to "MST" and POINT_TZ_UTCOFFSET to -7.
c  POINT_TEMP_UNITS:  Units to use for temp and dewpoint in point
c   forecast files.  Valid values are "F", "C", or "K".
c  POINT_WINDSPD_UNITS:  Sets the units for the windspeed values in
c   the tabular point forecasts.  Valid values are "KTS", "MPH", or 
c    "M/S"
c  POINT_VENT_UNITS: Sets the units for ventilation index in tabular
c    point forecasts.  Valid values are "KT-FT" or "M^2/S"
c  NUM_DOMAINS:  Tells lfmpost.nl how many entries are valid
c    for the namelist parameters that are domain dependent.
c  The following entries have up to 10 entries each, where each
c  entry corresponds to domains 1-10.  
c
c  LFM_NAME:  The name of the LAPS fsf/fua subdirectory
c    to use for LAPS output (MAKE_LAPS & WRITE_TO_LAPSDIR = .TRUE.)
c  MAKE_LAPS:  Set to .true. to create fsf/fua files 
c  WRITE_TO_LAPSDIR:  If .true. writes to
c   $LAPS_DATA_ROOT/lapsprd/fsf/{lfm_name} and fsf/{lfname}.  Otherwise,
c   fsf/fua go in MM5_DATA_ROOT/mm5prd/dxx/fsf and fua
c  MAKE_V5D:  .TRUE. cause vis5d format files to be written to
c    MM5_DATA_ROOT/mm5prd/dxx/v5d
c  GRIBUA:  Output isobaric GRIB upper-air data.  Output goes
c    in MM5_DATA_ROOT/mm5prd/dxx/grib
c  GRIBSFC: Output surface fields to GRIB files.  Output goes
c    in MM5_DATA_ROOT/mm5prd/dxx/grib
c  PROCESS_ID:  Process ID to use in PDS header when making GRIB
c  MAKE_POINTS:  Create tabular text point forecasts if .TRUE.
