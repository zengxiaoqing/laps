 &lfmpost_nl
 DOMAIN_NUM = 1,
 LFM_NAME = 'mm5', 
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
 MAKE_LAPS = .TRUE.,
 WRITE_TO_LAPSDIR = .FALSE.,
 MAKE_V5D = .FALSE.,
 V5D_COMPRESS = 1,
 DO_SMOOTHING = .false.,
 GRIBSFC = .TRUE.,
 GRIBUA = .TRUE.,
 TABLE_VERSION = 2,
 CENTER_ID = 59,
 SUBCENTER_ID = 2, 
 PROCESS_ID = 156,
 MAKE_POINTS = .TRUE.,
 POINT_TZ_UTCOFFSET = -7,
 POINT_TZ_LABEL = 'MST',
 POINT_TEMP_UNITS = 'F',
 POINT_WINDSPD_UNITS = 'MPH'
 POINT_VENT_UNITS = 'KT-FT'
&END
c 
c  DOMAIN_NUM:  Set to which domain you wish to process (for nested models)
c      1: Outer (main) grid, 2: first nest, etc.
c  MODEL_NAME:  Corresponds to fsf/fua subdirectory to be used
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
c  MAKE_LAPS:  Causes fsf/fua to be output
c  WRITE_TO_LAPSDIR:  If set to true, the fsf/fua data files will
c   be written to LAPS_DATA_ROOT/fxx/model_name.  If false, they
c   will be written to MODEL_DATA_ROOT/xxxprd/dxx/fxx
c  MAKE_POINTS:  Set to true to use lfmpost_points.txt file to create
c   tabular point forecasts.
c  MAKE_V5D:  Set to true to create Vis5D files.
c  VIS5D_COMPRESS:  Set to 1, 2, or 4 for maximum to no compression
c   for Vis5D files.
c  DO_SMOOTHING:  True causes most fields to be smoothed in output. 
c   Exceptions are precip and cloud fields.
c  GRIB_SFC:  Set to true to output the 2D fields to GRIB.
c  GRIB_UA:  Set to true to output the 3D fields to GRIB.
c  TABLE_VERSION/CENTER_ID/SUBCENTER_ID/PROCESS_ID:  Used to 
c    control these values in the PDS section of the output GRIB
c    files.  NOTE:  GRIB output uses FSL/LAPB (Center 59/Subcenter 2)
c    custom GRIB table.  Contact shaw@fsl.noaa.gov if you need a copy
c    of this GRIB table.  
c  DOMAIN_NUM:  Set to which domain you wish to process (for nested models)
c      1: Outer (main) grid, 2: first nest, etc.
c  MODEL_NAME:  Corresponds to fsf/fua subdirectory to be used
c  LEVELS_MB:  Pressure levels (in descending pressure order,
c    specified in mb) to which data will be interpolated.
c             pressure will use. (If using fsf/fua for LAPS
c             as in LAPS_DATA_ROOT/static/surface_analysis.nl
c  START_FILE_NUM/STOP_FILE_NUM/FILE_NUM_INC:  For MM5 split output
c  MAKE_LAPS:  Causes fsf/fua to be output
c  WRITE_TO_LAPSDIR:  If set to true, the fsf/fua data files will
c   be written to LAPS_DATA_ROOT/fxx/model_name.  If false, they
c   will be written to MODEL_DATA_ROOT/xxxprd/dxx/fxx
c  TABLE_VERSION/CENTER_ID/SUBCENTER_ID/PROCESS_ID:  Used to 
c    files.  NOTE:  GRIB output uses FSL/LAPB (Center 59/Subcenter 2)
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
