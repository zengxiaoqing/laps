 &verif_nl
 TYPE_OBS = 'R',
 PATH_TO_RAW_PROFILER = '/public/data/profiler/wind/noaanet/netcdf/',
 PATH_TO_RAW_SOUNDING = '/public/data/raob/netcdf/',
 RAOB_PROCESS_LAG = 7200,
 RAOB_PROCESS_LAG_BAL = 7200,
 VERIF_OUTPUT_DIR = 'noBal','Bal','Bkgd',
 VERIF_MISSING_DATA = -999.0
/
c VERIFICATION PARAMETERS
c type_obs = 'R' use raw obs to verify
c            'Q' use QC'd obs to verify
c path_to_raw_profiler = location of raw profiler data 
c path_to_raw_raob = location of raw raob data 
c raob_process_lag=how many seconds to look back for LAPS/raob files
c user chooses which of next two paths to use
c verif_output_dir = full path to put verif files in OR
c verif_output_ext = put files in $LAPS_DATA_ROOT+verif_output_ext
c verif_output_bal = if balance .eq. 1, put files in $LAPS_DATA_ROOT+verif_output_bal
c verif_missing_data = value in output files for missing data
c raob_process_lag=how many seconds to look back for LAPS/raob files
c verif_output_ext = put files in $LAPS_DATA_ROOT+verif_output_ext
c verif_output_bal = if balance .eq. 1, put files in $LAPS_DATA_ROOT+verif_output_bal
