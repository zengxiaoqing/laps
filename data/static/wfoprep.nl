 &wfoprep_nl
 fxa_data = '/data/fxa',
 model_name = 'CONUS211/Eta','CONUS212/MesoEta','CONUS215/MesoEta',
 model_code = 3,2,1,
 model_run_freq = 12,6,6,
 model_delay =  4,4,4,
 max_fcst_len = 36, 36, 36
 output_freq = 1, 1, 1,
 output_type = 'mm5',
 output_name = 'ETA211','ETA212','ETA215'
 min_vert_frac = 0.8,
 min_time_frac = 0.8,
/

! fxa_data = top directory containing the "Grid/SBN/netCDF" subdirectory.
!    If not specified, the FXA_DATA environment variable will be
!    queried.
!
! model_name:  One entry per source of data to process.  Each entry
!              consists of the subdirectories under
!              "FXA_DATA/Grid/SBN/netCDF where the source will be 
!              found.  Examples:  CONUS212/MesoEta, CONUS211/Eta
!
! model_code:  Integer code, one for each entry in model_name:
!              1 = get surface fields only from this source
!              2 = get upper-air (pressure level) data only
!              3 = get both surface and upper air
!
! model_run_freq:  Integer value for each entry in model_name.  The
!                  value is the frequency in hours at which the
!                  particular source is produced.  For the models
!                  from NCEP, this is typically 6 or 12.
!
! model_delay:  One entry per source, this is the number of hours
!               after the valid time corresponding to the 00hr 
!               forecast after which the entire run is typically
!               available.
!
! max_fcst_len: One entry per source, this tells the program the
!               maximum forecast hour to use from each source.
!
! output_freq:  Integer value for each source. This is the frequency
!               of the wfoprep output files.  For example, the
!               CONUS211/Eta typically has 6-hourly output.  Setting
!               output_freq = 1 causes wfoprep to time interpolate
!               between those 6-hourly output times to one-hourly
!               output.
!
! output_type:  One entry.  String value that specifies which
!                 modeling system to produce output for.
!                 'mm5':  Support MM5v3 Regridder intermediate format.
!                 'wrf':  Support WRF SI gribprep intermediate format.
!                 'rams': Support RAMS RALPH2 format.
!
! output_name:  One string entry per source.  This is the prefix that
!               will be used for the output files.  Typically it will
!               describe the model itself.
!
! min_vert_frac:  Minimum fraction of total pressure levels (0->1.0)
!                 required to be present before the particular forecast
!                 hour will be used.  The higher the setting (1 is max), the
!                 more stringent the requirements are for complete data sets.
!
! min_time_frac:  Like min_vert_frac, but pertains to the number of "good"
!                 forecast hours available compared to total requested.
!
