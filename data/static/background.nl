 &background_nl
 bgpaths='/public/data/grids/ruc2/40km_fsl-conus_hybb/netcdf',
'/public/data/grids/eta/48km_nat_isobaric/netcdf',
'',
 bgmodels=5,2,0,
 max_forecast_delta=6,
 oldest_forecast=12,
 forecast_length=1,
 use_analysis=.false.,
 cmodel='RUC40_NATIVE','ETA48_CONUS','LAPS',
 itime_inc =0,
 /
c
c bgpaths is a list of paths to background models in order of preference
c
c bgmodels describes the model type for the files found in each path
c          this variable works in conjunction with cmodel as indicated.
c 
c allowable values are:    (cmodel)
c        bgmodels = 0 ----> LAPS                                 (not tested)
c        bgmodels = 1 ----> RUC60_NATIVE                         (obsolete!)
c        bgmodels = 2 ----> ETA48_CONUS                          (tested)
c        bgmodels = 3 ----> CWB_20FA_LAMBERT (_NF or _RE)        (not tested)
c        bgmodels = 4 ----> RUC,          ETA,        AVN:       (all tested)
c                          {RUC40_NATIVE, ETA48_CONUS, AVN_SBN_CYLEQ}
c        bgmodels = 5 ----> RUC40_NATIVE                         (tested)
c        bgmodels = 6 ----> AVN (_AFWA_DEGRIB or _FSL_NETCDF)    (AVN_FSL_NETCDF tested)
c        bgmodels = 7 ----> ETA48_GRIB                           (not tested)
c        bgmodels = 8 ----> NOGAPS_AFWA_DEGRID                   (not tested)
c        bgmodels = 9 ----> NWS_CONUS                            (obsolete!)
c
c max_forecast_delta = model forecast output frequency (ie., hrly =1, 3hrly =3, etc)
c                      or controls the maximum allowable background model frequency
c                      for time interpolation.
c
c oldest_forecast = the oldest forecast to be used for making a background (= 12
c                   means a background will not be made from model output beyond
c                   12 hrs.
c
c forecast_length < 0 return only the file which matches i4time_now
c                 = 0 return files for i4time_now and the preceeding forecast
c                 > 0 return all files for i4time_now to >= i4time_now+forecast_length
c
c use_analysis =  forces backgrounds to be produced from model initial times.
c                 Note: this logical does not necessarily work as intended
c
c cmodel: this variable describes the specific type for a given value of bgmodel:
c         allowable values are included above with the allowable values of bgmodel.
c               
c
c itime_inc = controls time increment for model background.
c     itime_inc = 0   produce background at the analysis time (t)
c     this parameter is not neccessarily used and should be = 0 for now.
