 &background_nl
 bgpaths='/public/data/grids/ruc/hyb_A236/netcdf',
'/public/data/grids/eta/48km_nat_isobaric/netcdf',
'',
 bgmodels=5,2,0,
 max_forecast_delta=6,
 oldest_forecast=12,
 forecast_length=1,
 use_analysis=.false.,
 cmodel='RUC40_NATIVE','ETA48_CONUS','LAPS',
 itime_inc =0,
 smooth_fields = .false.,
 /
c
c bgpaths is a list of paths to background models in order of preference
c
c tested bgpaths for SBN grids:
c          '/data/fxa/Grid/SBN/netCDF/CONUS211/RUC/',
c          '/data/fxa/Grid/SBN/netCDF/CONUS211/Eta/',
c          '/data/fxa/Grid/SBN/netCDF/CONUS212/MesoEta/',
c          '/data/fxa/fd/Grid/SBN/netCDF/LATLON/AVN/',
c
c bgmodels describes the model type for the files found in each path
c          this variable works in conjunction with cmodel as indicated.
c 
c allowable values are:    (cmodel)
c        bgmodels = 0 ----> LAPS_FUA,MODEL_FUA,LAPS 
c                          (LAPS not tested)
c        bgmodels = 1 ----> RUC60_NATIVE                         (obsolete!)
c        bgmodels = 2 ----> ETA48_CONUS                          (tested)
c                           ORSM_HKO (Hong Kong Observ model     (tested) 
c                                     
c        bgmodels = 3 ----> CWB_20FA_LAMBERT (_NF or _RE)        (tested)
c        bgmodels = 4 ----> (SBN: RUC, ETA, AVN, MesoEta)        (all tested)
c                                 RUC40_NATIVE,
c                                      ETA48_CONUS,
c                                           AVN_SBN_CYLEQ
c                                                MesoEta_SBN
c        bgmodels = 5 ----> RUC40_NATIVE (50 level, 40km)        (tested)
c        bgmodels = 6 ----> AVN (AVN_AFWA_DEGRIB
c                             or AVN_FSL_NETCDF)(only AVN_FSL_NETCDF tested)
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
c         allowable names are included above with the allowable values of bgmodel.
c         new SBN (bgmodel=4) grid available 5-02 - MesoEta_SBN
c         For example, if bgmodel = 5, then cmodel = RUC40_NATIVE
c                      if bgmodel = 2, then cmodel = ETA48_CONUS
c         If bgmodel = 0 and cmodel = FUA_LAPS then corresponding bgpath must be
c                                     $LAPS_DATA_ROOT/lapsprd/fua/"model_type"; eg mm5.
c                                   = FUA_MODEL then lga will process fua/fsf from
c                                     a different domain. Set bgpath accordingly.
c               
c itime_inc = controls time increment for model background.
c     itime_inc = 0   produce background at the analysis time (t)
c     this parameter is not neccessarily used and should be = 0 for now.
c 
c smooth_fields  
c     Set to .true. to turn on 2dx smoother (should not normally be
c     required).
