 &background_nl
 bgpaths='/public/data/grids/ruc2/40km_fsl-conus_hybb/netcdf',
'/public/data/grids/eta/48km_nat_isobaric/netcdf',
'',
 bgmodels=5,2,0,
 max_forecast_delta=6,
 oldest_forecast=12,
 use_analysis=.false.,
 cmodel='RUC40_NATIVE','ETA48_CONUS','LAPS',
 itime_inc =0,
 /
c
c bgpaths is a list of paths to background models in order of preference
c bgsfcpaths is the path to the corresponding sfc data, 
c currently this is only used in the case of ruc2 (bgmodel 5)
c
c bgmodels describes the model type of the files found in each path
c allowable values are:
c        bgmodels = 0 ---> LAPS
c        bgmodels = 1 ---> RUC60_NATIVE
c        bgmodels = 2 ---> ETA48_CONUS
c        bgmodels = 3 ---> CWB_20FA_LAMBERT (_NF or _RE)
c        bgmodels = 4 ---> RUC_ or ETA_ (SBN_CONUS_211_GRID)
c        bgmodels = 5 ---> RUC40_NATIVE 
c        bgmodels = 6 ---> AVN (_AFWA_DEGRIB or _FSL_NETCDF)
c        bgmodels = 7 ---> ETA48_GRIB
c        bgmodels = 8 ---> NOGAPS_AFWA_DEGRID
c        bgmodels = 9 ---> NWS_CONUS
c
c If no forecast valid at the requested time and not older than oldest_forecast
c is available then go to the next bgmodel
c
c By default LAPS will not use an initial analysis file as it's background, 
c the use_analysis logical changes this behaivior (sometimes ... not
c totally debugged at this point ... 11-15-99) c
c
c A given value for bgmodels can indicate two different formats for that
c same data.  For example, bgmodels=6 is either "AFWA degrid AVN" or FSL's
c public "FSL Public AVN".  The software keys off these very specific names
c which must be properly defined in this namelist with 'cmodel' (currently
c only valid for Global AVN.
