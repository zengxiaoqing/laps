 &background_nl
 bgpaths='$FXA_DATA/Grid/SBN/netCDF/CONUS211/RUC/',
           '$FXA_DATA/Grid/SBN/netCDF/CONUS211/Eta/',
 bgmodels=4,4,
 max_forecast_delta=6,
 oldest_forecast=12,
 use_analysis=.false.,
 /
c
c bgpaths is a list of paths to background models in order of preference
c bgsfcpaths is the path to the corresponding sfc data, 
c currently this is only used in the case of ruc2 (bgmodel 5)
c
c bgmodels describes the model type of the files found in each path
c allowable values are:
c        bgmodels = 1 ---> RUC (60 km native grid)
c        bgmodels = 2 ---> ETA (48 km conus-c grid)
c        bgmodels = 3 ---> Taiwan FA Model (20km Lambert grid)
c        bgmodels = 4 ---> RUC/ETA (SBN CONUS 211 grid)
c        bgmodels = 5 ---> RUC (40 km native grid)
c        bgmodels = 6 ---> AFWA NOGAPS LL 1 deg grid
c        bgmodels = 7 ---> ETA48_GRIB
c        bgmodels = 8 ---> AFWA AVN LL 1 deg grid
c
c If no forecast valid at the requested time and not older than oldest_forecast
c is available then go to the next bgmodel
c
c By default LAPS will not use an analysis file as it's background, 
c the use_analysis logical changes this behaivior (sometimes ... not
c totally debugged at this point ... 11-15-99)
c

