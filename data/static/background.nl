 &background_nl
 bgpaths='/public/data/grids/ruc/60km_nat_hybb/netcdf',
         '/public/data/grids/maps/40km_fsl-conus_hybb/netcdf',
         '/public/data/grids/eta/48km_nat_isobaric/netcdf',
         '',
 bgmodels=1,5,2,0,
 /
c
c bgpaths is a list of paths to background models in order of preference
c bgmodels describes the model type of the files found in each path
c allowable values are:
c        bgmodels = 1 ---> RUC (60 km native grid)
c        bgmodels = 2 ---> ETA (48 km conus-c grid)
c        bgmodels = 3 ---> NOGAPS
c        bgmodels = 4 ---> RUC (60 km conus-c grid)
c        bgmodels = 5 ---> RUC (40 km native grid)
c
