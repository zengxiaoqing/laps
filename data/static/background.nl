 &background_nl
 bgpaths='/public/data/grids/ruc/hyb_A236/netcdf',
'/public/data/grids/nam/conus211/netcdf',
'',
 bgmodels=5,2,0,
 forecast_length=12,
 use_analysis=.false.,
 cmodel='RUC40_NATIVE','ETA48_CONUS','LAPS',
 itime_inc =0,
 smooth_fields = .false.,
 luse_sfc_bkgd = .false.,
 lgb_only = .false.,
 /
 &interp_controls
 max_delta_interp=6,
 oldest_fcst_interp=12,
 fcst_len_interp=1,
 use_anal_interp=.false.,
 itime_inc_interp=0,
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
c allowable values are:             (cmodel)
c -------------------------------------------------------------------------------------
c        bgmodels = 0 ----> LAPS_FUA,  MODEL_FUA, and  LAPS:
c                           - LAPS_FUA: uses lapsprd/"model"/fua and fsf for backgrounds;
c                             fua/fsf on same domain; veritcal and time interpolation possible.
c                           - MODEL_FUA: uses lapsprd/"model"/fua and fsf from a different
c                             laps domains; horiz and vertical interpolation required; time interp
c                             possible.
c                           - LAPS: not tested, but under construction 
c        bgmodels = 1 ----> RUC60_NATIVE                         (obsolete!)
c        bgmodels = 2 ----> ETA48_CONUS                          (tested)
c                           ORSM_HKO (Hong Kong Observ model     (tested) 
c                                     
c        bgmodels = 3 ----> CWB_20FA_LAMBERT
c                                "          _NF                  (tested, obsolete)
c                                "          _RE                  (tested, obsolete)
c                                "          _NF15                (tested, working)
c                                "          _GFS                 (tested, working)
c                                "          _NF45                (tested, working)
c                                "          _TFS                 (tested, working)
c
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
c        bgmodels = 10 ---> GFS_ISO or RUC_ISO.  This lets us ingest
c                     netCDF files created from GFS or RUC GRIB data run through
c                     Unidatas gribtonc decoders.
c        bgmodels = 11 ---> WRFARW.  This allows raw netcdf output files from
c                     WRF-ARW v2.1 and create the required lga/lgb files.
c                     Does not time interpolation yet, though, so if you don't
c                     have a raw WRF file with an exact match of the lga background
c                     time needed, it returns to the main lga program with a 0
c                     status to force lga to look for the next bgmodel source in your
c                     background.nl list. 
c        bgmodels = 12 ---> ECMWF - two options available as listed below:
c                           ESRL_NETCDF_LL---> global  area netCDF file from ESRL ITS. (Not tested 6-07)
c                           FMI_NETCDF_LL ---> limited area netCDF file with ECMWF data 
c                                              designed specifically for FMI project area.
c                                             (tested and working as of June-2007)
c        bgmodels = 13 ---> GFS 
c                           NAM 
c                     GRIB formatted files for which there is a Vtable.XXX found 
c                     in dir data/static/Variable_Tables. (ECMWF coming soon.)
c ----------------------------------------------------------------------------------------
c
c forecast_length = the length in hrs of the oldest forecast allowed (to be processed by lga)
c                   as a background in laps
c
c use_analysis =  .true. -> forces backgrounds to be produced from model initial times.
c                 .false.-> lga only uses model forecasts, no initial time allowed.

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
c
c luse_sfc_bkgd = true/false.  If true then lga uses 2m model variables for T, Td and P
c                 for lgb 2d fields. If false, lga interpolates hydrostatically
c                 using the 3d grids (laps) to derive sfc variables for lgb.  Currently
c                 only active for bgpath = /public/data/grids/eta/48km_nat_isobaric/netcdf
c                 (or identical type netcdf file) and cmodel = eta48conus (bgmodel = 2).
c                 When = true, lgb fields are only slightly better (ie., verify better with
c                 sfc obs).
c
c Namelist section 2: interp_controls
c     These variables are not active atm.
