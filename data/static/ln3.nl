 &ln3_nl
 PATH_TO_WSI_3D_RADAR='/public/data/radar/wsi/nexrad/netcdf/',
 MSNG_RADAR = 128,
 ICKINT = 15,
 ITOTWAIT = 120,
 IAGETH = 2701,
 /
c
c Author: J. Smart (9-22-98)
c ln3_driver.exe runtime parameters
c
c path_to_wsi_3d_radar = path to mosaic'ed (by WSI Corp.) nexrad products
c                        including layer composit reflectivities, echo top,
c                        and VIL. Only netCDF in FSL's public data base are
c                        valid currently.
c
c msng_radar = as defined by the netCDF fill value (signed)
c
c ckint = check interval used by wait_for_radar (seconds)
c
c itotwait = total wait time, in seconds
c
c iageth = age threshold above which ln3_driver will not wait for data
c          thus assumming that the nimbus data has stopped coming in.
