 &ln3_nl
 MSNG_RADAR = 128,
 ICKINT = 15,
 ITOTWAIT = 120,
 IAGETH = 2701,
 ISTART = 419,
 JSTART = 207,
 IEND = 863,
 JEND = 479,
 /
c
c ln3_driver.exe runtime parameters
c
c msng_radar = as defined by the netCDF fill value (signed)
c
c ckint = check interval used by wait_for_radar (seconds)
c
c itotwait = total wait time, in seconds
c
c iageth = age threshold above which ln3_driver will not wait for data
c          thus assumming that the nimbus data has stopped coming in.
c
c istart/jstart/iend/jend = estimates of the starting/ending i/j values
c                           for a given domain.  These are computed during
c                           localization by program genln3_bounds.exe.
