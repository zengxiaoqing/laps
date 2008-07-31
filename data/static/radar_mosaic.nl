 &radar_mosaic_nl
 n_radars_wideband=-1,
 n_radars_narrowband=0,
 i_window=900,
 mosaic_cycle_time=0,
 imosaic_3d=1,
 /

c
c
c  S. Albers / J. Smart
c 
c  n_radars_wideband: corresponds to the number of vxx files
c            if set to -1, then switch over to value of 'max_radars'
c 
c  n_radars_narrowband: corresponds to the number of narrowband radars
c                       in the 'rdr/xxx/vrc' directories
c
c  i_window: number of seconds to allow data to be mosaic'ed. Currently
c      set high because we only get data about once or twice per hour.
c
c  mosaic_cycle_time: Interval in seconds between multiple mosaic outputs 
c                     within the laps cycle. This should divide evenly into
c                     'laps_cycle_time'. If set to zero, then the mosaic 
c                     defaults to once per laps cycle.
c
c  imosaic_3d: = 0 for vrc output only; =1 for vrz output only; =2 for both.
c              Option 2 is being developed and may not yet be functional.
c              Option 2 is suggested only when n_radars_narrowband =0.
