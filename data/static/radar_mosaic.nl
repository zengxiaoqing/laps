 &radar_mosaic_nl
 n_radars_wideband=-1,
 n_radars_narrowband=0,
 c_radar_ext='v01','v02','v03','v04','v05','v06','v07','v08','v09','v10','v11','v12','v13','v14','v15',
 i_window=900,
 imosaic_3d=1,
 /

c
c
c  S. Albers / J. Smart
c 
c  n_radars_wideband: corresponds to the number of vxx files
c            if set to -1, then switch over to value of 'max_radars_cmn'
c 
c  n_radars_narrowband: corresponds to the number of rdr files
c
c  c_radar_ext: for example 'v01','v02','v03' (for vxx)  or '001','002' (for rdr)
c
c  i_window: number of seconds to allow data to be mosaic'ed. Currently
c      set high because we only get data about once or twice per hour.
c
c  imosaic_3d: = 0 for vrc output only; =1 for vrz output only; =2 for both.
c
