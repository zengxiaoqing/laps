 &radar_mosaic_nl
 c_radar_mosaic_type='vxx',
 n_radars=-1,
 c_radar_ext='v01','v02','v03','v04','v05','v06','v07','v08','v09','v10','v11','v12','v13','v14','v15',
 i_window=900,
 imosaic_3d=1,
 /

c
c
c  S. Albers / J. Smart
c 
c  c_radar_mosaic_type: input either 'vxx' or 'rdr' 
c      vxx files are in lapsprd/v01, v02, etc.         (wideband   / full volume / Level 2)
c      rdr files are in lapsprd/rdr/001,002, 003, etc. (narrowband / low-level   / Level 3)
c 
c      Pathways to data are automatically built in the program since it is lapsprd.
c
c  n_radars: corresponds to the number of vxx or rdr files
c            if set to -1, then switch over to value of 'max_radars_cmn'
c
c  c_radar_ext: for example 'v01','v02','v03' (for vxx)  or '001','002' (for rdr)
c
c  i_window: number of seconds to allow data to be mosaic'ed. Currently
c      set high because we only get data about once or twice per hour.
c
c  imosaic_3d: = 0 for vrc output only; =1 for vrz output only; =2 for both.
c
