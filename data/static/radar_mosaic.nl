 &radar_mosaic_nl
 c_radar_mosaic_type='vxx',
 n_radars=2,
 c_radar_ext='v01','v02','v04',
 i_window=2700,
 imosaic_3d=1,
 /

c
c
c  J. Smart (smart@fsl.noaa.gov 303-497-6597)
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
c  c_radar_ext: presently  'v01','v03','v04' (for vxx)  or '001','002' (for rdr)
c      For the current setup in FSL rrv files are being remapped to LAPS
c  domain for three radars: v01 = 'KCYS', v03 = 'KGLD', and v04 = 'KFTG'.
c
c  i_window: number of seconds to allow data to be mosaic'ed. Currently
c      set high because we only get data about once or twice per hour.
c
c  imosaic_3d: = 0 for vrc output only; =1 for vrz output only; =2 for both.
c
