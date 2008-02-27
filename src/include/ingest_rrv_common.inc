cc
c these come from the ingest_rrv_namelist:
c (data/static/ingest_rrv.nl)
cc
      character path_to_raw_rrv*200
      integer nxv01,nyv01,nzv01
      real    dxv01,dyv01
cc
c these come from the ingest_rrv_constants:
c (src/include/ingest_rrv_constants.dat)
cc
      character c_radar_name(max_radars_rrv)*4
      real    radar_la1(max_radars_rrv)
      real    radar_lo1(max_radars_rrv)
      real    radlat(max_radars_rrv)
      real    radlon(max_radars_rrv)

cc nxv01,nyv01 and nzv01 are not passed through
c  common since they are used to dimension arrays.
cc
      common /ingest_rrv_common/
     +path_to_raw_rrv,c_radar_name,dxv01,dyv01,
     +radar_la1,radar_lo1,radlat,radlon
