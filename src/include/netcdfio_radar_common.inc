
       real radialAzim(MAX_RAY_TILT)
       integer VCP, elevationNumber, numRadials

       character*5 radarName

       integer   Z(MAX_REF_GATES*MAX_RAY_TILT)
       integer   V(MAX_VEL_GATES*MAX_RAY_TILT)

       common/radar_tilt_data/
     1                        siteLat                        
     1                       ,siteLon                        
     1                       ,siteAlt                        
     1                       ,elevationAngle
     1                       ,elevationNumber
     1                       ,VCP
     1                       ,r_nyquist
     1                       ,radialAzim
     1                       ,Z
     1                       ,V
     1                       ,resolutionV
     1                       ,gateSizeV,gateSizeZ
     1                       ,firstGateRangeV,firstGateRangeZ
     1                       ,Z_scale, Z_offset
     1                       ,V_scale, V_offset
     1                       ,i4time_process
     1                       ,numRadials
     1                       ,ngates_ref_cdf
     1                       ,ngates_vel_cdf
     1                       ,radarName
