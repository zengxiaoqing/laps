
       real radialAzim(MAX_RAY_TILT)
       integer VCP, elevationNumber, numRadials

       integer   Z(MAX_REF_GATES*MAX_RAY_TILT)
       integer   V(MAX_VEL_GATES*MAX_RAY_TILT)

       common/radar_tilt_data/
     1                                siteLat                        
     1                               ,siteLon                        
     1                               ,siteAlt                        
     1                               ,elevationAngle
     1                               ,elevationNumber
     1                               ,VCP
     1                               ,r_nyquist
     1                               ,radialAzim
     1                               ,Z
     1                               ,V
     1                               ,resolutionV
     1                               ,i4time_process
     1                               ,numRadials
