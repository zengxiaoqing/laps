
       integer max_rays
       parameter (max_rays=380)

       real*4 radialAzim(max_rays)
       integer*2 VCP, elevationNumber, numRadials

       LOGICAL*1   Z                              ( 460*max_rays)
       LOGICAL*1   V                              ( 920*max_rays)

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
     1                               ,i4time_process
     1                               ,numRadials
