 &deriv_nl
 mode_evap=0,
 l_bogus_radar_w=.true.,
 /

c DERIV PARAMETERS
c
c mode_evap - flag for whether to evaporate radar echoes in the subcloud layer
c             (0) means no evaporation
c             (2) means do evaporation for 2D and 3D reflectivity data
c             (3) means do evaporation only for 3D reflectivity data
c
c             this is currently experimental while code is being developed 
c
c l_bogus_radar_w - flag for whether to call 'get_radar_deriv' to recalculate
c                   the cloud omega with consideration of radar data
c
