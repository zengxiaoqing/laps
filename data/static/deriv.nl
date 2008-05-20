 &deriv_nl
 mode_evap=0,
 l_bogus_radar_w=.true.,
 vv_to_height_ratio_Cu=0.5,
 vv_to_height_ratio_Sc=0.05,
 vv_for_St=.01,
 thresh_cvr_cty_vv=0.65,
 thresh_cvr_lwc=0.65,
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
c                   'get_radar_deriv' was contributed by Adan Teng from CWB
c
c vv_to_height_ratio_Cu - parameter for the cloud omega (vv.f/cloud_bogus_w)
c                         routine (units are inverse seconds)
c                         This is used in both cloud and radar bogusing
c
c vv_to_height_ratio_Sc - parameter for the cloud omega (vv.f/cloud_bogus_w)
c                         routine (units are 10^-3 inverse seconds)
c
c vv_for_St   - parameter for the cloud omega (vv.f/cloud_bogus_w) routine
c                         (units are meters/second)
c
c thresh_cvr_cty_vv - cloud cover threshold used for cloud type and cloud omega
c                     a lower value will increase the extent and magnitude of
c                     the cloud omega field
c
c thresh_cvr_lwc - cloud cover threshold used for cloud liquid/ice
