 &deriv_nl
 mode_evap=0,
 l_bogus_radar_w=.true.,
 l_deep_vv=.true.,
 vv_to_height_ratio_Cu=0.05,
 vv_to_height_ratio_Sc=0.05,
 vv_for_St=.01,
 hydrometeor_scale_pcp=1.,
 hydrometeor_scale_cldliq=0.2,
 hydrometeor_scale_cldice=0.002,
 c_z2m='rams'
 thresh_cvr_cty_vv=0.65,
 thresh_cvr_lwc=0.65,
 twet_snow=+1.3,
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
c l_deep_vv - flag that allows control of whether to use the newer method in 
c             'vv.f' that produces deep parabolic profiles spanning the
c             unstable and more stratiform regions of deep convective clouds
c             
c vv_to_height_ratio_Cu - parameter for the cloud omega (vv.f/cloud_bogus_w)
c                         routine (units are 10^-3 inverse seconds)
c                         This is used in both cloud and radar bogusing
c
c vv_to_height_ratio_Sc - parameter for the cloud omega (vv.f/cloud_bogus_w)
c                         routine (units are 10^-3 inverse seconds)
c
c vv_for_St   - parameter for the cloud omega (vv.f/cloud_bogus_w) routine
c                         (units are meters/second)
c
c hydrometeor_scale_cldliq & hydrometeor_scale_cldice & hydrometeor_scale_pcp:	
c    Factors that scale the hydrometeor concentrations for a grid
c    spacing. If set to a positive number it is treated as a simple scaling.
c    If set to a negative number, then (hydrometeor_scale = -hydrometeor_scale_cld[pcp]/dx)
c    Note that dx is in kilometers. This is the scale that LAPS analyzed values
c    are multiplied by to account for sub-grid scale effects.
c
c c_z2m - parameter for converting from radar reflectivity to precipitating
c         hydrometeor concentrations ('rams', 'albers', 'kessler')
c
c thresh_cvr_cty_vv - cloud cover threshold used for cloud type and cloud omega
c                     a lower value will increase the extent and magnitude of
c                     the cloud omega field
c
c thresh_cvr_lwc - cloud cover threshold used for cloud liquid/ice
c
c twet_snow - wet bulb snow melting threshold (degrees C)
