 &cloud_nl
 l_use_vis=.true.,
 l_use_vis_partial=.true.,
 l_use_vis_add=.true.,
 l_use_39=.true.,
 l_use_metars=.true.,
 l_use_radar=.true.,
 l_corr_parallax=.true.,
 latency_co2=4000,
 pct_req_lvd_s8a=75.,
 cld_weight_modelfg=.01,
 echotop_thr_a=4000.,4000.,276.15,
 i4_sat_window=1870,
 i4_sat_window_offset=-60,
 i_varadj=1,
 /

c CLOUD PARAMETERS
c
c l_use_vis - flag for whether to use visible satellite data
c
c l_use_vis_partial - flag for whether to use visible satellite data even
c                     if it covers just part of the domain
c
c l_use_vis_add - flag for whether to use visible satellite data for cloud
c                 building (instead of just cloud clearing)
c
c l_use_39 - flag for whether to use 3.9 micron satellite data, the
c            default is .true. This is still somewhat experimental, so this
c            should be set back to .false. if any problems are suspected.
c
c l_use_metars - flag for whether to use METARs (surface stations) in the
c                cloud fraction analysis. They are always used for internal 
c                verification.
c
c l_use_radar - flag for whether to use radar data in the cloud analysis. 
c
c l_corr_parallax - flag for whether to correct satellite parallax
c
c latency_co2 - Allowed time lag (in seconds) for using CO2-Slicing satellite 
c               data from the CTP file for cloud-top pressure information
c               Setting this to a negative value turns off CO2-Slicing usage
c
c pct_req_lvd_s8a - percent coverage required of IR LVD data for the cloud
c                   analysis to produce any output. Valid range is 0.-100. 
c                   A value of 0. means that IR data are not required and 
c                   the cloud analysis will produce output anyway using the 
c                   other data sources. 
c
c cld_weight_modelfg - Set weight for using model background clouds beyond a 
c                      certain effective radius of influence from the sfc 
c                      obs/pireps
c   cld_weight_modelfg = 0.    ! Model wt inactive, obs used to infinite radius
c   cld_weight_modelfg = 100.  ! Model used beyond  ~40km from nearest obs
c   cld_weight_modelfg = 1.    ! Model used beyond ~100km from nearest obs
c   cld_weight_modelfg = .01   ! Model used beyond ~250km from nearest obs
c   cld_weight_modelfg = .0001 ! Model used beyond ~630km from nearest obs
c 
c
c echotop_thr_a - 3 element array for the echo top ground clutter test
c                 1) height threshold for cold temperatures
c                 2) height threshold for warm temperatures
c                 3) surface air temperature cutoff (K)
c
c i4_sat_window - half-width of time window for satellite data
c
c i4_sat_window_offset - offset of sat time window center from systime
