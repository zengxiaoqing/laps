 &cloud_nl
 l_use_vis=.true.,
 l_use_vis_partial=.true.,
 l_use_vis_add=.true.,
 l_use_39=.true.,
 latency_co2=4000,
 pct_req_lvd_s8a=0.,
 i4_sat_window=1270,
 i4_sat_window_offset=-60,
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
c latency_co2 - Allowed time lag (in seconds) for using CO2-Slicing satellite 
c               data from the CTP file for cloud-top pressure information
c               Setting this to a negative value turns off CO2-Slicing usage
c
c pct_req_lvd_s8a - percent coverage required for IR LVD data
c                   Valid range is 0.-100. 
c
c i4_sat_window - half-width of time window for satellite data
c
c i4_sat_window_offset - offset of sat time window center from systime
