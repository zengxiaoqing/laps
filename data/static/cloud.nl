 &cloud_nl
 l_use_vis=.true.,
 l_use_vis_partial=.true.,
 l_use_vis_add=.true.,
 l_use_39=.true.,
 l_use_metars=.true.,
 l_use_radar=.true.,
 latency_co2=4000,
 pct_req_lvd_s8a=75.,
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
c l_use_metars - flag for whether to use METARs (surface stations) in the
c                cloud fraction analysis. They are always used for internal 
c                verification.
c
c l_use_radar - flag for whether to use radar data in the cloud analysis. 
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
c i4_sat_window - half-width of time window for satellite data
c
c i4_sat_window_offset - offset of sat time window center from systime
