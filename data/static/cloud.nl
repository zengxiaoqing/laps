 &cloud_nl
 l_use_vis=.true.,
 l_use_39=.false.,
 pct_req_lvd_s8a=0.,
 i4_sat_window=970,
 i4_sat_window_offset=-60,
 /

c CLOUD PARAMETERS
c
c l_use_vis - flag for whether to use visible satellite data
c
c l_use_39 - flag for whether to use 3.9 micron satellite data, the
c            default is .false. It is highly recommended to keep this set to 
c            .false. during the "under construction" software phase.
c
c pct_req_lvd_s8a - percent coverage required for IR LVD data
c                   Valid range is 0.-100. 
c
c i4_sat_window - half-width of time window for satellite data
c
c i4_sat_window_offset - offset of sat time window center from systime
