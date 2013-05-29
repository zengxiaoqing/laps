 &wind_nl
 l_use_raob=.true.,
 l_use_cdw=.true.,
 l_use_radial_vel=.true.,
 thresh_2_radarobs_lvl_unfltrd=40,
 thresh_4_radarobs_lvl_unfltrd=75,
 thresh_9_radarobs_lvl_unfltrd=999999,
 thresh_25_radarobs_lvl_unfltrd=110,
 stdev_thresh_radial=99999.,
 weight_bkg_const_wind=5e28,
 weight_radar=0.25,
 rms_thresh_wind=1.0,
 max_pr=3000,
 max_pr_levels=300,
 max_wind_obs=110000,
 r0_barnes_max_m=240000.,
 brns_conv_rate_wind = 0.8,
 qc_thresh_wind_def = 30.,
 qc_thresh_wind_pin = 10.,
 qc_thresh_wind_cdw = 10.,
 qc_thresh_wind_pro = 22.,
 /

c WIND PARAMETERS
c
c l_use_raob - Flag to determine whether to utilize RAOB data from the 'snd' 
c              file in the analysis. If RAOBs aren't analyzed they will be
c              still be utilized for independent data verification.
c
c l_use_cdw  - flag to determine whether to utilize cloud drift wind data from
c              the cdw file.
c
c l_use_radial_vel - flag to determine whether to utilize Doppler radial 
c                    velocity data
c
c thresh_2_radarobs_lvl_unfltrd - threshold number of Doppler obs per level
c                                 for subsampling by factor of 2
c
c thresh_4_radarobs_lvl_unfltrd - threshold number of Doppler obs per level
c                                 for subsampling by factor of 4
c
c thresh_9_radarobs_lvl_unfltrd - threshold number of Doppler obs per level
c                                 for subsampling by factor of 9
c
c thresh_25_radarobs_lvl_unfltrd - threshold number of Doppler obs per level
c                                  for subsampling by factor of 25
c
c stdev_thresh_radial - threshold standard deviation in radial velocity window
c                       (kernel) in meters per second. If measured stdev is 
c                       greater than this then subsampling will not be done for 
c                       this region
c
c weight_bkg_const_wind - Weight for Model Background. 
c                         Recommended values: 0. < value <= 1e+30.
c                         This controls how quickly the output values match the
c                         background if far from obs. Nominally this is equal
c                         to 1e30/err^2, where "1e30" is a scaling constant for
c                         the weights and "err" is the background error in m/s.
c
c weight_radar - weight for derived Doppler wind obs - equivalent to 1/err^2
c                where 'err' is the assumed radial velocity error in m/s. 
c
c rms_thresh_wind - Threshold for rms fit of analysis to obs (non-dimensional).
c                   Values are normalized relative to RMS instrument error of 
c                   the observations. This controls when to stop the 
c                   successive correction iterations at progressively smaller 
c                   radii of influence. Lower values tend to put more detail 
c                   in the analysis in the attempt to fit the obs.
c
c max_pr - Maximum number of wind profiles allowed for 'pro' + 'snd' files.
c
c max_pr_levels - Maximum number of levels per wind profile.
c
c max_wind_obs - Maximum total number of wind observations
c
c r0_barnes_max_m - length scale relating to where obs will start to blend in
c                   to the background
c
c brns_conv_rate_wind - rate of radius reduction for each telescoping Barnes 
c                       interation
c
c The following thresholds QC the winds according to the vector difference
c between observation and background.
c
c qc_thresh_wind_def - default threshold
c qc_thresh_wind_pin - threshold for acars and other point observations
c qc_thresh_wind_cdw - threshold for cloud drift winds
c qc_thresh_wind_pro - threshold for wind profilers
