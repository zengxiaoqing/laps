 &wind_nl
 l_use_raob=.true.,
 l_use_cdw=.true.,
 l_use_radial_vel=.true.,
 weight_bkg_const_wind=5e28,
 rms_thresh_wind=0.5,
 max_pr=1500,
 max_pr_levels=300,
 /

c WIND PARAMETERS
c
c l_use_raob - flag to determine whether to utilize RAOB data from the 'snd' 
c              file
c
c l_use_cdw  - flag to determine whether to utilize cloud drift wind data from
c              the cdw file.
c
c l_use_radial_vel - flag to determine whether to utilize Doppler radial 
c                    velocity data
c
c weight_bkg_const_wind - Weight for Model Background. 
c                         Recommended values: 0. < value <= 1e+30.
c                         This controls how quickly the output values match the
c                         background if far from obs. 
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