 &wind_nl
 l_use_raob=.true.,
 l_use_cdw=.true.,
 l_use_radial_vel=.true.,
 weight_bkg_const_wind=5e28,
 rms_thresh_wind=1.0,
 max_pr=1500,
 max_pr_levels=300,
 i_3d=0,
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
c
c i_3d - Valid values [-1,0,+1]
c        -1 sets l_3d always to .false. (old 3d weighting algorithm)
c         0 allows l_3d to be set automatically during runtime based on
c           estimated computing resources
c        +1 sets l_3d always to .true. (new 3d weighting algorithm)
c
c        Note: the new algorithm produces a more accurate analysis though it
c        takes more CPU and memory resources

