 &wind_nl
 l_use_raob=.true.,
 l_use_cdw=.true.,
 l_use_radial_vel=.true.,
 weight_bkg_const_wind=0.,
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
c                         Recommended values: 0. to 1e+30.
c                         This will make the output values match the background
c                         if far from obs. A value of zero means this parameter
c                         is not active.

