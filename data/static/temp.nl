 &temp_nl
 l_use_raob_t=.true.,
 weight_bkg_const_temp=0.,
 /
c TEMPERATURE (LT1/temp.exe) PARAMETERS
c
c l_use_raob_t - flag to determine whether to utilize RAOB data from the 'snd' 
c                file
c
c weight_bkg_const_temp - Weight for Model Background. Recommended values: 
c                         0. to 1e+30. Best might be 1e29 or 1e28. This will 
c                         make the output values match the background if far 
c                         from obs. A value of zero means this parameter is 
c                         not active.
