 &temp_nl
 l_use_raob_t=.true.,
 l_adjust_heights=.true.,
 weight_bkg_const_temp=5e28,
 /

c TEMPERATURE (LT1/temp.exe) PARAMETERS
c
c l_use_raob_t - flag to determine whether to utilize RAOB data from the 'snd' 
c                file
c
c l_adjust_heights - The height field is computed using a hydrostatic 
c                    integration of the temperatures. If the flag is .true.,
c                    the reference level for the integration is the model 
c                    background 500mb heights. If .false., the reference is 
c                    the surface pressures ('PS' field) from the LSX file.
c
c weight_bkg_const_temp - Weight for Model Background. 
c                         Recommended values: 0. < value <= 1e+30.
c                         This controls how quickly the output values match the
c                         background if far from obs. 
