 &temp_nl
 l_use_raob_t=.true.,
 l_adjust_heights=.true.,
 weight_bkg_const_temp=5e28,
 pres_mix_thresh=20000.,
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
c
c pres_mix_thresh - Depth of allowed mixed layer when the sfc temperature
c                   analysis is inserted and adiabatically propagated upward. 
c                   This is measured relative to the mean domain surface 
c                   pressure in pascals (i.e. relative to the reference or 
c                   average terrain). Default value is 20000. The surface theta
c                   is not allowed to exceed the pre-existing analyzed theta
c                   at the top of the mixed layer. If the terrain goes above
c                   the mixed layer top, the sfc analysis is not allowed to
c                   warm the column.

