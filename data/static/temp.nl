 &temp_nl
 l_read_raob_t=.true.,
 l_use_raob_t=.true.,
 mode_adjust_heights=1,
 weight_bkg_const_temp=5e28,
 pres_mix_thresh=20000.,
 rms_thresh_temp=1.0,
 radiometer_ht_temp=1000.,
 max_obs=100000,
 /

c TEMPERATURE (LT1/temp.exe) PARAMETERS
c
c l_read_raob_t - flag to determine whether RAOB data from the 'snd' file are 
c                 read in and potentially used
c
c l_use_raob_t - flag to determine whether to utilize RAOB data from the 'snd' 
c                file, given that the data are read in. This is active only 
c                when 'l_read_raob_t' is set to .true.
c
c mode_adjust_heights - The height field is computed using a hydrostatic 
c                       integration of the temperatures. Flag settings:
c
c                       0 - the reference is the surface pressures 
c                           ('PS' field) from the LSX file.
c
c                       1 - the reference level for the integration is the 
c                           model background 500mb heights. 
c
c                       2 - the reference is the mean sea level pressures
c                           ('MSL' field) from the LSX file.
c
c weight_bkg_const_temp - Weight for Model Background. 
c                         Recommended values: 0. < value <= 5e+28.
c                         This controls how quickly the output values match the
c                         background if far from obs. 
c
c rms_thresh_temp - Threshold for rms fit of analysis to obs (deg K). Values
c                   are normalized relative to RMS instrument error of the
c                   observations. This controls when to stop the successive 
c                   correction iterations at progressively smaller radii of 
c                   influence. Lower values tend to put more detail in the 
c                   analysis in the attempt to fit the obs.
c
c radiometer_ht_temp - Only use radiometer data below this altitude in meters
c                      (AGL)
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
c
c max_obs - Number of obs allowed in data structure

