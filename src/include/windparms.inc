

!       1 = NO Successive correction, > 1 indicates successive correction
        integer n_iter_wind
        parameter (n_iter_wind = 2)

        real*4 hor_radius
        parameter (hor_radius = 35.)

        integer MAX_PR,MAX_PR_LEVELS
        parameter (MAX_PR = 300)
        parameter (MAX_PR_LEVELS = 128)

!       Weights in wind analysis
        real*4 weight_meso
        parameter (weight_meso = 1.03)

        real*4 weight_sao
        parameter (weight_sao = 1.0)

        real*4 weight_pirep
        parameter (weight_pirep = 1.02)

        real*4 weight_prof
        parameter (weight_prof = 1.01)     ! Must be a unique value

        real*4 weight_radar
        parameter (weight_radar = 0.05)

        integer n_var
        parameter (n_var = 2)

        integer iradar_cycle_time
        parameter (iradar_cycle_time = 600)

        integer i4time_raob_window
        parameter (i4time_raob_window = 0)

!       Weight for Model Background. Recommended values: 0. to 1e+30.
!       This will make the output values match the background if far from obs.
!       A value of zero means this parameter is not active.
        real*4 weight_bkg_const                                 
        parameter (weight_bkg_const = 0.)
!       parameter (weight_bkg_const = 1.0e30)

!       If this is true, then we rotate the first guess and obs arrays to
!       Grid North U+V just prior to the call to 'laps_anl'. They are rotated
!       back to True North after the 'laps_anl' call and before the LW3 file
!       is written out. This flag would be especially useful if LAPS is 
!       being run in polar regions.
        logical l_grid_north
        parameter (l_grid_north = .true.)

