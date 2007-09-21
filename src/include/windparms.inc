

!       1 = NO Successive correction, > 1 indicates successive correction
        integer n_iter_wind
        parameter (n_iter_wind = 2)

        real hor_radius
        parameter (hor_radius = 35.)

!       Weights in wind analysis
        real weight_cdw
        parameter (weight_cdw = 0.04)      ! 5 m/s instrument error

        real weight_sfc
        parameter (weight_sfc = 1.0)

        real weight_pirep
        parameter (weight_pirep = 1.02)

        real weight_prof
        parameter (weight_prof = 1.01)     ! Must be a unique value

!       real weight_radar
!       parameter (weight_radar = 0.25)    ! 2 m/s radial velocity error
!                                          ! tangential velocity assume no
!                                          ! error with exact fit to 1st
!                                          ! pass background
        integer n_var
        parameter (n_var = 2)

        integer iradar_cycle_time
        parameter (iradar_cycle_time = 600)

