
!       LAPS Parameters
!       This version of lapsparms.for is different from the one in the
!       static directory and contains only those parameters that are
!       independent of which geographical domain you are running laps in.

!       LAPS Grid Dimensions

!       parameter VERTICAL_GRID = 'HEIGHT'

        character*8  VERTICAL_GRID
        parameter (VERTICAL_GRID = 'PRESSURE')

        real*4 height_interval   !  Not actively used
        parameter (height_interval = 1000.)

        integer i2_missing_data
        parameter (i2_missing_data = -99)
        real*4    r_missing_data
        parameter (r_missing_data = +1e37)

        real*4 ref_base
        parameter (ref_base = -10.)

        real*4 ref_base_useable
        parameter (ref_base_useable = 0.)   ! >= ref_base

