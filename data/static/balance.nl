 &balance_nl
 lrunbal = .true.,
 adv_anal_by_t_min = 0,
 cpads_type = 'snd',
 incl_clom = .true.,
 setdelo0  = .false.,
 c_erru=0.1,
 c_errub=1.5,
 c_errphi=.15,
 c_errphib=30.,
 c_delo=100.,
 comega_smooth=-1.,
/
c
c   adv_anal_by_t_min: time in minutes to advance the analyses
c and background inputs. Output in lapsprd/balance/lt1,lw3,lq3
c will be systime.dat + T. This namelist value is modified using
c sched.pl command line input -p T where T is in minutes past
c the systime value. The variable is only active when nest7grid.
c parms namelist variable c8_project_common = "Airdrop" in
c which case $DATAROOT/time/systime.dat is not rounded to the
c nearest "cycle_time".
c
c   cpads_type: For the Precision Airdrop Systems project that
c indicates the analysis method. snd: new method and recommended.
c pin: old method but available for backwards compatibility and
c comparison statistics.
c
c
c incl_clom: allows us to include or exclude cloud omega from 
c balance package. 
c and background inputs. Output in lapsprd/balance/lt1,lw3,lq3
c will be systime.dat + T. This namelist value is modified using
c which case $DATAROOT/time/systime.dat is not rounded to the
c
c setdelo0: When set = .true. then only continuity adjustments are made
c made and no dynamic balance.
c
c The following parameters are coefficients on the balance weights.
c The default values favor the wind observations. Also mentioned are
c suggested values to favor the height observations:
c
c c_erru    (default=0.1,  use 1.0 to favor heights instead of winds)
c c_errub   (default=1.5)
c c_errphi  (default=.15)  
c c_errphib (default=30.)
c c_delo    (default=100., use 1.0 to favor heights instead of winds)
c
c comega_smooth: When set to a whole number greater than zero, smooth the 
c                cloud omega field upon input with a box filter with a
c                half-width specified by the parameter value. Here are
c                some examples:
c
c                Value                Action
c                -----                ------
c                -1.                  Automated setting via grid spacing
c                 0.                  No smoothing (disabled)
c                 1.                  3x3 Box Filter
c                 2.                  5x5 Box Filter
