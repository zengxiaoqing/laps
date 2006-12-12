 &balance_nl
 lrunbal = .true.,
 adv_anal_by_t_min = 0,
 cpads_type = 'snd',
 incl_clom = .true.,
 setdelo0  = .false.,
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
