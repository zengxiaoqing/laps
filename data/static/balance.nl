 &balance_nl
 lrunbal = .true.
 adv_anal_by_t_min = 0,
 /
c
c adv_anal_by_t_min: time in minutes to advance the analyses
c and background inputs. Output in lapsprd/balance/lt1,lw3,lq3
c will be systime.dat + T. This namelist value is modified using
c sched.pl command line input -p T where T is in minutes past
c the systime value. The variable is only active when nest7grid.
c parms namelist variable c8_project_common = "Airdrop" in
c which case $DATAROOT/time/systime.dat is not rounded to the
c nearest "cycle_time".
