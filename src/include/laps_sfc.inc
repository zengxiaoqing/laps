c
c.....	LAPS Surface analysis include file.
c............................................
c.....  This file contains variables that are common to many
c.....  parts of the LAPS Surface analysis and that will need 
c.....  to be changed for porting LAPS to different domains.
c.....  Most, if not all, porting changes will be done automatically
c.....  when 'lapsparms.for' is changed.
c
	include 'lapsparms.for'
c
c.....	Define surface analysis parameters in terms of 'lapsparms'
c.....	parameters for easier porting.
c
	parameter(ni = NX_L, nj = NY_L, nk = NZ_L)
	parameter(mxstn = maxstns)
c
	real*4 fnorm(0:ni-1,0:nj-1)
	common/bimode2/ fnorm
c
	real*4 x1a(ni), x2a(nj), y2a(ni,nj)
	integer*4 i_loc(mxstn), j_loc(mxstn)
	common/verif/ x1a, x2a, y2a, i_loc, j_loc
c
	character*9 laps_domain
	parameter(laps_domain = LAPS_DOMAIN_FILE)
	integer*4 iden_ratio
	parameter(iden_ratio = IDEN_RATIO_TEMP)
	parameter(ni_maps = NX_M, nj_maps = NY_M)
c
	parameter(badflag = -99.9)
c
c..........................................................................
c..... Here are parameters to adjust FOR THE SURFACE ANALYSIS ONLY.
c..........................................................................
c
c.......... Quality control:  Set to 1 to skip the QC step, 0 to do the QC.
c
	parameter(iskip = 0)
c
c.......... Pressure reduction:  The main pressure analysis that LAPS produces
c..........     is a reduction to this elevation (m).  For example, the 
c..........     Colorado LAPS uses 1500 m, about the elevation of Denver, 
c..........     since it is representative of the elevations in the region 
c..........     of interest.
c
	parameter(redp_lvl = 1500.)
c
c.......... Surface Theta check:  Check to see that the surface potential
c..........     temperatures are not greater than the potential temperature
c..........     at an upper level.  Set this variable equal to the desired
c..........     upper level:
c
c..........		0 = No sfc theta check done
c..........		7 = Use 700 mb level
c..........	 	5 = Use 500 mb level
c
c..........     Recommended:  Use 700 mb most places, 500 mb over higher
c..........                   terrain areas (like Colorado).
c
	parameter(itheta = 5)
c
c