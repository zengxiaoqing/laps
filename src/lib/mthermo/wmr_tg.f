cdis    Forecast Systems Laboratory
cdis    NOAA/OAR/ERL/FSL
cdis    325 Broadway
cdis    Boulder, CO     80303
cdis 
cdis    Forecast Research Division
cdis    Local Analysis and Prediction Branch
cdis    LAPS 
cdis 
cdis    This software and its documentation are in the public domain and 
cdis    are furnished "as is."  The United States government, its 
cdis    instrumentalities, officers, employees, and agents make no 
cdis    warranty, express or implied, as to the usefulness of the software 
cdis    and documentation for any purpose.  They assume no responsibility 
cdis    (1) for the use of the software and documentation; or (2) to provide
cdis     technical support to users.
cdis    
cdis    Permission to use, copy, modify, and distribute this software is
cdis    hereby granted, provided that the entire disclaimer notice appears
cdis    in all copies.  All modifications to this software must be clearly
cdis    documented, and are solely the responsibility of the agent making 
cdis    the modifications.  If significant modifications or enhancements 
cdis    are made to this software, the FSL Software Policy Manager  
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
cdis 
	function wmr(p,t)

c   this function approximates the mixing ratio wmr (grams of water
c   vapor per kilogram of dry air) given the pressure p (mb) and the
c   temperature t (celsius). the formula used is given on p. 302 of the
c   smithsonian meteorological tables by roland list (6th edition).

c	baker, schlatter  17-may-1982	  original version.

c   eps = ratio of the mean molecular weight of water (18.016 g/mole)
c         to that of dry air (28.966 g/mole)

	data eps/0.62197/

c   the next two lines contain a formula by herman wobus for the 
c   correction factor wfw for the departure of the mixture of air
c   and water vapor from the ideal gas law. the formula fits values
c   in table 89, p. 340 of the smithsonian meteorological tables,
c   but only for temperatures and pressures normally encountered in
c   in the atmosphere.

	x = 0.02*(t-12.5+7500./p)
	wfw = 1.+4.5e-06*p+1.4e-03*x*x
	fwesw = wfw*esw(t)
	r = eps*fwesw/(p-fwesw)

c   convert r from a dimensionless ratio to grams/kilogram.

	wmr = 1000.*r
	return
	end
