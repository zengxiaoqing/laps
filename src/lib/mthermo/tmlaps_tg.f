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
	function tmlaps(thetae,p)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the temperature tmlaps (celsius) at pressure
c   p (millibars) along the moist adiabat corresponding to an equivalent
c   potential temperature thetae (celsius).
c   the algorithm was written by eric smith at colorado state
c   university.

	data crit/0.1/
c   cta = difference between kelvin and celsius temperatures.
c   crit = convergence criterion (degrees kelvin)

	eq0 = thetae

c   initial guess for solution

	tlev = 25.

c   compute the saturation equivalent potential temperature correspon-
c   ding to temperature tlev and pressure p.

	eq1 = ept(tlev,tlev,p)
	dif = abs(eq1-eq0)
	if (dif.lt.crit) go to 3
	if (eq1.gt.eq0) go to 1

c   dt is the initial stepping increment.

	dt = 10.
	i = -1
	go to 2
    1	dt = -10.
	i = 1
    2	tlev = tlev+dt
	eq1 = ept(tlev,tlev,p)
	dif = abs(eq1-eq0)
	if (dif.lt.crit) go to 3
	j = -1
	if (eq1.gt.eq0) j=1
	if (i.eq.j) go to 2

c   the solution has been passed. reverse the direction of search
c   and decrease the stepping increment.

	tlev = tlev-dt
	dt = dt/10.
	go to 2
    3	tmlaps = tlev
	return
	end
