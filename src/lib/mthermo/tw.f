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
	function tw(t,td,p)
c
c   this function returns the wet-bulb temperature tw (celsius)
c   given the temperature t (celsius), dew point td (celsius)
c   and pressure p (mb).  see p.13 in stipanuk (1973), referenced
c   above, for a description of the technique.
c
c	baker,schlatter	17-may-1982	original version
c
c   determine the mixing ratio line thru td and p.
	aw = w(td,p)
c
c   determine the dry adiabat thru t and p.
	ao = o(t,p)
	pi = p
c
c   iterate to locate pressure pi at the intersection of the two
c   curves .  pi has been set to p for the initial guess.
	do 4 i= 1,10
	   x= .02*(tmr(aw,pi)-tda(ao,pi))
	   if (abs(x).lt.0.01) go to 5
 4	   pi= pi*(2.**(x))
c   find the temperature on the dry adiabat ao at pressure pi.
 5	ti= tda(ao,pi)
c
c   the intersection has been located...now, find a saturation
c   adiabat thru this point. function os returns the equivalent 
c   potential temperature (c) of a parcel saturated at temperature
c   ti and pressure pi.
	aos= os(ti,pi)
c   function tsa returns the wet-bulb temperature (c) of a parcel at
c   pressure p whose equivalent potential temperature is aos.
	tw = tsa(aos,p)
	return
	end
