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
	function esrw(t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the saturation vapor pressure over liquid
c   water esrw (millibars) given the temperature t (celsius). the
c   formula used is due to richards, j.m., 1971: simple expression
c   for the saturation vapour pressure of water in the range -50 to
c   140c, british journal of applied physics, vol. 4, pp.l15-l18.
c   the formula was quoted more recently by wigley, t.m.l.,1974:
c   comments on 'a simple but accurate formula for the saturation
c   vapor pressure over liquid water,' journal of applied meteorology,
c   vol. 13, no. 5 (august) p.606.

	data cta,ts,ews/273.15,373.15,1013.25/

c   cta = difference between kelvin and celsius temperature
c   ts = temperature of the boiling point of water (k)
c   ews = saturation vapor pressure over liquid water at 100c

	data c1,     c2,     c3,     c4
     1	/ 13.3185,-1.9760,-0.6445,-0.1299 /
	tk = t+cta
	x = 1.-ts/tk
	px = x*(c1+x*(c2+x*(c3+c4*x)))
	vp = ews*exp(px)
	if (vp.lt.0) vp = 0.
	esrw = vp
	return
	end
