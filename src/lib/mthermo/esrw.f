cdis   
cdis    Open Source License/Disclaimer, Forecast Systems Laboratory
cdis    NOAA/OAR/FSL, 325 Broadway Boulder, CO 80305
cdis    
cdis    This software is distributed under the Open Source Definition,
cdis    which may be found at http://www.opensource.org/osd.html.
cdis    
cdis    In particular, redistribution and use in source and binary forms,
cdis    with or without modification, are permitted provided that the
cdis    following conditions are met:
cdis    
cdis    - Redistributions of source code must retain this notice, this
cdis    list of conditions and the following disclaimer.
cdis    
cdis    - Redistributions in binary form must provide access to this
cdis    notice, this list of conditions and the following disclaimer, and
cdis    the underlying source code.
cdis    
cdis    - All modifications to this software must be clearly documented,
cdis    and are solely the responsibility of the agent making the
cdis    modifications.
cdis    
cdis    - If significant modifications or enhancements are made to this
cdis    software, the FSL Software Policy Manager
cdis    (softwaremgr@fsl.noaa.gov) should be notified.
cdis    
cdis    THIS SOFTWARE AND ITS DOCUMENTATION ARE IN THE PUBLIC DOMAIN
cdis    AND ARE FURNISHED "AS IS."  THE AUTHORS, THE UNITED STATES
cdis    GOVERNMENT, ITS INSTRUMENTALITIES, OFFICERS, EMPLOYEES, AND
cdis    AGENTS MAKE NO WARRANTY, EXPRESS OR IMPLIED, AS TO THE USEFULNESS
cdis    OF THE SOFTWARE AND DOCUMENTATION FOR ANY PURPOSE.  THEY ASSUME
cdis    NO RESPONSIBILITY (1) FOR THE USE OF THE SOFTWARE AND
cdis    DOCUMENTATION; OR (2) TO PROVIDE TECHNICAL SUPPORT TO USERS.
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
