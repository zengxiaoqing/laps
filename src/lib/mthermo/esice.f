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

	function esice(t)
c
c   this function returns the saturation vapor pressure with respect to
c   ice esice (millibars) given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula used is based upon the integration of the clausius-
c   clapeyron equation by goff and gratch.  the formula appears on p.350
c   of the smithsonian meteorological tables, sixth revised edition,
c   1963.
c
	data cta,eis/273.16,6.1071/
c   cta = difference between kelvin and celsius temperature
c   eis = saturation vapor pressure (mb) over a water-ice mixture at 0c
	data c1,c2,c3/9.09718,3.56654,0.876793/
c   c1,c2,c3 = empirical coefficients in the goff-gratch formula
	if (t.le.0.) go to 5
	esice = 99999.
	write(6,3)esice
    3	format(' saturation vapor pressure for ice cannot be computed',
     1	       /' for temperature > 0c. esice =',f7.0)
	return
    5	continue
c   freezing point of water (k)
	tf = cta
	tk = t+cta
c   goff-gratch formula
	rhs = -c1*(tf/tk-1.)-c2*alog10(tf/tk)+c3*(1.-tk/tf)+alog10(eis)
	esi = 10.**rhs
	if (esi.lt.0.) esi = 0.
	esice = esi
	return
	end
