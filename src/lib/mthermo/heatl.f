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

	function heatl(key,t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the latent heat of
c		evaporation/condensation         for key=1
c		melting/freezing		 for key=2
c		sublimation/deposition		 for key=3
c   for water. the latent heat heatl (joules per kilogram) is a 
c   function of temperature t (celsius). the formulas are polynomial
c   approximations to the values in table 92, p. 343 of the smithsonian
c   meteorological tables, sixth revised edition, 1963 by roland list.
c   the approximations were developed by eric smith at colorado state
c   university.
c   polynomial coefficients

	data a0,a1,a2/ 3337118.5,-3642.8583, 2.1263947/
	data b0,b1,b2/-1161004.0, 9002.2648,-12.931292/
	data c0,c1,c2/ 2632536.8, 1726.9659,-3.6248111/
	hltnt = 0.
	tk = t+273.15
	if (key.eq.1) hltnt=a0+a1*tk+a2*tk*tk
	if (key.eq.2) hltnt=b0+b1*tk+b2*tk*tk
	if (key.eq.3) hltnt=c0+c1*tk+c2*tk*tk
	heatl = hltnt
	return
	end
