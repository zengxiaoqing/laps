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

	function esgg(t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the saturation vapor pressure over liquid
c   water esgg (millibars) given the temperature t (celsius). the
c   formula used, due to goff and gratch, appears on p. 350 of the
c   smithsonian meteorological tables, sixth revised edition, 1963,
c   by roland list.

	data cta,ews,ts/273.15,1013.246,373.15/

c   cta = difference between kelvin and celsius temperatures
c   ews = saturation vapor pressure (mb) over liquid water at 100c
c   ts = boiling point of water (k)

	data c1,      c2,      c3,      c4,       c5,       c6
     1	/ 7.90298, 5.02808, 1.3816e-7, 11.344, 8.1328e-3, 3.49149 /
	tk = t+cta

c   goff-gratch formula

	rhs = -c1*(ts/tk-1.)+c2*alog10(ts/tk)-c3*(10.**(c4*(1.-tk/ts))
     1	      -1.)+c5*(10.**(-c6*(ts/tk-1.))-1.)+alog10(ews)
	esw = 10.**rhs
	if (esw.lt.0.) esw = 0.
	esgg = esw
	return
	end
