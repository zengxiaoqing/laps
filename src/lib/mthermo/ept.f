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

	function ept(t,td,p)
c
c   this function returns the equivalent potential temperature ept
c   (celsius) for a parcel of air initially at temperature t (celsius),
c   dew point td (celsius) and pressure p (millibars).
c
c	baker,schlatter	17-may-1982	original version
c
c	the formula used
c   is eq.(43) in bolton, david, 1980: "the computation of equivalent
c   potential temperature," monthly weather review, vol. 108, no. 7
c   (july), pp. 1046-1053. the maximum error in ept in 0.3c.  in most
c   cases the error is less than 0.1c.
c
c   compute the mixing ratio (grams of water vapor per kilogram of
c   dry air).
	w = wmr(p,td)
c   compute the temperature (celsius) at the lifting condensation level.
	tlcl = tcon(t,td)
	tk = t+273.16
	tl = tlcl+273.16
	pt = tk*(1000./p)**(0.2854*(1.-0.00028*w))
	eptk = pt*exp((3.376/tl-0.00254)*w*(1.+0.00081*w))
	ept= eptk-273.16
	return
	end
