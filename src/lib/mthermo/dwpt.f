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

	function dwpt(t,rh)
c
c   this function returns the dew point (celsius) given the temperature
c   (celsius) and relative humidity (%). 
c
c	baker,schlatter	17-may-1982	original version
c
c	the formula is used in the
c   processing of u.s. rawinsonde data and is referenced in parry, h.
c   dean, 1969: "the semiautomatic computation of rawinsondes,"
c   technical memorandum wbtm edl 10, u.s. department of commerce,
c   environmental science services administration, weather bureau,
c   office of systems development, equipment development laboratory,
c   silver spring, md (october), page 9 and page ii-4, line 460.
	x = 1.-0.01*rh
c   compute dew point depression.
	dpd =(14.55+0.114*t)*x+((2.5+0.007*t)*x)**3+(15.9+0.117*t)*x**14
	dwpt = t-dpd
	return
	end
