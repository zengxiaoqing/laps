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

	function esw(t)

c   this function returns the saturation vapor pressure esw (millibars)
c   over liquid water given the temperature t (celsius). the polynomial
c   approximation below is due to herman wobus, a mathematician who
c   worked at the navy weather research facility, norfolk, virginia,
c   but who is now retired. the coefficients of the polynomial were
c   chosen to fit the values in table 94 on pp. 351-353 of the smith-
c   sonian meteorological tables by roland list (6th edition). the
c   approximation is valid for -50 < t < 100c.

c	baker, schlatter  17-may-1982	  original version.

c
c   es0 = saturation vapor ressure over liquid water at 0c

	data es0/6.1078/
	pol = 0.99999683       + t*(-0.90826951e-02 +
     1	   t*(0.78736169e-04   + t*(-0.61117958e-06 +
     2     t*(0.43884187e-08   + t*(-0.29883885e-10 +
     3     t*(0.21874425e-12   + t*(-0.17892321e-14 +
     4     t*(0.11112018e-16   + t*(-0.30994571e-19)))))))))
	esw = es0/pol**8
	return
	end
