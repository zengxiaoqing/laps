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

	function precpw(td,p,n)

c	baker, schlatter  17-may-1982	  original version.

c   this function computes total precipitable water precpw (cm) in a
c   vertical column of air based upon sounding data at n levels:
c	   td = dew point (celsius)
c	   p = pressure (millibars)
c   calculations are done in cgs units.

	dimension td(n),p(n)

c   g = acceleration due to the earth's gravity (cm/s**2)

	data g/980.616/

c   initialize value of precipitable water

	pw = 0.
	nl = n-1

c   calculate the mixing ratio at the lowest level.

	wbot = wmr(p(1),td(1))
	do 5 i=1,nl
	wtop = wmr(p(i+1),td(i+1))

c   calculate the layer-mean mixing ratio (g/kg).

	w = 0.5*(wtop+wbot)

c   make the mixing ratio dimensionless.

	wl = .001*w

c   calculate the specific humidity.

	ql = wl/(wl+1.)

c   the factor of 1000. below converts from millibars to dynes/cm**2.

	dp = 1000.*(p(i)-p(i+1))
	pw = pw+(ql/g)*dp
	wbot = wtop
    5	continue
	precpw = pw
	return
	end
