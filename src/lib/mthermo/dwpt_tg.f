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
	function dwpt(t,rh)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the dew point (celsius) given the temperature
c   (celsius) and relative humidity (%). the formula is used in the
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
