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
	function es(t)

c   this function returns the saturation vapor pressure es (mb) over
c   liquid water given the temperature t (celsius). the formula appears
c   in bolton, david, 1980: "the computation of equivalent potential
c   temperature," monthly weather review, vol. 108, no. 7 (july),
c   p. 1047, eq.(10). the quoted accuracy is 0.3% or better for
c   -35 < t < 35c.

c	baker, schlatter  17-may-1982	  original version.

c   es0 = saturation vapor pressure over liquid water at 0c

	data es0/6.1121/
	es = es0*exp(17.67*t/(t+243.5))
	return
	end
