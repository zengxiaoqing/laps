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
