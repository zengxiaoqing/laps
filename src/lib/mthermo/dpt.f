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
	function dpt(ew)
c
c   this function returns the dew point dpt (celsius), given the
c   water vapor pressure ew (millibars).
c
c	baker,schlatter	17-may-1982	original version
c
	data es0/6.1078/
c   es0 = saturation vapor pressure (mb) over water at 0c
c   return a flag value if the vapor pressure is out of range.
	if (ew.gt..06.and.ew.lt.1013.) go to 5
	dpt = 9999.
	return
    5	continue
c   approximate dew point by means of teten's formula.
c   the formula appears as eq.(8) in bolton, david, 1980:
c   "the computation of equivalent potential temperature,"
c   monthly weather review, vol 108, no. 7 (july), p.1047.
c   the formula is ew(t) = es0*10**(7.5*t/(t+237.3))
c            or    ew(t) = es0*exp(17.269388*t/(t+237.3))
c   the inverse formula is used below.
	x = alog(ew/es0)
	dnm = 17.269388-x
	t = 237.3*x/dnm
	fac = 1./(ew*dnm)
c   loop for iterative improvement of the estimate of dew point
   10	continue
c   get the precise vapor pressure corresponding to t.
	edp = esw(t)
c   estimate the change in temperature corresponding to (ew-edp)
c   assume that the derivative of temperature with respect to 
c   vapor pressure (dtdew) is given by the derivative of the
c   inverse teten formula.
	dtdew = (t+237.3)*fac
	dt = dtdew*(ew-edp)
	t = t+dt
	if (abs(dt).gt.1.e-04) go to 10
	dpt = t
	return
	end
