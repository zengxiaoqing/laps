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
	function ssh(p,t)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns saturation specific humidity ssh (grams of 
c   water vapor per kilogram of moist air) given the pressure p
c   (millibars) and the temperature t (celsius). the equation is given
c   in standard meteorological texts. if t is dew point (celsius), then
c   ssh returns the actual specific humidity.
c   compute the dimensionless mixing ratio.

	w = .001*wmr(p,t)

c   compute the dimensionless saturation specific humidity.

	q = w/(1.+w)
	ssh = 1000.*q
	return
	end
