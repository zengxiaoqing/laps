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
	function tcon(t,d)

c   this function returns the temperature tcon (celsius) at the lifting
c   condensation level, given the temperature t (celsius) and the
c   dew point d (celsius).

c	baker, schlatter  17-may-1982	  original version.

c
c   compute the dew point depression s.
	s = t-d
c   the approximation below, a third order polynomial in s and t,
c   is due to herman wobus. the source of data for fitting the
c   polynomial is unknown.

	dlt = s*(1.2185+1.278e-03*t+
     1        s*(-2.19e-03+1.173e-05*s-5.2e-06*t))
	tcon = t-dlt
	return
	end
