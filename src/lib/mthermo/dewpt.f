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
	function dewpt(ew)
c
c   this function yields the dew point dewpt (celsius), given the
c   water vapor pressure ew (millibars).
c
c	baker,schlatter	17-may-1982	original version
c
c   the empirical formula appears in bolton, david, 1980:
c   "the computation of equivalent potential temperature,"
c   monthly weather review, vol. 108, no. 7 (july), p. 1047, eq.(11).
c   the quoted accuracy is 0.03c or less for -35 < dewpt < 35c.
	enl = alog(ew)
	dewpt = (243.5*enl-440.8)/(19.48-enl)
	return
	end
