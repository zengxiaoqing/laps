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
	function tlcl(t,td)
c
c   this function yields the temperature tlcl (celsius) of the lifting
c   condensation level, given the temperature t (celsius) and the
c   dew point td (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula used is n bolton, david,
c   1980: "the computation of equivalent potential temperature,"
c   monthly weather review, vol. 108, no. 7 (july), p. 1048, eq.(15).
c
c   convert from celsius to kelvin degrees.
	tk = t+273.16
	tdk = td+273.16
	a = 1./(tdk-56.)
	b = alog(tk/tdk)/800.
	tc = 1./(a+b)+56.
	tlcl = tc-273.16
	return
	end
