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
	function tlcl1(t,td)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the temperature tlcl1 (celsius) of the lifting
c   condensation level (lcl) given the initial temperature t (celsius)
c   and dew point td (celsius) of a parcel of air.
c   eric smith at colorado state university has used the formula
c   below, but its origin is unknown.

	data cta/273.15/

c   cta = difference between kelvin and celsius temperature

	tk = t+cta

c   compute the parcel vapor pressure (mb).
	es = eslo(td)
	tlcl = 2840./(3.5*alog(tk)-alog(es)-4.805)+55.
	tlcl1 = tlcl-cta
	return
	end
