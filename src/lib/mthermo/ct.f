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
	function ct(wbar,pc,ps)
c
c   this function returns the convective temperature ct (celsius)
c   given the mean mixing ratio wbar (g/kg) in the surface layer,
c   the pressure pc (mb) at the convective condensation level (ccl)
c   and the surface pressure ps (mb).
c
c	baker,schlatter	17-may-1982	original version
c
c   compute the temperature (celsius) at the ccl.
	tc= tmr(wbar,pc)
c   compute the potential temperature (celsius), i.e., the dry
c   adiabat ao through the ccl.
	ao= o(tc,pc)
c	   compute the surface temperature on the same dry adiabat ao.
	ct= tda(ao,ps)
	return
	end
