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
	function thm(t,p)

c	baker, schlatter  17-may-1982	  original version.

c   this function returns the wet-bulb potential temperature thm
c   (celsius) corresponding to a parcel of air saturated at 
c   temperature t (celsius) and pressure p (millibars).

	f(x) =   1.8199427e+01+x*( 2.1640800e-01+x*( 3.0716310e-04+x*
     1	       (-3.8953660e-06+x*( 1.9618200e-08+x*( 5.2935570e-11+x*
     2	       ( 7.3995950e-14+x*(-4.1983500e-17)))))))
	thm = t
	if (p.eq.1000.) return

c   compute the potential temperature (celsius).

	thd = (t+273.15)*(1000./p)**.286-273.15
	thm = thd+6.071*(exp(t/f(t))-exp(thd/f(thd)))
	return
	end
