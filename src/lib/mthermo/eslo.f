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
	function eslo(t)
c
c   this function returns the saturation vapor pressure over liquid
c   water eslo (millibars) given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the formula is due to lowe, paul r.,1977: an approximating polynomial
c   for the computation of saturation vapor pressure, journal of applied
c   meteorology, vol 16, no. 1 (january), pp. 100-103.
c
c   the polynomial coefficients are a0 through a6.
	data a0,a1,a2,a3,a4,a5,a6
     1  /6.107799961,     4.436518521e-01, 1.428945805e-02,
     2   2.650648471e-04, 3.031240396e-06, 2.034080948e-08,
     3	 6.136820929e-11/
	es = a0+t*(a1+t*(a2+t*(a3+t*(a4+t*(a5+a6*t)))))
	if (es.lt.0.) es = 0.
	eslo = es
	return
	end
