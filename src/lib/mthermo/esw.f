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
	function esw(t)
c
c   this function returns the saturation vapor pressure esw (millibars)
c   over liquid water given the temperature t (celsius).
c
c	baker,schlatter	17-may-1982	original version
c
c   the polynomial approximation below is due to herman wobus, a mathematician who
c   worked at the navy weather research facility, norfolk, virginia,
c   but who is now retired. the coefficients of the polynomial were
c   chosen to fit the values in table 94 on pp. 351-353 of the smith-
c   sonian meteorological tables by roland list (6th edition). the
c   approximation is valid for -50 < t < 100c.
c
c   es0 = saturation vapor ressure over liquid water at 0c
	data es0/6.1078/
	pol = 0.99999683       + t*(-0.90826951e-02 +
     1	   t*(0.78736169e-04   + t*(-0.61117958e-06 +
     2     t*(0.43884187e-08   + t*(-0.29883885e-10 +
     3     t*(0.21874425e-12   + t*(-0.17892321e-14 +
     4     t*(0.11112018e-16   + t*(-0.30994571e-19)))))))))
	esw = es0/pol**8
	return
	end
